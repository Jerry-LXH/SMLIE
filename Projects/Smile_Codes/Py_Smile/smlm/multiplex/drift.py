"""Estimate one registered drift and apply it to every channel.

Conventions
-----------
* ``DriftTrace`` values are cumulative sample drift **in each channel's own
  pixel frame**. The same physical motion therefore appears with slightly
  different ``(dx, dy)`` per channel; compare or combine channels only after
  mapping through the registration (``DualViewGeometry``).
* Every trace is normalized so that frame 0 is zero; frame 0 is the drift
  reference frame in which crop windows are expressed.
* Images are only translated (no resampling of channel distortion); channel
  geometry stays in the coordinate layer.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

from ..datatypes import DriftTrace, Localizations, MultiImageStack
from ..imgprocess import DriftCorrector
from .coords import DualViewGeometry, map_channel_xy
from .emitters import channel_locs, channel_names


_CONCAT_FIELDS = ('frame', 'x', 'y', 'photons', 'bg', 'sigma', 'lp', 'z',
                  'z_lp')
_SIDES = ('top', 'bottom', 'left', 'right')


# ---------------------------------------------------------------------------
#  Small helpers
# ---------------------------------------------------------------------------
def _pixel_count(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(
            value, (int, np.integer)):
        raise TypeError(f'{name} must be an integer pixel count')
    if value < 0:
        raise ValueError(f'{name} must be non-negative')
    return int(value)


def _image_shape(shape):
    try:
        height, width = (int(v) for v in tuple(shape)[-2:])
    except (TypeError, ValueError):
        raise TypeError('image_shape must end with two integers '
                        '(height, width)') from None
    if height <= 0 or width <= 0:
        raise ValueError(f'invalid image shape {shape}')
    return height, width


def _drift_traces(drift):
    if not isinstance(drift, Mapping):
        raise TypeError('drift must be a {channel: DriftTrace} mapping, '
                        f'got {type(drift).__name__}')
    if not drift:
        raise ValueError('drift mapping is empty')
    for channel, trace in drift.items():
        if not isinstance(trace, DriftTrace):
            raise TypeError(f'drift[{channel!r}] must be a DriftTrace')
    return drift


def _channel_names(channels, available, name):
    return channel_names(available, channels, name=name)


def _normalize_locs(multi, locs):
    return channel_locs(multi, locs, allow_missing=True)


def _validate_locs(locs, channel, n_frames):
    if not isinstance(locs, Localizations):
        raise TypeError(f'locs[{channel!r}] must be Localizations')
    for name in ('x', 'y'):
        values = np.asarray(getattr(locs, name))
        if values.ndim != 1 or len(values) != len(locs):
            raise ValueError(f'locs[{channel!r}].{name} must be a 1D array '
                             'matching the localization count')
        if not np.all(np.isfinite(values)):
            raise ValueError(f'locs[{channel!r}].{name} contains non-finite '
                             'coordinates')
    frames = np.asarray(locs.frame)
    if frames.ndim != 1 or len(frames) != len(locs):
        raise ValueError(f'locs[{channel!r}].frame must be a 1D array '
                         'matching the localization count')
    if frames.size:
        if not np.all(np.isfinite(frames)):
            raise ValueError(f'locs[{channel!r}] contains non-finite frames')
        frame_idx = frames.astype(int)
        if not np.array_equal(frames, frame_idx):
            raise ValueError(f'locs[{channel!r}].frame must contain integers')
        if frame_idx.min() < 0 or frame_idx.max() >= n_frames:
            raise ValueError(f'locs[{channel!r}] frame indices must be in '
                             f'[0, {n_frames})')


def _validate_geometry(multi, geom):
    missing = [channel for channel in multi.channels
               if channel not in geom.windows]
    if missing:
        raise ValueError(f'geometry is missing image channels {missing}')
    if geom.reference not in multi.channels:
        raise ValueError(f'geometry reference {geom.reference!r} is not an '
                         'image channel')
    expected_shape = multi.data.shape[-2:]
    mismatched = {
        channel: geom.windows[channel][2:]
        for channel in multi.channels
        if geom.windows[channel][2:] != expected_shape
    }
    if mismatched:
        raise ValueError(f'geometry window sizes must match image shape '
                         f'{expected_shape}: {mismatched}')


def _map_channel_xy(xy, geom, source, target):
    """Deprecated shim: see :func:`smlm.multiplex.coords.map_channel_xy`."""
    return map_channel_xy(xy, geom, source, target)


def _channel_to_reference_jacobian(geom, channel, reference):
    """Linearize channel-to-reference coordinates at the window center."""
    if channel == reference:
        return np.eye(2)
    _, _, height, width = geom.windows[channel]
    center = np.array([(width - 1) / 2.0, (height - 1) / 2.0])
    step = 1e-3
    jacobian = np.empty((2, 2), dtype=float)
    for axis in range(2):
        delta = np.zeros(2, dtype=float)
        delta[axis] = step
        mapped = _map_channel_xy(
            np.vstack((center + delta, center - delta)), geom, channel,
            reference)
        jacobian[:, axis] = (mapped[0] - mapped[1]) / (2.0 * step)
    if not np.all(np.isfinite(jacobian)):
        raise ValueError(f'non-finite registration Jacobian for '
                         f'{channel!r} -> {reference!r}')
    if abs(float(np.linalg.det(jacobian))) < 1e-8:
        raise ValueError(f'singular registration Jacobian for '
                         f'{channel!r} -> {reference!r}')
    return jacobian


# ---------------------------------------------------------------------------
#  Registered point-cloud merge
# ---------------------------------------------------------------------------
def merge_registered(locs, geom: DualViewGeometry, *, reference='g',
                     channels=None, verbose=False):
    """Map selected channel localizations into ``reference`` and concatenate.

    Returns the merged localization set and one boolean provenance mask per
    selected channel. ``channels`` must include ``reference``.
    """
    available = list(geom.windows)
    selected = _channel_names(channels, available, 'channels')
    if reference not in selected:
        raise ValueError(f'estimation channels must include reference '
                         f'{reference!r}')
    ordered = [reference] + [channel for channel in selected
                             if channel != reference]

    chunks = []
    for channel in ordered:
        if channel not in locs or locs[channel] is None:
            raise ValueError(f'missing localizations for selected channel '
                             f'{channel!r}')
        current = locs[channel]
        xy = _map_channel_xy(np.column_stack((current.x, current.y)),
                             geom, channel, reference)
        chunks.append((channel, current, xy))

    kw = {}
    for name in _CONCAT_FIELDS:
        values = [
            (xy[:, 0] if name == 'x' else
             xy[:, 1] if name == 'y' else
             getattr(current, name))
            for _, current, xy in chunks
        ]
        kw[name] = (None if any(value is None for value in values) else
                    np.concatenate([np.asarray(value) for value in values]))

    merged = Localizations(**kw)

    parts = {}
    offset = 0
    for channel, current, _ in chunks:
        mask = np.zeros(len(merged), dtype=bool)
        mask[offset:offset + len(current)] = True
        parts[channel] = mask
        offset += len(current)
    if verbose:
        counts = ', '.join(f'{channel}:{int(mask.sum())}'
                           for channel, mask in parts.items())
        print(f'  merged registered cloud: {len(merged)} localizations '
              f'({counts}) in {reference!r} coordinates')
    return merged, parts


# ---------------------------------------------------------------------------
#  Estimation
# ---------------------------------------------------------------------------
def _default_uncertainty_px(locs, channels, geom, reference):
    """Median ``lp`` mapped from each channel's pixels into reference pixels.

    ``Localizations.lp`` carries the CRLB in **pixel** units (the Fisher matrix
    is assembled in pixel coordinates), so no nm/pixel_size conversion applies.
    """
    values = []
    for channel in channels:
        current = locs[channel]
        if current.lp is None:
            continue
        lp = np.asarray(current.lp, dtype=float)
        if lp.ndim != 1 or len(lp) != len(current):
            raise ValueError(f'locs[{channel!r}].lp must match the '
                             'localization count')
        valid = np.isfinite(lp) & (lp > 0)
        if not np.any(valid):
            continue
        jacobian = _channel_to_reference_jacobian(
            geom, channel, reference)
        scale = np.sqrt(np.sum(jacobian ** 2) / 2.0)
        values.append(lp[valid] * scale)
    if not values:
        raise ValueError('cannot infer drift uncertainty: selected '
                         'localizations need finite positive lp values; '
                         'otherwise pass uncertainty_px explicitly')
    return float(np.median(np.concatenate(values)))


def estimate_drift_multi(multi, locs, *, drift_reference=None,
                         estimate_channels=None, frames_per_seg=10,
                         do_stage2=True, uncertainty_px=None,
                         matched_pairs_threshold=30, r_scale=4.5,
                         max_iter=20, geom=None, verbose=True):
    """Estimate one registered drift and convert it to every channel.

    Parameters
    ----------
    multi : MultiImageStack
        Supplies channels, frame count, pixel size and (by default) geometry.
    locs : dict or sequence
        ``{channel: Localizations}`` in each channel's own analysis-window
        coordinates; a sequence ordered like ``multi.channels`` also works.
    drift_reference : str
        Channel whose pixel frame defines the estimated drift. Distinct from
        ``coord_reference`` (``geom.reference``), which defines registration and
        stays attached to the returned image stack.
    estimate_channels : sequence of str
        Channels contributing localizations to the merged cloud. Must include
        ``drift_reference``; other channels still receive a converted trace.

    Returns
    -------
    dict
        ``{channel: DriftTrace}`` for **every** image channel. Traces are in
        each channel's own pixel frame and normalized so frame 0 is zero.
        Provenance (drift/coord reference, estimation channels, uncertainty)
        lives in ``trace.metadata``.
    """
    if not isinstance(multi, MultiImageStack):
        raise TypeError('multi must be a MultiImageStack')
    if multi.n_frames <= 0:
        raise ValueError('multi must contain at least one frame')
    locs = _normalize_locs(multi, locs)
    geom = (multi.geometry(check_coverage=False) if geom is None else geom)
    if not isinstance(geom, DualViewGeometry):
        raise TypeError('geom must be a DualViewGeometry')
    _validate_geometry(multi, geom)

    drift_reference = (multi.ref_name if drift_reference is None
                       else drift_reference)
    if drift_reference not in multi.channels:
        raise ValueError(f'unknown drift reference channel '
                         f'{drift_reference!r}; available: {multi.channels}')
    channels = _channel_names(estimate_channels, multi.channels,
                              'estimate_channels')
    if drift_reference not in channels:
        raise ValueError('estimate_channels must include the drift reference '
                         f'{drift_reference!r}')

    for channel in channels:
        if locs[channel] is None:
            raise ValueError(f'missing localizations for channel {channel!r}')
        _validate_locs(locs[channel], channel, multi.n_frames)
    merged, _parts = merge_registered(
        locs, geom, reference=drift_reference, channels=channels,
        verbose=verbose)
    if len(merged) == 0:
        raise ValueError('no localizations available for drift estimation')

    uncertainty = (float(uncertainty_px) if uncertainty_px is not None
                   else _default_uncertainty_px(
                       locs, channels, geom, drift_reference))
    if not np.isfinite(uncertainty) or uncertainty <= 0:
        raise ValueError('uncertainty_px must be a positive finite value')

    from ..drift import CloudDriftEstimator
    estimator = CloudDriftEstimator(
        frames_per_seg=frames_per_seg,
        uncertainty=uncertainty,
        matched_pairs_threshold=matched_pairs_threshold,
        do_stage2=do_stage2,
        r_scale=r_scale,
        max_iter=max_iter,
        verbose=False,
    )
    estimated = estimator.estimate(merged, n_frames=multi.n_frames)
    if estimated.n_frames != multi.n_frames or \
            len(estimated.dy) != multi.n_frames:
        raise ValueError(f'estimator returned an invalid drift length for a '
                         f'{multi.n_frames}-frame stack')
    dx_ref = np.asarray(estimated.dx, dtype=float)
    dy_ref = np.asarray(estimated.dy, dtype=float)
    if not np.all(np.isfinite(dx_ref)) or not np.all(np.isfinite(dy_ref)):
        raise ValueError('drift estimator returned non-finite values')
    d_ref = np.column_stack((dx_ref - dx_ref[0], dy_ref - dy_ref[0]))

    provenance = {
        'estimator': 'registered-multi',
        'drift_reference': drift_reference,
        'coord_reference': geom.reference,
        'estimate_channels': list(channels),
        'uncertainty_px': uncertainty,
    }
    drift = {}
    for channel in multi.channels:
        jacobian = _channel_to_reference_jacobian(
            geom, channel, drift_reference)
        try:
            d_channel = np.linalg.solve(jacobian, d_ref.T).T
        except np.linalg.LinAlgError as exc:
            raise ValueError(f'cannot convert drift from {drift_reference!r} '
                             f'to {channel!r}') from exc
        drift[channel] = DriftTrace(
            dx=d_channel[:, 0],
            dy=d_channel[:, 1],
            pixel_size=multi.pixel_size,
            metadata={**provenance, 'channel': channel,
                      'frame': f'{channel}-analysis'},
        )
    if verbose:
        reference_total = float(np.hypot(*d_ref[-1]))
        summary = ', '.join(f'{channel}:{trace.total_drift:.2f}px'
                            for channel, trace in drift.items())
        print(f'  estimated {drift_reference!r} drift: '
              f'{reference_total:.2f} px from {channels}')
        print(f'  converted drift to channels: {summary}')
    return drift


# ---------------------------------------------------------------------------
#  Crop window
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CropWindow:
    """Crop window in frame-0 (drift reference) coordinates.

    ``y0``/``x0`` are the window origin in the corrected frame-0 frame, so the
    kept region is ``[y0, y0 + height) x [x0, x0 + width)``. ``margins`` records
    the pixels lost on each side to drift, padding and trimming, and
    ``requirements`` keeps the per-channel drift minimum for diagnostics.
    """
    y0: int
    x0: int
    height: int
    width: int
    margins: dict = field(default_factory=dict)
    requirements: dict = field(default_factory=dict)

    def __post_init__(self):
        for name in ('y0', 'x0', 'height', 'width'):
            object.__setattr__(self, name,
                               _pixel_count(getattr(self, name), name))
        if self.height <= 0 or self.width <= 0:
            raise ValueError('crop window must be non-empty')

    @property
    def shape(self):
        return (self.height, self.width)

    @property
    def y(self):
        return (self.y0, self.y0 + self.height)

    @property
    def x(self):
        return (self.x0, self.x0 + self.width)

    @property
    def slices(self):
        return (slice(*self.y), slice(*self.x))

    def __repr__(self):
        margins = ', '.join(f'{side}={value}'
                            for side, value in self.margins.items())
        return (f'CropWindow(y0={self.y0}, x0={self.x0}, '
                f'height={self.height}, width={self.width}'
                + (f'; margins: {margins}' if margins else '') + ')')


def drift_side_needs(trace):
    """Minimum pixels a single trace loses on each side of the frame.

    Each side uses ``max(0, ...)`` of that trace's drift range, so a
    one-directional drift does not waste the opposite edge.
    """
    if not isinstance(trace, DriftTrace):
        raise TypeError('trace must be a DriftTrace')
    dx = np.asarray(trace.dx, dtype=float)
    dy = np.asarray(trace.dy, dtype=float)
    if dx.ndim != 1 or dy.ndim != 1 or dx.size != dy.size or dx.size == 0:
        raise ValueError('drift trace must have matching non-empty dx and dy '
                         'arrays')
    if not np.all(np.isfinite(dx)) or not np.all(np.isfinite(dy)):
        raise ValueError('drift trace contains non-finite values')
    return {
        'top': int(np.ceil(max(0.0, -float(dy.min())))),
        'bottom': int(np.ceil(max(0.0, float(dy.max())))),
        'left': int(np.ceil(max(0.0, -float(dx.min())))),
        'right': int(np.ceil(max(0.0, float(dx.max())))),
    }


def _per_side(value, name):
    """Normalize ``int | {side: px} | 4-sequence`` into a per-side dict."""
    if isinstance(value, Mapping):
        unknown = [key for key in value if key not in _SIDES]
        if unknown:
            raise ValueError(f'{name} has unknown sides {unknown}; '
                             f'use {list(_SIDES)}')
        return {side: _pixel_count(value.get(side, 0), f'{name}[{side!r}]')
                for side in _SIDES}
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        count = _pixel_count(value, name)
        return {side: count for side in _SIDES}
    values = list(value)
    if len(values) != len(_SIDES):
        raise ValueError(f'{name} must be an int, a {{side: px}} mapping, or '
                         f'4 values in the order {list(_SIDES)}')
    return {side: _pixel_count(item, f'{name}[{side!r}]')
            for side, item in zip(_SIDES, values)}


def drift_requirements(drift, image_shape):
    """Per-channel minimum pixels lost on each side, plus the combined need."""
    height, width = _image_shape(image_shape)
    per_channel = {channel: drift_side_needs(trace)
                   for channel, trace in _drift_traces(drift).items()}
    combined = {side: max(need[side] for need in per_channel.values())
                for side in _SIDES}
    if combined['top'] + combined['bottom'] >= height or \
            combined['left'] + combined['right'] >= width:
        raise ValueError(
            f'drift needs margins {combined} for a {height}x{width} window; '
            'no common area remains')
    return per_channel, combined


def drift_crop_window(drift, image_shape, *, pad=0, extra_trim=0):
    """Smallest common window that keeps every channel's drifted frames valid.

    ``pad`` covers the subpixel interpolation boundary and ``extra_trim`` is an
    additional trim. Both accept an int applied to all sides, a ``{side: px}``
    mapping, or 4 values in the order ``(top, bottom, left, right)``. The drift
    requirement is the maximum across channels, so every channel keeps the same
    window. The result is expressed in frame-0 coordinates.
    """
    height, width = _image_shape(image_shape)
    pad = _per_side(pad, 'pad')
    extra = _per_side(extra_trim, 'extra_trim')
    per_channel, need = drift_requirements(drift, (height, width))
    margins = {side: need[side] + pad[side] + extra[side] for side in _SIDES}
    out_height = height - margins['top'] - margins['bottom']
    out_width = width - margins['left'] - margins['right']
    if out_height <= 0 or out_width <= 0:
        raise ValueError(
            f'margins {margins} (drift {need} + pad {pad} + trim {extra}) '
            f'leave no valid area in a {height}x{width} window')
    return CropWindow(
        y0=margins['top'], x0=margins['left'],
        height=out_height, width=out_width,
        margins=margins, requirements=per_channel,
    )


def _as_window(window):
    if isinstance(window, CropWindow):
        return window
    try:
        values = tuple(window)
    except TypeError:
        raise TypeError('window must be a CropWindow or '
                        '(y0, x0, height, width)') from None
    if len(values) != 4:
        raise ValueError('window must be (y0, x0, height, width)')
    return CropWindow(y0=values[0], x0=values[1],
                      height=values[2], width=values[3])


def _window_margins(window, shape):
    height, width = shape
    return {
        'top': window.y0,
        'left': window.x0,
        'bottom': height - window.y0 - window.height,
        'right': width - window.x0 - window.width,
    }


def _window_in_bounds(window, need, shape):
    """Ensure the window fits the region left by every channel's drift."""
    margins = _window_margins(window, shape)
    if min(margins.values()) < 0:
        raise ValueError(
            f'window {window.shape} at ({window.y0}, {window.x0}) does not fit '
            f'a {shape[0]}x{shape[1]} image')
    too_small = {side: (margins[side], need[side])
                 for side in _SIDES if margins[side] < need[side]}
    if too_small:
        detail = ', '.join(f'{side}: have {have} < need {required}'
                           for side, (have, required) in too_small.items())
        raise ValueError(
            f'window is too large for the drift ({detail}); build it with '
            'drift_crop_window(...) or trim more pixels')
    return margins


# ---------------------------------------------------------------------------
#  Correction
# ---------------------------------------------------------------------------
def correct_multi(multi, locs, drift, *, window, geom=None, edge=5,
                  interpolation_order=1, fill=0.0, verbose=True):
    """Translate every channel, then crop to one explicit common window.

    ``window`` is required and expressed in frame-0 coordinates; build it with
    :func:`drift_crop_window`. It is validated against every channel's drift so
    interpolation borders cannot silently enter the result.

    Returns ``(multi_corr, locs_corr)`` where ``locs_corr`` maps **every** image
    channel (including channels that did not contribute to estimation) to its
    corrected, cropped localizations.
    """
    if not isinstance(multi, MultiImageStack):
        raise TypeError('multi must be a MultiImageStack')
    drift = _drift_traces(drift)
    locs = _normalize_locs(multi, locs)
    geom = multi.geometry(check_coverage=False) if geom is None else geom
    if not isinstance(geom, DualViewGeometry):
        raise TypeError('geom must be a DualViewGeometry')
    _validate_geometry(multi, geom)

    missing = [channel for channel in multi.channels if channel not in drift]
    if missing:
        raise ValueError(f'missing drift traces for channels {missing}')
    for channel in multi.channels:
        trace = drift[channel]
        if trace.n_frames != multi.n_frames or \
                len(np.asarray(trace.dy)) != multi.n_frames:
            raise ValueError(f'{channel!r} drift has {trace.n_frames} frames, '
                             f'expected {multi.n_frames}')

    edge = _pixel_count(edge, 'edge')
    interpolation_order = _pixel_count(
        interpolation_order, 'interpolation_order')
    if interpolation_order > 5:
        raise ValueError('interpolation_order must be between 0 and 5')
    if not np.isfinite(fill):
        raise ValueError('fill must be finite')
    fill = float(fill)

    window = _as_window(window)
    height, width = multi.data.shape[-2:]
    channel_drift = {channel: drift[channel] for channel in multi.channels}
    _, need = drift_requirements(channel_drift, (height, width))
    margins = _window_in_bounds(window, need, (height, width))

    images = {}
    out_locs = {}
    new_windows = {}
    for channel in multi.channels:
        trace = drift[channel]
        current_locs = locs[channel]
        if current_locs is not None:
            _validate_locs(current_locs, channel, multi.n_frames)

        corrector = DriftCorrector(
            edge=edge, fill=fill, order=interpolation_order)
        image, crop = corrector.correct_image(
            multi.channel(channel), trace,
            window=(window.y0, window.x0, window.height, window.width))
        corrected_locs = None
        if current_locs is not None:
            corrected_locs, _ = corrector.correct_locs(
                current_locs, trace, window=crop)
        images[channel] = image
        out_locs[channel] = corrected_locs

        y0, x0, _, _ = geom.windows[channel]
        r0, c0, crop_height, crop_width = crop
        new_windows[channel] = (y0 + r0, x0 + c0, crop_height, crop_width)
        if verbose:
            kept = ('n/a' if corrected_locs is None
                    else f'{len(corrected_locs)}/{len(current_locs)}')
            print(f'  [{channel}] locs kept {kept}; '
                  f'total drift {trace.total_drift:.2f} px')

    data = np.stack(
        [np.asarray(images[channel].data) for channel in multi.channels],
        axis=0)
    offsets = np.asarray(
        [new_windows[channel][:2] for channel in multi.channels], dtype=int)
    provenance = drift[multi.channels[0]].metadata
    metadata = dict(multi.metadata)
    metadata['drift_crop'] = {
        'window': [window.y0, window.x0, window.height, window.width],
        'margins': {side: int(value) for side, value in margins.items()},
        'drift_requirements': need,
        'edge': edge,
        'interpolation_order': interpolation_order,
        'fill': fill,
        'coord_reference': geom.reference,
        'drift_reference': provenance.get('drift_reference'),
        'estimate_channels': list(provenance.get('estimate_channels', [])),
    }
    multi_corr = MultiImageStack(
        data=data,
        channels=list(multi.channels),
        reference=multi.index(geom.reference),
        reg=dict(geom.calibs),
        offsets=offsets,
        pixel_size=multi.pixel_size,
        metadata=metadata,
    )
    if verbose:
        print(f'  crop window: {window}')
        print(f'  corrected stack: {multi_corr}')
    return multi_corr, out_locs
