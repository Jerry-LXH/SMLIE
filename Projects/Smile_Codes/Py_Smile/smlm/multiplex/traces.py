"""Extract traces in several channels, per channel choosing the emitter source.

Why the sources are per channel
-------------------------------
A trace needs two sets of positions in the **output** channel's coordinates:

* the **aperture centres** (``centers_for_raw``) - where to integrate. Usually all
  channels use one physical reference channel's emitters (mapped through the
  registration), so both channels look at the same molecule even if one channel
  did not detect it in some frames;
* the **background mask centres** (``centers_for_bg``) - whose cores are cut out
  of the annulus. Here each channel normally wants **its own** emitters, because
  a molecule that only blinks in this channel never appears in the other
  channel's emitter list; leaving it in the ring inflates the background and
  biases the trace low.

Both are given as ``{output_channel: source_channel}``; a channel that does not
appear in ``centers_from`` produces no trace. Centres are exchanged in
``[x, y]`` order and mapped by :func:`smlm.multiplex.coords.map_channel_xy`.

Everything below :class:`smlm.trace.fill.TraceFillPipeline` is shared with the
single-channel flow, so the traces here carry the same fields and can be fed
straight into :class:`~smlm.trace.states.ThresholdStates` /
:class:`~smlm.trace.states.ChangepointStates`.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from ..datatypes import MultiImageStack
from ..trace.fill import TraceFillPipeline
from ..trace.signal import ApertureRaw
from ..trace.annulus import AnnulusRing
from ..trace.calculate import IntensityFromRawBg, PoissonSNR
from .coords import map_channel_xy
from .emitters import emitter_centers_xy


def _source_mapping(mapping, name, output_channels, available):
    if not isinstance(mapping, Mapping):
        raise TypeError(f'{name} must be a {{output_channel: source_channel}} '
                        'mapping')
    if not mapping:
        raise ValueError(f'{name} must contain at least one channel')
    unknown = [channel for channel in mapping if channel not in output_channels]
    if unknown:
        raise ValueError(f'{name} has unknown output channels {unknown}; '
                         f'available: {list(output_channels)}')
    bad_source = [channel for channel in mapping.values()
                  if channel not in available]
    if bad_source:
        raise ValueError(f'{name} refers to unknown source channels '
                         f'{bad_source}; available: {list(available)}')
    return dict(mapping)


def _centers_from(emitters, source):
    if source not in emitters:
        raise ValueError(f'emitters has no entry for source channel {source!r}; '
                         f'available: {list(emitters)}')
    xy = emitter_centers_xy(emitters[source])
    return np.asarray(xy, dtype=float).reshape(-1, 2)


def _per_channel(value, channels, name):
    """标量 -> 每通道同一个值; dict -> 按通道取值 (缺的通道报错)。"""
    if isinstance(value, Mapping):
        missing = [c for c in channels if c not in value]
        if missing:
            raise ValueError(f'{name} is missing channels {missing}')
        return {c: value[c] for c in channels}
    return {c: value for c in channels}


def _xy_to_yx(xy):
    xy = np.asarray(xy, float).reshape(-1, 2)
    return np.column_stack([xy[:, 1], xy[:, 0]])


def extract_multichannel_traces(multi, emitters, *, centers_from,
                                mask_from=None, sigma=1.5, n_frames=None,
                                mode='matched', radius=None, k_sigma=2.0,
                                bg_gap=2, bg_width=4, excl_radius=2,
                                bg_ref='frames', bg_smooth=0,
                                add_snr=True, verbose=False):
    """Extract one trace per aperture centre in every requested channel.

    Parameters
    ----------
    multi : MultiImageStack
        Drift-corrected stack; supplies images, channels and the geometry used
        to map emitter coordinates between channels.
    emitters : dict
        ``{channel: [Emitter]}`` from
        :func:`smlm.multiplex.analyze_multichannel_emitters` (``(N, 2) [x, y]``
        arrays are also accepted).
    centers_from : dict
        ``{output_channel: source_channel}`` aperture-centre source. Its keys
        define which channels get traces; other channels are skipped.
    mask_from : dict, optional
        ``{output_channel: source_channel}`` background-mask source. Defaults to
        ``centers_from`` (aperture centres mask themselves).
    sigma : float | dict
        PSF σ in px, scalar or ``{channel: px}``. Required for ``mode='matched'``
        (the weight kernel is built from it); with ``mode='hard'`` it only goes
        into the ``sigma`` field. **Not** estimated from localizations - in a
        dense field that estimate mixes in neighbours.
    n_frames : int, optional
        Frames per trace (default: all).
    mode : 'matched' | 'hard'
        ``matched`` needs an accurate σ; ``hard`` is a fixed disc and therefore
        the safer choice when the fitted σ is not trustworthy.
    radius : int, optional
        Fixed aperture radius (overrides ``k_sigma·σ``).
    k_sigma, bg_gap, bg_width, excl_radius
        Passed to :class:`~smlm.trace.signal.ApertureRaw` and
        :class:`~smlm.trace.annulus.AnnulusRing`.
    bg_ref, bg_smooth
        Background reference mode and centred moving-average width in frames
        (``bg_smooth=0`` = no smoothing).
    add_snr : bool
        Also write ``snr`` (needs ``intensity``; cheap).

    Returns
    -------
    dict
        ``{output_channel: [Trace]}`` only for channels listed in
        ``centers_from``; ``emitter_id`` follows the source emitter order.
    """
    if not isinstance(multi, MultiImageStack):
        raise TypeError('multi must be a MultiImageStack')
    if not isinstance(emitters, Mapping):
        raise TypeError('emitters must be a {channel: [Emitter]} mapping')

    centers_from = _source_mapping(centers_from, 'centers_from',
                                   multi.channels, multi.channels)
    if mask_from is None:
        mask_from = dict(centers_from)
    else:
        mask_from = _source_mapping(mask_from, 'mask_from', multi.channels,
                                    multi.channels)
    mask_from = {channel: mask_from.get(channel, centers_from[channel])
                 for channel in centers_from}

    n_frames = multi.n_frames if n_frames is None else int(n_frames)
    geom = multi.geometry(check_coverage=False)
    sigmas = _per_channel(sigma, list(centers_from), 'sigma')

    fillers = [ApertureRaw(mode=mode, radius=radius, k_sigma=k_sigma),
               AnnulusRing(gap=bg_gap, width=bg_width,
                           excl_radius=excl_radius, bg_ref=bg_ref,
                           bg_smooth=bg_smooth),
               IntensityFromRawBg()]
    if add_snr:
        fillers.append(PoissonSNR())

    traces = {}
    for channel, center_source in centers_from.items():
        center_xy = _centers_from(emitters, center_source)
        centers_xy = map_channel_xy(center_xy, geom, center_source, channel)

        mask_xy = map_channel_xy(
            _centers_from(emitters, mask_from[channel]), geom,
            mask_from[channel], channel)

        pipe = TraceFillPipeline(
            multi.data[multi.index(channel)],
            _xy_to_yx(centers_xy),
            centers_for_bg=_xy_to_yx(mask_xy),
            sigma=sigmas[channel],
            n_frames=n_frames,
            tag=channel,
        )
        traces[channel] = pipe.run(fillers, verbose=verbose)
        if verbose:
            print(f'  [{channel}] centres from {center_source!r}, '
                  f'background masked by {mask_from[channel]!r}: '
                  f'{len(centers_xy)} traces')
    return traces
