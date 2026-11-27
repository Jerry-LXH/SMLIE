"""Per-channel emitter clustering for multichannel data.

Clustering is done **inside each channel**, on that channel's own localizations
and in that channel's own coordinates. The result is a plain
``{channel: [Emitter]}`` mapping that later steps consume:

* one channel's emitters provide the physical aperture centres, mapped through
  the registration by
  :func:`smlm.multiplex.traces.extract_multichannel_traces`;
* every channel's own emitters provide the neighbour positions used to mask the
  background ring, because a molecule that only blinks in one channel never
  appears in another channel's emitter list.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from ..datatypes import MultiImageStack
from ..emitter import analyze_emitters_sm


def channel_names(available, channels=None, *, name='channels', required=()):
    """Resolve a channel selection against ``available``, keeping order."""
    if channels is None:
        names = list(available)
    elif isinstance(channels, str):
        names = [channels]
    else:
        names = list(channels)
    if not names:
        raise ValueError(f'{name} must contain at least one channel')
    if len(set(names)) != len(names):
        raise ValueError(f'{name} contains duplicate channels')
    unknown = [channel for channel in names if channel not in available]
    if unknown:
        raise ValueError(f'{name} contains unknown channels {unknown}; '
                         f'available: {list(available)}')
    missing = [channel for channel in required if channel not in names]
    if missing:
        raise ValueError(f'{name} must include {missing}')
    return names


def channel_locs(multi, locs, *, allow_missing=False):
    """Normalize ``locs`` to ``{channel: Localizations}`` in channel order.

    ``locs`` may be a mapping or a sequence aligned with ``multi.channels``.
    Missing channels become ``None`` when ``allow_missing`` is set, otherwise
    they raise.
    """
    if isinstance(locs, Mapping):
        unknown = [channel for channel in locs if channel not in multi.channels]
        if unknown:
            raise ValueError(f'locs contains unknown channels {unknown}; '
                             f'available: {multi.channels}')
        out = {channel: locs.get(channel) for channel in multi.channels}
    elif isinstance(locs, (list, tuple)):
        if len(locs) != multi.n_channels:
            raise ValueError(f'locs has {len(locs)} entries, expected '
                             f'{multi.n_channels} in channel order '
                             f'{multi.channels}')
        out = dict(zip(multi.channels, locs))
    else:
        raise TypeError('locs must be a {channel: Localizations} mapping or a '
                        'list aligned with multi.channels')
    if not allow_missing:
        missing = [channel for channel, value in out.items() if value is None]
        if missing:
            raise ValueError(f'missing localizations for channels {missing}')
    return out


def emitter_centers_xy(emitters) -> np.ndarray:
    """Mean position of each emitter as ``(N, 2)`` ``[x, y]``.

    A ``(N, 2)`` array is accepted as-is (already ``[x, y]``); otherwise an
    ``Emitter`` sequence is averaged through ``emitter.centroid``.
    """
    if isinstance(emitters, np.ndarray) and emitters.ndim == 2:
        if emitters.shape[1] != 2:
            raise ValueError('centers array must have shape (N, 2)')
        return emitters.astype(float, copy=True)
    items = list(emitters)
    if not items:
        return np.zeros((0, 2), dtype=float)
    return np.asarray([emitter.centroid for emitter in items], dtype=float)


def _per_channel_kw(common, overrides, channels, name='link_kw_per_channel'):
    """把公共参数与逐通道覆盖合成 ``{channel: kwargs}``。"""
    common = dict(common or {})
    if overrides is None:
        return {channel: dict(common) for channel in channels}
    if not isinstance(overrides, Mapping):
        raise TypeError(f'{name} must be a {{channel: dict}} mapping')
    unknown = [channel for channel in overrides if channel not in channels]
    if unknown:
        raise ValueError(f'{name} has unknown channels {unknown}; '
                         f'available: {list(channels)}')
    for channel, extra in overrides.items():
        if not isinstance(extra, Mapping):
            raise TypeError(f'{name}[{channel!r}] must be a dict, '
                            f'got {type(extra).__name__}')
    return {channel: {**common, **overrides.get(channel, {})}
            for channel in channels}


def analyze_multichannel_emitters(multi, locs, *, channels=None,
                                  n_frames=None, link_kw=None,
                                  link_kw_per_channel=None,
                                  return_all=False, verbose=False):
    """Cluster emitters independently in every selected channel.

    Parameters
    ----------
    multi : MultiImageStack
        Drift-corrected stack; supplies channel names and the frame count.
    locs : dict or sequence
        ``{channel: Localizations}`` in each channel's own local coordinates.
    channels : sequence of str, optional
        Subset to cluster. Defaults to every channel.
    link_kw : dict, optional
        Common parameters forwarded to
        :func:`smlm.emitter.analyze_emitters_sm`.
    link_kw_per_channel : dict, optional
        ``{channel: {...}}`` merged **on top of** ``link_kw`` for that channel.
        The two channels usually want opposite settings:

        * the channel that only supplies **background-mask positions** should
          keep as many emitters as possible (``min_frames=1``, ``ecc_min=None``,
          ``do_merge=False``) -- dropping one leaves its core unpunched and
          inflates every nearby trace's background;
        * the channel that supplies **aperture centres / traces** is the one to
          filter hard (long ``min_frames``, tight ``r_merge``, ``ecc_min`` on).
    return_all : bool
        Also return the pre-filter emitters (useful for diagnosing how many
        each stage drops).
    verbose : bool
        Print per-channel counts.

    Returns
    -------
    dict or (dict, dict)
        ``{channel: [Emitter]}`` of the **kept** emitters, in channel order;
        with ``return_all`` a second mapping of the same shape holding every
        emitter before filtering.
    """
    if not isinstance(multi, MultiImageStack):
        raise TypeError('multi must be a MultiImageStack')
    names = channel_names(multi.channels, channels)
    by_channel = channel_locs(multi, locs)
    n_frames = multi.n_frames if n_frames is None else int(n_frames)
    kwargs = _per_channel_kw(link_kw, link_kw_per_channel, names)

    emitters, everything = {}, {}
    for channel in names:
        kept, every = analyze_emitters_sm(by_channel[channel], n_frames,
                                          **kwargs[channel])
        emitters[channel] = list(kept)
        everything[channel] = list(every)
        if verbose:
            print(f'  [{channel}] kept {len(kept)} emitters from '
                  f'{len(every)} linked ({len(by_channel[channel])} '
                  f'localizations)')
    return (emitters, everything) if return_all else emitters

