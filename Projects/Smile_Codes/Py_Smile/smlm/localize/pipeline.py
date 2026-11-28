from __future__ import annotations

import numpy as np

from ..io.cache import result_path, run_cached, save_cached_locs
from .base import ClassicalPipeline
from .detect import LocalMaxDetector
from .filter import filter_locs
from .loc import MLELocalizer


def psf_sigma(locs, *, stat='median', clip=(0.5, 5.0)):
    """PSF σ (px) estimated from the **fitted** widths of one channel.

    Why measure it instead of using ``0.2165·λ/NA``
    -----------------------------------------------
    That formula assumes an ideal PSF in perfect focus. Real data has defocus
    and aberration and the fit comes out 20-50% wider (measured here: theory
    0.93 px vs fitted 1.35/1.43 px). σ is not a cosmetic parameter:

    * the aperture radius is ``ceil(k_sigma·σ)`` clamped to ``[2, 7]`` px, so
      0.93 -> 2 px and 1.43 -> 3 px - more than double the collecting area;
    * ``intensity`` is a matched-filter amplitude and a too-narrow kernel
      attenuates it by ``σ_k²/(σ_k² + σ_psf²)`` - about 0.6 for this data. Every
      SNR floor and every absolute threshold (``red_min``, ``min_drop_d``) then
      compares against a signal that was silently scaled down.

    Why the median rather than per-molecule
    ---------------------------------------
    Individual fits are noisy and, in a dense field, contaminated by neighbours.
    The median over tens of thousands of fits is stable. Per channel, not per
    file: chromatic defocus differs between the two views.

    Parameters
    ----------
    locs : Localizations
        Must carry ``.sigma`` (the fitted width in px).
    stat : 'median' | 'mean'
        Default median (outlier resistant).
    clip : (lo, hi), optional
        Drop fits outside this range first; ``None`` disables the clip.

    Returns
    -------
    float
        σ in px, ready for ``TraceFillPipeline(..., sigma=...)``.
    """
    if stat not in ('median', 'mean'):
        raise ValueError(f"stat must be 'median' or 'mean', got {stat!r}")
    raw = getattr(locs, 'sigma', None)
    s = np.asarray(raw, float) if raw is not None else np.zeros(0)
    s = s[np.isfinite(s)]
    if clip is not None:
        s = s[(s >= float(clip[0])) & (s <= float(clip[1]))]
    if s.size == 0:
        raise ValueError('no finite fitted sigma to estimate the PSF width '
                         f'from (clip={clip})')
    return float(s.mean() if stat == 'mean' else np.median(s))


def localize_stack(stack, *, k_sigma=1.5, det_kw=None, loc_kw=None,
                   filt_kw=None, src_path=None, read_cache=False,
                   verbose=False, tag=''):
    """Localize and filter one ``ImageStack`` with optional cache reuse.

    The same function handles a standalone stack and a channel returned by
    ``MultiImageStack.channel()``. Cache reads validate the image signature;
    detector, localizer, and filter settings are retained as provenance.
    """
    if not isinstance(read_cache, bool):
        raise TypeError('read_cache must be a bool')
    if read_cache:
        locs, _ = run_cached(
            stack, src_path, read_cache=True)
        if verbose:
            print(f'[cache] loaded {len(locs)} filtered locs -> '
                  f'{result_path(src_path)}')
            print('  cached localization parameters are in cache metadata')
        return locs

    if det_kw is None:
        det_kw = {}
    else:
        det_kw = dict(det_kw)
    det_kw.setdefault('verbose', verbose)
    loc_kw = {} if loc_kw is None else dict(loc_kw)
    filt_kw = {} if filt_kw is None else dict(filt_kw)

    method = ClassicalPipeline(
        detector=LocalMaxDetector(k_sigma=float(k_sigma), **det_kw),
        localizer=MLELocalizer(**loc_kw),
    )
    locs = method.run(stack)
    if verbose:
        print(f'  [{tag}] {len(locs)} raw locs')
    locs = filter_locs(locs, verbose=verbose, **filt_kw)

    if src_path is not None:
        save_cached_locs(
            locs,
            stack,
            src_path,
            method={
                'detector': method.detector,
                'localizer': method.localizer,
            },
            filter_kw=filt_kw,
            compressed=False,
        )
    if verbose:
        print(f'  [{tag}] {len(locs)} locs after filter '
              f'(k_sigma={k_sigma})')
    return locs
