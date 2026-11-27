from __future__ import annotations

from ..io.cache import result_path, run_cached, save_cached_locs
from .base import ClassicalPipeline
from .detect import LocalMaxDetector
from .filter import filter_locs
from .loc import MLELocalizer


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
