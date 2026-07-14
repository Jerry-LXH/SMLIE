"""Time-series engine: advance photophysics across frames.

SKELETON -- signatures and docstrings only.

Separated from the pipeline so the per-frame stepping logic (state reset, time
bookkeeping, exposure integration) lives in one place and is testable on its
own.
"""
from __future__ import annotations

import numpy as np


def run_movie(pipeline, n_frames: int) -> np.ndarray:
    """Drive a pipeline across ``n_frames``.

    Steps:
      1. Reset stateful stages once (photophysics initial state, noise maps).
      2. For each frame: set ctx.time / frame_index / exposure, run all stages,
         collect ctx.image.
      3. Stack frames -> (n_frames, H, W).

    Stateful stages persist their internal state (e.g. emitter ON/OFF) between
    frames; that persistence is what produces realistic blinking traces.
    """
    raise NotImplementedError
