"""Reference-based ROI extraction for multichannel image stacks."""
from __future__ import annotations

from collections.abc import Iterable, Sequence
from operator import index

import numpy as np

from ..datatypes import ImageStack, MultiImageStack


def _integer_pair(value: Iterable[int], name: str) -> tuple[int, int]:
    try:
        pair = tuple(value)
    except TypeError:
        raise ValueError(f'{name} must contain exactly two integers') from None
    if len(pair) != 2:
        raise ValueError(f'{name} must contain exactly two integers')
    try:
        return tuple(index(v) for v in pair)
    except TypeError:
        raise ValueError(f'{name} must contain exactly two integers') from None


class MultichannelROI:
    """Extract equally sized channel windows from an ``ImageStack``.

    ``y`` and ``x`` define channel 0 using half-open detector-coordinate ranges.
    Each ``offsets`` entry is a ``(dy, dx)`` displacement from channel 0; the
    first entry must be ``(0, 0)``. Channel windows retain the same shape.
    """

    def __init__(self, y: Iterable[int], x: Iterable[int],
                 offsets: Iterable[Iterable[int]] = ((0, 0),),
                 channels: Sequence[str] | None = None) -> None:
        self.y = _integer_pair(y, 'y')
        self.x = _integer_pair(x, 'x')
        if self.y[0] < 0 or self.x[0] < 0:
            raise ValueError('channel 0 ROI coordinates must be non-negative')
        if self.y[1] <= self.y[0] or self.x[1] <= self.x[0]:
            raise ValueError(f'channel 0 ROI must have positive size: '
                             f'y={self.y}, x={self.x}')

        try:
            self.offsets = tuple(_integer_pair(offset, f'offsets[{i}]')
                                 for i, offset in enumerate(offsets))
        except TypeError:
            raise ValueError('offsets must contain one or more (dy, dx) pairs') \
                from None
        if not self.offsets:
            raise ValueError('offsets must contain at least channel 0')
        if self.offsets[0] != (0, 0):
            raise ValueError('offsets[0] for channel 0 must be (0, 0)')

        if channels is None:
            self.channels = tuple(f'ch{i}' for i in range(len(self.offsets)))
        else:
            self.channels = tuple(channels)
            if len(self.channels) != len(self.offsets):
                raise ValueError(f'{len(self.offsets)} offsets require the same '
                                 f'number of channel names, got {len(self.channels)}')
            if any(not isinstance(channel, str) for channel in self.channels):
                raise TypeError('channel names must be strings')
            if len(set(self.channels)) != len(self.channels):
                raise ValueError('channel names must be unique')

    @property
    def shape(self) -> tuple[int, int]:
        return self.y[1] - self.y[0], self.x[1] - self.x[0]

    @property
    def windows(self) -> dict[str, tuple[int, int, int, int]]:
        """Return channel windows as ``{name: (y0, x0, height, width)}``."""
        height, width = self.shape
        return {
            channel: (self.y[0] + dy, self.x[0] + dx, height, width)
            for channel, (dy, dx) in zip(self.channels, self.offsets)
        }

    def _checked_windows(
            self, image: ImageStack
    ) -> tuple[np.ndarray, dict[str, tuple[int, int, int, int]]]:
        if not isinstance(image, ImageStack):
            raise TypeError('MultichannelROI expects an ImageStack')
        data = np.asarray(image.data)
        if data.ndim != 3:
            raise ValueError(f'ImageStack.data must be (T,H,W), got {data.shape}')

        _, image_height, image_width = data.shape
        windows = self.windows
        for channel, (y0, x0, height, width) in windows.items():
            if (y0 < 0 or x0 < 0 or y0 + height > image_height
                    or x0 + width > image_width):
                raise ValueError(
                    f'{channel} ROI y[{y0},{y0 + height}) '
                    f'x[{x0},{x0 + width}) is outside image '
                    f'{image_height}x{image_width}')
        return data, windows

    def __call__(self, image: ImageStack) -> MultiImageStack:
        """Convert an ``ImageStack`` into a ``MultiImageStack``."""
        data, windows = self._checked_windows(image)
        cropped = np.stack([
            data[:, y0:y0 + height, x0:x0 + width]
            for y0, x0, height, width in windows.values()
        ], axis=0)
        roi_metadata = {
            channel: {'y': (y0, y0 + height), 'x': (x0, x0 + width)}
            for channel, (y0, x0, height, width) in windows.items()
        }
        metadata = dict(image.metadata)
        metadata['multichannel_roi'] = {
            'reference': self.channels[0],
            'windows': roi_metadata,
        }
        return MultiImageStack(
            data=cropped,
            channels=list(self.channels),
            reference=0,
            offsets=np.array([window[:2] for window in windows.values()], dtype=int),
            pixel_size=image.pixel_size,
            metadata=metadata,
        )

    def show_roi(self, image: ImageStack,
                 figsize: tuple[float, float] = (10, 6)):
        """Return a figure showing every channel ROI on the source image."""
        _, windows = self._checked_windows(image)
        from .viz import plot_multichannel_roi
        return plot_multichannel_roi(image, windows, figsize=figsize)

    def __repr__(self) -> str:
        return (f'MultichannelROI(y={self.y}, x={self.x}, '
                f'offsets={self.offsets}, channels={self.channels})')
