"""Persistence for multichannel calibration bundles."""
from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO
import json
import os
from typing import Sequence

import numpy as np

from .roi import MultichannelROI
from .transform import DualViewCalibration

BUNDLE_FORMAT = 'smlm.multichannel_calibration'
BUNDLE_VERSION = 1


@dataclass
class MultiChannelCalibrationBundle:
    """A portable bundle of source-file provenance, ROI, and channel calibrations."""

    channels: tuple[str, ...]
    reference: str
    roi: dict
    source_files: list[dict[str, str]]
    calibrations: list[DualViewCalibration | None]

    def __post_init__(self):
        self.channels = tuple(self.channels)
        self.roi = dict(self.roi)
        self.source_files = [dict(item) for item in self.source_files]
        self.calibrations = list(self.calibrations)

        if not self.channels or len(set(self.channels)) != len(self.channels):
            raise ValueError('channels must contain unique channel names')
        if self.reference not in self.channels:
            raise ValueError('reference must be present in channels')
        if len(self.calibrations) != len(self.channels):
            raise ValueError('calibrations length must match channels')
        reference_index = self.channels.index(self.reference)
        for index, calibration in enumerate(self.calibrations):
            if index == reference_index:
                if calibration is not None:
                    raise ValueError('reference channel calibration must be None')
                continue
            if not isinstance(calibration, DualViewCalibration):
                raise TypeError(
                    f'calibration for {self.channels[index]!r} must be a '
                    'DualViewCalibration')
            if (calibration.reference != self.reference
                    or calibration.moving != self.channels[index]):
                raise ValueError(
                    f'calibration direction for {self.channels[index]!r} '
                    'does not match the bundle channels')
        for item in self.source_files:
            if not {'path', 'name'} <= item.keys():
                raise ValueError('each source file must contain path and name')

    @classmethod
    def from_fit(cls, calibrations, roi: MultichannelROI,
                 source_files: Sequence[os.PathLike | str]):
        """Create a bundle from fit results and the inputs used to produce them."""
        if not isinstance(roi, MultichannelROI):
            raise TypeError('roi must be a MultichannelROI')
        if isinstance(source_files, (str, os.PathLike)):
            source_files = [source_files]
        paths = [os.fspath(path) for path in source_files]
        if not paths:
            raise ValueError('source_files is empty')

        calibrations = list(calibrations)
        if len(calibrations) != len(roi.channels):
            raise ValueError('calibrations length must match roi.channels')
        reference_slots = [
            index for index, calibration in enumerate(calibrations)
            if calibration is None
        ]
        if len(reference_slots) != 1:
            raise ValueError(
                f'Expected one reference slot, got {len(reference_slots)}')
        reference = roi.channels[reference_slots[0]]
        roi_metadata = {
            'type': 'MultichannelROI',
            'coordinate_anchor': roi.channels[0],
            'reference': reference,
            'channels': list(roi.channels),
            'y': list(roi.y),
            'x': list(roi.x),
            'offsets': {
                channel: list(offset)
                for channel, offset in zip(roi.channels, roi.offsets)
            },
            'windows': {
                channel: list(window)
                for channel, window in roi.windows.items()
            },
        }
        source_metadata = [
            {'path': path, 'name': os.path.basename(path)}
            for path in paths
        ]
        return cls(
            channels=roi.channels,
            reference=reference,
            roi=roi_metadata,
            source_files=source_metadata,
            calibrations=calibrations,
        )

    def save(self, path):
        """Save the complete bundle as one compressed, versioned NPZ file."""
        calibration_entries = []
        archive = {}
        for index, calibration in enumerate(self.calibrations):
            if calibration is None:
                calibration_entries.append(None)
                continue
            key = f'calibration_{index}'
            stream = BytesIO()
            calibration.save(stream)
            archive[key] = np.frombuffer(stream.getvalue(), dtype=np.uint8)
            calibration_entries.append(key)

        manifest = {
            'format': BUNDLE_FORMAT,
            'version': BUNDLE_VERSION,
            'channels': list(self.channels),
            'reference': self.reference,
            'roi': self.roi,
            'source_files': self.source_files,
            'calibrations': calibration_entries,
        }
        archive['manifest'] = np.asarray(
            json.dumps(manifest, ensure_ascii=False))
        with open(path, 'wb') as file:
            np.savez_compressed(file, **archive)
        return path

    @classmethod
    def load(cls, path):
        """Load a bundle and reconstruct its ``DualViewCalibration`` objects."""
        with np.load(path, allow_pickle=False) as archive:
            manifest = json.loads(str(archive['manifest'].item()))
            if manifest.get('format') != BUNDLE_FORMAT:
                raise ValueError('file is not a multichannel calibration bundle')
            if manifest.get('version') != BUNDLE_VERSION:
                raise ValueError(
                    f'unsupported bundle version: {manifest.get("version")!r}')
            calibrations = []
            for key in manifest['calibrations']:
                if key is None:
                    calibrations.append(None)
                    continue
                calibration_data = archive[key].tobytes()
                calibrations.append(
                    DualViewCalibration.load(BytesIO(calibration_data)))

        return cls(
            channels=tuple(manifest['channels']),
            reference=manifest['reference'],
            roi=manifest['roi'],
            source_files=manifest['source_files'],
            calibrations=calibrations,
        )
