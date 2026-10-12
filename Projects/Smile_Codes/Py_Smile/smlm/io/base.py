# smlm/io/base.py
from abc import ABC, abstractmethod
from ..datatypes import RawDataset
import numpy as np


class Loader(ABC):
    @abstractmethod
    def load(self, path, interval=1) -> RawDataset:
        ...