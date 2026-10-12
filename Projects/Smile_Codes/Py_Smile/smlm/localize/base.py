# smlm/localize/base.py
from abc import ABC, abstractmethod
from ..datatypes import ImageStack, Candidates, Localizations


class Detector(ABC):
    @abstractmethod
    def detect(self, stack: ImageStack) -> Candidates: ...


class Localizer(ABC):
    @abstractmethod
    def localize(self, stack: ImageStack, cand: Candidates) -> Localizations: ...


class LocalizationMethod(ABC):
    """下游只认这一个接口: method.run(stack) -> Localizations。"""
    @abstractmethod
    def run(self, stack: ImageStack) -> Localizations: ...


class ClassicalPipeline(LocalizationMethod):
    """传统路线 = 检测器 + 定位器 组合。"""
    def __init__(self, detector: Detector, localizer: Localizer):
        self.detector = detector
        self.localizer = localizer

    def run(self, stack):
        cand = self.detector.detect(stack)
        return self.localizer.localize(stack, cand)