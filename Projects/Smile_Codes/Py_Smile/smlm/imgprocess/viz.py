# smlm/imgprocess/viz.py
"""单纯视觉处理: SOFI / RIED 结果的显示与对比辅助。"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt


def plot_sofi(img, ax=None, vmax_percentile=99.8, cmap="hot", title=None,
              show_colorbar=False):
    """显示一张 SOFI / RIED 图像 (自动裁到 vmax_percentile 分位)。

    img: (H, W) 或 (n, H, W) ndarray; 多张时取第一张。
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))
    img = np.asarray(img)
    if img.ndim == 3:
        img = img[0]
    v = np.percentile(img, vmax_percentile)
    im = ax.imshow(np.clip(img, 0, v), cmap=cmap,
                   interpolation="nearest", aspect="equal")
    if title:
        ax.set_title(title)
    ax.axis("off")
    if show_colorbar:
        plt.colorbar(im, ax=ax, fraction=0.046)
    return ax


def compare_widefield_sofi(stack, img, ax=None, vmax_percentile=99.8,
                           cmap_wide="gray", cmap_sofi="hot", titles=None):
    """并排: 宽场均值 vs SOFI 图像。stack: (T,H,W) 或 ImageStack。"""
    if ax is None:
        _, ax = plt.subplots(1, 2, figsize=(10, 5))
    data = stack.data if hasattr(stack, "data") else np.asarray(stack)
    ax[0].imshow(data.mean(0), cmap=cmap_wide,
                 interpolation="nearest", aspect="equal")
    ax[0].set_title((titles or ("Widefield (mean)", "SOFI"))[0])
    ax[0].axis("off")
    plot_sofi(img, ax=ax[1], vmax_percentile=vmax_percentile, cmap=cmap_sofi,
              title=(titles or ("Widefield (mean)", "SOFI"))[1])
    return ax


def line_profile(img, row=None, ax=None, normalize=True, label=None,
                 pixel_size=None):
    """画某一行(或列)的强度剖面, 用于对比 PSF 宽度。返回 ax。

    pixel_size: 有效像素尺寸 (nm/px)。给定时 x 轴以"距中心"的 nm 表示,
    这样不同分辨率(不同网格大小)的图像剖面才能在物理坐标上对齐对比
    (否则分辨率越高曲线被拉得越长, 无法直接比较)。
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    img = np.asarray(img)
    if img.ndim == 3:
        img = img[0]
    if row is None:
        row = img.shape[0] // 2
    p = img[row]
    if normalize and np.ptp(p) > 0:
        p = (p - p.min()) / (p.max() - p.min())
    if pixel_size:
        x = (np.arange(len(p)) - len(p) / 2) * pixel_size
        ax.plot(x, p, label=label)
        ax.set_xlabel("x (nm from center)")
    else:
        ax.plot(p, label=label)
        ax.set_xlabel("x (pixel)")
    ax.set_ylabel("intensity (normalized)")
    if label:
        ax.legend()
    return ax
