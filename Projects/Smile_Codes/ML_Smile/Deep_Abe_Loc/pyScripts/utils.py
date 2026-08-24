import matplotlib.pyplot as plt
import numpy as np

def show_two_image(X_train, y_train, num_images=3):
    
    fig, axes = plt.subplots(num_images, 2, figsize=(5, 2.5 * num_images))

    for i in range(num_images):
        # 左边显示 X_train
        axes[i, 0].imshow(X_train[i].squeeze(), cmap='hot', aspect='auto')
        axes[i, 0].set_title(f"Train Set X - Image {i}")
        axes[i, 0].axis("off")

        # 右边显示 y_train
        axes[i, 1].imshow(y_train[i].squeeze(), cmap='hot', aspect='auto')
        axes[i, 1].set_title(f"Train Set y - Image {i}")
        axes[i, 1].axis("off")

    plt.tight_layout()
    plt.show()

def show_three_image(X_train, y_train, y_pred, num_images=3):
    
    fig, axes = plt.subplots(num_images, 3, figsize=(7.5, 2.5 * num_images))

    for i in range(num_images):
        # 左边显示 X_train
        axes[i, 0].imshow(X_train[i].squeeze(), cmap='hot', aspect='auto')
        axes[i, 0].set_title(f"Train Set X - Image {i}")
        axes[i, 0].axis("off")

        # 右边显示 y_train
        axes[i, 1].imshow(y_train[i].squeeze(), cmap='hot', aspect='auto')
        axes[i, 1].set_title(f"Train Set y - Image {i}")
        axes[i, 1].axis("off")

        # 右边显示 y_pred
        axes[i, 2].imshow(y_pred[i].squeeze(), cmap='hot', aspect='auto')
        axes[i, 2].set_title(f"Predicted y - Image {i}")
        axes[i, 2].axis("off")

    plt.tight_layout()
    plt.show()


def render_3d_position(points):
    # 创建一个新的图形窗口
    fig = plt.figure(figsize=(3,3))
    # 添加一个三维坐标轴
    ax = fig.add_subplot(111, projection='3d')
    # 绘制散点图
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b', marker='o')
    # 设置坐标轴标签
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # 显示图形
    plt.show()
    # 提取 x 和 y 坐标
    x = points[:, 0]
    y = points[:, 1]

    # 创建二维散点图
    plt.figure(figsize=(2.5,2.5))
    plt.scatter(x, y, c='b', marker='o')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('XY')
    plt.grid(True)
    plt.show()
