import torch
import torch.nn as nn
import torch.nn.functional as F
from d2l import torch as d2l
from piq import ssim

class L1L2Loss(nn.Module):
    def __init__(self, alpha=0.6):  # α 控制 L1 和 L2 的比例
        super(L1L2Loss, self).__init__()
        self.alpha = alpha
        self.l1_loss = torch.nn.SmoothL1Loss()
        self.l2_loss = nn.MSELoss()
        self.ssim = ssim

    def forward(self, y_pred, y_true):
        l1 = self.l1_loss(y_pred, y_true) * 100 * 100 # for convinence of visualization
        l2 = self.l2_loss(y_pred, y_true) * 100 * 100
        # ssim = 1 - self.ssim(y_pred, y_true, data_range=1.0) 
        return  self.alpha * l1 + (1 - self.alpha) * l2

def  evaluate_accuracy_gpu(net, data_iter, device=None):
    if isinstance(net, nn.Module):
        net.eval()
        if not device:
            device = next(iter(net.parameters())).device
    metric = d2l.Accumulator(2)
    with torch.no_grad():
        for X, y in data_iter:
            if isinstance(X, list):
                X = [x.to(device) for x in X]
            else:
                X = X.to(device)
            y = y.to(device)
            ssim_value = ssim(net(X), y, data_range=1.0)  # 计算 SSIM
            metric.add(ssim_value.item() * X.shape[0] / 2, X.shape[0])
    return metric[0]/metric[1]

def train_cnn_denoise_deabe(net, train_iter, test_iter, num_epochs, lr, device, loss_function):
    def init_weights(m):
        if type(m) == nn.Conv2d or type(m) == nn.Linear:
            nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
    # print('Initializing model on CPU...')
    # net = net.to('cpu')
    net.apply(init_weights)  # **先初始化权重**
    net = net.to(device)  # **再转到 MPS/GPU**
    print('training on', device)
    optimizer = torch.optim.Adam(net.parameters(), lr=lr)
    loss = loss_function()
    animator = d2l.Animator(xlabel='epoch', xlim=[1, num_epochs], ylim=[0, 0.5], legend=['train loss', 'train acc', 'test acc'])
    timer, num_batches = d2l.Timer(), len(train_iter)
    for epoch in range(num_epochs):
        metric = d2l.Accumulator(3)
        net.train()
        for i, (X, y) in enumerate(train_iter):
            timer.start()
            X, y = X.to(device), y.to(device)
            optimizer.zero_grad()
            y_hat = net(X) # predicted y
            l = loss(y_hat, y)
            l.backward()
            optimizer.step()
            with torch.no_grad():
                metric.add(l * X.shape[0], ssim(y_hat, y) * X.shape[0]/2 , X.shape[0])
            timer.stop()
            train_l = metric[0] / metric[2]
            train_ssim = metric[1] / metric[2]
            if (i + 1) % (num_batches // 3) == 0 or i == num_batches - 1:
                animator.add(epoch + (i + 1) / num_batches, (train_l, train_ssim, None))
        test_ssim = evaluate_accuracy_gpu(net, test_iter)
        animator.add(epoch + 1, (None, None, test_ssim))
    print(f'loss {train_l:.3f}, train ssim {train_ssim:.3f}, test ssim {test_ssim:.3f}')
    print(f'{metric[2] * num_epochs / timer.sum():.1f} examples/sec on {str(device)}')