function coeFit=GF_DF(x,ima1)

clc;
ima1_1=ima1';
row_edge=round(x(1)-0.5*x(3):x(1)+1.5*x(3));
col_edge=round(x(2)-0.5*x(4):x(2)+1.5*x(4));
imgGrayData = ima1_1(row_edge,col_edge);
% 图像方向和数据方向一致
% 图像大小
% [row, col] = size(imgGrayData);
% 生成网格坐标
[X0, Y0] = meshgrid(col_edge, row_edge);
% 序号太大，拟合效果不好，建议缩小
x0 = X0 / 10;
y0 = Y0 / 10;
% uint8 转换为 double
imgData = double(imgGrayData);
% 画原图
% figure
% subplot(1, 2, 1)
% surf(X0, Y0, imgGrayData);
% imagesc(imgData)
% % smooth
% shading interp;
% 确定数据阈值
mean_img=mean(mean(imgGrayData));
thresh = 1.2*mean_img;
% 有效数据下标
[rows, cols] = find(imgData >= thresh);
% 有效数据个数
len = length(rows);
% 初始化
xfit = zeros(len, 1);
yfit = zeros(len, 1);
zfit = zeros(len, 1);
for i = 1 : len
    xfit(i) = x0(rows(i), cols(i));
    yfit(i) = y0(rows(i), cols(i));
    zfit(i) = imgData(rows(i), cols(i));
end
% 取对数
zfit = log(zfit);
% 系数矩阵
xfit2 = xfit .* xfit;
yfit2 = yfit .* yfit;
A = [xfit2, yfit2, xfit, yfit, ones(len, 1)];
val = A \ zfit;
% 将拟合结果还原
sigmaX = sqrt(-0.5 / val(1));
sigmaY = sqrt(-0.5 / val(2));
centerX = sigmaX * sigmaX * val(3);
centerY = sigmaY * sigmaY * val(4);
Amp = exp(val(5) + 0.5 * centerX * val(3) + 0.5 * centerY * val(4));
%q 确定拟合中心和σx&y
centerX = centerX * 10;
centerY = centerY * 10;
sigmaX = sigmaX * 10;
sigmaY = sigmaY * 10;
coeFit = [Amp; centerX; sigmaX; centerY; sigmaY];
% 画拟合图
% subplot(1, 2, 2)
% imgDataFit = coeFit(1) * exp(-(X0 - coeFit(2)).^2 / (2 * coeFit(3) * coeFit(3)) - (Y0 - coeFit(4)).^2 / (2 * coeFit(5) * coeFit(5)));
% surf(X0, Y0, imgDataFit);
% imagesc(imgDataFit)
% title('Gauss Fitting');
% shading interp