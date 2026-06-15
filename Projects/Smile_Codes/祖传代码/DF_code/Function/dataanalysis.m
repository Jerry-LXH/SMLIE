function [y,error_data]=dataanalysis(data_UCNP)
% data_UCNP=sortrows(data_UCNP,1,'descend');
data_UCNP(data_UCNP(:)<0)=[];
% diff_data=abs(diff(data_UCNP));
Q1 = prctile(data_UCNP, 25);  % 第一四分位数（25% 分位数）
Q3 = prctile(data_UCNP, 75);  % 第三四分位数（75% 分位数）
IQR = Q3 - Q1;          % 四分位距（IQR）
% 确定范围的大小scale
scale = 3;
upper_limit = Q3 + scale * IQR;  % 上限
lower_limit = Q1 - scale * IQR;  % 下限
data_UCNP = data_UCNP(data_UCNP <= upper_limit & data_UCNP >= lower_limit);
% 当方差大于均值时，数据不好看，因此，进一步缩小数据的范围，减小方差
if (std(data_UCNP) > 0.9 * mean(data_UCNP))
    scale=scale/1.5;
    upper_limit = Q3 + scale * IQR;  % 上限
    lower_limit = Q1 - scale * IQR;  % 下限
    data_UCNP = data_UCNP(data_UCNP <= upper_limit & data_UCNP >= lower_limit);
end

figure
h = histogram(data_UCNP);
% 每个histogram柱的y值
ycount = histcounts(data_UCNP,h.BinEdges);
% x轴每个柱的宽度
width=h.BinWidth;
% x轴范围
xlimit=h.BinLimits;
% 确定x值分布
xrange=xlimit(1)+width/2:width:xlimit(2)-width/2;
gauss_fittype = fittype('a*exp(-(x-b)^2/(2*c^2))+d', ...
    'independent', 'x', ...
    'dependent', 'y', ...
    'coefficients', {'a', 'b', 'c','d'});
% 设置startpoint
max1=max(ycount(:))-min(ycount(:));% 峰值
u1mean=mean(data_UCNP(:));% 均值
sigma=std(data_UCNP(:));% sigma
startpoint_G=[max1,u1mean,sigma,min(ycount(:))];% 一一对应fittype中的参数
% 用fit函数进行拟合
fitResult=fit(xrange',ycount',gauss_fittype,'startpoint',startpoint_G);
hold on
plot(fitResult)
hold off
error_data=std(data_UCNP);
xdata=coeffvalues(fitResult);
y=xdata(2);