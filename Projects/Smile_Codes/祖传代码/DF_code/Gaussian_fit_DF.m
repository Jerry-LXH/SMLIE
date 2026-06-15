function [centers,background,rsquare]= Gaussian_fit_DF(loc,ima1,edge)
ima1=ima1';
imgGrayData = ima1(max(round(loc(2))-edge,1):min(round(loc(2))+edge,256),max(round(loc(1))-edge,1):min(round(loc(1))+edge,256));
% 创建网格
[X,Y]=meshgrid(max(round(loc(2))-edge,1):min(round(loc(2))+edge,256),max(round(loc(1))-edge,1):min(round(loc(1))+edge,256));
% 拟合前准备
[xData, yData, zData] = prepareSurfaceData( X,Y, imgGrayData);
ft=fittype(@(u1,u2,ss1,ss2,A,th,C,x,y)A*exp(-(cos(th)*x+sin(th)*y-u1).^2./(2*ss1^2)- ...
    (-sin(th)*x+cos(th)*y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2)+C,'independent', ...
    {'x','y'},'dependent','z');
startpoint_G=[loc(2),loc(1),1.2,1.2,max(imgGrayData(:))*2*pi*1.2*1.2,0,mean(ima1(:))];
try
    [fitResult,gof]=fit([xData,yData],zData,ft,'startpoint',startpoint_G);
    rsquare=cell2mat(struct2cell(gof));
    fr_D=coeffvalues(fitResult);
    % 然后将拟合数据的每个参数分开
    numc=num2cell(fr_D,1);
    [u1,u2,ss1,ss2,A,th,background]=deal(numc{:});
    % 转换矩阵，将拟合的u1和u2转变为中心点
    cc=[cos(th),sin(th);-sin(th),cos(th)];
    xy=[u1;u2];
    xxyy=cc^(-1)*xy;
    % 计算得到的中心点分别存入center_x/y中
    centers=xxyy;
    center_x=xxyy(1);
    center_y=xxyy(2);
    [u1,u2]=deal(xxyy(1),xxyy(2));
    % 先确定积分函数的形式
    f=@(x,y)A*exp(-(x-u1).^2./(2*ss1^2)-(y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2);
    % 选取-10/+10σ的范围进行积分
    Integral=integral2(f,u1-5*ss1,u1+5*ss1,u2-5*ss2,u2+5*ss2);
catch Error
    Integral=-1;
    disp(Error.message);
end

% imagesc(imgGrayData)
% pause(1)
% figure
% f=surf(X,Y,imgGrayData);
% shading interp
% plot(fitResult,[xData,yData],zData)
