clear;clc
tic
%load data
path='E:\data\ZWR\20230718\29-15Er\';
name='29-15Er_980nm_0p1S_5frames_500mA_fov5.sif';%%%%%%%%%
seq=readsif([path name]);
% filename='E:\data\dye\data_dye.xlsx';

ima1=double(seq.imageData(:,:,1));
figure;imagesc(ima1);axis square
% 数据转置
ima1=ima1(100:355,100:355);
ima1=(ima1.*5.75)./300;
ima1_1=ima1';
% 二值化
S=rescale(ima1);
bg=imopen(S,strel('disk',6));
S=S-bg;
% find the threshold
level=graythresh(S);
bw=imbinarize(S,level);
[label,object]=bwlabel(bw,8);
% bw = bwmorph(bw,"clean");
% RGB_label=label2rgb(label);imshow(RGB_label);
% 将regionprops函数得到的数据存入到stats
stats=regionprops('table',label,ima1,'Centroid',...
    'Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');
num_d=find(stats.Area>25);% dye

% 找寻两个光斑连到一起的区域，并且分开
if (isempty(num_d)==0)
    for i=1:size(num_d,1)
        state_d=stats(num_d(i),:);
        x_d=state_d.BoundingBox;
        state_a=p_d(x_d,ima1_1);
        stats=[stats;state_a];
    end
end
% 删去光斑连在一起的数组和面积较小的区域
stats(num_d,:)=[];
stats(stats.MaxIntensity>2.5*mean(stats.MaxIntensity),:)=[];
stats(stats.Area<3,:)=[];
sizeg=size(stats,1);
fr_D=zeros(sizeg,7);
Integral=zeros(sizeg,1);
FWHM=zeros(sizeg,2);
% 图像的背景均值
bg_mean=mean(ima1(:));

% 每个PSF的半径和位置
centers = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2+1;

for i=1:sizeg
    % 确定需要拟合的区域
    x_edge=stats.BoundingBox(i,:);
    row_edge=round(x_edge(1)-1:x_edge(1)+x_edge(3)+1);
    col_edge=round(x_edge(2)-1:x_edge(2)+x_edge(4)+1);
    % 确定范围
    col_edge(col_edge>256)=[];
    col_edge(col_edge<0)=[];
    row_edge(row_edge>256)=[];
    row_edge(row_edge<0)=[];
    % 根据以上范围选取数据
    imgGrayData = ima1_1(row_edge,col_edge);
    % 创建网格
    [X,Y]=meshgrid(col_edge,row_edge);
    % 拟合前准备
    [xData, yData, zData] = prepareSurfaceData( X,Y, imgGrayData);
    % 设置fittype
    ft=fittype(@(u1,u2,ss1,ss2,A,th,C,x,y)A*exp(-(cos(th)*x+sin(th)*y-u1).^2./(2*ss1^2)- ...
        (-sin(th)*x+cos(th)*y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2)+C,'independent', ...
        {'x','y'},'dependent','z');
    % 设置startpoint
    max1=max(imgGrayData(:));% 峰值
    % sigmax=1.2;sigmay=1.2;
    startpoint_G=[centers(i,2),centers(i,1),1.2,1.2,...
        max1*2*pi*1.2*1.2,0,bg_mean];% 一一对应fittype中的参数
    % 用fit函数进行拟合
    fitResult=fit([xData,yData],zData,ft,'startpoint',startpoint_G);
    % 拟合得到的数据存入fr_D
        % figure
        % surf(X,Y,imgGrayData)
        % plot(fitResult,[xData,yData],zData)
    fr_D(i,:)=coeffvalues(fitResult);
    if (fr_D(i,5)>2* max1 * 2 * pi * 1.2 * 1.2)
        Integral(i)=-1;
        continue
    end
    % 对每个PSF进行积分
    Integral(i)=integration_DF(fr_D(i,:));
    % [fr_D(i,1),fr_D(i,2)]=gfcenter(fr_D(i,1),fr_D(i,2),fr_D(i,6));
end

% 绘图
figure
imagesc(ima1);colormap('gray');axis square;colorbar;
% 圈出每个数据点
hold on
viscircles(centers,radii,'color','#0072BD','LineWidth',1.5,'EnhanceVisibility',false);
plot(centers(18,1),centers(18,2),'+','Color','red')
% plot(fr_D(:,1),fr_D(:,2),'.','Color','red')
hold off

% 将数据导入到一个Excel表中

% 计算每个PSF的半峰宽FWHM_x & FWHM_Y
n=fr_D;
[ss1,ss2]=deal(n(:,3),n(:,4));
FWHM(:,1)=2.*ss1.*sqrt(-2*log(1/2));% X轴的半峰宽
FWHM(:,2)=2.*ss2.*sqrt(-2*log(1/2));% Y轴的半峰宽
stats.FWHM_X=FWHM(:,1);stats.FWHM_y=FWHM(:,2);
stats.Integral=Integral./seq.exposureTime;
% 计算偏转角度Theta
stats.Theta=fr_D(:,6).*360./pi;

% writematrix(Integral,filename,"WriteMode","append");
toc
figure
histogram(Integral)

%%
% seq=readsif('20230320_Hcy1_100nM_2mW_0pO2s_2000frame_fov2.sif');
seq=readsif('E:\data\20230401\Cy3-Cy5\220230401_cy3_cy5_532nm_0p1S_500frames_2mW_fov2.sif');
% 设置视频参数
outputVideo = VideoWriter('output_video.avi');
outputVideo.FrameRate = 1/seq.exposureTime;  % 每秒帧数
imgdata=seq.imageData;
open(outputVideo);
% 创建视频帧
numFrames = size(seq.imageData,3);  % 视频帧数
for i = 1:numFrames
    % 完全去除白边
    % imgdata(:,:,i)=rescale(imgdata(:,:,i));
    % imshow(imgdata(:,:,i),'border','tight','initialmagnification','fit');
    % 白边未完全去除
    imagesc(imgdata(:,:,i));
    set(gcf,'Position',[300,200,420,420]);
    axis square;
    colormap("gray");
    set(gca,'looseInset',[0 0 0 0])
    set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
    frame=getframe(gcf);
    % 写入当前帧到视频
    writeVideo(outputVideo, frame);  % 注意需要转换为uint8类型
    % 显示进度
    fprintf('Frame %d/%d\n', i, numFrames);
end

% 关闭视频
close(outputVideo);
disp('视频创建完成。');








