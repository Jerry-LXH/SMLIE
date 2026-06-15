%%
close all
clear;
% 100为进行100次模拟
time=100;
% 进行100次模拟，每次模拟10S，每个数据代表0.1S
x=zeros(time,time);
% 设置初始条件，每次模拟一A状态开始
x(:,1)=1;
% pa为k(A->B,1)*t(0.1S)，同理k(B->A)=0.5
pa=0.1;pb=0.05;
for j=1:time
% 因为第一个状态已确认为A，所以从第二个开始算
for i=2:time
    % 先判断状态，这里假设1为状态A，0为状态B
    if (x(j,i-1)==1)
        % rand(1) 在0~1取一个随机数，计算能否从状态A到B，假设随机数小于pa，则转化为B
        if (rand(1)<pa)
            x(j,i)=0;
        % rand(1) 大于pa，则无法转化为B
        else
            x(j,i)=1;
        end
    % 同上，先判断状态为B，则进一步进行概率计算
    else
        if (rand(1)<pb)
            x(j,i)=1;
        else
            x(j,i)=0;
        end
    end
end
end
% 通过sum函数进行总和计算处于A状态的数量，其中sum函数对每一列进行加和，因此y可以统计每个时间处于A状态的数量
y=sum(x);
t=0:0.1:time/10-0.1;
% 根据老师课上的公式，从而可以得到宏观状态下的分子状态随时间的变化
yy=time-(time*2/3*(1-exp(-1.5*t)));
% 之后，进行绘图，画出上述得到的两个数据
figure;plot(t,y,t,yy)
%%
close all
clc
% 创建一个分子的随时间的变化，在A，B态间的转变
zz=zeros(1,10);
zz(1)=1;
pa=0.1;pb=0.05;
for i=2:10000
    if (zz(i-1)==1)
        if (rand(1)<pa)
            zz(i)=0;
        else
            zz(i)=1;
        end
    else
        if (rand(1)<pb)
            zz(i)=1;
        else
            zz(i)=0;
        end
    end
end
% plot(zz);
% ylim([-0.5 1.5]);

% 计算A分子的dewll时间，即在A状态的持续时间
num=find(zz(:)==0);
num=[0;num];
dwell=diff(num)-1;
dwell(dwell(:)==0)=[];
dwell=dwell./10;
histogram(dwell)
mean(dwell)
%%
ax=0;
a1=0;
ll=[];
for i=1:10000
    if zz(i)==1
        a2=i;
    else
        ax=a2-a1;
        a2=i;
        a1=a2;
        if ax>0
            ll(end+1)=ax;
        end
    end
end

% if zz(i)==1
% ll(end+1)=a2-a1;
% end
sum(zz)

%%
clear;
clc
tic
% 多进程处理
pnumber=parpool(3);
for q=1:1
    qq={'590';'460';'372';'308';'260';'230'};% 808nm
%     qq={'500';'300';'188';'128';'100';'85'};% 980nm
    time={'0p1S';'0p1S';'0p2S';'0p4S';'0p8S';'2S'};
for x=1:1
%load data
path='E:\data\sm_l0308\20230302_808\UCNP_Nd_ZTL_1wX_808\';
name=['UCNP_Nd_ZTL_1wX_20230308_808nm_' char(qq(q)) 'mA_' char(time(q)) '_3frame_fov' char(string(x)) '.sif'];
% data_excel 保存位置
filename=[path,'data_UCNP_' char(qq(q)) '.xlsx'];
seq=readsif([path name]);

for m=1:3
    % 显示现在正在处理的数据
    name_ima=[ char(qq(q)) 'mA' ' fov' char(string(x))  ' frame' char(string(m)) ];
    disp(['现在是' name_ima '...'])
    ima1=double(seq.imageData(:,:,m));
    ima1=ima1(100:355,100:355);
    ima1=(ima1.*5.75)./300;
    % 数据转置
    ima1_1=ima1';
    % 二值化
    S=rescale(ima1);
    bg=imopen(S,strel('disk',4));
    S=S-bg;
    % find the threshold
    level=graythresh(S);
    bw=imbinarize(S,level);
%     bw = bwmorph(bw,"clean");
%     % 当图像比较乱时，可以通过该方法起到一定效果
%     bw = bwmorph(bw,"open");
    [label,object]=bwlabel(bw,4);
    % 将regionprops函数得到的数据存入到stats
    stats=regionprops('table',label,ima1,'Centroid',...
        'Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');
    num_d=find(stats.Area>30);% UCNP
    % 找寻两个光斑连到一起的区域，并且分开
    if (isempty(num_d)==0)
        for i=1:size(num_d,1)
            state_d=stats(num_d(i),:);
            x_d=state_d.BoundingBox;
            state_a=p_d(x_d,ima1_1);
            stats=[stats;state_a];
        end
    end
    % 删去光斑连在一起的数组
    stats(num_d,:)=[];
%     stats(stats.MaxIntensity<0.8*mean(stats.MaxIntensity),:)=[];
    stats(stats.MaxIntensity>2.5*mean(stats.MaxIntensity),:)=[];
    stats(stats.Area<5,:)=[];

    sizeg=size(stats,1);

    fr_D1=zeros(sizeg,7);

    FWHM=zeros(sizeg,2);
    
    center_x=zeros(sizeg,1);
    center_y=zeros(sizeg,1);
    % 图像的背景均值
    bg_mean=mean(ima1_1(:));
    x_edge=stats.BoundingBox(:,1:4);
    x_edge(:,1:2)=x_edge(:,1:2)+0.5;
    % 区域的row参数
    row_edge1=[x_edge(:,1)-1 x_edge(:,1)+x_edge(:,3)+1];
    % 区域的col参数
    col_edge1=[x_edge(:,2)-1 x_edge(:,2)+x_edge(:,4)+1];
    % 去除超出数组范围的参数
    col_edge1(col_edge1>256)=256;
    col_edge1(col_edge1<1)=1;
    row_edge1(row_edge1>256)=256;
    row_edge1(row_edge1<1)=1;
    % 每个PSF的初始位置
    centers = stats.Centroid;
    
    sto=ones(sizeg,1);
    % 设置初始位点
    startpoint_1=[centers(:,2),centers(:,1),1.2.*sto,1.2.*sto,stats.MaxIntensity*2*pi*1.2*1.2,zeros(sizeg,1),bg_mean.*sto];% 一一对应fittype中的参数
    % 每个点的积分值
    Integral1=zeros(sizeg,1);
    parfor i=1:sizeg
        imgGrayData = ima1_1(row_edge1(i,1):row_edge1(i,2),col_edge1(i,1):col_edge1(i,2));
        % 创建网格
        [X,Y]=meshgrid(col_edge1(i,1):col_edge1(i,2),row_edge1(i,1):row_edge1(i,2));
        % 拟合前准备
        [xData, yData, zData] = prepareSurfaceData(X,Y, imgGrayData);
        % 设置fittype
        ft=fittype(@(u1,u2,ss1,ss2,A,th,C,x,y)A*exp(-(cos(th)*x+sin(th)*y-u1).^2./(2*ss1^2)- ...
            (-sin(th)*x+cos(th)*y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2)+C,'independent', ...
            {'x','y'},'dependent','z');
        % 设置startpoint
        max1=max(imgGrayData(:));% 峰值
        % 用fit函数进行拟合
        fitResult=fit([xData,yData],zData,ft,'startpoint',startpoint_1(i,:));
        % 拟合得到的数据存入fr_D
        fr_D1(i,:)=coeffvalues(fitResult);
        % 将fr_D的数据分别存到七个变量中，与拟合一一对应
        numc=num2cell(fr_D1(i,:),1);
        [u1,u2,ss1,ss2,A,th,C]=deal(numc{:});
        % 计算每个数据点的积分强度
        Integral1(i)=integration_DF(fr_D1(i,:));
        % 转换矩阵，将拟合的u1和u2转变为中心点
        cc=[cos(th),sin(th);-sin(th),cos(th)];
        xy=[u1;u2];
        xxyy=cc^(-1)*xy;
        % 计算得到的中心点分别存入center_z/y中
        center_x(i)=xxyy(1);
        center_y(i)=xxyy(2);
        % figure
        % f=surf(X,Y,imgGrayData);
        % shading interp
        % plot(fitResult,[xData,yData],zData)
    end

    % 将数据导入到一个Excel表中
    % 计算每个PSF的半峰宽FWHM_x & FWHM_Y、中心点center_x/y

    fr_D1(:,1:2)=[center_x center_y];
    % 每个数据点的半径
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = diameters/2+1;
    % 只画每一组数据的第一帧
    if m==1
        % 绘图
        figure
        imagesc(ima1);colormap('gray');axis square;colorbar;
        % 圈出每个数据点
        hold on
        viscircles(centers,radii,'color','#0072BD','LineWidth',1.5,'EnhanceVisibility',false);
        % plot(row_edge(:,1),col_edge(:,1),'.','Color','red')
        hold off
        % 显示每个图像的名字
        title(name_ima);
    end
    
    n=fr_D1;  
    [ss1,ss2]=deal(n(:,3),n(:,4));
    FWHM(:,1)=2.*ss1.*sqrt(-2*log(1/2));% X轴的半峰宽
    FWHM(:,2)=2.*ss2.*sqrt(-2*log(1/2));% Y轴的半峰宽
    stats.FWHM_X=FWHM(:,1);
    stats.FWHM_y=FWHM(:,2);
    stats.Integral=Integral1;
    % 计算偏转角度Theta
    stats.Theta=fr_D1(:,6).*360./pi;
    % Excel_save(stats,filename);
    % writematrix(Integral,filename,"WriteMode","append");

end

end

end
toc
delete(gcp);


%%
path='E:\data\cy5_cy5p5_intensity\20230621\';
% 对数据进行拟合
files = dir(fullfile(path, '*.*'));
% 获取该文件下的每个数据文件
datanames = {files(~[files.isdir]).name}';
% 获取每个文件的格式
[~, ~, fileExtension] = fileparts(datanames);
expectedFormat='.xlsx';
% 选取.xlsx文件
datanames=datanames(strcmpi(fileExtension, expectedFormat));
% 获取数据的数量
datanumber=size(datanames,1);
% 初始化数据大小
intensity=zeros(datanumber,1);
error=zeros(datanumber,1);
% 对数据进行拟合
for i=1:1
    dataname=cell2mat(datanames(i,:));
    datapath=strcat(path,dataname)
    data_UCNP=readmatrix(datapath);
    [intensity(i),error(i)]=dataanalysis(data_UCNP);
end
datat=table(datanames,intensity,error);
data_ana_name=[path, 'data analysis.xlsx'];
% 删除进程池
%% 双边滤波
r = 3;% 滤波半径
a = 1;% 全局方差
b = 0.1;% 局部方差
[x,y] = meshgrid(-r:r);
w1 = exp(-(x.^2+y.^2)/(2*a^2));
% h = waitbar(0,'Applying bilateral filter...');
% set(h,'Name','Bilateral Filter Progress');
d10=d(:,:,10);
[m,n] = size(d10);
f_temp = padarray(d10,[r r],'symmetric');
g = zeros(m,n);
for i = r+1:m+r
    for j = r+1:n+r
        temp = f_temp(i-r:i+r,j-r:j+r);
        w2 = exp(-(temp-d10(i-r,j-r)).^2/(2*b^2));
        w = w1.*w2;
        s = temp.*w;
        g(i-r ...
            ,j-r) = sum(s(:))/sum(w(:));
    end
%     waitbar((i-r)/m);
end
% g = revertclass(g);

%%
parfor x=1:5
    % 转变路径格式为char
    dataname_mat=cell2mat(datanames(x,:))
    % 设置每个文件的
    data_save_name=[path,dataname_mat(1:end-4) '.xlsx'];
    fig_save_name=[path,dataname_mat(1:end-4) '.png'];
end





