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
% imagesc(seq.imageData(:,:,10))
for m=1:3
    % 显示现在正在处理的数据
    name_ima=[ char(qq(q)) 'mA' ' fov' char(string(x))  ' frame' char(string(m)) ];
    disp(['现在是' name_ima '...'])
    ima1=double(seq.imageData(:,:,m));
    % 选取数据范围
    ima1=ima1(100:355,100:355);
    % 根据CCD比例，对数据处理
    ima1=(ima1.*5.75)./300;
    % 数据转置
    ima1_1=ima1';
    % 二值化
    S=rescale(ima1);
    % 确定并去除背景值
    bg=imopen(S,strel('disk',4));
    S=S-bg;
    % find the threshold
    level=graythresh(S);
    bw=imbinarize(S,level);
    % 当图像比较乱时，可以通过该方法起到一定效果
    % bw = bwmorph(bw,"clean");
    % bw = bwmorph(bw,"open");
    [label,object]=bwlabel(bw,4);
    % 将regionprops函数得到的数据存入到stats
    stats=regionprops('table',label,ima1,'Centroid','Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');
    num_d=find(stats.Area>30);% UCNP
    % 找寻多个光斑连到一起的区域，重复上述运算，并将其分开
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
    % stats(stats.MaxIntensity<0.8*mean(stats.MaxIntensity),:)=[];
    stats(stats.MaxIntensity>2.5*mean(stats.MaxIntensity),:)=[];
    stats(stats.Area<5,:)=[];
    % 数据数量
    sizeg=size(stats,1);
    fr_D=zeros(sizeg,7);
    % 初始化FWHM
    FWHM=zeros(sizeg,2);
    % 初始化中心点位置
    center_x=zeros(sizeg,1);
    center_y=zeros(sizeg,1);
    % 图像的背景均值
    bg_mean=mean(ima1_1(:));
    % 确定每个数据点的范围
    x_edge=stats.BoundingBox(:,1:4);
    x_edge(:,1:2)=x_edge(:,1:2)+0.5;
    % 每个区域的row参数
    row_edge=[x_edge(:,1)-1 x_edge(:,1)+x_edge(:,3)+1];
    % 每个区域的col参数
    col_edge=[x_edge(:,2)-1 x_edge(:,2)+x_edge(:,4)+1];
    % 去除超出数组范围的参数
    col_edge(col_edge>256)=256;
    col_edge(col_edge<1)=1;
    row_edge(row_edge>256)=256;
    row_edge(row_edge<1)=1;
    % 每个PSF的初始位置
    centers = stats.Centroid;
    sto=ones(sizeg,1);
    % 设置拟合初始位点 sigma=1.2
    startpoint_G=[centers(:,2),centers(:,1),1.2.*sto,1.2.*sto,stats.MaxIntensity*2*pi*1.2*1.2,zeros(sizeg,1),bg_mean.*sto];% 一一对应fittype中的参数
    % 初始化每个点的积分值
    Integral=zeros(sizeg,1);
    parfor i=1:sizeg
        % 根据范围选取数据
        imgGrayData = ima1_1(row_edge(i,1):row_edge(i,2),col_edge(i,1):col_edge(i,2));
        % 创建网格
        [X,Y]=meshgrid(col_edge(i,1):col_edge(i,2),row_edge(i,1):row_edge(i,2));
        % 拟合前准备
        [xData, yData, zData] = prepareSurfaceData( X,Y, imgGrayData);
        % 设置fittype
        ft=fittype(@(u1,u2,ss1,ss2,A,th,C,x,y)A*exp(-(cos(th)*x+sin(th)*y-u1).^2./(2*ss1^2)- ...
            (-sin(th)*x+cos(th)*y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2)+C,'independent', ...
            {'x','y'},'dependent','z');
        % 设置startpoint
        % 用fit函数进行拟合
        fitResult=fit([xData,yData],zData,ft,'startpoint',startpoint_G(i,:));
        % 拟合得到的数据存入fr_D
        fr_D(i,:)=coeffvalues(fitResult);
        % 然后将拟合数据的每个参数分开
        numc=num2cell(fr_D(i,:),1);
        [u1,u2,ss1,ss2,A,th,~]=deal(numc{:});
        % 转换矩阵，将拟合的u1和u2转变为中心点
        cc=[cos(th),sin(th);-sin(th),cos(th)];
        xy=[u1;u2];
        xxyy=cc^(-1)*xy;
        % 计算得到的中心点分别存入center_x/y中
        center_x(i)=xxyy(1);
        center_y(i)=xxyy(2);
        % 计算每个数据点的积分强度
        % 先确定函数的中心点，并将其带入到函数中
        [u1,u2]=deal(xxyy(1),xxyy(2));
        % 先确定积分函数的形式
        f=@(x,y)A*exp(-(x-u1).^2./(2*ss1^2)-(y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2);
        % 选取-10/+10σ的范围进行积分
        Integral(i)=integral2(f,u1-10*ss1,u1+10*ss1,u2-10*ss2,u2+10*ss2);
        % 对每次拟合的原始数据和拟合结果进行绘图
        % figure
        % f=surf(X,Y,imgGrayData);
        % shading interp
        % plot(fitResult,[xData,yData],zData)
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
            % plot(fr_D(:,2),fr_D(:,1),'.','Color','red')
            hold off
            % 显示每个图像的名字
            title(name_ima);
        end
    end

    fr_D(:,1:2)=[center_x center_y];
    % 每个数据点的半径
    
    
    n=fr_D;  
    [ss1,ss2]=deal(n(:,3),n(:,4));
    FWHM(:,1)=2.*ss1.*sqrt(-2*log(1/2));% X轴的半峰宽
    FWHM(:,2)=2.*ss2.*sqrt(-2*log(1/2));% Y轴的半峰宽
    stats.FWHM_X=FWHM(:,1);
    stats.FWHM_y=FWHM(:,2);
    stats.Integral=Integral./seq.exposureTime;
    % 计算偏转角度Theta
    stats.Theta=fr_D(:,6).*360./pi;
    % Excel_save(stats,filename);
    % writematrix(Integral,filename,"WriteMode","append");
end

end

end
toc
delete(gcp);
%%
% laser power 980
% power=[500 300 188 128 100 85];
% laser power 808
power=[590 460 372 308 260 230];
% 初始化误差值数组
error_data=zeros(6,1);
y=zeros(6,1);
% load data
for i=1:6
    pathread='E:\data\sm_l0308\20230302_808\UCNP_Nd_ZTL_1wX_808\';
    nameread=['data_UCNP_' char(string(power(i))) '.xlsx'];
    filenamer=[pathread,nameread];
    data_UCNP=readmatrix(filenamer);
    data_UCNP(data_UCNP(:)>5*mean(data_UCNP(:)))=[];
    data_UCNP(data_UCNP(:)<0)=[];
    figure
    h=histogram(data_UCNP);
    % 每个柱的y值
    ycount=histcounts(data_UCNP,h.BinEdges);
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
    histogram(data_UCNP)
    hold on
    plot(fitResult)
    hold off
    error_data(i)=std(data_UCNP);
    xdata=coeffvalues(fitResult);
    y(i)=xdata(2);
end
% power intensity 980nm
% x=[21710 11613 5894 2835 1417 680];
% power intensity 808nm
x=[21710 14448 9611 6391 4252 2835];
figure
hold on
e=errorbar(x,y,error_data./2,'LineWidth',1);
e.Marker='.';
e.MarkerSize=10;
hold off
grid on;
set(gca,'xscale','log','yscale','log')
ax=gca;outerpos = ax.OuterPosition;% 去除白框
xlabel('Power Intensity (W/cm^2)');
ylabel('Photo luminescence (photons/s)');
legend('UCNP-Nd-ZTL-1wX-808',Location='northwest');
legend(Box="off")


%%



