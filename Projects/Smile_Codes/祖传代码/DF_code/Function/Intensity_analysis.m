function [stats_all]=Intensity_analysis(datapath, dataname_mat,paint)

seq = readsif([datapath dataname_mat]);
stats_all = table;
for m=1:size(seq.imageData,3)
    ima1=double(seq.imageData(:,:,m));
    % 选取数据范围
    if (size(ima1,1)==256)
        area_r=1;
        area_c=1;
    else
        area_r=128;
        area_c=128;
    end
    ima1=ima1(area_r:area_r+255,area_c:area_c+255);
    % 根据CCDfangda 比例，对数据处理
    ima1=(ima1.*5.75)./seq.gainDAC;
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
    [label,~]=bwlabel(bw,4);
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
    % 当识别的点比较多时，即识别错了
    if size(stats,1)>150
        stats(stats.MaxIntensity<1.5*mean(ima1(:)),:)=[];
        stats(stats.Area<5,:)=[];
    end
    % 去除重叠的值
    stats(stats.MaxIntensity>3*mean(stats.MaxIntensity),:)=[];
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
    % 记录每个点的半径
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = diameters/2+1;
    % 只画每一组数据的第一帧
    if (m==1 && paint==1)
        % 绘图
        fig=figure('Visible','off');
        imagesc(ima1);colormap('gray');axis square;colorbar;
        % 圈出每个数据点
        hold on
        viscircles(centers,radii,'color','#4F94CD','LineWidth',1.5,'EnhanceVisibility',false);
        % plot(fr_D(:,2),fr_D(:,1),'.','Color','red')
        hold off
        % 命名并保存
        title(dataname_mat(1:end-4),'Interpreter','none');
        % saveas(fig,fig_save_name);
        img_save_name=[datapath,dataname_mat(1:end-4) '.png'];
        saveas(fig,img_save_name);
        close(fig);
        % 保存原始图像
        fig=figure('Visible','off');
        imagesc(seq.imageData(:,:,m));colormap('gray');axis square;colorbar;
        hold on
        rectangle('Position',[area_r,area_c,255,255],'LineWidth',1,'Visible','on','EdgeColor','red');
        hold off
        % 显示每个图像的名字
        % 保存原始数据的图像
        title([dataname_mat(1:end-4) '_fit'],'Interpreter','none');
        img_save_name=[datapath,dataname_mat(1:end-4) '_orginal' '.png'];
        saveas(fig,img_save_name);
        close(fig)
    end

    for i=1:sizeg
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
        % 先确定函数的中心点，并将其代入到函数中
        [u1,u2]=deal(xxyy(1),xxyy(2));
        % 先确定积分函数的形式
        f=@(x,y)A*exp(-(x-u1).^2./(2*ss1^2)-(y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2);
        % 选取-10/+10σ的范围进行积分
        Integral(i)=integral2(f,u1-5*ss1,u1+5*ss1,u2-5*ss2,u2+5*ss2);
        % 对每次拟合的原始数据和拟合结果进行绘图
        % figure
        % f=surf(X,Y,imgGrayData);
        % shading interp
        % plot(fitResult,[xData,yData],zData)
    end

    fr_D(:,1:2)=[center_x center_y];
    % 每个数据点的半径
    n=fr_D;
    [ss1,ss2]=deal(n(:,3),n(:,4));
    FWHM(:,1)=2.*ss1.*sqrt(-2*log(1/2));% X轴的半峰宽
    FWHM(:,2)=2.*ss2.*sqrt(-2*log(1/2));% Y轴的半峰宽
    stats.FWHM_X=FWHM(:,1);
    stats.FWHM_y=FWHM(:,2);
    Integral=Integral./seq.exposureTime;
    stats.Integral=Integral;
    % 计算偏转角度Theta
    stats.Theta=fr_D(:,6).*360./pi;
    stats_all = [stats(:,[1,2,7,8,9,10]); stats_all];
end
