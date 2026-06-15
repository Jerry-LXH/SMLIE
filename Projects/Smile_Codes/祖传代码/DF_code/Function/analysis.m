function [finish] = analysis(namein)
pathin = 'F:\\20220331hcy\\';
pathout1 = 'F:\\20220331hcy_analysis\\fig1\\';
pathout2 = 'F:\\20220331hcy_analysis\\data1\\';

seq=readsif([pathin namein]);

% 多进程处理
parpool(6);
for m=614:2000
    tic
    % 显示现在正在处理的数据
    disp(['第',num2str(m),'帧处理中……']);
    ima1=double(seq.imageData(:,:,m));
    ima1=ima1(100:355,100:355);
    ima1=(ima1.*5.75)./300;
    % 数据转置
    ima1_1=ima1';
    % 二值化
    S=rescale(ima1);
    bg=imopen(S,strel('disk',6));
    S=S-bg;
    % find the threshold
    level=graythresh(S);
    bw=imbinarize(S,level);
    bw = bwmorph(bw,"clean");
%     % 当图像比较乱时，可以通过该方法起到一定效果
    bw = bwmorph(bw,"open");
    [label,object]=bwlabel(bw,8);
    % 将regionprops函数得到的数据存入到stats
    stats=regionprops('table',label,ima1,'Centroid',...
        'Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');
    num_d=find(stats.Area>30);% UCNP
    % 找寻两个光斑连到一起的区域，并且分开
    try
    if (isempty(num_d)==0)
        for i=1:size(num_d,1)
            state_d=stats(num_d(i),:);
            x_d=state_d.BoundingBox;
            state_a=p_d(x_d,ima1_1);
            stats=[stats;state_a];
        end
    end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:badsubscript')
        % 处理越界错误，例如输出错误信息并跳过当前迭代
            disp(['发生越界错误：' ME.message])
            continue
        else
        % 如果不是越界错误，则重新抛出异常
            rethrow(ME)
        end
    end
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
    fr_D=zeros(sizeg,7);
    FWHM=zeros(sizeg,2);
    % 图像的背景均值
    bg_mean=mean(ima1_1(:));
    x_edge=stats.BoundingBox(:,1:4);
    x_edge(:,1:2)=x_edge(:,1:2)+0.5;
    % 区域的row参数
    row_edge=[x_edge(:,1)-1 x_edge(:,1)+x_edge(:,3)+1];
    % 区域的col参数
    col_edge=[x_edge(:,2)-1 x_edge(:,2)+x_edge(:,4)+1];
    % 去除超出数组范围的参数
    col_edge(col_edge>256)=256;
    col_edge(col_edge<1)=1;
    row_edge(row_edge>256)=256;
    row_edge(row_edge<1)=1;
    % 每个PSF的初始位置
    centers = stats.Centroid;
    parfor i=1:sizeg
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
        max1=max(imgGrayData(:));% 峰值
        % sigmax=1.2;sigmay=1.2;
        startpoint_G=[centers(i,:),1.2,1.2,max1*2*pi*1.2*1.2,0,bg_mean];% 一一对应fittype中的参数
        % 用fit函数进行拟合
        fitResult=fit([xData,yData],zData,ft,'startpoint',startpoint_G);
        % 拟合得到的数据存入fr_D
        % figure
        % f=surf(X,Y,imgGrayData);
        % smooth(f)
        % shading interp
        % plot(fitResult,[xData,yData],zData)
        fr_D(i,:)=coeffvalues(fitResult);

    end

    % 将数据导入到一个Excel表中
    % 计算每个PSF的半峰宽FWHM_x & FWHM_Y、中心点center_x/y
    center_x=zeros(sizeg,1);
    center_y=zeros(sizeg,1);

    parfor i=1:sizeg
        n=fr_D; 
        th=n(i,6);
        cc=[cos(th),sin(th);-sin(th),cos(th)];
        xy=[n(i,1);n(i,2)];
        xxyy=cc^(-1)*xy;
        center_x(i)=xxyy(1);
        center_y(i)=xxyy(2);
        % 计算每组数据的积分
        Integral(i)=integration_DF(n(i,:));
    end
    % 
    fr_D(:,1:2)=[center_x center_y];
    % 每个数据点的半径
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = diameters/2+1;
    % 保存画圆前的图像
    %figure
    imagesc(ima1)
    colormap('gray')
    axis square
    colorbar
    filename1 = strcat(pathout1, 'before', num2str(m), '.png');
    saveas(gcf, filename1);

    % 保存画圆后的图像
    %figure
    imagesc(ima1)
    colormap('gray')
    axis square
    colorbar
    hold on
    viscircles(centers,radii,'color','#0072BD','LineWidth',1.5,'EnhanceVisibility',false);
    hold off
    filename2 = strcat(pathout1, 'after',num2str(m), '.png');
    saveas(gcf, filename2);
    end
    
    n=fr_D;  
    [ss1,ss2]=deal(n(:,3),n(:,4));
    FWHM(:,1)=2.*ss1.*sqrt(-2*log(1/2));% X轴的半峰宽
    FWHM(:,2)=2.*ss2.*sqrt(-2*log(1/2));% Y轴的半峰宽
    stats.FWHM_X=FWHM(:,1);
    stats.FWHM_y=FWHM(:,2);
    stats.Integral=Integral;
     % 计算偏转角度Theta
    stats.Theta=fr_D(:,6).*360./pi;
    filename3 = strcat(pathout2, ['data', num2str(m), '.xlsx']);
    Excel_save(stats,filename3,'writemode','append');
    disp(['完成第',num2str(m),'帧处理。']);
    end

% delete(gcp);
% finish= 'Finished';  

%end



% 
% parfor m=1:2000
%     tic
%     % 显示现在正在处理的数据
%     disp(['第',num2str(m),'帧处理中……']);
%     ima1=double(seq.imageData(:,:,m));
%     ima1=ima1(100:355,100:355);
%     ima1=(ima1.*5.75)./300;
%     % 数据转置
%     ima1_1=ima1';
%     % 二值化
%     S=rescale(ima1);
%     bg=imopen(S,strel('disk',6));  % 此处4可以调参
%     S=S-bg;
%     % find the threshold
%     level=graythresh(S);
%     bw=imbinarize(S,level);
% %     bw = bwmorph(bw,"clean");
% %     % 当图像比较乱时，可以通过该方法起到一定效果
% %     bw = bwmorph(bw,"open");
%     [label,object]=bwlabel(bw,4);
%     % 将regionprops函数得到的数据存入到stats
%     stats=regionprops('table',label,ima1,'Centroid',...
%         'Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');
%     num_d=find(stats.Area>30);           % 此处30可用于调参，指保存的光斑的最小面积
%     % 找寻两个光斑连到一起的区域，并且分开
%     try
%     if (isempty(num_d)==0)
%         for i=1:size(num_d,1)
%             state_d=stats(num_d(i),:);
%             x_d=state_d.BoundingBox;
%             state_a=p_d(x_d,ima1_1);
%             stats=[stats;state_a];
%         end
%     end
%     catch ME
%         if strcmp(ME.identifier, 'MATLAB:badsubscript')
%         % 处理越界错误，例如输出错误信息并跳过当前迭代
%             disp(['发生越界错误：' ME.message])
%             continue
%         else
%         % 如果不是越界错误，则重新抛出异常
%             rethrow(ME)
%         end
%     end
%     if (isempty(num_d)==0)
%         for i=1:size(num_d,1)
%             state_d=stats(num_d(i),:);
%             x_d=state_d.BoundingBox;
%             state_a=p_d(x_d,ima1_1);
%             stats=[stats;state_a];
%         end
%     end
%     % 删去光斑连在一起的数组
%     stats(num_d,:)=[];
% %     stats(stats.MaxIntensity<0.8*mean(stats.MaxIntensity),:)=[];
%     stats(stats.MaxIntensity>2.5*mean(stats.MaxIntensity),:)=[];
%     stats(stats.Area<5,:)=[];
%     sizeg=size(stats,1);
%     fr_D=zeros(sizeg,7);
%     FWHM=zeros(sizeg,2);
%     % 图像的背景均值
%     bg_mean=mean(ima1_1(:));
%     x_edge=stats.BoundingBox(:,1:4);
%     x_edge(:,1:2)=x_edge(:,1:2)+0.5;
%     % 区域的row参数
%     row_edge=[x_edge(:,1)-1 x_edge(:,1)+x_edge(:,3)+1];
%     % 区域的col参数
%     col_edge=[x_edge(:,2)-1 x_edge(:,2)+x_edge(:,4)+1];
%     % 去除超出数组范围的参数
%     col_edge(col_edge>256)=256;
%     col_edge(col_edge<1)=1;
%     row_edge(row_edge>256)=256;
%     row_edge(row_edge<1)=1;
%     % 每个PSF的初始位置
%     centers = stats.Centroid;
%     for i=1:sizeg
%         imgGrayData = ima1_1(row_edge(i,1):row_edge(i,2),col_edge(i,1):col_edge(i,2));
%         % 创建网格
%         [X,Y]=meshgrid(col_edge(i,1):col_edge(i,2),row_edge(i,1):row_edge(i,2));
%         % 拟合前准备
%         [xData, yData, zData] = prepareSurfaceData( X,Y, imgGrayData);
%         % 设置fittype
%         ft=fittype(@(u1,u2,ss1,ss2,A,th,C,x,y)A*exp(-(cos(th)*x+sin(th)*y-u1).^2./(2*ss1^2)- ...
%             (-sin(th)*x+cos(th)*y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2)+C,'independent', ...
%             {'x','y'},'dependent','z');
%         % 设置startpoint
%         max1=max(imgGrayData(:));% 峰值
%         % sigmax=1.2;sigmay=1.2;
%         startpoint_G=[centers(i,:),1.2,1.2,max1*2*pi*1.2*1.2,0,bg_mean];% 一一对应fittype中的参数
%         % 用fit函数进行拟合
%         fitResult=fit([xData,yData],zData,ft,'startpoint',startpoint_G);
%         % 拟合得到的数据存入fr_D
%         % figure
%         % f=surf(X,Y,imgGrayData);
%         % smooth(f)
%         % shading interp
%         % plot(fitResult,[xData,yData],zData)
%         fr_D(i,:)=coeffvalues(fitResult);
%         % 对每个PSF进行积分
%         [fr_D(i,1),fr_D(i,2)]=gfcenter(fr_D(i,1),fr_D(i,2),fr_D(i,6));
%         Integral(i)=integration_DF(fr_D(i,:));
% 
%     diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
%     radii = diameters/2+1;
% 
%     % 保存画圆前的图像
%     %figure
%     imagesc(ima1)
%     colormap('gray')
%     axis square
%     colorbar
%     filename1 = strcat(pathout1, 'before', num2str(m), '.png');
%     saveas(gcf, filename1);
% 
%     % 保存画圆后的图像
%     %figure
%     imagesc(ima1)
%     colormap('gray')
%     axis square
%     colorbar
%     hold on
%     viscircles(centers,radii,'color','#0072BD','LineWidth',1.5,'EnhanceVisibility',false);
%     hold off
%     filename2 = strcat(pathout1, 'after',num2str(m), '.png');
%     saveas(gcf, filename2);
% 
%     % 将数据导入到一个Excel表中
%     % 计算每个PSF的半峰宽FWHM_x & FWHM_Y
%     n=fr_D;
%     [ss1,ss2]=deal(n(:,3),n(:,4));
%     FWHM(:,1)=2.*ss1.*sqrt(-2*log(1/2));% X轴的半峰宽
%     FWHM(:,2)=2.*ss2.*sqrt(-2*log(1/2));% Y轴的半峰宽
%     stats.FWHM_X=FWHM(:,1);
%     stats.FWHM_y=FWHM(:,2);
%     stats.Integral=Integral;
%     % 计算偏转角度Theta
%     stats.Theta=fr_D(:,6).*360./pi;
%     filename3 = strcat(pathout2, ['data', num2str(m), '.xlsx']);
%     Excel_save(stats,filename3);
%     disp(['完成第',num2str(m),'帧处理。']);
% end
% 
% delete(gcp);
% finish= 'Finished';  
% 
% end