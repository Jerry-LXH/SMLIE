clear;
clc
tic
% 多进程处理
if isempty(gcp('nocreate'))
    pnumber=parpool(4);
end

filepath='E:\data\cy5_cy5p5_intensity\20230621\';
files = dir(fullfile(filepath, '*.*'));
% 获取该文件下的每个数据文件
filenames = {files([files.isdir]).name}';
% 去除自带的'.'和'..'的特殊文件夹
filenames(ismember(filenames,{'.','..'})) = [];

for p=1:size(filenames,1)
    % 读取每个文件夹下的数据文件
    path=cell2mat(filenames(p));
    path=strcat(filepath,path,'\');
    files = dir(fullfile(path, '*.*'));
    % 获取该文件下的每个数据文件
    datanames = {files(~[files.isdir]).name}';
    % 获取每个文件的格式
    [~, ~, fileExtension] = fileparts(datanames);
    expectedFormat='.sif';
    % 选取sif文件
    datanames=datanames(strcmpi(fileExtension, expectedFormat));
    % 获取数据的数量
    datanumber=size(datanames,1);
    for x=1:datanumber
        % 转变路径格式为char
        dataname_mat=cell2mat(datanames(x,:));
        % 设置每个文件的
        data_save_name=[path,dataname_mat(1:end-16) '.xlsx'];
        % 设置文件读取的路径
        pathname=strcat(path,dataname_mat);
        % 读取该文件
        seq=readsif(pathname);
        % 帧数大小
        framesize=size(seq.imageData,3);
        for m=startframe:startframe+3
            % 显示现在正在处理的数据
            powernum=strfind(dataname_mat, 'mA');
            name_ima=[ dataname_mat(powernum-3:powernum-1) 'mA' ' fov' dataname_mat(end-4)  ' frame' char(string(m)) ];
            disp(['现在是' name_ima '...'])
            ima1=double(seq.imageData(:,:,m));
            % 选取数据范围
            ima1=ima1(100:355,100:355);
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
            % 记录每个点的半径
            diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
            radii = diameters/2+1;
            % 只画每一组数据的第一帧
            if m==1
                % 绘图
                figure
                imagesc(ima1);colormap('gray');axis square;colorbar;
                % 圈出每个数据点
                hold on
                viscircles(centers,radii,'color','blue','LineWidth',1.5,'EnhanceVisibility',false);
                % plot(fr_D(:,2),fr_D(:,1),'.','Color','red')
                hold off
                % 显示每个图像的名字
                title(name_ima);
            end
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
            writematrix(Integral,data_save_name,"WriteMode","append");
        end

    end
end
toc
delete(gcp);

%%
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
for i=1:datanumber
    dataname=cell2mat(datanames(i,:));
    datapath=strcat(path,dataname);
    data_UCNP=readmatrix(datapath);
    [intensity(i),error(i)]=dataanalysis(data_UCNP);
end
datat=table(datanames,intensity,error);
data_ana_name=[path, 'data analysis.xlsx'];
% writetable(datat,data_ana_name)

