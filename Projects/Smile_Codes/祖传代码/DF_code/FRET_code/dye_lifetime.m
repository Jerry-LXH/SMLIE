%%
clear;
clc
tic
% 多进程处理
% if isempty(gcp('nocreate'))
%     pnumber=parpool(4);
% end
% 读取每个文件夹下的数据文件
path='E:\data\20240704\Cy3Cy5\';
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
mkdir([path,'Video_tif\']);
for x=1:datanumber
    % 转变路径格式为char
    dataname_mat=cell2mat(datanames(x,:));
    % 设置每个文件的
    data_save_name=[path,dataname_mat(1:end-9) '_lifetime' '.xlsx'];
    % 设置文件读取的路径
    pathname=strcat(path,dataname_mat);
    % 读取该文件
    seq=readsif(pathname);
    framesize=size(seq.imageData,3);
    % 确定开始采数据的位置
    for i=1:framesize
        intmean=mean(mean(seq.imageData(:,:,i)));
        if intmean>500
            startframe=i;
            break
        end
    end

    % 显示现在正在处理的数据
    disp(['现在是 ' dataname_mat(1:end-4) '...'])
    ima1=double(seq.imageData(:,:,startframe));
    % 选取数据范围
    area_r=128;
    area_c=128;
    dataedge=128;
    %ima1=ima1(dataedge:255+dataedge,1:256);% green
    ima1=ima1(dataedge:255+dataedge,257:512); % red
    %ima1=ima1(area_r:area_r+255,area_c:area_c+255);
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
    bw = bwmorph(bw,"clean");
    bw = bwmorph(bw,"open");
    [label,object]=bwlabel(bw,4);
    % 将regionprops函数得到的数据存入到stats
    stats=regionprops('table',label,ima1,'Centroid','Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');
    stats(stats.MaxIntensity<2*mean(ima1(:)),:)=[];
    stats(stats.Area<5,:)=[];

    % 数据数量
    sizeg=size(stats,1);
    % 初始化存活时间存储
    lifetime= zeros(sizeg,1);
    % 确定每个数据点的范围
    x_edge=stats.BoundingBox(:,1:4);
    x_edge(:,1:2)=x_edge(:,1:2)+0.5;
    % 每个区域的row参数
    row_edge=[x_edge(:,1)-1 , x_edge(:,1)+x_edge(:,3)+1];
    % 每个区域的col参数
    col_edge=[x_edge(:,2)-1 , x_edge(:,2)+x_edge(:,4)+1];
    % 去除超出数组范围的参数
    col_edge(col_edge>256)=256;
    col_edge(col_edge<1)=1;
    row_edge(row_edge>256)=256;
    row_edge(row_edge<1)=1;
    % 每个PSF的初始位置
    centers = stats.Centroid;
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = diameters/2+1;

    imagedata = double(seq.imageData);
    imagedata=(imagedata.*5.75)./seq.gainDAC;
    imagesc(ima1);colormap('gray');axis square;colorbar;
    % 圈出每个数据点
    hold on
    viscircles(centers,radii,'color','#4F94CD','LineWidth',1.5,'EnhanceVisibility',false);
    hold off
    pause(1)
    dataedge=128;
    data_bg = mean(imagedata(dataedge:255+dataedge,257:512,end),"all");
    Intensity_red=[];
    for i = 1 : sizeg
        % 根据范围选取数据
        for j = startframe : framesize
            % data_perframe=imagedata(dataedge:255+dataedge,1:256,j);% green
            data_perframe = imagedata(dataedge:255+dataedge,257:512,j);% red
            data_perframe = data_perframe';
            imgGrayData = data_perframe(row_edge(i,1):row_edge(i,2),col_edge(i,1):col_edge(i,2));
            MaxIntensity = max(imgGrayData(:));
            Intensity_red(i,j)=sum(imgGrayData,"all")-sum(data_bg.*ones(size(imgGrayData)),"all");
            % imagesc(imgGrayData);
            % colorbar
            % pause(1)
            if (MaxIntensity<2*data_bg)
                %break
            end
            lifetime(i,:) = lifetime(i,:)+1;
        end
    end
    lifetime=lifetime.*seq.exposureTime;
    writematrix(lifetime,data_save_name,"WriteMode","append");
end
toc
%%
% % 初始化存活时间存储
% % 确定每个数据点的范围
% result=[];
% %for i = 1 : size(centers,1)
%     % 根据范围选取数据
%     for j = 1 : size(green,3)
%         centers_round_green=round(centers_green);
%         % data_perframe = imagedata(dataedge:255+dataedge,257:512,j);% red
%         imgGrayData_green = green(centers_round_green(i,2)-3:centers_round_green(i,2)+3,centers_round_green(i,1)-3:centers_round_green(i,1)+3,j);
%         if max(imgGrayData_green(:))>2.5*mean(green(:,:,end), 'all') +6*std2(green(:,:,end))
%             a=[i j];
%             result=[result;a];
%             break
%         end
%     end
% %end

%%
for i = 15 : 15
    % 根据范围选取数据
    for j = 1 : size(green,3)
        centers_round_green=round(centers_green);
        % data_perframe = imagedata(dataedge:255+dataedge,257:512,j);% red
        imgGrayData_green = green(centers_round_green(i,2)-3:centers_round_green(i,2)+3,centers_round_green(i,1)-3:centers_round_green(i,1)+3,j);
        [X,Y]=meshgrid(1:7,1:7);
        surf(X,Y,imgGrayData_green)
        pause(1)
    end
end


%%
figure
num=[10,11,12,17,27,14];
for i=1:10
    imagesc(data(:,:,i));
    set (gca,'position',[0.1,0.1,0.8,0.8] );
    daspect([1 1 1])
%     hold on
%     viscircles([loc_red(:,1)+256,loc_red(:,2)],radii,'color','red','LineWidth',1.5,'EnhanceVisibility',false);
%     plot(loc_red(:,1)+256,loc_red(:,2),'Color','red','Marker','.','LineStyle','none')
%     viscircles(loc_green,radii,'color','green','LineWidth',1.5,'EnhanceVisibility',false);
%     plot(loc_green(:,1),loc_green(:,2),'Color','green','Marker','.','LineStyle','none')
%     hold off
%     pause(5)
hold on
viscircles([loc_red(num,1)+256,loc_red(num,2)],edge.*ones(size(num,2),1),'color','red','LineWidth',1.5,'EnhanceVisibility',false);
plot(loc_red(num,1)+256,loc_red(num,2),'Color','red','Marker','.','LineStyle','none')

viscircles(loc_green(num,:),edge.*ones(size(num,2),1),'color','green','LineWidth',1.5,'EnhanceVisibility',false);
plot(loc_green(num,1),loc_green(num,2),'Color','green','Marker','.','LineStyle','none')
for j=1:size(loc_red,1)
    text(loc_red(j,1)+256,loc_red(j,2),char(string(j)),'Color','red')
    text(loc_green(j,1),loc_green(j,2),char(string(j)),"Color",'green')
end
hold off
pause(3)
end

%%
Intensity_red=zeros(size(loc_green,1),size(data,3));
Intensity_green=zeros(size(loc_green,1),size(data,3));
edge=3;
for i=1:size(loc_green,1)
    for j=1:size(red,3)
        centers_round_green=round(loc_green);
        centers_round_red=round(loc_red);
        imgGrayData_green = green(centers_round_green(i,2)-edge:centers_round_green(i,2)+edge,centers_round_green(i,1)-edge:centers_round_green(i,1)+edge,j);
        imgGrayData_red = red(centers_round_red(i,2)-edge:centers_round_red(i,2)+edge,centers_round_red(i,1)-edge:centers_round_red(i,1)+edge,j);
%         [X,Y]=meshgrid(1:9,1:9);
%         surf(X,Y,imgGrayData_green)
%         figure
%         surf(X,Y,imgGrayData_red)
        Intensity_red(i,j)=sum(imgGrayData_red(:))-sum(mean(red(:,:,end),'all').*ones(1+2*edge),"all");
        Intensity_green(i,j)=sum(imgGrayData_green(:))-sum(mean(green(:,:,end),'all').*ones(1+2*edge),"all");
    end
end

%%
time=1000;
for i=1:56
    figure;
    %plot(0:time-1,Intensity_green(i,1:time),Color='green',LineWidth=1.5);
    hold on;
    plot(0:time-1,Intensity_red(i,1:time),Color='red',LineWidth=1.5)
    hold off;
end
