clear
clc
tic
path='G:\My work\data\trace\Cy3_Cy5_fret\20241115\';
%name='std_532nm_10frames_0p1s_8mW_fov14.sif';
calibration=DCcalibration(path);
toc


%%
% if isempty(gcp('nocreate'))
%     pnumber=parpool(4);
% end
clear;clc
path='G:\My work\data\trace\Cy3_Cy5_fret\20241115\';
files = dir(fullfile(path, '*532*.sif*'));
% 获取该文件下的每个数据文件
datanames = {files(~[files.isdir]).name}';
% 获取数据的数量
datanumber=size(datanames,1);

for num=3:3
    % 转变路径格式为char
    dataname_mat=cell2mat(datanames(num,:));
    % 设置每个文件的
    data_save_name=[path,dataname_mat(1:end-4)];
    % 设置文件读取的路径
    pathname=strcat(path,dataname_mat);
    % 读取该文件
    seq=readsif(pathname);
    framesize=size(seq.imageData,3);
    % 确定开始采数据的位置
    startframe=1;
    for i=1:framesize
        intmean=mean(mean(seq.imageData(:,:,i)));
        if intmean>500
            startframe=i;
            break
        end
    end
    % 显示现在正在处理的数据
    disp(['现在是 ' dataname_mat(1:end-4) '...'])
    dataedge=128;
    seq.imageData=seq.imageData.*5.75./300;
    % red=double(seq.imageData(dataedge:255+dataedge,257:512,startframe:end));
    % green=double(seq.imageData(dataedge:255+dataedge,1:256,startframe:end));
    green=double(seq.imageData(dataedge:255+dataedge,257:512,startframe:end));
    red=double(seq.imageData(dataedge:255+dataedge,1:256,startframe:end));
    loc_red=loc_db(red,50);
    filter=loc_red(:,1)>250|loc_red(:,1)<7|loc_red(:,2)>250|loc_red(:,2)<7;
    loc_red(filter,:)=[];
    figure;
    imagesc(red(:,:,1));axis square;colorbar;colormap('gray')
    % 圈出每个数据点
    hold on
    viscircles(loc_red,3.*ones(size(loc_red,1),1),'color','#4F94CD','LineWidth',1.5,'EnhanceVisibility',false);
    plot(loc_red(:,1),loc_red(:,2),'Color','red','Marker','.','LineStyle','none')
    % plot(centers(:,1),centers(:,2),'Color','red','Marker','.','LineStyle','none')
    hold off
    % 根据校准计算方法算出的三次校准曲线，计算绿光通道的信号位置
    data=seq.imageData(dataedge:255+dataedge,:,startframe:end);
    calibration=readmatrix([path 'calibration_result.xlsx']);
    loc_green=[];
    loc_red_x= loc_red(:,1);
    loc_red_y = loc_red(:,2);
    % 创建非线性转换矩阵
    % 构造设计矩阵
    X = [ones(size(loc_red_x)), loc_red_x, loc_red_y, loc_red_x.^2, loc_red_y.^2, loc_red_x.*loc_red_y, loc_red_x.^3, loc_red_y.^3, loc_red_x.^2.*loc_red_y, loc_red_x.*loc_red_y.^2];
    loc_green(:,1)=X*calibration(:,1);
    loc_green(:,2)=X*calibration(:,2);
    filter=loc_green(:,1)>250|loc_green(:,1)<7|loc_green(:,2)>250|loc_green(:,2)<7;
    loc_red(filter,:)=[];
    loc_green(filter,:)=[];
    % 初始化trace储存
  
            
            Intensity_red = time_trace(loc_red(:,1:2), red(:,:,1:end));
            Intensity_green = time_trace(loc_green(:,1:2), green(:,:,1:end));
    
            
    end
    Intensity_red=Intensity_red-median(median(Intensity_red(:,end-10,end),2));
    Intensity_green=Intensity_green-median(median(Intensity_green(:,end-10,end),2));
    writematrix(Intensity_green',[data_save_name,'_Trace_Green.xlsx'])
    writematrix(Intensity_red',[data_save_name,'_Trace_Red.xlsx'])
end


%%
clear;clc
path='E:\data\20240704\Cy3Cy55\';
files = dir(fullfile(path, '*532*.sif*'));
% 获取该文件下的每个数据文件
datanames = {files(~[files.isdir]).name}';
% 获取数据的数量
datanumber=size(datanames,1);

clear Intensity_green Intensity_red
dataname_mat=cell2mat(datanames(1,:));
% 设置每个文件的
data_save_name=[path,dataname_mat(1:end-4)];
pathname=strcat(path,dataname_mat)
seq=readsif(pathname);
dataedge=128;
for i=1:100
    intmean=mean(mean(seq.imageData(:,:,i)));
    if intmean>500
        startframe=i;
        break
    end
end
data=seq.imageData(dataedge:255+dataedge,:,startframe:end);
Intensity_red=readmatrix([data_save_name,'_Trace_Red.xlsx']);
Intensity_green=readmatrix([data_save_name,'_Trace_Green.xlsx']);
Intensity_red=Intensity_red';
Intensity_green=Intensity_green';
%FRET=Intensity_red./(Intensity_green+Intensity_green);
%%
clear state_red likest_red state_green likest_green
for i=1:size(Intensity_red,1)
    disp(['processing point: ',num2str(i),'/' num2str(size(Intensity_red,1))]);
    [state_red(i,:),likest_red(i,:),~]=trace_HMM(Intensity_red(i,1:492),seq.exposureTime);
    [state_green(i,:),likest_green(i,:),~]=trace_HMM(Intensity_green(i,1:492),seq.exposureTime);
end
% filter=find(state_red==2);
% Intensity_red=Intensity_red(filter,:);
% Intensity_green=Intensity_green(filter,:);
% likest_red=likest_red(filter,:);
% likest_green=likest_green(filter,:);

%% 先简单计算得到 Survival time
for i=1:size(Intensity_red,1)
    x=i;
    st(x,:)=size(Intensity_red,2)-sum(likest_red(x,:)==min(likest_red(x,:)));
end
st=st.*seq.exposureTime;
figure
histogram(st)
%% 手动筛选数据，设为0
st(st==0)=[];
histogram(st)
writematrix(st,[data_save_name,'_Survival Time.xlsx'])
%%
time=size(likest_red,2);
time=3000;
timescale=seq.exposureTime;
for i=1:69
    figure
    title(i)
    hold on;
    plot(0:timescale:(time-1)*timescale,Intensity_green(i,1:time),Color='green',LineWidth=1.5);
    %plot(0:timescale:(time-1)*timescale,likest_green(i,1:time),Color='black',LineWidth=1.5);
    plot(0:timescale:(time-1)*timescale,Intensity_red(i,1:time),Color='red',LineWidth=1.5)
    %plot(0:timescale:(time-1)*timescale,likest_red(i,1:time),Color='blue',LineWidth=1.5)
    hold off;
end
%%
for i=1:size(likest_red,1)
    if likest_green(i,1)<0
        num_green=find(likest_green(i,:)==min(likest_green(i,:)));
        likest_green(i,:)=likest_green(i,:)-median(Intensity_green(i,num_green));
        Intensity_green(i,:)=Intensity_green(i,:)-median(Intensity_green(i,num_green));
    end
    if min(likest_red(:))<0
        likest_red(i,:)=likest_red(i,:)-min(Intensity_red(i,:));
        Intensity_red(i,:)=Intensity_red(i,:)-min(Intensity_red(i,:));
    end
    if likest_green(i,1)<100
      FRET(i,3)=likest_red(i,1)/(likest_red(i,1)+likest_green(i,1));
      FRET(i,1)=likest_green(i,1);
      FRET(i,2)=likest_red(i,1);
    else
        try
        num=max(find(likest_red(i,:)>2*abs(min(likest_red(i,:)))));
        FRET(i,3)=likest_red(i,num)/(likest_red(i,num)+likest_green(i,num));
        FRET(i,2)=likest_red(i,num);
        FRET(i,1)=likest_green(i,num);
        catch
            i
        end
    end
end
figure
histogram(FRET(:,3))
%% 根据手筛的数据筛除数据
FRET(FRET(:,3)==0,:)=[];
FRET(:,3)=FRET(:,2)./(FRET(:,2)+FRET(:,1));
histogram(FRET(:,3));
mean(FRET(:,3))
% ff=1./(1+gama_df*(FRET(:,1)./FRET(:,2)));% 根据γ计算得到的FRET效率
%% 计算仪器的γ值
num_gama=intersect(find(state_green<4),find(state_red<4));
for m=1:size(num_gama,1)
    figure
    i=num_gama(m);
    title(i)
    hold on;
    plot(0:timescale:(time-1)*timescale,Intensity_green(i,1:time),Color='green',LineWidth=1.5);
    plot(0:timescale:(time-1)*timescale,likest_green(i,1:time),Color='black',LineWidth=1.5);
    plot(0:timescale:(time-1)*timescale,Intensity_red(i,1:time),Color='red',LineWidth=1.5)
    plot(0:timescale:(time-1)*timescale,likest_red(i,1:time),Color='blue',LineWidth=1.5)
    hold off;
end
for m=1:size(num_gama,1)
    i=num_gama(m);
    delta_donor(m,3)=max(likest_green(i,:))-min(likest_green(i,:));
    delta_donor(m,1)=max(likest_green(i,:));
    delta_donor(m,2)=min(likest_green(i,:));
    delta_acceptor(m,3)=max(likest_red(i,:))-min(likest_red(i,:));
    delta_acceptor(m,1)=max(likest_red(i,:));
    delta_acceptor(m,2)=min(likest_red(i,:));
end
gama_FRET=delta_acceptor(:,3)./delta_donor(:,3);
%%
histogram(FRET(:,3))
writematrix(FRET,[data_save_name,'_FRET.xlsx'])
%% 为绘图做准备，计算得到每个PSF坐标
for i=1:100
    intmean=mean(mean(seq.imageData(:,:,i)));
    if intmean>500
        startframe=i;
        break
    end
end
data=seq.imageData(dataedge:255+dataedge,:,startframe:end);
red=double(seq.imageData(dataedge:255+dataedge,257:512,startframe:end));
green=double(seq.imageData(dataedge:255+dataedge,1:256,startframe:end));
calibration=readmatrix([path 'calibration_result.xlsx']);
loc_red=loc_db(red,50);
filter=loc_red(:,1)>250|loc_red(:,1)<7|loc_red(:,2)>250|loc_red(:,2)<7;
loc_red(filter,:)=[];
loc_green=[];
loc_red_x= loc_red(:,1);
loc_red_y = loc_red(:,2);
% 创建非线性转换矩阵
% 构造设计矩阵
X = [ones(size(loc_red_x)), loc_red_x, loc_red_y, loc_red_x.^2, loc_red_y.^2, loc_red_x.*loc_red_y, loc_red_x.^3, loc_red_y.^3, loc_red_x.^2.*loc_red_y, loc_red_x.*loc_red_y.^2];
loc_green(:,1)=X*calibration(:,1);
loc_green(:,2)=X*calibration(:,2);
filter=loc_green(:,1)>250|loc_green(:,1)<7|loc_green(:,2)>250|loc_green(:,2)<7;
loc_red(filter,:)=[];
loc_green(filter,:)=[];

%% 绘图 查看具体数据
% edge=3;
% figure
% set(gcf,'outerposition',get(0,'screensize'));
% num=130;
% for i=106:106
%     %figure
%     subplot(122)
%     imagesc(red(:,:,i));
% %     hold on
% %     viscircles([loc_red(:,1)+256,loc_red(:,2)],radii,'color','red','LineWidth',1.5,'EnhanceVisibility',false);
% %     plot(loc_red(:,1)+256,loc_red(:,2),'Color','red','Marker','.','LineStyle','none')
% %     viscircles(loc_green,radii,'color','green','LineWidth',1.5,'EnhanceVisibility',false);
% %     plot(loc_green(:,1),loc_green(:,2),'Color','green','Marker','.','LineStyle','none')
% %     hold off
% %     pause(1)
%     hold on
%     axis square
%     title('Red')
%     viscircles([loc_red(num,1),loc_red(num,2)],edge.*ones(size(num,1),1),'color','red','LineWidth',1.5,'EnhanceVisibility',false);
%     plot(loc_red(num,1),loc_red(num,2),'Color','red','Marker','.','LineStyle','none')
%     text(loc_red(num,1),loc_red(num,2),char(string(num)),'Color','red')
%     hold off
% 
%     subplot(121)
%     imagesc(green(:,:,i))
%     title('Green')
%     axis square
%     hold on
%     viscircles(loc_green(num,:),edge.*ones(size(num,1),1),'color','green','LineWidth',1.5,'EnhanceVisibility',false);
%     plot(loc_green(num,1),loc_green(num,2),'Color','green','Marker','.','LineStyle','none')
%     text(loc_red(num,1),loc_red(num,2),char(string(num)),'Color','green')
%     hold off
%     pause(1)
% end
