clear;clc
% 
% if isempty(gcp('nocreate'))
%     pnumber=parpool(4);
% end

path='E:\data\20241020-pH\Cy5\';
files = dir(fullfile(path, '*.sif'));
% 获取该文件下的每个数据文件
datanames = {files(~[files.isdir]).name}';
% 获取数据的数量
datanumber=size(datanames,1);
for num=18:datanumber
    dataname_mat=cell2mat(datanames(num,:));
    % 设置每个文件的
    data_save_name=[path,dataname_mat(1:end-4),'Intensity_surival time.xlsx'];
    % 设置文件读取的路径
    pathname=strcat(path,dataname_mat);
    % 读取该文件
    seq=readsif(pathname);
    framesize=size(seq.imageData,3);
    startframe=1;
    % 确定开始采数据的位置
    for i=1:framesize
        intmean=mean(mean(seq.imageData(:,:,i)));
        if intmean>500
            startframe=i;
            break
        end
    end
    disp(['现在是 ' dataname_mat(1:end-4) '...'])
    dataedge=125;
    seq.imageData=seq.imageData.*5.75./300;
    red=double(seq.imageData(dataedge:255+dataedge,257:512,startframe:end));
    data=double(seq.imageData(dataedge:255+dataedge,dataedge:255+dataedge,startframe:end));
    %data=double(seq.imageData(:,:,startframe:end));
    red=data;
    %green=double(seq.imageData(dataedge:255+dataedge,1:256,startframe:end));
    %red=green;
    loc=loc_db(red,50);
    filter=loc(:,1)>250|loc(:,1)<7|loc(:,2)>250|loc(:,2)<7;
    loc(filter,:)=[];
    Intensity_red=zeros(size(loc,1),size(red,3));
%     writematrix(Intensity_red,[path,dataname_mat(1:end-4),'Intensity_red.xlsx'],'WriteMode','append');
    edge=3;
    for j=1:size(red,3)
        if j<1000
        %loc_tem=loc_molecule_DF(red(:,:,j),1);
        %filter=loc_tem(:,1)>250|loc_tem(:,1)<7|loc_tem(:,2)>250|loc_tem(:,2)<7;
        %loc_tem(filter,:)=[];
        %loc=loc_shift_correction(loc,loc_tem,2);
        end
        centers_round_red=round(loc);
        for i=1:size(loc,1)
            imgGrayData_red = red(centers_round_red(i,2)-edge:centers_round_red(i,2)+edge,centers_round_red(i,1)-edge:centers_round_red(i,1)+edge,j);
            Intensity_red(i,j)=sum(imgGrayData_red(:))-sum(mean(red(:,:,end),'all').*ones(1+2*edge),"all");
        end    
       
    end 
    Intensity_red=Intensity_red-median(median(Intensity_red(:,end-10,end),2));
    writematrix(Intensity_red,data_save_name,'WriteMode','append');
end
%%
clear;clc
path='E:\data\20241013\Cy5B\';
files = dir(fullfile(path, '*.sif*'));
% 获取该文件下的每个数据文件
datanames = {files(~[files.isdir]).name}';
% 获取数据的数量
datanumber=size(datanames,1);

clear  Intensity_red
dataname_mat=cell2mat(datanames(6,:));
% 设置每个文件的
data_save_name=[path,dataname_mat(1:end-4)];
pathname=strcat(path,dataname_mat)              
seq=readsif(pathname);
dataedge=150;
for i=1:100
    intmean=mean(mean(seq.imageData(:,:,i)));
    if intmean>500
        startframe=i;
        break
    end
end
red=double(seq.imageData(dataedge:255+dataedge,257:512,startframe:end));
green=double(seq.imageData(dataedge:255+dataedge,1:256,startframe:end));
loc_red=loc_db(red,50);
filter=loc_red(:,1)>250|loc_red(:,1)<7|loc_red(:,2)>250|loc_red(:,2)<7;
loc_red(filter,:)=[];
Intensity_red=readmatrix([data_save_name,'Intensity_surival time.xlsx']);
%Intensity_red(:,501:end)=[];
for i=1:size(Intensity_red,1)
    disp(['processing point: ',num2str(i),'/' num2str(size(Intensity_red,1))]);
    [state_red(i,:),likest_red(i,:),~]=trace_HMM(Intensity_red(i,:),seq.exposureTime);
end
%%
time=size(Intensity_red,2);
%num=find(state_red<5);
timescale=seq.exposureTime;
for i=51:122
    %x=num(i);
    x=i;
    figure;
    title(x)
    hold on;
    plot(0:timescale:(time-1)*timescale,Intensity_red(i,1:time),Color='red',LineWidth=1.5)
    plot(0:timescale:(time-1)*timescale,likest_red(i,1:time),Color='blue',LineWidth=1.5)
    hold off;
end

%%

%chosen_n=[];
for i=1:200
    if gcf().Number~=1
        gcf().Number
        chosen_n=horzcat(chosen_n,gcf().Number+50);
        close gcf
    end
end
%%
st=zeros(size(state_red,1),1);
%num=find(state_red<5);
for i=1:size(state_red,1)
    %x=num(i);
    x=i;
    st(x,:)=size(Intensity_red,2)-sum(likest_red(x,:)==min(likest_red(x,:)));
end
st=st.*seq.exposureTime;
for i=1:size(chosen_n,2)
    st(chosen_n(1,i),:)=0;
end

% for i=1:size(state_red,1)
%     st(i,:)=size(Intensity_red,2)-sum(likest_red(i,:)==min(likest_red(i,:)));
% end
%%
st(st==0)=[];
mean(st)
histogram(st)
writematrix(st,[path,dataname_mat(1:end-4),'_Surival Time.xlsx']);

%%
edge=3;
figure
%set(gcf,'outerposition',get(0,'screensize'));
num=64;
for i=1
    %figure
    imagesc(red(:,:,i));
    hold on
    viscircles([loc_red(:,1),loc_red(:,2)],edge.*ones(size(loc_red,1),1),'color','red','LineWidth',1.5,'EnhanceVisibility',false);
    plot(loc_red(:,1),loc_red(:,2),'Color','red','Marker','.','LineStyle','none')
    %viscircles(loc_green,radii,'color','green','LineWidth',1.5,'EnhanceVisibility',false);
    %plot(loc_green(:,1),loc_green(:,2),'Color','green','Marker','.','LineStyle','none')
    hold off
    % hold on
    % axis square
    % title('Red')
    % viscircles([loc_red(num,1),loc_red(num,2)],edge.*ones(size(num,1),1),'color','red','LineWidth',1.5,'EnhanceVisibility',false);
    % plot(loc_red(num,1),loc_red(num,2),'Color','red','Marker','.','LineStyle','none')
    % text(loc_red(num,1),loc_red(num,2),char(string(num)),'Color','red')
    % hold off
    % pause(1)
end
%% sif to tif
clear;clc
path='E:\data\20250113\';
files = dir(fullfile(path, '*.sif'));
% 获取该文件下的每个数据文件
datanames = {files(~[files.isdir]).name}';
% 获取数据的数量
datanumber=size(datanames,1);
mkdir([path,'Video_tif\']);
for num=3:3
    dataname_mat=cell2mat(datanames(num,:));
    % 设置每个文件的
    data_save_name=[path,'Video_tif\',dataname_mat(1:end-4),'.tif'];
    % 设置文件读取的路径
    pathname=strcat(path,dataname_mat);
    % 读取该文件
    seq=readsif(pathname);
    framesize=size(seq.imageData,3);
    startframe=1;
    % 确定开始采数据的位置
    for i=1:framesize
        intmean=mean(mean(seq.imageData(:,:,i)));
        if intmean>500
            startframe=i;
            break
        end
    end
    disp(['现在是 ' dataname_mat(1:end-4) '...'])
    %dataedge=128;
    seq.imageData=seq.imageData;
    %data=double(seq.imageData(:,:,startframe:end));
    data=double(seq.imageData(:,:,startframe:end));
    imagesc(data(:,:,1))
    % imagesc(sum(data,3))
    % colormap gray
    % axis square
    % set(gca,'looseInset',[0 0 0 0])
    % set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
    %red=double(seq.imageData(dataedge:255+dataedge,257:512,startframe:end));
    %green=double(seq.imageData(dataedge:255+dataedge,1:256,startframe:end));
    %red=green;
    red=data;
    % 指定输出 TIF 文件名
    outputFileName =data_save_name;
    numFrames=size(red,3);   
    % 写入剩余的帧
    for k = 1:numFrames
        fprintf('Frame %d/%d\n', k,numFrames);
        imwrite(uint16(red(:,:,k)), outputFileName, 'tif', 'WriteMode', 'append');
    end
    disp(['TIF 视频已保存为 ', outputFileName]);

end


%% trackit处理后的数据转为寿命
clear;clc
path='E:\data\20240704\Cy3Cy5\' ;
files = dir(fullfile(path, '*.mat'));
% 获取该文件下的每个数据文件
datanames = {files(~[files.isdir]).name}';
% 获取数据的数量
datanumber=size(datanames,1);
extime=0.1;
for num=1:datanumber
    dataname_mat=cell2mat(datanames(num,:));
    % 设置每个文件的
    data_save_name=[path,dataname_mat(1:end-15),'_surival time.xlsx'];
    % 设置文件读取的路径
    pathname=strcat(path,dataname_mat);
    load(pathname);
    time=zeros(size(trackedPar,2),1);
    for i=1:size(trackedPar,2)
        time(i)=size(trackedPar(i).xy,1)*extime;
    end
    figure
    histogram(time)
    writematrix(time,data_save_name);
end
%%
clear;clc
pathname="E:\data\20240923\Cy5H\Cy5H_638nm_1000frames_0p2s_4mW_fov1.sif";

seq=readsif(pathname);
framesize=size(seq.imageData,3);
startframe=1;
% 确定开始采数据的位置
for i=1:framesize
    intmean=mean(mean(seq.imageData(:,:,i)));
    if intmean>500
        startframe=i;
        break
    end
end
dataedge=130;
seq.imageData=seq.imageData.*5.75./300;
red=double(seq.imageData(dataedge:255+dataedge,257:512,startframe:end));
green=double(seq.imageData(dataedge:255+dataedge,1:256,startframe:end));
%%
figure
imagesc(red(:,:,1))
set(gcf,'Position',[300,200,420,420]);
axis square;
colormap("gray");
set(gca,'looseInset',[0 0 0 0])
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
crange = clim(gca);
figure
imagesc(red(:,:,100/seq.exposureTime))
set(gcf,'Position',[300,200,420,420]);
axis square;
colormap("gray");
set(gca,'looseInset',[0 0 0 0])
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
set(gca, 'CLim', crange);
