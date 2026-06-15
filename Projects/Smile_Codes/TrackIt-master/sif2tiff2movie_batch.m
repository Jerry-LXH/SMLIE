% 本代码作用是转换某文件夹下的所有或部分 .sif文件，转换为 .tif 文件并输出视频
% 本代码分为3节：
%   第1节：读取某路径下所有sif文件的文件名
%   第2节：查看256矩阵取值的区域
%   第3节：将sif文件转换为tiff文件、avi文件
% 注：256矩阵取值区域存储在 ...\256range\256range.csv 文件中，文件第1行第1列的元素代表行起始值，文件第1行第2列的元素列起始值

clear
clc
warning('off');

%%输入路径%%%
path='/Volumes/Lu_yk_mac/Data/SMLM/未命名文件夹/';%必须以 \ 结尾
%%输入路径%%%

Files=dir(strcat(path,'*.sif'));

% Files=dir(strcat('*.sif'));
LengthFiles=length(Files);
%% 查看256 * 256矩阵取值区域
% 鉴于单个染料单个条件下的数据是按照日期分类的，意味着这个文件夹下所有的数据光斑位置基本相同，故根据某一个.sif文件确定光斑位置即可

%%% 设置取值区域 %%%
row_start = 127; % 256矩阵第一行对应的512矩阵的某一行。即图像的y轴，默认：145。512矩阵中取值区域为：row_start : row_start+255
col_start = 127; % 256矩阵第一列对应的512矩阵的某一列。即图像的x轴，默认：150。512矩阵中取值区域为：col_start : col_start+255
%%% 设置取值区域 %%%

close all

for i=1%1:LengthFiles
    name=Files(i).name;
    path_name=strcat(path,name);
    sif_file = readsif(path_name);
    sif_data = sif_file.imageData;
end

% 显示前100帧里最亮的一帧，并将256矩阵取值范围显示在上面
[~,max_frame] = max(sum(sum(sif_data(:,:,1:10))));
imagesc(sif_data(:,:,max_frame))
colormap gray
set(gcf,'position',[400 50 600 500],'Color','w');
hold on

% 将256矩阵取值范围显示在上面
rectangle('Position',[col_start row_start 256 256],'EdgeColor','g')% rectangle参数含义：[矩形左上角点的x坐标(这是根据实践来的，matlab中的示例是左下角)，矩形左上角点的y坐标，x轴长度，y轴长度]

% 获取文件名
[~,file_name,~] = fileparts(name);

% 储存数据
mkdir(path,'256range')
saveas(gcf,[path,'256range\',file_name,'.png'])

% 存储256矩阵取值区域
range = [row_start,col_start];
writematrix(range,[path,'256range\','256range.csv'])

% close all
%% convert
% 若需转换所有文件，i的取值范围为 1:LengthFiles
% 若仅转换部分文件，i的取值范围为待转换文件在 Files 变量中的编号

%%% 输入参数 %%%
      tiff_file = 1; % 如果需要输出 .tif 文件，此项为1；若不需要，为0
       avi_file = 0; % 如果需要输出 .avi 文件，此项为1；若不需要，为0
convert512to256 = 1; % 如果需要将512矩阵转换为256矩阵，此项为1；若不需要，为0
%%% 输入参数 %%%

clear i name path_name sif_file sif_data max_frame file_name
close all

for i=1:LengthFiles
    try
    name=Files(i).name;
    path_name=strcat(path,name);
    s=readsif(path_name);
    
    %     names1=s.fileName;
    [~,file_name,~] = fileparts(name);

    si=s.imageData;
    DAC=s.gainDAC;
    [r,ss,t]=size(si);
    time=s.accumulateCycleTime;
    
    if convert512to256 == 1
        if r==512
            si=si(row_start:row_start+255,col_start:col_start+255,:);%取256x256分析
        end
    end
    
%     strs1=strfind(names1,'\');
%     simpNames1=names1((strs1(end)+1):(end-4));
    
    % tiff
    if tiff_file == 1
        
        for j=1:t
            %st=255*si(:,:,i)./max(max(si(:,:,i)));
            st=si(:,:,j);
            imwrite(uint16(st),[path,file_name,'_Tifstack_Uint16.tif'],'WriteMode','append')
        end
        clear si
    end

    %   %若无需输出moive，注释此行-最后的end前一行即可
    if avi_file == 1
        %movie
        s11=s.imageData;
        [r,ss,t]=size(s11);
        s1=s11.*5.75./s.gainDAC;
        if convert512to256 == 1
            if r==512
                s1=s1(row_start:row_start+255,col_start:col_start+255,:);%取256*256矩阵
            end
        end
        
        v=VideoWriter([path,file_name,'.avi','MPEG-4']);%储存到video中。采用'MPEG-4'是因为看到mathwork上提到'MPEG-4'可适用于mac系统
        open(v);
        for k=1:t
            fig=imagesc(s1(:,:,k));
            axis equal
            if size(s1,1) == 512
                xlim([1,512])
                ylim([1,512])
            else if size(s1,1) == 256
                    xlim([1,256])
                    ylim([1,256])
                end
            end
            set(gca,'FontSize',15)
            colormap gray
            title(['time=',num2str(k*time),'s'])
            pause(0.1)
            frame=getframe(gcf);
            writeVideo(v,frame);
        end
        
        set(gcf,'CloseRequestFcn','closereq');
        close(v)
        
        close all
    end
    
    catch ME
        disp(['Failed: 【',num2str(i),'】 ', path_name])
        disp(ME.message)
        disp('——————————————')
        continue
    end
end

% 分析完毕后，弹窗提醒
message_box = msgbox('Sif file convert completed');
set(message_box,'Units','normalized','OuterPosition',[0.3 0.4 0.3 0.3]);
set(findobj(message_box,'Type','Text'), 'FontSize', 42 );

