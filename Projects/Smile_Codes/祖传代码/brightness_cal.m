clear;
clc;
pathway = '/Volumes/Lu_yk_ExFat/20250224_RNO2/';
files = dir(fullfile(pathway, '*.sif*'));
datanames = {files(~[files.isdir]).name}';
% 获取数据的数量
datanumber=size(datanames,1);
brightness = zeros(1,datanumber);
for x=1:datanumber
    % 转变路径格式为char
    dataname_mat=cell2mat(datanames(x,:));
    % 设置每个文件的
    data_save_name=[path,'Figure/',dataname_mat(1:end-4) '.xlsx'];
    % 设置文件读取的路径
    pathname=strcat(pathway,dataname_mat);
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
    
    data = seq.imageData(128:383,128:383,startframe:end) .* 5.75 ./ seq.gainDAC;

    clear loc_total;
    loc_total = loc_molecule(data(:,:,1:50),seq.exposureTime);

    
    for i = 1:size(loc_total,1)
        brightness(i,x) = sum(data(loc_total(i,1)-2:loc_total(i,1)+2,loc_total(i,2)-2:loc_total(i,2)+2,loc_total(i,3)),'all')-mean(sum(sum(data(loc_total(i,1)-2:loc_total(i,1)+2,loc_total(i,2)-2:loc_total(i,2)+2,end-20:end))));
    end
end
