clear;clc
tic
%load data
% path='E:\data\ZWR\20230718\29-15Er\';
name='hcy_ph2p24_638nm_20220331_0p1s_2000f_10mw_fov2.sif';%%%%%%%%%
seq=readsif(name);

% filename='E:\data\dye\data_dye.xlsx';
%%
clc
imgdata=double(seq.imageData);
% figure;imagesc(imgdata);axis square
% 数据转置

imgdata=imgdata(128:383,128:383,:);
imgdata=(imgdata.*5.75)./300;
% ima1_1=imgdata';
% 二值化
locdata=[];
% outputString = sprintf('Parameter: ');
fprintf('%s', 'Parameter: ');
for i=1:1000
    fprintf('%d', i);
    S=rescale(imgdata(:,:,i));
    % imagesc(S);
    bg=imopen(S,strel('disk',4));
    S=S-bg;
    % find the threshold
    level=graythresh(S);
    bw=imbinarize(S,level);
    [label,object]=bwlabel(bw,8);
    % bw = bwmorph(bw,"clean");
    % RGB_label=label2rgb(label);imshow(RGB_label);
    % 将regionprops函数得到的数据存入到stats
    stats=regionprops('table',label,imgdata(:,:,i),'Centroid',...
        'Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');
    stats(stats.Area<7,:)=[];
    stats(stats.MaxIntensity<1.2*mean(stats.MaxIntensity),:)=[];
    center=round(stats.Centroid);
    locdata_temp=[center stats.Area i*ones(size(stats,1),1) (1:size(stats,1))'];
    locdata=[locdata;locdata_temp];
    if i ~= 1000                                                                                  
        fprintf(repmat('\b', 1, numel(num2str(i)))); % 退格符的数量与输出字符数量相同
    end 
end
fprintf('\n');
%%
scatter(locdata(:,1),locdata(:,2),Marker=".")
DBSCAN()
%%
sumdata=sum(imgdata(:,:,1:500),3);
bg=imopen(sumdata,strel('disk',4));
sumdata=sumdata-bg;
[X,Y]=meshgrid(1:256,1:256);
surf(X,Y,sumdata)
colorbar
ss=rescale(sumdata);
level=graythresh(ss);
bw=imbinarize(ss,2*level);
[label,object]=bwlabel(bw,4);
% bw = bwmorph(bw,"clean");
% RGB_label=label2rgb(label);imshow(RGB_label);
% 将regionprops函数得到的数据存入到stats
statsum=regionprops('table',label,sumdata,'Centroid',...
    'Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');



