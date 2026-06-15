function loc_filter_center=loc_db(data,length)
% 识别第一帧时，部分点没有被识别，尝试选取前100帧来进行DBSCAN，来选取点
loc_all=[];
framesize=size(data,3);

% 确定开始采数据的位置
for i=1:framesize
    intmean=mean(mean(data(:,:,i)));
    if intmean>500*5.75/300
        startframe=i;
        break
    end
end

for i=1:length
    loc=loc_molecule_DF(data(:,:,startframe+i-1),1);
    loc=[loc i.*ones(size(loc,1),1)];
    loc_all=[loc_all;loc];
end

% figure;
% imagesc(data(:,:,startframe+20-1));colormap('gray');axis square;colorbar;
% % 圈出每个数据点
% hold on
% viscircles(loc_all(loc_all(:,3)==20,1:2),3.*ones(size(loc_all(loc_all(:,3)==20,:),1),1),'color','#4F94CD','LineWidth',1.5,'EnhanceVisibility',false);
% %plot(fr_D(:,2),fr_D(:,1),'.','Color','red')
% hold off

% 
% figure
% imagesc(data(:,:,startframe));colormap('gray');axis square;colorbar;
% hold on
% plot(loc_all(:,1),loc_all(:,2),'.','Color','red')
% hold off

% 进行DBSCAN分析
clear loc_all_db

db=DBSCAN(loc_all(:,1:2),1,2);

loc_all_db=[loc_all db];
loc_all_db=sortrows(loc_all_db,4);

on_time=zeros(1);
loc_temp=zeros(1,2);
loc_filter_first=zeros(1,2);
loc_filter_last=zeros(1,2);
loc_filter_center=zeros(1,2);
j=1;k=1;
% 在排序后，数据根据点的顺序来排序
for i=1:(size(loc_all-db,1)-1)
    if loc_all_db(i,4)==loc_all_db(i+1,4)
        loc_temp(j,1)=loc_all_db(i,1);
        loc_temp(j,2)=loc_all_db(i,2);
        loc_temp(j,3)=loc_all_db(i,3);
        j=j+1;
    else
        loc_temp(j,1)=loc_all_db(i,1);
        loc_temp(j,2)=loc_all_db(i,2);
        loc_temp(j,3)=loc_all_db(i,3);
        loc_filter_first(k,1)=loc_temp(1,1);
        loc_filter_first(k,2)=loc_temp(1,2);
        loc_filter_last(k,1)=loc_temp(j,1);
        loc_filter_last(k,2)=loc_temp(j,2);
        loc_filter_center(k,1)=mean(loc_temp(:,1));
        loc_filter_center(k,2)=mean(loc_temp(:,2));
        on_time(k)=loc_temp(j,3);
        k=k+1;
        j=1;loc_temp=zeros(1,2);
    end
end

loc_filter_center=filter_duplicate_points(loc_filter_center,3);
loc=round(loc_filter_center);
filter=zeros(size(loc,1),1);
for i=1:size(loc,1)
    if data(loc(i,2),loc(i,1),1)<2*mean(data(:,:,end),"all")&&data(loc(i,2),loc(i,1),2)<1.5*mean(data(:,:,end),"all")
        filter(i,:)=1;
    end
end
loc_filter_center(filter==1,:)=[];

% figure;
% imagesc(sum(data(:,:,1:100),3));colormap('gray');axis square;colorbar;
% % 圈出每个数据点
% hold on
% viscircles(loc_filter_center,3.*ones(size(loc_filter_center,1),1),'color','#4F94CD','LineWidth',1.5,'EnhanceVisibility',false);
% %plot(fr_D(:,2),fr_D(:,1),'.','Color','red')
% hold off