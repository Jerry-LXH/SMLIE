%% load data
clear
file_name ="/Users/lu_yk/Desktop/Data/SMLM/20250210_R-4/R-4_0p1s_532nm_5000frs_2mw_fov1_sf1.sif";
s_frame =1000;  % start frame

[data1, ex_time, gainDAC] = pro_data(file_name, s_frame);    %ex_time: expousre time
for i=1:100
    intmean=mean(mean(data1(:,:,i)));
    if intmean>500/300*5.75
        startframe=i;
        break
    end
end
raw_start=120;
line_start=120;
data=data1(raw_start:raw_start+255,line_start:line_start+255,startframe:end);
imagesc(data1(:,:,startframe))
hold on
rectangle('Position',[line_start,raw_start,255,255],'LineWidth',1,'Visible','on','EdgeColor','red');
hold off
figure;imagesc(data(:,:,1));
%% snapshoot for checkframe
figure('visible','on')
imagesc(data1(:,:,2000))
axis square
% title('control 638ex');
colormap('gray')
%colorbar
set(gca,'looseInset',[0 0 0 0])
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
% set(gca, 'FontName', 'Arial');
% set(gca,'YDir','normal');
set(gcf,'Position',[300,200,420,420]);
%set(gcf,'position',[300 100 600 500]);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
% saveas(gcf, [file_name, '\control 638ex', 'tif']);
%% set parameter for loc_molecule
clear loc_total0
check_frame = 3;   % located frame of image
loc_total0 = loc_molecule(data(:,:,check_frame), ex_time);
check_localization(data(:,:,check_frame), loc_total0); 
% check_localization(data(:,:,check_frame), loc_red(loc_red(:,3)==check_frame,:)); 
axis square
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf,'-property','FontSize'), 'FontSize', 20);
title(['frame', ' ', num2str(check_frame)]);
colorbar
%% localing frame by frame
clear loc_total_o

loc_total_o = loc_molecule(data(:,:,:), ex_time);
%% set parameter for DBSCAN
check_frame = 500;   % located frame of image

check_localization(data(:,:,check_frame), loc_total_o(loc_total_o(:,3)==check_frame,:)); 
% check_localization(data(:,:,check_frame), loc_red(loc_red(:,3)==check_frame,:)); 
axis square
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf,'-property','FontSize'), 'FontSize', 20);
title(['frame', ' ', num2str(check_frame)]);
colorbar

figure('visible','on')
scatter(loc_total_o(:,2),loc_total_o(:,1),'.');
set(gcf,'position',[300 100 600 500]);
%% DBSCAN for localing
clear step r minn E on_time
t=1;E=1;step=500;r=1;minn=2;loc_total=zeros(1,3);

while t<=E
    clear db 
    i=size(loc_total,1);k=size(loc_total,1);
    for i=i:size(loc_total_o,1)
        if loc_total_o(i,3)>=((t-1)*step+1)&&loc_total_o(i,3)<=(t*step)
            loc_total(k,1:3)=loc_total_o(i,:);
            k=k+1;
        end
    end
%     disp('Load loc_total complete')
    db=dbscan(loc_total(:,1:2),r,minn);
    for i=1:size(loc_total,1)
        loc_total(i,3)=loc_total_o(i,3);
        loc_total(i,4)=db(i);
    end
%     disp('dbscan complete')
    loc_temp=zeros(1,2);
    loc_filter_first=zeros(1,2);
    loc_filter_last=zeros(1,2);
    loc_filter_center=zeros(1,2);
    loc_total=sortrows(loc_total,4);
    j=1;k=1;
    for i=1:(size(loc_total,1)-1)
        if loc_total(i,4)==loc_total(i+1,4)
            loc_temp(j,1)=loc_total(i,1);
            loc_temp(j,2)=loc_total(i,2);
            loc_temp(j,3)=loc_total(i,3);
            j=j+1;
        else
            loc_temp(j,1)=loc_total(i,1);
            loc_temp(j,2)=loc_total(i,2);
            loc_temp(j,3)=loc_total(i,3);
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
    loc_filter_last(1,:)=[];
    loc_filter_first(1,:)=[];
    loc_filter_center(1,:)=[];
    while loc_total(1,4)==-1
            loc_total(1,:)=[];
    end
   c{t}=loc_total; 
   t=t+1;
     on_time(1)=[];
     on_time=on_time.*ex_time;
     on_time=on_time';
     [loc_filter_single,inx,p]=unique(loc_total(:,4), 'rows');
     for i=1:size(inx,1)
         loc_filter_single(i,1)=loc_total(inx(i),1);
         loc_filter_single(i,2)=loc_total(inx(i),2);
         loc_filter_single(i,3)=loc_total(inx(i),3);
         loc_filter_single(i,4)=loc_total(inx(i),4);
     end
     loc_filter_single=sortrows(loc_filter_single,4);
%     scatter(loc_filter_first(:,2),loc_filter_first(:,1),'o');
%     set(gcf,'position',[300 100 600 500]);
%     title(['times=',t,'s'])
%     hold on
end

for i=1:size(c,2)
    scatter(c{i}(:,2),c{i}(:,1),'.')
    axis equal
    xlim([1,256])
    ylim([1,256])
    hold on
    viscircles([loc_filter_first(:,2),loc_filter_first(:,1)],3,'color','red','LineWidth',1.5,'EnhanceVisibility',false);
    pause(1)
end
    
% figure
% scatter(loc_total(:,2),loc_total(:,1),'.');
%     hold on
%     scatter(loc_filter_first(:,2),loc_filter_first(:,1),'o');
%     set(gcf,'position',[300 100 600 500]);
%     title(['times=',t,'s'])
%     check_frame = 1;   % located frame of image
%     check_localization(data(:,:,check_frame), loc_filter_first); 
%     % check_localization(data(:,:,check_frame), loc_red(loc_red(:,3)==check_frame,:)); 
%     axis square
%     set(gca, 'FontName', 'Times New Roman');
%     set(findall(gcf,'-property','FontSize'), 'FontSize', 20);
%     title(['frame', ' ', num2str(check_frame)]);
%     colorbar
%     %loc_total=loc_filter2;

data_temp=zeros(256,256);
for i=1:256
    for j=1:256
        for k=1:10
            data_temp(i,j)=data_temp(i,j)+data(i,j,k);
        end
    end
end
figure('name','intensity distribution','visible','on')
imagesc(data(:,:,1000));% a surface diagram of intensity distribution\
idlmapnames=transpose(h5read('idlcolormaps.h5','/name'));
mymapindex=6;%%%%%%%%%%%%%%%%
mycolormap=idlmapnames{mymapindex};
disp([string(mymapindex), mycolormap]);
idlrgbtables=h5read('idlcolormaps.h5', '/rgbt');
colormap(double(idlrgbtables(:,:,mymapindex))/255.);
colorbar off 
xlabel('x','FontSize',15);ylabel('y','FontSize',15);zlabel('intensity','FontSize',15);
%xlim([1,r]);ylim([1,ss])
shading interp
axis equal

%caxis([60 1300])
hold on
viscircles([loc_filter_first(:,2),loc_filter_first(:,1)],3,'color','red','LineWidth',1.5,'EnhanceVisibility',false);
%% extract fluorescence trace
% loc_single = unique(loc_total_o(:,1:2), 'rows');  % remove duplication of localizations
raw_series = time_trace(loc_filter_first(:,1:2), data(:,:,:));
[~, frames] = size(raw_series);
% loc_red = unique(loc_red(:,1:2), 'rows');  % remove duplication of localizations
% raw_series = time_trace(loc_red, data,2);
series = raw_series - median(raw_series(:,frames-10:frames), 2);

% cut_series = series(:,1:2000);
%% extract fluorescence trace(for drifting)
raw_series=[];
for j=1:size(loc_filter_first,1)
    loc_temp=[];
    jump=1;
    dis_max=6;
    a=[];
    [loc_max,a] = max(data1(raw_start+loc_filter_first(j,1)-1-3-jump:raw_start+loc_filter_first(j,1)-1+3+jump,line_start+loc_filter_first(j,2)-1-3-jump:line_start+loc_filter_first(j,2)-1+3+jump,1),[],'all');
    loc_temp=[raw_start+loc_filter_first(j,1)-1-3-jump-1+floor(a/((3+jump)*2+1))+1,line_start+loc_filter_first(j,2)-1-3-jump-1+a-floor(a/((3+jump)*2+1))*((3+jump)*2+1)];
    for i=1:size(data,3)-1000
        raw_series(j,i)=sum(data1(loc_temp(1)-3:loc_temp(1)+3,loc_temp(2)-3:loc_temp(2)+3,i+1000),'all');
        [loc_max,a] = max(data1(raw_start+loc_filter_first(j,1)-1-3-jump:raw_start+loc_filter_first(j,1)-1+3+jump,line_start+loc_filter_first(j,2)-1-3-jump:line_start+loc_filter_first(j,2)-1+3+jump,i+1),[],'all');
        raw_temp=raw_start+loc_filter_first(j,1)-1-3-jump-1+a-floor(a/((3+jump)*2+1))*((3+jump)*2+1);
        line_temp=line_start+loc_filter_first(j,2)-1-3-jump-1+floor(a/((3+jump)*2+1))+1;
        loc_mean=mean(data1(raw_temp-3:raw_temp+3,line_temp-3:line_temp+3,i+1),'all');
        if loc_max>2*loc_mean&&norm([raw_temp,line_temp]-loc_temp)<dis_max&&norm([raw_temp,line_temp]-[raw_start+loc_filter_first(j,1)-1,line_start+loc_filter_first(j,2)-1])<dis_max
            loc_temp=[raw_temp,line_temp];
        end
        raw_series(j,i+1)=sum(data1(loc_temp(1)-3:loc_temp(1)+3,loc_temp(2)-3:loc_temp(2)+3,i+1000),'all');
    end
    j
end
[~, frames] = size(raw_series);
% loc_red = unique(loc_red(:,1:2), 'rows');  % remove duplication of localizations
% raw_series = time_trace(loc_red, data,2);
series = raw_series - median(raw_series(:,frames-10:frames), 2);
%% plot picture(without fitting)
on_time=[];
off_time=[];
intesity=[];
for num=1:25
    %if blink_counts(num,1)>10
%for num=1:200
    % clear all
    % load('matlab.mat');
    intis1=series(num,1:end);
    intis1=intis1';
    %intis1=smooth(intis1,1);
    %intis1=medfilt1(intis1,5);
    %Set the time step
    timestep=ex_time;
    t=length(intis1);
    figure('visible','on')
    plot((1:t)*timestep,intis1)
    xlabel('time(s)')
    ylabel('Intensity')
    set(gca,'FontSize',15)
    %end
end
%% manual filtter
chosen_n=[];
for i=1:250
    if gcf().Number~=1
        chosen_n=horzcat(chosen_n,gcf().Number);
        close gcf
    elseif i==1
        chosen_n=horzcat(chosen_n,gcf().Number);
        close gcf
    end
end
%series_new=zeros(size(chosen_n,2),size(series,2));
for i=1:size(chosen_n,2)
    series_new(i,:)=series(chosen_n(i),:);
end
%series_new1=series_new';
%%
filter=find(blink_counts>20);
for i=size(filter,1):-1:1
    series(filter(i),:)=[];
end
%% manual analysis_blinking
%series=series_new;
split={};
dutycycle=[];
Intensity=[];
blink_counts=zeros(size(series,1),1);
for i=1:size(series,1)
    threshold=200;
    split_temp={};
    blink_tem=0;
    temp_on=[];
    temp_off=[];
    Intensity_tem=[];
    for j=1:size(series,2)
        if series(i,j)<=threshold
            if size(temp_on,1)~=0
                split_temp=horzcat(split_temp,temp_on);
                temp_off=horzcat(temp_off,series(i,j));
                temp_on=[];
            else
                temp_off=horzcat(temp_off,series(i,j));
            end
            % if j==size(series,2)
            %    split_temp=horzcat(split_temp,temp_off);
            % end
        else
            if size(temp_off,1)~=0
                split_temp=horzcat(split_temp,temp_off);
                temp_on=horzcat(temp_on,series(i,j));
                temp_off=[];
            else
                temp_on=horzcat(temp_on,series(i,j));
            end
            % if j==size(series,2)
            %    split_temp=horzcat(split_temp,temp_on);
            % end
        end
        on_time_temp=0;
        off_time_temp=0;
    end
    for k=1:size(split_temp,2)
        if split_temp{k}(1)<=threshold
            off_time_temp=off_time_temp+size(split_temp{k},2);
        else
            blink_tem=blink_tem+1;
            on_time_temp=on_time_temp+size(split_temp{k},2);
            Intensity_tem=horzcat(Intensity_tem,split_temp{k});
        end
    end
    blink_counts(i,1)=blink_tem;
    Intensity(i,1)=mean(Intensity_tem(:));
    dutycycle=horzcat(dutycycle,(on_time_temp/(on_time_temp+off_time_temp)));
    split=horzcat(split,split_temp);
end
dutycycle=dutycycle';
% 计算每帧亮度
Intensity_all=[];
for i=1:size(split,2)
    if split{i}(1)>threshold
        Intensity_all=horzcat(Intensity_all,split{i});
    end
end
on_time=[];
off_time=[];
for i=1:size(split,2)
    if split{i}(1)<=threshold
        off_time=horzcat(off_time,size(split{i},2)*ex_time);
    else
        on_time=horzcat(on_time,size(split{i},2)*ex_time);
    end
end
mean(on_time)
on_time=on_time';
off_time=off_time';
histogram(dutycycle);title('Dutycycle')
figure;histogram(on_time);title('On time')
figure;histogram(off_time);title("Off time")
figure;histogram(blink_counts);title("Blink")
figure;histogram(Intensity);title("Intensity")
figure;histogram(Intensity_all);title("Intensity All")
%% write data
mkdir([file_name{1}(1:end-9) '\']);
writematrix(series',[file_name{1}(1:end-9) '\' 'Trace' '.xlsx'])
%writematrix(dutycycle,[file_name{1}(1:end-4) '\' 'Dutycycle' '.xlsx'])
writematrix(on_time,[file_name{1}(1:end-9) '\' 'Ontime' '.xlsx'])
writematrix(off_time,[file_name{1}(1:end-9) '\' 'Offtime' '.xlsx'])
writematrix(blink_counts,[file_name{1}(1:end-9) '\' 'Blink counts' '.xlsx'])
%writematrix(Intensity_all',[file_name{1}(1:end-4) '\' 'Intensity' '.xlsx'])
%% plot trace(filter)
for num=1:size(series_new1,2)
    % clear all
    % load('matlab.mat');
    intis1=series_new1(:,num);
    %intis1=intis1';
    %intis1=smooth(intis1,1);
    %intis1=medfilt1(intis1,5);
    %Set the time step
    timestep=0.1;
    t=length(intis1);
    figure('visible','on')
    plot((1:t)*timestep,intis1)
    xlabel('time(s)')
    ylabel('Intensity')
    set(gca,'FontSize',15)   
end
%% plot trace(filter & with fitting)
on_time1=[];
off_time=[];
intesity=[];
total_photons=[];
STD_INTENSITY=[];
STD2_INTENSITY=[];

clear intensity_ontime on_time
series_new=series;
for num=1:size(series_new,1)
%for num=1:157
    % clear all
    % load('matlab.mat');
    clear Hint HMMint 
    intis1=series_new(num,:);
    intis1=intis1';
    %intis1=smooth(intis1,1);
    %intis1=medfilt1(intis1,5);
    %Set the time step
    timestep=0.1;

    t=length(intis1);%The length of observation sequence

    ChooseMethod=1;
    % Method 1: Using hierarchical clustering
    if ChooseMethod==1
        X=intis1;
        %Plot the observed data curve
        % figure('visible','on')
        % I=plot((1:t)*timestep,intis1)
        % xlabel('time(s)')
        % ylabel('Intensity')
        % set(gca,'FontSize',15)

        Y=pdist(X,'euclidean');%Compute the pairwise distances
        Z=linkage(Y,'average');%Generate clustering hierarchical tree according to distance information
        %Calculate cophenetic correlation coefficient for the hierarchical cluster tree, a larger value indicates that the tree fits the distance well
        C=cophenet(Z,Y)

        %%%%%%Natural Divisions
        %Set cutoff threshold according from standard deviation and cophenetic correlation coefficient
        I1=find(intis1~=0);
        cutoff=20*std(intis1(I1((end-50):end)));%If the number of categories
        %is small,try increasing the cutoff,otherwise reducing the cutoffD
        %clustering result
        seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
        seq=seq';

        %Determine emission value
        F1=length(unique(seq));
        intensity=zeros(1,F1);
        for i=1:F1
            intensity(i)=mean(intis1(find(seq==i)));
        end

        Ctrace=zeros(1,length(intis1));
        for i=1:length(intensity)
            Ctrace(find(seq==i))=intensity(i);
        end
        R=corrcoef(Ctrace',intis1);
        i=0.1;
        R(1,2)

        while R(1,2)<0.88
            cutoff=(20-i)*std(intis1(end-50:end));
            i=i+0.1;
            seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
            seq=seq';
            F1=length(unique(seq));
            intensity=zeros(1,F1);
            for j=1:F1
                intensity(j)=mean(intis1(find(seq==j)));
            end
            for j=1:F1
                Ctrace(find(seq==j))=intensity(j);
            end
            R=corrcoef(Ctrace,intis1);
        end

        %%Method 2: Finding the peak position of the histogram of the observation series
    elseif ChooseMethod==2
        %Set peak threshold, number of bins, minimum peak separation, method to determine emission, method to determine dividing line

        %Set Minimum peak height
        PeakThreshold=1;
        %Set number of bins
        Nhist=round(t/20);
        %Set Minimum peak separation
        MinPeakdistance=round(Nhist/20);
        %The emission is determined by median(method=0),peak(method=1), or mean(method=2) value
        method=0;
        %The deviding line is determined by peak width(bmethod=0) or mean value of two neighbouring peaks(bmethod=1)
        bmethod=1;

        %Seq:Emission sequence
        %location: Emission intensity
        [seq,intensity]=HMMEmitTrc(intis1,timestep,PeakThreshold,Nhist,MinPeakdistance,method,bmethod);
    end

    %%%%%% Baum-welch algorithm is used to estimate the transfer probability, emission probability 
    lloc=length(intensity);
    TRANS_GUESS = eye(lloc)*1+rand(lloc,lloc)*0.0000001;
    EMIS_GUESS = eye(lloc)*0.99+rand(lloc,lloc)*0.01;
    [TRANS_EST2, EMIS_EST2] = hmmtrain(seq, TRANS_GUESS, EMIS_GUESS);

    %%%%%% The state sequence is estimated by Viterbi algorithm
    likelystates = hmmviterbi(seq, TRANS_EST2, EMIS_EST2);
    likest=likelystates;
    likest=sort(likest);
    likest(size(likest,2)+1)=-123.45;
    k=1;j=0;intensity_time=[0,0];
    for i=1:(size(likest,2)-1)
        if likest(i)==likest(i+1)
            k=k+1;
        else
            j=j+1;
            intensity_time(j,1)=k*ex_time;
            intensity_time(j,2)=intensity(j);
            k=1;
        end
    end
    try
    intensity_time=sortrows(intensity_time,2);
    on_time=sum(intensity_time(:,1))-intensity_time(1,1);
    intensity_ontime(num,1)=on_time;
    intensity_time(1,:)=[];
    intensity_time=sortrows(intensity_time,1);
    intensity1=intensity_time(size(intensity_time,1),2);
    intensity_ontime(num,2)=intensity1;
    likest=likelystates ;
    for i=1:length(intensity)
        likest(find(likelystates==i))=intensity(i);
    end
    split={};
    split_data={};
    temp=[];
    temp_data=[];
    count=1;
    count2=1;
    base_state=likelystates(size(likelystates,2));
    for i=1:size(likelystates,2)-1
        if likelystates(1,i)==likelystates(1,i+1)
            temp(count)=likelystates(1,i);
            temp_data(count)=intis1(i,1);
            count=count+1;
        else
            temp(count)=likelystates(1,i);
            temp_data(count)=intis1(i,1);
            split{count2}=temp;
            split_data{count2}=temp_data;
            count2=count2+1;
            temp=[];
            temp_data=[];
            count=1;
        end
    end
    total_photons(num)=0;
    blink_counts(num,1)=0;
    for i=1:size(split,2)
        if split{i}(1)==base_state
            off_time=horzcat(off_time,[size(split{i},2)]);
        else
            blink_counts(num,1)=blink_counts(num,1)+1;
            on_time1=horzcat(on_time1,[size(split{i},2)]);
            intesity=horzcat(intesity,[sum(split_data{i})/size(split{i},2)]);
            total_photons(num)=total_photons(num)+sum(split_data{i});
            STD2_INTENSITY=horzcat(STD2_INTENSITY,[(std(split_data{i}))^2/(sum(split_data{i})/size(split{i},2))]);
            STD_INTENSITY=horzcat(STD_INTENSITY,[std(split_data{i})/(sum(split_data{i})/size(split{i},2))]);
        end
    end
    hold on
    plot((1:t)*timestep,likest);
    disp(['The number of state is ',num2str(length(intensity))])
    disp(['The probability transition matrix is',newline])
    disp(num2str(TRANS_EST2))
    end
    %saveas(I, ['E:\My work\data\trace\TEMP\', 'trace', int2str(num), '-intensity-',num2str(intensity1),'-on_time-',num2str(on_time),'.fig']);
end
on_time1=on_time1'.*ex_time;
off_time=off_time'.*ex_time;
intesity=intesity';
total_photons=total_photons';
STD_INTENSITY=STD_INTENSITY';
STD2_INTENSITY=STD2_INTENSITY';

%% data for plotting;
clear dataX dataY
dataY=series(117,:)';
dataX=(1:size(series,2)).*ex_time;
%dataY=dataY';
dataX=dataX';
%% check trace
clear Hint HMMint num state_num
good_num=1;
for num=1:size(series,1)
    timestep=0.1;
    t=size(series,2);%The length of observation sequence
    X=series(num,:)';
    % figure('visible','on')
    % plot((1:t)*timestep,X)
    % xlabel('time(s)')
    % ylabel('Intensity')
    % set(gca,'FontSize',15)
%     X=smooth(X,3);
%     X=medfilt1(X,5);
    %Plot the observed data curve
    % figure('visible','on')
    % plot((1:t)*timestep,X)
    % xlabel('time(s)')
    % ylabel('Intensity')
    % set(gca,'FontSize',15)

    Y=pdist(X,'euclidean');%Compute the pairwise distances
    Z=linkage(Y,'average');%Generate clustering hierarchical tree according to distance information
    %Calculate cophenetic correlation coefficient for the hierarchical cluster tree, a larger value indicates that the tree fits the distance well
    C=cophenet(Z,Y)

    %%%%%%Natural Divisions
    %Set cutoff threshold according from standard deviation and cophenetic correlation coefficient
    I1=find(X~=0);
    cutoff=20*std(X(I1((end-50):end)));%If the number of categories
    %is small,try increasing the cutoff,otherwise reducing the cutoffD
    %clustering result
    seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
    seq=seq';

    %     Determine emission value
        F1=length(unique(seq));
        intensity=zeros(1,F1);
        for i=1:F1
            intensity(i)=mean(X(find(seq==i)));
        end

        Ctrace=zeros(1,length(X));
        for i=1:length(intensity)
            Ctrace(find(seq==i))=intensity(i);
        end
        R=corrcoef(Ctrace',X);
        i=0.1;
        R(1,2);

        while R(1,2)<0.90
            cutoff=(20-i)*std(X(end-50:end));
            i=i+0.1;
            seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
            seq=seq';
            F1=length(unique(seq));
            intensity=zeros(1,F1);
            for j=1:F1
                intensity(j)=mean(X(find(seq==j)));
            end
            for j=1:F1
                Ctrace(find(seq==j))=intensity(j);
            end
            R=corrcoef(Ctrace,X);
        end
     clear lloc
    % Baum-welch algorithm is used to estimate the transfer probability, emission probability 
    lloc=length(intensity);
    TRANS_GUESS = eye(lloc)*1+rand(lloc,lloc)*0.01;
    EMIS_GUESS = eye(lloc)*1+rand(lloc,lloc)*0.01;
    [TRANS_EST2, EMIS_EST2] = hmmtrain(seq, TRANS_GUESS, EMIS_GUESS);

    %%%%%% The state sequence is estimated by Viterbi algorithm
    likelystates = hmmviterbi(seq, TRANS_EST2, EMIS_EST2);

    %Draw analysis results
    hold on
    likest=likelystates ;
    for i=1:length(intensity)
    likest(find(likelystates==i))=intensity(i);
    end 
    % plot((1:t)*timestep,likest);
    % disp(['The number of state is ',num2str(length(intensity))])
    % disp(['The probability transition matrix is',newline])
    % disp(num2str(TRANS_EST2))

    Hint=likelystates;
    for i=length(intensity):(-1):1
    Hint(find(Hint==i))=intensity(i);
    end

    HMMint=Hint';%State sequence

    for i=1:size(HMMint,1)
        HMMint(i,2)=i;
    end

    j=1;k=1;q=0;
    clear HMMint_temp Hint_new
    for i=1:(size(HMMint,1)-1)
        if HMMint(i,1)==HMMint(i+1,1)
            HMMint_temp(j)=X(i);
            j=j+1;
        else
            HMMint_temp(j)=X(i);
            Hint_new(k,1)=std(HMMint_temp);
            Hint_new(k,2)=mean(HMMint_temp);
            Hint_new(k,3)=j;
            q=q+j;
            j=1;
            k=k+1;
            clear HMMint_temp
        end
        if i==size(HMMint,1)-1
            HMMint_temp(j)=X(i);
            Hint_new(k,1)=std(HMMint_temp);
            Hint_new(k,2)=mean(HMMint_temp);
            Hint_new(k,3)=j;
            q=q+j;
        end     
    end
    state_num(num)=length(intensity);
    if length(intensity)<=3 && length(intensity)~=1
        good_num=good_num+1;
    end
end
state_num=state_num';
clear a
k=1;
for i=1:size(state_num,1)
    if state_num(i)>5 | state_num(i)==1 | state_num(i)==0
        a(k)=i;
        k=k+1;
    end
end
%% filter trace
clear loc_filter_single1 
loc_single_filt=[0,0];
loc_filter_single1=loc_filter_first;

for i=1:size(a,2)
    loc_filter_single1(a(i),1)=-1;
end
k=1;
for i=1:size(loc_filter_single1,1)
    if loc_filter_single1(i,1)~=-1
        loc_single_filt(k,1:2)=loc_filter_single1(i,1:2);
        k=k+1;
    end
end
loc_single_filt(1,:)=[];
raw_series = time_trace(loc_single_filt, data(:,:,:),2);
[~, frames] = size(raw_series);
% loc_red = unique(loc_red(:,1:2), 'rows');  % remove duplication of localizations
% raw_series = time_trace(loc_red, data,2);
series = raw_series - median(raw_series(:,frames-10:frames), 2);

% cut_series = series(:,1:2000);
[~, frames] = size(series);
threshold_on = mean(series(:, frames-10:frames), 2) + 10*std(series(:, frames-8:frames), 0, 2);   % threshold for fluorescent state
%% plot picture(with fitting)
on_time1=[];
off_time=[];
intesity=[];
total_photons=[];
STD_INTENSITY=[];
STD2_INTENSITY=[];
mkdir([file_name{1}(1:end-9) '\']);
clear intensity_ontime on_time
for num=1:size(series,1)
%for num=1:size(series,1)
    % clear all
    % load('matlab.mat');
    clear Hint HMMint 
    intis1=series(num,:);
    intis1=intis1';
    %intis1=smooth(intis1,1);
    %intis1=medfilt1(intis1,5);
    %Set the time step
    timestep=0.1;

    t=length(intis1);%The length of observation sequence

    ChooseMethod=1;
    % Method 1: Using hierarchical clustering
    if ChooseMethod==1
        X=intis1;
        %Plot the observed data curve
        fig=figure('visible','off');
        I=plot((1:t)*timestep,intis1);
        xlabel('time(s)')
        ylabel('Intensity')
        set(gca,'FontSize',15)

        Y=pdist(X,'euclidean');%Compute the pairwise distances
        Z=linkage(Y,'average');%Generate clustering hierarchical tree according to distance information
        %Calculate cophenetic correlation coefficient for the hierarchical cluster tree, a larger value indicates that the tree fits the distance well
        C=cophenet(Z,Y)

        %%%%%%Natural Divisions
        %Set cutoff threshold according from standard deviation and cophenetic correlation coefficient
        I1=find(intis1~=0);
        cutoff=20*std(intis1(I1((end-50):end)));%If the number of categories
        %is small,try increasing the cutoff,otherwise reducing the cutoffD
        %clustering result
        seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
        seq=seq';

        %Determine emission value
        F1=length(unique(seq));
        intensity=zeros(1,F1);
        for i=1:F1
            intensity(i)=mean(intis1(find(seq==i)));
        end

        Ctrace=zeros(1,length(intis1));
        for i=1:length(intensity)
            Ctrace(find(seq==i))=intensity(i);
        end
        R=corrcoef(Ctrace',intis1);
        i=0.1;
        R(1,2)

        while R(1,2)<0.88
            cutoff=(20-i)*std(intis1(end-50:end));
            i=i+0.1;
            seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
            seq=seq';
            F1=length(unique(seq));
            intensity=zeros(1,F1);
            for j=1:F1
                intensity(j)=mean(intis1(find(seq==j)));
            end
            for j=1:F1
                Ctrace(find(seq==j))=intensity(j);
            end
            R=corrcoef(Ctrace,intis1);
        end

        %%Method 2: Finding the peak position of the histogram of the observation series
    elseif ChooseMethod==2
        %Set peak threshold, number of bins, minimum peak separation, method to determine emission, method to determine dividing line

        %Set Minimum peak height
        PeakThreshold=1;
        %Set number of bins
        Nhist=round(t/20);
        %Set Minimum peak separation
        MinPeakdistance=round(Nhist/20);
        %The emission is determined by median(method=0),peak(method=1), or mean(method=2) value
        method=0;
        %The deviding line is determined by peak width(bmethod=0) or mean value of two neighbouring peaks(bmethod=1)
        bmethod=1;

        %Seq:Emission sequence
        %location: Emission intensity
        [seq,intensity]=HMMEmitTrc(intis1,timestep,PeakThreshold,Nhist,MinPeakdistance,method,bmethod);
    end

    %%%%%% Baum-welch algorithm is used to estimate the transfer probability, emission probability 
    lloc=length(intensity);
    TRANS_GUESS = eye(lloc)*1+rand(lloc,lloc)*0.0000001;
    EMIS_GUESS = eye(lloc)*0.99+rand(lloc,lloc)*0.01;
    [TRANS_EST2, EMIS_EST2] = hmmtrain(seq, TRANS_GUESS, EMIS_GUESS);

    %%%%%% The state sequence is estimated by Viterbi algorithm
    likelystates = hmmviterbi(seq, TRANS_EST2, EMIS_EST2);
    likest=likelystates;
    likest=sort(likest);
    likest(size(likest,2)+1)=-123.45;
    k=1;j=0;intensity_time=[0,0];
    for i=1:(size(likest,2)-1)
        if likest(i)==likest(i+1)
            k=k+1;
        else
            j=j+1;
            intensity_time(j,1)=k*ex_time;
            intensity_time(j,2)=intensity(j);
            k=1;
        end
    end
    try
    intensity_time=sortrows(intensity_time,2);
    on_time=sum(intensity_time(:,1))-intensity_time(1,1);
    intensity_ontime(num,1)=on_time;
    intensity_time(1,:)=[];
    intensity_time=sortrows(intensity_time,1);
    intensity1=intensity_time(size(intensity_time,1),2);
    intensity_ontime(num,2)=intensity1;
    likest=likelystates ;
    for i=1:length(intensity)
        likest(find(likelystates==i))=intensity(i);
    end
    split={};
    split_data={};
    temp=[];
    temp_data=[];
    count=1;
    count2=1;
    base_state=likelystates(size(likelystates,2));
    for i=1:size(likelystates,2)-1
        if likelystates(1,i)==likelystates(1,i+1)
            temp(count)=likelystates(1,i);
            temp_data(count)=intis1(i,1);
            count=count+1;
        else
            temp(count)=likelystates(1,i);
            temp_data(count)=intis1(i,1);
            split{count2}=temp;
            split_data{count2}=temp_data;
            count2=count2+1;
            temp=[];
            temp_data=[];
            count=1;
        end
    end
    total_photons(num)=0;
    blink_counts(num,1)=0;
    for i=1:size(split,2)
        if split{i}(1)==base_state
            off_time=horzcat(off_time,[size(split{i},2)]);
        else
            blink_counts(num,1)=blink_counts(num,1)+1;
            on_time1=horzcat(on_time1,[size(split{i},2)]);
            intesity=horzcat(intesity,[sum(split_data{i})/size(split{i},2)]);
            total_photons(num)=total_photons(num)+sum(split_data{i});
            STD2_INTENSITY=horzcat(STD2_INTENSITY,[(std(split_data{i}))^2/(sum(split_data{i})/size(split{i},2))]);
            STD_INTENSITY=horzcat(STD_INTENSITY,[std(split_data{i})/(sum(split_data{i})/size(split{i},2))]);
        end
    end
    hold on
    plot((1:t)*timestep,likest);
    hold off
    img_save_name=[file_name{1}(1:end-9) '\' '_Trace_' char(string(num))  '.png'];
    saveas(fig,img_save_name);
    close(fig)
    disp(['The number of state is ',num2str(length(intensity))])
    disp(['The probability transition matrix is',newline])
    disp(num2str(TRANS_EST2))
    end
    %saveas(I, ['E:\My work\data\trace\TEMP\', 'trace', int2str(num), '-intensity-',num2str(intensity1),'-on_time-',num2str(on_time),'.fig']);
end
on_time1=on_time1'.*ex_time;
off_time=off_time'.*ex_time;
intesity=intesity';
total_photons=total_photons';
STD_INTENSITY=STD_INTENSITY';
STD2_INTENSITY=STD2_INTENSITY';
% save data
writematrix(series',[file_name{1}(1:end-9) '\' 'Trace' '.xlsx'])
writematrix(on_time1,[file_name{1}(1:end-9) '\' 'Ontime' '.xlsx'])
writematrix(off_time,[file_name{1}(1:end-9) '\' 'Offtime' '.xlsx'])
writematrix(blink_counts,[file_name{1}(1:end-9) '\' 'Blink counts' '.xlsx'])
%% manual filtter
chosen_n=[];
for i=1:250
    if gcf().Number~=1
        chosen_n=horzcat(chosen_n,gcf().Number);
        close gcf
    elseif i==1
        chosen_n=horzcat(chosen_n,gcf().Number);
        close gcf
    end
end
%series_new=zeros(size(chosen_n,2),size(series,2));
for i=1:size(chosen_n,2)
    series_new(i,:)=series(chosen_n(i),:);
end
%series_new1=series_new';


%% manual analysis
%series_new=series_new';
split={};
for i=1:size(series_new,1)
    threshold=200;
    temp=[];
    for j=1:size(series_new,2)
        if series_new(i,j)<=threshold
            series_new(i,j)=0;
        else
            temp=horzcat(temp,series_new(i,j));
        end
    end
    split{i}=temp;
end
on_time=[];
intensity=[];
std_intensity=[];
std2_intensity=[];
for i=1:size(split,2)
    on_time(i)=size(split{1,i},2)*ex_time;
    intensity(i)=mean(split{1,i});
    std_num(i)=std(split{1,i});
    total_photon(i)=sum(split{1,i});
    std_intensity(i)=std(split{1,i})/mean(split{1,i});
    std2_intensity(i)=std(split{1,i})^2/mean(split{1,i});
end
on_time=on_time';
intensity=intensity';
total_photon=total_photon';
std_intensity=std_intensity';
std2_intensity=std2_intensity';


%% data for plotting;
clear dataX dataY
dataY=series(134,:);
dataX=(1:size(series,2)).*ex_time;
dataY=dataY';
dataX=dataX';
%% set parameter for HiddenMarhov_analysis
% clear all
% load('matlab.mat');
clear Hint HMMint
intis1=series(113,:);
intis1=intis1';
%intis1=smooth(intis1,1);
%intis1=medfilt1(intis1,5);
%Set the time step
timestep=0.1;

t=length(intis1);%The length of observation sequence

ChooseMethod=1;
% Method 1: Using hierarchical clustering
if ChooseMethod==1
    X=intis1;
    %Plot the observed data curve
    figure('visible','on')
    plot((1:t)*timestep,intis1)
    xlabel('time(s)')
    ylabel('Intensity')
    set(gca,'FontSize',15)
    
    Y=pdist(X,'euclidean');%Compute the pairwise distances
    Z=linkage(Y,'average');%Generate clustering hierarchical tree according to distance information
    %Calculate cophenetic correlation coefficient for the hierarchical cluster tree, a larger value indicates that the tree fits the distance well
    C=cophenet(Z,Y)
    
    %%%%%%Natural Divisions
    %Set cutoff threshold according from standard deviation and cophenetic correlation coefficient
    I1=find(intis1~=0);
    cutoff=20*std(intis1(I1((end-50):end)));%If the number of categories
    %is small,try increasing the cutoff,otherwise reducing the cutoffD
    %clustering result
    seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
    seq=seq';
    
    %Determine emission value
    F1=length(unique(seq));
    intensity=zeros(1,F1);
    for i=1:F1
        intensity(i)=mean(intis1(find(seq==i)));
    end
    
    Ctrace=zeros(1,length(intis1));
    for i=1:length(intensity)
        Ctrace(find(seq==i))=intensity(i);
    end
    R=corrcoef(Ctrace',intis1);
    i=0.1;
    R(1,2)

    while R(1,2)<0.88
        cutoff=(20-i)*std(intis1(end-50:end));
        i=i+0.1;
        seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
        seq=seq';
        F1=length(unique(seq));
        intensity=zeros(1,F1);
        for j=1:F1
            intensity(j)=mean(intis1(find(seq==j)));
        end
        for j=1:F1
            Ctrace(find(seq==j))=intensity(j);
        end
        R=corrcoef(Ctrace,intis1);
    end

    %%Method 2: Finding the peak position of the histogram of the observation series
elseif ChooseMethod==2
    %Set peak threshold, number of bins, minimum peak separation, method to determine emission, method to determine dividing line
    
    %Set Minimum peak height
    PeakThreshold=1;
    %Set number of bins
    Nhist=round(t/20);
    %Set Minimum peak separation
    MinPeakdistance=round(Nhist/20);
    %The emission is determined by median(method=0),peak(method=1), or mean(method=2) value
    method=0;
    %The deviding line is determined by peak width(bmethod=0) or mean value of two neighbouring peaks(bmethod=1)
    bmethod=1;
    
    %Seq:Emission sequence
    %location: Emission intensity
    [seq,intensity]=HMMEmitTrc(intis1,timestep,PeakThreshold,Nhist,MinPeakdistance,method,bmethod);
end

%%%%%% Baum-welch algorithm is used to estimate the transfer probability, emission probability 
lloc=length(intensity);
TRANS_GUESS = eye(lloc)*0.8+rand(lloc,lloc)*0.2;
EMIS_GUESS = eye(lloc)*0.8+rand(lloc,lloc)*0.2;
[TRANS_EST2, EMIS_EST2] = hmmtrain(seq, TRANS_GUESS, EMIS_GUESS);

%%%%%% The state sequence is estimated by Viterbi algorithm
likelystates = hmmviterbi(seq, TRANS_EST2, EMIS_EST2);

Draw analysis results
hold on
likest=likelystates ;
for i=1:length(intensity)
likest(find(likelystates==i))=intensity(i);
end

plot((1:t)*timestep,likest);
disp(['The number of state is ',num2str(length(intensity))])
disp(['The probability transition matrix is',newline])
disp(num2str(TRANS_EST2))

Hint=likelystates;
for i=length(intensity):(-1):1
Hint(find(Hint==i))=intensity(i);
end

HMMint=Hint';%State sequence
%% Intensity
raw_series_I=[];series_I=[];
loc_total_I = loc_molecule(data(:,:,1), ex_time);
% loc_single = unique(loc_total_o(:,1:2), 'rows');  % remove duplication of localizations
raw_series_I = time_trace(loc_total_I, data(:,:,:),2);
[~, frames] = size(raw_series_I);
% loc_red = unique(loc_red(:,1:2), 'rows');  % remove duplication of localizations
% raw_series = time_trace(loc_red, data,2);
series_I = raw_series_I - median(raw_series_I(:,frames-10:frames), 2);
% cut_series = series(:,1:2000);
%% rate convert
series=series';
rate=2;
tst=[];
for k=1:size(series,1)
    ts=series(k,:);
    n=1;
    for i=1:floor(size(series,2)/rate)
    tst(k,n)=sum(ts(i:i+rate-1));
    n=n+1;
    end
end
seriest=tst;
for i=1:size(seriest,2)
    seriest(2,i)=seriest(1,i);
    seriest(1,i)=0.1*i*rate;
end
seriest=seriest';
%% load data
clear all
file_name = 'D:\My work\data\trace\OSS+T\20240116\Cy5_300x3dilu_OSS+T_20240116_6mW_1000f_0p1s_startf1_fov1.sif';
s_frame =1;  % start frame

[data, ex_time, gainDAC] = pro_data(file_name, s_frame);    %ex_time: expousure time
data=data(:,:,:);
data0=data(:,:,10);
%% plot frame
data0=data(:,:,50);
% for i=1:size(data0,1)
%     for j=1:size(data0,2)
%         if data0(i,j)>200
%             data0(i,j)=data0(i,j)-90;
%         end
%         if data0(i,j)>300
%              data0(i,j)=data0(i,j)-190;
%         end
%     end
% end
figure
fig=imagesc(data0);
%hold on;plot([110 144],[ 99 147],'o','MarkerSize',8,'Color','Yellow')
%axis(theaxis)
%zlim([0,140]);
axis equal
axis off
xlim([1,size(data0,1)])
ylim([1,size(data0,1)])
%fig.Parent.CLim=[30,150]
% caxis([100,8000])
set(gca,'FontSize',15)
colormap gray
%% get picture of frame
data0=data(:,:,1);
% for i=1:size(data0,1)
%     for j=1:size(data0,2)
%         if data0(i,j)>200
%             data0(i,j)=data0(i,j)-90;
%         end
%         if data0(i,j)>300
%              data0(i,j)=data0(i,j)-190;
%         end
%     end
% end
figure
fig=imagesc(data0);
%hold on;plot([110 144],[ 99 147],'o','MarkerSize',8,'Color','Yellow')
%axis(theaxis)
%zlim([0,140]);
axis equal
axis off
xlim([1,256])
ylim([1,256])
%fig.Parent.CLim=[30,150]
% caxis([100,8000])
set(gca,'FontSize',15)
colormap gray
loc_frame=loc_molecule(data0, ex_time);
check_localization(data0, loc_frame); 
% check_localization(data(:,:,check_frame), loc_red(loc_red(:,3)==check_frame,:)); 
axis square
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf,'-property','FontSize'), 'FontSize', 20);
colorbar
hold on
figure('name','intensity distribution','visible','on')
% fig=imagesc(data0);
% idlmapnames=transpose(h5read('idlcolormaps.h5','/name'));
% mymapindex=6;%%%%%%%%%%%%%%%%
% mycolormap=idlmapnames{mymapindex};
% disp([string(mymapindex), mycolormap]);
% fig.Parent.CLim=[30,170]
% idlrgbtables=h5read('idlcolormaps.h5', '/rgbt');
% colormap(double(idlrgbtables(:,:,mymapindex))/255.);
% colorbar 
% xlabel('x','FontSize',15);ylabel('y','FontSize',15);zlabel('intensity','FontSize',15);
% %xlim([1,r]);ylim([1,ss])
% shading interp
% axis equal

%caxis([60 1300])
brightness=[];
for i=1:size(loc_frame,1)
    brightness(i)=sum(data0(loc_frame(i,1)-3:loc_frame(i,1)+3, loc_frame(i,2)-3:loc_frame(i,2)+3),'all');
end
brightness=brightness';
%% SLB (filter)
tracenum=[];
k=1;
for i=1:size(trackedPar,2)
    if size(trackedPar(i).xy,1)>30 && size(trackedPar(i).xy,1)<1000
        tracenum(k)=i;
        k=k+1;
    end
end
%% SLB (plot1)
maxframe=300;
figure
for j=3:3
    cmap=jet;
    plot(trackedPar(tracenum(j)).xy(:,1),trackedPar(tracenum(j)).xy(:,2))
    hold on
    for i=1:size(trackedPar(tracenum(j)).xy,1)
        scatter(trackedPar(tracenum(j)).xy(i,1),trackedPar(tracenum(j)).xy(i,2),100,'MarkerEdgeColor',[.4 .4 .4],'MarkerFaceColor',cmap(ceil(i/maxframe*256),:),'LineWidth',1)
    end
end
%% SLB (plot2)
maxframe=400;
cmap=jet;
figure
for j=1:1
    plot(trackedPar(tracenum(j)).xy(:,1),trackedPar(tracenum(j)).xy(:,2),'color',cmap(ceil(size(trackedPar(tracenum(j)).xy,1)/maxframe*256),:),'LineWidth',0.5)
    hold on
end
%% MSD (all trajectories in one plot)
MSD=[];
frame=300;
%dis_t=[];
for i=1:frame
    MSD(i,1)=i;
    MSD(i,2)=0;
    k=0;
    for j=1:size(tracenum,2)
        if size(trackedPar(tracenum(j)).xy,1)>=i
            MSD(i,2)=MSD(i,2)+norm(trackedPar(tracenum(j)).xy(i,1:2)-trackedPar(tracenum(j)).xy(1,1:2))^2;
            k=k+1;
        end
    end
    MSD(i,2)=MSD(i,2)/k;
end
%% MSD (single trajectory)
tracenum=[365];
i=1;
X_data=[];
Y_data=[];
temp=[];
for j=1:size(trackedPar(tracenum(i)).xy,1)
   for k=1:size(trackedPar(tracenum(i)).xy,1)
       if j<=k
        temp(j,k-j+1)=norm(trackedPar(tracenum(i)).xy(j,1:2)-trackedPar(tracenum(i)).xy(k,1:2))^2;
       end
   end
   % X_data(j)=(j-1)*0.1;
   % Y_data(j)=norm(trackedPar(tracenum(i)).xy(j,:)-trackedPar(tracenum(i)).xy(1,:))^2;
end
for j=1:size(temp,1)
    X_data(j)=(j-1)*0.1;
    Y_temp=0;
    for k=1:size(temp,2)
        if temp(k,j)~=0
            Y_temp=Y_temp+temp(k,j);
        end
    end
    Y_data(j)=Y_temp/(size(temp,2)-j+1);
end
X_data=X_data';
Y_data=Y_data';
% data_MSD=[];
% for i=1:size(X_data,1)
%     data_MSD(i,1)=X_data(i);
%     data_MSD(i,2)=Y_data(i);
% end
plot(X_data,Y_data)
%% MSD&Diffusion constant (statics)
temp=[];
A=[];
D=[];
for i=1:size(tracenum,2)
     X_data=[];
     Y_data=[];
    for j=1:size(trackedPar(tracenum(i)).xy,1)
       for k=1:size(trackedPar(tracenum(i)).xy,1)
           if j<=k
            temp(j,k-j+1)=norm(trackedPar(tracenum(i)).xy(j,1:2)-trackedPar(tracenum(i)).xy(k,1:2))^2;
           end
       end
       % X_data(j)=(j-1)*0.1;
       % Y_data(j)=norm(trackedPar(tracenum(i)).xy(j,:)-trackedPar(tracenum(i)).xy(1,:))^2;
    end
    for j=1:size(temp,1)
        X_data(j)=(j-1)*0.1;
        Y_temp=0;
        for k=1:size(temp,2)
            if temp(k,j)~=0
                Y_temp=Y_temp+temp(k,j);
            end
        end
        Y_data(j)=Y_temp/(size(temp,2)-j+1);
    end
    X_data=X_data';
    Y_data=Y_data';
    f=fittype('4*D*t^a','independent','t','coefficients',{'D','a'});
    try
        cfun=fit(X_data(1:30),Y_data(1:30),f,'StartPoint',[0,1])
        A_temp=cfun.a;
        D_temp=cfun.D;
    catch
        A_temp=-1;
        D_temp=-1;
    end
    if A_temp~=-1&&D_temp~=-1
        A=horzcat(A,[A_temp]);
        D=horzcat(D,[D_temp]);
    end
end
A=A';
D=D';
%% get intensity for SLB
data=[];data1=[];
file_name = 'F:\data\trace\SLB\20240221\Cy5-tracking_OSS_1000f_1p5x_fov2.sif';
s_frame =1;  % start frame

[data1, ex_time, gainDAC] = pro_data(file_name, s_frame);    %ex_time: expousre time
raw_start=110;
line_start=150;
data=data1(raw_start:raw_start+255,line_start:line_start+255,:);
for i=1:size(trackedPar,2)
    for j=1:size(trackedPar(i).xy,1)
        try
            trackedPar(i).xy(j,3)=sum(data1(raw_start+floor(trackedPar(i).xy(j,1)-3)-1:raw_start+floor(trackedPar(i).xy(j,1)+3)-1,line_start+floor(trackedPar(i).xy(j,2)-3)-1:line_start+floor(trackedPar(i).xy(j,2)+3)-1,trackedPar(i).Frame(j)),'all');
        end
    end
end
%%
SNR=[];
k=0;
for i=1:size(trackedPar,2)
    if size(trackedPar(i).xy,2)==3 && size(trackedPar(i).xy,1)>=30
        k=k+1;
        SNR(k,1)=sum(trackedPar(i).xy(:,3))/size(trackedPar(i).xy,1)/std(trackedPar(i).xy(:,3));
    end
end
%% load data (cell)
clear all
file_name = 'H:\My work\data\trace\tublin\U2OS_HESi5b_200x_808nm_0p1s_500f_vol_1.sif';
s_frame =1;  % start frame

[data1, ex_time, gainDAC] = pro_data(file_name, s_frame);    %ex_time: expousre time
%% localing frame by frame (cell)
clear loc_total_o
loc_total_o = loc_molecule(data1(:,:,1:1000), ex_time);

%% extract fluorescence trace
raw_series=[];
j=43;
raw_series = time_trace(loc_filter_first(j,1:2), data(:,:,:),2);
loc_temp=[];
jump=2;
a=[];
[loc_max,a] = max(data1(raw_start+loc_filter_first(j,1)-1-3-jump:raw_start+loc_filter_first(j,1)-1+3+jump,line_start+loc_filter_first(j,2)-1-3-jump:line_start+loc_filter_first(j,2)-1+3+jump,1),[],'all');
loc_temp=[raw_start+loc_filter_first(j,1)-1-3-jump-1+floor(a/((3+jump)*2+1))+1,line_start+loc_filter_first(j,2)-1-3-jump-1+a-floor(a/((3+jump)*2+1))*((3+jump)*2+1)];
for i=1:999
    raw_series(2,i)=sum(data1(loc_temp(1)-3:loc_temp(1)+3,loc_temp(2)-3:loc_temp(2)+3,i),'all');
    [loc_max,a] = max(data1(raw_start+loc_filter_first(j,1)-1-3-jump:raw_start+loc_filter_first(j,1)-1+3+jump,line_start+loc_filter_first(j,2)-1-3-jump:line_start+loc_filter_first(j,2)-1+3+jump,i+1),[],'all');
    raw_temp=raw_start+loc_filter_first(j,1)-1-3-jump-1+a-floor(a/((3+jump)*2+1))*((3+jump)*2+1);
    line_temp=line_start+loc_filter_first(j,2)-1-3-jump-1+floor(a/((3+jump)*2+1))+1;
    loc_mean=mean(data1(raw_temp-3:raw_temp+3,line_temp-3:line_temp+3,i+1),'all');
    if loc_max>2*loc_mean&&norm([raw_temp,line_temp]-loc_temp)<4&&norm([raw_temp,line_temp]-[raw_start+loc_filter_first(j,1)-1,line_start+loc_filter_first(j,2)-1])<4
            loc_temp=[raw_temp,line_temp];
    end
    raw_series(2,i+1)=sum(data1(loc_temp(1)-3:loc_temp(1)+3,loc_temp(2)-3:loc_temp(2)+3,i+1),'all');
    loc_temp
end
figure
plot(0.1:0.1:0.1*size(raw_series,2),raw_series(1,:))
figure
plot(0.1:0.1:0.1*size(raw_series,2),raw_series(2,:))
%% snapshoot for checkframe
figure('visible','on')
imagesc(data1(:,:,1))
axis square
% title('control 638ex');
colormap('hot')
colorbar

set(gca, 'FontName', 'Arial');
set(gca,'YDir','normal');
set(gcf,'position',[300 100 600 500]);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
% saveas(gcf, [file_name, '\control 638ex', 'tif']);
%%
SNR=[];
for i=1:size(trackedPar,2)
    try
        SNR(i,1)=sum(trackedPar(i).xy(:,3))/size(trackedPar(i).xy,1)/std(trackedPar(i).xy(:,3));
    end
end
