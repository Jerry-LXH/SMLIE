%% load data
clear all;clc
file_name = "E:\LHX041-PH5-TWEEN20-0709\041_1p0s_638nm_3mw_fov1_1000frs_sf1.sif";
s_frame =22;   % start frame

[data, ex_time, gainDAC] = pro_data(file_name, s_frame);    %ex_time: expousre time
data=data(128:383,128:383,:);
%% snapshoot for checkframe
figure('visible','on')
imagesc(data(:,:,34))
axis square
% title('control 638ex');
colorbar

set(gca, 'FontName', 'Arial');
set(gca,'YDir','normal');
set(gcf,'position',[300 100 600 500]);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
% saveas(gcf, [file_name, '\control 638ex', 'tif']);W
%% set parameter for loc_molecule
clear loc_total0
check_frame = 3;   % located frame of image
loc_total0 = loc_molecule(data(:,:,check_frame), ex_time);%start_frame
check_localization(data(:,:,check_frame), loc_total0); 
% check_localization(data(:,:,check_frame), loc_red(loc_red(:,3)==check_frame,:)); 
axis square
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf,'-property','FontSize'), 'FontSize', 20);
title(['frame', ' ', num2str(check_frame)]);
colorbar
%% localing frame by frame
clear loc_total_o
loc_total_o = loc_molecule(data(:,:,1:300), ex_time);
%% set parameter for DBSCAN
check_frame = 78;   % located frame of image

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
    db=DBSCAN(loc_total(:,1:2),r,minn);
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
imagesc(data(:,:,1));% a surface diagram of intensity distribution\
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
scatter(loc_filter_first(:,2),loc_filter_first(:,1),'o');


%cut_series = series(:,1:20);
%% extract fluorescence trace
raw_series = time_trace(loc_filter_first(:,1:2), data(:,:,:));
[~, frames] = size(raw_series);
series = raw_series - median(raw_series(:,frames-10:frames), 2);%remove background

%% check trace
clear Hint HMMint num state_num
good_num=1;
for num=1:300%size(series,1)
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

        while R(1,2)<0.88
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
    TRANS_GUESS = eye(lloc)*1+rand(lloc,lloc)*0.0000001;
    EMIS_GUESS = eye(lloc)*0.99+rand(lloc,lloc)*0.01;
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
    disp(['processing',num2str(num),'/',num2str(size(series,1))]);
end
state_num=state_num';
clear a
k=1;
for i=1:size(state_num,1)
    if state_num(i)>2 | state_num(i)==1 | state_num(i)==0
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
raw_series = time_trace(loc_single_filt, data(:,:,:));
[~, frames] = size(raw_series);
% loc_red = unique(loc_red(:,1:2), 'rows');  % remove duplication of localizations
% raw_series = time_trace(loc_red, data);
series = raw_series - median(raw_series(:,frames-10:frames), 2);

% cut_series = series(:,1:2000);
[~, frames] = size(series);
threshold_on = mean(series(:, frames-10:frames), 2) + 10*std(series(:, frames-8:frames), 0, 2);   % threshold for fluorescent state
% BGBGBG = mean(data,[1 2]);
% BGBG = reshape(BGBGBG,[1 size(BGBGBG,3)]);
% for i = 1:size(series,1)
%     series(i,:) = series(i,:)-BGBG;
% end
%% plot picture(with fitting)
on_time1=[];
off_time=[];
intesity=[];
clear intensity_ontime on_time
for num=1:size(series,1)
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
        figure('visible','on')
        I=plot((1:t)*timestep,intis1)
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
    for i=1:size(split,2)
        if split{i}(1)==base_state
            off_time=horzcat(off_time,[size(split{i},2)]);
        else
            on_time1=horzcat(on_time1,[size(split{i},2)]);
            intesity=horzcat(intesity,[sum(split_data{i})/size(split{i},2)]);
        end
    end
    % hold on
    % plot((1:t)*timestep,likest);
    % disp(['The number of state is ',num2str(length(intensity))])
    % disp(['The probability transition matrix is',newline])
    % disp(num2str(TRANS_EST2))
    end
    %saveas(I, ['E:\My work\data\trace\TEMP\', 'trace', int2str(num), '-intensity-',num2str(intensity1),'-on_time-',num2str(on_time),'.fig']);
end
on_time1=on_time1'.*0.1;
off_time=off_time'.*0.1;
intesity=intesity';
%% plot picture(without fitting)
on_time=[];
off_time=[];
intesity=[];
for num=100:140%size(series,1)
    % clear all
    % load('matlab.mat');
    intis1=series(num,:);
    intis1=intis1';
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
%% data for plotting
clear dataX dataY
dataY=series(99,:);
dataX=(1:size(series,2)).*ex_time;
dataY=dataY';
dataX=dataX';
%% set parameter for HiddenMarhov_analysis
% clear all
% load('matlab.mat');
clear Hint HMMint
intis1=series(1,:);
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
%%
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
on_time=[];
off_time=[];
intesity=[];
for i=1:size(split,2)
    if split{i}(1)==base_state
        off_time=horzcat(off_time,[size(split{i},2)]);
    else
        on_time=horzcat(on_time,[size(split{i},2)]);
        intesity=horzcat(intesity,[sum(split_data{i})/size(split{i},2)]);
    end
end
on_time=on_time.*0.1;
off_time=off_time.*0.1;
%%
merge_data=zeros(257,257);
for i=1:size(data,3)
    merge_data=merge_data+data(:,:,i);
end
loc_total0 = loc_molecule(merge_data, ex_time);
check_localization(merge_data, loc_total0); 
% check_localization(data(:,:,check_frame), loc_red(loc_red(:,3)==check_frame,:)); 
axis square
set(gca, 'FontName', 'Times New Roman');
set(findall(gcf,'-property','FontSize'), 'FontSize', 20);
title(['frame', ' ', num2str(check_frame)]);
colorbar
%%
intensity=[];
baseline=[];
timestep=0.1;
t=length(intis1);
intis1=series(79,:);
intis1=intis1';
for i=1:size(intis1,2)
    if intis1(i)>400
        intis2(i)=1;
        intensity=horzcat(intensity,[intis1(i)]);
    else
        intis2(i)=0;
        baseline=horzcat(baseline,[intis1(i)]);
    end
end
figure('visible','on')
plot((1:t)*timestep,intis1)
xlabel('time(s)')
ylabel('Intensity')
set(gca,'FontSize',15)
%%
rate=5;
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