%% load data (Note that data is in absolute e- unit)
clear all;clc;
file_name = '/Volumes/SMILeSSD/PacBio/UnknownSample_Dilute_ForLoading/20251219/3-singlemolecule_532nm_0p1S_2000frames_0p5mW_fov2.sif';
s_frame = 1;
[data,ex_time,gainDAC] = pro_data(file_name,s_frame);
data = data(128:384,128:384,:);

%% calculate avg.brightness vs time
avg_intensity = squeeze(mean(data,[1,2]))';
max_intensity = max(avg_intensity);
avg_intensity = avg_intensity/max_intensity;
figure;
plot(1:length(avg_intensity), avg_intensity, 'LineWidth', 1.8);
xlabel('Frame', 'FontSize', 14);
ylabel('Brightness', 'FontSize', 14);
title('Brightness Over Frames', 'FontSize', 16);
grid on;



%% snapshot of single frame
% BG = mean(data(:,:,end-10:end),'all');
% a = sum(data(:,:,210:230) ,3);
a = data(:,:,100);
figure
imagesc(a);
axis square; colorbar;
% set figure
set(gca,'FontName','Times New Roman');
set(gca,'YDir','Normal');
set(gcf,'Position',[500 300 500 400]);
set(gca,'FontSize',15);
colormap gray;
clim([0 300]);
%saveas(gcf,'FrameX.tif');% save figure

%% locate molecule
clear loc_total;
loc_total = loc_molecule(data(:,:,:),ex_time);

%% localization total figure
% figure('name','intensity distribution','visible','on')
imagesc(data(:,:,10));% a surface diagram of intensity distribution\
xlabel('x','FontSize',15);ylabel('y','FontSize',15);zlabel('intensity','FontSize',15);
shading interp
axis equal
xlim([1,256])
ylim([1,256])
hold on
scatter(loc_total(:,2),loc_total(:,1),'o','r');


%% check field (catch specific field)
loc_x = 108;
loc_y = 78;
data_1 = data(loc_x-4:loc_x+4,loc_y-4:loc_y+4,:);


%% field frame (check single frame of specific field)
check_frame = 1200; 
imagesc(data_1(:,:,check_frame));
txt = ['t = ',num2str(check_frame*ex_time),'s'];
text(0.05, 0.93, txt, 'Units', 'normalized', 'FontSize', 60,'Color','w','FontWeight','bold','FontName','Times New Roman');
colormap gray;
axis square;
axis off;
clim([0 150]);
set(gcf,'position',[256 128 512 512]);
set(gca,'position',[0 0 1 1]);


%% check single frame (with local maxium marked)
check_frame = 1300; %check single frame data
check_localization(data(:,:,check_frame),loc_total(loc_total(:,3)==check_frame,:));
title(['Frame',' ',num2str(check_frame)]);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',15);
colormap gray;
clim([0 80]);


%% point view (integrate over 7*7 field)
point = [47 51];
point_raw_series = time_trace(point, data(:,:,:));

[~, frames] = size(point_raw_series);
point_raw_series = point_raw_series - mean(point_raw_series(:,frames-20:frames)); % The Last 20 frames are used to substract baseline from each trace.

% bg_point = [225 68];
% 
% bg_point_raw_series = time_trace(bg_point, data(:,:,:));
% [~, frames] = size(bg_point_raw_series);
% bg_point_series = bg_point_raw_series - mean(bg_point_raw_series(:,frames-20:frames));
% point_raw_series = point_raw_series-bg_point_series;


figure;
plot(point_raw_series);


%% movie
clear v;
figure('Name','Intensity distribution by time-movie');
set(gcf,'Position',[800 600 600 500]);
    v = VideoWriter('RNO-COOH-10mW.avi');
open(v);
firstfr = 1;
total_fr = 1000;
set(gca,'nextplot','replace');
%title('Intensity distribution by time');
%fig = imagesc(data(:,:,1));
%colormap gray;
for i = firstfr:total_fr
    fig = imagesc(data(:,:,i));
    colorbar;
    axis square;
    xlim([1,256]);
    ylim([1,256]);
    set(gca,'FontSize',15);
    set(gca,'YDir','normal');
    %fig.Parent.Clim = [50,300];
    colormap gray;
    clim([0 150]);
    title(['time = ',num2str((i-firstfr)*ex_time),'s']);
    pause(0.04);
    frame = getframe(gcf);
    writeVideo(v,frame);
end
%set(gcf,'CloseRequstFcn','closereq');
close(v);

%% partial movie
clear v;
figure('Name','Intensity distribution by time-movie');
set(gcf,'Position',[800 600 600 500]);
v = VideoWriter('RNO-COOH-10mW.avi');
open(v);
firstfr = 541;
total_fr = 541;
set(gca,'nextplot','replace');
%title('Intensity distribution by time');
%fig = imagesc(data(:,:,1));
%colormap gray;
index = 31;
x =loc_filter_first(index,1);
y =loc_filter_first(index,2);
filed_location = data(x-6:x+6,y-6:y+6,:);
for i = firstfr:total_fr
    fig = imagesc(filed_location(:,:,i));
    colorbar;
    axis square;
    xlim([1,13]);
    ylim([1,13]);
    set(gca,'FontSize',15);
    set(gca,'YDir','normal');
    %fig.Parent.Clim = [50,300];
    colormap gray;
    clim([0 80]);
    title(['time = ',num2str((i-firstfr)*ex_time),'s']);
    pause(0.04);
    frame = getframe(gcf);
    writeVideo(v,frame);
end
%set(gcf,'CloseRequstFcn','closereq');
close(v);

%% on-state molecules counts (Showing Photon-bleaching)
loc_frame = loc_total(:,3);
freq = tabulate(loc_frame);% turn the amount of on-state molecules per frame to frequency
plot(freq(:,1).*ex_time,freq(:,2));
axis square;
xlabel('Time (s)','FontName','Times New Roman','FontSize',17);
ylabel('Population','FontName','Times New Roman','FontSize',17);
title('On-state molecules counts at 15mW','FontName','Times New Roman','FontSize',17);


%% DBSCAN of localized points
clear step r minn E on_intime
t=1;E=1;step=size(data,3);r=1;minn=2;loc_total_1=zeros(1,3);

while t<=E 
    clear db
    i=size(loc_total_1,1);k=size(loc_total_1,1);
    for i=i:size(loc_total,1)
        if loc_total(i,3)>=((t-1)*step+1)&&loc_total(i,3)<=(t*step)
            loc_total_1(k,1:3)=loc_total(i,:);
            k=k+1;
        end
    end
    disp('Load loc_total for dbscan complete') % load loc_total to a new matrix for dbscan
    db = dbscan(loc_total_1(:,1:2),r,minn);
    loc_total_db = loc_total_1;
    for i=1:size(db)
        loc_total_db(i,4)=db(i); % Append to each loc cluster information, db=[-1,1,2,...]
    end
    disp('DBSCAN complete')

    loc_db=sortrows(loc_total_db,4);  % Sort accroading to clusters (if a loc belongs to the same emitter)
    j=1;
    k=1;
    loc_temp=zeros(1,3);
    loc_filter_first=zeros(1,2);
    loc_filter_last=zeros(1,2);
    loc_filter_center=zeros(1,2);
    on_intime=zeros;
    for i = 1:size(loc_db,1)-1
        if loc_db(i,4)==loc_db(i+1,4)
            loc_temp(j,1:3)=loc_db(i,1:3);
            j=j+1;
        else
            loc_temp(j,1:3)=loc_db(i,1:3);
            loc_filter_first(k,:)=loc_temp(1,1:2);
            loc_filter_last(k,:)=loc_temp(end,1:2);
            loc_filter_center(k,:)=mean(loc_temp(:,1:2));
            on_intime(j)=loc_temp(1,3);
            k=k+1;
            j=1;
            loc_temp=zeros(1,3);
        end
    end
    loc_filter_last(1,:)=[];
    loc_filter_first(1,:)=[];
    loc_filter_center(1,:)=[];
    while loc_db(1,4)==-1
        loc_db(1,:)=[];
    end
    c{t}=loc_db;
    t=t+1;
    on_intime(1)=[];
    on_intime=on_intime.*ex_time;
    on_intime=on_intime';
    [loc_single,inx,p]=unique(loc_db(:,4), 'rows');
    for i = 1:size(inx,1)
        loc_single(i,1:4)=loc_db(inx(i),1:4);
    end
    loc_single=sortrows(loc_single,4);
end

for i=1:size(c,2)
    scatter(c{i}(:,2),c{i}(:,1),'.')
    axis equal
    xlim([1,256])
    ylim([1,256])
    hold on
    pause(1)
end

data_temp=zeros(256,256);
for i=1:256
    for j=1:256
        for k=1:10
            data_temp(i,j)=data_temp(i,j)+data(i,j,k);
        end
    end
end
figure('name','intensity distribution','visible','on')
% imagesc(data(:,:,1));% a surface diagram of intensity distribution\
xlabel('x','FontSize',15);ylabel('y','FontSize',15);zlabel('intensity','FontSize',15);
shading interp
axis equal
xlim([1,256])
ylim([1,256])
hold on
scatter(loc_filter_first(:,2),loc_filter_first(:,1),'o','r');

%% extract fluorescence trace
raw_series = time_trace(loc_filter_first(:,1:2), data(:,:,:));
[~, frames] = size(raw_series);
% 
series = raw_series - median(raw_series(:,frames-10:frames), 2);%remove background

% bg_point = [36 89];
% point_raw_series = time_trace(bg_point, data(:,:,:));
% [~, frames] = size(point_raw_series);
% point_series = point_raw_series - mean(point_raw_series(:,frames-20:frames));
% series = series-point_series;





%% manual filt trace

% fake_spot = [1;2;4;6;10;11;15;19;20;21;23;24;25;28;31;34;35;37;38;42;44;45;47;48;49;51;54;57;59;61;62;63;65;71;72;73;74;75;77;79;81;82;83;85;86;87;88;89;93;97;104;105;106;110;111;113;114;118;120;127;128;132;135;136;137;141;143;144;145;146;149;150;152;154;155;156;157;158;163;164;169;170;171;172;175;177;178;180;182;183;184;185;188;191;193;194;195;199;200;202;203;205;208;210;211;213;215;216;217;218;222;225;226;229;231;233;237;238;241;245;246;249;252;253;256;258;260;262;263;264;265;267;271;272;274;279;280;281;282;283;284;286;288;289;290;293;294;295;296;297;300;303;304;306;307;308;309;312;314;316;320;323;326;330;331;333;334;335;338;339;342;343;346;347;348;350;352;354;355;363;365;368;370;372;373;376;377;378;379;382;383;386;388;389;390;392;395;398;399;400;401;402;403;404;407;408;409;412;414;415;416;417;419;421;424;425;426;429;432;434;435;436;438;440;443;447;448;450;451;458;464;465;466;470;473;474;475;476;477;480;481;483;484;485;490;491;494;499;503;505;509;513;515;517;520;521;523;524;526;529;533;537;538;539;540;541;542;543;544;546;549;550;552;554;556;561;563;566;568];%input the index number of fake spots
% select_spot = [];
% 
% n = isempty(fake_spot);
% if n == 0
%     for i = 1:size(fake_spot,1)
%         series(fake_spot(i),1) = 9999.9999;
%     end
% 
% elseif n == 1
%     for i = 1:size(series,1)
%         if ~ismember(i,select_spot)
%             series(i,1) = 9999.9999;
%         end
%     end
% end

% 阈值筛去过暗的点
for i = 1:size(series,1)
    if series(i,1) < 240
        series(i,1) = 9999.9999;
    end
end

% 阈值筛去过亮的点
% for i = 1:size(series,1)
%     if any(series(i,:) > 650)
%         series(i,1) = 9999.9999;
%     end
% end


filt = series(:,1) == 9999.9999;
series(filt,:) = [];


%% filter by blinking times
for i = 1:length(switchtime)
    if 1 > switchtime(i) || switchtime(i) > 15
        series(i,1) = 9999.9999;
    end
end

filt = series(:,1) == 9999.9999;
series(filt,:) = [];

%% analysis of unblinking brightness & survival time

clear parameter

intensity_total = zeros(size(series,1),1);
survival_time = zeros(size(series,1),1);

for num = 1:size(series,1)

    for bleach_frame = 1:size(series,2)
        if series(num,bleach_frame)<150
            break
        end
    end
    intensity_total(num) = mean(series(num,1:bleach_frame));
    survival_time(num) = bleach_frame * ex_time;
end

max_num = max([length(survival_time),length(intensity_total)]);

survival_time = [survival_time; NaN(max_num-length(survival_time),1)];
intensity_total = [intensity_total; NaN(max_num-length(intensity_total),1)];

parameter = [survival_time, intensity_total];



%% analysis of blinking by threshold
on_time=[];
off_time=[];
intensity_total=[];
switchtime=[];
bleach_time=[];
clear parameter

lower_limit = 250;%!!!!!!!!!
% upper_limit = 600;

duty_cycle = zeros(size(series,1),1);
bleach_time = zeros(size(series,1),1);

for num = 1:size(series,1)
    valid_signal = series(num,:) > lower_limit;
    non_signal = ~valid_signal;
    [labels_1,num_segments_1] = bwlabel(valid_signal);
    [labels_2,num_segments_2] = bwlabel(non_signal);

    on_time_single = zeros(num_segments_1,1);
    off_time_single = zeros(num_segments_2-1,1);
    intensity_single = zeros(num_segments_1,1);

    for i = 1:num_segments_1
        signal_index = (labels_1 == i);
        on_time_single(i) = sum(signal_index) * ex_time;
        intensity_single(i) = mean(series(num,signal_index));
    end

    for i = 1:num_segments_2-1
        signal_index = (labels_2 == i);
        off_time_single(i) = sum(signal_index) * ex_time;
    end
    bleach_index = (labels_2 == num_segments_2);
    bleach_time(num) =  (size(data,3)-sum(bleach_index)) * ex_time;

    duty_cycle(num) = sum(on_time_single)/(sum(on_time_single) + sum(off_time_single));
    on_time = [on_time ; on_time_single];
    off_time = [off_time ; off_time_single];
    intensity_total = [intensity_total ; intensity_single];
    switchtime = [switchtime ; num_segments_1];

end

max_num = max([length(on_time),length(off_time),length(intensity_total),length(duty_cycle),length(switchtime),length(bleach_time)]);

on_time = [on_time; NaN(max_num-length(on_time),1)];
off_time = [off_time; NaN(max_num-length(off_time),1)];
intensity_total = [intensity_total; NaN(max_num-length(intensity_total),1)];
duty_cycle = [duty_cycle; NaN(max_num-length(duty_cycle),1)];
switchtime = [switchtime; NaN(max_num-length(switchtime),1)];
bleach_time = [bleach_time; NaN(max_num-length(bleach_time),1)];

parameter = [on_time, off_time, intensity_total,duty_cycle,switchtime,bleach_time];


%% filter fake_parameter
for i = 1:size(parameter,1)
    if parameter(i,2) < 0.2
        parameter(i,2) = NaN;
    end
    if parameter(i,1) > 15
        parameter(i,1) = NaN;
    end
    if parameter(i,4) == 1
        parameter(i,4) = NaN;
    end
end

%% check trace by HMM
clear Hint HMMint num state_num
good_num=1;
state_num = zeros(size(series,1),1);
survival_time = [];
for num=1:size(series,1)
    timestep=0.1;
    t=size(series,2);%The length of observation sequence
    X=series(num,:)';
    Y=pdist(X,'euclidean');%Compute the pairwise distances
    Z=linkage(Y,'average');%Generate clustering hierarchical tree according to distance information
    %Calculate cophenetic correlation coefficient for the hierarchical cluster tree, a larger value indicates that the tree fits the distance well
    C=cophenet(Z,Y)

    I1=find(X~=0);
    cutoff=20*std(X(I1((end-50):end)));

    seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
    seq=seq';
    intensity=zeros(1,length(unique(seq)));
    for i=1:length(unique(seq))
        intensity(i)=mean(X(find(seq==i)));
    end

    %intergrate trace
    Ctrace=zeros(1,length(X));
    for i=1:length(intensity)
        Ctrace(find(seq==i))=intensity(i);
    end

    R=corrcoef(Ctrace',X);
    i=0.1;

    while R(1,2)<0.88 %editable!!!!!!!!!!!!!!!!
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
    likest=zeros(size(likelystates));
    for i=1:length(intensity)
    likest(find(likelystates==i))=intensity(i);
    end
    % plot((1:t)*timestep,likest);
    % disp(['The number of state is ',num2str(length(intensity))]);
    % disp(['The probability transition matrix is',newline]);
    % disp(num2str(TRANS_EST2));

    Hint=likelystates;
    for i=length(intensity):(-1):1
    Hint(find(Hint==i))=intensity(i);
    end

    % HMMint=Hint';%State sequence
    % for i=1:size(HMMint,1)
    %     HMMint(i,2)=i;
    % end
    % 
    % j=1;k=1;q=0;
    % clear HMMint_temp Hint_new
    % for i=1:(size(HMMint,1)-1)
    %     if HMMint(i,1)==HMMint(i+1,1)
    %         HMMint_temp(j)=X(i);
    %         j=j+1;
    %     else
    %         HMMint_temp(j)=X(i);
    %         Hint_new(k,1)=std(HMMint_temp);
    %         Hint_new(k,2)=mean(HMMint_temp);
    %         Hint_new(k,3)=j;
    %         q=q+j;
    %         j=1;
    %         k=k+1;
    %         clear HMMint_temp
    %     end
    %     if i==size(HMMint,1)-1
    %         HMMint_temp(j)=X(i);
    %         Hint_new(k,1)=std(HMMint_temp);
    %         Hint_new(k,2)=mean(HMMint_temp);
    %         Hint_new(k,3)=j;
    %         q=q+j;
    %     end     
    % end

    % state_num(num)=length(intensity);
    % if length(intensity)<=2 | length(intensity)~=1
    %     good_num=good_num+1;
    % end

    state_num(num)=length(intensity);
            survival_frame = 0;
            if length(intensity) == 2 %| length(intensity)~=1
                good_num=good_num+1;

                %%读取状态序列分析survival time
                if likelystates(5) ~= 1
                    for bright_frame = 1:5
                        if likelystates(bright_frame) == 2
                            break
                        end
                    end
                    for bleach_frame = bright_frame:size(data,3)
                        if likelystates(bleach_frame) == 1
                            survival_frame = (bleach_frame-bright_frame);
                            break
                        end
                    end
                end
            end
            if survival_frame ~= 0
                survival_time = [survival_time;survival_frame*ex_time];
            end
    disp(['processing',num2str(num),'/',num2str(size(series,1))]);

end

state_num=state_num';
clear good_trace
k=1;
for i=1:size(state_num,2)
    if state_num(i) == 2
        good_trace(k)=i;
        k=k+1;
    end
end

%% filter trace
clear loc_filter_single
loc_HMM_filter = zeros(size(good_trace,1),2);
for i=1:size(good_trace,2)
    loc_HMM_filter(i,1:2)=loc_filter_first(good_trace(i),1:2);
end

raw_series = time_trace(loc_HMM_filter, data(:,:,:));
[~, frames] = size(raw_series);
series = raw_series - median(raw_series(:,frames-10:frames), 2);


%% plot trace with fitting
on_time=[];
off_time=[];
intensity_total=[];
clear parameter on_intime

for num=1:size(series,1)
    clear Hint HMMint 
    intis1=series(num,:);
    intis1=intis1';
  
    timestep=0.1;%Set the time step
    
    t=length(intis1);%The length of observation sequence

    ChooseMethod=1;
    % Method 1: Using hierarchical clustering
    %%Method 2: Finding the peak position of the histogram of the observation series

    if ChooseMethod == 1
        X=intis1;
        % plot the observed data curve
        figure('visible','on')
        I = plot((1:t)*timestep,intis1);
        xlabel('time(s)')
        ylabel('Intensity')
        set(gca,'FontSize',15)

        Y=pdist(X,'euclidean');%Compute the pairwise distances
        Z=linkage(Y,'average');%Generate clustering hierarchical tree according to distance information
        %Calculate cophenetic correlation coefficient for the hierarchical cluster tree, a larger value indicates that the tree fits the distance well
        C=cophenet(Z,Y)

        I1=find(X~=0);
        cutoff=20*std(X(I1((end-50):end)));

        seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
        seq=seq';
        intensity=zeros(1,length(unique(seq)));
        for i=1:length(unique(seq))
            intensity(i)=mean(X(find(seq==i)));
        end

        %intergrate trace
        Ctrace=zeros(1,length(X));
        for i=1:length(intensity)
            Ctrace(find(seq==i))=intensity(i);
        end

        R=corrcoef(Ctrace',X);
        i=0.1;

        while R(1,2)<0.88 %editable!!!!!!!!!!!!!!!!
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

    elseif ChooseMethod == 2
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

    k=1;j=0;
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
    if size(intensity_time,1) >1

    intensity_time=sortrows(intensity_time,2);
    on_time_1=sum(intensity_time(:,1))-intensity_time(1,1);
    on_time=horzcat(on_time,on_time_1);
    intensity_time(1,:)=[];
    intensity_time=sortrows(intensity_time,1);
    intensity1=intensity_time(size(intensity_time,1),2);
    intensity_total=horzcat(intensity_total,intensity1);
    likest=likelystates;
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
            on_time=horzcat(on_time,[size(split{i},2)]);
            intensity=horzcat(intensity,[sum(split_data{i})/size(split{i},2)]);
        end
    end
    % hold on
    % plot((1:t)*timestep,likest);
    % disp(['The number of state is ',num2str(length(intensity))])
    % disp(['The probability transition matrix is',newline])
    % disp(num2str(TRANS_EST2))
    else
    end
end

on_time=on_time'.*ex_time;
off_time=off_time'.*ex_time;
intensity_total=intensity_total';

max_num = max([length(on_time),length(off_time),length(intensity_total)]);

on_time = [on_time; NaN(max_num-length(on_time),1)];
off_time = [off_time; NaN(max_num-length(off_time),1)];
intensity_total = [intensity_total; NaN(max_num-length(intensity_total),1)];

parameter = [on_time, off_time, intensity_total];


%% plot picture(without fitting)
on_time=[];
off_time=[];
intesity=[];
for num=251:300%size(series,1)
    % clear all
    % load('matlab.mat');
    intis1=series(num,:);
    intis1=intis1';
    %intis1=smooth(intis1,1);
    %intis1=mentis1,5);
    %Set the time step
    timestep=0.1;
    t=length(intis1);
    figure('visible','on')
    plot((1:t)*timestep,intis1)
    xlabel('time(s)')
    ylabel('Intensity')
    set(gca,'FontSize',15)   
end


%% activation time
raw_series = time_trace(loc_filter_first(:,1:2), data(:,:,:));
[~, frames] = size(raw_series);
series = raw_series - mean(raw_series(:,frames-20:frames),2);
ac_time = zeros(size(series,1),1);

for spotnum = 1:size(series,1)
    for framenum = 1:frames
        if series(spotnum,framenum) >= 400 %threshold of on-state
            ac_time(spotnum) = framenum*ex_time;
            break
        else 
        end
    end
end


%% merge TIFF
tiffile = '/Volumes/Lu_yk_mac/Data/SMLM/20250220_RNO2-TPP_Hela/RNO2-TPP_5uM_Hela_100ms_532nm_5000frs_500uw_fov1_sf1_Tifstack_Uint16.tif';

info = imfinfo(tiffile);
numframes = numel(info);

img = imread(tiffile,1);
imgsum = double(img);
for i = 2:500
    frame = imread(tiffile,i);
    imgsum = imgsum+double(frame);
    clc
    disp([num2str(i),'/',num2str(numframes)])
end
%%
imagesc(imgsum);
colormap gray;
axis square;
axis off;
clim([500000 9000000]);
set(gcf,'position',[256 128 512 512]);
set(gca,'position',[0 0 1 1]);

%% merge sif and to tiff
clear all;clc;
file_name = '/Volumes/Lu_yk_mac/Data/SMLM/20250220_RNO2-TPP_Hela/RNO2-TPP_5uM_Hela_100ms_532nm_5000frs_500uw_fov1_sf1.sif';
s_frame = 11;
[data,ex_time,gainDAC] = pro_data(file_name,s_frame);
%%
merge_sif = sum(data(:,:,2:100),3);
t = size(merge_sif,3);
for i=1:t
    %st=255*si(:,:,i)./max(max(si(:,:,i)));
    st=merge_sif(:,:,i);
    imwrite(uint16(st),[file_name(1:end-4) '_Tifstack.tif'],'WriteMode','append')
end
%%
merge_sif = merge_sif - 10000;
imagesc(merge_sif);
colormap gray;
axis square;
axis off;
clim([10000 200000]);
set(gcf,'position',[256 128 512 512]);
set(gca,'position',[0 0 1 1]);

%% bright total
brightness = zeros(size(loc_total,1),1);
for i = 1:size(loc_total,1)
    brightness(i) = sum(data(loc_total(i,1)-2:loc_total(i,1)+2,loc_total(i,2)-2:loc_total(i,2)+2,loc_total(i,3)),'all')-mean(sum(sum(data(loc_total(i,1)-2:loc_total(i,1)+2,loc_total(i,2)-2:loc_total(i,2)+2,end-20:end))));
end































