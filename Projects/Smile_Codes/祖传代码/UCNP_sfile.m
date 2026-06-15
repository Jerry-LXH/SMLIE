clear;
clc
tic
% 多进程处理
% if isempty(gcp('nocreate'))
%     pnumber=parpool(4);
% end
% 读取每个文件夹下的数据文件
path= '/Volumes/Lu_yk_ExFat/未命名文件夹/';
SMbrightness(path);
% Data2Video(path);
% SMdatafit(path);
toc

%%
function SMdatafit(path)
path=[path 'Figure/'];
% 对数据进行拟合
files = dir(fullfile(path, '*.*'));
% 获取该文件下的每个数据文件
datanames = {files(~[files.isdir]).name}';
% 获取每个文件的格式
[~, ~, fileExtension] = fileparts(datanames);
expectedFormat='.xlsx';
% 选取.xlsx文件
datanames=datanames(strcmpi(fileExtension, expectedFormat));
% 获取数据的数量
datanumber=size(datanames,1);
% 初始化数据大小
intensity=zeros(datanumber,1);
error=zeros(datanumber,1);
% 对数据进行拟合
for i=1:datanumber
    dataname=cell2mat(datanames(i,:));
    datapath=strcat(path,dataname);
    data_UCNP=readmatrix(datapath);
    [intensity(i),error(i)]=dataanalysis(data_UCNP);
end
datat=table(datanames,intensity,error);
data_ana_name=[path, 'data analysis.xlsx'];
writetable(datat,data_ana_name)
end


%%
function Data2Video(path)
files = dir(fullfile(path, '*0p1s*.sif*'));
% 获取该文件下的每个数据文件
datanames = {files(~[files.isdir]).name}';
% 获取数据的数量
datanumber=size(datanames,1);
mkdir([path,'Video/']);
for i=1:datanumber
    dataname_mat=cell2mat(datanames(i,:));
    data_save_name=[path,'Video/',dataname_mat(1:end-4) '.avi'];
    outputVideo = VideoWriter(data_save_name);
    seq=readsif(strcat(path,dataname_mat));
    outputVideo.FrameRate = 10/seq.exposureTime;  % 每秒帧数
    imgdata=seq.imageData;
    open(outputVideo);
    % 创建视频帧
    numFrames = size(seq.imageData,3);  % 视频帧数
    %numFrames=500;
    startframe=1;
    framesize=size(seq.imageData,3);
    for k =1:framesize
        intmean=mean(mean(seq.imageData(:,:,k)));
        if intmean>500
            startframe=k;
            break
        end
    end
    for j = startframe:min(numFrames,framesize)
        % 完全去除白边
        % imgdata(:,:,i)=rescale(imgdata(:,:,i));
        % imshow(imgdata(:,:,i),'border','tight','initialmagnification','fit');
        % 白边未完全去除
        imagesc(imgdata(:,:,j));
        set(gcf,'Position',[300,200,420,420]);
        axis square;
        colormap("gray");
        set(gca,'looseInset',[0 0 0 0])
        set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
        if j==startframe
           crange = clim(gca);
        end
        set(gca, 'CLim', crange);
        %clim([0 10000])
        frame=getframe(gcf);
        % 写入当前帧到视频
        writeVideo(outputVideo, frame);  % 注意需要转换为uint8类型
        % 显示进度
        fprintf('Frame %d/%d\n', j, min(numFrames,framesize));
    end
    % 关闭视频
    close(outputVideo);
    disp('视频创建完成。');

end
end




%%
function SMbrightness(path)
files = dir(fullfile(path, '*.sif*'));
% 获取该文件下的每个数据文件
datanames = {files(~[files.isdir]).name}';
% 获取数据的数量
datanumber=size(datanames,1);
mkdir([path,'Figure/']);
%parfor x=1:3
for x=1:datanumber
    % 转变路径格式为char
    dataname_mat=cell2mat(datanames(x,:));
    % 设置每个文件的
    data_save_name=[path,'Figure/',dataname_mat(1:end-4) '.xlsx'];
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
    count=1;
    for m=startframe:min(startframe+25,framesize)
        % 显示现在正在处理的数据
        % powernum=strfind(dataname_mat, 'mA');
        disp(['现在是 ' dataname_mat(1:end-4) ' 第' char(string(m)) '帧 ' '...'])
        ima1=double(seq.imageData(:,:,m));
        % 选取数据范围
        area_r=127;
        dataedge=50;
        area_c=127;
        %ima1=ima1(dataedge:255+dataedge,1:256);%green
        % ima1=ima1(dataedge:255+dataedge,area_c:128+255);%red
        ima1=ima1(area_r:area_r+255,area_c:area_c+255);%单通道
        % 根据CCDfangda 比例，对数据处理
        ima1=(ima1.*5.75)./seq.gainDAC;
        % 数据转置
        ima1_1=ima1';
        % 二值化
        S=rescale(ima1);
        % 确定并去除背景值
        % bg=imopen(S,strel('disk',4));
        % S=S-bg;
        % find the threshold
        level=graythresh(S);
        bw=imbinarize(S,level*0.7);
        % 当图像比较乱时，可以通过该方法起到一定效果
        % bw = bwmorph(bw,"clean");
        % bw = bwmorph(bw,"open");
        [label,~]=bwlabel(bw,4);
        % 将regionprops函数得到的数据存入到stats
        stats=regionprops('table',label,ima1,'Centroid','Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');
        num_d=find(stats.Area>25);% UCNP
        % 判断找到的点的数目，如果点的数目太少，减小level，重新识别
        if size(stats,1)<2
            bg=imopen(S,strel('disk',4));
            S=S-bg;
            bw=imbinarize(S,level*0.7);
            % bw = bwmorph(bw,"open");
            [label,~]=bwlabel(bw,4);
            % 将regionprops函数得到的数据存入到stats
            stats=regionprops('table',label,ima1,'Centroid','Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');
            num_d=find(stats.Area>25);
        end
        % 找寻多个光斑连到一起的区域，重复上述运算，并将其分开
        if (isempty(num_d)==0)
            for i=1:size(num_d,1)
                state_d=stats(num_d(i),:);
                x_d=state_d.BoundingBox;
                state_a=p_d(x_d,ima1_1);
                stats=[stats;state_a];
            end
        end
        % 删去光斑连在一起的数组
        stats(num_d,:)=[];
        % 当识别的点比较多时，即识别错了
        if size(stats,1)>50
            stats(stats.MaxIntensity<1.5*mean(ima1(:)),:)=[];
            stats(stats.Area<5,:)=[];
        end
        stats(stats.MaxIntensity<2*mean(ima1(:)),:)=[];
        stats(stats.Area<5,:)=[];
        % 去除重叠的值
        stats(stats.MaxIntensity>2.5*mean(stats.MaxIntensity),:)=[];
        % 数据数量
        sizeg=size(stats,1);
        fr_D=zeros(sizeg,7);
        % 初始化FWHM
        FWHM=zeros(sizeg,2);
        % 初始化中心点位置
        center_x=zeros(sizeg,1);
        center_y=zeros(sizeg,1);
        % 图像的背景均值
        bg_mean=mean(ima1_1(:));
        % 确定每个数据点的范围
        x_edge=stats.BoundingBox(:,1:4);
        x_edge(:,1:2)=x_edge(:,1:2)+0.5;
        % 每个区域的row参数
        row_edge=[x_edge(:,1)-1 x_edge(:,1)+x_edge(:,3)+1];
        % 每个区域的col参数
        col_edge=[x_edge(:,2)-1 x_edge(:,2)+x_edge(:,4)+1];
        % 去除超出数组范围的参数
        col_edge(col_edge>256)=256;
        col_edge(col_edge<1)=1;
        row_edge(row_edge>256)=256;
        row_edge(row_edge<1)=1;
        % 每个PSF的初始位置
        centers = stats.Centroid;
        sto=ones(sizeg,1);
        % 设置拟合初始位点 sigma=1.2
        startpoint_G=[centers(:,2),centers(:,1),1.2.*sto,1.2.*sto,stats.MaxIntensity*2*pi*1.2*1.2,zeros(sizeg,1),bg_mean.*sto];% 一一对应fittype中的参数
        % 初始化每个点的积分值
        Integral=zeros(sizeg,1);
        % 记录每个点的半径
        diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
        radii = diameters/2+1;
        fitResult=[];
        % 只画每一组数据的第一帧
        if count==1
            % 绘图
            fig=figure('Visible','off');
            imagesc(ima1);colormap('gray');axis square;colorbar;
            clim([0 500]);
            % 圈出每个数据点
            hold on
            viscircles(centers,radii,'color','#4F94CD','LineWidth',1.5,'EnhanceVisibility',false);
            % plot(fr_D(:,2),fr_D(:,1),'.','Color','red')
            hold off
            % 命名并保存
            title(dataname_mat(1:end-4));
            % saveas(fig,fig_save_name);
            img_save_name=[path 'Figure/' dataname_mat(1:end-4) '_StartFra_' char(string(m)) '.png'];
            saveas(fig,img_save_name);
            close(fig);
            % 保存原始图像
            fig=figure('Visible','off');
            imagesc(seq.imageData(:,:,m));colormap('gray');axis square;colorbar;
            hold on
            rectangle('Position',[area_r,area_c,255,255],'LineWidth',1,'Visible','on','EdgeColor','red');
            hold off
            % 显示每个图像的名字
            % 保存原始数据的图像
            title([dataname_mat(1:end-4) '_fit']);
            img_save_name=[path 'Figure/' dataname_mat(1:end-4) '_StartFra_' char(string(m)) '_orginal' '.png'];
            saveas(fig,img_save_name);
            close(fig)
            count=count+1;
        end

        for i=1:sizeg
            % 根据范围选取数据
            imgGrayData = ima1_1(row_edge(i,1):row_edge(i,2),col_edge(i,1):col_edge(i,2));
            % 创建网格
            [X,Y]=meshgrid(col_edge(i,1):col_edge(i,2),row_edge(i,1):row_edge(i,2));
            % 拟合前准备
            [xData, yData, zData] = prepareSurfaceData( X,Y, imgGrayData);
            % 设置fittype
            ft=fittype(@(u1,u2,ss1,ss2,A,th,C,x,y)A*exp(-(cos(th)*x+sin(th)*y-u1).^2./(2*ss1^2)- ...
                (-sin(th)*x+cos(th)*y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2)+C,'independent', ...
                {'x','y'},'dependent','z');
            % 设置startpoint
            % 用fit函数进行拟合
            try
                fitResult=fit([xData,yData],zData,ft,'startpoint',startpoint_G(i,:));
            catch Error
                disp(['现在是 ' dataname_mat(1:end-4) ' 第' char(string(m)) '帧 ' '第' char(string(i)) '拟合点报错' ]);
                disp(Error.message);
            end
            % 拟合得到的数据存入fr_D
            fr_D(i,:)=coeffvalues(fitResult);
            % 然后将拟合数据的每个参数分开
            numc=num2cell(fr_D(i,:),1);
            [u1,u2,ss1,ss2,A,th,~]=deal(numc{:});
            % 转换矩阵，将拟合的u1和u2转变为中心点
            cc=[cos(th),sin(th);-sin(th),cos(th)];
            xy=[u1;u2];
            xxyy=cc^(-1)*xy;
            % 计算得到的中心点分别存入center_x/y中
            center_x(i)=xxyy(1);
            center_y(i)=xxyy(2);
            % 计算每个数据点的积分强度
            % 先确定函数的中心点，并将其代入到函数中
            [u1,u2]=deal(xxyy(1),xxyy(2));
            % 先确定积分函数的形式
            f=@(x,y)A*exp(-(x-u1).^2./(2*ss1^2)-(y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2);
            % 选取-10/+10σ的范围进行积分
            Integral(i)=integral2(f,u1-5*ss1,u1+5*ss1,u2-5*ss2,u2+5*ss2);
            % 对每次拟合的原始数据和拟合结果进行绘图
%             figure
%             f=surf(X,Y,imgGrayData);
%             shading interp
%             plot(fitResult,[xData,yData],zData)
        end

        fr_D(:,1:2)=[center_x center_y];
        % 每个数据点的半径
        n=fr_D;
        [ss1,ss2]=deal(n(:,3),n(:,4));
        FWHM(:,1)=2.*ss1.*sqrt(-2*log(1/2));% X轴的半峰宽
        FWHM(:,2)=2.*ss2.*sqrt(-2*log(1/2));% Y轴的半峰宽
        stats.FWHM_X=FWHM(:,1);
        stats.FWHM_y=FWHM(:,2);
        % Integral=Integral./seq.exposureTime;
        stats.Integral=Integral;
        % 计算偏转角度Theta
        stats.Theta=fr_D(:,6).*360./pi;
        % Excel_save(stats,filename);
        % if isempty(Integral)==0
        %     try
        %         writematrix(Integral,data_save_name,"WriteMode","append");
        %         %writematrix(fr_D,data_save_name,"WriteMode","append");
        %     catch
        %         pause(1);
        %         writematrix(Integral,data_save_name,"WriteMode","append");
        %         %writematrix(Integral,data_save_name, 'Sheet', 1, 'Range', 'A1',"WriteMode","append");
        %         %writematrix(fr_D,data_save_name,"WriteMode","append");
        %     end
        % end
    end

    %通过HMM计算survival_time

        %载入数据
        data = seq.imageData(128:384,128:384,startframe+1:end);
        ex_time = seq.exposureTime;
        survival_time = [];

        %定位数据点
        clear loc_total;
        loc_total = loc_molecule(data(:,:,:),ex_time);

        %DBSCAN
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
                loc_total_db(i,4)=db(i);
            end
            disp('DBSCAN complete')

            loc_db=sortrows(loc_total_db,4);
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

        % for i=1:size(c,2)
        %     scatter(c{i}(:,2),c{i}(:,1),'.')
        %     axis equal
        %     xlim([1,256])
        %     ylim([1,256])
        %     hold on
        %     pause(1)
        % end

        data_temp=zeros(256,256);
        for i=1:256
            for j=1:256
                for k=1:10
                    data_temp(i,j)=data_temp(i,j)+data(i,j,k);
                end
            end
        end

        %提取trace
        raw_series = time_trace(loc_filter_first(:,1:2), data(:,:,:));
        [~, frames] = size(raw_series);

        series = raw_series - median(raw_series(:,frames-10:frames), 2);%remove background

        %HMM
        clear Hint HMMint num state_num
        good_num=1;
        state_num = zeros(size(series,1),1);
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
            % hold on
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
        max_num = max([length(Integral),length(survival_time)]);

        Integral = [Integral; NaN(max_num-length(Integral),1)];
        survival_time = [survival_time; NaN(max_num-length(survival_time),1)];
        

        parameter = [Integral, survival_time];
        writematrix(parameter,data_save_name,"WriteMode","append");
end
end


