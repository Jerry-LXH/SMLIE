function simpleVideoPlayer ()
    % 创建视频播放器GUI
    h.fig = figure('Position', [100 200 600 500]);
    axis equal
    set(gca().XAxis,'Visible','off')
    set(gca().YAxis,'Visible','off')
    h.selectfile=uicontrol('Style', 'pushbutton', 'String', '选择sif文件','Position', [10 480 100 20]);
    h.moleculeloc=uicontrol('Style', 'pushbutton', 'String', '单通道定位','Position', [120 480 80 20]);
    h.slider = uicontrol('Style', 'slider', 'Position', [175 20 200 20]);
    h.play = uicontrol('Style', 'pushbutton', 'String', '播放', 'Position', [100 20 50 20]);
    h.pause = uicontrol('Style', 'pushbutton', 'String', '暂停', 'Position', [400 20 50 20], 'Enable', 'off');
    h.selectedp = uicontrol('Style', 'pushbutton', 'String', '选择一个分子', 'Position', [210 480 100 20]);
    h.DCcalibration=uicontrol('Style', 'pushbutton', 'String', '双通道校准', 'Position', [320 480 100 20]);
    h.DClocalization=uicontrol('Style', 'pushbutton', 'String', '双通道定位', 'Position', [430 480 100 20]);
    h.chanel=1;
    h.setxlim=5000;

    
    % 设置文件选择回调函数
    set(h.selectfile,'Callback', @selectfileCallback);

    % 设置分子定位回调函数
    set(h.moleculeloc,'Callback', @moleculelocCallback);

    % 设置滑块回调函数
    set(h.slider, 'Callback', @sliderCallback);

    % 设置播放按钮回调函数
    set(h.play, 'Callback', @playCallback); 

    % 设置暂停按钮回调函数
    set(h.pause, 'Callback', @pauseCallback);

    % 设置分子选择按钮回调函数
    set(h.selectedp,'Callback', @selectedpCallback);

    % 设置双通道校准按钮回调函数
    set(h.DCcalibration,'Callback',@DCcalibrationCallback);

    % 设置双通道定位按钮回调函数
    set(h.DClocalization,'Callback',@DClocalizationCallback);

    % 存储GUI句柄
    guidata(h.fig, h);
end

function selectfileCallback(hObject,eventdata,handles)
    [filename,path]=uigetfile('041_0p1s_638nm_10mw_fov3_5000frs_sf1.sif', 'Select a MATLAB file',"D:\LHX\lhx-1-041-0703");
    s_frame=input('请设置起始帧数：','s');  % start frame
    s_frame=str2num(s_frame);
    [data1, ex_time, gainDAC] = pro_data(strcat(path,filename), s_frame);    %ex_time: expousre time
    for i=1:100
        intmean=mean(mean(data1(:,:,i)));
        if intmean>500/300*5.75
            startframe=i;
            break
        end
    end
    raw_start=128;
    red=data1(raw_start:raw_start+255,257:512,startframe:end);
    green=data1(raw_start:raw_start+255,1:256,startframe:end);
    imagesc(data1(:,:,1))
    axis square
    set(gca().XAxis,'Visible','off')
    set(gca().YAxis,'Visible','off')
    %fig.Parent.CLim=[20,280]
    colormap gray
    msgbox('数据载入完成', '提示');
    handles=guidata(get(hObject, 'Parent'));
    handles.data=data1;   
    handles.ex_time=ex_time;
    assignin('base','data1',handles.data);
    handles.dataexportDC=struct('loc_red',[],'loc_green',[],'trace_data_red',[],'trace_data_green',[]);
    handles.dataexportSC=struct('loc',[],'trace_data',[]);
    guidata(hObject,handles)    
end

function sliderCallback(hObject, eventdata,handles)
    % 滑块回调函数，更新视频帧
    handles=guidata(hObject);
    frameNum = round(get(hObject, 'Value')*size(handles.data,3))+1;
    imagesc(handles.data(:,:,frameNum));
    try
        viscircles([handles.loc_result(:,2),handles.loc_result(:,1)],3,'color','red','LineWidth',1.5,'EnhanceVisibility',false);
        catch
    end
    axis square
    set(gca().XAxis,'Visible','off')
    set(gca().YAxis,'Visible','off')
    colormap gray
    drawnow;
end

function playCallback(hObject, eventdata, handles)
    % 播放按钮回调函数
    handles=guidata(hObject); 
    set(hObject, 'Enable', 'off');
    set(handles.pause, 'Enable', 'on');
    handles.porc=1;
    guidata(hObject,handles)
    for i = round(get(handles.slider, 'Value')*size(handles.data,3))+1:size(handles.data,3)
        handles=guidata(hObject); 
        if handles.porc==1
            imagesc(handles.data(:,:,i));
            axis square
            set(gca().XAxis,'Visible','off')
            set(gca().YAxis,'Visible','off')
            colormap gray
            set(handles.slider, 'Value', i/size(handles.data,3));
            drawnow;
            pause(handles.ex_time);
        else
            break
        end
    end
    set(handles.pause, 'Enable', 'off');
    set(hObject, 'Enable', 'on');
    guidata(hObject,handles)
end

function pauseCallback(hObject, eventdata, handles)
    % 暂停按钮回调函数
    handles=guidata(hObject);
    set(handles.play, 'Enable', 'on');
    handles.porc=0;
    guidata(hObject,handles)
end

function moleculelocCallback(hObject, eventdata, handles)
    clear loc_total_green
    num = input('请输入定位帧数: ', 's');
    num=str2num(num)
    handles=guidata(hObject); 
    loc_total_green = loc_molecule(handles.data(:,:,1:num), handles.ex_time);
    clear step r minn E on_time
    t=1;E=1;step=size(handles.data,3);r=1;minn=2;loc_total=zeros(1,3);
    while t<=E
        clear db 
        i=size(loc_total,1);k=size(loc_total,1);
        for i=i:size(loc_total_green,1)
            if loc_total_green(i,3)>=((t-1)*step+1)&&loc_total_green(i,3)<=(t*step)
                loc_total(k,1:3)=loc_total_green(i,:);
                k=k+1;
            end
        end
    %     disp('Load loc_total complete')
        db=dbscan(loc_total(:,1:2),r,minn);
        for i=1:size(loc_total,1)
            loc_total(i,3)=loc_total_green(i,3);
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
         on_time=on_time.*handles.ex_time;
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
                data_temp(i,j)=data_temp(i,j)+handles.data(i,j,k);
            end
        end
    end
    viscircles([loc_filter_first(:,2),loc_filter_first(:,1)],3,'color','red','LineWidth',1.5,'EnhanceVisibility',false);
    handles.loc_result=loc_filter_first;
    intensity_SC=[];
    for num=1:size(loc_total_green,1)
        m1=mean(mean(handles.data(:,:,1)));
        data_signal=handles.data(loc_total_green(num,1)-3:loc_total_green(num,1)+3, loc_total_green(num,2)-3:loc_total_green(num,2)+3,loc_total_green(num,3));
        % cov1=cov(sum(data_signal));
        % cov2=cov(sum(data_signal,2));
        % [~,xm1]=max(max(data_signal));
        % [~,xm2]=max(data_signal(:,xm1));
        % ft=fittype(@(u1,u2,ss1,ss2,A,th,C,x,y)A*exp(-(cos(th)*x+sin(th)*y-u1).^2./(2*ss1^2)-(-sin(th)*x+cos(th)*y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2)+C,'independent',{'x','y'},'dependent','z');%加入扭转角th和背景C
        % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        % opts.Display = 'Off';
        % opts.Lower = [-Inf -Inf 0];
        % height=data_signal(4,4);
        % sss1=(cov1/(2*pi*height))^(1/4)/2.15;
        % sss2=(cov2/(2*pi*height))^(1/4)/2.15;
        % startPoints = [4 4 sss1 sss2 (height-m1)*2*pi*sss1*sss2  0 m1];%Set the initial fitting parameters
        % xdata=[];
        % ydata=[];
        % zdata=[];
        % k=0;
        % for i=1:7
        %     for j=1:7
        %         i;
        %         j;
        %         k=k+1;
        %         xdata(k)=loc_total_green(num,1)-3+i-1;
        %         ydata(k)=loc_total_green(num,2)-3+j-1;
        %         zdata(k)=data_signal(i,j);
        %     end
        % end
        % try
        %     [fitresult, ~] = fit( [xdata', ydata'],zdata', ft,'Start',startPoints)
        %     intensity_SC(num)=fitresult.A;
        % catch
        % intensity_SC_2(num)=sum(data_signal,"all")-sum(handles.data(loc_total_green(num,1)-3:loc_total_green(num,1)+3, loc_total_green(num,2)-3:loc_total_green(num,2)+3,end-10:end)/10,"all");
        % end
        intensity_SC_2(num)=sum(data_signal,"all")-sum(handles.data(loc_total_green(num,1)-3:loc_total_green(num,1)+3, loc_total_green(num,2)-3:loc_total_green(num,2)+3,end-10:end)/10,"all");
    end
    trace_all = time_trace(loc_filter_first(:,1:2), handles.data(:,:,1:end));
    trace_all = trace_all-median(trace_all(:,end-10,end),2);
    assignin('base', 'loc_total_o', loc_total_green);
    assignin('base', 'intensity_SC', intensity_SC');
    assignin('base', 'series', trace_all);
    assignin('base','ex_time',handles.ex_time);
    assignin('base', 'loc_filter_first', loc_filter_first);
    assignin('base', 'intensity_SC_2', intensity_SC_2');
    guidata(hObject,handles);
end

function selectedpCallback(hObject, eventdata, handles)
    handles=guidata(hObject);
    setxlim=handles.setxlim;
    if handles.chanel==1
        radius=3;
        time=size(handles.data,3);
        [x,y]=ginput(1);
        axis square
        set(gca().XAxis,'Visible','off')
        set(gca().YAxis,'Visible','off')
        colormap gray
        drawnow;
        hold on
        scatter(x,y,36,'yellow','x')
        hold off
        handles.x=x;
        handles.y=y;
        handles.data_point=handles.data(round(y)-radius:round(y)+radius,round(x)-radius:round(x)+radius,:);
        handles.data_point_trace = time_trace([round(y),round(x)], handles.data(:,:,1:end));
        handles.data_point_trace=handles.data_point_trace-median(handles.data_point_trace(:,end-10,end),2);
        handles.hfig_point = figure('Position', [100 200 600 500]);
        subplot(1,2,2)
        plot(0:0.1:(time-1)/10,handles.data_point_trace(1:time),Color='red',LineWidth=1.5)
        axis square
        xlim([0,setxlim])
        subplot(1,2,1)
        handles.figpoint=imagesc(handles.data_point(:,:,1));
        axis square
        set(gca().XAxis,'Visible','off')
        set(gca().YAxis,'Visible','off')
        %handles.figpoint.Parent.CLim=[20,200];
        colormap gray
        handles.sliderpoint = uicontrol('Style', 'slider', 'Position', [175 20 200 20]);
        handles.playpoint = uicontrol('Style', 'pushbutton', 'String', '播放', 'Position', [100 20 50 20]);
        handles.pausepoint = uicontrol('Style', 'pushbutton', 'String', '暂停', 'Position', [400 20 50 20], 'Enable', 'off');
        handles.SCdataexport = uicontrol('Style', 'pushbutton', 'String', '导出数据', 'Position', [20 20 50 20]);
        set(handles.SCdataexport,'Callback',@SCdataexportCallback);
        set(handles.sliderpoint, 'Callback', @sliderpointCallback);
        set(handles.playpoint,'Callback',@playpointCallback);
        set(handles.pausepoint,'Callback',@pausepointCallback);
        guidata(hObject,handles);
        guidata(handles.hfig_point,handles);
        
    else
        radius=3;
        time=size(handles.data,3);
        [x,y]=ginput(1);
        x=round(x)
        y=round(y)
        axis square
        set(gca().XAxis,'Visible','off')
        set(gca().YAxis,'Visible','off')
        colormap gray
        drawnow;
        hold on
        scatter(x,y,36,'yellow','x')
        hold off
        handles.data_point=handles.data(round(y)-radius:round(y)+radius,round(x)-radius:round(x)+radius,:);
        handles.data_point_trace = time_trace([round(y),round(x)], handles.data(:,:,1:end));
        handles.data_point_trace=handles.data_point_trace-median(handles.data_point_trace(:,end-10,end),2);
        handles.hfig_point = figure('Position', [100 200 600 500]);
        subplot(2,2,2)
        if  handles.chanelselect=='g'
            colorp{1}='green';
            colorp{2}='red';
        else
            colorp{1}='red';
            colorp{2}='green';
            x=x-256;
        end
        handles.x=x;
        handles.y=y;
        plot(0:0.1:(time-1)/10,handles.data_point_trace(1:time),Color=colorp{1},LineWidth=1.5)
        axis square
        xlim([0,setxlim])
        subplot(2,2,1)
        handles.figpoint=imagesc(handles.data_point(:,:,1));
        axis square
        set(gca().XAxis,'Visible','off')
        set(gca().YAxis,'Visible','off')
        %handles.figpoint.Parent.CLim=[20,200];
        colormap gray
        X = [1, x, y, x.^2, y.^2, x.*y, x.^3, y.^3, x.^2.*y, x.*y.^2];
        if handles.chanelselect=='g'
            x_other=X*handles.calibrationgreen2red(:,1);
            y_other=X*handles.calibrationgreen2red(:,2);
            handles.x_other=round(x_other);
            handles.y_other=round(y_other);
        else
            x_other=X*handles.calibrationred2green(:,1);
            y_other=X*handles.calibrationred2green(:,2);
            handles.x_other=round(x_other);
            handles.y_other=round(y_other);
        end
        if handles.chanelselect=='g'
            handles.data_point_other=handles.data(round(y_other)-radius:round(y_other)+radius,256+round(x_other)-radius:256+round(x_other)+radius,:);
            handles.data_point_trace_other = time_trace([round(y_other),256+round(x_other)], handles.data(:,:,1:end));
            handles.data_point_trace_other=handles.data_point_trace_other-median(handles.data_point_trace_other(:,end-10,end),2);
        else
            handles.data_point_other=handles.data(round(y_other)-radius:round(y_other)+radius,round(x_other)-radius:round(x_other)+radius,:);
            handles.data_point_trace_other = time_trace([round(y_other),round(x_other)], handles.data(:,:,1:end));
            handles.data_point_trace_other=handles.data_point_trace_other-median(handles.data_point_trace_other(:,end-10,end),2);
        end
        subplot(2,2,4)
        plot(0:0.1:(time-1)/10,handles.data_point_trace_other(1:time),Color=colorp{2},LineWidth=1.5)
        axis square
        xlim([0,setxlim])
        subplot(2,2,3)
        handles.figpoint=imagesc(handles.data_point_other(:,:,1));
        axis square
        set(gca().XAxis,'Visible','off')
        set(gca().YAxis,'Visible','off')
        %handles.figpoint.Parent.CLim=[20,200];
        colormap gray

        % % 检查定位位置
        % figure
        % imagesc(handles.data(:,:,1))
        % colormap gray
        % hold on
        % scatter(x,y,'green','o')
        % scatter(x_other,y_other,'red','o')

        handles.DCsliderpoint = uicontrol('Style', 'slider', 'Position', [175 10 200 20]);
        handles.DCplaypoint = uicontrol('Style', 'pushbutton', 'String', '播放', 'Position', [100 10 50 20]);
        handles.DCpausepoint = uicontrol('Style', 'pushbutton', 'String', '暂停', 'Position', [400 10 50 20], 'Enable', 'off');
        handles.DCdataexport = uicontrol('Style', 'pushbutton', 'String', '导出数据', 'Position', [20 10 50 20]);
        set(handles.DCdataexport, 'Callback', @DCdataexportCallback);
        set(handles.DCsliderpoint, 'Callback', @DCsliderpointCallback);
        set(handles.DCplaypoint,'Callback',@DCplaypointCallback);
        set(handles.DCpausepoint,'Callback',@DCpausepointCallback);
        guidata(hObject,handles);
        guidata(handles.hfig_point,handles);

        handles.fig_red_green=figure;

        hold on;
        plot(0:0.1:(time-1)/10,handles.data_point_trace(1:time),Color=colorp{1},LineWidth=1.5);
        %plot(0:0.1:(time-1)/10,likest_green(i,1:time),Color='green',LineWidth=1.5);
        plot(0:0.1:(time-1)/10,handles.data_point_trace_other(1:time),Color=colorp{2},LineWidth=1.5);
        %plot(0:0.1:(time-1)/10,likest_red(i,1:time),Color='red',LineWidth=1.5)
        hold off;
        guidata(hObject,handles)
        %assignin('base','data_trace', handles.dat1
        % aexportDC);

    end
end

function sliderpointCallback(hObject, eventdata, handles)
    handles=guidata(hObject);
    setxlim=handles.setxlim;
    figure(handles.hfig_point)
    time=size(handles.data,3);
    subplot(1,2,1)
    get(handles.sliderpoint, 'Value');
    frameNum = round(get(handles.sliderpoint, 'Value')*size(handles.data_point,3))+1;
    fig=imagesc(handles.data_point(:,:,frameNum));
    axis square
    set(gca().XAxis,'Visible','off')
    set(gca().YAxis,'Visible','off')
    fig.Parent.CLim=[20,200];
    colormap gray
    subplot(1,2,2)
    subplot(1,2,2)
    plot(0:0.1:(time-1)/10,handles.data_point_trace(1:time),Color='green',LineWidth=1.5)
    xlim([0,setxlim])
    axis square
    hold on
    scatter((frameNum-1)*handles.ex_time,handles.data_point_trace(frameNum),30,"red",'o')
    hold off
    % hold on
    % scatter((frameNum-1)*handles.ex_time,handles.data_point_trace(frameNum),30,"yellow",'filled')
    % try
    %     scatter((frameNum-1-1)*handles.ex_time,handles.data_point_trace(frameNum-1),0)
    % catch
    % end
    % hold off
    drawnow;
end

function playpointCallback(hObject, eventdata, handles)
    % 播放按钮回调函数
    handles=guidata(hObject); 
    setxlim=handles.setxlim;
    figure(handles.hfig_point)
    time=size(handles.data,3);
    set(hObject, 'Enable', 'off');
    set(handles.pausepoint, 'Enable', 'on');
    handles.porcpoint=1;
    guidata(hObject,handles)
    for i = round(get(handles.sliderpoint, 'Value')*size(handles.data_point,3))+1:size(handles.data_point,3)
        subplot(1,2,1)
        handles=guidata(hObject); 
        if handles.porcpoint==1 
            fig=imagesc(handles.data_point(:,:,i));
            axis square
            set(gca().XAxis,'Visible','off')
            set(gca().YAxis,'Visible','off')
            fig.Parent.CLim=[20,200]
            colormap gray
            set(handles.sliderpoint, 'Value', i/size(handles.data_point,3));
            subplot(1,2,2)
            plot(0:0.1:(time-1)/10,handles.data_point_trace(1:time),Color='green',LineWidth=1.5)
            axis square
            xlim([0,setxlim])
            hold on
            scatter((i-1)*handles.ex_time,handles.data_point_trace(i),30,"red",'o')
            hold off
            % hold on
            % scatter((i-1)*handles.ex_time,handles.data_point_trace(i),30,"yellow",'filled')
            % try
            %     scatter((frameNum-1-1)*handles.ex_time,handles.data_point_trace(frameNum-1),0)
            % catch
            % end
            % hold off
            drawnow;
            pause(handles.ex_time/10);
        else
            break
        end
    end
    set(handles.pausepoint, 'Enable', 'off');
    set(hObject, 'Enable', 'on');
    guidata(hObject,handles);
end

function pausepointCallback(hObject, eventdata, handles)
    % 暂停按钮回调函数
    handles=guidata(hObject);
    set(handles.playpoint, 'Enable', 'on');
    handles.porcpoint=0;
    guidata(hObject,handles)
end

function DCcalibrationCallback(hObject,eventdata,handles)
    handles=guidata(hObject);
    handles.chanel=2;
    folderPath = uigetdir("H:\My work\data\trace\Cy3_Cy5_fret\20241115");
    folderPath = horzcat(folderPath,'\')
    % 定义对话框的标题和问题
    dlgTitle = '通道选择';
    dlgMessage = '请选择优先定位的通道';
    
    % 显示询问对话框
    userSelection = questdlg(dlgMessage, dlgTitle, '红光通道','绿光通道','取消');
    
 
    % 根据用户的选择执行不同的操作
    switch userSelection
        case '红光通道'
            handles.chanelselect='r';
            handles.calibrationred2green=DCcalibrationred2green(folderPath);
        case '绿光通道'
            handles.chanelselect='g';
            handles.calibrationgreen2red=DCcalibrationgreen2red(folderPath);
    end
    guidata(hObject,handles)
end

function DClocalizationCallback(hObject,eventdata,handles)
    handles=guidata(hObject);
    % 定义对话框的标题和问题
    dlgTitle = '通道选择';
    dlgMessage = '请选择优先定位的通道';
    
    % 显示询问对话框
    userSelection = questdlg(dlgMessage, dlgTitle, '红光通道','绿光通道','取消','红光通道');
    
    % 根据用户的选择执行不同的操作
    switch userSelection
        case '红光通道'
            num = input('请输入红光定位帧数: ', 's');
            num=str2num(num)
            loc_total_red = loc_molecule(handles.data(:,257:512,1:num), handles.ex_time);
            clear step r minn E on_time
            t=1;E=1;step=500;r=1;minn=2;loc_total=zeros(1,3);
            while t<=E
                clear db 
                i=size(loc_total,1);k=size(loc_total,1);
                for i=i:size(loc_total_red,1)
                    if loc_total_red(i,3)>=((t-1)*step+1)&&loc_total_red(i,3)<=(t*step)
                        loc_total(k,1:3)=loc_total_red(i,:);
                        k=k+1;
                    end
                end
            %     disp('Load loc_total complete')
                db=dbscan(loc_total(:,1:2),r,minn);
                for i=1:size(loc_total,1)
                    loc_total(i,3)=loc_total_red(i,3);
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
                 on_time=on_time.*handles.ex_time;
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
                        data_temp(i,j)=data_temp(i,j)+handles.data(i,j,k);
                    end
                end
            end
            viscircles([loc_filter_first(:,2),loc_filter_first(:,1)],3,'color','red','LineWidth',1.5,'EnhanceVisibility',false);
            loc_red=[loc_filter_first(:,2),loc_filter_first(:,1)];
            loc_green_x= loc_red(:,1);
            loc_green_y = loc_red(:,2);
            X = [ones(size(loc_red_x)), loc_red_x, loc_red_y, loc_red_x.^2, loc_red_y.^2, loc_red_x.*loc_red_y, loc_red_x.^3, loc_red_y.^3, loc_red_x.^2.*loc_red_y, loc_red_x.*loc_red_y.^2];
            loc_green(:,1)=X*handles.calibrationred2green(:,1);
            loc_green(:,2)=X*handles.calibrationred2green(:,2);
            loc_green=[loc_green(:,2),loc_green(:,1)];
            loc_green=round(loc_green);
            viscircles([loc_green(:,2),loc_green(:,1)],3,'color','green','LineWidth',1.5,'EnhanceVisibility',false);
            guidata(hObject,handles);
        case '绿光通道'
            num = input('请输入绿光定位帧数: ', 's');
            num=str2num(num); 
            loc_total_green = loc_molecule(handles.data(:,1:256,1:num), handles.ex_time);
            clear step r minn E on_time
            t=1;E=1;step=500;r=1;minn=2;loc_total=zeros(1,3);
            while t<=E
                clear db 
                i=size(loc_total,1);k=size(loc_total,1);
                for i=i:size(loc_total_green,1)
                    if loc_total_green(i,3)>=((t-1)*step+1)&&loc_total_green(i,3)<=(t*step)
                        loc_total(k,1:3)=loc_total_green(i,:);
                        k=k+1;
                    end
                end
            %     disp('Load loc_total complete')
                db=dbscan(loc_total(:,1:2),r,minn);
                for i=1:size(loc_total,1)
                    loc_total(i,3)=loc_total_green(i,3);
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
                 on_time=on_time.*handles.ex_time;
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
                        data_temp(i,j)=data_temp(i,j)+handles.data(i,j,k);
                    end
                end
            end
            viscircles([loc_filter_first(:,2),loc_filter_first(:,1)],3,'color','green','LineWidth',1.5,'EnhanceVisibility',false);
            loc_green=[loc_filter_first(:,2),loc_filter_first(:,1)];
            loc_green_x= loc_green(:,1);
            loc_green_y = loc_green(:,2);
            X = [ones(size(loc_green_x)), loc_green_x, loc_green_y, loc_green_x.^2, loc_green_y.^2, loc_green_x.*loc_green_y, loc_green_x.^3, loc_green_y.^3, loc_green_x.^2.*loc_green_y, loc_green_x.*loc_green_y.^2];
            loc_red(:,1)=X*handles.calibrationgreen2red(:,1);
            loc_red(:,2)=X*handles.calibrationgreen2red(:,2);
            loc_red=[loc_red(:,2),loc_red(:,1)];
            loc_red=round(loc_green);
            viscircles([256+loc_red(:,1),loc_red(:,2)],3,'color','red','LineWidth',1.5,'EnhanceVisibility',false);
            guidata(hObject,handles);
        case '取消'
            disp('已取消');
    end
end

function DCsliderpointCallback(hObject, eventdata, handles)
    handles=guidata(hObject);
    setxlim=handles.setxlim;
    figure(handles.hfig_point)
    time=size(handles.data,3);
    if  handles.chanelselect=='g'
        colorp{1}='green';
        colorp{2}='red';
    else
        colorp{1}='red';
        colorp{2}='green';
    end
    subplot(2,2,1)
    get(handles.DCsliderpoint, 'Value');
    frameNum = round(get(handles.DCsliderpoint, 'Value')*size(handles.data_point,3))+1;
    fig=imagesc(handles.data_point(:,:,frameNum));
    axis square
    set(gca().XAxis,'Visible','off')
    set(gca().YAxis,'Visible','off')
    %fig.Parent.CLim=[20,200];
    colormap gray
    subplot(2,2,2)
    plot(0:0.1:(time-1)/10,handles.data_point_trace(1:time),Color=colorp{1},LineWidth=1.5)
    xlim([0,setxlim])
    axis square
    hold on
    scatter((frameNum-1)*handles.ex_time,handles.data_point_trace(frameNum),30,"black",'o')
    hold off
    subplot(2,2,3)
    get(handles.DCsliderpoint, 'Value');
    frameNum = round(get(handles.DCsliderpoint, 'Value')*size(handles.data_point_other,3))+1;
    fig=imagesc(handles.data_point_other(:,:,frameNum));
    axis square
    set(gca().XAxis,'Visible','off')
    set(gca().YAxis,'Visible','off')
    %fig.Parent.CLim=[20,200];
    colormap gray
    subplot(2,2,4)
    plot(0:0.1:(time-1)/10,handles.data_point_trace_other(1:time),Color=colorp{2},LineWidth=1.5)
    xlim([0,setxlim])
    axis square
    hold on
    scatter((frameNum-1)*handles.ex_time,handles.data_point_trace_other(frameNum),30,"black",'o')
    hold off
    % hold on
    % scatter((frameNum-1)*handles.ex_time,handles.data_point_trace(frameNum),30,"yellow",'filled')
    % try
    %     scatter((frameNum-1-1)*handles.ex_time,handles.data_point_trace(frameNum-1),0)
    % catch
    % end
    % hold off
    drawnow;
end
% 
function DCplaypointCallback(hObject, eventdata, handles)
    % 播放按钮回调函数
    handles=guidata(hObject); 
    setxlim=handles.setxlim;
    figure(handles.hfig_point)
    time=size(handles.data,3);
    if  handles.chanelselect=='g'
        colorp{1}='green';
        colorp{2}='red';
    else
        colorp{1}='red';
        colorp{2}='green';
    end
    set(hObject, 'Enable', 'off');
    set(handles.DCpausepoint, 'Enable', 'on');
    handles.porcpoint=1;
    guidata(hObject,handles)
    for i = round(get(handles.DCsliderpoint, 'Value')*size(handles.data_point,3))+1:size(handles.data_point,3)
         handles=guidata(hObject);
        if handles.porcpoint==1
            subplot(2,2,1)
            fig=imagesc(handles.data_point(:,:,i));
            axis square
            set(gca().XAxis,'Visible','off')
            set(gca().YAxis,'Visible','off')
            %fig.Parent.CLim=[20,200]
            colormap gray
            set(handles.DCsliderpoint, 'Value', i/size(handles.data_point,3));
            subplot(2,2,2)
            plot(0:0.1:(time-1)/10,handles.data_point_trace(1:time),Color=colorp{1},LineWidth=1.5)
            axis square
            xlim([0,setxlim])
            hold on
            scatter((i-1)*handles.ex_time,handles.data_point_trace(i),30,"black",'o')
            hold off
            subplot(2,2,3)
            fig=imagesc(handles.data_point_other(:,:,i));
            axis square
            set(gca().XAxis,'Visible','off')
            set(gca().YAxis,'Visible','off')
            %fig.Parent.CLim=[20,200]
            colormap gray
            %set(handles.DCsliderpoint, 'Value', i/size(handles.data_point_other,3));
            subplot(2,2,4)
            plot(0:0.1:(time-1)/10,handles.data_point_trace_other(1:time),Color=colorp{2},LineWidth=1.5)
            axis square
            xlim([0,setxlim])
            hold on
            scatter((i-1)*handles.ex_time,handles.data_point_trace_other(i),30,"black",'o')
            hold off
            % hold on
            % scatter((i-1)*handles.ex_time,handles.data_point_trace(i),30,"yellow",'filled')
            % try
            %     scatter((frameNum-1-1)*handles.ex_time,handles.data_point_trace(frameNum-1),0)
            % catch
            % end
            % hold off
            drawnow;
            pause(handles.ex_time/10);
        else
            break
        end
    end
    set(handles.DCpausepoint, 'Enable', 'off');
    set(hObject, 'Enable', 'on');
    guidata(hObject,handles);
end

function DCpausepointCallback(hObject, eventdata, handles)
    % 暂停按钮回调函数
    handles=guidata(hObject)
    set(handles.DCplaypoint, 'Enable', 'on');
    handles.porcpoint=0;
    guidata(hObject,handles)
    handles.porcpoint
end

function DCdataexportCallback(hObject, eventdata, handles)
   handles=guidata(hObject);
   if handles.chanelselect=='g'
        handles.dataexportDC.loc_green=horzcat(handles.dataexportDC.loc_green,[handles.y;handles.x]);
        handles.dataexportDC.loc_red=horzcat(handles.dataexportDC.loc_red,[handles.y_other;handles.x_other]);
    else
        handles.dataexportDC.loc_green=horzcat(handles.dataexportDC.loc_green,[handles.y_other;handles.x_other]);
        handles.dataexportDC.loc_red=horzcat(handles.dataexportDC.loc_red,[handles.y;handles.x+256]);
    end
    if handles.chanelselect=='g'
        handles.dataexportDC.trace_data_green=horzcat(handles.dataexportDC.trace_data_green,handles.data_point_trace');
        handles.dataexportDC.trace_data_red=horzcat(handles.dataexportDC.trace_data_red,handles.data_point_trace_other');
    else
        handles.dataexportDC.trace_data_red=horzcat(handles.dataexportDC.trace_data_red,handles.data_point_trace');
        handles.dataexportDC.trace_data_green=horzcat(handles.dataexportDC.trace_data_green,handles.data_point_trace_other');
    end
    assignin('base','data_trace', handles.dataexportDC);
    guidata(handles.fig,handles)
end

function SCdataexportCallback(hObject, eventdata, handles)
   handles=guidata(hObject);
   handles.dataexportSC.trace_data=horzcat(handles.dataexportSC.trace_data, handles.data_point_trace');
   handles.dataexportSC.loc=horzcat(handles.dataexportSC.loc,[round(handles.y);round(handles.x)]);
   guidata(handles.fig,handles)
   assignin('base','data_trace', handles.dataexportSC);
end