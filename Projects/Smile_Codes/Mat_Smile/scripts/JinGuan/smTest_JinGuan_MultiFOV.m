%% Parameters
    clear all; clc;
    file_name = '/Volumes/SMILeSSD/PacBio/UnknownSample_Dilute_ForLoading/20260415/532nm(0p8mW)_1.5x_0p1S_5000frames_num3_alb_MultiFOV.sif';

%% Load data
    fprintf('Loading data: %s\n', file_name);
    if strcmpi(file_name(end-2:end), 'sif')
        [raw_data, ex_time, gainDAC] = io.readSIFData(file_name);
    else
        raw_data = io.readTiffStack(file_name);
        ex_time  = 0.1;     % <-- set manually for TIFF
    end

    interval = 0; % time interval between two adjacent exposure. In unit of seconds.
    oneFrameTime = ex_time+interval;

%% Windowing data
    windowed_raw_data = raw_data(100:412,100:412,:); % Should be a bit larger then 257*257 (128:384) for drift correction. Should be CENTERED.
    frames = size(windowed_raw_data,3);

%% save to tif
    outFile = [file_name(1:end-4) '_windowed.tif']; % Get rid of '.sif'.
    io.writeTiffStack(windowed_raw_data, outFile);

%% clear data
    clear raw_data

%% Calculate Correlation
    corr_array = zeros(frames-1,1);
    for b = 1:frames -1 
        corr_array(b) = corr2(windowed_raw_data(:,:,b), windowed_raw_data(:,:,b+1));
    end
    figure('Name','Correlation');
    plot(1:frames-1, corr_array, 'Color',[0.3 0.5 0.8], 'LineWidth', 0.8);

%% Calculate L1 Norm
    l1_norm = zeros(frames-1,1);
    for b = 1:frames -1 
        l1_norm(b) = sum(abs(windowed_raw_data(:,:,b) - windowed_raw_data(:,:,b+1)), 'all')/sum(windowed_raw_data(:,:,b),'all');
    end
    figure('Name','Delta L1 Norm');
    plot(1:frames-1, l1_norm, 'Color',[0.3 0.5 0.8], 'LineWidth', 0.8);

%% Calculate Laplacian variance
    H = fspecial('laplacian', 0.2);
    lap_var = zeros(frames,1);
    for b = 1:frames
        L = imfilter(windowed_raw_data(:,:,b), H, 'replicate');
        lap_var(b) = var(L(:));
    end
    figure('Name','Laplacian variance');
    plot(1:frames, lap_var, 'Color',[0.3 0.5 0.8], 'LineWidth', 0.8);

%% Check FOV
    l1_th = 0.21;

    mask = l1_norm(:) < l1_th ;
    d = diff([false; mask; false]);
    start_idx = find(d == 1);
    end_idx   = find(d == -1) - 1;
    len = end_idx - start_idx + 1;
    valid = len >= 100; % threshold
    sections = [start_idx(valid), end_idx(valid)];

%% Plot L1 norm
    figure('Name','L1 Norm');
    plot(1:length(l1_norm), l1_norm, 'Color', [0.3 0.5 0.8], 'LineWidth', 0.8);
    hold on;

    % threshold line
    yline(l1_th , '--r', 'Threshold', 'LineWidth', 1);

    % highlight valid sections
    for i = 1:size(sections,1)
        idx = sections(i,1):sections(i,2);
        plot(idx, l1_norm(idx), 'g', 'LineWidth', 2.5);
    end

    xlabel('Frame difference index');
    ylabel('Normalized L1 difference');
    title('Stable Sections in L1 Norm');
    grid on;
    hold off;

%% Manually set section
    % sections = [1 10; 413 423; 472 482; 595 605;763 773;1230 1240;1857 1867;2402 2412;2604 2614;2851 2861;3084 3094;3306 3316;4076 4086;4287 4297;4440 4450;4760 4770]; % a11_MULTI1
    %x = [192;466;610;809;960;1115;1253;1426;1607;1748;1880;2035;2214;2363;2599;2759;2910;3057;3248;3421;3560;3760;3938;4136;4302;4449;4626;4760];%532nm(1p5mW)_1.5x_0p1S_5000frames_num2_a11_MultiFOV.sif
    %x = [1;125;340;515;698;893;1112;1283;1428;1554;1694;1820;1973;2103;2264;2456;2625;2830;2995;3320;3489;3649;3778;3944;4246;4442;4571;4717;4863];%532nm(1p5mW)_1.5x_0p1S_5000frames_num3_alb_MultiFOV.sif
    x = [522;714;936;1170;1364;1542;1681;1823;1965;2098;2236;2593;2755;2886;2921;3086;3213;3351;3501;3643;3791;4157;4429;4610;4794];%532nm(0p8mW)_1.5x_0p1S_5000frames_num3_alb_MultiFOV.sif
    % x = [132;359;509;708;868;1012;1196;1354;1509;1710;1901;2028;2180;2320;2460;2654;2844;3075;3304;3473;3647;3790;3953;4149;4307;4530;4686;4812];%532nm(0p8mW)_1.5x_0p1S_5000frames_num2_a11_MultiFOV
    sections = [x zeros(size(x))];

%% check data
    start_frames = sections(:,1);
    for i = 1:numel(start_frames)
        start_frame = start_frames(i);
        figure;
        viz.plotImage(windowed_raw_data(:,:,start_frame), 1, 'hot',['Frame ',num2str(start_frame)],'max');
    end

%% localize
    lasting_frame = 100;
    roi = 201;
    k_sigma = 1.8;
    [rows, cols, ~] = size(windowed_raw_data);
    start_row = (rows-roi)/2+1;
    end_row = (rows+roi)/2;
    start_col = (cols-roi)/2+1;
    end_col = (cols+roi)/2;
    data = windowed_raw_data(start_row:end_row,start_col:end_col,:);
    searching_radius = 2; %pixel
    emitters_all = struct([]);
    brightness_all = [];
    background_all = [];
    sigma_all      = [];
    sigma_loc_all  = [];
    selected_data = [];
    delta_sum = zeros(lasting_frame, 2);

    for i = 1:numel(start_frames)
        start_frame = start_frames(i);
        end_frame = start_frame + lasting_frame-1; 
        data_temp = data(:,:,start_frame:end_frame);
        detected = detect.findMaxima(data_temp,k_sigma,5);
        localized = localize.LSQ1.locMolecules(data_temp,detected,5);
        [~, ~, loc_idx]=postproc.drift.correctDrift(data_temp,localized(:,[1,2,7]),delta_sum,5,roi);
        [localized] = localized (loc_idx,:);
        brightness = localized(:,3);
        background = localized(:,4);
        sigma      = localized(:,5);
        sigma_loc  = localized(:,6);

        emitters = postproc.emitter.findEmitters(data_temp, localized(:,[1,2,7]), 10, searching_radius);
        emitters = postproc.emitter.mergeEmitters(emitters, searching_radius);
        
        emitters_filt = postproc.emitter.filterEmitters_firstframe(emitters);
        % emitters_filt = postproc.emitter.filterEmitters_end(emitters_filt); 
        emitters_filt = postproc.emitter.filterEmitters_short(emitters_filt, 2,'consecutive');
        [emitters_filt,stats_jump] = postproc.emitter.filterEmitters_jumping(emitters_filt,3.5);

        % 当前循环之前已经累计了多少个 localization
        loc_offset = numel(brightness_all);

        % 修正 emitters_filt 中的 loc_idx，使其对应到全局数组
        for n = 1:numel(emitters_filt)
            emitters_filt(n).loc_idx = emitters_filt(n).loc_idx + loc_offset;
            emitters_filt(n).start_frame = start_frame;
        end
        
        % 合并emitters
        if isempty(emitters_all)
            emitters_all = emitters_filt;
        else
            emitters_all = [emitters_all, emitters_filt];
        end

        % 合并 localization 对应的属性数组
        brightness_all = [brightness_all; brightness];
        background_all = [background_all; background];
        sigma_all      = [sigma_all; sigma];
        sigma_loc_all  = [sigma_loc_all; sigma_loc];
        selected_data = cat(3, selected_data, data_temp);
    end

%% Filter emitters by brightness stability (short trace version)
   emitters_all_filt = emitters_all;
 
%% Filter by brightness
    bm = [];
    N = numel(emitters_all);
    for k = 1:N
        idx = emitters_all(k).loc_idx;
        bm(k) = mean(brightness_all(idx), 'omitnan') / ex_time;
    end
    thresh = prctile(bm, 99);
    keep_idx = bm <= thresh;
    emitters_all_filt = emitters_all(keep_idx);

%% Statistics
    stats = postproc.emitter.collectEmitterStatistics(emitters_all_filt,lasting_frame,ex_time,interval,brightness_all,sigma_all,sigma_loc_all,background_all);
    postproc.emitter.plotEmitterStatistics(stats,100);

%% Brightness Viz
    pd = analysis.photophys.fitBrightness_lognormal(stats.brightness_mean);
    n = numel(stats.brightness_mean);

    figure;
    histogram(stats.brightness_mean, 35);
    xlabel('Mean Intensity (photon/sec)');
    ylabel('Counts');
    title(sprintf('Histogram of mean intensities, n=%d', n))
    box on;

    medIntensity = median(stats.brightness_mean);
    fprintf('  median intensity = %.6f photon/sec\n', medIntensity);

%% Check all emitters
    figure;
    viz.plotImage(selected_data, 1:10, 'gray','Drift-corrected Image and Emitters','max');
    hold on;
    for k = 1:60 %numel(emitters_all)
        plot(emitters_all(k).col+0.5, emitters_all(k).row+0.5, '-', 'LineWidth', 1);
    end
    title('Emitter trajectories overlaid on image');

%% Extract trace 1
    raw_series = analysis.photophys.extractTrace(stats.pos_mean_px, data(:,:,:));
    bg = prctile(raw_series, 3, 2);  
    stats.trace = raw_series - bg;

%% Extract trace 2
    int_range = 4;
    bg_range = 7;
    raw_series = analysis.photophys.extractTrace(stats.pos_mean_px, data(:,:,:),int_range,true,bg_range,'median'); 
    stats.trace = raw_series;

%% Calculate and Check single emitter trace
    start_frames_em = [emitters_all_filt.start_frame]';
    time_window = [start_frames_em, start_frames_em + lasting_frame-1];
    N = numel(emitters_all_filt);
    stats.trace_temp = zeros(N,lasting_frame);
    for index = 1:N %61
        start_frame = time_window(index, 1);
        end_frame   = time_window(index, 2);
        trace_temp = stats.trace(index,start_frame:end_frame);
        stats.trace_temp(index,:) = trace_temp;
        %postproc.emitter.checkTrace(start_frame,end_frame,ex_time,trace_temp,stats.brightness_em(index,:),stats.pos_matrix(:,1,index),stats.pos_matrix(:,2,index));
    end

%% Check emitter video
    index = 2;
    start_frame = time_window(index, 1);
    end_frame   = time_window(index, 2);
    data_temp = data(:,:,start_frame:end_frame);
    postproc.emitter.checkMovie(index, data_temp, emitters_all_filt, stats.pos_mean_px, stats.pos_matrix, oneFrameTime, 'frameStep',1,'pauseTime', 0.2,'clim', [0 150],'preFrames', 0,'postFrames',0);

%% PhotoAnalysis of Trace
    % emitters_filt = analysis.photophys.analyzeStates(emitters_filt, 'HMM', 'bleachTail',50,'bgThreshold',100,'bicPenalty', 0.1, 'minStateSep', 1.0);
    states = analysis.photophys.analyzeStates(stats.trace_temp, 'CHANGEPOINT','penalty',3.0,'minSegLen',4,'mergeThr', 2.5, 'bleachTail',20,'bgThreshold',60, 'minStateSep', 2.5);
    %states = analysis.photophys.analyzeStates(stats.trace_temp, 'CHANGEPOINT','penalty',4.0,'minSegLen',4,'mergeThr', 3.0, 'bleachTail',50,'bgThreshold',60, 'minStateSep', 3.0);

%% Visualization
    % analysis.photophys.plotStates(emitters_filt,1:15);
    for index = 1:20
        analysis.photophys.plotState(states(index));
    end

%% Analysis brightest-state intensity in brightness_em
    N = numel(states);
    brightestEmMat = nan(N, 10);   % 第10列新增：前5帧平均亮度

    for n = 1:N
        brightestEmMat(n,1) = n;
        brightestEmMat(n,6) = states(n).nStates;
        brightestEmMat(n,7) = states(n).bleachFrame;
        brightestEmMat(n,8) = states(n).lifetime;

        sInfo = states(n).stateInfo;

        % brightness_em 第 n 条轨迹
        yb = stats.brightness_em(n,:);

        % ========= 新增：前5帧亮度统计 =========
        yb_valid = yb(~isnan(yb));
        if ~isempty(yb_valid)
            brightestEmMat(n,10) = mean(yb_valid(1:min(5,numel(yb_valid))), 'omitnan');
        end

        % 跳过没有有效状态的 emitter
        if isempty(sInfo)
            continue;
        end

        % 找最亮态：meanIntensity 最大的那个 state
        meanList = [sInfo.meanIntensity];
        [~, idxBright] = max(meanList);

        brightLabel = sInfo(idxBright).label;

        brightestEmMat(n,2) = brightLabel;
        brightestEmMat(n,3) = sInfo(idxBright).meanIntensity;
        brightestEmMat(n,4) = sInfo(idxBright).stdIntensity;
        brightestEmMat(n,5) = sInfo(idxBright).occupancy;

        % 在 sequence 中找到最亮态对应的 frame
        seq = states(n).sequence(:)';

        % 长度保护
        F = min(numel(seq), numel(yb));
        seq = seq(1:F);
        yb  = yb(1:F);

        idxFrames = (seq == brightLabel);

        if any(idxFrames)
            brightestEmMat(n,9) = mean(yb(idxFrames), 'omitnan');
        end
    end

%% MLE拟合的最亮态亮度
    brightestIntensity = brightestEmMat(:,9) / ex_time;
    brightestIntensity = brightestIntensity(~isnan(brightestIntensity) & brightestIntensity > 0);

    pd = analysis.photophys.fitBrightness_lognormal(brightestIntensity);
    n = numel(brightestIntensity);

    figure;
    histogram(brightestIntensity, 80);
    xlabel('Brightest State Intensity (photon/sec)');
    ylabel('Counts');
    title(sprintf('Histogram of brightest state intensities, n=%d', n))
    box on;

    medIntensity = median(brightestIntensity);
    fprintf('  median intensity = %.6f photon/sec\n', medIntensity);

%% 新增：前5帧亮度统计
    first5Intensity = brightestEmMat(:,10) / ex_time;
    first5Intensity = first5Intensity(~isnan(first5Intensity) & first5Intensity > 0);

    pd_first5 = analysis.photophys.fitBrightness_lognormal(first5Intensity);
    n_first5 = numel(first5Intensity);

    figure;
    histogram(first5Intensity, 40);
    xlabel('First 5-frame Mean Intensity (photon/sec)');
    ylabel('Counts');
    title(sprintf('Histogram of first 5-frame mean intensities, n=%d', n_first5))
    box on;

    medIntensity = median(first5Intensity);
    fprintf('  median intensity = %.6f photon/sec\n', medIntensity);