%% Parameter settings
    clear all;clc;

    % file_name = '/Volumes/SMILeSSD/PacBio/UnknownSample_Dilute_ForLoading/20260308/532nm_1.5x_0p1S_5000frames_0p64mW_num5_alb_fov1.sif';
    % file_name = '/Volumes/SMILeSSD/PacBio/UnknownSample_Dilute_ForLoading/20260306/532nm_1.5x_0p1S_5000frames_0p46mW_num1_alb_fov2.sif';
    % file_name = '/Volumes/SMILeSSD/PacBio/UnknownSample_Dilute_ForLoading/20260306/532nm_1.5x_0p1S_5000frames_0p46mW_num3_a11_fov1.sif';
    file_name = '/Volumes/SMILeSSD/PacBio/UnknownSample_Dilute_ForLoading/20260415/532nm(0p8mW)_1.5x_0p1S_5000frames_num2_a11_fov2.sif';
    %file_name = '/Volumes/SMILeSSD/PacBio/UnknownSample_Dilute_ForLoading/20260320/532nm_1.5x_0p1S_5000frames_073mW_num6_bazhuayu_fov3.sif';
    

%% Load data (Note that data is in absolute e- unit)
    [raw_data,ex_time,gainDAC] = io.readSIFData(file_name);
    interval = 0; % time interval between two adjacent exposure. In unit of seconds.
    oneFrameTime = ex_time+interval;

%% Windowing data
    start_fm = 1;
    end_fm = 5000;
    time_window = start_fm:end_fm;
    windowed_raw_data = raw_data(100:412,100:412,time_window); % Should be a bit larger then 257*257 (128:384) for drift correction. Should be CENTERED.
    frames = size(windowed_raw_data,3);

%% clear data
    clear raw_data

%% save to tif
    outFile = [file_name(1:end-4) '_Tifstack.tif']; % Get rid of '.sif'.
    io.writeTiffStack(windowed_raw_data, outFile);
    
%% Detect molecule using local maxium
    clear loc_total;
    k_sigma = 2.0;
    uncorrected_detected_total = detect.findMaxima(windowed_raw_data,k_sigma); % Loc*3 (y,x,frame); Note x corresponds to columns, y corresponds to rows.

%% Visualize detections
    num_show = min(1, frames); % Select frames to visualize
    frames_to_show = 1:num_show; 
    figure;
    viz.plotImage(windowed_raw_data, frames_to_show, 'hot','Drift-uncorrected Image and Detections');
    hold on;
    viz.overlayLocs(uncorrected_detected_total,frames_to_show);

%% Localize using MLE estimator (slow)
    uncorrected_super_loc_total_raw = localize.MLE1.locMolecules(windowed_raw_data,uncorrected_detected_total,5);
    uncorrected_super_loc_total = localize.filterBadLocs(uncorrected_super_loc_total_raw);
    uncorrected_loc_total = uncorrected_super_loc_total(:,[1,2,7]);

%% Visualize the uncorrected subpixel locs
    num_show = min(5000, frames);
    frames_to_show = 1:num_show;
    figure;
    viz.plotImage(windowed_raw_data, frames_to_show, 'hot','Drift-uncorrected Image and Locs','max');
    hold on;
    viz.overlayLocs(uncorrected_loc_total,frames_to_show,true);

%% Save localization result
    [file_dir, file_base, ~] = fileparts(file_name);
    save_name = fullfile(file_dir, ...
        [file_base '_uncorrected_super_loc_total.mat']);
    save(save_name, 'uncorrected_super_loc_total', '-v7.3');
    fprintf('Detection saved to:\n%s\n', save_name);

%% Load localization result 
    [file_dir, file_base, ~] = fileparts(file_name);
    load_name = fullfile(file_dir, ...
        [file_base '_uncorrected_super_loc_total.mat']);
    S = load(load_name);
    uncorrected_super_loc_total = S.uncorrected_super_loc_total;
    uncorrected_super_loc_total = uncorrected_super_loc_total(ismember(uncorrected_super_loc_total(:,7), time_window), :);
    uncorrected_super_loc_total(:,7) = uncorrected_super_loc_total(:,7) - start_fm+1;
    uncorrected_loc_total = uncorrected_super_loc_total(:,[1,2,7]);
    fprintf('Detection loaded from:\n%s\n', load_name);

%% Drift estimation using NP-cloud
    drift_corr = true;
    uncertainty = mean(uncorrected_super_loc_total(:,6));
    if drift_corr == true
        delta_sum = postproc.drift.estimateDrift_cloudxy(windowed_raw_data, uncorrected_loc_total, 15, uncertainty, 30); % Note that loc_total have only 3 columns of (r,c,f).
    else
        delta_sum = zeros(frames, 2);
    end

%% Correct the Drift
    [data, loc_total, loc_idx] = postproc.drift.correctDrift(windowed_raw_data,uncorrected_loc_total,delta_sum,5,257);

    [corrected_super_loc_total] = uncorrected_super_loc_total(loc_idx,:); % this contains all 7 columns

%% Visualize the drift
    viz.plotTracking(delta_sum, oneFrameTime, [], 'Drift Tracking'); % visualize

%% save dirft-corrected data
    outFile = [file_name(1:end-4) '_Tifstack_corrected.tif'];
    io.writeTiffStack(data, outFile);

%% read tif stack
    inFile = [file_name(1:end-4) '_Tifstack_corrected.tif'];
    data = io.readTiffStack(inFile);

%% clear data
    clear windowed_raw_data

%% Other paras of MLE result
    brightness = corrected_super_loc_total(:,3);   % photon
    background = corrected_super_loc_total(:,4);   % background
    sigma      = corrected_super_loc_total(:,5);   % PSF width
    sigma_loc  = corrected_super_loc_total(:,6);   % localization precision

    figure;
    subplot(4,1,1);
    histogram(brightness, 180);
    xlabel('Photons (N)');
    ylabel('Count');
    title('Brightness distribution(by locs)');

    subplot(4,1,2);
    histogram(background, 140);
    xlabel('Photons (N)');
    ylabel('Count');
    title('Background distribution(by locs)');

    subplot(4,1,3);
    histogram(sigma, 120);
    xlabel('Sigma (px)');
    ylabel('Count');
    title('Fitted PSF width (\sigma)');

    subplot(4,1,4);
    histogram(sigma_loc, 100);
    xlabel('Localization precision (px)');
    ylabel('Count');
    title('\sigma_{loc} (precision, CRLB)');


%% Visualize the corrected subpixel locs
    num_show = min(7000, frames);
    frames_to_show = 1:num_show;
    figure;
    viz.plotImage(data, frames_to_show, 'hot','Drift-corrected Image and Locs','max');
    hold on;
    viz.overlayLocs(loc_total,frames_to_show,true);


%% Emitter Analysis of localized points
    bleach_time = 50; %seconds
    bleach_frames = bleach_time/(oneFrameTime); 
    searching_radius = 2; %pixel
    emitters = postproc.emitter.findEmitters(data, loc_total, bleach_frames, searching_radius);
    emitters = postproc.emitter.mergeEmitters(emitters, searching_radius);

%% Filter short/long-lived emitters
    livetime_th = 0.3; %seconds
    emitters_filt = postproc.emitter.filterEmitters_short(emitters, round(livetime_th/ex_time),'consecutive');
    emitters_filt =  postproc.emitter.filterEmitters_firstframe(emitters_filt);
    emitters_filt = postproc.emitter.filterEmitters_end(emitters_filt); 
    [emitters_filt,stats_jump] = postproc.emitter.filterEmitters_jumping(emitters_filt,3.5);
    %emitters_filt = emitters;
    postproc.emitter.plotJumpStats(stats_jump);

%% Check all emitters
    figure;
    viz.plotImage(data, 1:5000, 'gray','Drift-corrected Image and Emitters','max');
    hold on;
    for k = 1:numel(emitters_filt)
        plot(emitters_filt(k).col+0.5, emitters_filt(k).row+0.5, '-', 'LineWidth', 1);
    end
    title('Emitter trajectories overlaid on image');

%% Important Data Collection (position, mean position, survival time)
    stats = postproc.emitter.collectEmitterStatistics(emitters_filt,frames,ex_time,interval,brightness,sigma,sigma_loc,background);
    postproc.emitter.plotEmitterStatistics(stats,90);

%% Check Brightness in First 10 Frames
    F = 10; 
    B = stats.brightness_em(:, 1:F);
    avg_all = mean(B,2,'omitnan')/ex_time;
    figure;
    histogram(avg_all(:),50);   % 对前 F 列所有元素画直方图
    title('Brightness distribution in first few frames');
    xlabel('Brightness(Photon/sce)');
    ylabel('Occurancce');

%% Fit lifetime and brightness
    tau = analysis.photophys.fitLifetime_exp(stats.survival_sec);
    pd = analysis.photophys.fitBrightness_lognormal(stats.brightness_mean);
    n = numel(stats.brightness_mean);
    figure;
    histogram(stats.brightness_mean, 30);
    xlabel('Mean Intensity (photon/sec)');
    ylabel('Counts');
    title(sprintf('Histogram of mean intensities, n=%d', n))
    box on;
    medIntensity = median(stats.brightness_mean);
    fprintf('  median intensity = %.6f photon/sec\n', medIntensity);

%% Test of First Frame Emitters
    ONframe = 1:2;
    pos_ff = postproc.emitter.checkFrameEmitters(stats.pos_matrix, data, ONframe);

%% Extract trace 1
    raw_series = analysis.photophys.extractTrace(stats.pos_mean_px, data(:,:,:));
    bg = prctile(raw_series, 5, 2);
    stats.trace = raw_series - bg;

%% Extract trace 2
    int_range = 4;
    int_area = (2*int_range+1)^2;
    raw_series = analysis.photophys.extractTrace(stats.pos_mean_px, data(:,:,:),int_range);
    % bg = prctile(raw_series, 5, 2);  
    stats.trace = raw_series - stats.bg_mean*int_area;


%% Check single emitter and its trace
    start_frame = 1;
    end_frame = frames;
    for index = 5 %61
        postproc.emitter.checkTrace(start_frame,end_frame,ex_time,stats.trace(index,:),stats.brightness_em(index,:),stats.pos_matrix(:,1,index),stats.pos_matrix(:,2,index));
    end

%% Check emitter video
    postproc.emitter.checkMovie(5, data, emitters_filt, stats.pos_mean_px, stats.pos_matrix, oneFrameTime, 'frameStep',3,'pauseTime', 0.03,'clim', [0 150],'preFrames', 100,'postFrames',300);

%% PhotoAnalysis of Trace
    % states = analysis.photophys.analyzeStates(stats.trace, 'HMM', 'bleachTail',50,'bgThreshold',100,'bicPenalty', 2.0, 'minStateSep', 0);
    % states = analysis.photophys.analyzeStates(stats.trace, 'CHANGEPOINT','penalty',2.0,'minSegLen',4,'mergeThr', 3.0, 'bleachTail',50,'bgThreshold',60, 'minStateSep', 3.0);
    states = analysis.photophys.analyzeStates(stats.trace, 'CHANGEPOINT','penalty',2.5,'minSegLen',4,'mergeThr', 3.3, 'bleachTail',50,'bgThreshold',60, 'minStateSep', 3.3);

%% Visualization
    % analysis.photophys.plotStates(emitters_filt,1:15);
    for index = 1:15
        analysis.photophys.plotState( ...
            states(index), ...
            'useBgSub', true, ...
            'saveFig', false, ...
            'saveDir', '/Users/jerryling/Documents/Uni_graduate/SmilE_matter/DiscussionAndWork/DNA_jinguan/20260411_bazhuayu_trace', ...
            'saveName', sprintf('state_%03d', index) ...
        );
        %close(gcf);
    end

%% Analysis Brightest steps
    N = numel(states);
    brightestMat = nan(N, 8);

    for n = 1:N
        brightestMat(n,1) = n;
        brightestMat(n,6) = states(n).nStates;
        brightestMat(n,7) = states(n).bleachFrame;
        brightestMat(n,8) = states(n).lifetime;
        sInfo = states(n).stateInfo;
        % 跳过没有有效状态的 emitter
        if isempty(sInfo)
            continue;
        end
        % 背景值
        bgVal = 0;
        if isfield(states(n), 'background') && ~isempty(states(n).background) && ~isnan(states(n).background)
            bgVal = states(n).background;
        end
        % 最亮态（背景扣除与否不影响排序，因为只是减同一个常数）
        meanList = [sInfo.meanIntensity];
        [~, idxBright] = max(meanList);
        brightestMat(n,2) = sInfo(idxBright).label;
        brightestMat(n,3) = sInfo(idxBright).meanIntensity - bgVal;  % 扣背景后的 mean
        brightestMat(n,4) = sInfo(idxBright).stdIntensity;
        brightestMat(n,5) = sInfo(idxBright).occupancy;
    end

%% 显示结果矩阵
    disp('Columns: [emitterIdx, brightestLabel, brightestMean, brightestStd, brightestOccupancy, nStates, bleachFrame, lifetime]');
    disp(brightestMat);

%% 提取最亮态平均强度，去掉 NaN
    brightestIntensity = brightestMat(:,3);
    brightestIntensity = brightestIntensity(~isnan(brightestIntensity));

%% 作直方图
    figure;
    histogram(brightestIntensity/ex_time, 30);
    xlabel('Brightest state mean intensity(photon/sec)');
    ylabel('Counts');
    title('Histogram of brightest state intensity across emitters');
    box on;

%% 寿命
    lifetimes = brightestMat(:,8)*ex_time;
    lifetimes = lifetimes(~isnan(lifetimes));
    tau = analysis.photophys.fitLifetime_exp(lifetimes);

%% 所有态亮度分布
    allStateIntensity = [];

    for n = 1:numel(states)
        sInfo = states(n).stateInfo;
        if isempty(sInfo)
            continue;
        end

        % 背景值
        bgVal = 0;
        if isfield(states(n), 'background') && ~isempty(states(n).background) && ~isnan(states(n).background)
            bgVal = states(n).background;
        end

        % 只保留非背景态
        labels = [sInfo.label];
        validIdx = labels > 0;

        if any(validIdx)
            meanVals = [sInfo(validIdx).meanIntensity] - bgVal;
            allStateIntensity = [allStateIntensity; meanVals(:)];
        end
    end

    figure;
    histogram(allStateIntensity / ex_time, 35);
    xlabel('State mean intensity (photon/sec)');
    ylabel('Counts');
    title('Histogram of all non-background state intensities');
    box on;