%% ======================================================================
    %  pipeline_main.m — 单分子定位分析主流程
    %  四个模块: detectAndLocalize → correctDrift → analyzeEmitters → analyzeTraces
    % =======================================================================
    clear; clc;
    parameters = struct();

%% ===================== Parameters =====================

    % --- 文件与 I/O ---
    parameters.file_name = ...
        '/Volumes/SMILeSSD/Optics/TE2000UTests/Cy3/20260603/532nm_1x_0p1S_1000frames_1p00mW_cy3.sif';
    parameters.interval = 0;

    % --- Windowing ---
    parameters.row_range = 100:412;
    parameters.col_range = 100:412;

    % --- Detection / Localization ---
    parameters.k_sigma        = 2.0;
    parameters.edge           = 5;
    parameters.viz_enabled    = true;
    parameters.viz_max_frames = 5000;

    % --- Drift Estimation ---
    parameters.drift_corr           = true;
    parameters.drift_frames_per_seg = 15;
    parameters.drift_min_locs       = 30;

    % --- Drift Correction ---
    parameters.row_width = 227;
    parameters.col_width = 454;

    % --- Emitter Analysis ---
    parameters.bleach_time      = 50; % beyond this time, undetected emitters will be marked as bleached. Note that accurate detection is necessary. If highly blinking molecules exists, one may set large bleach_time, then all locs will be added into searching set.
    parameters.searching_radius = 2;

    % --- Emitter Filtering (nan = skip) ---
    parameters.livetime_th      = 0.2; % Filter short-lived emitters
    parameters.filter_firstframe = nan; % Filter non-first-frame emitters
    parameters.filter_end = true; % Filter unbleached emitters
    parameters.jump_threshold   = 3.5; % Filter jumping/unstable emitters 

    % --- Trace Extraction ---
    parameters.int_range = 4;
    parameters.bg_range  = 7;
    parameters.bg_extracted = true;
    parameters.bg_extraction_method = 'median';

    % --- State Analysis ---
    parameters.emittersFIT_enabled = false;
    parameters.state_method      = 'CHANGEPOINT';
    parameters.state_penalty     = 4.0;  % penalty of adding steps
    parameters.state_min_seg_len = 4; % minium length of a step, in frame units
    parameters.state_merge_thr   = 2.5; % merging threthold in noise-std units

    parameters.state_bicPenalty = 2.0; % HMM BIC penalty multiplier

    parameters.state_bleach_tail = 50; % considered bleached
    parameters.state_bg_threshold = 60; % bg threshold, under which is set as bg
    parameters.state_min_state_sep = 3.0; % a second merge, for HMM 0 is recommanded
    



%% ===================== Read & Window =====================

    [raw_data, ex_time, ~] = io.readSIFData(parameters.file_name);
    parameters.ex_time        = ex_time;          % ← 写入 parameters，全程透传
    parameters.one_frame_time = ex_time + parameters.interval;

    windowed_raw_data = raw_data(parameters.row_range, parameters.col_range, :);
    parameters.frames = size(windowed_raw_data, 3);
    clear raw_data;

    fprintf('Windowed data: %d × %d × %d\n', ...
        size(windowed_raw_data,1), size(windowed_raw_data,2), parameters.frames);

    % Export Uncorrected TIFF (optional)
    % out_tif = [parameters.file_name(1:end-4), '_Tifstack.tif'];
    % io.writeTiffStack(windowed_raw_data, out_tif);

%% ===================== Module 1: Detect + Localize =====================

    [file_dir, file_base, ~] = fileparts(parameters.file_name);
    cacheFile = fullfile(file_dir, [file_base, '_locResult.mat']);

    locResult = pipeline.detectAndLocalize(windowed_raw_data, parameters, ...
        'cacheFile', cacheFIle);

%% ===================== Module 2: Drift Correction =====================

    [data, driftResult] = pipeline.correctDrift( ...
        windowed_raw_data, locResult, parameters);
    % clear windowed_raw_data;

    % Export Corrected TIFF (optional)
    % out_corr_tif = [parameters.file_name(1:end-4), '_Tifstack_corrected.tif'];
    % io.writeTiffStack(data, out_corr_tif);

%% ===================== Module 3: Emitter Analysis =====================

    emitterResult = pipeline.analyzeEmitters(data, driftResult, parameters);

    % 说明
    % 输出中包含最重要两个结构体：
    % emitterResult.emitters_filt(k) 代表第k个经过筛选的发光点，包含长度为定位次数的数列：行坐标.row、列坐标.col、帧数.frames、在所有定位的矩阵loc_total中的行数.loc_idx、首次detect到的帧数.on_frame、漂白的帧数.bleach_frame、存活状态.alive。该变量以emitter为中心储存信息，可兼容tracking追踪。
    % emitterResult.stats 包含一系列重要的统计数据，大部份以矩阵或单值储存，以方便后续可视化/进一步处理。其中.pos_matrix储存了所有的emitters的位置信息，未检测到的帧以Nan记录、.pos_mean_px计算了emitter平均位置用于可视化、.survival_sec计算了基于detect和emitter分析的存活时间、.brightness_em/mean/sum计算了拟合亮度信息、.bg_mean为拟合的平均背景、.sigma_mean为拟合的平均标准差、.sigma_loc_mean则为估计的定位精度。

%% ===================== Module 4: Trace & State Analysis =====================

    traceResult = pipeline.analyzeTraces(data, emitterResult, parameters);

%% ===================== Save All Results =====================

    save_path = fullfile(file_dir, [file_base, '_pipeline_result.mat']);
    save(save_path, 'parameters', 'locResult', 'driftResult', ...
        'emitterResult', 'traceResult', '-v7.3');
    fprintf('All results saved: %s\n', save_path);

%% ======================================================================
%  以下为交互式检查，按需取消注释运行
% =======================================================================


%% [Emitters] Check some
    ON_frames = 1:2;
    pos_ff = postproc.emitter.checkFrameEmitters( ...
        emitterResult.stats.pos_matrix, data, ON_frames);

%% [Emitters] Check all
    figure;
    viz.plotImage(data, 1:parameters.frames, 'gray', ...
        'Drift-corrected Image and Emitters', 'max');
    hold on;
    for k = 1:numel(emitterResult.emitters_filt)
        plot(emitterResult.emitters_filt(k).col + 0.5, ...
            emitterResult.emitters_filt(k).row + 0.5, '-', 'LineWidth', 1);
    end
    title('Emitter trajectories overlaid on image');

%% [Emitters] Check Local Movie
    postproc.emitter.checkMovie(5, data, emitterResult.emitters_filt, ...
        emitterResult.stats.pos_mean_px, emitterResult.stats.pos_matrix, ...
        parameters.one_frame_time, ...
        'frameStep', 3, 'pauseTime', 0.005, 'clim', [0 150], 'preFrames', 100);

%% [Trace] Inspect Trace
    for idx = 1:10
        postproc.emitter.checkTrace( ...
            1, parameters.frames, parameters.ex_time, ...
            traceResult.trace(idx, :), ...
            emitterResult.stats.brightness_em(idx, :), ...
            emitterResult.stats.pos_matrix(:, 1, idx), ...
            emitterResult.stats.pos_matrix(:, 2, idx));
    end

%% [Trace] Check States
    for idx = 11:30
        analysis.photophys.plotState(traceResult.states(idx));
    end

%% Display Brightest State Matrix
% disp('Columns: [idx, label, mean, std, occ, nStates, bleachFr, lifetime]');
% disp(traceResult.brightestMat);