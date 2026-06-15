%% 初始化+重要中间变量说明 
    clear all; clc;
    parameters = struct();
% [uncorrected_detected_total]
%   3列矩阵 [row, col, frame]，detect 阶段的粗定位结果，后续送入 localization。
%
% [uncorrected_super_loc_total]
%   7列矩阵，MLE localization 的完整输出。后续用于 drift correction。
%   该步骤耗时较长，建议保存结果以避免重复计算。若需手动构造，必须遵循7列格式。
%
% [delta_sum]
%   drift estimation 的输出，描述各帧相对于首帧的累积漂移。
%
% [loc_total, brightness, background, sigma, sigma_loc]
%   drift correction 后，corrected_super_loc_total 的各列拆分。
%
% [emitters_filt]
%   结构体，包含筛选后合格 emitter 的全部信息。
%
% [stats]
%   在 emitters_filt 基础上汇总的统计信息：平均位置、存活时间、亮度等。
%
% [series]
%   背景提取后的强度 trace，行索引与 emitters_filt / stats 一一对应。

%% Parameters Setting
    % 所有可调参数统一放入 parameters 结构体，方便后续封装为 pipeline 函数。
    % 使用时通过 parameters.xxx 访问；传入子函数只需传一个结构体即可。


    % --- 文件与 I/O ---
    parameters.file_name = ...
        '/Volumes/SMILeSSD/Optics/Ti2Tests/TirfModule/Cy3/20260603/532_cy3_0.1s_1.79A/TUC-001.tif';
    %parameters.file_name = ...
      %  '/Volumes/SMILeSSD/Optics/Ti2Tests/TirfModule/Beads_532/20260514_532_beads_0.1s_1.59A/TUC-001.tif';
    parameters.interval = 0; % 相邻帧之间的死时间（秒）。曝光时间从 .sif 文件头读取，总帧间隔 = exposure + interval
    parameters.gain = 33.448; % DN/e-
    parameters.ex_time = 0.1;
    parameters.offset = 1893.68;

    % --- Windowing 参数 --- 
    parameters.row_range = 624:1424;
    parameters.col_range = 624:1424; % 裁切的行/列范围。应略大于最终裁剪的范围（如257*257）以保证 drift correction 有足够重叠，且尽量居中

    % --- Detection/Localization 参数 ---
    psf_estimated = 1.2; % pixel size
    parameters.k_sigma = 6.5; % local maximum 检测阈值：像素值需超过局部背景 + k_sigma × 局部标准差才被标记为候选点
    parameters.edge = ceil(psf_estimated*3.3-0.5); % 边缘防护宽度（像素）。距图像边缘 < edge 的 detection 将被跳过，避免 PSF 截断导致拟合失败
    parameters.viz_enabled = true; % 是否绘制 detection / localization 的叠加图
    parameters.viz_max_frames = 5000; % 可视化时最多显示的帧数（避免帧数过多导致绘图缓慢）


    % --- Drift Estimation 参数 ---
    parameters.drift_corr = true;           % 是否执行漂移校正；false 则跳过，delta_sum 置零
    parameters.drift_frames_per_seg = 15;   % 每段帧数
    parameters.drift_min_locs = 30;         % 每段最少定位数

    % --- Drift Correction 参数 ---
    parameters.row_width = 501;
    parameters.col_width = 501;

    % --- Emitter Analysis 参数 ---
    parameters.bleach_time = 4;            % 光漂白搜索窗口（秒），决定 emitter 聚类的时间跨度
    parameters.searching_radius = 2;        % emitter 聚类搜索半径（像素）
    parameters.livetime_th = 0.2;           % 最短存活时间阈值（秒），低于此值的 emitter 被过滤
    parameters.jump_threshold = 3.5;        % 位置跳跃过滤阈值（σ 倍数）

    % --- Trace Extraction 参数 ---
    parameters.int_range = ceil(psf_estimated*2.5-0.5);               % 积分半径（像素），积分面积 = (2*int_range+1)^2
    parameters.bg_range = parameters.int_range+3;                % 扣除背景积分半径（像素）

    % --- State Analysis 参数 ---
    parameters.state_method = 'CHANGEPOINT';
    parameters.state_penalty = 4.0;         % changepoint 惩罚系数
    parameters.state_min_seg_len = 4;       % 最短 segment 长度（帧）
    parameters.state_merge_thr = 3.0;       % 状态合并阈值
    parameters.state_bleach_tail = 50;      % 漂白尾部截断（帧）
    parameters.state_bg_threshold = 130;     % 背景态判定阈值
    parameters.state_min_state_sep = 3.0;   % 最小态间距

%% Read Raw Data
    raw_data = io.readTiffStack(parameters.file_name);
    raw_data = double(raw_data) - parameters.offset;
    raw_data = raw_data / parameters.gain;
    parameters.one_frame_time = parameters.ex_time + parameters.interval; % 计算单帧总时长（曝光 + 死时间），供后续使用

%% Windowing
    windowed_raw_data = raw_data(parameters.row_range, parameters.col_range, 1:100);
    parameters.frames = size(windowed_raw_data, 3);   % 记录总帧数，后续多处复用
    fprintf('Windowed Data: %d × %d × %d\n', ...
        size(windowed_raw_data,1), size(windowed_raw_data,2), parameters.frames);

%% Release raw data
    clear raw_data;

%% Export TIFF stack（optional）
    % 将裁切后的数据保存为 .tif，便于在 ImageJ / Fiji 等软件中浏览
    out_tif = [parameters.file_name(1:end-4), '_Tifstack.tif'];
    io.writeTiffStack(windowed_raw_data, out_tif);
    fprintf('TIFF stack stored: %s\n', out_tif);

%% Viz
    figure;
    viz.plotImage(windowed_raw_data, 1:10, 'hot','Drift-uncorrected Image and Detections');
    hold on;

%% Hot Pixel Detection for Non-bleaching Particles
    %  Key idea: hot pixel is single-pixel sharp; real particle is PSF-broadened
    %  Insert after Windowing, before MLE / trace extraction

    % --- Parameters ---
    hot_pixel_params.spatial_kernel = 5;       % neighborhood size for local comparison
    hot_pixel_params.sharpness_thresh = 0.2;   % fraction of energy in center pixel vs ROI
    hot_pixel_params.persist_frac = 0.95;      % fraction of frames pixel must be "sharp" to be hot
    hot_pixel_params.intensity_percentile = 90; % only examine bright pixels (avoid noise floor)

    [nRows, nCols, nFrames] = size(windowed_raw_data);

    % --- Step 1: Time-averaged image to find candidate bright pixels ---
    avg_img = mean(windowed_raw_data, 3);

    % Only consider pixels above a brightness threshold (hot pixels & particles are both bright)
    intensity_thresh = prctile(avg_img(:), hot_pixel_params.intensity_percentile);
    bright_mask = avg_img > intensity_thresh;

    % --- Step 2: Sharpness ratio per pixel per frame ---
    % For each pixel, compute: (center pixel value) / (sum of local 3x3 patch)
    % Hot pixel: ratio ~ 1 (all signal in one pixel)
    % Real PSF:  ratio ~ 0.15-0.4 (signal spread over neighbors)

    half_k = floor(hot_pixel_params.spatial_kernel / 2);
    pad_data = padarray(windowed_raw_data, [half_k half_k 0], 'symmetric');

    % Compute local sum using convolution (fast)
    local_kernel = ones(hot_pixel_params.spatial_kernel);
    kernel_area = hot_pixel_params.spatial_kernel^2;

    % Preallocate sharpness count
    sharp_count = zeros(nRows, nCols);

    % Process in batches to save memory
    batch_size = min(100, nFrames);
    for f_start = 1:batch_size:nFrames
        f_end = min(f_start + batch_size - 1, nFrames);
        batch = windowed_raw_data(:,:,f_start:f_end);
        
        for f_idx = 1:size(batch, 3)
            frame = batch(:,:,f_idx);
            % Local sum via conv2
            local_sum = conv2(frame, local_kernel, 'same');
            % Avoid division by zero
            local_sum(local_sum <= 0) = 1;
            % Sharpness: center value / local sum
            sharpness = frame ./ local_sum;
            % Count frames where pixel is "sharp" AND bright
            is_sharp = (sharpness > hot_pixel_params.sharpness_thresh) & bright_mask;
            sharp_count = sharp_count + double(is_sharp);
        end
    end

    % --- Step 3: Hot pixel = persistently sharp across most frames ---
    persistence = sharp_count / nFrames;
    hot_pixel_mask = persistence > hot_pixel_params.persist_frac;

    % --- Step 4: Safety check - exclude anything that looks like a PSF ---
    % Dilate hot pixel candidates and check if neighbors are also bright
    % (if yes, it's likely a real particle, not a hot pixel)
    se = strel('square', 3);
    neighbor_bright = imdilate(bright_mask, se) & ~bright_mask; % ring of neighbors

    % For each hot pixel candidate, check if its neighbors are also significantly bright
    hot_pixel_mask_refined = hot_pixel_mask;
    [hot_r, hot_c] = find(hot_pixel_mask);
    for idx = 1:length(hot_r)
        r = hot_r(idx); c = hot_c(idx);
        r_min = max(1, r-1); r_max = min(nRows, r+1);
        c_min = max(1, c-1); c_max = min(nCols, c+1);
        patch = avg_img(r_min:r_max, c_min:c_max);
        center_val = avg_img(r, c);
        neighbor_vals = patch(:);
        neighbor_vals(neighbor_vals == center_val) = []; % remove center
        % If neighbors have significant signal (>30% of center), likely a real particle
        if mean(neighbor_vals) > 0.3 * center_val
            hot_pixel_mask_refined(r, c) = false;
        end
    end

    hot_pixel_mask = hot_pixel_mask_refined;

    num_hot = sum(hot_pixel_mask(:));
    fprintf('Hot pixels detected: %d / %d (%.4f%%)\n', num_hot,nRows*nCols, 100*num_hot/numel(hot_pixel_mask));

    % --- Step 5: Replace hot pixels frame-by-frame ---
    corrected_data = windowed_raw_data;
    k = hot_pixel_params.spatial_kernel;
    for f = 1:nFrames
        frame = windowed_raw_data(:,:,f);
        frame_filtered = medfilt2(frame, [k k]);
        frame(hot_pixel_mask) = frame_filtered(hot_pixel_mask);
        corrected_data(:,:,f) = frame;
    end
    windowed_raw_data = corrected_data;
    fprintf('Hot pixels replaced with %d×%d spatial median.\n', k, k);

    % --- Step 6: Diagnostic figure ---
    figure('Name', 'Hot Pixel Diagnostics (Non-bleaching)', 'Position', [100 100 1400 400]);

    subplot(1,4,1);
    imagesc(avg_img); axis image; colorbar;
    title('Time-averaged image');

    subplot(1,4,2);
    imagesc(persistence); axis image; colorbar;
    title('Sharpness persistence map');

    subplot(1,4,3);
    imagesc(hot_pixel_mask); axis image;
    title(sprintf('Hot Pixel Mask (N=%d)', num_hot));
    colormap(gca, [0 0 0; 1 0 0]);

    subplot(1,4,4);
    corrected_avg = mean(windowed_raw_data, 3);
    imagesc(corrected_avg); axis image; colorbar;
    title('After correction');

%% Detection (Roughly find locs)
    uncorrected_detected_total = detect.findMaxima(windowed_raw_data, parameters.k_sigma, parameters.edge); % 输出: [row, col, frame] 三列矩阵

%% Visualize Detection Results
    if parameters.viz_enabled
        frames_to_show = 1:min(1, parameters.frames);

        figure;
        viz.plotImage(windowed_raw_data, frames_to_show, 'hot', ...
            'Drift-uncorrected Image and Detections');
        hold on;
        viz.overlayLocs(uncorrected_detected_total, frames_to_show);
    end

%% Localization
    % 注：此步骤计算量大，建议完成后保存结果
    uncorrected_super_loc_total_raw = ...
        localize.MLE1.locMolecules(windowed_raw_data, uncorrected_detected_total, parameters.edge,1); % 输出 7 列: [row, col, brightness, background, sigma, sigma_loc, frame], sigma_loc 即 CRLB 估计的定位不确定度（像素）

    % 过滤不合理的 localization（如负亮度、sigma 异常、位置跑出 ROI 等）
    uncorrected_super_loc_total = localize.filterBadLocs(uncorrected_super_loc_total_raw);

    % 提取精简三列 [row, col, frame]，供可视化和 drift correction 使用
    uncorrected_loc_total = uncorrected_super_loc_total(:, [1, 2, 7]);

    fprintf('Localization Complete: %d → %d（%d bad locs filtered）\n', ...
        size(uncorrected_super_loc_total_raw, 1), ...
        size(uncorrected_super_loc_total, 1), ...
        size(uncorrected_super_loc_total_raw, 1) - size(uncorrected_super_loc_total, 1));

%% Visualize Localization Results
    if parameters.viz_enabled
        frames_to_show = 1:min(parameters.viz_max_frames, parameters.frames);

        figure;
        viz.plotImage(windowed_raw_data, frames_to_show, 'hot', ...
            'Drift-uncorrected Image and Locs', 'max');
        hold on;
        viz.overlayLocs(uncorrected_loc_total, frames_to_show, true);
    end

%% Save Localization Results
    % MLE localization 耗时最长，保存后可跳过重复计算直接加载
    [file_dir, file_base, ~] = fileparts(parameters.file_name);
    save_name = fullfile(file_dir, [file_base, '_uncorrected_super_loc_total.mat']);
    save(save_name, 'uncorrected_super_loc_total', '-v7.3');
    fprintf('Localization result saved: %s\n', save_name);

%% Load Localization Results
    % 若此前已保存，可直接从此 section 开始运行，跳过 detection + localization
    [file_dir, file_base, ~] = fileparts(parameters.file_name);
    load_name = fullfile(file_dir, [file_base, '_uncorrected_super_loc_total.mat']);
    S = load(load_name);
    uncorrected_super_loc_total = S.uncorrected_super_loc_total;
    uncorrected_loc_total = uncorrected_super_loc_total(:, [1, 2, 7]);
    fprintf('Localization result loaded: %s\n', load_name);

%% Drift Estimation (NP-cloud)
    % 用 localization 的中位定位精度作为 kernel bandwidth
    uncertainty = median(uncorrected_super_loc_total(:,6));

    if parameters.drift_corr
        delta_sum = postproc.drift.estimateDrift_cloudxy( ...
            windowed_raw_data, uncorrected_loc_total, ...
            parameters.drift_frames_per_seg, uncertainty, ...
            parameters.drift_min_locs, false);
        % uncorrected_loc_total 仅含 [row, col, frame] 三列
    else
        delta_sum = zeros(parameters.frames, 2);
        fprintf('Drift correction disabled, delta_sum set to zero.\n');
    end

%% Correct Drift
    % 根据 delta_sum 对原始数据和定位进行漂移校正，同时裁切到目标尺寸
    [data, loc_total, loc_idx] = postproc.drift.correctDrift( ...
        windowed_raw_data, uncorrected_loc_total, delta_sum, ...
        parameters.edge, parameters.row_width, parameters.col_width);

    % loc_idx 为保留下来的定位索引，用于同步提取完整 7 列信息
    corrected_super_loc_total = uncorrected_super_loc_total(loc_idx, :);

%% Visualize Drift
    viz.plotTracking(delta_sum, parameters.one_frame_time, [], 'Drift Tracking');

%% Save Drift-corrected TIFF Stack (optional)
    out_corr_tif = [parameters.file_name(1:end-4), '_Tifstack_corrected.tif'];
    io.writeTiffStack(data, out_corr_tif);
    fprintf('Corrected TIFF stack saved: %s\n', out_corr_tif);

%% Load Drift-corrected TIFF Stack (optional)
    % 若此前已保存校正后数据，可直接从此处开始运行
    in_corr_tif = [parameters.file_name(1:end-4), '_Tifstack_corrected.tif'];
    data = io.readTiffStack(in_corr_tif);
    fprintf('Corrected TIFF stack loaded: %s\n', in_corr_tif);

%% Unpack MLE Columns
    % 将校正后的 7 列 localization 拆分为独立变量，便于后续引用
    brightness = corrected_super_loc_total(:,3);   % 光子数
    background = corrected_super_loc_total(:,4);   % 背景
    sigma      = corrected_super_loc_total(:,5);   % PSF 宽度（像素）
    sigma_loc  = corrected_super_loc_total(:,6);   % 定位精度 / CRLB（像素）

%% Visualize MLE Distributions
    figure;
    subplot(4,1,1);
    histogram(brightness, 180);
    xlabel('Photons (N)'); ylabel('Count');
    title('Brightness distribution (by locs)');

    subplot(4,1,2);
    histogram(background, 140);
    xlabel('Photons (N)'); ylabel('Count');
    title('Background distribution (by locs)');

    subplot(4,1,3);
    histogram(sigma, 120);
    xlabel('Sigma (px)'); ylabel('Count');
    title('Fitted PSF width (\sigma)');

    subplot(4,1,4);
    histogram(sigma_loc, 100);
    xlabel('Localization precision (px)'); ylabel('Count');
    title('\sigma_{loc} (precision, CRLB)');

%% Visualize Corrected Localizations
    if parameters.viz_enabled
        frames_to_show = 1:min(parameters.viz_max_frames, parameters.frames);

        figure;
        viz.plotImage(data, frames_to_show, 'hot', ...
            'Drift-corrected Image and Locs', 'max');
        hold on;
        viz.overlayLocs(loc_total, frames_to_show, true);
    end

%% Emitter Analysis
    % 将逐帧定位按空间和时间邻近关系聚类为 emitter
    bleach_frames = parameters.bleach_time / parameters.one_frame_time;
    emitters = postproc.emitter.findEmitters(data, loc_total, ...
        bleach_frames, parameters.searching_radius);
    emitters = postproc.emitter.mergeEmitters(emitters, parameters.searching_radius);

%% Filter Emitters
    % 依次过滤：短寿命 → 末帧残留 → 位置跳跃异常
    min_frames = round(parameters.livetime_th / parameters.ex_time);
    emitters_filt = emitters;
    emitters_filt = postproc.emitter.filterEmitters_short(emitters_filt, min_frames, 'consecutive');
    emitters_filt = postproc.emitter.filterEmitters_firstframe(emitters_filt);
    %emitters_filt = postproc.emitter.filterEmitters_end(emitters_filt);
    %[emitters_filt, stats_jump] = postproc.emitter.filterEmitters_jumping( ...
    %    emitters_filt, parameters.jump_threshold);
    %postproc.emitter.plotJumpStats(stats_jump);

    fprintf('Emitters after filtering: %d / %d\n', numel(emitters_filt), numel(emitters));

%% Visualize All Emitter Trajectories
    figure;
    viz.plotImage(data, 1:parameters.frames, 'gray', ...
        'Drift-corrected Image and Emitters', 'max');
    hold on;
    for k = 1:numel(emitters_filt)
        plot(emitters_filt(k).col + 0.5, emitters_filt(k).row + 0.5, '-', 'LineWidth', 1);
    end
    title('Emitter trajectories overlaid on image');

%% Collect Emitter Statistics
    % 汇总每个 emitter 的平均位置、存活时间、亮度等
    stats = postproc.emitter.collectEmitterStatistics( ...
        emitters_filt, parameters.frames, parameters.ex_time, parameters.interval, ...
        brightness, sigma, sigma_loc, background);
    postproc.emitter.plotEmitterStatistics(stats, 50);

%% Check Brightness in First Few Frames
    F = 10;
    B = stats.brightness_em(:, 1:F);
    avg_all = mean(B, 2, 'omitnan') / parameters.ex_time;

    figure;
    histogram(avg_all(:), 50);
    title('Brightness distribution in first few frames');
    xlabel('Brightness (photon/s)');
    ylabel('Occurrence');

%% Fit Lifetime and Brightness (Emitter Analysis Results)
    tau = analysis.photophys.fitLifetime_exp(stats.survival_sec);
    pd  = analysis.photophys.fitBrightness_lognormal(stats.brightness_mean);
    % 追踪寿命转化为计数存活曲线
    t_max = max(stats.survival_sec);
    t_axis = 0:t_max;
    N_surviving = arrayfun(@(t) sum(stats.survival_sec > t), t_axis);

    % 拟合单指数衰减
    f = fit(t_axis(:), N_surviving(:), 'exp1');  % N0 * exp(-k*t)
    tau_counting = -1 / f.b;  % 特征寿命

    %可视化：计数存活曲线
    figure('Name', 'Counting Survival Curve', 'Color', 'w');

    % 主图：存活曲线 + 拟合
    subplot(2,1,1);
    hold on;
    stairs(t_axis, N_surviving, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Data');
    plot(t_axis, f(t_axis), 'r-', 'LineWidth', 2, 'DisplayName', ...
        sprintf('Exp fit: \\tau = %.2f s', tau_counting));
    hold off;
    xlabel('Time (s)');
    ylabel('N surviving');
    title(sprintf('Counting Survival Curve (N_0 = %d, \\tau = %.2f s)', ...
        N_surviving(1), tau_counting));
    legend('Location', 'northeast');
    set(gca, 'FontSize', 12, 'Box', 'on');

    % 半对数图：检验单指数性
    subplot(2,1,2);
    hold on;
    stairs(t_axis, N_surviving, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Data');
    plot(t_axis, f(t_axis), 'r-', 'LineWidth', 2, 'DisplayName', 'Exp fit');
    hold off;
    set(gca, 'YScale', 'log');
    xlabel('Time (s)');
    ylabel('N surviving (log)');
    title('Semi-log: single exponential \rightarrow straight line');
    legend('Location', 'northeast');
    set(gca, 'FontSize', 12, 'Box', 'on');
    ylim([max(1, min(N_surviving(N_surviving>0))), N_surviving(1)*1.1]);

    sgtitle(sprintf('Tracking \\tau = %.2f s | Counting \\tau = %.2f s', tau, tau_counting), ...
        'FontSize', 14, 'FontWeight', 'bold');



%% Check First-frame Emitters
    ON_frames = 1:2;
    pos_ff = postproc.emitter.checkFrameEmitters(stats.pos_matrix, data, ON_frames);

%% Extract Trace 
    [stats.trace,stats.snr ,~,~] = analysis.photophys.extractTrace(stats.pos_mean_px, data(:,:,:),parameters.int_range,true, parameters.bg_range,'median'); 

%% Inspect Single Emitter Traces
    start_frame = 1;
    end_frame   = parameters.frames;
    for index = 70:86
        postproc.emitter.checkTrace(start_frame, end_frame, parameters.ex_time, ...
            stats.trace(index,:), stats.brightness_em(index,:), ...
            stats.pos_matrix(:,1,index), stats.pos_matrix(:,2,index));
    end

%% Inspect Emitter Movie
    postproc.emitter.checkMovie(5, data, emitters_filt, ...
        stats.pos_mean_px, stats.pos_matrix, parameters.one_frame_time, ...
        'frameStep', 3, 'pauseTime', 0.005, 'clim', [0 130], 'preFrames', 100);

%% State Analysis (Changepoint)
    states = analysis.photophys.analyzeStates(stats.trace, ...
        parameters.state_method, ...
        'penalty',      parameters.state_penalty, ...
        'minSegLen',    parameters.state_min_seg_len, ...
        'mergeThr',     parameters.state_merge_thr, ...
        'bleachTail',   parameters.state_bleach_tail, ...
        'bgThreshold',  parameters.state_bg_threshold, ...
        'minStateSep',  parameters.state_min_state_sep);

%% Visualize States
    for index = 70:86
        analysis.photophys.plotState(states(index));
    end
%% Brightness
    for index = 70:86
        fprintf('%d\n',mean(states(index).trace)/1000);
    end

%% SNR calculation
    % 假设 stats.trace 和 stats.snr 已经通过 extractTrace 提取
    % [stats.trace, stats.snr, ~, ~] = analysis.photophys.extractTrace(stats.pos_mean_px, data(:,:,:), parameters.int_range, true, parameters.bg_range, 'median'); 

    % 假设 states 是 analyzeStates 函数的输出结果，且与 stats.trace 行数一一对应
    % states = analysis.photophys.analyzeStates(...)

    % 1. 筛选出只有 1 个 step（即 nStates == 1）的合规 traces
    isValidTrace = [states.nStates] == 1; 
    valid_indices = find(isValidTrace);
    num_valid = length(valid_indices);

    % 初始化 Temporal SNR 数组 (每个 trace 1 个值)
    temporal_snr_values = zeros(num_valid, 1);

    % 初始化一个 cell 数组用于收集所有 trace 的 ON 态 Spatial SNR (每帧 1 个值)
    all_spatial_snr_cells = cell(num_valid, 1);

    % 2. 遍历所有合规 trace 计算两种 SNR
    for i = 1:num_valid
        orig_idx = valid_indices(i); % 获取该分子在原始矩阵中的行号
        s = states(orig_idx);
        
        % ==========================================
        % 计算 Temporal SNR (时间信噪比) - 每个分子1个值
        % ==========================================
        Signal_emitter = s.stateInfo(1).meanIntensity;
        sigma_emitter  = s.stateInfo(1).stdIntensity;
        
        bg_mask = (s.sequence == 0);
    
        Signal_bg = mean(s.trace(bg_mask));
        sigma_bg  = std(s.trace(bg_mask));
    
        denominator = sqrt(sigma_emitter^2 + sigma_bg^2);
        if denominator < 1e-6
            denominator = 1e-6;
        end
        
        temporal_snr_values(i) = (Signal_emitter - Signal_bg) / denominator;
        
        % ==========================================
        % 收集 Spatial SNR (空间信噪比) - 收集所有 ON 态帧
        % ==========================================
        current_spatial_snr_trace = stats.snr(orig_idx, :);
        
        % ON 态的 mask (假设非 0 即为有信号的 ON 态)
        on_mask = (s.sequence > 0);
        
        if sum(on_mask) > 0
            % 将该分子所有 ON 态帧的 Spatial SNR 保存为列向量
            all_spatial_snr_cells{i} = current_spatial_snr_trace(on_mask)';
        end
    end

    % ========== 调试：揪出 Spatial SNR < 5 的 trace ==========
    snr_threshold = 5;

    fprintf('\n============================================\n');
    fprintf('  DEBUG: Traces with median Spatial SNR < %.1f\n', snr_threshold);
    fprintf('============================================\n');

    low_snr_count = 0;

    for i = 1:num_valid
        if isempty(all_spatial_snr_cells{i})
            continue;
        end
        
        median_spatial_snr = median(all_spatial_snr_cells{i});
        
        if median_spatial_snr < snr_threshold
            low_snr_count = low_snr_count + 1;
            orig_idx = valid_indices(i);
            s = states(orig_idx);
            
            % 基本信息
            pos_r = round(stats.pos_mean_px(orig_idx, 1));
            pos_c = round(stats.pos_mean_px(orig_idx, 2));
            
            % ON 态帧数
            on_mask = (s.sequence > 0);
            n_on = sum(on_mask);
            
            % 背景统计
            bg_vals = stats.bg_mean(orig_idx);
            mean_bg = mean(bg_vals);
            std_bg  = std(bg_vals);
            
            % 信号统计
            raw_on = stats.trace(orig_idx, on_mask);
            trace_on = stats.trace(orig_idx, on_mask);
            
            % Spatial SNR 分布
            spsnr_on = all_spatial_snr_cells{i};
            
            fprintf('\n--- [%d] orig_idx=%d, pos=(%d,%d) ---\n', low_snr_count, orig_idx, pos_r, pos_c);
            fprintf('  ON frames:        %d / %d total\n', n_on, length(s.sequence));
            fprintf('  Temporal SNR:     %.2f\n', temporal_snr_values(i));
            fprintf('  Spatial SNR:      median=%.2f, mean=%.2f, min=%.2f, max=%.2f\n', ...
                median(spsnr_on), mean(spsnr_on), min(spsnr_on), max(spsnr_on));
            fprintf('  Raw intensity:    mean=%.1f, std=%.1f\n', mean(raw_on), std(raw_on));
            fprintf('  BG-sub signal:    mean=%.1f, std=%.1f\n', mean(trace_on), std(trace_on));
            fprintf('  Background/px:    mean=%.2f, std=%.2f\n', mean_bg, std_bg);
            fprintf('  State intensity:  mean=%.1f, std=%.1f\n', s.stateInfo(1).meanIntensity, s.stateInfo(1).stdIntensity);
            
            % 检查是否靠近边缘
            [imgR, imgC, ~] = size(data);
            edge_dist = min([pos_r-1, imgR-pos_r, pos_c-1, imgC-pos_c]);
            fprintf('  Distance to edge: %d px\n', edge_dist);
            
            % 检查邻近分子数量
            all_pos = stats.pos_mean_px;
            dists = sqrt((all_pos(:,1) - pos_r).^2 + (all_pos(:,2) - pos_c).^2);
            dists(orig_idx) = Inf; % 排除自身
            n_neighbors_close = sum(dists < 10);
            nearest_dist = min(dists);
            fprintf('  Neighbors <10px:  %d (nearest=%.1f px)\n', n_neighbors_close, nearest_dist);
            
            % 只详细打印前 20 个，避免刷屏
            if low_snr_count >= 20
                fprintf('\n  ... (showing first 20 only, more exist)\n');
                break;
            end
        end
    end

    % 汇总统计
    median_all = cellfun(@(x) median(x), all_spatial_snr_cells(~cellfun(@isempty, all_spatial_snr_cells)));
    n_below = sum(median_all < snr_threshold);

    fprintf('\n============================================\n');
    fprintf('  SUMMARY\n');
    fprintf('  Total valid traces:       %d\n', num_valid);
    fprintf('  Median Spatial SNR < %.1f: %d  (%.1f%%)\n', snr_threshold, n_below, 100*n_below/num_valid);
    fprintf('  Median Spatial SNR >= %.1f: %d\n', snr_threshold, num_valid - n_below);
    fprintf('  Global median Spatial SNR: %.2f\n', median(median_all));
    fprintf('============================================\n');


    % 3. 拼接并剔除可能的 NaN 或 Inf 异常值
    % Temporal SNR 清洗
    valid_temp_mask = ~isnan(temporal_snr_values) & ~isinf(temporal_snr_values);
    clean_temporal_snr = temporal_snr_values(valid_temp_mask);

    % Spatial SNR 拼接与清洗
    all_spatial_snr_values = vertcat(all_spatial_snr_cells{:}); % 拼接成一个超长列向量
    valid_spat_mask = ~isnan(all_spatial_snr_values) & ~isinf(all_spatial_snr_values);
    clean_spatial_snr = all_spatial_snr_values(valid_spat_mask);

    % 4. 绘制 SNR 分布对比图
    figure('Name', 'Single-Molecule SNR Distributions', 'Color', 'w', 'Position', [100, 100, 1000, 450]);

    % --- 子图 1: Temporal SNR (Trace-level) ---
    subplot(1, 2, 1);
    histogram(clean_temporal_snr, 'BinMethod', 'auto', 'FaceColor', '#0072BD', 'EdgeColor', 'w', 'FaceAlpha', 0.8);
    xlabel('Temporal SNR', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Frequency / Count (Traces)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Temporal SNR\n(1-Step, %d Traces)', length(clean_temporal_snr)), 'FontSize', 14);
    grid on; set(gca, 'FontSize', 11, 'LineWidth', 1.2);

    temp_mean = mean(clean_temporal_snr);
    temp_median = median(clean_temporal_snr);
    xline(temp_mean, '--r', sprintf('Mean: %.2f', temp_mean), 'LineWidth', 2, 'LabelVerticalAlignment', 'top');
    xline(temp_median, '-.g', sprintf('Median: %.2f', temp_median), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');

    % --- 子图 2: Spatial SNR (Frame-level) ---
    subplot(1, 2, 2);
    histogram(clean_spatial_snr, 'BinMethod', 'auto', 'FaceColor', '#D95319', 'EdgeColor', 'w', 'FaceAlpha', 0.8);
    xlabel('Spatial SNR (All ON-state frames)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Frequency / Count (Frames)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Spatial SNR\n(All ON-Frames, N = %d)', length(clean_spatial_snr)), 'FontSize', 14);
    grid on; set(gca, 'FontSize', 11, 'LineWidth', 1.2);

    spat_mean = mean(clean_spatial_snr);
    spat_median = median(clean_spatial_snr);
    xline(spat_mean, '--r', sprintf('Mean: %.2f', spat_mean), 'LineWidth', 2, 'LabelVerticalAlignment', 'top');
    xline(spat_median, '-.g', sprintf('Median: %.2f', spat_median), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');

%% ========== Plot: Spatial SNR Histogram ==========
    % 汇总所有 ON 态帧的 Spatial SNR
    clean_spatial_snr = vertcat(all_spatial_snr_cells{:});
    avg_spatial_snr = mean(clean_spatial_snr);

    figure('Color', 'w', 'Position', [200 200 560 420]);
    histogram(clean_spatial_snr, 'BinWidth', 1, ...
        'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'w', 'LineWidth', 0.5);

    xlabel('SNR', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Counts', 'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.2, 'Box', 'on');

    % 标注均值（模仿示例图的样式）
    xlims = xlim; ylims = ylim;
    text(xlims(2)*0.95, ylims(2)*0.92, ...
        sprintf('\\bfAverage: \\color[rgb]{0.6,0.1,0.1}%.2f', avg_spatial_snr), ...
        'FontSize', 20, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    title('Spatial SNR Distribution (ON frames)', 'FontSize', 13);
    drawnow;

    % 保存图片
    [file_dir, file_base, ~] = fileparts(parameters.file_name);
    fig_name = fullfile(file_dir, [file_base, '_spatial_snr_histogram.png']);
    exportgraphics(gcf, fig_name, 'Resolution', 300);
    fprintf('Spatial SNR histogram saved: %s\n', fig_name);


%% ========== Export: Temporal & Spatial SNR to Excel ==========
    % --- Sheet 1: Per-molecule summary ---
    % 每个合规分子一行：原始索引、位置、Temporal SNR、Spatial SNR (median of ON frames)
    median_spatial_per_mol = zeros(num_valid, 1);
    mean_spatial_per_mol   = zeros(num_valid, 1);
    for i = 1:num_valid
        if ~isempty(all_spatial_snr_cells{i})
            median_spatial_per_mol(i) = median(all_spatial_snr_cells{i});
            mean_spatial_per_mol(i)   = mean(all_spatial_snr_cells{i});
        end
    end

    T_summary = table( ...
        valid_indices(:), ...
        round(stats.pos_mean_px(valid_indices, 1)), ...
        round(stats.pos_mean_px(valid_indices, 2)), ...
        temporal_snr_values(:), ...
        median_spatial_per_mol(:), ...
        mean_spatial_per_mol(:), ...
        'VariableNames', {'MoleculeIndex', 'Row_px', 'Col_px', ...
                        'Temporal_SNR', 'Spatial_SNR_median', 'Spatial_SNR_mean'});

    % --- Sheet 2: Frame-by-frame Spatial SNR (每个分子所有 ON 帧) ---
    % 整理为长表格：MoleculeIndex, Frame, Spatial_SNR
    mol_idx_col   = [];
    frame_col     = [];
    spsnr_col     = [];

    for i = 1:num_valid
        orig_idx = valid_indices(i);
        s = states(orig_idx);
        on_mask = (s.sequence > 0);
        on_frames = find(on_mask);
        spsnr_on  = stats.snr(orig_idx, on_mask)';
        
        n = length(on_frames);
        mol_idx_col = [mol_idx_col; repmat(orig_idx, n, 1)];
        frame_col   = [frame_col; on_frames(:)];
        spsnr_col   = [spsnr_col; spsnr_on(:)];
    end

    T_framewise = table(mol_idx_col, frame_col, spsnr_col, ...
        'VariableNames', {'MoleculeIndex', 'Frame', 'Spatial_SNR'});

    % --- 写入 Excel ---
    excel_name = fullfile(file_dir, [file_base, '_SNR_analysis.xlsx']);

    writetable(T_summary,   excel_name, 'Sheet', 'Per_Molecule_Summary');
    writetable(T_framewise, excel_name, 'Sheet', 'Framewise_Spatial_SNR');

    fprintf('SNR analysis exported: %s\n', excel_name);
    fprintf('  Sheet 1 "Per_Molecule_Summary": %d molecules\n', height(T_summary));
    fprintf('  Sheet 2 "Framewise_Spatial_SNR": %d rows\n', height(T_framewise));

%% Compute from States （已修复：仅计算非背景态的平均亮度，并新增总光子数估计）parameters.bg_range
    %  Compute Summary Statistics
    % ================================================================
    N = numel(states);

    % 1) 态数量分布（每个 emitter 的 step 数）
    nStatesVec = nan(N, 1);
    for n = 1:N
        nStatesVec(n) = states(n).nStates;
    end

    % 2) 每个 emitter 的平均亮度（仅计算非背景态/发光态的均值）
    meanBrightnessVec = nan(N, 1);
    for n = 1:N
        y = states(n).trace;
        seq = states(n).sequence;
        
        % 仅提取 sequence > 0 的帧（排除背景和漂白态）
        on_mask = (seq > 0);
        if any(on_mask)
            meanBrightnessVec(n) = mean(y(on_mask), 'omitnan');
        end
    end

    % 3 & 4 & 5) 仅包含 1 个亮态（排除漂白/暗态后）的分子的亮度、寿命与总光子数
    singleStateLevels = [];
    singleStateLifetimes = []; 
    singleStateTotalPhotons = []; % 【新增】：用于存储单态分子的总光子数
    
    for n = 1:N
        sInfo = states(n).stateInfo;
        if isempty(sInfo), continue; end
        lvls = [sInfo.meanIntensity];
        
        % 使用逐分子的背景估计排除暗态/漂白态
        bg_n = states(n).background;
        bright_lvls = lvls(lvls > bg_n);
        
        % 仅当该分子只有一个亮态时，才将其亮度与寿命记录下来
        if numel(bright_lvls) == 1
            lvl_frame = bright_lvls(1); % 单位: counts/frame
            lt_frames = states(n).lifetime; % 单位: frames
            
            singleStateLevels = [singleStateLevels, lvl_frame];  %#ok<AGROW>
            singleStateLifetimes = [singleStateLifetimes, lt_frames]; %#ok<AGROW>
            
            % 【新增】：计算总光子数 = 亮度(counts/frame) * 寿命(frames)
            singleStateTotalPhotons = [singleStateTotalPhotons, lvl_frame * lt_frames]; %#ok<AGROW>
        end
    end
    singleStateLevels = singleStateLevels(:);
    singleStateLifetimes = singleStateLifetimes(:);
    singleStateTotalPhotons = singleStateTotalPhotons(:); % 【新增】
    
    % 将帧数转换为时间 (秒)
    singleStateLifetimes_sec = singleStateLifetimes * parameters.ex_time;

    fprintf('[analyzeTraces] Summary: %d emitters, %d single-state emitters collected.\n', ...
        N, numel(singleStateLevels));

    %================================================================
    %  Optional: Fitting & Visualization (controlled by viz_enabled)
    % ================================================================
    tau_fit = [];
    photon_fit_mu = []; % 【新增】：总光子数拟合均值
    brightness_fit = struct();
    single_fit = struct(); 
        
    if parameters.viz_enabled
        % --- Lifetime fit (仅使用单态分子) ---
        validLT = singleStateLifetimes_sec(isfinite(singleStateLifetimes_sec) & singleStateLifetimes_sec > 0);
        if numel(validLT) >= 10
            tau_fit = analysis.photophys.fitLifetime_exp(validLT);
            fprintf('[analyzeTraces] 1-State Lifetime fit: tau = %.2f s\n', tau_fit);
        else
            warning('[analyzeTraces] 单态有效寿命数据不足 (%d), 跳过拟合。', numel(validLT));
        end

        % --- Mean brightness fit (lognormal) ---
        validMB = meanBrightnessVec(isfinite(meanBrightnessVec) & meanBrightnessVec > 0);
        if numel(validMB) >= 10
            pd = fitdist(validMB, 'Lognormal');
            brightness_fit.distribution = 'Lognormal';
            brightness_fit.mu    = pd.mu;
            brightness_fit.sigma = pd.sigma;
            brightness_fit.mean  = exp(pd.mu + pd.sigma^2 / 2);
            brightness_fit.median = exp(pd.mu);
            fprintf('[analyzeTraces] Brightness fit (lognormal): median = %.1f, mean = %.1f counts/frame\n', ...
                brightness_fit.median, brightness_fit.mean);
        else
            warning('[analyzeTraces] 有效亮度数据不足 (%d), 跳过拟合。', numel(validMB));
        end
        
        % --- Single State Brightness fit (Normal) ---
        validSingle = singleStateLevels(isfinite(singleStateLevels) & singleStateLevels > 0);
        validSingle_plot = validSingle / parameters.ex_time; % 转换为 counts/s
        if numel(validSingle_plot) >= 10
            pd_single = fitdist(validSingle_plot, 'Normal');
            single_fit.mu = pd_single.mu;
            single_fit.sigma = pd_single.sigma;
            single_fit.median = median(validSingle_plot);
        end
        
        % --- 【新增】：Total Photon Budget fit (Exponential) ---
        validPhotons = singleStateTotalPhotons(isfinite(singleStateTotalPhotons) & singleStateTotalPhotons > 0);
        if numel(validPhotons) >= 10
            pd_photons = fitdist(validPhotons, 'Exponential');
            photon_fit_mu = pd_photons.mu;
            fprintf('[analyzeTraces] Total Photon Budget fit: mean = %.0f counts\n', photon_fit_mu);
        end

        % --- Visualization ---
        % 【修改点】：将画板加宽，并准备使用 2x3 的布局以容纳 5 张图
        figure('Name', 'analyzeTraces Summary', 'Position', [80 80 1400 700]);

        % 1) 态数量分布
        subplot(2, 3, 1);
        validNS = nStatesVec(isfinite(nStatesVec));
        edges = (0:max(validNS)) + 0.5;
        histogram(validNS, edges, 'FaceColor', [0.3 0.6 0.9]);
        xlabel('Number of States');
        ylabel('# Emitters');
        title(sprintf('State Count Distribution (N=%d)', numel(validNS)));
        box on;

        % 2) 平均亮度分布 + lognormal 拟合曲线
        subplot(2, 3, 2);
        validMB_plot = validMB / parameters.ex_time;
        histogram(validMB_plot, 40, 'Normalization', 'pdf', 'FaceColor', [0.9 0.5 0.2]);
        hold on;
        if ~isempty(fieldnames(brightness_fit))
            x_fit = linspace(min(validMB_plot), max(validMB_plot), 200);
            scale = 1 / parameters.ex_time;
            y_fit = pdf('Lognormal', x_fit, ...
                brightness_fit.mu + log(scale), brightness_fit.sigma);
            plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
            legend('Data', sprintf('Lognormal (med=%.0f)', brightness_fit.median / parameters.ex_time));
        end
        xlabel('Mean Brightness (counts/s)');
        ylabel('PDF');
        title(sprintf('Mean Brightness (ON-States, N=%d)', numel(validMB)));
        hold off; box on;

        % 3) 仅1个态的亮度分布 + 中位数 + 正态拟合
        subplot(2, 3, 3);
        if ~isempty(validSingle_plot)
            histogram(validSingle_plot, 40, 'Normalization', 'pdf', 'FaceColor', [0.2 0.8 0.4]);
            hold on;
            
            if ~isempty(fieldnames(single_fit))
                x_fit3 = linspace(min(validSingle_plot), max(validSingle_plot), 200);
                y_fit3 = pdf('Normal', x_fit3, single_fit.mu, single_fit.sigma);
                plot(x_fit3, y_fit3, 'r-', 'LineWidth', 2);
                legend('Data', sprintf('Gaussian (Med=%.0f, \\mu=%.0f)', ...
                    single_fit.median, single_fit.mu));
            else
                legend(sprintf('Median=%.0f', median(validSingle_plot)));
            end
            
            xlabel('Single State Intensity (counts/s)');
            ylabel('PDF');
            title(sprintf('1-State Emitters Brightness (N=%d)', numel(validSingle_plot)));
            hold off; 
        else
            title('No 1-State Emitters Found');
        end
        box on;

        % 4) 仅1个态的寿命分布 + 指数拟合曲线
        subplot(2, 3, 4);
        if ~isempty(validLT)
            histogram(validLT, 40, 'Normalization', 'pdf', 'FaceColor', [0.7 0.3 0.7]);
            hold on;
            if ~isempty(tau_fit)
                x_fit = linspace(0, max(validLT), 200);
                y_fit = (1 / tau_fit) * exp(-x_fit / tau_fit);
                plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
                legend('Data', sprintf('Exp (\\tau=%.2f s)', tau_fit));
            end
            xlabel('Lifetime (s)');
            ylabel('PDF');
            title(sprintf('1-State Emitters Lifetime (N=%d)', numel(validLT)));
            hold off; 
        else
            title('No 1-State Emitters Found');
        end
        box on;
        
        % 5) 【新增】：总光子数分布 + 指数拟合曲线
        subplot(2, 3, 5);
        if ~isempty(validPhotons)
            histogram(validPhotons, 40, 'Normalization', 'pdf', 'FaceColor', [0.8 0.4 0.4]);
            hold on;
            if ~isempty(photon_fit_mu)
                x_fit5 = linspace(0, max(validPhotons), 200);
                % 指数分布的 PDF: f(x) = (1/mu) * exp(-x/mu)
                y_fit5 = (1 / photon_fit_mu) * exp(-x_fit5 / photon_fit_mu);
                plot(x_fit5, y_fit5, 'r-', 'LineWidth', 2);
                legend('Data', sprintf('Exp (Mean=%.0f)', photon_fit_mu));
            end
            xlabel('Total Photons (counts)');
            ylabel('PDF');
            title(sprintf('Total Photon Budget (N=%d)', numel(validPhotons)));
            hold off; 
        else
            title('No Data for Total Photons');
        end
        box on;

        sgtitle('analyzeTraces: Photophysics Summary');
    end