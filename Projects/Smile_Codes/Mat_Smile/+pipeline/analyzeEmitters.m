function result = analyzeEmitters(data, driftResult, parameters)
%PIPELINE.ANALYZEEMITTERS  Emitter 聚类、过滤、统计.
%
%   result = pipeline.analyzeEmitters(data, driftResult, parameters)
%
% INPUTS
%   data         — [H × W × F] 校正后图像栈
%   driftResult  — pipeline.correctDrift 返回的 struct（非 data）
%   parameters   — 参数结构体
%
% OUTPUT  (struct)
%   result.emitters_raw   — 聚类后、过滤前的 emitter 数组
%   result.emitters_filt  — 过滤后的 emitter 数组
%   result.stats          — collectEmitterStatistics 输出
%   result.stats_jump     — jumping filter 诊断信息
%   result.n_before       — 过滤前 emitter 数
%   result.n_after        — 过滤后 emitter 数

    arguments
        data
        driftResult  (1,1) struct
        parameters   (1,1) struct
    end

    % ================================================================
    %  Unpack
    % ================================================================
    loc_total  = driftResult.loc_total;
    corrected  = driftResult.corrected_super_loc_total;
    brightness = corrected(:, 3);
    background = corrected(:, 4);
    sigma      = corrected(:, 5);
    sigma_loc  = corrected(:, 6);

    % ================================================================
    %  Emitter Clustering
    % ================================================================
    fprintf('[analyzeEmitters] Finding emitters ...\n');
    bleach_frames = parameters.bleach_time / parameters.one_frame_time;
    emitters = postproc.emitter.findEmitters(data, loc_total, ...
        bleach_frames, parameters.searching_radius);
    emitters = postproc.emitter.mergeEmitters(emitters, parameters.searching_radius);

    n_before = numel(emitters);
    fprintf('[analyzeEmitters] Found %d emitters (before filtering).\n', n_before);

    % ================================================================
    %  Filtering
    %  顺序: short-lived → first-frame → last-frame residual → position-jumping
    %  每一步通过对应参数是否为 NaN 来决定是否执行；NaN 则跳过。
    % ================================================================
    fprintf('[analyzeEmitters] Filtering ...\n');
    emitters_filt = emitters;

    % --- 1. Short-lived filtering (参数: livetime_th) ---
    if ~isnan(parameters.livetime_th)
        min_frames    = round(parameters.livetime_th / parameters.ex_time);
        emitters_filt = postproc.emitter.filterEmitters_short( ...
                            emitters_filt, min_frames, 'consecutive');
        fprintf('[analyzeEmitters]   short-lived filter applied (livetime_th = %.3g)\n', ...
            parameters.livetime_th);
    else
        fprintf('[analyzeEmitters]   short-lived filter skipped (livetime_th = NaN)\n');
    end

    % --- 2. First-frame filtering (参数: filter_firstframe) ---
    if ~isnan(parameters.filter_firstframe)
        emitters_filt = postproc.emitter.filterEmitters_firstframe(emitters_filt);
        fprintf('[analyzeEmitters]   first-frame filter applied\n');
    else
        fprintf('[analyzeEmitters]   first-frame filter skipped (filter_firstframe = NaN)\n');
    end

    % --- 3. Last-frame (end) filtering (参数: filter_end) ---
    if ~isnan(parameters.filter_end)
        emitters_filt = postproc.emitter.filterEmitters_end(emitters_filt);
        fprintf('[analyzeEmitters]   end-frame filter applied\n');
    else
        fprintf('[analyzeEmitters]   end-frame filter skipped (filter_end = NaN)\n');
    end

    % --- 4. Position-jumping filtering (参数: jump_threshold) ---
    stats_jump = [];
    if ~isnan(parameters.jump_threshold)
        [emitters_filt, stats_jump] = postproc.emitter.filterEmitters_jumping( ...
            emitters_filt, parameters.jump_threshold);
        fprintf('[analyzeEmitters]   jumping filter applied (jump_threshold = %.3g)\n', ...
            parameters.jump_threshold);
    else
        fprintf('[analyzeEmitters]   jumping filter skipped (jump_threshold = NaN)\n');
    end

    n_after = numel(emitters_filt);
    fprintf('[analyzeEmitters] After filtering: %d / %d\n', n_after, n_before);

    if parameters.viz_enabled && ~isempty(stats_jump)
        postproc.emitter.plotJumpStats(stats_jump);
    end

    % ================================================================
    %  Visualization: Emitter Trajectories
    % ================================================================
    if parameters.viz_enabled
        figure;
        viz.plotImage(data, 1:parameters.frames, 'gray', ...
            'Drift-corrected Image and Emitters', 'max');
        hold on;
        for k = 1:numel(emitters_filt)
            plot(emitters_filt(k).col + 0.5, ...
                 emitters_filt(k).row + 0.5, '-', 'LineWidth', 1);
        end
        title('Emitter trajectories overlaid on image');
    end

    % ================================================================
    %  Collect Statistics
    % ================================================================
    fprintf('[analyzeEmitters] Collecting statistics ...\n');
    stats = postproc.emitter.collectEmitterStatistics( ...
        emitters_filt, parameters.frames, parameters.ex_time, ...
        parameters.interval, brightness, sigma, sigma_loc, background);

    if parameters.viz_enabled
        postproc.emitter.plotEmitterStatistics(stats, 50);
    end

    % ================================================================
    %  Pack output
    % ================================================================
    result.emitters_raw  = emitters;
    result.emitters_filt = emitters_filt;
    result.stats         = stats;
    result.stats_jump    = stats_jump;
    result.n_before      = n_before;
    result.n_after       = n_after;
end