function fileResult = pipeline_electric(file_name, parameters)
% pipeline_charge  处理单个 SIF 文件的电响应发光分析全流程
%
% 输入:
%   file_name   - SIF 文件路径
%   parameters  - 参数结构体（可选字段见下）
%
% 输出:
%   fileResult  - 结构体，包含所有中间结果和最终 metrics

    arguments
        file_name   char
        parameters  struct = struct()
    end

%% ===================== Default Parameters =====================
    parameters.file_name = file_name;
    parameters = setDefault(parameters, 'interval', 0);

    % --- Windowing ---
    parameters = setDefault(parameters, 'row_range', 100:412);
    parameters = setDefault(parameters, 'col_range', 100:412);

    % --- Detection / Localization ---
    parameters = setDefault(parameters, 'k_sigma',        2.0);
    parameters = setDefault(parameters, 'edge',           5);
    parameters = setDefault(parameters, 'viz_enabled',    false);
    parameters = setDefault(parameters, 'viz_max_frames', 5000);

    % --- Drift Estimation ---
    parameters = setDefault(parameters, 'drift_corr',           true);
    parameters = setDefault(parameters, 'drift_frames_per_seg', 15);
    parameters = setDefault(parameters, 'drift_min_locs',       30);

    % --- Drift Correction ---
    parameters = setDefault(parameters, 'row_width', 257);
    parameters = setDefault(parameters, 'col_width', 257);

    % --- Emitter Analysis ---
    parameters = setDefault(parameters, 'bleach_time',      50);
    parameters = setDefault(parameters, 'searching_radius', 2);

    % --- Emitter Filtering ---
    parameters = setDefault(parameters, 'livetime_th',       0.2);
    parameters = setDefault(parameters, 'filter_firstframe', true);
    parameters = setDefault(parameters, 'filter_end',        nan);   % nan = 不过滤（电致实验特有）
    parameters = setDefault(parameters, 'jump_threshold',    3.5);

    % --- Trace Extraction ---
    parameters = setDefault(parameters, 'int_range',             4);
    parameters = setDefault(parameters, 'bg_range',              7);
    parameters = setDefault(parameters, 'bg_extracted',          true);
    parameters = setDefault(parameters, 'bg_extraction_method',  'median');

    % --- State Analysis ---
    parameters = setDefault(parameters, 'emittersFIT_enabled',    true);
    parameters = setDefault(parameters, 'state_method',       'CHANGEPOINT');
    parameters = setDefault(parameters, 'state_penalty',      2.5);
    parameters = setDefault(parameters, 'state_min_seg_len',  2.0);
    parameters = setDefault(parameters, 'state_merge_thr',    1.3);
    parameters = setDefault(parameters, 'state_bicPenalty',   2.0);
    parameters = setDefault(parameters, 'state_bleach_tail',  50);
    parameters = setDefault(parameters, 'state_bg_threshold', 60);
    parameters = setDefault(parameters, 'state_min_state_sep', 1.4);

    % --- Charge-specific ---
    parameters = setDefault(parameters, 'initial_n_frames', 50);
    % parameters = setDefault(parameters, 'delta_min_occ',    4);
    parameters = setDefault(parameters, 'voltage_period',  5.0);   % 电压完整周期 (秒)
    parameters = setDefault(parameters, 'n_periods_fit', 4);

%% ===================== Read & Window =====================

    [raw_data, ex_time, ~] = io.readSIFData(file_name);
    parameters.ex_time        = ex_time;
    parameters.one_frame_time = ex_time + parameters.interval;

    windowed_raw_data = raw_data(parameters.row_range, parameters.col_range, :);
    parameters.frames = size(windowed_raw_data, 3);
    clear raw_data;

    fprintf('[pipeline_charge] Windowed data: %d × %d × %d\n', ...
        size(windowed_raw_data,1), size(windowed_raw_data,2), parameters.frames);

%% ===================== Module 1: Detect & Localize =====================

    [file_dir, file_base, ~] = fileparts(file_name);
    cache_loc = fullfile(file_dir, [file_base, '_locResult.mat']);

    locResult = pipeline.detectAndLocalize(windowed_raw_data, parameters, ...
        'cacheFile', cache_loc);

%% ===================== Module 2: Drift Correction =====================

    [data, driftResult] = pipeline.correctDrift( ...
        windowed_raw_data, locResult, parameters);

    clear windowed_raw_data;

%% ===================== Module 3: Emitter Analysis =====================

    emitterResult = pipeline.analyzeEmitters(data, driftResult, parameters);

%% ===================== Module 4: Trace & State Analysis =====================

    traceResult = pipeline.analyzeTraces(data, emitterResult, parameters);

%% ===================== Module 5: Charge-specific Metrics =====================

    metrics = computeChargeMetrics(traceResult.trace, parameters);

%% ===================== Pack Output =====================

    fileResult = struct( ...
        'file_name',     file_name, ...
        'parameters',    parameters, ...
        'emitterResult', emitterResult, ...
        'traceResult',   traceResult, ...
        'metrics',       metrics, ...
        'n_emitters',    size(traceResult.trace, 1));

%% ===================== Auto-save =====================

    save_name = fullfile(file_dir, [file_base '_chargeResult.mat']);
    save(save_name, 'fileResult', '-v7.3');
    fprintf('[pipeline_charge] Result saved: %s\n', save_name);

end

%% ======================== Local Functions ========================
function [metrics, debug] = computeChargeMetrics(traces, params)
    % computeChargeMetrics - Detect voltage onset and compute charge metrics
    %
    % Inputs:
    %   traces  - [n_emitters x n_frames] raw intensity traces
    %   params  - struct with optional fields:
    %       .voltage_period   - fallback full period if auto-detection fails
    %       .n_periods_fit    - number of full periods to fit (default: 4)
    %       .min_onset_frame  - minimum onset frame to search (default: auto)
    %       .viz_enabled      - plot debug figures (default: false)
    %
    % Outputs:
    %   metrics - struct with onset, brightness, stateMax, stateMin, etc.
    %   debug   - struct with intermediate data for diagnostics

    %% Parse parameters
    if ~isfield(params, 'n_periods_fit'),   params.n_periods_fit = 4;     end
    if ~isfield(params, 'viz_enabled'),     params.viz_enabled = false;   end
    if ~isfield(params, 'min_onset_frame'), params.min_onset_frame = [];  end

    [n_emitters, n_frames] = size(traces);

    %% 1. Ensemble average
    trace_means = mean(traces, 2);
    trace_means(trace_means == 0) = 1;
    norm_traces = traces ./ trace_means;
    avg_trace = mean(norm_traces, 1);

    %% 2. Auto-detect period from autocorrelation
    sig = avg_trace - mean(avg_trace);
    [acf, lags] = xcorr(sig, 'coeff');
    mid = find(lags == 0);
    acf_pos = acf(mid+1:end);

    min_search_lag = 20;
    [~, pk_locs] = findpeaks(acf_pos(min_search_lag:end), 'MinPeakProminence', 0.05);

    if ~isempty(pk_locs)
        full_period = pk_locs(1) + min_search_lag - 1;
    elseif isfield(params, 'voltage_period')
        full_period = params.voltage_period;
    else
        error('Cannot auto-detect period. Provide params.voltage_period');
    end
    half_period = round(full_period / 2);

    %% 3. Onset detection via square wave template matching
    N_fit = params.n_periods_fit * full_period;

    if isempty(params.min_onset_frame)
        min_onset = max(10, round(half_period / 2));
    else
        min_onset = params.min_onset_frame;
    end
    max_onset = n_frames - N_fit;

    if max_onset <= min_onset
        error('Trace too short: need at least %d frames after onset for %d periods.', ...
            N_fit, params.n_periods_fit);
    end

    % Square wave template: +1 for first half-period, -1 for second, repeat
    sq_template = ones(1, N_fit);
    for k = 0:N_fit-1
        if mod(floor(k / half_period), 2) == 1
            sq_template(k+1) = -1;
        end
    end

    cost = nan(1, n_frames);
    fit_b = nan(1, n_frames);

    for t = min_onset:max_onset
        % Before onset: constant model
        before = avg_trace(1:t-1);
        mu_b = mean(before);
        rss_before = sum((before - mu_b).^2);
        
        % After onset: y = a + b * sq_template
        after = avg_trace(t : t + N_fit - 1);
        S1 = sum(sq_template);
        Sy = sum(after);
        Syq = sum(after .* sq_template);
        n_a = N_fit;
        
        det_val = n_a * n_a - S1 * S1;
        if abs(det_val) < 1e-10, continue; end
        
        a_val = (n_a * Sy - S1 * Syq) / det_val;
        b_val = (n_a * Syq - S1 * Sy) / det_val;
        
        rss_after = sum((after - a_val - b_val * sq_template).^2);
        cost(t) = rss_before + rss_after;
        fit_b(t) = b_val;
    end

    % Global minimum
    [~, idx] = min(cost(min_onset:max_onset));
    onset = min_onset + idx - 1;

    %% 4. Phase assignment (bright first or dark first)
    first_half_bright = (fit_b(onset) > 0);

    %% 5. Label frames within the fitted region only
    fit_end = onset + N_fit - 1;
    frames_fit = onset:fit_end;
    blk_idx = floor((frames_fit - onset) / half_period);
    phase = mod(blk_idx, 2);  % 0 = first half, 1 = second half

    if first_half_bright
        idx_bright = frames_fit(phase == 0);
        idx_dark   = frames_fit(phase == 1);
    else
        idx_bright = frames_fit(phase == 1);
        idx_dark   = frames_fit(phase == 0);
    end

    %% 6. Per-emitter metrics
    initial_brightness = mean(traces(:, 1:onset-1), 2);
    stateMax = mean(traces(:, idx_bright), 2);
    stateMin = mean(traces(:, idx_dark), 2);

    metrics.voltage_onset_frame = onset;
    metrics.half_period = half_period;
    metrics.full_period = full_period;
    metrics.first_half_bright = first_half_bright;
    metrics.initial_brightness = initial_brightness;
    metrics.stateMax = stateMax;
    metrics.stateMin = stateMin;
    metrics.normed_delta = (stateMax - stateMin) ./ initial_brightness;
    metrics.bright_dark_ratio = (stateMax - initial_brightness) ./ (initial_brightness - stateMin);
    metrics.bright_frames = idx_bright;
    metrics.dark_frames = idx_dark;
    metrics.fit_end_frame = fit_end;

    %% 7. Debug struct
    debug.avg_trace = avg_trace;
    debug.cost = cost;
    debug.onset = onset;
    debug.half_period = half_period;
    debug.full_period = full_period;
    debug.N_fit = N_fit;
    debug.fit_end = fit_end;

    %% 8. Visualization
    if params.viz_enabled
        plotOnsetDebug(avg_trace, cost, onset, half_period, full_period, ...
            N_fit, first_half_bright, idx_bright, idx_dark, acf_pos, min_search_lag);
    end

end

%% ========== Local plotting function ==========
function plotOnsetDebug(avg_trace, cost, onset, half_period, full_period, ...
        N_fit, first_half_bright, idx_bright, idx_dark, acf_pos, min_search_lag)

    n_frames = length(avg_trace);
    fit_end = onset + N_fit - 1;

    % Build display fit line
    fit_line = nan(1, n_frames);
    fit_line(1:onset-1) = mean(avg_trace(1:onset-1));

    mu_bright = mean(avg_trace(idx_bright));
    mu_dark = mean(avg_trace(idx_dark));

    for f = onset:fit_end
        blk = floor((f - onset) / half_period);
        is_first_half = (mod(blk, 2) == 0);
        if first_half_bright
            fit_line(f) = mu_bright * is_first_half + mu_dark * (~is_first_half);
        else
            fit_line(f) = mu_dark * is_first_half + mu_bright * (~is_first_half);
        end
    end

    figure('Name', 'Onset Detection', 'Position', [100 100 1200 800]);

    % Panel 1: Autocorrelation
    subplot(3,1,1);
    max_lag_plot = min(length(acf_pos), round(n_frames/2));
    plot(1:max_lag_plot, acf_pos(1:max_lag_plot), 'b-', 'LineWidth', 1);
    hold on;
    xline(full_period, 'r--', 'LineWidth', 2);
    xline(half_period, 'g--', 'LineWidth', 1.5);
    xlabel('Lag (frames)'); ylabel('ACF');
    title(sprintf('Autocorrelation — full period = %d, half period = %d', full_period, half_period));
    legend('ACF', 'Full period', 'Half period', 'Location', 'northeast');
    grid on;

    % Panel 2: Ensemble average + fit
    subplot(3,1,2);
    plot(avg_trace, 'b-', 'LineWidth', 0.5); hold on;
    plot(fit_line, 'r-', 'LineWidth', 1.5);
    xline(onset, 'r--', 'LineWidth', 2);
    xline(fit_end, 'k:', 'LineWidth', 1.5);
    xlabel('Frame'); ylabel('Norm. intensity');
    title(sprintf('Ensemble Average — onset = %d, fit end = %d', onset, fit_end));
    legend('Data', 'Square wave fit', 'Onset', 'Fit end', 'Location', 'northeast');
    grid on;

    % Panel 3: Cost curve
    subplot(3,1,3);
    plot(cost, 'b-', 'LineWidth', 1); hold on;
    xline(onset, 'r--', 'LineWidth', 2);
    xlabel('Frame'); ylabel('RSS');
    title('Onset cost curve (RSS_{before} + RSS_{after})');
    grid on;

end

function s = setDefault(s, field, value)
    if ~isfield(s, field)
        s.(field) = value;
    end
end