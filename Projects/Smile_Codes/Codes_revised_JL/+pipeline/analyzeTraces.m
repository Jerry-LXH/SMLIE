function result = analyzeTraces(data, emitterResult, parameters)
%PIPELINE.ANALYZETRACES  Trace 提取 + State 分析 + 光物理拟合.
%
%   result = pipeline.analyzeTraces(data, emitterResult, parameters)
%
% INPUTS
%   data           — [H × W × F] 校正后图像栈
%   emitterResult  — pipeline.analyzeEmitters 的输出
%   parameters     — 参数结构体
%
% OUTPUT  (struct)
%   result.trace          — [N × F] 背景扣除后的强度 trace
%   result.states         — [1 × N] state analysis 结构体数组
%   result.brightestMat   — [N × 8] 每个 emitter 最亮态汇总
%                           列: [idx, label, mean, std, occupancy,
%                                nStates, bleachFrame, lifetime_frames]
%   result.tau_survival   — 基于 emitter 存活时间的指数拟合
%   result.tau_states     — 基于 state lifetime 的指数拟合
%   result.pd_brightness  — 亮度分布的对数正态拟合

    arguments
        data
        emitterResult  (1,1) struct
        parameters     (1,1) struct
    end

    stats = emitterResult.stats;

    % ================================================================
    %  Preliminary Photophysics Fits（仅依赖 emitter 统计量）
    % ================================================================
    if parameters.emittersFIT_enabled
        fprintf('[analyzeTraces] Fitting emitter survival lifetime ...\n');
        tau_survival = analysis.photophys.fitLifetime_exp(stats.survival_sec);

        fprintf('[analyzeTraces] Fitting emitter brightness distribution ...\n');
        pd_brightness = analysis.photophys.fitBrightness_lognormal(stats.brightness_mean);
    end
    % ================================================================
    %  Trace Extraction
    % ================================================================
    fprintf('[analyzeTraces] Extracting intensity traces ...\n');
    trace = analysis.photophys.extractTrace( ...
        stats.pos_mean_px, data, ...
        parameters.int_range, parameters.bg_extracted, parameters.bg_range, parameters.bg_extraction_method);

    % ================================================================
    %  State Analysis (Changepoint)
    % ================================================================
    fprintf('[analyzeTraces] Running state analysis (%s) ...\n', ...
        parameters.state_method);
    states = analysis.photophys.analyzeStates(trace, ...
        parameters.state_method, ...
        'bicPenalty',  parameters.state_bicPenalty, ...
        'penalty',      parameters.state_penalty, ...
        'minSegLen',    parameters.state_min_seg_len, ...
        'mergeThr',     parameters.state_merge_thr, ...
        'bleachTail',   parameters.state_bleach_tail, ...
        'bgThreshold',  parameters.state_bg_threshold, ...
        'minStateSep',  parameters.state_min_state_sep);

    

    % ================================================================
    %  Compute Summary Statistics
    % ================================================================
    N = numel(states);

    % 1) 态数量分布（每个 emitter 的 step 数）
    nStatesVec = nan(N, 1);
    for n = 1:N
        nStatesVec(n) = states(n).nStates;
    end

    % 2) 每个 emitter 的平均亮度（整条 trace 的均值，排除 NaN）
    meanBrightnessVec = mean(trace, 2, 'omitnan');

    % 3) 所有拟合到的态的亮度（排除暗态，合并所有 emitter）
    allStateLevels = [];
    for n = 1:N
        sInfo = states(n).stateInfo;
        if isempty(sInfo), continue; end
        lvls = [sInfo.meanIntensity];
        % 使用逐分子的背景估计排除暗态
        bg_n = states(n).background;
        lvls = lvls(lvls > bg_n);
        allStateLevels = [allStateLevels, lvls(:)'];  %#ok<AGROW>
    end
    allStateLevels = allStateLevels(:);

    % 4) 寿命分布（帧数 → 时间）
    lifetimeVec = nan(N, 1);
    for n = 1:N
        lifetimeVec(n) = states(n).lifetime;
    end
    lifetimeVec_sec = lifetimeVec * parameters.ex_time;

    fprintf('[analyzeTraces] Summary: %d emitters, %d total state levels collected.\n', ...
        N, numel(allStateLevels));

    % ================================================================
    %  Optional: Fitting & Visualization (controlled by viz_enabled)
    % ================================================================
        tau_fit = [];
        brightness_fit = struct();
        
    if parameters.viz_enabled
            % --- Lifetime fit ---
        validLT = lifetimeVec_sec(isfinite(lifetimeVec_sec) & lifetimeVec_sec > 0);
        if numel(validLT) >= 10
            tau_fit = analysis.photophys.fitLifetime_exp(validLT);
            fprintf('[analyzeTraces] Lifetime fit: tau = %.2f s\n', tau_fit);
        else
            warning('[analyzeTraces] 有效寿命数据不足 (%d), 跳过拟合。', numel(validLT));
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
        % --- Visualization ---
        figure('Name', 'analyzeTraces Summary', 'Position', [80 80 1100 700]);

        % 1) 态数量分布
        subplot(2, 2, 1);
        validNS = nStatesVec(isfinite(nStatesVec));
        edges = (0:max(validNS)) + 0.5;
        histogram(validNS, edges, 'FaceColor', [0.3 0.6 0.9]);
        xlabel('Number of States');
        ylabel('# Emitters');
        title(sprintf('State Count Distribution (N=%d)', numel(validNS)));
        box on;

        % 2) 平均亮度分布 + lognormal 拟合曲线
        subplot(2, 2, 2);
        validMB_plot = validMB / parameters.ex_time;
        h = histogram(validMB_plot, 40, 'Normalization', 'pdf', 'FaceColor', [0.9 0.5 0.2]);
        hold on;
        if ~isempty(fieldnames(brightness_fit))
            x_fit = linspace(min(validMB_plot), max(validMB_plot), 200);
            % 变换到 counts/s 域的 lognormal pdf
            scale = 1 / parameters.ex_time;
            y_fit = pdf('Lognormal', x_fit, ...
                brightness_fit.mu + log(scale), brightness_fit.sigma);
            plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
            legend('Data', sprintf('Lognormal (med=%.0f)', brightness_fit.median / parameters.ex_time));
        end
        xlabel('Mean Brightness (counts/s)');
        ylabel('PDF');
        title(sprintf('Mean Brightness (N=%d)', numel(validMB)));
        hold off; box on;

        % 3) 所有态亮度分布
        subplot(2, 2, 3);
        histogram(allStateLevels / parameters.ex_time, 50, 'FaceColor', [0.2 0.8 0.4]);
        xlabel('State Intensity (counts/s)');
        ylabel('# States');
        title(sprintf('All State Levels (N=%d, excl. dark)', numel(allStateLevels)));
        box on;

        % 4) 寿命分布 + 指数拟合曲线
        subplot(2, 2, 4);
        h = histogram(validLT, 40, 'Normalization', 'pdf', 'FaceColor', [0.7 0.3 0.7]);
        hold on;
        if ~isempty(tau_fit)
            x_fit = linspace(0, max(validLT), 200);
            y_fit = (1 / tau_fit) * exp(-x_fit / tau_fit);
            plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
            legend('Data', sprintf('Exp (\\tau=%.2f s)', tau_fit));
        end
        xlabel('Lifetime (s)');
        ylabel('PDF');
        title(sprintf('Lifetime Distribution (N=%d)', numel(validLT)));
        hold off; box on;

        sgtitle('analyzeTraces: Photophysics Summary');
    end

    % ================================================================
    %  Pack Output
    % ================================================================
    result.trace              = trace;
    result.states             = states;

    result.nStatesVec         = nStatesVec;          % (N×1) 每个 emitter 的态数量
    result.meanBrightnessVec  = meanBrightnessVec;   % (N×1) 每个 emitter 的平均亮度
    result.allStateLevels     = allStateLevels;      % (M×1) 所有非暗态的亮度值
    result.lifetimeVec_sec    = lifetimeVec_sec;     % (N×1) 寿命 (s)
    result.tau_fit            = tau_fit;             % 指数拟合寿命常数 (s), 或 []
    result.brightness_fit    = brightness_fit;      % 亮度 lognormal 拟合参数, 或空 struct
end