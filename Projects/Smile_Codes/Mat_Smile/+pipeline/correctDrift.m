function [data, result] = correctDrift(windowed_raw_data, locResult, parameters, options)
%PIPELINE.CORRECTDRIFT  漂移估计 + 图像/定位校正 + 裁切.
%
%   [data, result] = pipeline.correctDrift(windowed_raw_data, locResult, parameters)
%   [data, result] = pipeline.correctDrift(..., 'cacheFile', filepath)
%
% INPUTS
%   windowed_raw_data  — [H × W × F] 未校正图像栈
%   locResult          — pipeline.detectAndLocalize 的输出
%   parameters         — 参数结构体
%
% OUTPUTS
%   data   — [row_width × col_width × F] 漂移校正后图像栈
%            （独立输出，避免打包进 struct 造成内存拷贝）
%   result — struct:
%     .delta_sum                 — [F × 2] 各帧累积漂移 (row, col)
%     .corrected_super_loc_total — [M' × 7] 校正后保留的完整定位
%     .loc_total                 — [M' × 3] 校正后 [row, col, frame]
%     .loc_idx                   — 保留定位在 locResult 中的索引
%     .uncertainty               — 中位定位精度 (px)

    arguments
        windowed_raw_data
        locResult         (1,1) struct
        parameters        (1,1) struct
        options.cacheFile       = ''
    end

    % ================================================================
    %  Cache — 加载
    %  注意：data 可能很大（数 GB），cacheFile 会相应较大。
    %        若仅需缓存轻量结构体，可改为只存 result 并用 TIFF 存 data。
    % ================================================================
    if ~isempty(options.cacheFile) && isfile(options.cacheFile)
        S = load(options.cacheFile, 'data', 'result');
        data   = S.data;
        result = S.result;
        fprintf('[correctDrift] Loaded from cache: %s\n', options.cacheFile);
        return
    end

    % ================================================================
    %  Drift Estimation
    % ================================================================
    uncertainty = median(locResult.super_loc_total(:, 6));

    if parameters.drift_corr
        fprintf('[correctDrift] Estimating drift (cloud-XY) ...\n');
        delta_sum = postproc.drift.estimateDrift_cloudxy( ...
            windowed_raw_data, locResult.loc_total, ...
            parameters.drift_frames_per_seg, uncertainty, ...
            parameters.drift_min_locs, false);
    else
        delta_sum = zeros(parameters.frames, 2);
        fprintf('[correctDrift] Drift correction disabled; delta_sum = 0.\n');
    end

    % ================================================================
    %  Apply Drift Correction + Crop
    % ================================================================
    fprintf('[correctDrift] Applying correction & cropping to %d×%d ...\n', ...
        parameters.row_width, parameters.col_width);
    [data, loc_total, loc_idx] = postproc.drift.correctDrift( ...
        windowed_raw_data, locResult.loc_total, delta_sum, ...
        parameters.edge, parameters.row_width, parameters.col_width);

    corrected_super_loc_total = locResult.super_loc_total(loc_idx, :);
    % 注：corrected_super_loc_total 的 row/col（第1-2列）仍为未校正值，
    %     校正后的位置由 loc_total 提供。第3-7列（亮度/背景/sigma 等）保持不变。

    % ================================================================
    %  Visualization: Drift Tracking
    % ================================================================
    if parameters.viz_enabled
        viz.plotTracking(delta_sum, parameters.one_frame_time, [], 'Drift Tracking');
    end

    % ================================================================
    %  Visualization: MLE Distribution
    % ================================================================
    if parameters.viz_enabled
        brightness_v = corrected_super_loc_total(:, 3);
        background_v = corrected_super_loc_total(:, 4);
        sigma_v      = corrected_super_loc_total(:, 5);
        sigma_loc_v  = corrected_super_loc_total(:, 6);

        figure;
        subplot(4,1,1);
        histogram(brightness_v, 180);
        xlabel('Photons (N)'); ylabel('Count');
        title('Brightness distribution (by locs)');

        subplot(4,1,2);
        histogram(background_v, 140);
        xlabel('Photons (N)'); ylabel('Count');
        title('Background distribution (by locs)');

        subplot(4,1,3);
        histogram(sigma_v, 120);
        xlabel('Sigma (px)'); ylabel('Count');
        title('Fitted PSF width (\sigma)');

        subplot(4,1,4);
        histogram(sigma_loc_v, 100);
        xlabel('Localization precision (px)'); ylabel('Count');
        title('\sigma_{loc} (precision, CRLB)');
    end

    % ================================================================
    %  Visualization: Corrected Localizations
    % ================================================================
    if parameters.viz_enabled
        frames_to_show = 1:min(parameters.viz_max_frames, parameters.frames);
        figure;
        viz.plotImage(data, frames_to_show, 'hot', ...
            'Drift-corrected Image and Locs', 'max');
        hold on;
        viz.overlayLocs(loc_total, frames_to_show, true);
    end

    % ================================================================
    %  Pack output
    % ================================================================
    result.delta_sum                 = delta_sum;
    result.corrected_super_loc_total = corrected_super_loc_total;
    result.loc_total                 = loc_total;
    result.loc_idx                   = loc_idx;
    result.uncertainty               = uncertainty;

    % ================================================================
    %  Cache — 保存
    % ================================================================
    if ~isempty(options.cacheFile)
        save(options.cacheFile, 'data', 'result', '-v7.3');
        fprintf('[correctDrift] Cached: %s\n', options.cacheFile);
    end
end