function result = detectAndLocalize(windowed_raw_data, parameters, options)
%PIPELINE.DETECTANDLOCALIZE  粗定位 (detection) + MLE 精定位 (localization).
%
%   result = pipeline.detectAndLocalize(windowed_raw_data, parameters)
%   result = pipeline.detectAndLocalize(..., 'cacheFile', filepath)
%
% INPUTS
%   windowed_raw_data  — [H × W × F] 裁切后的原始图像栈
%   parameters         — 参数结构体（需含 k_sigma, edge, viz_enabled,
%                         viz_max_frames, frames, file_name）
%
% NAME-VALUE
%   cacheFile — .mat 缓存路径。若文件存在则直接加载并返回；
%               若不存在则正常计算后自动保存。
%
% OUTPUT  (struct)
%   result.detected_total   — [N × 3] 粗定位 [row, col, frame]
%   result.super_loc_total  — [M × 7] MLE 精定位（已过滤坏点）
%   result.loc_total        — [M × 3] 精简三列 [row, col, frame]
%   result.n_raw            — 过滤前定位数
%   result.n_filtered       — 过滤后定位数

    arguments
        windowed_raw_data
        parameters        (1,1) struct
        options.cacheFile       = ''
    end

    % ================================================================
    %  Cache — 加载
    % ================================================================
    if ~isempty(options.cacheFile) && isfile(options.cacheFile)
        S = load(options.cacheFile, 'result');
        result = S.result;
        fprintf('[detectAndLocalize] Loaded from cache: %s\n', options.cacheFile);

        % --- 提示性输出：显示缓存中的结果情况 ---
        fprintf('[detectAndLocalize] Cache contents:\n');
        fprintf('[detectAndLocalize]   detected_total:  %d detections\n', size(result.detected_total, 1));
        fprintf('[detectAndLocalize]   super_loc_total: %d locs, %d columns\n', ...
            size(result.super_loc_total, 1), size(result.super_loc_total, 2));
        fprintf('[detectAndLocalize]   loc_total:       %d locs\n', size(result.loc_total, 1));
        fprintf('[detectAndLocalize]   n_raw = %d, n_filtered = %d\n', ...
            result.n_raw, result.n_filtered);
        if ~isempty(result.loc_total)
            fprintf('[detectAndLocalize]   Row range:   [%.2f, %.2f]\n', ...
                min(result.loc_total(:,1)), max(result.loc_total(:,1)));
            fprintf('[detectAndLocalize]   Col range:   [%.2f, %.2f]\n', ...
                min(result.loc_total(:,2)), max(result.loc_total(:,2)));
            fprintf('[detectAndLocalize]   Frame range: [%d, %d]\n', ...
                min(result.loc_total(:,3)), max(result.loc_total(:,3)));
        end
        return
    end

    % ================================================================
    %  Legacy File Check
    %  若存在旧版 _uncorrected_super_loc_total.mat，直接加载作为定位结果，
    %  跳过 detection 和 MLE localization。
    % ================================================================
    %[file_dir, file_base, ~] = fileparts(parameters.file_name);
    %legacy_file = fullfile(file_dir, [file_base, '_uncorrected_super_loc_total.mat']);

    if isfile(legacy_file)
        fprintf('[detectAndLocalize] Found legacy file: %s\n', legacy_file);
        fprintf('[detectAndLocalize] Skipping detection & MLE, using legacy localizations.\n');

        S_legacy = load(legacy_file);

        % --- 提示性输出：显示文件中的变量情况 ---
        fnames = fieldnames(S_legacy);
        fprintf('[detectAndLocalize] Legacy file contains %d variable(s):\n', numel(fnames));
        for ii = 1:numel(fnames)
            var_size = size(S_legacy.(fnames{ii}));
            fprintf('[detectAndLocalize]   %-40s  size = [%s]\n', ...
                fnames{ii}, strjoin(string(var_size), ' x '));
        end

        % 兼容不同变量名：优先查找 super_loc_total，其次 uncorrected_super_loc_total
        if isfield(S_legacy, 'super_loc_total')
            super_loc_total = S_legacy.super_loc_total;
            fprintf('[detectAndLocalize] Using variable: super_loc_total\n');
        elseif isfield(S_legacy, 'uncorrected_super_loc_total')
            super_loc_total = S_legacy.uncorrected_super_loc_total;
            fprintf('[detectAndLocalize] Using variable: uncorrected_super_loc_total\n');
        else
            % 若只有一个变量，直接取出
            super_loc_total = S_legacy.(fnames{1});
            fprintf('[detectAndLocalize] Using fallback variable: %s\n', fnames{1});
        end

        loc_total = super_loc_total(:, [1, 2, 7]);
        n_filt    = size(super_loc_total, 1);

        fprintf('[detectAndLocalize] Legacy localizations loaded: %d locs, %d columns.\n', ...
            n_filt, size(super_loc_total, 2));
        fprintf('[detectAndLocalize]   Row range:   [%.2f, %.2f]\n', ...
            min(super_loc_total(:,1)), max(super_loc_total(:,1)));
        fprintf('[detectAndLocalize]   Col range:   [%.2f, %.2f]\n', ...
            min(super_loc_total(:,2)), max(super_loc_total(:,2)));
        fprintf('[detectAndLocalize]   Frame range: [%d, %d]\n', ...
            min(super_loc_total(:,7)), max(super_loc_total(:,7)));

        if parameters.viz_enabled
            frames_to_show = 1:min(parameters.viz_max_frames, parameters.frames);
            figure;
            viz.plotImage(windowed_raw_data, frames_to_show, 'hot', ...
                'Drift-uncorrected Image and Locs (legacy)', 'max');
            hold on;
            viz.overlayLocs(loc_total, frames_to_show, true);
        end

        % Pack output（detected_total 不可用，置空）
        result.detected_total  = [];
        result.super_loc_total = super_loc_total;
        result.loc_total       = loc_total;
        result.n_raw           = n_filt;   % 无法追溯原始数目，设为与过滤后相同
        result.n_filtered      = n_filt;

        % 存入 cache 供后续直接加载
        if ~isempty(options.cacheFile)
            save(options.cacheFile, 'result', '-v7.3');
            fprintf('[detectAndLocalize] Cached (from legacy): %s\n', options.cacheFile);
        end
        return
    end

    % ================================================================
    %  Detection（粗定位）
    % ================================================================
    fprintf('[detectAndLocalize] Running detection ...\n');
    detected_total = detect.findMaxima(windowed_raw_data, parameters.k_sigma);

    if parameters.viz_enabled
        frames_to_show = 1:min(parameters.viz_max_frames, parameters.frames);
        figure;
        viz.plotImage(windowed_raw_data, frames_to_show, 'hot', ...
            'Drift-uncorrected Image and Detections');
        hold on;
        viz.overlayLocs(detected_total, frames_to_show);
    end

    % ================================================================
    %  MLE Localization（精定位 + 坏点过滤）
    % ================================================================
    fprintf('[detectAndLocalize] Running MLE localization ...\n');
    super_loc_raw   = localize.MLE1.locMolecules( ...
                        windowed_raw_data, detected_total, parameters.edge);
    super_loc_total = localize.filterBadLocs(super_loc_raw);
    loc_total       = super_loc_total(:, [1, 2, 7]);

    n_raw  = size(super_loc_raw, 1);
    n_filt = size(super_loc_total, 1);
    fprintf('[detectAndLocalize] Localizations: %d → %d (%d removed)\n', ...
        n_raw, n_filt, n_raw - n_filt);

    if parameters.viz_enabled
        frames_to_show = 1:min(parameters.viz_max_frames, parameters.frames);
        figure;
        viz.plotImage(windowed_raw_data, frames_to_show, 'hot', ...
            'Drift-uncorrected Image and Locs', 'max');
        hold on;
        viz.overlayLocs(loc_total, frames_to_show, true);
    end

    % ================================================================
    %  Pack output
    % ================================================================
    result.detected_total  = detected_total;
    result.super_loc_total = super_loc_total;
    result.loc_total       = loc_total;
    result.n_raw           = n_raw;
    result.n_filtered      = n_filt;

    % ================================================================
    %  Cache — 保存
    % ================================================================
    if ~isempty(options.cacheFile)
        save(options.cacheFile, 'result', '-v7.3');
        fprintf('[detectAndLocalize] Cached: %s\n', options.cacheFile);
    end
end