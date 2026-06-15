%% buildModel.m — 手动配置极性 & 粒径，汇总 metrics 并可视化
clear; clc;

%% ===== 手动配置区（修改这里） =====
    % 极性: 'plus_to_minus' 表示 +1→-1, 'minus_to_plus' 表示 -1→+1
    polarity =  'plus_to_minus';

    % 粒径选择: 从 [1, 2, 3, 4] 中选 2~4 个
    %   1 = 2.5 nm,  2 = 6 nm,  3 = 10 nm,  4 = small
    size_idx = [1, 2, 3, 4];

    % ===== 关键词 =====
    switch polarity
        case 'plus_to_minus'
            pol_keyword = 'first no E 1V and -1V';
            pol_label   = '+1 to -1';
        case 'minus_to_plus'
            pol_keyword = 'first no E -1V and 1V';
            pol_label   = '-1 to +1';
        otherwise
            error('无效极性，请使用 ''plus_to_minus'' 或 ''minus_to_plus''');
    end
    fprintf('极性: %s\n', pol_label);
    fprintf('粒径: %s\n\n', strjoin(arrayfun(@num2str, size_idx, 'Uni', false), ', '));

    % ===== 粒径 → 目录映射 =====
    size_info = struct();
    size_info(1).tag   = 'nm2p5';
    size_info(1).label = '2.5 nm';
    size_info(1).color = [0.2  0.6  0.8];
    size_info(1).dirs  = {
        '/Volumes/SMILeSSD/UCNP charge data/20251028 98Yb2Er@Lu-2.5nm';
    };

    size_info(2).tag   = 'nm6';
    size_info(2).label = '6 nm';
    size_info(2).color = [0.9  0.4  0.4];
    size_info(2).dirs  = {
        '/Volumes/SMILeSSD/UCNP charge data/20251028 98Yb2Er@Lu-6nm';
    };

    size_info(3).tag   = 'nm10';
    size_info(3).label = '10 nm';
    size_info(3).color = [0.4  0.8  0.4];
    size_info(3).dirs  = {
        '/Volumes/SMILeSSD/UCNP charge data/20251028 98Yb2Er@Lu-10nm';
        '/Volumes/SMILeSSD/UCNP charge data/20251029 98Yb2Er@Lu-10nm';
    };

    size_info(4).tag   = 'small';
    size_info(4).label = 'small';
    size_info(4).color = [0.7  0.5  0.9];
    size_info(4).dirs  = {
        '/Volumes/SMILeSSD/UCNP charge data/20251222+1223 small-98Yb2Er no NaYF4 300mA-laser/1222';
        '/Volumes/SMILeSSD/UCNP charge data/20251222+1223 small-98Yb2Er no NaYF4 300mA-laser/1223';
    };

%% ===== 搜索 locResult 文件 & 提取 metrics =====
    result_pattern = '**/*_chargeResult.mat';

    metrics_data = struct();

    for s = 1:numel(size_idx)
        idx = size_idx(s);
        tag = size_info(idx).tag;
        
        init_bright = [];
        norm_delta  = [];
        bd_ratio    = [];
        
        dirs = size_info(idx).dirs;
        fprintf('--- 搜索 %s ---\n', size_info(idx).label);
        for r = 1:numel(dirs)
            if ~isfolder(dirs{r})
                fprintf('  [MISS] 目录不存在: %s\n', dirs{r});
                continue;
            end
            d = dir(fullfile(dirs{r}, result_pattern));
            for k = 1:numel(d)
                fpath = fullfile(d(k).folder, d(k).name);
                if startsWith(d(k).name, '._')
                    continue;
                end
                % 极性过滤
                if ~contains(fpath, pol_keyword)
                    continue;
                end
                
                try
                    S = load(fpath, 'fileResult');
                    m = S.fileResult.metrics;
                    init_bright = [init_bright; m.initial_brightness(:)]; %#ok<AGROW>
                    norm_delta  = [norm_delta;  m.normed_delta(:)];       %#ok<AGROW>
                    bd_ratio    = [bd_ratio;    m.bright_dark_ratio(:)];  %#ok<AGROW>
                    fprintf('  [OK] %s  (N_em = %d)\n', d(k).name, numel(m.normed_delta));
                catch ME
                    fprintf('  [SKIP] %s : %s\n', d(k).name, ME.message);
                end
            end
        end
        
        metrics_data.(tag).initial_brightness = init_bright;
        metrics_data.(tag).normed_delta       = norm_delta;
        metrics_data.(tag).bright_dark_ratio  = bd_ratio;
        
        fprintf('>>> %s 总计: %d 个发光体\n\n', size_info(idx).label, numel(norm_delta));
    end

%% ===== 可视化：三个 metrics 各画一张叠加直方图 =====
    metric_fields  = {'normed_delta', 'initial_brightness', 'bright_dark_ratio'};
    metric_titles  = {
        ['Normalized \Delta brightness  (' pol_label ')'];
        ['Initial brightness  (' pol_label ')'];
        ['Bright/Dark ratio  (' pol_label ')'];
    };
    metric_xlabels = {
        '\Delta brightness / initial brightness';
        'Initial brightness (counts)';
        'Bright-dark ratio';
    };

    nbins = 50;

    for mi = 1:3
        field = metric_fields{mi};
        
        figure('Position', [50 + (mi-1)*50, 100, 900, 500]);
        hold on;
        
        all_data = [];
        dsets  = {};
        labels = {};
        colors = [];
        
        for s = 1:numel(size_idx)
            idx = size_idx(s);
            tag = size_info(idx).tag;
            d_raw = metrics_data.(tag).(field);
            d_clean = d_raw(isfinite(d_raw));
            
            dsets{end+1}  = d_clean; %#ok<SAGROW>
            labels{end+1} = size_info(idx).label; %#ok<SAGROW>
            colors = [colors; size_info(idx).color]; %#ok<AGROW>
            all_data = [all_data; d_clean(:)]; %#ok<AGROW>
        end
        
        if isempty(all_data)
            title([metric_titles{mi} ' — 无数据']);
            hold off;
            continue;
        end
        
        edges = linspace(min(all_data), max(all_data), nbins + 1);
        xfit  = linspace(min(all_data), max(all_data), 500);
        
        h_objs = gobjects(numel(dsets), 1);
        for s = 1:numel(dsets)
            d = dsets{s};
            disp_name = sprintf('%s  (N=%d, \\mu=%.3f, \\sigma=%.3f)', ...
                labels{s}, numel(d), mean(d), std(d));
            h_objs(s) = histogram(d, edges, ...
                'FaceColor', colors(s,:), 'EdgeColor', 'none', ...
                'FaceAlpha', 0.35, 'Normalization', 'probability', ...
                'DisplayName', disp_name);
        end
        
        % KDE 曲线
        for s = 1:numel(dsets)
            d = dsets{s};
            if numel(d) > 2
                [fk, xk] = ksdensity(d, xfit);
                plot(xk, fk / max(fk) * max(h_objs(s).Values), '-', ...
                    'Color', colors(s,:)*0.7, 'LineWidth', 2.2, 'HandleVisibility', 'off');
            end
        end
        
        xlabel(metric_xlabels{mi});
        ylabel('Probability');
        title(metric_titles{mi});
        legend('Location', 'best', 'FontSize', 9);
        box on;
        hold off;
    end

%% ===== 计算 max-mean 和 mean-min =====

    for s = 1:numel(size_idx)
        idx = size_idx(s);
        tag = size_info(idx).tag;
        
        IB    = metrics_data.(tag).initial_brightness;
        nd    = metrics_data.(tag).normed_delta;
        ratio = metrics_data.(tag).bright_dark_ratio;
        
        valid = isfinite(IB) & isfinite(nd) & isfinite(ratio) ;
        
        IB_v    = IB(valid);
        nd_v    = nd(valid);
        ratio_v = ratio(valid);
        
        mean_minus_min = nd_v .* IB_v ./ (ratio_v + 1);            % IB - D
        max_minus_mean = ratio_v .* nd_v .* IB_v ./ (ratio_v + 1); % B - IB
        
        metrics_data.(tag).max_minus_mean = max_minus_mean;
        metrics_data.(tag).mean_minus_min = mean_minus_min;
        
        fprintf('[%s] max-mean: μ=%.2f, σ=%.2f | mean-min: μ=%.2f, σ=%.2f  (N_valid=%d)\n', ...
        size_info(idx).label, ...
        mean(max_minus_mean), std(max_minus_mean), ...
        mean(mean_minus_min), std(mean_minus_min), sum(valid));
    end

%% ===== 可视化: max-mean 和 mean-min 叠加直方图 =====
    derived_fields  = {'max_minus_mean', 'mean_minus_min'};
    derived_titles  = {
        ['Max − Mean  (' pol_label ')'];
        ['Mean − Min  (' pol_label ')'];
    };
    derived_xlabels = {
        'Max − Mean (counts)';
        'Mean − Min (counts)';
    };

    for mi = 1:2
        field = derived_fields{mi};
        
        figure('Position', [100 + (mi-1)*60, 650, 900, 500]);
        hold on;
        
        all_data = [];
        dsets    = {};
        labels   = {};
        colors   = [];
        
        for s = 1:numel(size_idx)
            idx = size_idx(s);
            tag = size_info(idx).tag;
            
            if ~isfield(metrics_data.(tag), field)
                continue;
            end
            d_raw   = metrics_data.(tag).(field);
            d_clean = d_raw(isfinite(d_raw));
            
            % 可选: 去除极端离群值 (超出 ±5σ)
            mu_d = mean(d_clean);
            sg_d = std(d_clean);
            d_clean = d_clean(abs(d_clean - mu_d) < 5*sg_d);
            
            dsets{end+1}  = d_clean; %#ok<SAGROW>
            labels{end+1} = size_info(idx).label; %#ok<SAGROW>
            colors = [colors; size_info(idx).color]; %#ok<AGROW>
            all_data = [all_data; d_clean(:)]; %#ok<AGROW>
        end
        
        if isempty(all_data)
            title([derived_titles{mi} ' — 无数据']);
            hold off;
            continue;
        end
        
        edges = linspace(min(all_data), max(all_data), nbins + 1);
        xfit  = linspace(min(all_data), max(all_data), 500);
        
        h_objs = gobjects(numel(dsets), 1);
        for s = 1:numel(dsets)
            d = dsets{s};
            disp_name = sprintf('%s  (N=%d, \\mu=%.2f, \\sigma=%.2f)', ...
                labels{s}, numel(d), mean(d), std(d));
            h_objs(s) = histogram(d, edges, ...
                'FaceColor', colors(s,:), 'EdgeColor', 'none', ...
                'FaceAlpha', 0.35, 'Normalization', 'probability', ...
                'DisplayName', disp_name);
        end
        
        % KDE 曲线
        for s = 1:numel(dsets)
            d = dsets{s};
            if numel(d) > 2
                [fk, xk] = ksdensity(d, xfit);
                plot(xk, fk / max(fk) * max(h_objs(s).Values), '-', ...
                    'Color', colors(s,:)*0.7, 'LineWidth', 2.2, 'HandleVisibility', 'off');
            end
        end
        
        xlabel(derived_xlabels{mi}, 'FontSize', 11);
        ylabel('Probability', 'FontSize', 11);
        title(derived_titles{mi}, 'FontSize', 13);
        legend('Location', 'best', 'FontSize', 9);
        box on; grid on;
        hold off;
    end

%% ===== 额外: 归一化版本 (除以 IB，无量纲) =====
    derived_norm_fields  = {'max_minus_mean_norm', 'mean_minus_min_norm'};
    derived_norm_titles  = {
        ['(Max − Mean) / IB  (' pol_label ')'];
        ['(Mean − Min) / IB  (' pol_label ')'];
    };
    derived_norm_xlabels = {
        '(Max − Mean) / Initial Brightness';
        '(Mean − Min) / Initial Brightness';
    };

    % 计算归一化版本
    for s = 1:numel(size_idx)
        idx = size_idx(s);
        tag = size_info(idx).tag;
        
        IB    = metrics_data.(tag).initial_brightness;
        nd    = metrics_data.(tag).normed_delta;
        ratio = metrics_data.(tag).bright_dark_ratio;
        
        valid = isfinite(IB) & isfinite(nd) & isfinite(ratio)  & (IB > 0);
        
        nd_v    = nd(valid);
        ratio_v = ratio(valid);
    
        % ===== 正确的归一化公式 (除以 IB 后) =====
        metrics_data.(tag).max_minus_mean_norm = ratio_v .* nd_v ./ (ratio_v + 1);
        metrics_data.(tag).mean_minus_min_norm = nd_v ./ (ratio_v + 1);
    end

    % 可视化归一化版本
    for mi = 1:2
        field = derived_norm_fields{mi};
        
        figure('Position', [200 + (mi-1)*60, 650, 900, 500]);
        hold on;
        
        all_data = [];
        dsets    = {};
        labels   = {};
        colors   = [];
        
        for s = 1:numel(size_idx)
            idx = size_idx(s);
            tag = size_info(idx).tag;
            
            if ~isfield(metrics_data.(tag), field)
                continue;
            end
            d_raw   = metrics_data.(tag).(field);
            d_clean = d_raw(isfinite(d_raw));
            
            mu_d = mean(d_clean);
            sg_d = std(d_clean);
            d_clean = d_clean(abs(d_clean - mu_d) < 5*sg_d);
            
            dsets{end+1}  = d_clean; %#ok<SAGROW>
            labels{end+1} = size_info(idx).label; %#ok<SAGROW>
            colors = [colors; size_info(idx).color]; %#ok<AGROW>
            all_data = [all_data; d_clean(:)]; %#ok<AGROW>
        end
        
        if isempty(all_data)
            title([derived_norm_titles{mi} ' — 无数据']);
            hold off;
            continue;
        end
        
        edges = linspace(min(all_data), max(all_data), nbins + 1);
        xfit  = linspace(min(all_data), max(all_data), 500);
        
        h_objs = gobjects(numel(dsets), 1);
        for s = 1:numel(dsets)
            d = dsets{s};
            disp_name = sprintf('%s  (N=%d, \\mu=%.3f, \\sigma=%.3f)', ...
                labels{s}, numel(d), mean(d), std(d));
            h_objs(s) = histogram(d, edges, ...
                'FaceColor', colors(s,:), 'EdgeColor', 'none', ...
                'FaceAlpha', 0.35, 'Normalization', 'probability', ...
                'DisplayName', disp_name);
        end
        
        for s = 1:numel(dsets)
            d = dsets{s};
            if numel(d) > 2
                [fk, xk] = ksdensity(d, xfit);
                plot(xk, fk / max(fk) * max(h_objs(s).Values), '-', ...
                    'Color', colors(s,:)*0.7, 'LineWidth', 2.2, 'HandleVisibility', 'off');
            end
        end
        
        xlabel(derived_norm_xlabels{mi}, 'FontSize', 11);
        ylabel('Probability', 'FontSize', 11);
        title(derived_norm_titles{mi}, 'FontSize', 13);
        legend('Location', 'best', 'FontSize', 9);
        box on; grid on;
        hold off;
    end    

%% ===== 统计量汇总（normed_delta）=====
    fprintf('\n===== Distribution Summary: normed_delta (%s) =====\n', pol_label);
    for s = 1:numel(size_idx)
        idx = size_idx(s);
        tag = size_info(idx).tag;
        d = metrics_data.(tag).normed_delta;
        d = d(isfinite(d));
        fprintf('%8s : N=%4d  mean=%.4f  std=%.4f  median=%.4f  IQR=[%.4f, %.4f]\n', ...
            size_info(idx).label, numel(d), mean(d), std(d), median(d), ...
            prctile(d,25), prctile(d,75));
    end

    % Cohen's d
    fprintf('\n===== Pairwise Cohen''s d (normed_delta) =====\n');
    for i = 1:numel(size_idx)
        for j = i+1:numel(size_idx)
            di = metrics_data.(size_info(size_idx(i)).tag).normed_delta;
            dj = metrics_data.(size_info(size_idx(j)).tag).normed_delta;
            di = di(isfinite(di)); dj = dj(isfinite(dj));
            pooled_std = sqrt(((numel(di)-1)*var(di) + (numel(dj)-1)*var(dj)) / ...
                            max(numel(di)+numel(dj)-2, 1));
            cohen_d = abs(mean(di) - mean(dj)) / max(pooled_std, eps);
            fprintf('  %s vs %s : d = %.3f\n', ...
                size_info(size_idx(i)).label, size_info(size_idx(j)).label, cohen_d);
        end
    end

    % Overlap Coefficient
    fprintf('\n===== Pairwise Overlap Coefficient (normed_delta) =====\n');
    all_nd = [];
    for s = 1:numel(size_idx)
        tag = size_info(size_idx(s)).tag;
        d = metrics_data.(tag).normed_delta;
        all_nd = [all_nd; d(isfinite(d))]; %#ok<AGROW>
    end
    if ~isempty(all_nd)
        xcommon = linspace(min(all_nd), max(all_nd), 2000);
        kdes = cell(numel(size_idx), 1);
        for s = 1:numel(size_idx)
            tag = size_info(size_idx(s)).tag;
            d = metrics_data.(tag).normed_delta;
            d = d(isfinite(d));
            if numel(d) > 2
                kdes{s} = ksdensity(d, xcommon);
            else
                kdes{s} = zeros(size(xcommon));
            end
        end
        for i = 1:numel(size_idx)
            for j = i+1:numel(size_idx)
                ovl = trapz(xcommon, min(kdes{i}, kdes{j}));
                fprintf('  %s vs %s : OVL = %.3f\n', ...
                    size_info(size_idx(i)).label, size_info(size_idx(j)).label, ovl);
            end
        end
    end
    n_types = numel(size_idx);
    fprintf('\n===== buildModel 完成 =====\n');



%% ===== 1D Classification: normed_delta only K-fold =====
    % 自动适应 2~4 种粒径组合
    % 使用分层 K-Fold 交叉验证（测试集不参与模型构建）

    fprintf('\n========== 1D Gaussian Classification (normed_delta) ==========\n');

    n_types = numel(size_idx);

    % --- 全量训练: 从纯样品提取高斯参数（用于后续混合样品分类）---
    model1D = struct();
    for s = 1:n_types
        idx = size_idx(s);
        tag = size_info(idx).tag;
        dat = metrics_data.(tag).normed_delta;
        %dat = metrics_data.(tag).initial_brightness;
        dat = dat(isfinite(dat));
        
        model1D(s).mu    = mean(dat);
        model1D(s).sigma = max(std(dat), eps);
        model1D(s).tag   = tag;
        model1D(s).label = size_info(idx).label;
        model1D(s).color = size_info(idx).color;
        model1D(s).N_train = numel(dat);
        model1D(s).data  = dat;
        
        fprintf('  Type %d [%s]: mu=%.4f, sigma=%.4f, N=%d\n', ...
            s, model1D(s).label, model1D(s).mu, model1D(s).sigma, model1D(s).N_train);
    end

    % --- 分类参数 ---
    ratio_th_1D = 3;          % 置信度阈值
    priors_1D   = ones(1, n_types) / n_types;   % 均匀先验

    % --- 合并所有纯样品数据, 真实标签已知 ---
    all_train_data  = [];
    all_train_label = [];
    for s = 1:n_types
        dat = model1D(s).data;
        all_train_data  = [all_train_data;  dat(:)];          %#ok<AGROW>
        all_train_label = [all_train_label; s*ones(numel(dat),1)]; %#ok<AGROW>
    end
    N_train_total = numel(all_train_data);

    % =====================================================
    % ===== 分层 K-Fold 交叉验证 (Stratified K-Fold CV) =====
    % =====================================================
    K_fold = 5;
    rng_seed = 42;   % 固定随机种子, 保证可重复
    fprintf('\n--- Stratified %d-Fold Cross-Validation ---\n', K_fold);

    % 为每个类别分别分配 fold 索引 (分层: 保证每折各类比例一致)
    fold_assignment = zeros(N_train_total, 1);
    rng(rng_seed);
    for s = 1:n_types
        class_idx = find(all_train_label == s);
        n_s = numel(class_idx);
        perm = class_idx(randperm(n_s));   % 类内随机打乱
        for i = 1:n_s
            fold_assignment(perm(i)) = mod(i-1, K_fold) + 1;
        end
    end

    % 存储 CV 预测结果
    pred_label_1D_train = zeros(N_train_total, 1);   % CV 预测标签
    confidence_1D_train = zeros(N_train_total, 1);   % CV 置信度 (ratio)

    % 逐折验证
    for k = 1:K_fold
        test_mask  = (fold_assignment == k);
        train_mask = ~test_mask;
        
        % --- 在当前训练折上重新估计高斯参数 ---
        cv_model_mu    = zeros(1, n_types);
        cv_model_sigma = zeros(1, n_types);
        for s = 1:n_types
            dat_s = all_train_data(all_train_label == s & train_mask);
            cv_model_mu(s)    = mean(dat_s);
            cv_model_sigma(s) = max(std(dat_s), eps);
        end
        
        % --- 对当前测试折做预测 ---
        test_indices = find(test_mask);
        for i = 1:numel(test_indices)
            ii = test_indices(i);
            x = all_train_data(ii);
            
            if ~isfinite(x)
                pred_label_1D_train(ii) = 0;
                confidence_1D_train(ii) = 0;
                continue;
            end
            
            % 计算各类后验 (未归一化)
            p = zeros(1, n_types);
            for s = 1:n_types
                p(s) = priors_1D(s) * normpdf(x, cv_model_mu(s), cv_model_sigma(s));
            end
            
            % 排序, 取最高与次高的比值
            [p_sorted, idx_sorted] = sort(p, 'descend');
            ratio = p_sorted(1) / max(p_sorted(2), realmin);
            confidence_1D_train(ii) = ratio;
            
            if ratio < ratio_th_1D
                pred_label_1D_train(ii) = 0;   % ambiguous → unassigned
            else
                pred_label_1D_train(ii) = idx_sorted(1);
            end
        end
        
        % 报告每折结果
        fold_test_idx = find(test_mask);
        fold_correct = sum(pred_label_1D_train(fold_test_idx) == all_train_label(fold_test_idx) & ...
                           pred_label_1D_train(fold_test_idx) > 0);
        fold_assigned = sum(pred_label_1D_train(fold_test_idx) > 0);
        fprintf('    Fold %d: N_test=%d, Assigned=%d, Correct=%d\n', ...
            k, numel(fold_test_idx), fold_assigned, fold_correct);
    end

    % =====================================================
    % ===== 评价指标 (基于 CV held-out 预测) =====
    % =====================================================
    fprintf('\n--- 1D Model Evaluation (%d-Fold CV, test ≠ train) ---\n', K_fold);
    fprintf('  ratio_th = %.1f\n', ratio_th_1D);

    % 混淆矩阵 (不含 unassigned)
    assigned_mask = (pred_label_1D_train > 0);
    conf_mat_1D = zeros(n_types, n_types);
    for true_s = 1:n_types
        for pred_s = 1:n_types
            conf_mat_1D(true_s, pred_s) = sum( ...
                all_train_label == true_s & pred_label_1D_train == pred_s);
        end
    end

    % 准确率 (assigned only)
    correct_1D = sum(diag(conf_mat_1D));
    total_assigned_1D = sum(conf_mat_1D(:));
    accuracy_1D = correct_1D / max(total_assigned_1D, 1);
    unassigned_rate_1D = sum(~assigned_mask) / N_train_total;

    fprintf('  Overall Accuracy (assigned): %.2f%%\n', accuracy_1D*100);
    fprintf('  Unassigned Rate: %.2f%%\n', unassigned_rate_1D*100);
    fprintf('  Confusion Matrix (rows=true, cols=pred):\n');

    % 打印表头
    fprintf('%12s', '');
    for s = 1:n_types
        fprintf('%12s', model1D(s).label);
    end
    fprintf('%12s\n', 'Unassigned');
    for true_s = 1:n_types
        fprintf('%12s', model1D(true_s).label);
        for pred_s = 1:n_types
            fprintf('%12d', conf_mat_1D(true_s, pred_s));
        end
        n_unassigned_s = sum(all_train_label == true_s & pred_label_1D_train == 0);
        fprintf('%12d\n', n_unassigned_s);
    end

    % Per-class precision, recall, F1
    fprintf('\n  Per-class metrics:\n');
    fprintf('%12s %10s %10s %10s\n', 'Class', 'Precision', 'Recall', 'F1');
    precision_1D = zeros(n_types,1);
    recall_1D    = zeros(n_types,1);
    f1_1D        = zeros(n_types,1);
    for s = 1:n_types
        TP = conf_mat_1D(s,s);
        FP = sum(conf_mat_1D(:,s)) - TP;
        FN = sum(conf_mat_1D(s,:)) - TP;
        precision_1D(s) = TP / max(TP+FP, 1);
        recall_1D(s)    = TP / max(TP+FN, 1);
        f1_1D(s) = 2*precision_1D(s)*recall_1D(s) / max(precision_1D(s)+recall_1D(s), eps);
        fprintf('%12s %10.3f %10.3f %10.3f\n', ...
            model1D(s).label, precision_1D(s), recall_1D(s), f1_1D(s));
    end

    % --- 额外: 逐折准确率的均值 ± 标准差 ---
    fold_accuracies = zeros(K_fold, 1);
    for k = 1:K_fold
        fold_idx = find(fold_assignment == k);
        fold_assigned_mask = pred_label_1D_train(fold_idx) > 0;
        fold_correct = sum(pred_label_1D_train(fold_idx) == all_train_label(fold_idx) & fold_assigned_mask);
        fold_accuracies(k) = fold_correct / max(sum(fold_assigned_mask), 1);
    end
    fprintf('\n  CV Accuracy (per-fold): %.2f%% ± %.2f%%\n', ...
        mean(fold_accuracies)*100, std(fold_accuracies)*100);

    % Log-likelihood & BIC (基于全量模型, 衡量模型拟合质量)
    LL_1D = 0;
    for i = 1:N_train_total
        x = all_train_data(i);
        if ~isfinite(x), continue; end
        p_total = 0;
        for s = 1:n_types
            p_total = p_total + priors_1D(s) * normpdf(x, model1D(s).mu, model1D(s).sigma);
        end
        LL_1D = LL_1D + log(max(p_total, realmin));
    end
    BIC_1D = -2*LL_1D + (3*n_types - 1)*log(N_train_total);
    fprintf('\n  Log-Likelihood (full model): %.2f\n', LL_1D);
    fprintf('  BIC (full model): %.2f\n', BIC_1D);

    % --- 1D 可视化 ---
    figure('Position', [100 100 1200 500]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % 子图1: 训练数据分布 + 高斯拟合 (PDF)
    nexttile;
    hold on;
    all_dat_range = all_train_data(isfinite(all_train_data));
    edges_1D = linspace(min(all_dat_range), max(all_dat_range), 60);
    xfit_1D = linspace(edges_1D(1), edges_1D(end), 500);

    for s = 1:n_types
        histogram(model1D(s).data, edges_1D, ...
            'FaceColor', model1D(s).color, 'EdgeColor', 'none', ...
            'FaceAlpha', 0.5, 'Normalization', 'probability', ...
            'DisplayName', sprintf('%s (N=%d)', model1D(s).label, model1D(s).N_train));
        
        % 高斯拟合也要对应缩放: probability = pdf * bin_width
        bw_1D = edges_1D(2) - edges_1D(1);
        yfit = bw_1D * normpdf(xfit_1D, model1D(s).mu, model1D(s).sigma);
        plot(xfit_1D, yfit, '-', 'Color', model1D(s).color*0.7, 'LineWidth', 2.5, ...
            'HandleVisibility', 'off');
    end
    
    xlabel('\Delta brightness / initial brightness');
    ylabel('Probability');
    title(sprintf('1D Model Fit Curve — Full Training Data [%s]', pol_label));
    legend('Location', 'best');
    box on; hold off;

    % 子图2: 混淆矩阵热力图 (含 Unassigned 列)
    nexttile;

    % 构建 n_types × (n_types+1) 矩阵, 最后一列为 Unassigned
    conf_mat_1D_ext = zeros(n_types, n_types + 1);
    conf_mat_1D_ext(:, 1:n_types) = conf_mat_1D;
    for true_s = 1:n_types
        conf_mat_1D_ext(true_s, end) = sum(all_train_label == true_s & pred_label_1D_train == 0);
    end

    % 归一化: 分母为每个真实类别的总粒子数 (assigned + unassigned)
    row_totals = sum(conf_mat_1D_ext, 2);  % 每行总数 = 该类全部粒子
    conf_mat_1D_pct = conf_mat_1D_ext ./ max(row_totals, 1) * 100;

    imagesc(conf_mat_1D_pct);
    colormap(gca, parula);
    colorbar;
    caxis([0 100]);

    xlabels_ext = [{model1D.label}, {'Unassigned'}];
    set(gca, 'XTick', 1:(n_types+1), 'XTickLabel', xlabels_ext);
    set(gca, 'YTick', 1:n_types, 'YTickLabel', {model1D.label});
    xlabel('Predicted'); ylabel('True');
    title(sprintf('CV Confusion Matrix (%%) | Acc=%.1f%% ± %.1f%%', ...
        mean(fold_accuracies)*100, std(fold_accuracies)*100));

    for r = 1:n_types
        for c = 1:(n_types + 1)
            text(c, r, sprintf('%.1f%%', conf_mat_1D_pct(r,c)), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
    box on;

    sgtitle(sprintf('1D Gaussian Naive Bayes — Stratified %d-Fold CV', K_fold), 'FontSize', 14);


%% ===== 2D Classification: normed_delta + initial_brightness K-fold =====
    % 使用分层 K-Fold 交叉验证（测试集不参与模型构建）

    fprintf('\n\n========== 2D Gaussian Classification (normed_delta + initial_brightness) ==========\n');

    % --- 全量训练: 从纯样品提取 2D 高斯参数（用于后续混合样品分类）---
    model2D = struct();
    for s = 1:n_types
        idx = size_idx(s);
        tag = size_info(idx).tag;
        
        nd = metrics_data.(tag).normed_delta;
        ib = metrics_data.(tag).initial_brightness;
        
        % 对齐 & 过滤
        valid_s = isfinite(nd) & isfinite(ib);
        nd = nd(valid_s);
        ib = ib(valid_s);
        
        X_s = [nd(:), ib(:)];   % N x 2
        
        model2D(s).mu    = mean(X_s, 1);          % 1x2
        model2D(s).Sigma = cov(X_s);              % 2x2
        model2D(s).tag   = tag;
        model2D(s).label = size_info(idx).label;
        model2D(s).color = size_info(idx).color;
        model2D(s).N_train = size(X_s, 1);
        model2D(s).data  = X_s;
        
        % 正则化: 确保协方差矩阵正定
        [V, D] = eig(model2D(s).Sigma);
        D(D < 1e-6) = 1e-6;
        model2D(s).Sigma = V * D / V;
        
        fprintf('  Type %d [%s]: mu=[%.4f, %.1f], N=%d\n', ...
            s, model2D(s).label, model2D(s).mu(1), model2D(s).mu(2), model2D(s).N_train);
        fprintf('    Sigma = [%.4f, %.2f; %.2f, %.1f]\n', ...
            model2D(s).Sigma(1,1), model2D(s).Sigma(1,2), ...
            model2D(s).Sigma(2,1), model2D(s).Sigma(2,2));
    end

    % --- 分类参数 ---
    priors_2D   = ones(1, n_types) / n_types;
    ratio_th_2D  = 3;

    % --- 合并所有纯样品 2D 数据, 真实标签已知 ---
    all_train_2D = [];
    all_train_label_2D = [];
    for s = 1:n_types
        X_s = model2D(s).data;
        all_train_2D = [all_train_2D; X_s];                          %#ok<AGROW>
        all_train_label_2D = [all_train_label_2D; s*ones(size(X_s,1),1)]; %#ok<AGROW>
    end
    N_train_2D = size(all_train_2D, 1);

    % =====================================================
    % ===== 分层 K-Fold 交叉验证 (Stratified K-Fold CV) =====
    % =====================================================
    K_fold_2D = 5;
    rng_seed_2D = 42;   % 固定随机种子, 保证可重复
    fprintf('\n--- Stratified %d-Fold Cross-Validation (2D) ---\n', K_fold_2D);

    % 为每个类别分别分配 fold 索引 (分层: 保证每折各类比例一致)
    fold_assignment_2D = zeros(N_train_2D, 1);
    rng(rng_seed_2D);
    for s = 1:n_types
        class_idx_s = find(all_train_label_2D == s);
        n_s = numel(class_idx_s);
        perm = class_idx_s(randperm(n_s));   % 类内随机打乱
        for i = 1:n_s
            fold_assignment_2D(perm(i)) = mod(i-1, K_fold_2D) + 1;
        end
    end

    % 存储 CV 预测结果
    pred_label_2D_train = zeros(N_train_2D, 1);   % CV 预测标签
    confidence_2D_train = zeros(N_train_2D, 1);   % CV 置信度 (ratio)
    posterior_2D_train  = zeros(N_train_2D, n_types);  % CV 后验概率

    % 逐折验证
    for k = 1:K_fold_2D
        test_mask  = (fold_assignment_2D == k);
        train_mask = ~test_mask;
        
        % --- 在当前训练折上重新估计 2D 高斯参数 ---
        cv_model_mu_2D    = zeros(n_types, 2);
        cv_model_Sigma_2D = zeros(2, 2, n_types);
        for s = 1:n_types
            dat_s = all_train_2D(all_train_label_2D == s & train_mask, :);
            cv_model_mu_2D(s, :) = mean(dat_s, 1);
            S = cov(dat_s);
            % 正则化
            [V, D_eig] = eig(S);
            D_eig(D_eig < 1e-6) = 1e-6;
            S = V * D_eig / V;
            cv_model_Sigma_2D(:,:,s) = S;
        end
        
        % --- 对当前测试折做预测 ---
        test_indices = find(test_mask);
        for i = 1:numel(test_indices)
            ii = test_indices(i);
            xi = all_train_2D(ii, :);   % 1x2
            
            if any(~isfinite(xi))
                pred_label_2D_train(ii) = 0;
                confidence_2D_train(ii) = 0;
                continue;
            end
            
            % 计算各类后验 (未归一化)
            p = zeros(1, n_types);
            for s = 1:n_types
                p(s) = priors_2D(s) * mvnpdf(xi, cv_model_mu_2D(s,:), cv_model_Sigma_2D(:,:,s));
            end
            
            % 归一化后验
            p_sum = sum(p);
            if p_sum > 0
                posterior_2D_train(ii,:) = p / p_sum;
            end
            
            % 排序, 取最高与次高的比值
            [p_sorted, idx_sorted] = sort(p, 'descend');
            ratio = p_sorted(1) / max(p_sorted(2), realmin);
            confidence_2D_train(ii) = ratio;
            
            if ratio < ratio_th_2D
                pred_label_2D_train(ii) = 0;   % ambiguous → unassigned
            else
                pred_label_2D_train(ii) = idx_sorted(1);
            end
        end
        
        % 报告每折结果
        fold_test_idx = find(test_mask);
        fold_correct = sum(pred_label_2D_train(fold_test_idx) == all_train_label_2D(fold_test_idx) & ...
                           pred_label_2D_train(fold_test_idx) > 0);
        fold_assigned = sum(pred_label_2D_train(fold_test_idx) > 0);
        fprintf('    Fold %d: N_test=%d, Assigned=%d, Correct=%d\n', ...
            k, numel(fold_test_idx), fold_assigned, fold_correct);
    end

    % =====================================================
    % ===== 评价指标 (基于 CV held-out 预测) =====
    % =====================================================
    fprintf('\n--- 2D Model Evaluation (%d-Fold CV, test ≠ train) ---\n', K_fold_2D);
    fprintf('  ratio_th = %.1f\n', ratio_th_2D);

    % 混淆矩阵 (不含 unassigned)
    conf_mat_2D = zeros(n_types, n_types);
    for true_s = 1:n_types
        for pred_s = 1:n_types
            conf_mat_2D(true_s, pred_s) = sum( ...
                all_train_label_2D == true_s & pred_label_2D_train == pred_s);
        end
    end

    correct_2D = sum(diag(conf_mat_2D));
    total_assigned_2D = sum(conf_mat_2D(:));
    accuracy_2D = correct_2D / max(total_assigned_2D, 1);
    unassigned_rate_2D = sum(pred_label_2D_train == 0) / N_train_2D;

    fprintf('  Overall Accuracy (assigned): %.2f%%\n', accuracy_2D*100);
    fprintf('  Unassigned Rate: %.2f%%\n', unassigned_rate_2D*100);
    fprintf('  Confusion Matrix (rows=true, cols=pred):\n');

    fprintf('%12s', '');
    for s = 1:n_types
        fprintf('%12s', model2D(s).label);
    end
    fprintf('%12s\n', 'Unassigned');
    for true_s = 1:n_types
        fprintf('%12s', model2D(true_s).label);
        for pred_s = 1:n_types
            fprintf('%12d', conf_mat_2D(true_s, pred_s));
        end
        n_unassigned_s = sum(all_train_label_2D == true_s & pred_label_2D_train == 0);
        fprintf('%12d\n', n_unassigned_s);
    end

    fprintf('\n  Per-class metrics:\n');
    fprintf('%12s %10s %10s %10s\n', 'Class', 'Precision', 'Recall', 'F1');
    precision_2D = zeros(n_types,1);
    recall_2D    = zeros(n_types,1);
    f1_2D        = zeros(n_types,1);
    for s = 1:n_types
        TP = conf_mat_2D(s,s);
        FP = sum(conf_mat_2D(:,s)) - TP;
        FN = sum(conf_mat_2D(s,:)) - TP;
        precision_2D(s) = TP / max(TP+FP, 1);
        recall_2D(s)    = TP / max(TP+FN, 1);
        f1_2D(s) = 2*precision_2D(s)*recall_2D(s) / max(precision_2D(s)+recall_2D(s), eps);
        fprintf('%12s %10.3f %10.3f %10.3f\n', ...
            model2D(s).label, precision_2D(s), recall_2D(s), f1_2D(s));
    end

    % --- 逐折准确率的均值 ± 标准差 ---
    fold_accuracies_2D = zeros(K_fold_2D, 1);
    for k = 1:K_fold_2D
        fold_idx = find(fold_assignment_2D == k);
        fold_assigned_mask = pred_label_2D_train(fold_idx) > 0;
        fold_correct = sum(pred_label_2D_train(fold_idx) == all_train_label_2D(fold_idx) & fold_assigned_mask);
        fold_accuracies_2D(k) = fold_correct / max(sum(fold_assigned_mask), 1);
    end
    fprintf('\n  CV Accuracy (per-fold): %.2f%% ± %.2f%%\n', ...
        mean(fold_accuracies_2D)*100, std(fold_accuracies_2D)*100);

    % Log-likelihood & BIC (基于全量模型, 衡量模型拟合质量)
    LL_2D = 0;
    for i = 1:N_train_2D
        xi = all_train_2D(i,:);
        if any(~isfinite(xi)), continue; end
        p_total = 0;
        for s = 1:n_types
            p_total = p_total + priors_2D(s) * mvnpdf(xi, model2D(s).mu, model2D(s).Sigma);
        end
        LL_2D = LL_2D + log(max(p_total, realmin));
    end
    % 2D高斯: 每类5参数 (2 mu + 3 unique in Sigma) + n_types-1 priors
    n_params_2D = n_types * 5 + (n_types - 1);
    BIC_2D = -2*LL_2D + n_params_2D * log(N_train_2D);
    fprintf('\n  Log-Likelihood (full model): %.2f\n', LL_2D);
    fprintf('  BIC (full model): %.2f\n', BIC_2D);

    % --- 1D vs 2D 对比 ---
    fprintf('\n--- Model Comparison ---\n');
    fprintf('  %-20s %12s %12s\n', '', '1D', '2D');
    fprintf('  %-20s %11.2f%% %11.2f%%\n', 'CV Accuracy', accuracy_1D*100, accuracy_2D*100);
    fprintf('  %-20s %11.2f%% %11.2f%%\n', 'CV Acc ± std', ...
        mean(fold_accuracies)*100, mean(fold_accuracies_2D)*100);
    fprintf('  %-20s %11.2f%% %11.2f%%\n', 'Unassigned Rate', unassigned_rate_1D*100, unassigned_rate_2D*100);
    fprintf('  %-20s %12.1f %12.1f\n', 'Log-Likelihood', LL_1D, LL_2D);
    fprintf('  %-20s %12.1f %12.1f\n', 'BIC', BIC_1D, BIC_2D);

    % --- 2D 可视化 ---
    figure('Position', [50 50 1400 900]);
    tl2 = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % 颜色设置
    clr_fill = [0.75 0.75 0.75];   % unassigned
    clr_line_un = [0.5 0.5 0.5];
    for s = 1:n_types
        clr_fill = [clr_fill; model2D(s).color];           %#ok<AGROW>
        clr_line_un = [clr_line_un; model2D(s).color*0.6]; %#ok<AGROW>
    end
    names_all = [{'Unassigned'}, {model2D.label}];

    % 分组索引 (基于 CV 预测)
    lab_2D = pred_label_2D_train;
    x_2D = all_train_2D(:,1);   % normed_delta
    y_2D = all_train_2D(:,2);   % initial_brightness

    idx_groups = cell(n_types+1, 1);
    idx_groups{1} = (lab_2D == 0);
    for s = 1:n_types
        idx_groups{s+1} = (lab_2D == s);
    end

    % ----- Tile 1: Normed delta marginal (按真实标签, probability + Gaussian fit) -----
    nexttile(1);
    hold on;
    edges_nd = linspace(min(x_2D), max(x_2D), 50);
    bw_nd = edges_nd(2) - edges_nd(1);
    xfit_nd = linspace(edges_nd(1), edges_nd(end), 500);

    for s = 1:n_types
        d_s = x_2D(all_train_label_2D == s);
        histogram(d_s, edges_nd, ...
            'FaceColor', model2D(s).color, 'EdgeColor', 'none', ...
            'FaceAlpha', 0.5, 'Normalization', 'probability', ...
            'DisplayName', sprintf('%s (N=%d)', model2D(s).label, numel(d_s)));
        % 高斯拟合
        mu_nd_s = model2D(s).mu(1);
        sigma_nd_s = sqrt(model2D(s).Sigma(1,1));
        yfit_nd = bw_nd * normpdf(xfit_nd, mu_nd_s, sigma_nd_s);
        plot(xfit_nd, yfit_nd, '-', 'Color', model2D(s).color*0.6, 'LineWidth', 2.2, ...
            'HandleVisibility', 'off');
    end
    xlabel('\Delta brightness / initial brightness');
    ylabel('Probability');
    title('Normed \Delta Brightness (marginal)');
    legend('Location', 'best', 'FontSize', 8); box on; hold off;

    % ----- Tile 2: Initial brightness marginal (按真实标签) -----
    nexttile(2);
    hold on;
    edges_ib = linspace(min(y_2D), max(y_2D), 50);
    bw_ib = edges_ib(2) - edges_ib(1);
    xfit_ib = linspace(edges_ib(1), edges_ib(end), 500);

    for s = 1:n_types
        d_s = y_2D(all_train_label_2D == s);
        histogram(d_s, edges_ib, ...
            'FaceColor', model2D(s).color, 'EdgeColor', 'none', ...
            'FaceAlpha', 0.5, 'Normalization', 'probability', ...
            'DisplayName', sprintf('%s (N=%d)', model2D(s).label, numel(d_s)));
        mu_ib_s = model2D(s).mu(2);
        sigma_ib_s = sqrt(model2D(s).Sigma(2,2));
        yfit_ib = bw_ib * normpdf(xfit_ib, mu_ib_s, sigma_ib_s);
        plot(xfit_ib, yfit_ib, '-', 'Color', model2D(s).color*0.6, 'LineWidth', 2.2, ...
            'HandleVisibility', 'off');
    end
    xlabel('Initial Brightness (counts)');
    ylabel('Probability');
    title('Initial Brightness (marginal)');
    legend('Location', 'best', 'FontSize', 8); box on; hold off;

    % ----- Tile 3: Summary text -----
    nexttile(3);
    axis off;
    str = {sprintf('\\bf 2D Classification Summary (%d-Fold CV)', K_fold_2D), ''};
    str{end+1} = sprintf('  Polarity: %s', pol_label);
    str{end+1} = sprintf('  ratio\\_th = %.1f', ratio_th_2D);
    str{end+1} = '';
    for k = 0:n_types
        n_k = sum(lab_2D == k);
        str{end+1} = sprintf('  %s :  N = %d  (%.1f%%)', ...
            names_all{k+1}, n_k, 100*n_k/N_train_2D);             %#ok<AGROW>
    end
    str{end+1} = '';
    str{end+1} = sprintf('  CV Accuracy: %.1f%% \\pm %.1f%%', ...
        mean(fold_accuracies_2D)*100, std(fold_accuracies_2D)*100);
    str{end+1} = sprintf('  Unassigned: %.1f%%', unassigned_rate_2D*100);
    str{end+1} = sprintf('  BIC: %.1f', BIC_2D);
    text(0.05, 0.95, str, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 11, ...
        'FontName', 'FixedWidth', 'Interpreter', 'tex');

    % ----- Tile 4: Joint scatter (true label) -----
    nexttile(4);
    hold on;
    true_colors = zeros(N_train_2D, 3);
    for s = 1:n_types
        mask_s = (all_train_label_2D == s);
        true_colors(mask_s, :) = repmat(model2D(s).color, sum(mask_s), 1);
    end
    scatter(x_2D, y_2D, 15, true_colors, 'filled', 'MarkerFaceAlpha', 0.4);

    % 绘制 2D 高斯等高线 (全量模型)
    theta_ell = linspace(0, 2*pi, 100);
    for s = 1:n_types
        mu_s = model2D(s).mu;
        [V, D_eig] = eig(model2D(s).Sigma);
        % 1-sigma 和 2-sigma 椭圆
        for n_sig = [1, 2]
            ell = n_sig * V * sqrt(D_eig) * [cos(theta_ell); sin(theta_ell)];
            plot(mu_s(1) + ell(1,:), mu_s(2) + ell(2,:), '-', ...
                'Color', model2D(s).color*0.6, 'LineWidth', 1.5);
        end
        plot(mu_s(1), mu_s(2), 'x', 'Color', model2D(s).color*0.5, ...
            'MarkerSize', 12, 'LineWidth', 3);
    end
    xlabel('\Delta brightness / initial brightness');
    ylabel('Initial Brightness (counts)');
    title('Joint Distribution (true labels + model ellipses)');
    grid on; box on; hold off;

    % ----- Tile 5: Joint scatter (CV predicted label) -----
    nexttile(5);
    hold on;
    % Unassigned first
    scatter(x_2D(idx_groups{1}), y_2D(idx_groups{1}), 10, clr_fill(1,:), 'filled', ...
        'MarkerFaceAlpha', 0.2, 'DisplayName', 'Unassigned');
    for s = 1:n_types
        scatter(x_2D(idx_groups{s+1}), y_2D(idx_groups{s+1}), 20, clr_fill(s+1,:), 'filled', ...
            'MarkerFaceAlpha', 0.5, 'DisplayName', names_all{s+1});
    end

    % 绘制决策边界 (基于全量模型, 网格法)
    n_grid = 200;
    xg = linspace(min(x_2D), max(x_2D), n_grid);
    yg = linspace(min(y_2D), max(y_2D), n_grid);
    [Xg, Yg] = meshgrid(xg, yg);
    Zg = zeros(n_grid, n_grid);

    for gi = 1:n_grid
        for gj = 1:n_grid
            xi_g = [Xg(gi,gj), Yg(gi,gj)];
            p_g = zeros(1, n_types);
            for s = 1:n_types
                p_g(s) = priors_2D(s) * mvnpdf(xi_g, model2D(s).mu, model2D(s).Sigma);
            end
            [~, Zg(gi,gj)] = max(p_g);
        end
    end
    contour(Xg, Yg, Zg, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    xlabel('\Delta brightness / initial brightness');
    ylabel('Initial Brightness (counts)');
    title('2D CV Predictions + Decision Boundaries');
    legend('Location', 'best');
    grid on; box on; hold off;

    % ----- Tile 6: Confusion matrix heatmap (CV, 含 Unassigned 列) -----
    nexttile(6);
    
    % 构建 n_types × (n_types+1) 矩阵, 最后一列为 Unassigned
    conf_mat_2D_ext = zeros(n_types, n_types + 1);
    conf_mat_2D_ext(:, 1:n_types) = conf_mat_2D;
    for true_s = 1:n_types
        conf_mat_2D_ext(true_s, end) = sum(all_train_label_2D == true_s & pred_label_2D_train == 0);
    end
    
    % 归一化: 分母为每个真实类别的总粒子数 (assigned + unassigned)
    row_totals_2D = sum(conf_mat_2D_ext, 2);
    conf_mat_2D_pct = conf_mat_2D_ext ./ max(row_totals_2D, 1) * 100;
    
    imagesc(conf_mat_2D_pct);
    colormap(gca, parula);
    colorbar;
    caxis([0 100]);
    
    xlabels_ext_2D = [{model2D.label}, {'Unassigned'}];
    set(gca, 'XTick', 1:(n_types+1), 'XTickLabel', xlabels_ext_2D);
    set(gca, 'YTick', 1:n_types, 'YTickLabel', {model2D.label});
    xlabel('Predicted'); ylabel('True');
    title(sprintf('CV Confusion Matrix (%%) | Acc=%.1f%% \\pm %.1f%%', ...
        mean(fold_accuracies_2D)*100, std(fold_accuracies_2D)*100));
    
    for r = 1:n_types
        for c = 1:(n_types + 1)
            text(c, r, sprintf('%.1f%%', conf_mat_2D_pct(r,c)), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
    box on;

    sgtitle(sprintf('2D Multivariate Gaussian — Stratified %d-Fold CV [%s]', K_fold_2D, pol_label), 'FontSize', 14);

%% 加载混合样数据

    %file_name = '/Volumes/SMILeSSD/UCNP charge data/20260317 98Yb2Er mix 2.5nm 6nm (shell) 300mA-laser/7 98Yb2Er-ZJ first no E 1V and -1V every5s 5cycles then no E.sif';
    %file_name = '/Volumes/SMILeSSD/UCNP charge data/20260317 98Yb2Er mix 2.5nm 6nm (shell) 300mA-laser/4 98Yb2Er-ZJ first no E -1V and 1V every5s 5cycles then no E.sif';
    
    %file_name = '/Volumes/SMILeSSD/UCNP charge data/20260401 Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and 98Yb2Er@Lu(without Y) 300mA-laser/3mix Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu)/5 98Yb2Er-ZJ 3mix first no E 1V and -1V every5s 5cycles then no E.sif';
    %file_name = '/Volumes/SMILeSSD/UCNP charge data/20260401 Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and 98Yb2Er@Lu(without Y) 300mA-laser/3mix Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu)/1 98Yb2Er-ZJ 3mix first no E -1V and 1V every5s 5cycles then no E.sif';

    file_name = '/Volumes/SMILeSSD/UCNP charge data/20260401 Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and 98Yb2Er@Lu(without Y) 300mA-laser/4mix Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and small98Yb2Er@Lu(without Y)/4 98Yb2Er-ZJ 4mix first no E 1V and -1V every5s 5cycles then no E.sif';
    %file_name = '/Volumes/SMILeSSD/UCNP charge data/20260401 Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and 98Yb2Er@Lu(without Y) 300mA-laser/4mix Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and small98Yb2Er@Lu(without Y)/1 98Yb2Er-ZJ 4mix first no E -1V and 1V every5s 5cycles then no E.sif';
    [file_dir, file_base, ~] = fileparts(file_name);

    cacheFile = fullfile(file_dir, [file_base, '_chargeResult.mat']);
    S = load(cacheFile, 'fileResult');
    result = S.fileResult;
    fprintf('Loaded from cache: %s\n', cacheFile);
    parameters = result.parameters;
    [raw_data, ex_time, ~] = io.readSIFData(parameters.file_name);
    windowed_raw_data = raw_data(parameters.row_range, parameters.col_range, :);
    clear raw_data;
    fprintf('Windowed data: %d × %d × %d\n', ...
        size(windowed_raw_data,1), size(windowed_raw_data,2), parameters.frames);
    [file_dir, file_base, ~] = fileparts(parameters.file_name);
    cache_loc = fullfile(file_dir, [file_base, '_locResult.mat']);
    locResult = pipeline.detectAndLocalize(windowed_raw_data, parameters, ...
        'cacheFile', cache_loc);
    [data, driftResult] = pipeline.correctDrift( ...
        windowed_raw_data, locResult, parameters);
    clear windowed_raw_data;
    emitterResult = result.emitterResult;
    metrics = result.metrics;
    normed_delta_mix = metrics.normed_delta(:);
    initial_brightness_mix = metrics.initial_brightness(:);
    N_mix = numel(normed_delta_mix);
    pos_mix = emitterResult.stats.pos_mean_px;
    fprintf('  Mixed sample: N = %d emitters\n', N_mix);


%% ===== 分类 (混合样品) =====
    fprintf('\n--- 1D Classification on Mixed Sample ---\n');

    label_1D_mix = zeros(N_mix, 1);
    confidence_1D_mix = zeros(N_mix, 1);
    posterior_1D_mix = zeros(N_mix, n_types);

    for i = 1:N_mix
        x = normed_delta_mix(i);
        
        if ~isfinite(x)
            label_1D_mix(i) = 0;
            continue;
        end
        
        p = zeros(1, n_types);
        for s = 1:n_types
            p(s) = priors_1D(s) * normpdf(x, model1D(s).mu, model1D(s).sigma);
        end
        
        % 归一化后验
        p_sum = sum(p);
        if p_sum > 0
            posterior_1D_mix(i,:) = p / p_sum;
        end
        
        [p_sorted, idx_sorted] = sort(p, 'descend');
        ratio = p_sorted(1) / max(p_sorted(2), realmin);
        confidence_1D_mix(i) = ratio;
        
        if ratio < ratio_th_1D
            label_1D_mix(i) = 0;
        else
            label_1D_mix(i) = idx_sorted(1);
        end
    end

    % 统计
    fprintf('  ratio_th = %.1f\n', ratio_th_1D);
    fprintf('  Unassigned: N = %d  (%.1f%%)\n', sum(label_1D_mix==0), 100*mean(label_1D_mix==0));
    for s = 1:n_types
        fprintf('  %8s:   N = %d  (%.1f%%)\n', model1D(s).label, ...
            sum(label_1D_mix==s), 100*mean(label_1D_mix==s));
    end

    % ===== 2D 分类 (混合样品) =====
    fprintf('\n--- 2D Classification on Mixed Sample ---\n');

    label_2D_mix = zeros(N_mix, 1);
    confidence_2D_mix = zeros(N_mix, 1);
    posterior_2D_mix = zeros(N_mix, n_types);

    for i = 1:N_mix
        xi = [normed_delta_mix(i), initial_brightness_mix(i)];
        
        if any(~isfinite(xi))
            label_2D_mix(i) = 0;
            continue;
        end
        
        p = zeros(1, n_types);
        for s = 1:n_types
            p(s) = priors_2D(s) * mvnpdf(xi, model2D(s).mu, model2D(s).Sigma);
        end
        
        % 归一化后验
        p_sum = sum(p);
        if p_sum > 0
            posterior_2D_mix(i,:) = p / p_sum;
        end
        
        [p_sorted, idx_sorted] = sort(p, 'descend');
        ratio = p_sorted(1) / max(p_sorted(2), realmin);
        confidence_2D_mix(i) = ratio;
        
        if ratio < ratio_th_2D
            label_2D_mix(i) = 0;
        else
            label_2D_mix(i) = idx_sorted(1);
        end
    end

    % 统计
    fprintf('  ratio_th = %.1f\n', ratio_th_2D);
    fprintf('  Unassigned: N = %d  (%.1f%%)\n', sum(label_2D_mix==0), 100*mean(label_2D_mix==0));
    for s = 1:n_types
        fprintf('  %8s:   N = %d  (%.1f%%)\n', model2D(s).label, ...
            sum(label_2D_mix==s), 100*mean(label_2D_mix==s));
    end

%% ===== 1D vs 2D 比较 =====
    fprintf('\n--- 1D vs 2D Agreement (Mixed Sample) ---\n');
    agree_mix = (label_1D_mix == label_2D_mix);
    disagree_mix = ~agree_mix & (label_1D_mix > 0 | label_2D_mix > 0);
    fprintf('  Total emitters: %d\n', N_mix);
    fprintf('  Agreement:      %d (%.1f%%)\n', sum(agree_mix), 100*mean(agree_mix));
    fprintf('  Disagreement:   %d (%.1f%%)\n', sum(disagree_mix), 100*mean(disagree_mix));

    % 交叉表
    fprintf('\n  Cross-tabulation (1D rows × 2D cols):\n');
    fprintf('%12s', '1D\\2D');
    fprintf('%10s', 'Unasgn');
    for s = 1:n_types
        fprintf('%10s', model1D(s).label);
    end
    fprintf('\n');
    for r = 0:n_types
        if r == 0
            fprintf('%12s', 'Unasgn');
        else
            fprintf('%12s', model1D(r).label);
        end
        for c = 0:n_types
            fprintf('%10d', sum(label_1D_mix == r & label_2D_mix == c));
        end
        fprintf('\n');
    end

%% ===== 可视化 1: 图像标注 (1D & 2D 并排) =====
    figure('Position', [50 100 1400 600]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % — 1D 分类结果叠加在图像上 —
    nexttile;
    mean_img = mean(data(:,:,1:min(700, size(data,3))), 3);
    imagesc(mean_img); colormap(gca, gray); axis image; hold on;

    % Unassigned
    idx_0 = (label_1D_mix == 0);
    scatter(pos_mix(idx_0,2), pos_mix(idx_0,1), 80, ...
        [0.6 0.6 0.6], 'o', 'LineWidth', 0.8);

    h_leg = gobjects(n_types, 1);
    for s = 1:n_types
        idx_s = (label_1D_mix == s);
        h_leg(s) = scatter(pos_mix(idx_s,2), pos_mix(idx_s,1), 150, ...
            model1D(s).color, 'o', 'LineWidth', 1.5);
    end
    legend(h_leg, {model1D.label}, 'Location', 'best', 'FontSize', 11);
    title(sprintf('1D Classification | ratio_{th}=%.1f | Unassigned=%.1f%%', ...
        ratio_th_1D, 100*mean(label_1D_mix==0)));
    hold off;

    % — 2D 分类结果叠加在图像上 —
    nexttile;
    imagesc(mean_img); colormap(gca, gray); axis image; hold on;

    idx_0 = (label_2D_mix == 0);
    scatter(pos_mix(idx_0,2), pos_mix(idx_0,1), 80, ...
        [0.6 0.6 0.6], 'o', 'LineWidth', 0.8);

    h_leg = gobjects(n_types, 1);
    for s = 1:n_types
        idx_s = (label_2D_mix == s);
        h_leg(s) = scatter(pos_mix(idx_s,2), pos_mix(idx_s,1), 150, ...
            model2D(s).color, 'o', 'LineWidth', 1.5);
    end
    legend(h_leg, {model2D.label}, 'Location', 'best', 'FontSize', 11);
    title(sprintf('2D Classification | ratio_{th}=%.1f | Unassigned=%.1f%%', ...
        ratio_th_2D, 100*mean(label_2D_mix==0)));
    hold off;

    sgtitle('Mixed Sample Classification (Image Overlay)', 'FontSize', 14);

%% ===== 可视化 2: 1D 直方图 + 模型高斯曲线 =====
    figure('Position', [100 100 900 500]);
    hold on;

    valid_1D = isfinite(normed_delta_mix);
    nd_valid = normed_delta_mix(valid_1D);
    lab_valid_1D = label_1D_mix(valid_1D);

    nbins_mix = 60;
    edges_mix = linspace(min(nd_valid), max(nd_valid), nbins_mix+1);
    bw_mix = edges_mix(2) - edges_mix(1);

    % Unassigned
    histogram(nd_valid(lab_valid_1D==0), edges_mix, ...
        'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'k', ...
        'FaceAlpha', 0.5, 'DisplayName', 'Unassigned');

    % 各类
    for s = 1:n_types
        histogram(nd_valid(lab_valid_1D==s), edges_mix, ...
            'FaceColor', model1D(s).color, 'EdgeColor', 'k', ...
            'FaceAlpha', 0.6, 'DisplayName', model1D(s).label);
    end

    % 高斯拟合曲线 (来自训练模型)
    xfit_mix = linspace(edges_mix(1), edges_mix(end), 500);
    for s = 1:n_types
        N_s = sum(lab_valid_1D == s);
        yfit_s = N_s * bw_mix * normpdf(xfit_mix, model1D(s).mu, model1D(s).sigma);
        plot(xfit_mix, yfit_s, '-', 'Color', model1D(s).color*0.6, ...
            'LineWidth', 2.5, 'DisplayName', sprintf('Fit %s', model1D(s).label));
    end

    xlabel('\Delta brightness / initial brightness');
    ylabel('Counts');
    title(sprintf('Mixed Sample — 1D Classification [%s]', pol_label));
    legend('Location', 'best');
    box on; hold off;

%% ===== 可视化 3: 2D 散点图 + 决策边界 =====
    figure('Position', [100 50 1000 700]);
    hold on;

    valid_2D = isfinite(normed_delta_mix) & isfinite(initial_brightness_mix);
    nd_2d = normed_delta_mix(valid_2D);
    ib_2d = initial_brightness_mix(valid_2D);
    lab_2d = label_2D_mix(valid_2D);

    % Unassigned
    scatter(nd_2d(lab_2d==0), ib_2d(lab_2d==0), 15, [0.7 0.7 0.7], 'filled', ...
        'MarkerFaceAlpha', 0.3, 'DisplayName', 'Unassigned');

    % 各类
    for s = 1:n_types
        mask_s = (lab_2d == s);
        scatter(nd_2d(mask_s), ib_2d(mask_s), 30, model2D(s).color, 'filled', ...
            'MarkerFaceAlpha', 0.6, 'DisplayName', model2D(s).label);
    end

    % 2D 高斯椭圆 (1σ, 2σ)
    theta_ell = linspace(0, 2*pi, 200);
    for s = 1:n_types
        mu_s = model2D(s).mu;
        [V, D] = eig(model2D(s).Sigma);
        for n_sig = [1, 2]
            ell = n_sig * V * sqrt(D) * [cos(theta_ell); sin(theta_ell)];
            plot(mu_s(1) + ell(1,:), mu_s(2) + ell(2,:), '-', ...
                'Color', model2D(s).color*0.6, 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
        plot(mu_s(1), mu_s(2), 'x', 'Color', model2D(s).color*0.4, ...
            'MarkerSize', 14, 'LineWidth', 3, 'HandleVisibility', 'off');
    end

    % 决策边界
    n_grid = 200;
    xg = linspace(min(nd_2d), max(nd_2d), n_grid);
    yg = linspace(min(ib_2d), max(ib_2d), n_grid);
    [Xg, Yg] = meshgrid(xg, yg);
    Zg = zeros(n_grid, n_grid);
    for gi = 1:n_grid
        for gj = 1:n_grid
            xi_g = [Xg(gi,gj), Yg(gi,gj)];
            p_g = zeros(1, n_types);
            for s = 1:n_types
                p_g(s) = priors_2D(s) * mvnpdf(xi_g, model2D(s).mu, model2D(s).Sigma);
            end
            [~, Zg(gi,gj)] = max(p_g);
        end
    end
    contour(Xg, Yg, Zg, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    xlabel('\Delta brightness / initial brightness');
    ylabel('Initial Brightness (counts)');
    title(sprintf('Mixed Sample — 2D Classification + Decision Boundaries [%s]', pol_label));
    legend('Location', 'best');
    grid on; box on; hold off;

%% ===== 可视化 4: 不一致标注高亮 =====
    figure('Position', [150 150 700 600]);
    imagesc(mean_img); colormap(gray); axis image; hold on;

    % 一致的 (assigned)
    agree_assigned = agree_mix & (label_1D_mix > 0);
    scatter(pos_mix(agree_assigned,2), pos_mix(agree_assigned,1), 100, ...
        'g', 'o', 'LineWidth', 1, 'DisplayName', 'Agree (assigned)');

    % 不一致的
    scatter(pos_mix(disagree_mix,2), pos_mix(disagree_mix,1), 200, ...
        'm', 'x', 'LineWidth', 2.5, 'DisplayName', 'Disagree (1D≠2D)');

    % 均 unassigned
    both_unassigned = (label_1D_mix == 0) & (label_2D_mix == 0);
    scatter(pos_mix(both_unassigned,2), pos_mix(both_unassigned,1), 60, ...
        [0.5 0.5 0.5], '.', 'DisplayName', 'Both unassigned');

    legend('Location', 'best', 'FontSize', 11);
    title(sprintf('1D vs 2D Agreement Map | Agree=%.1f%% | Disagree=%.1f%%', ...
        100*mean(agree_mix), 100*mean(disagree_mix)));
    hold off;

%% ===== 可视化 5: 分类汇总统计图 =====
    figure('Position', [200 200 800 400]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % 饼图/条形图 - 1D
    nexttile;
    counts_1D = zeros(n_types+1, 1);
    counts_1D(1) = sum(label_1D_mix == 0);
    for s = 1:n_types
        counts_1D(s+1) = sum(label_1D_mix == s);
    end
    bar_colors = [0.7 0.7 0.7];
    for s = 1:n_types
        bar_colors = [bar_colors; model1D(s).color]; %#ok<AGROW>
    end
    b1 = bar(counts_1D, 'FaceColor', 'flat');
    b1.CData = bar_colors;
    set(gca, 'XTickLabel', [{'Unassigned'}, {model1D.label}]);
    ylabel('Count');
    title('1D Classification Summary');
    for k = 1:(n_types+1)
        text(k, counts_1D(k)+0.5, sprintf('%d\n(%.1f%%)', counts_1D(k), 100*counts_1D(k)/N_mix), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    box on;

    % 饼图/条形图 - 2D
    nexttile;
    counts_2D = zeros(n_types+1, 1);
    counts_2D(1) = sum(label_2D_mix == 0);
    for s = 1:n_types
        counts_2D(s+1) = sum(label_2D_mix == s);
    end
    b2 = bar(counts_2D, 'FaceColor', 'flat');
    b2.CData = bar_colors;
    set(gca, 'XTickLabel', [{'Unassigned'}, {model2D.label}]);
    ylabel('Count');
    title('2D Classification Summary');
    for k = 1:(n_types+1)
        text(k, counts_2D(k)+0.5, sprintf('%d\n(%.1f%%)', counts_2D(k), 100*counts_2D(k)/N_mix), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    box on;

    sgtitle(sprintf('Mixed Sample Classification Summary (N=%d) [%s]', N_mix, pol_label), 'FontSize', 13);

    fprintf('\n===== Mixed Sample Classification Complete =====\n');


    

%% ===== 保存模型 =====
    classificationModel.model1D = model1D;
    classificationModel.model2D = model2D;
    classificationModel.params.ratio_th_1D = ratio_th_1D;
    classificationModel.params.ratio_th_2D = ratio_th_2D;
    classificationModel.params.priors_1D = priors_1D;
    classificationModel.params.priors_2D = priors_2D;
    classificationModel.params.polarity = polarity;
    classificationModel.params.size_idx = size_idx;
    classificationModel.eval.accuracy_1D = accuracy_1D;
    classificationModel.eval.accuracy_2D = accuracy_2D;
    classificationModel.eval.unassigned_1D = unassigned_rate_1D;
    classificationModel.eval.unassigned_2D = unassigned_rate_2D;
    classificationModel.eval.conf_mat_1D = conf_mat_1D;
    classificationModel.eval.conf_mat_2D = conf_mat_2D;
    classificationModel.eval.LL_1D = LL_1D;
    classificationModel.eval.LL_2D = LL_2D;
    classificationModel.eval.BIC_1D = BIC_1D;
    classificationModel.eval.BIC_2D = BIC_2D;
    classificationModel.eval.precision_1D = precision_1D;
    classificationModel.eval.recall_1D = recall_1D;
    classificationModel.eval.f1_1D = f1_1D;
    classificationModel.eval.precision_2D = precision_2D;
    classificationModel.eval.recall_2D = recall_2D;
    classificationModel.eval.f1_2D = f1_2D;

    save('classificationModel.mat', 'classificationModel');
    fprintf('\n模型已保存至 classificationModel.mat\n');