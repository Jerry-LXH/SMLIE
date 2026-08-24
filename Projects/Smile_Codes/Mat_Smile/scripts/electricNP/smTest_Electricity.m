
%% Parameter settings
    clear all;clc;

    % file_name = '/Volumes/SMILeSSD/UCNP charge data/20251028 98Yb2Er@Lu-2.5nm/6 Y@98Yb2Er@Lu-2.5nm-ZJ first no E 1V and -1V every5s 5cycles then no E.sif';
    % file_name = '/Volumes/SMILeSSD/UCNP charge data/20251028 98Yb2Er@Lu-6nm/9 Y@98Yb2Er@Lu-6nm-ZJ first no E -1V and 1V every5s 5cycles then no E.sif';
    % file_name = '/Volumes/SMILeSSD/UCNP charge data/20251222+1223 small-98Yb2Er no NaYF4 300mA-laser/1223/6 s98Yb2Er-ZJ first no E 1V and -1V every5s 5cycles then no E.sif';
    % file_name = '/Volumes/SMILeSSD/UCNP charge data/20251029 98Yb2Er@Lu-10nm/3 Y@98Yb2Er@Lu-10nm-ZJ first no E 1V and -1V every5s 5cycles then no E.sif';


    % file_name = '/Volumes/SMILeSSD/UCNP charge data/20260317 98Yb2Er mix 2.5nm 6nm (shell) 300mA-laser/4 98Yb2Er-ZJ first no E -1V and 1V every5s 5cycles then no E.sif';
    % file_name = '/Volumes/SMILeSSD/UCNP charge data/20260401 Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and 98Yb2Er@Lu(without Y) 300mA-laser/3mix Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu)/1 98Yb2Er-ZJ 3mix first no E -1V and 1V every5s 5cycles then no E.sif';
    file_name = '/Volumes/SMILeSSD/UCNP charge data/20260401 Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and 98Yb2Er@Lu(without Y) 300mA-laser/4mix Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and small98Yb2Er@Lu(without Y)/1 98Yb2Er-ZJ 4mix first no E -1V and 1V every5s 5cycles then no E.sif';

%% Load Pipeline
    [raw_data,ex_time,gainDAC] = io.readSIFData(file_name);
    interval = 0; % time interval between two adjacent exposure. In unit of seconds.
    oneFrameTime = ex_time+interval;
    windowed_raw_data = raw_data(100:412,100:412,:); % Should be a bit larger then 257*257 (128:384) for drift correction. Should be CENTERED.
    frames = size(windowed_raw_data,3);
    clear raw_data;

%% Detect and Loc Pipeline
    clear loc_total;

    % detect local max, k_sigma is the threshold control parameter
    k_sigma = 2;
    uncorrected_detected_total = detect.findMaxima(windowed_raw_data,k_sigma); 

    % localize using MLE. LSQ is also avaliable. The output has 7 columns: row_loc, col_loc, Brightness N, background b, width sigma, uncertainty sigma_loc and frame number f. 
    uncorrected_super_loc_total_raw = localize.MLE1.locMolecules(windowed_raw_data,uncorrected_detected_total,5);
    uncorrected_super_loc_total = localize.filterBadLocs(uncorrected_super_loc_total_raw);
    uncorrected_loc_total = uncorrected_super_loc_total(:,[1,2,7]);

    % save localization results
    [file_dir, file_base, ~] = fileparts(file_name);
    save_name = fullfile(file_dir, ...
        [file_base '_uncorrected_super_loc_total.mat']);
    save(save_name, 'uncorrected_super_loc_total', '-v7.3');
    fprintf('Detection saved to:\n%s\n', save_name);

%% Drift and Emitters Pipeline 。
    % read saved localizations (N*7 matrix)
    [file_dir, file_base, ~] = fileparts(file_name);
    load_name = fullfile(file_dir, ...
        [file_base '_uncorrected_super_loc_total.mat']);
    S = load(load_name);
    uncorrected_super_loc_total = S.uncorrected_super_loc_total;
    uncorrected_loc_total = uncorrected_super_loc_total(:,[1,2,7]);
    fprintf('Detection loaded from:\n%s\n', load_name);

    % drift correction, note that the input has to be N*3 matrix
    drift_corr = true;
    uncertainty = mean(uncorrected_super_loc_total(:,6));
    if drift_corr == true
        delta_sum = postproc.drift.estimateDrift_cloudxy(windowed_raw_data, uncorrected_loc_total, 15, uncertainty,30); % Note that loc_total have only 3 columns of (r,c,f).
    else
        delta_sum = zeros(frames, 2);
    end

    % correct image and localizations by fitted drift. Note that loc_idx represents the row-index of drift-corrected localizations in drift-uncorrected localization matrix. Note that the last input of correctDrift() controls the size of drift-corrected image (should be odd number, noramlly set to 257px).
    [data, loc_total, loc_idx] = postproc.drift.correctDrift(windowed_raw_data,uncorrected_loc_total,delta_sum,5,257);
    [corrected_super_loc_total] = uncorrected_super_loc_total(loc_idx,:); % this contains all 7 columns
    viz.plotTracking(delta_sum, oneFrameTime, [], 'Drift Tracking'); % visualize drift

    % extract other localization information of drift-corrected locs
    brightness = corrected_super_loc_total(:,3);   % photon
    background = corrected_super_loc_total(:,4);   % background
    sigma      = corrected_super_loc_total(:,5);   % PSF width
    sigma_loc  = corrected_super_loc_total(:,6);   % localization 

    % Visualize the corrected subpixel locs
    %num_show = min(7000, frames);
    %frames_to_show = 1:num_show;
    %figure;
    %viz.plotImage(data, frames_to_show, 'hot','Drift-corrected Image and Locs','max');
    %hold on;
    %viz.overlayLocs(loc_total,frames_to_show,true);

    % Emittters analysis: put nearby locs into one emitter. bleach_time sets the max waitting time beyond which the emitter will be set as bleached if not detected.
    bleach_time = 10; %seconds
    bleach_frames = bleach_time/(oneFrameTime); 
    searching_radius = 2; %pixel
    emitters = postproc.emitter.findEmitters(data, loc_total, bleach_frames, searching_radius);
    emitters = postproc.emitter.mergeEmitters(emitters, searching_radius);

    % Filter short-lived/not-in-first-frame/unbleached-until-last-frame/large-jumping(potentially 2 overlapping emitters) emitters
    livetime_th = 0.2; %seconds
    emitters_filt = postproc.emitter.filterEmitters_short(emitters, round(livetime_th/ex_time),'consecutive');
    emitters_filt =  postproc.emitter.filterEmitters_firstframe(emitters_filt);
    %emitters_filt = postproc.emitter.filterEmitters_end(emitters_filt); 
    [emitters_filt,stats_jump] = postproc.emitter.filterEmitters_jumping(emitters_filt,3.5);
    %emitters_filt = emitters;
    % postproc.emitter.plotJumpStats(stats_jump);

    % Check all emitters
    %figure;
    %viz.plotImage(data, 1:700, 'gray','Drift-corrected Image and Emitters','max');
    %hold on;
    %for k = 1:numel(emitters_filt)
    %    plot(emitters_filt(k).col+0.5, emitters_filt(k).row+0.5, '-', 'LineWidth', 1);
    %end
    %title('Emitter trajectories overlaid on image');

    % Important Data Collection (position, mean position, survival time)
    stats = postproc.emitter.collectEmitterStatistics(emitters_filt,frames,ex_time,interval,brightness,sigma,sigma_loc,background);
    % postproc.emitter.plotEmitterStatistics(stats);


    % Extract trace (background calculated by surrounding pixels)
    int_range = 4;
    bg_range = 7;
    raw_series = analysis.photophys.extractTrace(stats.pos_mean_px, data(:,:,:),int_range,true,bg_range,'median'); 
    stats.trace = raw_series;

%% Check single emitter and its trace
    start_frame = 1;
    end_frame = frames;
    for index = 21:60 %61
        postproc.emitter.checkTrace(start_frame,end_frame,ex_time,stats.trace(index,:),stats.brightness_em(index,:),stats.pos_matrix(:,1,index),stats.pos_matrix(:,2,index));
    end

%% PhotoAnalysis of Trace
    % emitters_filt = analysis.photophys.analyzeStates(emitters_filt, 'HMM', 'bleachTail',50,'bgThreshold',100,'bicPenalty', 0.1, 'minStateSep', 1.0);
    states = analysis.photophys.analyzeStates(stats.trace, 'CHANGEPOINT','penalty',2.5,'minSegLen',20,'mergeThr', 1.2, 'bleachTail',50,'bgThreshold',60, 'minStateSep', 1.3);

    % Delta of max and min (note that the second input sets the minium occuring times, excluding badly-fitted traces)
    delta_brightness = analysis.photophys.computeDelta(states, 4);

    % initial brightness from first 50 frames
    nInit = min(50, size(stats.trace, 2));
    %initial_brightness = mean(stats.trace(:, 1:nInit), 2, 'omitnan'); % calculated by bg-extracted trace
    initial_brightness = mean(stats.brightness_em(:, 1:nInit), 2, 'omitnan'); % calculated by MLE result (should be indentical/very close to trace result)

    % normalizd delta_brightness
    normed_delta_brightness = delta_brightness./initial_brightness;

    % Viz
    delta_valid = normed_delta_brightness(~isnan(normed_delta_brightness));
    figure;
    histogram(delta_valid, 90, ...
        'FaceColor', [0.2 0.6 0.8], ...
        'EdgeColor', 'k');
    xlabel('\Delta brightness / initial brightness');
    ylabel('Counts');
    title('Distribution of normalized \Delta brightness');
    box on;

%% Visualization of Fitted States
    % analysis.photophys.plotStates(emitters_filt,1:15);
    for index = 1:30
        analysis.photophys.plotState(states(index));
    end

%% Check emitter video 
    % Note: The first input is index.
    postproc.emitter.checkMovie(33, data, emitters_filt, stats.pos_mean_px, stats.pos_matrix, oneFrameTime, 'frameStep',2,'pauseTime', 0.005,'clim', [0 250],'preFrames', 100);

%% Save Single Particle Result
    [file_dir, file_base, ~] = fileparts(file_name);
    save_name = fullfile(file_dir, ...
        [file_base '_normed_delta_brightness.mat']);
    save(save_name, 'normed_delta_brightness', '-v7.3');

    % Save Single Particle Result (brightness)
    [file_dir, file_base, ~] = fileparts(file_name);
    save_name = fullfile(file_dir, ...
        [file_base '_initial_brightness.mat']);
    save(save_name, 'initial_brightness', '-v7.3');

%% Read mat files of saved brightness/delta data
    dataDir_2p5 = '/Volumes/SMILeSSD/UCNP charge data/2.5nm';   
    dataDir_6   = '/Volumes/SMILeSSD/UCNP charge data/6nm';   
    dataDir_10 = '/Volumes/SMILeSSD/UCNP charge data/10nm';
    dataDir_small = '/Volumes/SMILeSSD/UCNP charge data/small';

    pattern = '*-1V and 1V every5s 5cycles then no E_normed_delta_brightness.mat';

    normed_delta_brightness_all_2p5nm = utils.read_mat_1d_vectors(dataDir_2p5, pattern);
    normed_delta_brightness_all_6nm   = utils.read_mat_1d_vectors(dataDir_6, pattern);
    normed_delta_brightness_all_10nm   = utils.read_mat_1d_vectors(dataDir_10, pattern);
    normed_delta_brightness_all_small   = utils.read_mat_1d_vectors(dataDir_small, pattern);

%% 四种粒子 normalized delta brightness 分布对比（单图）
    figure('Position', [100 100 900 500]);
    hold on;

    % 准备数据（去除 NaN）
    d_2p5  = normed_delta_brightness_all_2p5nm(isfinite(normed_delta_brightness_all_2p5nm));
    d_6    = normed_delta_brightness_all_6nm(isfinite(normed_delta_brightness_all_6nm));
    d_10   = normed_delta_brightness_all_10nm(isfinite(normed_delta_brightness_all_10nm));
    d_sm   = normed_delta_brightness_all_small(isfinite(normed_delta_brightness_all_small));

    % 统一 bin edges
    all_data = [d_2p5(:); d_6(:); d_10(:); d_sm(:)];
    nbins = 50;
    edges = linspace(min(all_data), max(all_data), nbins + 1);

    % 颜色定义
    c1 = [0.2  0.6  0.8];   % 2.5 nm
    c2 = [0.9  0.4  0.4];   % 6 nm
    c3 = [0.4  0.8  0.4];   % 10 nm
    c4 = [0.7  0.5  0.9];   % small

    % ---- 方法 A：用 histogram 半透明叠加 ----
    h1 = histogram(d_2p5, edges, 'FaceColor', c1, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.35, 'Normalization', 'probability', ...
        'DisplayName', sprintf('2.5 nm  (N=%d, \\mu=%.3f, \\sigma=%.3f)', numel(d_2p5), mean(d_2p5), std(d_2p5)));
    h2 = histogram(d_6,   edges, 'FaceColor', c2, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.35, 'Normalization', 'probability', ...
        'DisplayName', sprintf('6 nm    (N=%d, \\mu=%.3f, \\sigma=%.3f)', numel(d_6), mean(d_6), std(d_6)));
    h3 = histogram(d_10,  edges, 'FaceColor', c3, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.35, 'Normalization', 'probability', ...
        'DisplayName', sprintf('10 nm   (N=%d, \\mu=%.3f, \\sigma=%.3f)', numel(d_10), mean(d_10), std(d_10)));
    h4 = histogram(d_sm,  edges, 'FaceColor', c4, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.35, 'Normalization', 'probability', ...
        'DisplayName', sprintf('small   (N=%d, \\mu=%.3f, \\sigma=%.3f)', numel(d_sm), mean(d_sm), std(d_sm)));

    % ---- 方法 B：叠加核密度估计曲线 (KDE) ----
    xfit = linspace(min(all_data), max(all_data), 500);

    if numel(d_2p5) > 2
        [f1, xi1] = ksdensity(d_2p5, xfit);
        plot(xi1, f1 / max(f1) * max(h1.Values), '-', 'Color', c1*0.7, 'LineWidth', 2.2, 'HandleVisibility', 'off');
    end
    if numel(d_6) > 2
        [f2, xi2] = ksdensity(d_6, xfit);
        plot(xi2, f2 / max(f2) * max(h2.Values), '-', 'Color', c2*0.7, 'LineWidth', 2.2, 'HandleVisibility', 'off');
    end
    if numel(d_10) > 2
        [f3, xi3] = ksdensity(d_10, xfit);
        plot(xi3, f3 / max(f3) * max(h3.Values), '-', 'Color', c3*0.7, 'LineWidth', 2.2, 'HandleVisibility', 'off');
    end
    if numel(d_sm) > 2
        [f4, xi4] = ksdensity(d_sm, xfit);
        plot(xi4, f4 / max(f4) * max(h4.Values), '-', 'Color', c4*0.7, 'LineWidth', 2.2, 'HandleVisibility', 'off');
    end

    xlabel('\Delta brightness / initial brightness');
    % xlabel('Initial brightness(pps)');
    ylabel('Probability');
    title('Comparison of normalized \Delta brightness across particle types');
    % title('Comparison of Initial brightness across particle types');
    legend('Location', 'best', 'FontSize', 9);
    box on;
    hold off;

%% ---- 补充统计量打印输出 ----
    fprintf('\n===== Distribution Summary =====\n');
    types  = {'2.5 nm', '6 nm', '10 nm', 'small'};
    dsets  = {d_2p5, d_6, d_10, d_sm};
    for k = 1:4
        d = dsets{k};
        fprintf('%s : N=%4d  mean=%.4f  std=%.4f  median=%.4f  IQR=[%.4f, %.4f]\n', ...
            types{k}, numel(d), mean(d), std(d), median(d), ...
            prctile(d,25), prctile(d,75));
    end

    % ---- 两两 Cohen's d ----
    fprintf('\n===== Pairwise Cohen''s d =====\n');
    for i = 1:4
        for j = i+1:4
            di = dsets{i}; dj = dsets{j};
            pooled_std = sqrt(((numel(di)-1)*var(di) + (numel(dj)-1)*var(dj)) / ...
                            (numel(di)+numel(dj)-2));
            cohen_d = abs(mean(di) - mean(dj)) / max(pooled_std, eps);
            fprintf('  %s vs %s : d = %.3f\n', types{i}, types{j}, cohen_d);
        end
    end

    % ---- 两两分布重叠面积 (Overlap Coefficient via KDE) ----
    fprintf('\n===== Pairwise Overlap Coefficient =====\n');
    xcommon = linspace(min(all_data), max(all_data), 2000);
    kdes = cell(4,1);
    for k = 1:4
        if numel(dsets{k}) > 2
            kdes{k} = ksdensity(dsets{k}, xcommon);
        else
            kdes{k} = zeros(size(xcommon));
        end
    end
    for i = 1:4
        for j = i+1:4
            ovl = trapz(xcommon, min(kdes{i}, kdes{j}));
            fprintf('  %s vs %s : OVL = %.3f  (1=complete overlap)\n', types{i}, types{j}, ovl);
        end
    end

%% Build Model (from single-kind beads, 2 kind)
    type1 = normed_delta_brightness_all_2p5nm(:);
    type2 = normed_delta_brightness_all_6nm(:);
    mix   = normed_delta_brightness(:);
    % remove NaN
    type1 = type1(isfinite(type1));
    type2 = type2(isfinite(type2));
    N = numel(mix);
    label = zeros(N,1);   % 0 = unassigned, 1 = type1, 2 = type2
    % fit Gaussian model
    mu1 = mean(type1);
    sg1 = std(type1);
    mu2 = mean(type2);
    sg2 = std(type2);
    % avoid zero std
    sg1 = max(sg1, eps);
    sg2 = max(sg2, eps);
    % prior probabilities, can also use 0.5 / 0.5
    prior1 = 0.5;
    prior2 = 0.5;
    for i = 1:N
        x = mix(i);
        if ~isfinite(x)
            label(i) = 0;
            continue;
        end
        p1 = prior1 * normpdf(x, mu1, sg1);
        p2 = prior2 * normpdf(x, mu2, sg2);
        % confidence criterion
        ratio = max(p1, p2) / max(min(p1, p2), realmin);
        if ratio < 3
            label(i) = 0;   % ambiguous
        elseif p1 > p2
            label(i) = 1;
        else
            label(i) = 2;
        end
    end
    % Viz
    figure;
    viz.plotImage(data, 1:700, 'gray', 'Classification');
    hold on;
    % positions
    pos = stats.pos_mean_px;

    idx_2p5 = (label == 1);
    idx_6   = (label == 2);

    h1 = scatter(pos(idx_2p5,2), pos(idx_2p5,1), 150, 'b', 'o', 'LineWidth', 1);
    h2 = scatter(pos(idx_6,2), pos(idx_6,1), 150, 'o', 'o', 'LineWidth', 1);

    legend([h1, h2], {'2.5 nm', '6 nm'}, 'Location', 'best');
    title('2-type classification');

    % Viz 2
    delta_all = normed_delta_brightness(:);

    % keep only finite values and matched labels
    valid = isfinite(delta_all) & isfinite(label(:));

    delta_all = delta_all(valid);
    label_valid = label(valid);

    delta_0 = delta_all(label_valid == 0);
    delta_1 = delta_all(label_valid == 1);
    delta_2 = delta_all(label_valid == 2);

    % ===== unified bin edges =====
    nbins = 60;
    edges = linspace(min(delta_all), max(delta_all), nbins+1);
    binWidth = edges(2) - edges(1);

    xfit = linspace(min(delta_all), max(delta_all), 500);

    figure;
    hold on;

    % histograms
    histogram(delta_0, edges, ...
        'FaceColor', [0.85 0.85 0.85], ...
        'EdgeColor', 'k', ...
        'DisplayName', 'Unassigned');

    histogram(delta_1, edges, ...
        'FaceColor', [0.2 0.6 0.8], ...
        'EdgeColor', 'k', ...
        'DisplayName', '2.5 nm');

    histogram(delta_2, edges, ...
        'FaceColor', [0.9 0.4 0.4], ...
        'EdgeColor', 'k', ...
        'DisplayName', '6 nm');

    % ===== fitted Gaussian curves =====
    % scale PDF to histogram counts
    y1 = numel(delta_1) * binWidth * normpdf(xfit, mu1, sg1);
    y2 = numel(delta_2) * binWidth * normpdf(xfit, mu2, sg2);

    plot(xfit, y1, '-', ...
        'Color', [0.0 0.3 0.8], ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Fit 2.5 nm');

    plot(xfit, y2, '-', ...
        'Color', [0.8 0.1 0.1], ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Fit 6 nm');

    xlabel('\Delta brightness / initial brightness');
    ylabel('Counts');
    title('Distribution of normalized \Delta brightness with Gaussian fits');
    legend('Location', 'best');
    box on;
    hold off;

%% Build Model (from single-kind beads, 3 kind)
    type1 = normed_delta_brightness_all_2p5nm(:);
    type2 = normed_delta_brightness_all_6nm(:);
    type3 = normed_delta_brightness_all_10nm(:);
    mix   = normed_delta_brightness(:);

    % remove NaN / Inf
    type1 = type1(isfinite(type1));
    type2 = type2(isfinite(type2));
    type3 = type3(isfinite(type3));

    N = numel(mix);
    label = zeros(N,1);   % 0 = unassigned, 1 = type1, 2 = type2, 3 = type3

    % fit Gaussian model
    mu1 = mean(type1);   sg1 = std(type1);
    mu2 = mean(type2);   sg2 = std(type2);
    mu3 = mean(type3);   sg3 = std(type3);

    % avoid zero std
    sg1 = max(sg1, eps);
    sg2 = max(sg2, eps);
    sg3 = max(sg3, eps);

    % prior probabilities
    prior1 = 1/3;
    prior2 = 1/3;
    prior3 = 1/3;

    % confidence threshold
    ratio_th = 3;

    for i = 1:N
        x = mix(i);

        if ~isfinite(x)
            label(i) = 0;
            continue;
        end

        p1 = prior1 * normpdf(x, mu1, sg1);
        p2 = prior2 * normpdf(x, mu2, sg2);
        p3 = prior3 * normpdf(x, mu3, sg3);

        p_all = [p1, p2, p3];

        % sort probabilities
        [p_sorted, idx_sorted] = sort(p_all, 'descend');

        % compare best vs second best
        ratio = p_sorted(1) / max(p_sorted(2), realmin);

        if ratio < ratio_th
            label(i) = 0;   % ambiguous
        else
            label(i) = idx_sorted(1);
        end
    end

    % Viz 1
    figure;
    viz.plotImage(data, 1:700, 'gray', 'Classification');
    hold on;

    % positions
    pos = stats.pos_mean_px;

    idx_1 = (label == 1);
    idx_2 = (label == 2);
    idx_3 = (label == 3);

    h1 = scatter(pos(idx_1,2), pos(idx_1,1), 150, 'b', 'o', 'LineWidth', 1);
    h2 = scatter(pos(idx_2,2), pos(idx_2,1), 150, 'r', 'o', 'LineWidth', 1);
    h3 = scatter(pos(idx_3,2), pos(idx_3,1), 150, 'g', 'o', 'LineWidth', 1);

    legend([h1, h2, h3], {'2.5 nm', '6 nm', '10 nm'}, 'Location', 'best');
    title('3-type classification');

    % Viz 2
    delta_all = normed_delta_brightness(:);

    % keep only finite values and matched labels
    valid = isfinite(delta_all) & isfinite(label(:));

    delta_all = delta_all(valid);
    label_valid = label(valid);

    delta_0 = delta_all(label_valid == 0);
    delta_1 = delta_all(label_valid == 1);
    delta_2 = delta_all(label_valid == 2);
    delta_3 = delta_all(label_valid == 3);

    % ===== unified bin edges =====
    nbins = 60;
    edges = linspace(min(delta_all), max(delta_all), nbins+1);
    binWidth = edges(2) - edges(1);

    figure;
    hold on;

    % histograms
    h0 = histogram(delta_0, edges, ...
        'FaceColor', [0.85 0.85 0.85], ...
        'EdgeColor', 'k', ...
        'DisplayName', 'Unassigned');

    h1 = histogram(delta_1, edges, ...
        'FaceColor', [0.2 0.6 0.8], ...
        'EdgeColor', 'k', ...
        'DisplayName', '2.5nm');

    h2 = histogram(delta_2, edges, ...
        'FaceColor', [0.9 0.4 0.4], ...
        'EdgeColor', 'k', ...
        'DisplayName', '6nm');

    h3 = histogram(delta_3, edges, ...
        'FaceColor', [0.4 0.8 0.4], ...
        'EdgeColor', 'k', ...
        'DisplayName', '10nm');

    % ===== fitted Gaussian curves =====
    xfit = linspace(min(delta_all), max(delta_all), 500);

    % scale PDF to histogram counts
    y1 = numel(delta_1) * binWidth * normpdf(xfit, mu1, sg1);
    y2 = numel(delta_2) * binWidth * normpdf(xfit, mu2, sg2);
    y3 = numel(delta_3) * binWidth * normpdf(xfit, mu3, sg3);

    p1 = plot(xfit, y1, '-', 'Color', [0.0 0.3 0.8], 'LineWidth', 2.5, ...
        'DisplayName', 'Fit 2.5nm');
    p2 = plot(xfit, y2, '-', 'Color', [0.8 0.1 0.1], 'LineWidth', 2.5, ...
        'DisplayName', 'Fit 6nm');
    p3 = plot(xfit, y3, '-', 'Color', [0.1 0.6 0.1], 'LineWidth', 2.5, ...
        'DisplayName', 'Fit 10nm');

    xlabel('\Delta brightness / initial brightness');
    ylabel('Counts');
    title('Distribution of normalized \Delta brightness with Gaussian fits');
    legend('Location', 'best');
    box on;
    hold off;

%% Combined visualization: marginal + joint distribution

    normed_delta_brightness_all = normed_delta_brightness;
    initial_brightness_all = initial_brightness;

    x = initial_brightness_all / ex_time;
    y = normed_delta_brightness_all;

    valid = ~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y);
    x_valid = x(valid);
    y_valid = y(valid);

    figure;
    tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Top-left: initial brightness histogram
    nexttile(1);
    histogram(x_valid, 30, ...
        'FaceColor', [0.2 0.6 0.8], ...
        'EdgeColor', 'k');
    xlabel('Initial Brightness (Photons/sec)');
    ylabel('Counts');
    title('Initial Brightness');
    box on;

    % Top-right: empty
    nexttile(2);
    axis off;

    % Bottom-left: joint scatter
    nexttile(3);
    scatter(x_valid, y_valid, 15, ...
        'MarkerFaceColor', [0.85 0.33 0.10], ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.4);
    xlabel('Initial Brightness (Photons/sec)');
    ylabel('\Delta brightness / initial brightness');
    title('Joint Distribution');
    grid on;
    box on;

    % Bottom-right: normalized delta histogram
    nexttile(4);
    histogram(y_valid, 30, ...
        'FaceColor', [0.4 0.7 0.4], ...
        'EdgeColor', 'k');
    xlabel('\Delta brightness / initial brightness');
    ylabel('Counts');
    title('Normalized \Delta Brightness');
    box on;

%% Combined visualization: marginal + joint distribution (color-coded by model)

    x = initial_brightness(:) / ex_time;
    y = normed_delta_brightness(:);

    valid = isfinite(x) & isfinite(y) & isfinite(label);
    x_v = x(valid);
    y_v = y(valid);
    lab = label(valid);

    % color & name
    clr = [0.75 0.75 0.75;    % 0 = unassigned
        0.2  0.6  0.8;     % 1 = 2.5 nm
        0.9  0.4  0.4;     % 2 = 6 nm
        0.4  0.8  0.4];    % 3 = 10 nm

    clr_line = [0.5 0.5 0.5;
                0.0 0.3 0.8;
                0.8 0.1 0.1;
                0.1 0.6 0.1];

    names = {'Unassigned', '2.5 nm', '6 nm', '10 nm'};

    idx = cell(4,1);
    for k = 0:3
        idx{k+1} = (lab == k);
    end

    figure('Position', [100 100 1000 850]);
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % ===== Top-left: Initial brightness marginal =====
    nexttile(1);
    hold on;
    edges_x = linspace(min(x_v), max(x_v), 40);
    for k = [1 2 3 4]   % plot unassigned first (behind)
        histogram(x_v(idx{k}), edges_x, ...
            'FaceColor', clr(k,:), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.6, 'DisplayName', names{k});
    end
    xlabel('Initial Brightness (Photons/sec)');
    ylabel('Counts');
    title('Initial Brightness (marginal)');
    legend('Location', 'best');
    box on; hold off;

    % ===== Top-right: Classification summary =====
    nexttile(2);
    axis off;
    str = {sprintf('\\bf Classification Summary'), ''};
    for k = 0:3
        n_k = sum(lab == k);
        str{end+1} = sprintf('  %s :  N = %d  (%.1f%%)', ...
            names{k+1}, n_k, 100*n_k/numel(lab));
    end
    str{end+1} = '';
    str{end+1} = sprintf('  Total :  N = %d', numel(lab));
    str{end+1} = '';
    str{end+1} = sprintf('  ratio\\_th = %.1f', ratio_th);
    str{end+1} = sprintf('  prior = [1/3, 1/3, 1/3]');
    text(0.05, 0.95, str, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 11, ...
        'FontName', 'FixedWidth', 'Interpreter', 'tex');

    % ===== Bottom-left: Joint scatter =====
    nexttile(3);
    hold on;
    % unassigned behind
    scatter(x_v(idx{1}), y_v(idx{1}), 12, clr(1,:), 'filled', ...
        'MarkerFaceAlpha', 0.25, 'DisplayName', names{1});
    for k = 2:4
        scatter(x_v(idx{k}), y_v(idx{k}), 25, clr(k,:), 'filled', ...
            'MarkerFaceAlpha', 0.55, 'DisplayName', names{k});
    end

    % add horizontal decision boundaries (Gaussian intersection points)
    pairs = [1 2; 2 3];   % mu1>mu2>mu3
    mus = [mu1 mu2 mu3];
    sgs = [sg1 sg2 sg3];
    xb = linspace(min(y_v), max(y_v), 1e4);
    for pp = 1:size(pairs,1)
        a = pairs(pp,1); b = pairs(pp,2);
        lp_a = log(1/3) - 0.5*((xb - mus(a))/sgs(a)).^2 - log(sgs(a));
        lp_b = log(1/3) - 0.5*((xb - mus(b))/sgs(b)).^2 - log(sgs(b));
        [~, ci] = min(abs(lp_a - lp_b));
        yline(xb(ci), '--', sprintf('%.3f', xb(ci)), ...
            'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, ...
            'LabelHorizontalAlignment', 'left', ...
            'HandleVisibility', 'off');
    end

    xlabel('Initial Brightness (Photons/sec)');
    ylabel('\Delta brightness / initial brightness');
    title('Joint Distribution (colored by model)');
    legend('Location', 'best');
    grid on; box on; hold off;

    % ===== Bottom-right: Normalized delta brightness marginal + Gaussian fits =====
    nexttile(4);
    hold on;
    edges_y = linspace(min(y_v), max(y_v), 50);
    bw = edges_y(2) - edges_y(1);

    for k = [1 2 3 4]
        histogram(y_v(idx{k}), edges_y, ...
            'FaceColor', clr(k,:), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.6, 'DisplayName', names{k});
    end

    % Gaussian model curves (scaled to assigned counts)
    xfit = linspace(min(y_v), max(y_v), 500);
    mus_all = [mu1 mu2 mu3];
    sgs_all = [sg1 sg2 sg3];
    for k = 1:3
        n_k = sum(idx{k+1});
        yfit = n_k * bw * normpdf(xfit, mus_all(k), sgs_all(k));
        plot(xfit, yfit, '-', 'Color', clr_line(k+1,:), 'LineWidth', 2.5, ...
            'DisplayName', sprintf('Gauss %s', names{k+1}));
    end

    xlabel('\Delta brightness / initial brightness');
    ylabel('Counts');
    title('Normalized \Delta Brightness (marginal)');
    legend('Location', 'best');
    box on; hold off;

    sgtitle('Mixed-bead classification: model-based coloring', 'FontSize', 14);