%% Parameter settings
    clear all;clc;

    file_name = '/Volumes/SMILeSSD/Optics/Ti2Tests/TirfModule/Cy3_polylys/20260616/poly_cy3_0.1s_200fm_1.66A_var/20260616-115551604/TUC-001.tif';
    %file_name = '/Volumes/SMILeSSD/Optics/TE2000UTests/Cy3_polylys/20260616/532nm_1x_0p1S_100frames_0p33mW_poly_cy3_varill_fov1.sif';
    %file_name_filt = '/Volumes/SMILeSSD/Optics/Ti2Tests/TirfModule/Cy3_polylys/20260227_532_cy3_0.2s_1.82A_filter/TUC-002.tif';
    gain = 33.448; % DN/e-
    ex_time = 0.1;
    offset = 1893.68;
    interval = 0; % time interval between two adjacent exposure
    oneFrameTime = ex_time+interval;
    
%% read folder into tiff
    folderPath = '/Volumes/SMILeSSD/Optics/Ti2Tests/TirfModule/Cy3_polylys/20260206_532_cy3Test_0.2s_1.72A/2_filt';
    files = dir(fullfile(folderPath, '*.tif'));
    [~,idx] = sort({files.name});
    files = files(idx);
    % 读第一张确定尺寸
    img = imread(fullfile(folderPath,files(1).name));
    [ny,nx] = size(img);

    n = numel(files);
    raw_data = zeros(ny,nx,n,class(img));

    for i = 1:n
        raw_data(:,:,i) = imread(fullfile(folderPath,files(i).name));
    end

%% write tiff
    outFile = [file_name(1:end-4) '.tif'];
    io.writeTiffStack(raw_data, outFile);

%% read tif stack
    raw_data = io.readTiffStack(file_name);
    % raw_data_filt = io.readTiffStack(file_name_filt);

%% Convert to electrons unit
    raw_data = double(raw_data) - offset;
    raw_data = raw_data / gain;
    %raw_data_filt = double(raw_data_filt) - offset;
    %raw_data_filt = raw_data_filt / gain;
    windowed_raw_data = raw_data;

%% Load sif data (Note that data is in absolute e- unit)
    [raw_data,ex_time,gainDAC] = io.readSIFData(file_name);
    interval = 0; % time interval between two adjacent exposure. In unit of seconds.
    oneFrameTime = ex_time+interval;
    windowed_raw_data = raw_data;
    
%% Hot Pixel Detection for Non-bleaching Particles
    %  Key idea: hot pixel is single-pixel sharp; real particle is PSF-broadened
    %  Insert after Windowing, before MLE / trace extraction

    % --- Parameters ---
    hot_pixel_params.spatial_kernel = 5;       % neighborhood size for local comparison
    hot_pixel_params.sharpness_thresh = 0.2;   % fraction of energy in center pixel vs ROI
    hot_pixel_params.persist_frac = 0.95;      % fraction of frames pixel must be "sharp" to be hot
    hot_pixel_params.intensity_percentile = 90; % only examine bright pixels (avoid noise floor)
    windowed_raw_data = raw_data;
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

%% Show time-domain var
    % 假设 raw_data 大小为 2048 x 2048 x T
    [H, W, T] = size(windowed_raw_data);
    frame_rate = 10;
    time_axis = (1:T) / frame_rate;

    % ==========================================
    % 1. 空间低通滤波 (消除荧光分子的疏密噪声，保留整体照明形态)
    % ==========================================
    sigma = 60.0;  % 与手动拟合保持一致
    smoothed_data = zeros(H, W, T, class(windowed_raw_data));

    if (H * W) > 1e6  % 大图像分块处理
        block_size = 512;
        for t = 1:T
            frame = windowed_raw_data(:, :, t);
            for i = 1:block_size:H
                i_end = min(i + block_size - 1, H);
                for j = 1:block_size:W
                    j_end = min(j + block_size - 1, W);
                    smoothed_data(i:i_end, j:j_end, t) = ...
                        imgaussfilt(frame(i:i_end, j:j_end), sigma);
                end
            end
        end
    else
        for t = 1:T
            smoothed_data(:, :, t) = imgaussfilt(windowed_raw_data(:, :, t), sigma);
        end
    end

    % ==========================================
    % 2. 每帧高斯拟合 - 提取光斑亮度和背景
    % ==========================================
    % 为了加速，对平滑后的数据进行下采样
    downsample_factor = 6;

    signal_intensity = zeros(1, T);
    background_intensity = zeros(1, T);
    signal_std = zeros(1, T);
    background_std = zeros(1, T);
    spot_fwhm_x = zeros(1, T);
    spot_fwhm_y = zeros(1, T);
    fit_quality = zeros(1, T);

    fprintf('Fitting Gaussian for each frame...\n');
    for t = 1:T
        % 下采样平滑数据
        frame_smooth_ds = smoothed_data(1:downsample_factor:H, 1:downsample_factor:W, t);
        [H_ds_actual, W_ds_actual] = size(frame_smooth_ds);
        [xx_ds, yy_ds] = meshgrid(1:W_ds_actual, 1:H_ds_actual);

        try
            est = optics.illumi.fitIlluminationGaussian(frame_smooth_ds);

            % 严格的参数检查：排除不合理的拟合
            % 下采样后的 FWHM（需要乘以下采样因子还原到原始空间）
            fwhm_x_check = 2.355 * est.sx * downsample_factor;
            fwhm_y_check = 2.355 * est.sy * downsample_factor;
            fwhm_ratio = max(fwhm_x_check, fwhm_y_check) / (min(fwhm_x_check, fwhm_y_check) + eps);

            valid_signal = est.A > 0 && est.B >= 0;
            valid_sigma = est.sx > 0 && est.sy > 0;
            % FWHM 应该在合理范围：600-3000 px（对应 σ ~250-1300）
            valid_fwhm = fwhm_x_check > 100 && fwhm_x_check < 3500 && ...
                         fwhm_y_check > 100 && fwhm_y_check < 3500;
            valid_shape = fwhm_ratio < 5;  % 光斑不要太椭圆
            valid_snr = est.A > est.B * 1.1;  % 信号至少比背景高 10%

            if valid_signal && valid_sigma && valid_fwhm && valid_shape && valid_snr
                signal_intensity(t) = est.A;
                background_intensity(t) = est.B;
                spot_fwhm_x(t) = fwhm_x_check;
                spot_fwhm_y(t) = fwhm_y_check;

                fit_img = est.A * exp(-((xx_ds - est.x0).^2) / (2 * est.sx^2) ...
                                     -((yy_ds - est.y0).^2) / (2 * est.sy^2)) + est.B;
                residual = frame_smooth_ds - fit_img;
                fit_quality(t) = sqrt(mean(residual(:).^2)) / (est.A + eps);
                signal_std(t) = sqrt(mean(residual(:).^2));
                background_std(t) = sqrt(mean(residual(:).^2)) * 0.5;
            else
                signal_intensity(t) = NaN;
                background_intensity(t) = NaN;
                spot_fwhm_x(t) = NaN;
                spot_fwhm_y(t) = NaN;
                fit_quality(t) = NaN;
                signal_std(t) = NaN;
                background_std(t) = NaN;
            end
        catch
            signal_intensity(t) = NaN;
            background_intensity(t) = NaN;
            spot_fwhm_x(t) = NaN;
            spot_fwhm_y(t) = NaN;
            fit_quality(t) = NaN;
            signal_std(t) = NaN;
            background_std(t) = NaN;
        end

        if mod(t, max(1, T/10)) == 0
            fprintf('  Frame %d / %d\n', t, T);
        end
    end

    % 计算激发强度（est.A 已经是去背景的峰值）
    excitation_intensity = signal_intensity;  % A = 高斯峰振幅（已去背景）
    excitation_std = signal_std;

    % 数据平滑用于展示趋势
    window_size = 5;
    excitation_smooth = movmean(excitation_intensity, window_size, 'omitnan');
    background_smooth = movmean(background_intensity, window_size, 'omitnan');

    % ==========================================
    % 4. 可视化 - Figure 1: 亮度、背景和光斑大小
    % ==========================================
    figure(1);
    set(gcf, 'Position', [100, 200, 1200, 700], 'Name', 'TIRF: Gaussian Fit Analysis');

    valid_points = ~isnan(excitation_intensity);

    % 上：亮度和背景
    subplot(2, 1, 1);
    hold on;
    errorbar(time_axis(valid_points), excitation_intensity(valid_points), excitation_std(valid_points), 'o', ...
        'Color', [0.8 0.2 0.2], 'LineWidth', 2, 'MarkerSize', 4, 'CapSize', 3, ...
        'DisplayName', 'Excitation Intensity');
    plot(time_axis(valid_points), excitation_smooth(valid_points), '-', ...
        'Color', [0.8 0.2 0.2], 'LineWidth', 3.5, 'DisplayName', 'Excitation Trend');

    errorbar(time_axis(valid_points), background_intensity(valid_points), background_std(valid_points), 's', ...
        'Color', [0.2 0.2 0.8], 'LineWidth', 2, 'MarkerSize', 4, 'CapSize', 3, ...
        'DisplayName', 'Background');
    plot(time_axis(valid_points), background_smooth(valid_points), '-', ...
        'Color', [0.2 0.2 0.8], 'LineWidth', 3.5, 'DisplayName', 'Background Trend');
    hold off;

    ylabel('Intensity (A.U.)', 'FontSize', 12, 'FontWeight', 'bold');
    title('TIRF: Excitation Intensity & Background vs Tilting Angle', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 11);
    grid on; box on;
    set(gca, 'FontSize', 11, 'LineWidth', 1.2);

    % 下：光斑大小（FWHM X 和 Y，还原至原始像素空间）
    subplot(2, 1, 2);
    hold on;
    plot(time_axis(valid_points), spot_fwhm_x(valid_points), 'o-', ...
        'Color', [0.8 0.4 0.1], 'LineWidth', 2, 'MarkerSize', 4, ...
        'DisplayName', 'FWHM X');
    plot(time_axis(valid_points), spot_fwhm_y(valid_points), 's-', ...
        'Color', [0.1 0.4 0.8], 'LineWidth', 2, 'MarkerSize', 4, ...
        'DisplayName', 'FWHM Y');
    hold off;

    xlabel('Time (seconds)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('FWHM (pixels)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Illumination Spot Size (σ) Evolution', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 11);
    grid on; box on;
    set(gca, 'FontSize', 11, 'LineWidth', 1.2);


%% ROI-based intensity measurement (alternative strategy - avoid edge truncation)
    % ==========================================
    % 策略2：中心 ROI 积分亮度（避免边界截断影响）
    % ==========================================
    roi_size = round(200);  % 中心 200x200 像素 ROI
    roi_intensity = zeros(1, T);
    roi_std = zeros(1, T);

    fprintf('Computing center ROI integrated intensity...\n');
    for t = 1:T
        frame = windowed_raw_data(:, :, t);
        [H_frame, W_frame] = size(frame);

        % 定义中心 ROI
        roi_start_y = max(1, floor(H_frame/2) - roi_size/2);
        roi_end_y = min(H_frame, floor(H_frame/2) + roi_size/2);
        roi_start_x = max(1, floor(W_frame/2) - roi_size/2);
        roi_end_x = min(W_frame, floor(W_frame/2) + roi_size/2);

        roi_data = frame(roi_start_y:roi_end_y, roi_start_x:roi_end_x);

        % 积分亮度（总和除以面积，得到平均值）
        roi_intensity(t) = mean(roi_data(:));
        roi_std(t) = std(roi_data(:));
    end

    % 平滑 ROI 亮度趋势
    roi_smooth = movmean(roi_intensity, window_size, 'omitnan');

    % ==========================================
    % 可视化 - Figure 2: ROI 积分亮度
    % ==========================================
    figure(2);
    set(gcf, 'Position', [100, 950, 1100, 500], 'Name', 'TIRF: Center ROI Intensity');

    valid_roi = ~isnan(roi_intensity);

    hold on;
    errorbar(time_axis(valid_roi), roi_intensity(valid_roi), roi_std(valid_roi), 'o', ...
        'Color', [0.2 0.6 0.2], 'LineWidth', 2, 'MarkerSize', 4, 'CapSize', 3, ...
        'DisplayName', sprintf('Center ROI (%dx%d) Mean Intensity', roi_size, roi_size));
    plot(time_axis(valid_roi), roi_smooth(valid_roi), '-', ...
        'Color', [0.2 0.6 0.2], 'LineWidth', 3.5, 'DisplayName', 'Trend');
    hold off;

    xlabel('Time (seconds)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Intensity (A.U.)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('TIRF: Center ROI (%dx%d) Integrated Intensity vs Tilting Angle', roi_size, roi_size), ...
        'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 11);
    grid on; box on;

    set(gca, 'FontSize', 11, 'LineWidth', 1.2);
%% avg
    mean_raw_data = mean(windowed_raw_data(:,:,155:160),3);
    %mean_raw_data_filt = mean(raw_data_filt(:,:,:),3);

%% clear data
    clear raw_data

%% Viz
    figure;
    viz.plotImage(mean_raw_data, 1, 'hot','Mean Image (without filter)');
    %figure;
    %viz.plotImage(mean_raw_data_filt, 1, 'hot','Mean Image (with filtert)');

%% Freq Domain
    F = fftshift(abs(fft2(mean_raw_data)));
    figure;
    viz.plotImage(log(F+1), 1, 'gray','Freq Domain');

%% Lowpass Filter
    L_est = imgaussfilt(mean_raw_data, 10);
    %L_est_filt = imgaussfilt(mean_raw_data_filt, 50);
    figure;
    viz.plotImage(L_est, 1, 'hot','Lowpassed Image (without filter)');
    %figure;
    %viz.plotImage(L_est_filt, 1, 'hot','Lowpassed Image (with filtert)');
    %figure;
    %viz.plotImage(log10(L_est./L_est_filt), 1, 'hot','Lowpassed Image (ratio)');

%% Gaussian Fit
    est = optics.illumi.fitIlluminationGaussian(L_est);
    
%% Fit model
    [ny,nx] = size(L_est);
    [xx,yy] = meshgrid(1:nx,1:ny);

    A  = est.A;
    x0 = est.x0;
    y0 = est.y0;
    sx = est.sx;
    sy = est.sy;
    B  = est.B;

    fit_img = A * exp(-((xx-x0).^2)/(2*sx^2) ...
                 -((yy-y0).^2)/(2*sy^2)) + B;
    figure;
    viz.plotImage(fit_img, 1, 'hot','Fit Image');

%% Residual
    residual_map = L_est - fit_img;
    figure;
    viz.plotImage(residual_map, 1, 'hot','Residual Image');


