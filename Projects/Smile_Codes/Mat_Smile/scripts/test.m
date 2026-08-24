clear;
clc;

% 设置路径
pathway = '/Volumes/SMILeSSD/Optics/TE2000UTests/Cy3/20260407_lhx';
files = dir(fullfile(pathway, '9-cy3_532nm_1p0mW_0P1s_100frames_fov1.sif'));
datanames = {files(~[files.isdir]).name}';

%% 获取数据的数量
datanumber = size(datanames, 1);

% 使用cell数组存储亮度，因为每个文件可能不同数量的分子
brightness = cell(datanumber, 1);
gaussian_params = cell(datanumber, 1); % 存储高斯拟合参数

for x = 1:datanumber
    fprintf('处理文件 %d/%d: %s\n', x, datanumber, datanames{x});
    
    % 设置每个文件的保存名称
    dataname_mat = datanames{x};
    data_save_name = [pathway, 'Brightness\', dataname_mat(1:end-4), '.xlsx'];
    
    % 设置文件读取的路径
    pathname = fullfile(pathway, dataname_mat);
    
    % 读取该文件
    [raw_data, parameters.ex_time, gainDAC] = io.readSIFData(pathname);
    framesize =size(raw_data, 3);
    
    % 确定开始采数据的位置
    startframe = 1;
    for i = 1:framesize
        intmean = mean(mean(raw_data(:,:,i)));
        if intmean > 500
            startframe = i;
            break
        end
    end
    
    % 处理数据
    data = raw_data(128:383, 128:383, startframe:end);
    
    % 调用定位函数，在前xx帧中定位分子
    clear loc_total;
    loc_total = loc_molecule_fast(data(:,:,1:5));
    
    if isempty(loc_total)
        fprintf('  文件中没有找到分子\n');
        continue;
    end
    
    % 为每个分子初始化存储
    num_molecules = size(loc_total, 1);
    temp_brightness = zeros(num_molecules, 1);
    temp_params = cell(num_molecules, 1);
    
    % 对每个分子进行高斯拟合
    for i = 1:num_molecules
        try
            % 获取分子位置
            row = loc_total(i, 1);
            col = loc_total(i, 2);
            frame = loc_total(i, 3);
            
            % 提取分子周围的区域（使用7x7区域以获得更好的拟合）
            roi_size = 7;
            half_size = floor(roi_size/2);
            
            % 确保ROI在图像范围内
            row_start = max(row - half_size, 1);
            row_end = min(row + half_size, size(data, 1));
            col_start = max(col - half_size, 1);
            col_end = min(col + half_size, size(data, 2));
            
            % 提取ROI
            roi = double(data(row_start:row_end, col_start:col_end, frame));
            
            % 如果ROI太小，跳过
            if numel(roi) < 9
                temp_brightness(i) = sum(roi(:));
                continue;
            end
            
            % 准备用于高斯拟合的数据
            [X, Y] = meshgrid(1:size(roi, 2), 1:size(roi, 1));
            Z = roi(:);
            
            % 初始参数估计
            [max_val, max_idx] = max(roi(:));
            [max_row, max_col] = ind2sub(size(roi), max_idx);
            
            % 估计背景（使用ROI边缘的值）
            bg = mean([roi(1,:), roi(end,:), roi(:,1)', roi(:,end)']);
            
            % 设置初始参数: [A, x0, y0, sigma, bg]
            % A: 振幅, x0,y0: 中心位置, sigma: 高斯宽度, bg: 背景
            initial_params = [max_val - bg, max_col, max_row, 1.5, bg];
            
            % 设置参数边界
            lb = [0, 1, 1, 0.5, 0]; % 下限
            ub = [Inf, size(roi, 2), size(roi, 1), 5, max_val]; % 上限
            
            % 定义高斯函数
            gaussian2d = @(params, X, Y) params(5) + params(1) * ...
                exp(-((X - params(2)).^2 + (Y - params(3)).^2) / (2 * params(4)^2));
            
            % 定义拟合函数
            fit_func = @(params) gaussian2d(params, X, Y) - roi;
            
            % 使用lsqnonlin进行非线性最小二乘拟合
            options = optimoptions('lsqnonlin', ...
                'Display', 'off', ...
                'MaxIterations', 100, ...
                'FunctionTolerance', 1e-8);
            
            % 执行拟合
            params_fit = lsqnonlin(fit_func, initial_params, lb, ub, options);
            
            % 保存拟合参数
            temp_params{i} = params_fit;
            
            % 计算高斯拟合的积分作为亮度（2D高斯积分 = A * 2πσ²）
            brightness_fit = params_fit(1) * 2 * pi * params_fit(4)^2;
            
            % 保存亮度
            temp_brightness(i) = brightness_fit;
            
            % 可选：显示拟合结果（调试用）
            % if mod(i, 10) == 0
            %     figure;
            %     subplot(1,2,1);
            %     imagesc(roi);
            %     title('原始ROI');
            %     subplot(1,2,2);
            %     fitted_roi = gaussian2d(params_fit, X, Y);
            %     imagesc(fitted_roi);
            %     title('高斯拟合');
            %     sgtitle(sprintf('分子 %d/%d', i, num_molecules));
            % end
            
        catch ME
            fprintf('  分子 %d 拟合失败: %s\n', i, ME.message);
            % 如果拟合失败，使用简单求和作为备份
            temp_brightness(i) = sum(roi(:));
        end
    end
    
    % 保存到cell数组
    brightness{x} = temp_brightness;
    gaussian_params{x} = temp_params;
    
    % 可选：保存结果到文件
    if ~exist([pathway, 'Brightness\'], 'dir')
        mkdir([pathway, 'Brightness']);
    end
    
    % 保存亮度数据
    brightness_table = table(temp_brightness, 'VariableNames', {'Brightness'});
    writetable(brightness_table, data_save_name, 'Sheet', 'Brightness');
    
    % 可选：保存高斯参数
    param_names = {'Amplitude', 'X_center', 'Y_center', 'Sigma', 'Background'};
    param_data = zeros(num_molecules, 5);
    for i = 1:num_molecules
        if ~isempty(temp_params{i})
            param_data(i, :) = temp_params{i};
        end
    end
    param_table = array2table(param_data, 'VariableNames', param_names);
    writetable(param_table, data_save_name, 'Sheet', 'Gaussian_Params');
    
    fprintf('  处理完成，找到 %d 个分子\n', num_molecules);
end

%% Viz
    frames_to_show = 1:min(1, 5);

    figure;
    viz.plotImage(data, frames_to_show, 'hot', ...
        'Drift-uncorrected Image and Detections');
    hold on;
    viz.overlayLocs(loc_total, frames_to_show);

%% 分析所有文件的亮度分布
all_brightness = [];
for x = 1:datanumber
    if ~isempty(brightness{x})
        all_brightness = [all_brightness; brightness{x}];
    end
end

if ~isempty(all_brightness)
    figure;
    plot(1);
    histogram(all_brightness, 50);
    xlabel('亮度');
    ylabel('计数');
    title('所有分子亮度分布');
    
    %subplot(2,2,2);
    %boxplot(all_brightness);
    %ylabel('亮度');
    %title('亮度箱线图');
    
    %subplot(2,2,3);
    %plot(all_brightness, 'o');
    %xlabel('分子序号');
    %ylabel('亮度');
    %title('亮度散点图');
    
    %subplot(2,2,4);
    %[f, xi] = ksdensity(all_brightness);
    %plot(xi, f);
    %xlabel('亮度');
    %ylabel('概率密度');
    %title('亮度概率密度估计');
    
    % 显示统计信息
    fprintf('\n亮度统计信息:\n');
    fprintf('  均值: %.2f\n', mean(all_brightness));
    fprintf('  中位数: %.2f\n', median(all_brightness));
    fprintf('  标准差: %.2f\n', std(all_brightness));
    fprintf('  最小值: %.2f\n', min(all_brightness));
    fprintf('  最大值: %.2f\n', max(all_brightness));
end