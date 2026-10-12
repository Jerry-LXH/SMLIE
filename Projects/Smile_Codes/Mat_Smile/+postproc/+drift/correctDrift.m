function [corr_data, corr_locs, corr_idx] = correctDrift(raw_data, loc_total, delta_sum, edge, row_width, col_width)
% correctDrift  根据估计的漂移量，校正图像栈和分子定位坐标
%
% ============ 功能简介 ============
%
%   超分辨成像过程中，样品会随时间发生缓慢漂移。本函数做三件事：
%     1) 逐帧平移图像以抵消漂移（亚像素精度，双线性插值）；
%     2) 从平移后的图像中，居中裁切出 row_width × col_width 大小的区域；
%     3) 对分子定位坐标做相同的漂移校正和坐标转换，
%        只保留落在裁切区域内（且远离边缘 edge 像素）的点。
%
%   因此，raw_data 的尺寸必须 > 输出尺寸，多出的部分作为平移缓冲区。
%
% ============ 输入参数 ============
%
%   raw_data   — 原始图像栈 [rows × cols × frames]，需预先居中裁切，
%                尺寸必须大于 [row_width × col_width]。
%   loc_total  — 分子定位结果 [N × 3]，各列为 (row, col, frame)。
%   delta_sum  — 各帧累积漂移 [frames × 2]，各列为 (row_drift, col_drift)。
%   edge       — 边缘保护距离（像素），默认 5。
%   row_width  — 输出图像行数，默认 257。
%   col_width  — 输出图像列数，默认与 row_width 相同。
%
% ============ 输出参数 ============
%
%   corr_data  — 校正后图像栈 [row_width × col_width × frames]。
%   corr_locs  — 校正后定位坐标 [M × 3]（局部坐标系）。
%   corr_idx   — corr_locs 各点在 loc_total 中的原始行号，便于回溯。

% ---- 默认参数 ----
if nargin < 4 || isempty(edge),      edge      = 5;          end
if nargin < 5 || isempty(row_width), row_width = 257;        end
if nargin < 6 || isempty(col_width), col_width = row_width;  end

% ---- 读取输入尺寸 ----
[rows, cols, frames] = size(raw_data);

% ---- 计算居中裁切范围（兼容非整除） ----
margin_row = rows - row_width;
margin_col = cols - col_width;
if margin_row < 0
    error('raw_data rows (%d) < row_width (%d), unable to window data.', rows, row_width);
end
if margin_col < 0
    error('raw_data cols (%d) < col_width (%d), unable to window data.', cols, col_width);
end

start_row = floor(margin_row / 2) + 1;
end_row   = start_row + row_width - 1;

start_col = floor(margin_col / 2) + 1;
end_col   = start_col + col_width - 1;

% ---- Odd margin warning ----
if mod(margin_row, 2) ~= 0
    fprintf('[Warn] Row margin %d is odd: top = %d, bottom = %d\n', ...
        margin_row, start_row - 1, rows - end_row);
end
if mod(margin_col, 2) ~= 0
    fprintf('[Warn] Col margin %d is odd: left = %d, right = %d\n', ...
        margin_col, start_col - 1, cols - end_col);
end

% ---- Print layout info ----
fprintf('Input  : %d x %d x %d frames\n', rows, cols, frames);
fprintf('Output : rows [%d:%d], cols [%d:%d] -> %d x %d\n', ...
    start_row, end_row, start_col, end_col, row_width, col_width);

% ---- 漂移缓冲区充足性检查 ----
% 四个方向的可用缓冲区宽度（像素）
top_margin    = start_row - 1;
bottom_margin = rows - end_row;
left_margin   = start_col - 1;
right_margin  = cols - end_col;

% imtranslate 的平移量为 -delta_sum，输出裁切窗口内像素实际采样自
%   原图的 [start_row + drift_row, end_row + drift_row] 区间（列方向同理）。
% 因此需要: drift_row >= -top_margin  且 drift_row <= bottom_margin
%           drift_col >= -left_margin 且 drift_col <= right_margin
drift_row = delta_sum(:, 1);
drift_col = delta_sum(:, 2);

% 各方向漂移极值（正值表示确实存在该方向的漂移需求）
max_drift_up    = max(0, -min(drift_row));   % 最大上漂（消耗 top_margin）
max_drift_down  = max(0,  max(drift_row));   % 最大下漂（消耗 bottom_margin）
max_drift_left  = max(0, -min(drift_col));   % 最大左漂（消耗 left_margin）
max_drift_right = max(0,  max(drift_col));   % 最大右漂（消耗 right_margin）

% 逐帧判断是否溢出
overflow_top    = drift_row < -top_margin;
overflow_bottom = drift_row >  bottom_margin;
overflow_left   = drift_col < -left_margin;
overflow_right  = drift_col >  right_margin;
overflow_any    = overflow_top | overflow_bottom | overflow_left | overflow_right;

n_overflow = sum(overflow_any);

if n_overflow > 0
    % 各方向的溢出量（超出缓冲的像素数）
    excess_top    = max(0, max_drift_up    - top_margin);
    excess_bottom = max(0, max_drift_down  - bottom_margin);
    excess_left   = max(0, max_drift_left  - left_margin);
    excess_right  = max(0, max_drift_right - right_margin);

    % 推荐的最小输入尺寸：输出尺寸 + 两侧漂移极值之和（向上取整留余量）
    min_rows = row_width + ceil(max_drift_up) + ceil(max_drift_down);
    min_cols = col_width + ceil(max_drift_left) + ceil(max_drift_right);

    warning('correctDrift:bufferInsufficient', ...
        ['\n', ...
         '========== 漂移缓冲区不足 ==========\n', ...
         '共 %d / %d 帧（%.1f%%）的漂移超出缓冲范围，\n', ...
         '这些帧的裁切区域将包含零填充像素（黑边），可能影响后续分析。\n', ...
         '\n', ...
         '  方向      可用缓冲    实际需求    溢出量     溢出帧数\n', ...
         '  上 (↑)    %6.1f px    %6.2f px    %5.2f px    %d\n', ...
         '  下 (↓)    %6.1f px    %6.2f px    %5.2f px    %d\n', ...
         '  左 (←)    %6.1f px    %6.2f px    %5.2f px    %d\n', ...
         '  右 (→)    %6.1f px    %6.2f px    %5.2f px    %d\n', ...
         '\n', ...
         '建议操作（任选其一）:\n', ...
         '  (1) 增大 raw_data 的裁切范围（至少 %d × %d）\n', ...
         '  (2) 缩小 row_width / col_width\n', ...
         '  (3) 若溢出帧集中在首尾，可截断该部分帧后重新估计漂移\n', ...
         '==========================================\n'], ...
        n_overflow, frames, 100 * n_overflow / frames, ...
        top_margin,    max_drift_up,    excess_top,    sum(overflow_top), ...
        bottom_margin, max_drift_down,  excess_bottom, sum(overflow_bottom), ...
        left_margin,   max_drift_left,  excess_left,   sum(overflow_left), ...
        right_margin,  max_drift_right, excess_right,  sum(overflow_right), ...
        min_rows, min_cols);
else
    fprintf(['Buffer check passed:\n', ...
         '  Max drift   U %.2f / D %.2f / L %.2f / R %.2f px\n', ...
         '  Available   U %d / D %d / L %d / R %d px\n'], ...
    max_drift_up, max_drift_down, max_drift_left, max_drift_right, ...
    top_margin, bottom_margin, left_margin, right_margin);
end

% ---- 初始化输出 ----
corr_data = zeros(row_width, col_width, frames, 'like', raw_data);
corr_locs = [];
corr_idx  = [];

% ---- 逐帧校正 ----
for i = 1:frames
    % 亚像素平移（imtranslate 格式为 [x, y] = [col, row]）
    img_corr = imtranslate(raw_data(:,:,i), -delta_sum(i,[2 1]), 'linear', 'FillValues', 0);

    % 居中裁切
    corr_data(:,:,i) = img_corr(start_row:end_row, start_col:end_col);

    % ---- 校正该帧的定位坐标 ----
    idx = loc_total(:,3) == i;
    if ~any(idx), continue; end
    idx_global = find(idx);

    % 减去漂移
    locs = loc_total(idx, 1:2) - delta_sum(i,:);

    % 转换为裁切窗口的局部坐标
    locs(:,1) = locs(:,1) - start_row + 1;
    locs(:,2) = locs(:,2) - start_col + 1;

    % 只保留远离边缘的点
    in = locs(:,1) >= 1 + edge & locs(:,1) <= row_width - edge & ...
         locs(:,2) >= 1 + edge & locs(:,2) <= col_width - edge;

    if any(in)
        corr_locs = [corr_locs; locs(in,:), i * ones(sum(in),1)];
        corr_idx  = [corr_idx;  idx_global(in)];
    end
end

fprintf('Drift Correction complete. Output size: [%d × %d × %d], remaining locs: %d / %d。\n', ...
    row_width, col_width, frames, size(corr_locs,1), size(loc_total,1));

end