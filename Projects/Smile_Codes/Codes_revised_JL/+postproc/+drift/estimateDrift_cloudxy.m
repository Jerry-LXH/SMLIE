function [delta_sum] = estimateDrift_cloudxy(raw_data, loc_total, frames_per_seg, uncertainty, matched_pairs_threshold, do_stage2)
% estimateDrift_cloudxy  基于 NP-Cloud 算法估计 XY 方向的逐帧漂移
%
% ============ 功能简介 ============
%
%   在超分辨显微成像中，样品会在采集过程中缓慢漂移（通常几十到几百纳米）。
%   如果不校正，重构图像会变模糊。本函数将全部帧分成若干"时间段"（segment），
%   利用相邻段之间分子定位坐标的最近邻匹配来估计漂移量，最终输出每一帧的
%   累积漂移。
%
%   算法基于 Nature Communications (2025) 16:9031 的 NP-Cloud 简化版，
%   核心思路如下：
%     1) 将全部帧按 frames_per_seg 分成若干段；
%     2) 以第一段的定位为参考，对后续每段在搜索半径 r = 3×uncertainty 内
%        做最近邻匹配，构建"位移向量云"（displacement cloud）；
%     3) 对位移向量取均值得到该段漂移估计，并迭代至收敛；
%     4) 将段级漂移线性分配到每帧。
%
%   可选的 Stage 2（RR-NP Cloud）进一步利用 Stage 1 校正后的全部定位构建
%   一个信息量更丰富的"增强参考点云"，再对每段重新估计漂移，结果更鲁棒。
%
% ============ 输入参数 ============
%
%   raw_data                — 图像栈 [rows × cols × frames]，仅用于获取帧数。
%   loc_total               — 分子定位 [N × 3]，各列为 (row, col, frame)。
%   frames_per_seg          — 每段帧数，推荐 10；闪烁多时可设 15–30。
%   uncertainty             — 定位不确定度（像素），搜索半径 r = 3 × uncertainty。
%   matched_pairs_threshold — 最少匹配对数（默认 30），不足则跳过该段。
%   do_stage2               — 是否执行 Stage 2（RR-NP Cloud），默认 false。
%
% ============ 输出参数 ============
%
%   delta_sum — [frames × 2]，每帧累积漂移 (row_drift, col_drift)。
%              在 correctDrift 中用此量对图像和定位做校正。
%
% Written by Jerry Ling, 2026.1.20

% ---- 默认参数 ----
if nargin < 5 || isempty(matched_pairs_threshold)
    matched_pairs_threshold = 30;
end
if nargin < 6 || isempty(do_stage2)
    do_stage2 = false;
end

% ---- 基本参数 ----
[~, ~, frames] = size(raw_data);
delta_sum = zeros(frames, 2);
r         = 1.5 * 3 * uncertainty;                    % 最近邻搜索半径
min_locs  = 2 * matched_pairs_threshold;        % 段内最少定位数
end_seg   = ceil(frames / frames_per_seg);       % 总段数
resample_factor = 12;                            % Stage 2 重采样倍数（论文推荐值）

tic;

% ======================================================================
%  Stage 1：以第一段为参考，逐段估计累积漂移
% ======================================================================
fprintf('\n=============== Stage 1: Use first segment as reference ===============\n');

% ---- 提取第一段定位作为参考点云 Q ----
idx_q = ismember(loc_total(:,3), 1:frames_per_seg);
Q     = loc_total(idx_q, 1:2);
fprintf('Reference segment (seg 1, frames 1-%d): %d locs\n', frames_per_seg, size(Q,1));

if size(Q,1) < min_locs
    fprintf('[Warn] Seg 1 has only %d locs < threshold %d, estimate may be unstable.\n', ...
        size(Q,1), min_locs);
end

% ---- 逐段迭代估计 ----
sum_trace        = [0, 0];   % 段级累积漂移（跨段持续累加）
iter_sum_s1      = 0;        % 总迭代次数（用于统计）
skip_count_s1    = 0;        % 跳过的段数
precision_list_s1 = [];      % 各有效段的理论校正精度

for seg = 2:end_seg

    start_frame  = (seg-1) * frames_per_seg + 1;
    end_frame    = min(seg * frames_per_seg, frames);
    n_frames_seg = end_frame - start_frame + 1;        % 兼容最后一个不完整段
    idx_r        = ismember(loc_total(:,3), start_frame:end_frame);

    next_seg = false;   % 标记是否需要跳过该段

    for iter = 1:20
        % 用当前累积漂移预校正该段定位
        R = loc_total(idx_r, 1:2) - sum_trace;

        if size(R,1) < min_locs
            fprintf('[Warn] Segment %d（Frame %d–%d）locs %d < %d, skipped.\n', ...
                seg, start_frame, end_frame, size(R,1), min_locs);
            next_seg = true;
            break
        end

        % 对参考 Q 的每个点，在 R 中搜索半径 r 内的最近邻
        [nnIdx, nnDist] = knnsearch(R, Q, 'K', 1);
        mask          = nnDist <= r;
        matched_pairs = sum(mask);

        if matched_pairs < matched_pairs_threshold
            fprintf('[Warn] Segment %d（Frame %d–%d）matched pairs %d < %d, skipped.\n', ...
                seg, start_frame, end_frame, matched_pairs, matched_pairs_threshold);
            next_seg = true;
            break
        end

        % 构建位移向量云：R(matched) − Q(matched) = 残余漂移
        dxy = R(nnIdx(mask), :) - Q(mask, :);

        % 收敛判据：期望精度 ≈ uncertainty / sqrt(matched_pairs)
        if iter == 1
            judge_radius = uncertainty / sqrt(matched_pairs);
        end

        delta_y   = mean(dxy(:,1));
        delta_x   = mean(dxy(:,2));
        sum_trace = sum_trace + [delta_y, delta_x];
        centeroff = hypot(delta_y, delta_x);

        if centeroff < judge_radius
            fprintf('  Segment %2d (Frame %4d–%4d): %2d iterations, %4d pairs, residual %.4f px\n', ...
                seg, start_frame, end_frame, iter-1, matched_pairs, centeroff);
            break
        end
    end

    % ---- 跳过的段：沿用上一帧漂移（零增量外推） ----
    if next_seg
        skip_count_s1 = skip_count_s1 + 1;
        for f = start_frame:end_frame
            delta_sum(f,:) = delta_sum(f-1,:);
        end
        continue
    end

    iter_sum_s1 = iter_sum_s1 + (iter - 1);
    precision_list_s1 = [precision_list_s1; uncertainty / sqrt(matched_pairs)];

    % ---- 将段级漂移线性分配到该段各帧 ----
    frame_drift = (sum_trace - delta_sum(start_frame-1, :)) / n_frames_seg;
    for f = start_frame:end_frame
        delta_sum(f,:) = delta_sum(f-1,:) + frame_drift;
    end
end

% ---- Stage 1 报告 ----
t_s1 = toc;
valid_segs_s1  = (end_seg - 1) - skip_count_s1;
total_drift_s1 = hypot(delta_sum(end,1), delta_sum(end,2));

fprintf('\n--- Stage 1 Summary ---\n');
fprintf('Segments: %d total (%d valid, %d skipped)\n', end_seg-1, valid_segs_s1, skip_count_s1);
if valid_segs_s1 > 0
    fprintf('Avg iterations: %.2f per seg\n', iter_sum_s1 / valid_segs_s1);
    fprintf('Precision (sigma/sqrt(n)): median %.4f px, mean %.4f px\n', ...
        median(precision_list_s1), mean(precision_list_s1));
end
fprintf('Total drift: %.2f px (dRow %.2f, dCol %.2f)\n', ...
    total_drift_s1, delta_sum(end,1), delta_sum(end,2));
fprintf('Elapsed: %.3f s (%.3f ms/seg)\n', t_s1, t_s1/max(end_seg-1,1)*1000);

% ======================================================================
%  Stage 2（可选）：RR-NP Cloud — 以全局重采样参考重新估计
% ======================================================================
if do_stage2
    fprintf('\n========= Stage 2: Use Full Loc Set as Reference =========\n');
    tic_s2 = tic;

    % 备份 Stage 1 结果，以便 Stage 2 失败时回退
    delta_sum_s1 = delta_sum;

    % ---- 2.1 用 Stage 1 结果校正全部定位 ----
    %   对每个定位点减去其所在帧的累积漂移，得到初步校正坐标。
    frame_list   = loc_total(:, 3);
    corrected_xy = loc_total(:, 1:2) - delta_sum_s1(frame_list, :);

    % ---- 2.2 从各段均匀采样，构建增强参考点云 Q2 ----
    %   为保证 Q2 在时间上均匀分布，从每段取相同数量的定位点。
    %   目标总数 ≈ resample_factor × 第一段定位数。
    n_target  = resample_factor * size(Q, 1);
    n_per_seg = ceil(n_target / end_seg);
    Q2 = [];
    for seg = 1:end_seg
        sf = (seg-1)*frames_per_seg + 1;
        ef = min(seg*frames_per_seg, frames);
        idx_seg = ismember(frame_list, sf:ef);
        seg_xy  = corrected_xy(idx_seg, :);
        if size(seg_xy,1) <= n_per_seg
            Q2 = [Q2; seg_xy];                          %#ok<AGROW>
        else
            sel = randperm(size(seg_xy,1), n_per_seg);
            Q2  = [Q2; seg_xy(sel,:)];                  %#ok<AGROW>
        end
    end
    fprintf('Enhanced ref Q2: %d locs (approx %d x %d = target %d)\n', ...
        size(Q2,1), size(Q,1), resample_factor, n_target);

    % ---- 2.3 对每段独立估计绝对漂移 ----
    %   与 Stage 1 不同，这里每段独立估计（不跨段累加），
    %   并以 Stage 1 的估计值为迭代起点加速收敛。
    seg_drift          = zeros(end_seg, 2);   % 各段绝对漂移
    seg_valid          = true(end_seg, 1);    % 是否有效
    skip_count_s2      = 0;
    iter_sum_s2        = 0;
    precision_list_s2  = [];

    for seg = 1:end_seg
        start_frame  = (seg-1)*frames_per_seg + 1;
        end_frame    = min(seg*frames_per_seg, frames);
        idx_r        = ismember(frame_list, start_frame:end_frame);

        % 以 Stage 1 中该段中心帧的漂移值为初始猜测
        mid_frame  = round((start_frame + end_frame) / 2);
        sum_trace2 = delta_sum_s1(mid_frame, :);

        next_seg = false;

        for iter = 1:20
            % 用当前猜测预校正该段定位
            R = loc_total(idx_r, 1:2) - sum_trace2;

            if size(R,1) < min_locs
                fprintf('[Warn] Segment %d（Frame %d–%d）locs %d < %d, skipped.\n', ...
                    seg, start_frame, end_frame, size(R,1), min_locs);
                next_seg = true;
                break
            end

            % 对 R 中每个点，在 Q2 中搜索最近邻
            %   注意方向：knnsearch(Q2, R) 以 Q2 为数据库、R 为查询
            [nnIdx2, nnDist2] = knnsearch(Q2, R, 'K', 1);
            mask2         = nnDist2 <= r;
            matched_pairs = sum(mask2);

            if matched_pairs < matched_pairs_threshold
                fprintf('[Warn] Segment %d（Frame %d–%d）matched pairs %d < %d, skipped.\n', ...
                    seg, start_frame, end_frame, matched_pairs, matched_pairs_threshold);
                next_seg = true;
                break
            end

            % 位移向量：R(matched) − Q2(nearest) = 残余漂移
            dxy = R(mask2, :) - Q2(nnIdx2(mask2), :);

            if iter == 1
                judge_radius2 = uncertainty / sqrt(matched_pairs);
            end

            delta_y    = mean(dxy(:,1));
            delta_x    = mean(dxy(:,2));
            sum_trace2 = sum_trace2 + [delta_y, delta_x];
            centeroff  = hypot(delta_y, delta_x);

            if centeroff < judge_radius2
                fprintf(' Segment %2d (Frame %4d–%4d): %2d iterations, %4d pairs, residual %.4f px\n', ...
                    seg, start_frame, end_frame, iter-1, matched_pairs, centeroff);
                break
            end
        end

        if next_seg
            skip_count_s2    = skip_count_s2 + 1;
            seg_valid(seg)   = false;
            seg_drift(seg,:) = delta_sum_s1(mid_frame, :);   % 回退到 Stage 1
        else
            iter_sum_s2      = iter_sum_s2 + (iter - 1);
            seg_drift(seg,:) = sum_trace2;
            precision_list_s2 = [precision_list_s2; uncertainty / sqrt(matched_pairs)]; %#ok<AGROW>
        end
    end

    % ---- 2.4 将段级漂移插值到每帧 ----
    %   用各段中心帧位置做线性插值，仅使用有效段，跳过段由插值自动覆盖。
    seg_mids = zeros(end_seg, 1);
    for seg = 1:end_seg
        sf = (seg-1)*frames_per_seg + 1;
        ef = min(seg*frames_per_seg, frames);
        seg_mids(seg) = (sf + ef) / 2;
    end

    if sum(seg_valid) >= 2
        delta_sum(:,1) = interp1(seg_mids(seg_valid), seg_drift(seg_valid,1), (1:frames)', 'linear', 'extrap');
        delta_sum(:,2) = interp1(seg_mids(seg_valid), seg_drift(seg_valid,2), (1:frames)', 'linear', 'extrap');
    elseif sum(seg_valid) == 1
        delta_sum = repmat(seg_drift(seg_valid,:), frames, 1);
        fprintf('[Warn] Only 1 valid segment, cannot interpolate; using constant drift for all frames.\n');
    else
        delta_sum = delta_sum_s1;
        fprintf('[Warn] No valid segments, Stage 2 failed; keeping Stage 1 result.\n');
    end

    % 归一化：令第一帧漂移为 [0,0]，保持与 Stage 1 一致的参考系
    drift_offset = delta_sum(1, :);
    delta_sum    = delta_sum - drift_offset;

    % ---- Stage 2 报告 ----
    t_s2 = toc(tic_s2);
    valid_segs_s2  = sum(seg_valid);
    total_drift_s2 = hypot(delta_sum(end,1), delta_sum(end,2));
    fprintf('\n--- Stage 2 Summary ---\n');
    fprintf('Segments: %d total (%d valid, %d skipped)\n', end_seg, valid_segs_s2, skip_count_s2);
    if valid_segs_s2 > 0 && ~isempty(precision_list_s2)
        fprintf('Avg iterations: %.2f per seg\n', iter_sum_s2 / valid_segs_s2);
        fprintf('Precision (sigma/sqrt(n)): median %.4f px, mean %.4f px\n', ...
            median(precision_list_s2), mean(precision_list_s2));
    end
    fprintf('Total drift: %.2f px (dRow %.2f, dCol %.2f)\n', ...
        total_drift_s2, delta_sum(end,1), delta_sum(end,2));
    fprintf('Elapsed: %.3f s (%.3f ms/seg)\n', t_s2, t_s2/max(end_seg,1)*1000);

    % ---- Stage comparison ----
    if ~isempty(precision_list_s1) && ~isempty(precision_list_s2)
        p1 = median(precision_list_s1);
        p2 = median(precision_list_s2);
        fprintf('\n--- Stage 1 vs Stage 2 ---\n');
        fprintf('Skipped segs: %d -> %d\n', skip_count_s1, skip_count_s2);
        fprintf('Median precision: %.4f px -> %.4f px', p1, p2);
        if p1 > 0
            fprintf(' (%+.1f%%)', (p2/p1 - 1)*100);
        end
        fprintf('\n');
    end
end

% ---- 总体输出 ----
t_total = toc;
fprintf('\n=============== Drift Estimation Complete ===============\n');
fprintf('Frames: %d, Segments: %d (%d frames/seg)\n', frames, end_seg, frames_per_seg);
fprintf('Final drift: %.2f px (dRow %.2f, dCol %.2f)\n', ...
    hypot(delta_sum(end,1), delta_sum(end,2)), delta_sum(end,1), delta_sum(end,2));
fprintf('Total time: %.3f s\n\n', t_total);

end