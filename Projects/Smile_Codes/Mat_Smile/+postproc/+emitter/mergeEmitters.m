function emitters = mergeEmitters(emitters, r_merge)
% 合并空间上重叠的 emitter
% r_merge: 合并半径，建议取 r_track 或略大（如 1.5 * r_track）

if nargin < 2
    r_merge = 1.5;
end

N = numel(emitters);
merged = false(1, N);

% 计算每个 emitter 的平均位置
mean_pos = zeros(N, 2);
for k = 1:N
    mean_pos(k,:) = [mean(emitters(k).row), mean(emitters(k).col)];
end

% 用 pdist2 找空间邻近对
D = pdist2(mean_pos, mean_pos);
D(logical(eye(N))) = Inf;  % 排除自身

for k = 1:N
    if merged(k)
        continue
    end
    
    % 找所有与 k 空间重叠且未被合并的 emitter
    neighbors = find(D(k,:) < r_merge & ~merged);
    neighbors(neighbors == k) = [];
    
    for j = neighbors
        if merged(j)
            continue
        end
        
        % 检查时间上是否有大量重叠帧（同帧同 emitter 不合理）
        overlap_frames = intersect(emitters(k).frames, emitters(j).frames);
        total_frames = union(emitters(k).frames, emitters(j).frames);
        overlap_ratio = numel(overlap_frames) / numel(total_frames);
        
        if overlap_ratio > 0.3
            % 重叠帧太多 → 可能真的是两个紧邻 emitter，不合并
            continue
        end
        
        % 合并 j → k
        emitters(k).row    = [emitters(k).row,    emitters(j).row];
        emitters(k).col    = [emitters(k).col,    emitters(j).col];
        emitters(k).frames = [emitters(k).frames, emitters(j).frames];
        emitters(k).loc_idx = [emitters(k).loc_idx, emitters(j).loc_idx];
        
        % 按帧排序
        [emitters(k).frames, sort_i] = sort(emitters(k).frames);
        emitters(k).row    = emitters(k).row(sort_i);
        emitters(k).col    = emitters(k).col(sort_i);
        emitters(k).loc_idx = emitters(k).loc_idx(sort_i);
        
        % 更新时间信息
        emitters(k).on_frame     = emitters(k).frames(1);
        emitters(k).last_frame   = emitters(k).frames(end);
        emitters(k).bleach_frame = emitters(k).last_frame;
        emitters(k).alive        = emitters(k).alive || emitters(j).alive;
        
        merged(j) = true;
        
        % 更新 k 的平均位置供后续比较
        mean_pos(k,:) = [mean(emitters(k).row), mean(emitters(k).col)];
        D(k,:) = pdist2(mean_pos(k,:), mean_pos);
        D(:,k) = D(k,:)';
        D(k,k) = Inf;
    end
end

% 删除被合并掉的条目
emitters(merged) = [];
fprintf('Left %d emitters(merging %d）\n', numel(emitters), sum(merged));
end