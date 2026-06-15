function emitters = findEmitters(data, loc_total, max_gap, r_track)
% FINDEMITTERS  Extract emitters from a localization series.
%
% Clusters localizations into emitters based on time-gap and spatial
% distance.  Uses global optimal assignment (matchpairs) instead of
% greedy nearest-neighbor to link localizations each frame.

% INPUTS
%   data      - image stack [rows x cols x frames], only frame count is used.
%               Alternatively, a scalar integer giving total number of frames.
%   loc_total - [N x 3] matrix: (row_loc, col_loc, frame)
%   max_gap   - max dark frames before declaring bleach (default: 30)
%   r_track   - spatial linking radius, same units as locs (default: 1)
%
% OUTPUTS
%   emitters(k):
%     .row          - [double array]  row positions
%     .col          - [double array]  col positions
%     .frames       - [int array]     frame numbers
%     .loc_idx      - [int array]     row indices into loc_total
%     .on_frame     - int, first detected frame
%     .bleach_frame - int, = last_frame (for alive emitters too)
%     .last_frame   - int, last detected frame
%     .alive        - logical, true if never declared bleached

% ---- Default parameters ----
if nargin < 3 || isempty(max_gap)
    max_gap = 30;
end
if nargin < 4 || isempty(r_track)
    r_track = 1;
end

% ---- Determine total number of frames ----
if isscalar(data)
    n_frames = double(data);
else
    [~, ~, n_frames] = size(data);
end

% ---- Pre-group localizations by frame ----
% This replaces the per-frame "loc_total(:,3) == f" scan in the original,
% reducing overall complexity from O(N * F) to O(N log N) for this step.
loc_total(:,3) = round(loc_total(:,3));
[locs_by_frame, loc_idx_by_frame] = groupLocsByFrame(loc_total, n_frames);

% ---- Initialize emitter storage ----
emitters     = struct([]);
num_emitters = 0;
alive_list   = [];       % indices of currently alive emitters

tic;
for f = 1:n_frames

    % ---- Close dead emitters (gap exceeded max_gap) ----
    dead_mask = false(size(alive_list));
    for j = 1:numel(alive_list)
        k = alive_list(j);
        if f - emitters(k).last_frame > max_gap
            emitters(k).alive = false;
            emitters(k).bleach_frame = emitters(k).last_frame;
            dead_mask(j) = true;
        end
    end
    alive_list(dead_mask) = [];

    % ---- Get localizations for current frame ----
    locs      = locs_by_frame{f};       % [Nf x 2] positions
    loc_idx_f = loc_idx_by_frame{f};    % [1 x Nf]  indices into loc_total
    if isempty(locs)
        if mod(f, 10) == 0 || f == n_frames
            fprintf('Progress: %.1f%%\n', (f / n_frames) * 100);
        end
        continue
    end

    matched = false(size(locs, 1), 1);  % which locs are matched

    % ---- Match localizations to alive emitters via matchpairs ----
    if ~isempty(alive_list)
        n_alive = numel(alive_list);

        % Collect last-known positions of all alive emitters
        prev_pos = zeros(n_alive, 2);
        for j = 1:n_alive
            k = alive_list(j);
            prev_pos(j, :) = [emitters(k).row(end), emitters(k).col(end)];
        end

        % Pairwise distance matrix [n_locs_f x n_alive]
        d = sqrt((locs(:,1) - prev_pos(:,1)').^2 + ...
                 (locs(:,2) - prev_pos(:,2)').^2);

        % Disallow links beyond r_track
        d(d >= r_track) = Inf;

        % Global optimal 1-to-1 assignment (Jonker-Volgenant via matchpairs)
        if any(isfinite(d(:)))
            pairs = matchpairs(d, r_track);   % [M x 2]: (loc_row, alive_col)

            for j = 1:size(pairs, 1)
                i          = pairs(j, 1);             % index into locs
                k          = alive_list(pairs(j, 2)); % emitter index
                idx_global = loc_idx_f(i);

                % Append to existing emitter
                emitters(k).row(end+1)     = locs(i, 1);
                emitters(k).col(end+1)     = locs(i, 2);
                emitters(k).frames(end+1)  = f;
                emitters(k).loc_idx(end+1) = idx_global;
                emitters(k).last_frame     = f;

                matched(i) = true;
            end
        end
    end

    % ---- Create new emitters for unmatched localizations ----
    for i = find(~matched)'
        pos        = locs(i, :);
        idx_global = loc_idx_f(i);

        k = num_emitters + 1;
        emitters(k).row         = pos(1);
        emitters(k).col         = pos(2);
        emitters(k).frames      = f;
        emitters(k).loc_idx     = idx_global;
        emitters(k).on_frame    = f;
        emitters(k).bleach_frame = NaN;
        emitters(k).last_frame  = f;
        emitters(k).alive       = true;

        num_emitters      = k;
        alive_list(end+1) = k;                    %#ok<AGROW>
    end

    % ---- Progress report ----
    if mod(f, 10) == 0 || f == n_frames
        fprintf('Progress: %.1f%%\n', (f / n_frames) * 100);
    end
end

% ---- Finalize: set bleach_frame for emitters still alive at the end ----
N = numel(emitters);
for k = 1:N
    if emitters(k).alive
        emitters(k).bleach_frame = emitters(k).last_frame;
    end
end

t_total = toc;
fprintf('Total time = %.3f s\n', t_total);
fprintf('Totally %d emitters found.\n', N);

end


%% =====================================================================
%  Helper function: pre-group localizations by frame
%  =====================================================================
function [locs_by_frame, idx_by_frame] = groupLocsByFrame(loc_total, n_frames)
% GROUPLOCSBYFRAME  Sort-then-split localizations into per-frame cell arrays.
%   This is done once before the main loop, replacing the repeated
%   "loc_total(:,3) == f" scan that made the original O(N * F).
%
%   locs_by_frame{f} = [Nf x 2]  subpixel positions for frame f
%   idx_by_frame{f}  = [1  x Nf] row indices into loc_total

locs_by_frame = cell(n_frames, 1);
idx_by_frame  = cell(n_frames, 1);

[sorted_frames, sort_order] = sort(loc_total(:, 3));
n_locs = numel(sorted_frames);
if n_locs == 0, return; end

% Find boundaries where frame number changes
breaks = [1; find(diff(sorted_frames) > 0) + 1];
ends   = [breaks(2:end) - 1; n_locs];

for j = 1:numel(breaks)
    f = sorted_frames(breaks(j));
    if f >= 1 && f <= n_frames
        sel = sort_order(breaks(j):ends(j));
        locs_by_frame{f} = loc_total(sel, 1:2);
        idx_by_frame{f}  = sel(:)';
    end
end

end