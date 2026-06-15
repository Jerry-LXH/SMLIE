function emitters = findEmitters_copy(data, loc_total, max_gap, r_track)

% This function is used to extract emitters from localizations series. Based on time-gap and space-distance, new localizations will be either put into an existing emitter or marked as a new one. This allows for new-coming emitters to be recorded. The complexity approaches O(Loc_number * local_density)

% ---- inputs ----
% [data] is the image stack. It's used only for extracting total frames.
% [loc_total] is the localizations series. It should be a 2D [N*3] matrix carrying subpixel-level locs. N = total number of locolizations, with the 3 columns representing (row_loc, col_loc, frames). 
% [max_gap] sets the time-threshold beyond which an undetected emitter will be marked as bleached. Note it's in frame-unit. If duty cycle is small, one may set a very large number.
% [r_track] sets the space-threshold within which localizations can be clustered into one emitter. For drift-corrected fixed emitters, this value should typically be 3-4 sigma.

% ---- outputs ----
% [emitters] is a strucure array. emitter{k} contains infomation of i-th detected emitter. For each emitter, [.row]/[.col]/[.frames]/[.idx] are arrays of length N_on, N_on is the times this emitter is detected. Specially, [.idx] records the row-index of initial loc_total. [.on_frame] and [.bleach_frame] records when the emitter is first detected and bleached. For those unbleached till the end, [.bleach_frame] is recorded as [.last_frame]. [.last_frame] is the last time an emitter is detected. It's mainly a functional variable. [.alive] records the state of an emitter. Those living till end will have "true" value. One may use this to filter out those emitters.
% emitters(k):
% row: [double array]
% col: [double array]
% frames: [int array]
% loc_idx: [int array]
% on_frame: int
% bleach_frame: int
% last_frame: int
% alive: logical

if nargin < 3 || isempty(max_gap)
    max_gap = 30;
end
if nargin < 4 || isempty(r_track)
    r_track = 1;   % in same units as loc_total
end

[~,~,frames] = size(data);
emitters = struct([]);
tic;
for f = 1:frames
    idx = loc_total(:,3) == f;
    locs = loc_total(idx,1:2);
    loc_idx_f = find(idx);   % indices of locs in loc_total
    if isempty(locs)
        continue
    end

    num_emitters = numel(emitters);
    matched = false(1, num_emitters);

    % ---- loop through localizations in this frame ----
    for i = 1:size(locs,1)
        pos = locs(i,:);
        idx_global = loc_idx_f(i);   % global index in loc_total
        best_k = 0;
        best_d = inf;

        % search through currently alive emitters
        for k = 1:num_emitters
            if ~emitters(k).alive || matched(k)
                continue
            end
            if f - emitters(k).last_frame > max_gap
                continue
            end

            prev_pos = [emitters(k).row(end), emitters(k).col(end)];
            d = hypot(pos(1)-prev_pos(1), pos(2)-prev_pos(2));

            if d < r_track && d < best_d
                best_d = d;
                best_k = k;
            end
        end

        if best_k > 0 
            % append to existing emitter
            emitters(best_k).row(end+1) = pos(1);
            emitters(best_k).col(end+1) = pos(2);
            emitters(best_k).frames(end+1) = f;
            emitters(best_k).last_frame = f;
            matched(best_k) = true;
            emitters(best_k).loc_idx(end+1) = idx_global;

        else
            % create new emitter
            k = num_emitters + 1;
            emitters(k).row = pos(1);
            emitters(k).col = pos(2);
            emitters(k).frames = f;
            emitters(k).on_frame = f;
            emitters(k).last_frame = f;
            emitters(k).bleach_frame = NaN;
            emitters(k).alive = true;
            emitters(k).loc_idx = idx_global;   % store index

            matched(k) = true;
            num_emitters = k;
        end


    end

    % ---- close dead emitters ----
    for k = 1:num_emitters
        if emitters(k).alive && f - emitters(k).last_frame > max_gap
            emitters(k).alive = false;
            emitters(k).bleach_frame = emitters(k).last_frame;
        end
    end

    if mod(f,10) == 0 || f == frames
        fprintf('Progress: %.1f%%\n', (f/frames)*100);
    end
end

% ---- Deal with those still alive by the end of detection. Their bleach_frame are recorded as the last_frame. The alive tag may be used later to filter these out for photophysical analysis. Meanwhile add some extra information.
N = numel(emitters);
for k = 1:numel(emitters)
    if emitters(k).alive
        emitters(k).bleach_frame = emitters(k).last_frame; % or = frames
    end
end

t_total = toc;
fprintf('Total time = %.3f s\n', t_total);
fprintf('Totally %d emitters found.\n', N);