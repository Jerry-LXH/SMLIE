function emitters_filt = filterEmitters_short( ...
    emitters, min_on_frames, mode)

% mode:
%   'consecutive' (default) - 连续ON帧数过滤
%   'total'                - 总ON帧数过滤

if nargin < 3
    mode = 'consecutive';
end

keep = false(1, numel(emitters));

for k = 1:numel(emitters)

    fr = emitters(k).frames(:)';
    fr = sort(fr);

    if isempty(fr)
        continue
    end

    switch lower(mode)

        % ===== 连续 ON 帧 =====
        case 'consecutive'

            d = diff(fr);

            run_len = 1;
            max_run = 1;

            for i = 1:numel(d)
                if d(i) == 1
                    run_len = run_len + 1;
                else
                    max_run = max(max_run, run_len);
                    run_len = 1;
                end
            end

            max_run = max(max_run, run_len);

            if max_run >= min_on_frames
                keep(k) = true;
            end

        % ===== 总 ON 帧 =====
        case 'total'

            if numel(fr) >= min_on_frames
                keep(k) = true;
            end

        otherwise
            error('Unknown mode. Use "consecutive" or "total".');
    end
end

emitters_filt = emitters(keep);

num = numel(emitters);
num_filt = numel(emitters_filt);
filtered = (num-num_filt)/num*100;

fprintf(['Totally %d emitters left. %.2f%% emitters without ' ...
         '%d ON frames (%s mode) are filtered.\n'], ...
         num_filt, filtered, min_on_frames, mode);

end
