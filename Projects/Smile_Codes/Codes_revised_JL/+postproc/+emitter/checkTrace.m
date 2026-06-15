function checkTrace(start_frame, end_frame, oneFrameTime, trace, brightness_em, row_trace, col_trace)
%CHECKTRACE Visualize trace, fitted brightness, blinking pattern, and track
%
% Inputs:
%   start_frame    - global start frame of this window
%   end_frame      - global end frame of this window
%   oneFrameTime   - exposure time per frame, in seconds
%   trace          - raw intensity trace, column or row vector, length = window length
%   brightness_em  - fitted brightness trace, column or row vector, length = window length
%                    NaN indicates undetected/off frames
%   row_trace      - row position trace, column or row vector, length = window length
%                    use NaN for undetected frames if needed
%   col_trace      - col position trace, column or row vector, length = window length
%                    use NaN for undetected frames if needed
%
% Notes:
%   1. All input traces must have the same length.
%   2. Blinking pattern is inferred directly from brightness_em:
%         ON  -> ~isnan(brightness_em)
%         OFF -> isnan(brightness_em)

    % ---- reshape to column vectors ----
    trace = trace(:);
    brightness_em = brightness_em(:);
    row_trace = row_trace(:);
    col_trace = col_trace(:);

    % ---- basic checks ----
    if ~isscalar(start_frame) || ~isscalar(end_frame) || ~isscalar(oneFrameTime)
        error('start_frame, end_frame, and oneFrameTime must be scalars.');
    end

    if start_frame < 1 || end_frame < start_frame
        error('Invalid frame range: require 1 <= start_frame <= end_frame.');
    end

    win_len = end_frame - start_frame + 1;

    if length(trace) ~= win_len
        error('Length of trace (%d) does not match window length (%d).', ...
            length(trace), win_len);
    end

    if length(brightness_em) ~= win_len
        error('Length of brightness_em (%d) does not match window length (%d).', ...
            length(brightness_em), win_len);
    end

    if length(row_trace) ~= win_len
        error('Length of row_trace (%d) does not match window length (%d).', ...
            length(row_trace), win_len);
    end

    if length(col_trace) ~= win_len
        error('Length of col_trace (%d) does not match window length (%d).', ...
            length(col_trace), win_len);
    end

    % ---- time / frame axis ----
    frame_axis = (start_frame:end_frame).';
    time_axis = (frame_axis - 1) * oneFrameTime;   % seconds

    % ---- blinking pattern inferred from brightness_em ----
    is_on = ~isnan(brightness_em);

    % ---- track validity ----
    has_pos = ~isnan(row_trace) & ~isnan(col_trace);

    % ---- plotting ----
    figure

    % 1) Raw trace
    subplot(3,1,1)
    plot(time_axis, trace, 'k-', 'LineWidth', 1.2)
    hold on
    plot([time_axis(1), time_axis(end)], [0 0], '--', ...
        'Color', [0 0 0], 'LineWidth', 1)
    hold off
    xlim([time_axis(1), time_axis(end)])
    xlabel('Time (s)')
    ylabel('Intensity')
    title(sprintf('Raw trace (frames %d-%d)', start_frame, end_frame))
    set(gca, 'FontSize', 12)

   % 2) Fitted brightness
    subplot(3,1,2)
    brightness_plot = brightness_em;
    brightness_plot(isnan(brightness_plot)) = 0;   % connect to baseline

    plot(time_axis, brightness_plot, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5)
    hold on
    plot([time_axis(1), time_axis(end)], [0 0], '--', ...
        'Color', [0 0.4470 0.7410], 'LineWidth', 1)
    hold off

    xlim([time_axis(1), time_axis(end)])
    xlabel('Time (s)')
    ylabel('Brightness')
    title('Fitted brightness')
    set(gca, 'FontSize', 12)

    % 3) Blinking pattern
    subplot(3,1,3)
    stem(time_axis(is_on), ones(sum(is_on),1), 'filled', 'MarkerSize', 4)
    xlim([time_axis(1), time_axis(end)])
    ylim([0, 1.2])
    xlabel('Time (s)')
    ylabel('Emitter ON')
    title('Blinking pattern')
    set(gca, 'FontSize', 12)

    % 4) Track (row-col trajectory)
    figure
    plot(col_trace(has_pos), row_trace(has_pos), '-o', 'LineWidth', 1.2, 'MarkerSize', 4)
    set(gca, 'YDir', 'reverse')   % image coordinates: row increases downward
    axis equal
    grid on
    xlabel('Column (pixel)')
    ylabel('Row (pixel)')
    title('Emitter track')
    set(gca, 'FontSize', 12)

    % Optional: mark start/end points if there are valid positions
    if any(has_pos)
        hold on
        valid_idx = find(has_pos);
        plot(col_trace(valid_idx(1)), row_trace(valid_idx(1)), 'go', ...
            'MarkerFaceColor', 'g', 'MarkerSize', 7)
        plot(col_trace(valid_idx(end)), row_trace(valid_idx(end)), 'ro', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 7)
        legend({'Track', 'Start', 'End'}, 'Location', 'best')
        hold off
    end
end