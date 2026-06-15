function plotTracking(delta, ex_time, varargin)

% VISUALIZE_DRIFT Visualize XY drift trajectory colored by time. Note that the unit is in pixels.
%
%   visualize_drift(delta, ex_time)
%   visualize_drift(delta, ex_time, figHandle)
%   visualize_drift(delta, ex_time, figHandle, titleStr)
%
% INPUTS:
%   delta     : [F x 2] matrix, delta(:,1)=dy, delta(:,2)=dx (pixel)
%   ex_time   : exposure time per frame (seconds)
%   figHandle : This one is normally [] so that tracking is plotted in a new figure.
%   titleStr  : Title string (default: 'tracking')


    % ---------------- parse inputs ----------------
    figHandle = [];
    titleStr  = 'Tracking';

    if numel(varargin) >= 1
        figHandle = varargin{1};
    end
    if numel(varargin) >= 2
        titleStr = varargin{2};
    end

    % ---------------- figure ----------------
    if ~isempty(figHandle) && ishandle(figHandle)
        figure(figHandle);
    else
        figure;
    end

    % ---------------- data ----------------
    y = delta(:,1);   % Δy
    x = delta(:,2);   % Δx
    F = size(delta,1);
    t = (0:F-1) * ex_time;

    % ---------------- plot ----------------
    hold on
    patch( ...
        'XData', x, ...
        'YData', y, ...
        'ZData', zeros(F,1), ...
        'CData', t(:), ...
        'EdgeColor','interp', ...
        'LineWidth',2, ...
        'FaceColor','none' ...
    );

    axis equal
    axis auto
    view(2)

    colormap(flipud(jet))
    cb = colorbar;
    cb.Label.String = 'Time (s)';
    caxis([min(t) max(t)])

    xlabel('x/col (pixel)')
    ylabel('y/row (pixel)')
    title(titleStr)

    set(gca, 'FontSize', 12)
end
