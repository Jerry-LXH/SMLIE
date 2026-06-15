function checkMovie(index, data, emitters, pos_mean_px, pos_matrix, oneFrameTime, varargin)
% Play local ROI movie of a single emitter. One may adjust ROI, speed, time-window by adding parameters.

p = inputParser;
addParameter(p, 'roiRadius', 7);
addParameter(p, 'frameStep', 3);
addParameter(p, 'preFrames', 100);
addParameter(p, 'postFrames', 300);
addParameter(p, 'clim', [0 250]);
addParameter(p, 'pauseTime', 0.01);  
parse(p, varargin{:});
opt = p.Results;

[H, W, totalFrames] = size(data);

firstfr = max(emitters(index).on_frame - opt.preFrames, 1);
endfr   = min(emitters(index).last_frame + opt.postFrames, totalFrames);

roi_r = opt.roiRadius;
r = pos_mean_px(index,1);
c = pos_mean_px(index,2);

r = max(roi_r+1, min(H-roi_r, r));
c = max(roi_r+1, min(W-roi_r, c));

r0 = r - roi_r;
c0 = c - roi_r;

local_field = data(r-roi_r:r+roi_r, ...
                    c-roi_r:c+roi_r, ...
                    1:endfr);

fig = figure( ...
    'Name', sprintf('Emitter %d local movie', index), ...
    'NumberTitle','off', ...
    'Position',[900 200 500 520]);

ax = axes('Parent', fig);

for i = firstfr:opt.frameStep:endfr

    imagesc(ax, local_field(:,:,i));
    axis(ax,'image');
    colormap(ax, gray);
    colorbar(ax);
    clim(ax, opt.clim);
    set(ax,'YDir','normal','FontSize',12);

    hold(ax,'on');

    if ~isnan(pos_matrix(i,1,index))
        rr = pos_matrix(i,1,index) + 0.5;
        cc = pos_matrix(i,2,index) + 0.5;

        scatter(ax, ...
            cc - c0 + 1, ...
            rr - r0 + 1, ...
            40, 'g', 'x');
    end

    hold(ax,'off');

    title(ax, sprintf('Emitter %d, t = %.2f s', ...
        index, (i-1)*oneFrameTime));

    drawnow;
    pause(opt.pauseTime);   % <-- 使用控制参数
end

end