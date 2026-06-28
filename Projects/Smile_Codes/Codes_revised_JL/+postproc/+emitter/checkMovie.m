function checkMovie(index, data, emitters, pos_mean_px, pos_matrix, oneFrameTime, varargin)
% Play local ROI movie of a single emitter. One may adjust ROI, speed, time-window by adding parameters.

p = inputParser;
addParameter(p, 'roiRadius', 7);
addParameter(p, 'frameStep', 3);
addParameter(p, 'preFrames', 100);
addParameter(p, 'postFrames', 300);
addParameter(p, 'clim', [0 250]);
addParameter(p, 'pauseTime', 0.01);  
% 新增可选参数 snr，默认值为空
addParameter(p, 'snr', []);  
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

    % 判断是否传入了 snr 参数
    if ~isempty(opt.snr)
        % 假设 snr 形状与 pos_matrix 相同 [frames x 2 x emitters]
        % 我们取第1列的值作为当前帧的 SNR
        current_snr = opt.snr(index,i); 
        
        if ~isnan(current_snr)
            title_str = sprintf('Emitter %d, t = %.2f s, SNR = %.2f', ...
                index, (i-1)*oneFrameTime, current_snr);
        else
            title_str = sprintf('Emitter %d, t = %.2f s, SNR = NaN', ...
                index, (i-1)*oneFrameTime);
        end
    else
        % 如果没有传入 snr，保持原来的标题格式
        title_str = sprintf('Emitter %d, t = %.2f s', ...
            index, (i-1)*oneFrameTime);
    end
    
    title(ax, title_str);

    drawnow;
    pause(opt.pauseTime);   % <-- 使用控制参数
end

end