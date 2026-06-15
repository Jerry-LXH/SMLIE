function plotImage(data, frames_to_show, cmap, tit, proj)
% This function show max projection of specified frames.
%
% data   : R x C x F image stack
% frames : subset of frames to project
% cmap   : optional colormap string (default: 'gray')
% titleStr : the title

    % ----- defaults -----
    if nargin < 3 || isempty(cmap)
        cmap = 'gray';
    end
    if nargin < 4 || isempty(tit)
        tit = 'Image';
    end
    if nargin < 5 || isempty(proj)
        proj = 'mean';
    end
    switch lower(proj)
        case 'mean'
            img_proj = mean(data(:,:,frames_to_show), 3);
            titleStr = sprintf('Mean Projection of %s', tit);
        case 'max'
            img_proj = max(data(:,:,frames_to_show),[], 3);
            titleStr = sprintf('Max Projection of %s', tit);
        otherwise
            error('Undefined projection method.')
    end
    imagesc(img_proj);
    axis image;
    colormap(cmap);
    colorbar;
    title(sprintf('%s (%d frames shown)', titleStr, numel(frames_to_show)));
end
