function overlayLocs(locs, frames_to_show, margin, varargin)
% This function overlay localizations on current axes.

% ---- inputs ----
% [locs] : Nx3 matrix (row, col, frame)
% [frames_to_show]: frames to show.
% [margin]: when sets true, coordinated will be added 0.5 to convert to pixel-center frame.

% Optional Name-Value pairs:
%   'Color'  : marker color (default 'g')
%   'Marker' : marker style (default 'x')
%   'Size'   : marker size (default 20)
%   'LineWidth' : width for markers
    if nargin < 3 || isempty(margin) 
        margin = false;
    end
    % Extract localizations belonging to these frames
    
    if size(locs,2)==3
        mask = ismember(locs(:,3), frames_to_show);
        locs = locs(mask, :);
    elseif size(locs,2)==2
        fprintf('Whole locs set will be used.\n')
    else
        error('Wrong columns.')
    end

    p = inputParser;
    addParameter(p,'Color','g');
    addParameter(p,'Marker','x');
    addParameter(p,'Size',20);
    addParameter(p,'LineWidth',1.2);
    parse(p,varargin{:});
    opts = p.Results;

    if margin == false
        scatter(locs(:,2), locs(:,1), opts.Size, opts.Color, opts.Marker, ...
            'LineWidth', opts.LineWidth);
    elseif margin == true
        scatter(locs(:,2)+0.5, locs(:,1)+0.5, opts.Size, opts.Color, opts.Marker, ...
            'LineWidth', opts.LineWidth);
    else
        fprint("Error! Undefined input on third position.");
    end

end
