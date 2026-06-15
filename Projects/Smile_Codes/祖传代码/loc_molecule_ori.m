function loc_total = loc_molecule(data, ex_time)
%{
adjustable parameters:
1. threshold
2. filter bright point
%}

% set threshold
[rows, columns, frames] = size(data);
% bg_data = data(:, :, frames-round(10/ex_time)+1:frames);    % value of last second imaging
% bg_down = prctile(bg_data, 5, 'all');
% bg_up = prctile(bg_data, 95, 'all');
% bg = bg_data(bg_data>bg_down & bg_data<bg_up);
% threshold = mean(bg) + 9*std2(bg);
threshold = mean(data(:,:,end), 'all') +7*std2(data(:,:,end)); % adjust threshold by std's coefficient

% locating frame by frame (scan radius: 3 pixels)
loc_total = [];
for frame = 1:frames
    loc = [];
    for row = 4:rows-3
        for column = 4:columns-3
            temp_max = max(data(row-3:row+3,column-3:column+3, frame), [], 'all');
            if data(row, column, frame) == temp_max && data(row, column, frame) >= threshold
                if data(row, column, frame) / min(data(row-1:row+1, column-1:column+1, frame), [], 'all') < 10   % filter bright point
                    loc = [loc; [row column frame]];
                end
            end
        end
    end
    if isempty(loc)     % exclude frames without localizations
        continue
    else
       loc_total = [loc_total; loc]; %#ok<*AGROW>
    end
    clc;
    disp(['processing frame',num2str(frame)]);
%     else
%         loc_fi = loc_filter(loc);   % filter located single molecule
%     end
%     loc_total = [loc_total; loc_fi]; %#ok<*AGROW>
end
clc;
disp('processing complete');
end


