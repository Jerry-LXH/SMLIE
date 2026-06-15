function [data, ex_time, gainDAC] = pro_data(file_name, s_frame)

% read data from sif file
file_name = strcat(file_name);
data_set = readsif(file_name);
ex_time = data_set.exposureTime;
gainDAC = data_set.gainDAC;
raw_data = data_set.imageData;

% data conversion
% [~, ~, frames] = size(raw_data);
% data = raw_data(:,:,s_frame:frames/4+s_frame-1) .* 5.75 ./ gainDAC;
data = raw_data(:,:,s_frame:end) .* 5.75 ./ gainDAC;
% data = raw_data(:,:,s_frame:end);

end

