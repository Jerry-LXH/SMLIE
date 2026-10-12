function filtered_data = filter_hotpixel(raw_data, thr)
% FILTER_HOTPIXEL Remove isolated hot pixels in image stack
%
% raw_data : H x W x T stack (single/double/uint16)
% thr      : brightness factor above local neighborhood (default 1.8)
%
% Method:
%   - compare each pixel with 3x3 neighborhood median
%   - remove only isolated spikes
%   - preserve multi-pixel PSFs

if nargin < 2
    thr = 1.8;   % good starting value for 1.5–3x spikes
end

filtered_data = raw_data;
[~,~,T] = size(raw_data);

for k = 1:T
    frame = raw_data(:,:,k);

    % local median (background estimate)
    medf = medfilt2(frame,[3 3]);

    % spike detection
    spike = frame > thr .* medf;

    % remove spikes
    frame(spike) = medf(spike);

    filtered_data(:,:,k) = frame;
end

end
