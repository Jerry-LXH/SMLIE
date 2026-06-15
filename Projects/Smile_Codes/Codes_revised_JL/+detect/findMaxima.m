function detected_total = findMaxima(data, k_sigma, edge, F_EM)

% This function find local maximum that satisfy some given conditions.

% ---- inputs ----
% [data] is the image stack shaped [row*col*frames].
% [k_sigma] sets the lower bound of Signal-to-Noise Ratio. Valid maximum value must exceed k_sigma * local_sigma. Typically set as 4-6.
% [edge] will mask out those edge maxmium. Set carefully so that the integration later in pipeline won't exceed boundary.
% [F_EM] compensate for Electronic Magnify noise. If not using EMCCD, set this to 1.

% ---- outputs ----
% [detected_local] is the pixel-level detection of local maximums shaped [N*3], the 3 cols reading: row coordinates, col coordinates and frame number.

if nargin < 2 || isempty(k_sigma)
    k_sigma = 4;  % detection sensitivity (3-5 is reasonable)
end
if nargin < 3 || isempty(edge)
    edge = 5;
end
if nargin < 4 || isempty(F_EM)
    F_EM = sqrt(2); 
end

[rows, cols, frames] = size(data);

% ---- structuring elements ----
se7 = true(7);      % for local max detection

detected_total = zeros(1e6, 3);
cnt = 0;
win = 15;
tic;
for frame = 1:frames
    img = data(:,:,frame);

    % --- fast background estimation ---
    bg = imboxfilt(img, win);

    % --- fast noise estimate (shot noise dominated) ---
    sigma = F_EM * sqrt(abs(bg));
    T = bg + k_sigma * sigma;

    % --- local maxima ---
    local_max = img == imdilate(img, se7);

    % ---- edge mask ----
    edge_mask = false(rows, cols);
    edge_mask(1+edge:rows-edge, 1+edge:cols-edge) = true;

    % ---- neighbor mask ----
    kernel = ones(3);
    above_T = img >= T;
    neighbor_count = conv2(double(above_T), kernel, 'same');

    % ---- final mask ----
    mask = local_max & img >= T & edge_mask & (neighbor_count >= 2);


    [r, c] = find(mask);
    n = numel(r);

    if n > 0
        detected_total(cnt+1:cnt+n, :) = [r, c, frame*ones(n,1)];
        cnt = cnt + n;
    end

    if mod(frame,50)==0
        fprintf('processing frame %d / %d\n', frame, frames);
    end
end
t_total = toc;
t_per = t_total / frames;
fprintf('Total time = %.3f s, Per-frame = %.3f ms\n', t_total, t_per*1000);
detected_total = detected_total(1:cnt,:);
fprintf('Detection complete, %.2f detects/frame.\n',size(detected_total,1)/frames);
end