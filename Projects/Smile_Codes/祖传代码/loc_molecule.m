function loc_total = loc_molecule(data, ex_time)

[rows, cols, frames] = size(data);

% ---- threshold estimation ----
threshold = mean(data(:,:,end), 'all') + 7*std2(data(:,:,end));
contrast_th = 10;

se7 = true(7);   % for local maxima
se3 = true(3);   % for local contrast

loc_total = zeros(1e6, 3);
cnt = 0;

for frame = 1:frames
    
    img = data(:,:,frame);

    % --- local maximum (7x7) ---
    local_max = img == imdilate(img, se7);

    % --- local minimum (3x3), avoid zero division ---
    local_min = imerode(img, se3);
    local_min(local_min <= 0) = eps;

    % --- contrast ---
    contrast = img ./ local_min;

    % --- edge mask ---
    edge_mask = false(rows, cols);
    edge_mask(4:rows-3, 4:cols-3) = true;

    % --- final mask ---
    mask = local_max & ...
           (img >= threshold) & ...
           (contrast < contrast_th) & ...
           edge_mask;

    [r, c] = find(mask);
    n = numel(r);

    if n > 0
        loc_total(cnt+1:cnt+n, :) = [r, c, frame*ones(n,1)];
        cnt = cnt + n;
    end

    if mod(frame,50)==0
        fprintf('processing frame %d / %d\n', frame, frames);
    end
end

loc_total = loc_total(1:cnt, :);

disp('processing complete');
end
