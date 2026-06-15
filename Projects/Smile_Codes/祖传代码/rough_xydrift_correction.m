function [corr_data, delta_sum] = rough_xydrift_correction(raw_data, loc_total)

    % This function is intended for a quick xy correction so that trace/intensity analysis would be more accurate. The correction is based on r-radius search between 2 adjacent point clusters. Pixel-level locolization input have been tested, though in principle sub-pixel(SMLM) data shall work.

    % This is a simplified version of Nature Communications | (2025) 16:9031
    % Written by Jerry Ling, 2026.1.9
    
    % ---- inputs ----
    % [raw_data] should be a 3D [rows*cols*frames] matrix. It must be LARGER then corr_data to correct. Note that it should be CENTERED, for example raw_data = data(100:412,100:412,:).
    % [loc_total] should be a 2D [N*3] matrix. N = total number of locolizations, with the 3 columns representing (x_loc, y_loc, frame). 
    
    % ---- parameters ----
    % [corr_rows/corr_cols] decides the size of the output corrected image.
    % [num_clusters] decides the frames in one cluster. If there are blinking then one may set larger number as 15-30. However, note that the drift should be very small compared to cluster time.
    % [r] is the searching radius. r should in principle be larger then 3-sigma to ensure the same emitter be counted in adjacent clusters.

    % ---- outputs ----
    % [corr_data] is the corrected data. It's totally compatible to other codes.
    % [delta_sum] is the size-(n_frames,2) matrix encoding drift over time. 


    [rows, cols, frames] = size(raw_data);
    delta_trace = zeros(frames, 2); 
    corr_rows = 257;
    corr_cols = 257;
    start_row = (rows-corr_rows)/2+1;
    end_row = (rows+corr_rows)/2;
    start_col = (cols-corr_cols)/2+1;
    end_col = (cols+corr_cols)/2;
    corr_data = zeros(corr_rows, corr_cols, frames, 'like', raw_data);
    num_cluster = 10; % frames that count as a time-cluster
    r = 3; % cloud-searching 


    for f=1:num_cluster:(floor(frames/10)-1)*num_cluster

        f_q = f:f+num_cluster-1;
        f_r = f_q + num_cluster; 
        idx_q  = ismember(loc_total(:,3), f_q); 
        idx_r = ismember(loc_total(:,3), f_r); 
        Q = loc_total(idx_q(:) ,1:2);  % query
        R = loc_total(idx_r(:),1:2);  % reference
        
        if isempty(R) || isempty(Q)
            continue
        end

        idx_cell = rangesearch(R, Q, r); % search points in next cluster around the r-circle of points of the first frame. 

        dx_all = []; % collection of x-displacements
        dy_all = []; % collection of y-displacements

        for i = 1:numel(idx_cell)
            idxs = idx_cell{i};
            if isempty(idxs)
                continue
            end

            dxy = R(idxs,:) - Q(i,:);
            dx_all = [dx_all; dxy(:,1)];  % accumulate in rows
            dy_all = [dy_all; dxy(:,2)];
        end

        if isempty(dx_all)
            continue
        end

        dr = hypot(dx_all, dy_all);
        valid = dr < r;
        if ~any(valid)
            continue
        end
        % ----use mean----
        delta_x = mean(dx_all(valid));
        delta_y = mean(dy_all(valid));   
        % ----distribute correction to each frame----
        delta_trace(f_r,1) = delta_x / num_cluster; 
        delta_trace(f_r,2) = delta_y / num_cluster;
        
    end

    % ---- delta_trace time-averaging ----
    %delta_trace_smooth = movmean(delta_trace, 5, 1);
    %delta_sum = cumsum(delta_trace_smooth, 1);

    delta_sum = cumsum(delta_trace, 1);
    
    % ---- correct raw data ----

    for i = 1:frames
        % ---- Subpixel correction with interpolation ----
        img_corr = imtranslate( ...
            raw_data(:,:,i), ...
            -delta_sum(i,[2 1]), ...   
            'linear', ...
            'FillValues', 0);

        % ---- Ensure the 257*257 window----
        corr_data(:,:,i) = img_corr( ...
            start_row:end_row, ...
            start_col:end_col);

        idx = loc_total(:,3) == i;
        if ~any(idx)
            continue
        end

    end

    

