function stats = collectEmitterStatistics(emitters,frames,ex_time,interval, ...
                                          brightness, ...
                                          sigma, ...
                                          sigma_loc,bg)

    N = numel(emitters);
    oneFrameTime = ex_time+interval;
    % --- Position matrix ---
    pos_matrix = nan(frames, 2, N);

    for k = 1:N
        f = emitters(k).frames;
        pos_matrix(f,1,k) = emitters(k).row;
        pos_matrix(f,2,k) = emitters(k).col;
    end

    pos_mean = squeeze(mean(pos_matrix,1,'omitnan')).';
    pos_mean_px = round(pos_mean + 0.5);

    % --- Survival ---
    survival = zeros(N,1);
    for k = 1:N
        t_first = emitters(k).on_frame;
        t_last  = emitters(k).bleach_frame;
        if isnan(t_last)
            t_last = emitters(k).last_frame;
        end
        survival(k) = t_last - t_first + 1;
    end
    survival_sec = survival * oneFrameTime;

    % --- Brightness & sigma ---
    brightness_em = nan(N,frames);
    sigma_em      = nan(N,frames);
    sigma_loc_em  = nan(N,frames);
    bg_em = nan(N,frames);

    for k = 1:N 
        idx = emitters(k).loc_idx;
        f   = emitters(k).frames;
        brightness_em(k,f) = brightness(idx);
        sigma_em(k,f)      = sigma(idx);
        sigma_loc_em(k,f)  = sigma_loc(idx);
        bg_em(k,f) = bg(idx);
    end

    stats.pos_matrix        = pos_matrix;
    stats.pos_mean_px       = pos_mean_px;
    stats.survival_sec      = survival_sec;
    stats.brightness_em = brightness_em;
    stats.bg_em = bg_em;
    stats.bg_mean = mean(bg_em,2,'omitnan');
    stats.brightness_mean   = mean(brightness_em,2,'omitnan')/ex_time; % photon/sec
    stats.brightness_sum    = sum(brightness_em,2,'omitnan');
    stats.sigma_mean        = mean(sigma_em,2,'omitnan');
    stats.sigma_loc_mean    = mean(sigma_loc_em,2,'omitnan');
end