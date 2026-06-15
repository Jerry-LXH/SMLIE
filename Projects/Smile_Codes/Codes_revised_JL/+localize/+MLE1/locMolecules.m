function super_loc_total = locMolecules(img, detected_total, edge, F_EM, use_pixel)

% This function locates the subpixel coordinates based on Maxium Likelihood Estimation in the ROI near the rough detections. 

% ---- inputs ----
% [img] is the whole image stack shaped [row*col*frames].
% [detected_total] is the rough detection of local maximums shaped [N*3], the 3 cols reading: row coordinates, col coordinates and frame number.

% ---- outputs ----
% The 7 outputs columns read: row_loc, col_loc, Brightness N, background b, width sigma, uncertainty sigma_loc and frame number f. 

% Note that detected_total must not appear at the edge of the image! (Or the roi will exceed the boundary)

if nargin < 5 || isempty(F_EM)
    F_EM = sqrt(2); % default excess noise factor for EMCCD
end
if nargin < 6 || isempty(use_pixel)
    use_pixel = false; 
end

N = size(detected_total,1);
super_loc_total = nan(N, 7);
tic;
% profile on
for i = 1:size(detected_total,1)
    frame = detected_total(i,3);
    rough_row = detected_total(i,1);
    rough_col = detected_total(i,2);
    r1 = rough_row - edge;
    r2 = rough_row + edge;
    c1 = rough_col - edge;
    c2 = rough_col + edge;
    if r1 < 1 || c1 < 1 || r2 > size(img,1) || c2 > size(img,2)
        fprintf('ROI out of range !');
        continue
    end
    roi = double(img(r1:r2, c1:c2, frame));

    % ---- MLE fitting ----
    try
        est = localize.MLE1.locSingle(roi, F_EM, use_pixel);   % returns [x, y, N, b, sigma, sigma_loc]
    catch
        fprintf("Wrong! \n");
        continue
    end

    r     = est(1);
    c     = est(2);
    Np    = est(3);     % photons
    b     = est(4);     % background
    sigma = est(5);     % PSF width
    slog  = est(6);     % localization precision

    % Filter conditions (NaN or nonphysical)
    if ~isfinite(Np)    || Np <= 0     || ...
       ~isfinite(sigma) || sigma <= 0 || ...
       ~isfinite(slog)  || slog <= 0
        continue
    end

    % ---- Convert to absolute super-res coordinates ----
    row = r + rough_row - 0.5;
    col = c + rough_col - 0.5;
    super_loc_total(i,:) = [row, col, Np, b, sigma, slog, frame];

    if mod(i,2000) == 0 || i == N
        fprintf('Progress: %.1f%%\n', (i/N)*100);
    end
end
t_total = toc;
t_per = t_total / N;
fprintf('Total time = %.3f s, Per-loc = %.3f ms\n', t_total, t_per*1000);

% ---- Remove failed rows ----
super_loc_total = super_loc_total(all(isfinite(super_loc_total),2),:);

fprintf('MLE estimation complete. %d localizations kept.\n', ...
        size(super_loc_total,1));

% profile viewer