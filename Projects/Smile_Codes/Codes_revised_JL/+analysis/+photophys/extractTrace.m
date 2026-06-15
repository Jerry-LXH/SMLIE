function [series, snr_series, bg_series, raw_series] = extractTrace(loc, data, edge, doBgSubtract, bgEdge, bgMethod, exclRadius)
% extractTrace  Integrated intensity trace with robust local background.
%
% ---- inputs ----
% [loc]           ALL emitter coordinates, [N x 2], [row, col]
% [data]          image stack, [rows x cols x frames]
% [edge]          half-size of signal integration square (default 3)
% [doBgSubtract]  true/false (default false)
% [bgEdge]        half-size of outer background square (default edge+3)
% [bgMethod]      'median' | 'trimmean' | 'mean' (default 'median')
% [exclRadius]    exclusion radius around each neighbor (default = edge)
%                 Pixels within exclRadius of any OTHER emitter are excluded
%                 from the background annulus.
%
% ---- outputs ----
% [series]        background-corrected trace (or raw if doBgSubtract=false)
% [snr_series]    spatial SNR per frame
% [bg_series]     estimated background per pixel per frame
% [raw_series]    raw integrated trace

% ======== Defaults ========
if nargin < 3 || isempty(edge),         edge = 3;            end
if nargin < 4 || isempty(doBgSubtract), doBgSubtract = false; end
if nargin < 5 || isempty(bgEdge),       bgEdge = edge + 3;   end
if nargin < 6 || isempty(bgMethod),     bgMethod = 'median';  end
if nargin < 7 || isempty(exclRadius),   exclRadius = edge;    end

if bgEdge <= edge
    error('bgEdge must be larger than edge.');
end

[loc_num, ~] = size(loc);
[rows, cols, frames] = size(data);

intSize = 2*edge + 1;
intArea = intSize^2;

% ========== Info output ==========
fprintf('========================================\n');
fprintf('  extractTrace: Configuration Summary\n');
fprintf('========================================\n');
fprintf('  Emitters:          %d\n', loc_num);
fprintf('  Image size:        %d x %d x %d frames\n', rows, cols, frames);
fprintf('  Signal ROI:        %d x %d  (edge = %d)\n', intSize, intSize, edge);
if doBgSubtract
    bgOutSize = 2*bgEdge + 1;
    bgPixelCount = bgOutSize^2 - intArea;
    fprintf('  Background:        ON\n');
    fprintf('  Background ROI:    %d x %d  (bgEdge = %d)\n', bgOutSize, bgOutSize, bgEdge);
    fprintf('  Background pixels: ~%d per frame (before neighbor exclusion)\n', bgPixelCount);
    fprintf('  Background method: ''%s''\n', bgMethod);
    fprintf('  Neighbor masking:  ON  (exclRadius = %d)\n', exclRadius);
    fprintf('  Sigma-clipping:    3-sigma, 3 iterations\n');
else
    fprintf('  Background:        OFF  (raw integration only)\n');
end
fprintf('========================================\n\n');

% ========== Preallocate ==========
raw_series = nan(loc_num, frames);
series     = nan(loc_num, frames);
snr_series = nan(loc_num, frames);
bg_series  = nan(loc_num, frames);

% ========== Precompute neighbor exclusion masks (frame-independent) ==========
neighborMasks = cell(loc_num, 1);

for num = 1:loc_num
    r = round(loc(num,1));
    c = round(loc(num,2));

    % Outer background patch boundaries (clipped to image)
    r_out1 = max(1, r - bgEdge);
    r_out2 = min(rows, r + bgEdge);
    c_out1 = max(1, c - bgEdge);
    c_out2 = min(cols, c + bgEdge);

    nR = r_out2 - r_out1 + 1;
    nC = c_out2 - c_out1 + 1;

    % Start with all true (all pixels available)
    validMask = true(nR, nC);

    % Exclude signal box of the current emitter
    in_r1 = (r - edge) - r_out1 + 1;
    in_r2 = (r + edge) - r_out1 + 1;
    in_c1 = (c - edge) - c_out1 + 1;
    in_c2 = (c + edge) - c_out1 + 1;
    validMask(in_r1:in_r2, in_c1:in_c2) = false;

    % Exclude regions around ALL OTHER emitters
    for m = 1:loc_num
        if m == num
            continue;
        end

        rn = round(loc(m, 1));
        cn = round(loc(m, 2));

        % Quick check: is this neighbor anywhere near our bg patch?
        if rn < r_out1 - exclRadius || rn > r_out2 + exclRadius || ...
           cn < c_out1 - exclRadius || cn > c_out2 + exclRadius
            continue;
        end

        % Exclusion box of neighbor in patch coordinates
        excl_r1 = max(1, (rn - exclRadius) - r_out1 + 1);
        excl_r2 = min(nR, (rn + exclRadius) - r_out1 + 1);
        excl_c1 = max(1, (cn - exclRadius) - c_out1 + 1);
        excl_c2 = min(nC, (cn + exclRadius) - c_out1 + 1);

        validMask(excl_r1:excl_r2, excl_c1:excl_c2) = false;
    end

    neighborMasks{num} = validMask;
end

% ========== Main extraction loop ==========
for num = 1:loc_num
    r = round(loc(num,1));
    c = round(loc(num,2));

    % Signal region indices
    r_in = (r-edge):(r+edge);
    c_in = (c-edge):(c+edge);

    % Background patch boundaries
    r_out1 = max(1, r - bgEdge);
    r_out2 = min(rows, r + bgEdge);
    c_out1 = max(1, c - bgEdge);
    c_out2 = min(cols, c + bgEdge);

    validMask = neighborMasks{num};
    nValidPixels = sum(validMask, 'all');

    for frame = 1:frames
        img = data(:,:,frame);

        % Raw signal integration
        patch_in = img(r_in, c_in);
        raw_val = sum(double(patch_in), 'all');
        raw_series(num, frame) = raw_val;

        if doBgSubtract
            % Extract background patch, apply precomputed mask
            patch_out = double(img(r_out1:r_out2, c_out1:c_out2));
            bg_pixels = patch_out(validMask);
            bg_pixels = bg_pixels(isfinite(bg_pixels));

            % Iterative sigma-clipping (3 iterations, 3-sigma)
            bg_pixels = sigmaClip(bg_pixels, 3.0, 3);

            if numel(bg_pixels) < 5
                % Fallback: use unclipped values with median
                bg_pixels_raw = patch_out(validMask);
                bg_pixels_raw = bg_pixels_raw(isfinite(bg_pixels_raw));
                if isempty(bg_pixels_raw)
                    bg_per_pixel = 0;
                    bg_std = 0;
                else
                    bg_per_pixel = median(bg_pixels_raw);
                    bg_std = std(bg_pixels_raw);
                end
            else
                bg_std = std(bg_pixels);
                switch lower(bgMethod)
                    case 'median'
                        bg_per_pixel = median(bg_pixels);
                    case 'trimmean'
                        if numel(bg_pixels) >= 5
                            bg_per_pixel = trimmean(bg_pixels, 20);
                        else
                            bg_per_pixel = median(bg_pixels);
                        end
                    case 'mean'
                        bg_per_pixel = mean(bg_pixels);
                    otherwise
                        error('Unknown bgMethod: %s', bgMethod);
                end
            end

            bg_val = bg_per_pixel * intArea;
            total_signal = raw_val - bg_val;

            % Spatial SNR (shot noise + background noise model)
            if total_signal > 0 && bg_std > 0
                snr_val = total_signal / sqrt(total_signal + intArea * (bg_std^2));
            else
                snr_val = 0;
            end

            series(num, frame)     = total_signal;
            snr_series(num, frame) = snr_val;
            bg_series(num, frame)  = bg_per_pixel;
        else
            series(num, frame)     = raw_val;
            snr_series(num, frame) = 0;
            bg_series(num, frame)  = 0;
        end
    end
end

fprintf('Done. %d traces x %d frames extracted.\n', loc_num, frames);

% ========== Report masking statistics ==========
validCounts = cellfun(@(m) sum(m, 'all'), neighborMasks);
fprintf('  Background pixels per emitter: min=%d, median=%d, max=%d\n', ...
    min(validCounts), round(median(validCounts)), max(validCounts));
end


%% ========== Helper: Iterative Sigma-Clipping ==========
function data_clean = sigmaClip(data, nSigma, nIter)
% Remove outliers iteratively: reject points > nSigma * std from median.
    data_clean = data(:);
    for iter = 1:nIter
        if numel(data_clean) < 3
            break;
        end
        med = median(data_clean);
        s = std(data_clean);
        if s < 1e-10
            break;
        end
        keep = abs(data_clean - med) <= nSigma * s;
        if all(keep)
            break;
        end
        data_clean = data_clean(keep);
    end
end