function super_loc_total = locMolecules(img, detected_total, edge, F_EM, use_pixel)

if nargin < 5 || isempty(F_EM)
    F_EM = sqrt(2); % default excess noise factor for EMCCD
end
if nargin < 6 || isempty(use_pixel)
    use_pixel = false; 
end

N = size(detected_total,1);
super_loc_total = nan(N,7);

tic;
for i = 1:N

    frame = detected_total(i,3);
    rough_row = detected_total(i,1);
    rough_col = detected_total(i,2);

    r1 = rough_row - edge;
    r2 = rough_row + edge;
    c1 = rough_col - edge;
    c2 = rough_col + edge;

    if r1<1 || c1<1 || r2>size(img,1) || c2>size(img,2)
        continue
    end

    roi = double(img(r1:r2, c1:c2, frame));

    try
        est = localize.LSQ1.locSingle(roi, F_EM, use_pixel);
    catch
        continue
    end

    r     = est(1);
    c     = est(2);
    Np    = est(3);     % photons
    b     = est(4);     % background
    sigma = est(5);     % PSF width
    slog  = est(6);     % localization precision

    if ~isfinite(Np) || Np<=0 || ...
       ~isfinite(sigma) || sigma<=0 || ...
       ~isfinite(slog) || slog<=0
        continue
    end

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

super_loc_total = super_loc_total(all(isfinite(super_loc_total),2),:);

fprintf('LS fitting complete. %d localizations kept.\n', ...
        size(super_loc_total,1));
end