function [emitters_filt, stats] = filterEmitters_jumping(emitters, ecc_th, min_pts)
% FILTEREMITTERS_JUMPING  Remove double-PSF artifacts by detecting
% anisotropic (ellipsoidal) localization clouds.
%
% Physical basis:
%   Single emitter  → isotropic localization noise  → circular scatter
%   Double-PSF fit  → jumps along inter-center axis  → elongated scatter
%
% The eigenvalue ratio of the 2D position covariance matrix quantifies
% this elongation. This metric is scale-invariant, so it does NOT bias
% against dim emitters — only shape matters, not absolute spread.
%
% Inputs:
%   emitters  - struct array with fields .row, .col, .frames, ...
%   ecc_th    - base eigenvalue ratio threshold (default: 3.0)
%               internally adjusted for sample size
%   min_pts   - minimum localizations to evaluate (default: 8)
%
% Outputs:
%   emitters_filt - filtered struct array
%   stats         - (optional) diagnostic struct for threshold tuning

if nargin < 2 || isempty(ecc_th),   ecc_th  = 3.0; end
if nargin < 3 || isempty(min_pts),  min_pts = 8;    end

N    = numel(emitters);
keep = true(1, N);

ecc_all      = nan(1, N);   % eigenvalue ratio per emitter
angle_all    = nan(1, N);   % major axis orientation (deg)
sig_maj_all  = nan(1, N);   % std along major axis
sig_min_all  = nan(1, N);   % std along minor axis
npts_all     = zeros(1, N);

for k = 1:N

    r = emitters(k).row(:);
    c = emitters(k).col(:);
    n = numel(r);
    npts_all(k) = n;

    if n < min_pts
        continue          % too few points — keep by default
    end

    % ---- covariance eigen-decomposition ----
    pos = [r - mean(r), c - mean(c)];
    C   = (pos' * pos) / (n - 1);        % 2x2 covariance
    [V, D] = eig(C);
    lambda  = sort(diag(D), 'descend');   % λ1 >= λ2

    lambda_max = lambda(1);
    lambda_min = max(lambda(2), eps);

    ecc = lambda_max / lambda_min;        % elongation ratio

    % major-axis angle (for diagnostics)
    [~, idx] = max(diag(D));
    angle_all(k)   = atan2d(V(2,idx), V(1,idx));
    ecc_all(k)     = ecc;
    sig_maj_all(k) = sqrt(lambda_max);
    sig_min_all(k) = sqrt(lambda_min);

    % ---- sample-size-adaptive threshold ----
    % For isotropic 2D Gaussian with n samples, the expected
    % eigenvalue ratio is inflated at small n.  Correction term
    % derived from Wishart distribution asymptotics:
    %   E[λ_max/λ_min] ≈ 1 + 2·sqrt(2/n)  for isotropic case
    % We add a margin on top of that.
    ecc_null = 1 + 2 * sqrt(2 / n);       % expected ratio under H0
    ecc_threshold = ecc_th * ecc_null;     % scale user threshold

    if ecc > ecc_threshold
        keep(k) = false;
    end
end

emitters_filt = emitters(keep);

% ---- console summary ----
n_kept = sum(keep);
fprintf('Kept %d / %d emitters (filtered %.1f%% ellipsoidal traces)\n', ...
    n_kept, N, (1 - n_kept/N) * 100);

% ---- optional diagnostics struct ----
if nargout >= 2
    stats.eccentricity   = ecc_all;
    stats.sigma_major    = sig_maj_all;
    stats.sigma_minor    = sig_min_all;
    stats.major_angle    = angle_all;
    stats.n_localizations = npts_all;
    stats.keep           = keep;
    stats.ecc_threshold  = ecc_th;
end

end