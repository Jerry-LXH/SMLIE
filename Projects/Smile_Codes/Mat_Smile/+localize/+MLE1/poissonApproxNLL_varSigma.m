function L = poissonApproxNLL_varSigma(theta, roi, xx, yy , F_EM, use_pixel)

% This function gives negative likelihood function -l(experiment|estimatation) as a function of estimated parameters x0, y0, N, b, s. Symmetric Gaussian model with Poisson noise is used.

% ---- inputs ----
% [theta] contains the estimated variables.
% [roi] is the partial image, containing the 'experiment results' in likelihood function.
% [xx] and [yy] is the 2D coordinate grids of roi. 
% [F_EM] compensate for Electronic Magnify noise. If not using EMCCD, set this to 1.

% ---- output ----
% [L] is the negative likelihood function with terms containing variables in theta.

if nargin < 5 || isempty(F_EM)
    F_EM = sqrt(2); % default excess noise factor for EMCCD
end
if nargin < 6 || isempty(use_pixel)
    use_pixel = false; 
end
x0 = theta(1);
y0 = theta(2);
N  = exp(theta(3));
b  = exp(theta(4));
s = exp(theta(5));

if use_pixel == true
    % ----- pixel-integrated Gaussian -----
    Ex = 0.5 * ( ...
        erf((xx - x0 + 0.5)/(sqrt(2)*s)) - ...
        erf((xx - x0 - 0.5)/(sqrt(2)*s)));
    Ey = 0.5 * ( ...
        erf((yy - y0 + 0.5)/(sqrt(2)*s)) - ...
        erf((yy - y0 - 0.5)/(sqrt(2)*s)));
    model = b + N * Ex .* Ey;

elseif use_pixel == false
    % ----- Continuous Gaussian -----
    G = exp(-((xx-x0).^2+(yy-y0).^2)/(2*s^2)) ...
        /(2*pi*s^2);
    model = b + N * G;
else
    error('Undefined use_pixel input.')
end

var_eff = F_EM^2 * (model+1e-12); % EMCCD effective variance

L = 0.5 * sum( (roi(:) - model(:)).^2 ./ var_eff(:) + log(var_eff(:)) ); % For poisson distribution, log P = - n_hat + n*log(n_hat) - Constant
end
