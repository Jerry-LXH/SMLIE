function r = residual(theta, roi, xx, yy, F_EM, use_pixel)

% Weighted least-squares residual
% Compatible with MLE noise model (Poisson + EMCCD option)

if nargin < 5 || isempty(F_EM)
    F_EM = sqrt(2);      % default excess noise factor for EMCCD
end
if nargin < 6 || isempty(use_pixel)
    use_pixel = false;  
end

x0 = theta(1);
y0 = theta(2);
N  = exp(theta(3));
b  = exp(theta(4));
s  = exp(theta(5));

% ----- Gaussian model -----
if use_pixel
    Ex = 0.5 * ( ...
        erf((xx - x0 + 0.5)/(sqrt(2)*s)) - ...
        erf((xx - x0 - 0.5)/(sqrt(2)*s)));

    Ey = 0.5 * ( ...
        erf((yy - y0 + 0.5)/(sqrt(2)*s)) - ...
        erf((yy - y0 - 0.5)/(sqrt(2)*s)));

    model = b + N * Ex .* Ey;
else
    G = exp(-((xx-x0).^2+(yy-y0).^2)/(2*s^2)) ...
        /(2*pi*s^2);
    model = b + N * G;
end

% ----- Effective variance (Poisson + EMCCD) -----
var_eff = F_EM^2 * (model + 1e-12);

% ----- Weighted residual -----
r = (model - roi) ./ sqrt(var_eff);
r = r(:);

end