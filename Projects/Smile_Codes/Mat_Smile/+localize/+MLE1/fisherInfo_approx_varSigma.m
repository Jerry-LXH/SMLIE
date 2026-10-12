function F = fisherInfo_approx_varSigma(theta, xx, yy, F_EM)

% This function calculates the Fisher Information of the model at estimated parameters [theta].

% ---- inputs ----
% [theta] contains the estimated variables.
% [xx] and [yy] is the 2D coordinate grids of roi. 
% [F_EM] compensate for Electronic Magnify noise. If not using EMCCD, set this to 1.

% ---- output ----
% [F] is the Fisher Information matrix. Diagonal terms encodes uncertainty of each parameter estimated. 

if nargin < 5 || isempty(F_EM)
    F_EM = sqrt(2); % default EM excess noise
end

x0 = theta(1); y0 = theta(2);
N  = exp(theta(3));
b  = exp(theta(4));
s = exp(theta(5));

% ----- pixel-integrated Gaussian -----
Ex = 0.5 * ( ...
    erf((xx - x0 + 0.5)/(sqrt(2)*s)) - ...
    erf((xx - x0 - 0.5)/(sqrt(2)*s)));

Ey = 0.5 * ( ...
    erf((yy - y0 + 0.5)/(sqrt(2)*s)) - ...
    erf((yy - y0 - 0.5)/(sqrt(2)*s)));

E = Ex .* Ey;
model = b + N * E;

% ----- derivatives -----
% derivatives of integrated Gaussian
gx_p = exp(-((xx - x0 + 0.5).^2)/(2*s^2));
gx_m = exp(-((xx - x0 - 0.5).^2)/(2*s^2));
gy_p = exp(-((yy - y0 + 0.5).^2)/(2*s^2));
gy_m = exp(-((yy - y0 - 0.5).^2)/(2*s^2));

dEx_dx = (gx_m - gx_p) / (sqrt(2*pi)*s);
dEy_dy = (gy_m - gy_p) / (sqrt(2*pi)*s);

d_dx = N * dEx_dx .* Ey;
d_dy = N * Ex .* dEy_dy;

d_dN = E;
d_db = ones(size(model));

% σ derivative (stable approximation)
r2 = (xx-x0).^2 + (yy-y0).^2;
A = N/(2*pi*s^2);
Ecenter = exp(-r2/(2*s^2));
d_ds = A * Ecenter .* (r2/s^3 - 2/s);

% ----- log-parameter chain rule -----
d_dlogN = d_dN * N;
d_dlogb = d_db * b;
d_dlogs = d_ds * s;

J = [d_dx(:), d_dy(:), d_dlogN(:), d_dlogb(:), d_dlogs(:)]; % Jacobian

mu = model(:) + 1e-12;
W = 1 ./ (F_EM^2 * mu) + 1 ./ (2 * mu.^2);
F = J' * (W .* J);

end
