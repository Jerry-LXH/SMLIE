function [locs_refined] = locSingle(roi, F_EM, use_pixel,sigma0)

% Least-squares 2D Gaussian fitting 
% Output format identical to MLE version:
% [y0, x0, N, b, sigma, sigma_loc]

if nargin < 4 || isempty(sigma0)
    sigma0 = 1.2; 
end

sizeROI = size(roi);
ny = (sizeROI(1)-1)/2;
nx = (sizeROI(2)-1)/2;
[yy,xx] = ndgrid(-ny:ny, -nx:nx);

roi = double(roi);

opts = optimoptions('lsqnonlin',...
    'Algorithm','levenberg-marquardt',...
    'Display','off');

% ---- initial guesses ----
b0 = max(median(roi(:)),1);
N0 = max(sum(max(roi(:)-b0,0)),10);
y0 = 0;
x0 = 0;

theta0 = [x0; y0; log(N0); log(b0); log(sigma0)];

% ---- LS residual function ----
fun = @(theta) localize.LSQ1.residual(theta, roi, xx, yy, F_EM, use_pixel);

% ---- optimization ----
theta_hat = lsqnonlin(fun, theta0, [], [], opts);

% ---- unpack ----
x0 = theta_hat(1);
y0 = theta_hat(2);
Np = exp(theta_hat(3));
b  = exp(theta_hat(4));
sigma = exp(theta_hat(5));

% ---- covariance estimation (Gauss-Newton approximation) ----
[J,res] = localize.LSQ1.jacobian(theta_hat, roi, xx, yy, F_EM, use_pixel);

sigma_r2 = sum(res.^2) / (numel(res) - numel(theta_hat));

C = inv(J.'*J) * sigma_r2;

sigma_loc = sqrt(0.5*(C(1,1)+C(2,2)));

locs_refined = [y0, x0, Np, b, sigma, sigma_loc];

end