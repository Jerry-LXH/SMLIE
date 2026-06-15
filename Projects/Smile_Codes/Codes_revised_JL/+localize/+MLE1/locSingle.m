function [locs_refined] = locSingle(roi, F_EM, use_pixel)

% This subpixel localization is beasd on MLE esitmator. Note the model here is a var-adjustable 2D Gaussian with Poisson noise.

% ---- input ---- 
% [roi] should be a part of image. Ideally it sholud only contain one emitter near the center to be estimated. 

% ---- output ----
% [locs_refined] is a 6-number array including: (row,col)(in cordinates relative to the center), Brightness N, background b, width sigma and uncertainty sigma_loc.

sizeROI = size(roi);
ny = (sizeROI(1)-1)/2;
nx = (sizeROI(2)-1)/2;
[yy,xx] = ndgrid(-ny:ny, -nx:nx); % grid matrice
opts = optimoptions('fminunc','Algorithm','quasi-newton','Display','off'); % Quasi-Newton algorithm;
% opts = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'Display','off');

% ---- initial parameter guesses ---- 
b0 = max(median(roi(:)),1);         % background initial guess
N0 = max(sum(max(roi(:)-b0,0)),10); % photons, >=10 to avoid degeneracy
y0 = 0;
x0 = 0;
sigma0 = 1.2;                % initial PSF sigma guess (pixels)
theta0 = [x0; y0; log(N0); log(b0); log(sigma0)];


% ---- negative log-likelihood function ----
fun = @(theta) localize.MLE1.poissonApproxNLL_varSigma(theta, roi, xx, yy, F_EM, use_pixel);
theta_hat = fminunc(fun, theta0, opts);

% ---- extract optimized parameters ----
x0 = theta_hat(1);
y0 = theta_hat(2);
Np = exp(theta_hat(3));
b  = exp(theta_hat(4));
sigma = exp(theta_hat(5));

% ---- Fisher Information for CRLB ----
F = localize.MLE1.fisherInfo_approx_varSigma(theta_hat, xx, yy, F_EM);
C = inv(F); % inverse
sigma_loc = sqrt(0.5*(C(1,1)+C(2,2)));

locs_refined = [y0, x0, Np, b, sigma, sigma_loc];