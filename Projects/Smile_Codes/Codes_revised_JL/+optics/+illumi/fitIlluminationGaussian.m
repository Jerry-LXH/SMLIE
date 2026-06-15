function params = fitIlluminationGaussian(L_est)

L_est = double(L_est);

[ny,nx] = size(L_est);
[xx,yy] = meshgrid(1:nx,1:ny);

% 初始值
A0 = max(L_est(:)) - min(L_est(:));
B0 = min(L_est(:));
x0 = nx/2;
y0 = ny/2;
sx0 = nx/4;
sy0 = ny/4;

theta0 = [A0, x0, y0, sx0, sy0, B0];

opts = optimoptions('lsqnonlin','Display','off');

theta = lsqnonlin(@(p) optics.illumi.residual(p,L_est,xx,yy), theta0,[],[],opts);

params.A  = theta(1);
params.x0 = theta(2);
params.y0 = theta(3);
params.sx = theta(4);
params.sy = theta(5);
params.B  = theta(6);

end