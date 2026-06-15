addpath('mFunc');

lambda = 0.532;
pIn = 3:15;
NA = 1.1;
RI = 1.33;
upsample = 3;
xSize = 128;
xSizeUp = xSize*upsample;
xPixelSize = 0.130; % xPixelSize = 0.130; % um
xPixelSizeUp = xPixelSize/upsample;
zSize = 16;
zSizeUp = zSize*upsample;
zPixelSize = 0.130; 
zPixelSizeUp = zPixelSize/upsample;
aValue = 0.1;
bgValue = 95;
flagBg = 0;% 0.8NA: 1; 1.1NA: 0; 0.71NA: 0
zernIndex = 3:15;
zernType = 'random1';
% conType = 'ANSI';
zernNum = length(zernIndex);
density = 0.0005;
mean = 6000;
std = 2000;
para = 3;

[pointsBlur2D_noabe, pointsBlur2D, coeffs, points, emitterNumber] = gen_blurred_points_img(lambda, pIn, NA, RI, upsample, xSize, xPixelSize, zSize, zPixelSize, aValue, zernIndex, zernType, zernNum, density, mean, std);

% 显示相位掩膜
figure;
imagesc(pointsBlur2D_noabe);
colormap gray; colorbar;
axis equal; axis off;
figure;
imagesc(pointsBlur2D);
colormap gray; colorbar;
axis equal; axis off;
