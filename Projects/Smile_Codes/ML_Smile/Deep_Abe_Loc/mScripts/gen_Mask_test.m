
addpath('mFunc');

lambda = 0.532;
pIn = 3:15;
NA = 1.1;
RI = 1.33;
upsample = 6;
xSize = 128;
xSizeUp = xSize*upsample;
xPixelSize = 0.130; % xPixelSize = 0.130; % um
xPixelSizeUp = xPixelSize/upsample;
zSize = 48;
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
density = 0.00005;
mean = 6000;
std = 2000;
para = 3;

% eng_wavefront = generate_DH_PSF(r, theta, idx, [1,2,3,4,5], [1,3,5,7,9], 0.15, 30);
% eng_wavefront = generate_DH_PSF(r, theta, idx, [5], [9], 0.15, 20);
mod = imread([fileFolderOut, fileNameBase,'Mask_mod', '.tif']);
ang = imread([fileFolderOut, fileNameBase,'Mask_angle', '.tif']);
eng_wavefront = exp(1i.*ang);
% disp(eng_wavefront);

figure;
imagesc(mod);
colormap hsv; colorbar;
axis equal; axis off;

% 显示相位掩膜
figure;
imagesc(angle(eng_wavefront));
colormap hsv; colorbar;
axis equal; axis off;

% WriteTifStack(eng_wavefront,[fileFolderOut, fileNameBase, '.tif'], 32);
% WriteTifStack(PSF3d,[fileFolderOut, fileNameBase, '3d_01', '.tif'], 32); 


