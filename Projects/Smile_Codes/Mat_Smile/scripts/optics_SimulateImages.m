% This script is used for generating images.

% this script is aimed to:
% 1) simulate the 3D PSF based on NA, RI, pixel size
% 2) simulate the aberated 3D PSF based on Zernike Functions

%% parameters setting
% ---- physical parameters ----
lambda = 0.570; % um
NA = 1.49; % numerical aperture
RI = 1.515; % refraction index
wavenumber = 2*pi/lambda;

% ---- simulation parameters ----
upsample = 10; % sub-pixel grids defining.
mag = 150;
xSize = 257;
xSizeUp = xSize*upsample;
xPixelSize = 16/mag; % xPixelSize = 0.130; % um
xPixelSizeUp = xPixelSize/upsample;
zSize = 8; % how many defocus is allowed
zSizeUp = zSize*upsample;
zPixelSize = xPixelSize; 
zPixelSizeUp = zPixelSize/upsample;

% ---- emitters parameters ----
density = 0.0010;
photon_mean = 600;
photon_std = 20; % allow variance of brightness

% ---- background & noise ----
gn_mean = 2; % gaussian
gn_std = 0.5;
pn_mean = 4; % poisson
flagBg = 0;

para = 'test01';

%% save metadata
fileFolderOut = ['gendata/', num2str(para), '/'];
fileNameBase = [num2str(para), '_'];
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end
metadata = struct();
metadata.title = 'Simulation Metadata';
metadata.author = 'Jerry Ling';
metadata.date = datetime("now");
metadata.params = struct('lambda',lambda,  'NA', NA ,'RI', RI ,'unsample', upsample,'xSize', xSize ,'xPixelSize', xPixelSize ,'zSize', zSize ,'zPixelSize', zPixelSize, 'pointsdensity', density, 'photonmean', photon_mean, 'photonstd', photon_std, 'gn_mean',gn_mean,'gn_std',gn_std,'pn_mean',pn_mean);
jsonStr = jsonencode(metadata); 
fid = fopen([fileFolderOut, fileNameBase, 'metadata', '.json'], 'w');
fprintf(fid, '%s', jsonStr);
fclose(fid);

%% generate mask
fprintf('Generating pupil function... \n');
[r, theta, idx] =  gen.def_pupilcoor(xSizeUp, xPixelSizeUp, lambda, NA);
r0 = r(idx); 
theta0 = theta(idx);
pupilMask = zeros(xSizeUp, xSizeUp, 'single');
pupilMask(idx) = 1;
fprintf('Yes !\n');

fprintf('Generating/Loading Musk function... \n');
Mask = 0; % No mask
fprintf('Yes !\n');
if Mask == 0
    fprintf('NO MASK GIVEN. Cicular pupil assumed.\n');
    pupilFun = pupilMask;
else
    fprintf('USING GIVEN MUSK.\n')
    if ~isequal(size(pupilMask), size(Mask))
        error('Error: Mask and Grids must be the same size!\n');
    end
    pupilFun = pupilMask.*Mask;
end

%% calculate 3d psf according to mask
fprintf('Generating noabe PSF... \n');
PSF3d_noabe = gen.gen_PSF3D(xPixelSizeUp, xSizeUp, zPixelSizeUp, zSizeUp, lambda, RI, pupilFun, pupilMask); 
fprintf('Yes !\n');
% WriteTifStack(PSF3d_noabe,[fileFolderOut, fileNameBase,'PSF3d_noabe' ,'.tif'], 16);

%% generate image
emitters = [];
y = [];
tic;
PSF3d = PSF3d_noabe;

% Generate the random points with Gaussian brightness
[positions, emitterNumber] = gen.gen_random_points(xSizeUp, zSizeUp, upsample, density, photon_mean, photon_std);

% Calculete 2D image based on these points 
pointsBlur2D_noabe = gen.gen_blurred2D(xSizeUp,zSizeUp,positions,emitterNumber,PSF3d_noabe);
fprintf('Max pixel: %.2f',max(pointsBlur2D_noabe(:)));

% Downsampling 
pointsBlur2D_noabe = gen.downsample2D(pointsBlur2D_noabe, upsample); 

% Add noise
pointsBlur2D_noabe = gen.addnoise(pointsBlur2D_noabe, gn_mean, gn_std, pn_mean);

y = cat(3, y, pointsBlur2D_noabe);
t_total = toc;
fprintf('Size of image: ')
disp(size(y))
fprintf('%d points has been generated. \n',emitterNumber);
fprintf('Total time = %.3f s\n', t_total);
% fprintf('Total time = %.3f s, Per-frame = %.3f ms\n', t_total, t_per*1000);

%% save data
save([fileFolderOut, fileNameBase, 'data_', num2str(upsample),'UP','.mat'], 'pointsBlur2D_noabe','emitters');

%% Visualize
viz.max_project(y, 1, 'hot','Simulated_single_frame');

