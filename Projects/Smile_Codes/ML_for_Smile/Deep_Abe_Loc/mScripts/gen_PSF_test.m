addpath('mFunc');

lambda = 0.532;
pIn = 3:15;
NA = 1.1;
RI = 1.33;
upsample = 1;
xSize = 128;
xSizeUp = xSize*upsample;
xPixelSize = 0.130; % xPixelSize = 0.130; % um
xPixelSizeUp = xPixelSize/upsample;
zSize = 32*4;
zSizeUp = zSize*upsample;
zPixelSize = 0.130; 
zPixelSizeUp = zPixelSize/upsample;
aValue = 0.4;
bgValue = 95;
flagBg = 0;% 0.8NA: 1; 1.1NA: 0; 0.71NA: 0
zernIndex = 3:15;
zernType = 'random1';
% conType = 'ANSI';
zernNum = length(zernIndex);
density = 0.0005;
mean = 6000;
std = 2000;
para = 'sphr';
wavenumber = 2*pi/lambda;

fileFolderOut = ['simuAndModel/PSF', '/'];
fileNameBase = ['PSF_', para, '_'];
if isequal(exist(fileFolderOut, 'dir'),7)
    disp(['output folder:' fileFolderOut]);
else
    mkdir(fileFolderOut);
    disp(['output folder created:' fileFolderOut]);
end

fprintf('Generating pupil function... \n');
[r, theta, idx] =  def_pupilcoor(xSizeUp, xPixelSizeUp, lambda, NA);
r0 = r(idx); 
theta0 = theta(idx);
pupilMask = zeros(xSizeUp, xSizeUp, 'single');
pupilMask(idx) = 1;
fprintf('Yes !\n');

pupilFun = pupilMask;

% Generate zernike coefficients
coeffs = zeros(1, zernNum, 'single'); 
coeffs(10)= aValue ;
% coeffs = gen_zern_coeffs(pIn,aValue,zernType);
coeffs = wavenumber * coeffs; % Make into wavenumber units

% Generate aberrated pupil function
phi = zeros(xSizeUp, xSizeUp, 'single');
phi(idx) = create_wavefront(zernIndex, coeffs, r0, theta0); 
pupilFunAbe = pupilFun.*exp(1i*phi);
PSF3d = gen_PSF3D(xPixelSizeUp, xSizeUp, zPixelSizeUp, zSizeUp, lambda, RI, pupilFunAbe, pupilMask);
% PSF3d = coeffPSF_modified(zernIndex, coeffs_zero, xSizeUp, xPixelSizeUp, lambda, NA, zSizeUp, zPixelSizeUp, RI);  % 生成不同GL模态的PSF

WriteTifStack(PSF3d,[fileFolderOut, fileNameBase, '.tif'], 16);


