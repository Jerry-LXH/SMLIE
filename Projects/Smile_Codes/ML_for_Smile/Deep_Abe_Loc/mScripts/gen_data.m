% this script is aimed to:
% 1) simulate the 3D PSF based on NA, RI, pixel size
% 2) simulate the aberated 3D PSF based on Zernike Functions

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
zSize = 32;
zSizeUp = zSize*upsample;
zPixelSize = 0.130; 
zPixelSizeUp = zPixelSize/upsample;
aValue = 0.1; % std(phase), aValue/2*12=rms
gn_mean = 2;
gn_std = 0.5;
pn_mean = 5;
flagBg = 0;% 0.8NA: 1; 1.1NA: 0; 0.71NA: 0
zernIndex = 3:15;
zernType = 'random1';
% conType = 'ANSI';
zernNum = length(zernIndex);
density = 0.0001;
photon_mean = 6000;
photon_std = 500;
para = '5_very_low_density_not_downsampled';

%%% SAVE METADATA
fileFolderOut = ['simuAndModel/para', num2str(para), '/'];
fileNameBase = ['para', num2str(para), '_'];
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
metadata.params = struct('lambda',lambda, 'pIn', pIn , 'NA', NA ,'RI', RI ,'unsample', upsample,'xSize', xSize ,'xPixelSize', xPixelSize ,'zSize', zSize ,'zPixelSize', zPixelSize, 'aValue', aValue , 'photonmean', photon_mean, 'photonstd', photon_std, 'pointsdensity', density, 'zernIndex', zernIndex, 'zernType', zernType, 'gn_mean',gn_mean,'gn_std',gn_std,'pn_mean',pn_mean);
jsonStr = jsonencode(metadata); 
fid = fopen([fileFolderOut, fileNameBase, 'metadata', '.json'], 'w');
fprintf(fid, '%s', jsonStr);
fclose(fid);


y = [];
X= [];
X_coeffs = [];
emitters = [];
rounds = 5; 

%%% GENERATING MASK

fprintf('Generating pupil function... \n');
[r, theta, idx] =  def_pupilcoor(xSizeUp, xPixelSizeUp, lambda, NA);
r0 = r(idx); 
theta0 = theta(idx);
pupilMask = zeros(xSizeUp, xSizeUp, 'single');
pupilMask(idx) = 1;
fprintf('Yes !\n');

fprintf('Generating/Loading Musk function... \n');
% Mask = generate_DH_PSF(r, theta, idx, [1,2,3,4,5], [1,3,5,7,9], 0.15, 1);
% mod = abs(Mask);
% WriteTifStack(mod,[fileFolderOut, fileNameBase,'Mask_mod', '.tif'], 32);
% angle = angle(Mask);
% WriteTifStack(angle,[fileFolderOut, fileNameBase,'Mask_angle', '.tif'], 32);

% mod = imread(['simuAndModel/Mask/DH_mod.tif']);
% ang = imread(['simuAndModel/Mask/DH_angle.tif']);
% Mask = mod .* exp(1i.*ang);

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

% This is the same across rounds
fprintf('Generating noabe PSF... \n');
PSF3d_noabe = gen_PSF3D(xPixelSizeUp, xSizeUp, zPixelSizeUp, zSizeUp, lambda, RI, pupilFun, pupilMask); 
WriteTifStack(PSF3d_noabe,[fileFolderOut, fileNameBase,'PSF3d_noabe' ,'.tif'], 16);
fprintf('Yes !\n');


%%% RANDOM ABE AND NOISE

wavenumber = 2*pi/lambda;
for round = 1:rounds

    fprintf('Round %d ...\n',round);

    % Generate random zernike coefficients
    coeffs = gen_zern_coeffs(pIn,aValue,zernType);
    coeffs = wavenumber * coeffs; % Make into wavenumber units

    % Generate aberrated pupil function
    phi = zeros(xSizeUp, xSizeUp, 'single');
    phi(idx) = create_wavefront(zernIndex, coeffs, r0, theta0); 
    pupilFunAbe = pupilFun.*exp(1i*phi);
    
    % Generate 3D PSF (aberrated)
    PSF3d = gen_PSF3D(xPixelSizeUp, xSizeUp, zPixelSizeUp, zSizeUp, lambda, RI, pupilFunAbe, pupilMask);

    % Generate the random points with Gaussian brightness
    [positions, emitterNumber] = gen_random_points(xSizeUp, zSizeUp, upsample, density, photon_mean, photon_std);

    % Calculete 2D image based on these points 
    pointsBlur2D_noabe = gen_blurred2D(xSizeUp,zSizeUp,positions,emitterNumber,PSF3d_noabe);
    pointsBlur2D = gen_blurred2D(xSizeUp,zSizeUp,positions,emitterNumber,PSF3d);
    if round == 1 
        disp(max(pointsBlur2D_noabe(:)));
        disp(max(pointsBlur2D(:)));
        disp(sqrt(mean(phi(idx).^2)));
    end
    % Add noise
    pointsBlur2D = addnoise(pointsBlur2D, gn_mean, gn_std, pn_mean);

    % Downsampling 
    pointsBlur2D_noabe = downsample2D(pointsBlur2D_noabe, upsample); 
    pointsBlur2D = downsample2D(pointsBlur2D, upsample); 

    % Normalization and save to tif stack
    pointsBlur2D_noabe = pointsBlur2D_noabe/max(pointsBlur2D_noabe(:)); % normalize by max value
    pointsBlur2D = pointsBlur2D/max(pointsBlur2D(:));  
    
    % Save a typical blurred/unblurred image (low res)
    if round == 1 
        WriteTifStack(pointsBlur2D_noabe,[fileFolderOut, fileNameBase, '2d_points_noabe_','rounds', num2str(round), '.tif'], 32);
        WriteTifStack(pointsBlur2D,[fileFolderOut, fileNameBase, '2d_points_', 'rounds', num2str(round), '.tif'], 32);
        WriteTifStack(PSF3d,[fileFolderOut, fileNameBase,'PSF3d_abe' ,'.tif'], 16);
    end

    % Store data
    y = cat(3, y, pointsBlur2D_noabe);
    X = cat(3, X, pointsBlur2D);
    coeffs = coeffs(:); 
    X_coeffs = cat(2, X_coeffs, coeffs); % num_coeffs * rounds
    emitters = cat(3, emitters, positions);
end

fprintf('%d points has been generated. \n',emitterNumber);
fprintf('%d rounds has been generated.\n',rounds)
fprintf('Size of X: ')
disp(size(X))
fprintf('Size of emitters:')
disp(size(emitters))

% Save 1. blurred and unblurred image stakes 2. Zernike coeffs 3. emitters
save([fileFolderOut, fileNameBase, 'dataset_', num2str(rounds), 'rounds_',num2str(upsample),'UP','.mat'], 'y', 'X', 'X_coeffs','emitters');



