function [PSF3d] = coeffPSF(zerIndex, coeffs, xSize, xPixelSize, lambda, NA, zSize, zPixelSize, RI)
% calculate the PSF
% coeffs: um
% all pixel size: um
% lambda: um
% This calculation rely on the effct that the impulse response function h(x,y) is the Fourier transform of the calibrated pupil function P(u,v), whose coordinates are in the frequency domain. Given NA, lambda, and pixel size, we can calculate the pupil function in frequency domain and then the PSF. The 3D PSF is calculated by the convolution of the 2D PSF with the defocus wavefront, i.e., the multiplication of the pupil function with the phase term of the defocus wavefront followed then by fft.
% Note that each pixel in the pupil function is related to the real pupil size by (freq_x_pixel)*(freq_pixel_size)*(lambda*zi/n), while the def_pupilcoor function gives the normalized frequency.
 

ySize = xSize;
wavenumber = 2*pi/lambda;
coeffs = wavenumber * coeffs; % turn into wavenumber units

% pupil coordinates, they share x-y indice as card-coordinates
[r, theta, idx] =  def_pupilcoor(xSize, xPixelSize, lambda, NA);
r0 = r(idx); % 
theta0 = theta(idx);

% generate distorted wavefront in phase units
phi = zeros(xSize, ySize, 'single');
phi(idx) = create_wavefront(zerIndex, coeffs, r0, theta0); % in phase unit: pi


% add artifical mask

% eng_wavefront = mod(2 * theta + 5 * r1.^2, 2 * pi);
% disp(theta);
eng_wavefront = generate_DH_PSF(r, theta, idx, [1,2,3,4,5], [1,3,5,7,9], 0.15, 30);
% phi = phi + eng_wavefront;
% disp(eng_wavefront);

% express distortion in complex units
pupilMask = zeros(xSize, ySize, 'single');
pupilMask(idx) = 1;
% pupilFun = pupilMask.*exp(1i*phi).*eng_wavefront;
pupilFun = pupilMask.*exp(1i*(phi+eng_wavefront));

% disp(abs(eng_wavefront));
% 2D PSF at focal plane
% prf = fftshift(ifft2(ifftshift(pupilFun)));
% PSF2d = abs(prf).^2;

% 3D PSF
dConst = defocus_wavefront(xPixelSize, xSize, ySize, lambda, RI, pupilMask);
PSF3d = zeros(xSize, ySize, zSize, 'single');
Soz = (zSize+1)/2;
for i = 1:zSize
    zPos = (i - Soz) * zPixelSize;
    pupilFun2D = pupilFun .*exp(1i*zPos*dConst);
    prf = fftshift(ifft2(ifftshift(pupilFun2D)));
    PSF3d(:,:,i) = abs(prf).^2;
end

% normalization
% PSF2d = PSF2d/max(PSF2d(:));
PSF3d = PSF3d/max(PSF3d(:));

