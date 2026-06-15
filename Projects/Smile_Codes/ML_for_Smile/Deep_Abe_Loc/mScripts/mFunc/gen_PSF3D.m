function PSF3d = gen_PSF3D(xPixelSize, xSize, zPixelSize, zSize, lambda, RI, pupilFun, pupilMask)
% Generate 3D PSF based on pupil function
dConst = defocus_wavefront(xPixelSize, xSize, lambda, RI, pupilMask); 
PSF3d = zeros(xSize, xSize, zSize, 'single');
Soz = (zSize+1)/2; 
for i = 1:zSize
    zPos = (i - Soz) * zPixelSize;
    pupilFun2D = pupilFun .*exp(1i*zPos*dConst);
    prf = fftshift(ifft2(ifftshift(pupilFun2D)));
    PSF3d(:,:,i) = abs(prf).^2;
end
PSF3d = PSF3d/max(PSF3d(:)); % normalization