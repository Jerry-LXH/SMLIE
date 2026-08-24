function blurredImg = conv_with_PSF(img, PSF3d)
% convolve the image with the PSF (3d)
PSF3d = PSF3d/sum(PSF3d(:,:,10), "all");
OTF = fftn(ifftshift(PSF3d));
blurredImg = real(ifftn(fftn(img).*OTF));
