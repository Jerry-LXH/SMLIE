function dConst = defocus_wavefront(xPixelSize, Sx, lambda, RI, pupilMask)
% generate the constant term of defocus wavefront

freSampling = 1/xPixelSize; % length^-1
pixelSizePhase = freSampling/Sx;

% calculate defocus function: dConst
dConst = zeros(Sx,Sx, 'single');
Sox = (Sx+1)/2; % Sox == Soy
for i = 1:Sx
    for j = 1:Sx
        if(pupilMask(i,j)==1)
            rSQ = (i-Sox)^2 + (j-Sox)^2;
            rSQ = rSQ * pixelSizePhase^2;
            dConst(i,j) = sqrt(1-(lambda/RI)^2*rSQ); 
        end
    end
end
dConst = 2*pi*RI/lambda*dConst; % in wavenumber unit
