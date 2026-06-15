function [SNRspots, intPeak] = calc_snr(spots, originalIm)

% [SNRspots, intPeak] = calc_snr(spots, originalIm)
%
% Calculate the signal-to-noise ratio (SNR) of single-molecule detections
%
% Input: 
%   spots       -   2d-array containing the xy-coordinates of spots. x-coordinates 
%                   (horizontal direction) in first row and y-coordinates (vertical direction) in seconds row.
%   originalIm  -   Grayscale image where the SNR is calculated on
%
% Output:
%   SNRspots    -   list of SNRs
%   intPeak     -   list containing the maximum pixel value in a window of 17x17 pixels around
%                   the spot position

%Radius in which the original image pixel values are set to zero around the
%spot positions to calculate the background noise
spotRadius = 5;

%Half size of the window which is cut out of the original image where the SNR is be calculated
halfWindowSize = 8;

%Convert original image to double
originalIm = double(originalIm);

%Get number of spots
nSpots = size(spots,1);

%Get image dimensions
imageSize = size(originalIm);

[columnsInImage, rowsInImage] = meshgrid(1:imageSize(2), 1:imageSize(1));

%Create an image mask containing the spot positions plus a disc around each
%spot position with the radius defined by spotRadius
spotMask = false(imageSize);

%Iterate through all spots
for spotIdx = 1:nSpots
    curSpotY = round(spots(spotIdx, 2));
    curSpotX = round(spots(spotIdx, 1));
    spotMask =  spotMask | (rowsInImage - curSpotY).^2 ...
        + (columnsInImage - curSpotX).^2 <= spotRadius.^2;
end

%Create image containing only background by setting all values in a radius
%around each spot position to 0.
bgMask = spotMask == 0;
bgIm = originalIm.*bgMask;

%Initialize array of peak intensities and SNRs
intPeak = zeros(1,nSpots);
SNRspots = zeros(1,nSpots);

%Iterate through all spots
for spotIdx = 1:nSpots
    
    %Get coordinates of current spot
    curSpotY        = round(spots(spotIdx, 2));
    curSpotX        = round(spots(spotIdx, 1));
    
    %Get the boundaries of the window where the SNR is calculated
    xMin            = max(curSpotX - halfWindowSize, 1);
    xMax            = min(curSpotX + halfWindowSize, size(originalIm, 2));
    yMin            = max(curSpotY - halfWindowSize, 1);
    yMax            = min(curSpotY + halfWindowSize, size(originalIm, 1));
    
    %Get background image of current spot
    curSpotBgIm     = bgIm(yMin:yMax, xMin:xMax);
    
    %Define a radius around the spot position in which the pixels are
    %considered for getting the maximum and mean intensity of the spot
    spotRadius = 1;
    
    %Create a spot mask for cutting the spot out of the original image
    [columnsInSpotImage, rowsInSpotImage] = meshgrid(1:xMax-xMin+1, 1:yMax-yMin+1);
    spotMask =  (rowsInSpotImage - (yMax-yMin)/2-1).^2 ...
        + (columnsInSpotImage - (xMax-xMin)/2-1).^2 <= spotRadius.^2;
    
    %Cut out small image from the original image where only pixels inside a
    %radius aroung the spot position are nonzero
    curSpotIm = spotMask.*originalIm(yMin:yMax, xMin:xMax);
        
    %Get pixel values of background image
    bgPixelValues       = curSpotBgIm(curSpotBgIm~=0);
    
    %Get mean background intensity
    meanBgI             = mean(bgPixelValues);
    
    %Get standard deviation of background 
    stdBg               = std(bgPixelValues);
    
    %Get values of the pixels around the spot position
    spotPixelValues     = curSpotIm(curSpotIm ~= 0);  
    
    %Calculate mean intensity of the pixels around the spot position
    meanSpotI           = mean(spotPixelValues);
    
    %Get the pixel with the highest intensity
    intPeak(spotIdx)    = max(spotPixelValues);
    
    %Calculate the SNR
    SNRspots(spotIdx)   = round((meanSpotI - meanBgI)/stdBg,2);
end

end








