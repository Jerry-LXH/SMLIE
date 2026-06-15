function [spots,stdWaveletAll,threshold] = find_spots_wavelet(stack, thresholdFactor, frameRange, ROI, ui)
%
% [spots,stdWaveletAll,threshold] = find_spots_wavelet(stack, para, ROI, ui)
%
% Detect spots in a stack of single-molecule image by first filtering the image
% stack with a wavelet filter and finding local maxima in the filtered stack
% through image dilation. Only local maxima with values above a specified threshold
% are kept.
%
% Input:
%   stack           -   3d-array of pixel values. 
%                           1st dimension (column): y-coordinate of the image plane 
%                           2nd dimension (row):    x-coordinate of the image plane 
%                           3rd dimension: frame number
%   thresholdFactor:    used to calculate the detection threshold.
%   frameRange:         array containing first and last frame
%                       of the image stack that should be
%                      	analysed
%   ROI             -   2-column vector containing coordinates of a region of interest where                    
%                       the spots should be detected. 1st column x-coordinates (horizontal),
%                       2nd column y-coordinates (vertical).
%   ui(optional)    -   Structure array containing all ui handles of
%                       TrackIt. Used to display tracking progress.
%
% Output:
%   spots           -   Cell array that has as many cells as there are
%                       frames in the image stack. Each cell contains a
%                       2-column array with x-coordinates in the 1st column
%                       and y-coordinates in the 2nd column.
%                       Eg. spots{10} = [5.1,10.2; 20.5,30.1] implies that
%                       there are two spots in the 10th frame of the stack,
%                        one at (x=5.1, y=10.2) and one at(x=20.5, y=30.1). 
%   stdWaveletAll   -   Estimated background noise which was used to
%                       calculate the detection threshold
%   threshold       -   Threshold above which spots are detected





if nargin == 4
    ui.editFeedbackWin.String = '';
end


origFeedbackWin = ui.editFeedbackWin.String(1:end,:);

%Get movie dimensions
stackSize = size(stack);

%Get number of frames in movie
nFrames = stackSize(3);

%Make sure last frame to analyze is not larger than the number of frames in the image stack 
if frameRange(2) > nFrames
    frameRange(2) = nFrames;
end

%Save ROI x and y coordinates separately
ROIx = ROI(:,1);
ROIy = ROI(:,2);

%----Cut out ROI of original stack to increase speed-----------------------

%Increase ROI to reduce border effects when filtering
additionalPixels = 5;

roiMaxX = ceil(min(max(ROIx)+additionalPixels,stackSize(2)));
roiMinX = floor(max(min(ROIx)-additionalPixels,1));
roiMaxY = ceil(min(max(ROIy)+additionalPixels,stackSize(1)));
roiMinY = floor(max(min(ROIy)-additionalPixels,1));

if ~isa(stack,'double') && ~isa(stack,'single')
    filtered = single(stack(roiMinY:roiMaxY,roiMinX:roiMaxX,:));
else
    filtered = stack(roiMinY:roiMaxY,roiMinX:roiMaxX,:);
end
%--------------------------------------------------------------------------

%Initialize cell array for spots
spots = repmat({zeros(0,2)},nFrames,1);

%Initialize variable to save the framewise standard deviations of the first
%wavelet filtered images
stdWavelet1 = zeros(nFrames,1);

ui.editFeedbackWin.String = char('Filtering...', origFeedbackWin);
drawnow

if double(get(gcf,'CurrentCharacter')) == 24
    %User pressed Strg+x so stop and return
    return
end

%Filter stack with wavelet filter
parfor k = frameRange(1):frameRange(2)
    [filtered(:,:,k), stdWavelet1(k)] = wavelet_filter(filtered(:,:,k));
end

%Stdev of the first order wavelet map. Used to estimate background noise.
stdWaveletAll = round(mean(stdWavelet1(stdWavelet1 ~= 0)),2);

%Calculate threshold for spot finding based on the estimated background
%noise of the first order wavelet image.
threshold = stdWaveletAll*thresholdFactor;

if double(get(gcf,'CurrentCharacter')) == 24
    %User pressed Strg+x so stop and return
    return
end

ui.editFeedbackWin.String = char('Finding Spots...', origFeedbackWin);
drawnow


parfor k = frameRange(1):frameRange(2)
    curFrameFiltered = filtered(:,:,k);
    
    %----------Local maxima finding through image dilation-----------------
    %Image dilation mask
    se = strel('square',4);
    
    %Identify local maxima as pixels which have the same value before and
    %after image dilation
    binary = curFrameFiltered == imdilate(curFrameFiltered,se);
    
    %Set spots below the user defined threshold to zero
    binary(curFrameFiltered < threshold) = 0;
    
    %Find spot positions
    [spotsY, spotsX] = find(binary);
    
    %Calculate position on non-cropped image
    spotsX = spotsX + roiMinX-1;
    spotsY = spotsY + roiMinY-1;
    
    %Look up which spots are inside ROI
    inROI = inpolygon(spotsX,spotsY,ROIx,ROIy);
    spotsCurFrame = [spotsX spotsY];
    
    %Save only spots which are inside ROI
    spots{k} = spotsCurFrame(inROI == 1,:);
end

if double(get(gcf,'CurrentCharacter')) == 24
    %User pressed Strg+x so stop and return
    return
end

end






