function [pupilFun, r0, theta0, idx] = gen_grids(xSize, xPixelSize, lambda, NA, Mask)
% Generate pupil grids based either on a circle

[r, theta, idx] =  def_pupilcoor(xSize, xPixelSize, lambda, NA);
r0 = r(idx); 
theta0 = theta(idx);

pupilMask = zeros(xSize, xSize, 'single');
pupilMask(idx) = 1;

if ~exist('Mask', 'var')
    fprintf('NO MASK GIVEN. Cicular pupil assumed.\n');
    pupilFun = pupilMask;
else
    fprintf('USING GIVEN MUSK.\n')
    if ~isequal(size(pupilMask), size(Mask))
        error('Error: Mask and Grids must be the same size!\n');
    end
    pupilFun = pupilMask.*Mask;
end