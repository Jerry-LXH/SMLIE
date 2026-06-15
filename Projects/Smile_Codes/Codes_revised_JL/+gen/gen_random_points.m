function [positions, emitterNumber] = gen_random_points(xSizeUp, zSizeUp, upsample, density, mean, std)
% This function is aimed to generate random points in 3D space based on upsampling factor
% density: number of points per pixel^3
% fprintf('Generating random points... \n');
emitterNumber = round(density*(xSizeUp/upsample)*(xSizeUp/upsample)*(zSizeUp/upsample));
Sox = xSizeUp/2;
Soz = zSizeUp/2;
x = randi([1, xSizeUp], emitterNumber, 1)- Sox;
y = randi([1, xSizeUp], emitterNumber, 1)- Sox;
z = randi([1, zSizeUp], emitterNumber, 1)- Soz;

brightness = mean + std * randn(emitterNumber, 1);
brightness = max(0, min(2147483647, brightness));
positions = [x, y, z, brightness];
% fprintf('%d points has been generated. \n',emitterNumber);
