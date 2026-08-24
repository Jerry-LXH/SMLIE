function downFig = downsample(fig, upsample)
% down sample using pooling
poolSize = [upsample, upsample, upsample];
kernel = ones(poolSize) / prod(poolSize); 
downFig = convn(fig, kernel, 'valid');   
downFig = downFig(1:upsample:end, 1:upsample:end, 1:upsample:end); 
