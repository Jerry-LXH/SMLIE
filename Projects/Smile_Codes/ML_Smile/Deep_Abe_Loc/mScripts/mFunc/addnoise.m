function imageNoise = addnoise(image, average, sigma, lambda)
% Add shot noise and Gaussian noise to the image
imageNoise = zeros(size(image));
imageNoise = imageNoise + poissrnd(image) + normrnd(average, sigma, size(image));% + poissrnd(lambda, size(image));
end