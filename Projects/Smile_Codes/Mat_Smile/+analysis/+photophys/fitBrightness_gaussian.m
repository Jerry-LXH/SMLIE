function [mu, sigma] = fitBrightness_gaussian(brightness)
% Fit brightness with normal distribution
% [brightness] should be a list in unit photons/sec.

    % Remove invalid values
    brightness = brightness(brightness > 0 & ~isnan(brightness));

    % MLE for Gaussian
    mu = mean(brightness);
    sigma = std(brightness);

    pd = makedist('Normal','mu',mu,'sigma',sigma);

    % Plot
    figure
    histogram(brightness,40,'Normalization','pdf')
    hold on

    x = linspace(min(brightness),max(brightness),200);
    y = pdf(pd,x);

    plot(x,y,'r','LineWidth',2)
    xlabel('Brightness (photons/frame)')
    ylabel('PDF')
    title(sprintf('Gaussian Fit: \\mu = %.2f, \\sigma = %.2f', mu, sigma))
    legend('Data','Gaussian fit')
    grid on

    fprintf('Brightness Gaussian fit: mu = %.6f, sigma = %.6f\n', mu, sigma);
end