function pd = fitBrightness_lognormal(brightness)
% Fit brightness distribution with lognormal model

    brightness = brightness(brightness > 0 & ~isnan(brightness));
    n = numel(brightness);

    pd = fitdist(brightness,'Lognormal');

    mu = pd.mu;
    sigma = pd.sigma;

    figure
    h = histogram(brightness,40);
    hold on

    x = linspace(min(brightness),max(brightness),300);
    binWidth = h.BinWidth;
    y = pdf(pd,x) * n * binWidth;

    plot(x,y,'r','LineWidth',2)
    xlabel('Brightness (photons/sec)')
    ylabel('Count')
    title(sprintf('Lognormal Fit: \\mu=%.3f, \\sigma=%.3f, n=%d', mu, sigma, n))
    legend('Data','Lognormal fit')
    grid on

    fprintf('Lognormal fit:\n');
    fprintf('  n = %d\n', n);
    fprintf('  mu_log = %.6f\n', mu);
    fprintf('  sigma_log = %.6f\n', sigma);
end