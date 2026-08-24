function tau = fitLifetime_exp(lifetime)
% Fit lifetime with exponential distribution
% [lifetime] should be a list, in seconds unit.

    % Remove invalid values
    lifetime = lifetime(lifetime > 0 & ~isnan(lifetime));
    n = numel(lifetime);

    % MLE for exponential
    tau = mean(lifetime);

    % Create PDF object
    pd = makedist('Exponential','mu',tau);

    % Plot
    figure
    nBins = 25;
    h = histogram(lifetime, nBins);
    hold on

    % 将 PDF 缩放到 counts：scale = N * binWidth
    binWidth = h.BinWidth;
    N = numel(lifetime);
    scaleFactor = N * binWidth;

    x = linspace(0, max(lifetime), 200);
    y = pdf(pd, x) * scaleFactor;

    plot(x, y, 'r', 'LineWidth', 2)
    xlabel('Lifetime (s)')
    ylabel('Counts')
    title(sprintf('Exponential Fit: \\tau = %.4f s, n=%d', tau, n))
    legend('Data', 'Exponential fit')
    grid on

    fprintf('Lifetime exponential fit:\n');
    fprintf('  n = %d\n', n);
    fprintf('  tau = %.6f s\n', tau);
end