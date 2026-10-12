function plotEmitterStatistics(stats,bin_num)

    if nargin < 2
        bin_num = 40;
    end
    n = numel(stats.survival_sec);
    figure

    subplot(5,1,1)
    histogram(stats.survival_sec,bin_num)
    xlabel('Survival time (s)')
    ylabel('Occurrence')
    title(sprintf('Lifetime distribution (emitters) [n=%d]', n))

    subplot(5,1,2)
    histogram(stats.brightness_mean,bin_num)
    xlabel('Photons/sec')
    ylabel('Occurrence')
    title(sprintf('Brightness distribution (emitters) [n=%d]', n))

    subplot(5,1,3)
    histogram(stats.brightness_sum,bin_num)
    xlabel('Photons')
    ylabel('Occurrence')
    title(sprintf('Total photon distribution (emitters) [n=%d]', n))

    subplot(5,1,4)
    histogram(stats.sigma_mean,bin_num)
    xlabel('Sigma (Pixel)')
    ylabel('Occurrence')
    title(sprintf('Width distribution (emitters) [n=%d]', n))

    subplot(5,1,5)
    histogram(stats.sigma_loc_mean,bin_num)
    xlabel('Sigma loc')
    ylabel('Occurrence')
    title(sprintf('Uncertainty distribution (emitters) [n=%d]', n))
end