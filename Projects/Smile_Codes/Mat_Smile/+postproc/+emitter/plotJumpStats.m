function plotJumpStats(stats)
% Quick diagnostic plots for filterEmitters_jumping

figure('Position', [100 100 1200 400]);

% 1) Eccentricity histogram
subplot(1,3,1);
ecc = stats.eccentricity(~isnan(stats.eccentricity));
kp  = stats.keep(~isnan(stats.eccentricity));
histogram(ecc(kp),  30, 'FaceColor', [0.3 0.6 1], 'EdgeColor', 'none'); hold on;
histogram(ecc(~kp), 30, 'FaceColor', [1 0.3 0.3], 'EdgeColor', 'none');
xlabel('\lambda_{max} / \lambda_{min}');
ylabel('Count');
legend('Kept', 'Filtered');
title('Eccentricity distribution');
set(gca, 'YScale', 'log');

% 2) σ_major vs eccentricity
subplot(1,3,2);
valid = ~isnan(stats.eccentricity);
scatter(stats.sigma_major(valid & stats.keep), ...
        stats.eccentricity(valid & stats.keep), ...
        8, [0.3 0.6 1], 'filled', 'MarkerFaceAlpha', 0.4); hold on;
scatter(stats.sigma_major(valid & ~stats.keep), ...
        stats.eccentricity(valid & ~stats.keep), ...
        12, [1 0.3 0.3], 'filled');
xlabel('\sigma_{major} (px)');
ylabel('Eccentricity');
title('\sigma_{major} vs eccentricity');

% 3) Major axis orientation — double-PSF should cluster
subplot(1,3,3);
filt_angle = stats.major_angle(valid & ~stats.keep);
if ~isempty(filt_angle)
    polarhistogram(deg2rad(filt_angle), 36, ...
        'FaceColor', [1 0.3 0.3], 'EdgeColor', 'none');
    title('Major axis angle (filtered)');
else
    text(0.5, 0.5, 'No filtered emitters', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
end

sgtitle('filterEmitters\_jumping diagnostics');
end