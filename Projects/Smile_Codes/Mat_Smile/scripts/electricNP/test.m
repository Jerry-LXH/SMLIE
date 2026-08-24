clear all; clc;
file_name = '/Volumes/SMILeSSD/UCNP charge data/20251028 98Yb2Er@Lu-10nm/2 Y@98Yb2Er@Lu-10nm-ZJ first no E -1V and 1V every5s 5cycles then no E.sif';
params = struct();
params.k_sigma          = 2.0;
params.bleach_time      = 10;
params.state_penalty    = 2.5;
params.state_min_seg_len = 20;
params.filter_end       = nan;
params.viz_enabled    = false;
params.emittersFIT_enabled = false;
% ... 其他需要覆盖默认值的参数

%% 运行 pipeline
fileResult = pipeline.special.pipeline_electric(file_name,params);

%% 画 20 条 trace 看看状态分配
for idx = 1:10
    analysis.photophys.plotState(fileResult.traceResult.states(idx));
    title(sprintf('Emitter #%d  (nStates = %d)', idx, fileResult.traceResult.nStatesVec(idx)));
end

%% 提取 metrics
metrics = fileResult.metrics;   % 根据你 pipeline 的实际输出结构调整字段名
fprintf('===== computeChargeMetrics 测试 =====\n');
fprintf('总 emitter 数: %d\n', numel(metrics.initial_brightness));

assert(isfield(metrics, 'initial_brightness'), '缺少字段: initial_brightness');
assert(isfield(metrics, 'normed_delta'),       '缺少字段: normed_delta');
assert(isfield(metrics, 'bright_dark_ratio'),  '缺少字段: bright_dark_ratio');
fprintf('✓ 三个输出字段均存在\n');

nEm = numel(metrics.initial_brightness);
assert(numel(metrics.normed_delta)      == nEm, '维度不一致: normed_delta');
assert(numel(metrics.bright_dark_ratio) == nEm, '维度不一致: bright_dark_ratio');
fprintf('✓ 维度一致 (%d × 1)\n', nEm);

valid_ib    = sum(isfinite(metrics.initial_brightness) & metrics.initial_brightness > 0);
valid_nd    = sum(isfinite(metrics.normed_delta));
valid_ratio = sum(isfinite(metrics.bright_dark_ratio));

fprintf('\n--- 有效数据统计 ---\n');
fprintf('initial_brightness > 0 : %d / %d\n', valid_ib, nEm);
fprintf('normed_delta 有效      : %d / %d\n', valid_nd, nEm);
fprintf('bright_dark_ratio 有效 : %d / %d\n', valid_ratio, nEm);

fprintf('\n--- 数值摘要 ---\n');
fprintf('initial_brightness : median = %.2f, range = [%.2f, %.2f]\n', ...
    median(metrics.initial_brightness, 'omitnan'), ...
    min(metrics.initial_brightness), max(metrics.initial_brightness));

fprintf('normed_delta       : median = %.4f, range = [%.4f, %.4f]\n', ...
    median(metrics.normed_delta, 'omitnan'), ...
    min(metrics.normed_delta), max(metrics.normed_delta));

fprintf('bright_dark_ratio  : median = %.4f, range = [%.4f, %.4f]\n', ...
    median(metrics.bright_dark_ratio, 'omitnan'), ...
    min(metrics.bright_dark_ratio), max(metrics.bright_dark_ratio));

% normed_delta 应 >= 0（亮度差不会为负）
if all(metrics.normed_delta(isfinite(metrics.normed_delta)) >= 0)
    fprintf('✓ normed_delta 全部 >= 0\n');
else
    warning('normed_delta 存在负值，请检查');
end

% bright_dark_ratio 应 > 0
if all(metrics.bright_dark_ratio(isfinite(metrics.bright_dark_ratio)) > 0)
    fprintf('✓ bright_dark_ratio 全部 > 0\n');
else
    warning('bright_dark_ratio 存在非正值，请检查');
end

%% 可视化
figure('Name', 'Charge Metrics Test', 'Position', [100 100 1200 400]);

subplot(1,3,1);
histogram(metrics.initial_brightness(isfinite(metrics.initial_brightness)), 30);
xlabel('Initial Brightness (counts)');
ylabel('# Emitters');
title('Initial Brightness');
xline(median(metrics.initial_brightness, 'omitnan'), 'r--', 'Median');

subplot(1,3,2);
histogram(metrics.normed_delta(isfinite(metrics.normed_delta)), 30);
xlabel('Normed \Delta');
ylabel('# Emitters');
title('Normalized Delta');
xline(median(metrics.normed_delta, 'omitnan'), 'r--', 'Median');

subplot(1,3,3);
histogram(metrics.bright_dark_ratio(isfinite(metrics.bright_dark_ratio)), 30);
xlabel('Bright/Dark Ratio');
ylabel('# Emitters');
title('上均差 / 下均差');
xline(1, 'k--', 'Symmetric');
xline(median(metrics.bright_dark_ratio, 'omitnan'), 'r--', 'Median');

sgtitle('computeChargeMetrics 输出验证');

%% 散点图: normed_delta vs bright_dark_ratio
figure('Name', 'Scatter', 'Position', [100 550 500 450]);
valid = isfinite(metrics.normed_delta) & isfinite(metrics.bright_dark_ratio);
scatter(metrics.normed_delta(valid), metrics.bright_dark_ratio(valid), 20, ...
    metrics.initial_brightness(valid), 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Normed \Delta');
ylabel('Bright/Dark Ratio');
colorbar; 
ylabel(colorbar, 'Initial Brightness');
yline(1, 'k--');
title('normed\_delta vs bright\_dark\_ratio');

fprintf('\n===== 测试完成 =====\n');