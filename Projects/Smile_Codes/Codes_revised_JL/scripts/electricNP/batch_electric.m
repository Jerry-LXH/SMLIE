%% batch_charge.m — 批量处理多个 SIF 文件
clear; clc;

%% ===== 文件列表 =====
    file_list = {
        '/Volumes/SMILeSSD/UCNP charge data/20260401.../1 98Yb2Er-ZJ 4mix....sif';
        '/Volumes/SMILeSSD/UCNP charge data/20260401.../2 98Yb2Er-ZJ 4mix....sif';
        % ... 添加更多文件
    };

%% ===== 自动搜索多个目录 =====
    root_dirs = {
        '/Volumes/SMILeSSD/UCNP charge data/20251028 98Yb2Er@Lu-2.5nm';
        '/Volumes/SMILeSSD/UCNP charge data/20251028 98Yb2Er@Lu-6nm';
        '/Volumes/SMILeSSD/UCNP charge data/20251028 98Yb2Er@Lu-10nm';
        '/Volumes/SMILeSSD/UCNP charge data/20251029 98Yb2Er@Lu-10nm';
        '/Volumes/SMILeSSD/UCNP charge data/small';
        '/Volumes/SMILeSSD/UCNP charge data/20251222+1223 small-98Yb2Er no NaYF4 300mA-laser/1222';
        '/Volumes/SMILeSSD/UCNP charge data/20251222+1223 small-98Yb2Er no NaYF4 300mA-laser/1223';
        '/Volumes/SMILeSSD/UCNP charge data/20260317 98Yb2Er mix 2.5nm 6nm (shell) 300mA-laser';
        '/Volumes/SMILeSSD/UCNP charge data/20260401 Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and 98Yb2Er@Lu(without Y) 300mA-laser/3mix Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu)';
        '/Volumes/SMILeSSD/UCNP charge data/20260401 Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and 98Yb2Er@Lu(without Y) 300mA-laser/4mix Y@98Yb2Er@Lu mix 2.5nm 6nm 10nm(Lu) and small98Yb2Er@Lu(without Y)';
    };

    sif_pattern = '**/*every5s 5cycles then no E.sif';

    file_list = {};
    for r = 1:numel(root_dirs)
        d = dir(fullfile(root_dirs{r}, sif_pattern));
        file_list = [file_list; fullfile({d.folder}, {d.name})']; %#ok<AGROW>
    end

%% ===== 统一参数 =====
    params = struct();
    params.k_sigma          = 2.0;
    params.bleach_time      = 10;
    params.state_penalty    = 2.5;
    params.state_min_seg_len = 20;
    params.filter_end       = nan;
    params.viz_enabled    = false;
    params.emittersFIT_enabled = false;
    % ... 其他需要覆盖默认值的参数

%% ===== 逐文件处理 =====
    nFiles = numel(file_list);
    results = cell(nFiles, 1);

    for f = 1:nFiles
        fprintf('\n========== Processing file %d / %d ==========\n', f, nFiles);
        try
            results{f} = pipeline.special.pipeline_electric(file_list{f}, params);
        catch ME
            fprintf('[ERROR] File %d failed: %s\n', f, ME.message);
            results{f} = struct('file_name', file_list{f}, 'error', ME.message);
        end
    end

    fprintf('\n===== All processing complete =====\n');
