function plotState(state, varargin)
%PLOTSTATE Visualize state analysis result for a single trace.
%
%   plotState(state)
%   plotState(state, 'useBgSub', true)
%   plotState(state, 'saveFig', true, 'saveDir', 'D:\result')
%   plotState(state, 'saveFig', true, 'saveDir', 'D:\result', 'saveName', 'trace_01')
%
% Input:
%   state - struct with fields:
%       trace
%       sequence
%       fitTrace
%       stateInfo
%       bleachFrame
%       lifetime
%
% Optional Name-Value pairs:
%   'useBgSub'   logical, whether to plot background-subtracted version
%                using state.background (default false)
%   'saveFig'    logical, whether to save the figure (default false)
%   'saveDir'    char/string, base directory where a subfolder "fig" will be created
%   'saveName'   char/string, file name for saving (default 'state_plot')
%
% Example:
%   plotState(state)
%   plotState(state, 'useBgSub', true)
%   plotState(state, 'saveFig', true, 'saveDir', 'D:\result', 'saveName', 'trace_01')

    % ---- parse options ------------------------------------------------
    p = inputParser;
    addParameter(p, 'useBgSub', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'saveFig', false, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'saveDir', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'saveName', 'state_plot', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});

    useBgSub = logical(p.Results.useBgSub);
    saveFigFlag = logical(p.Results.saveFig);
    saveDir = char(p.Results.saveDir);
    saveName = char(p.Results.saveName);

    % ---- check required fields ---------------------------------------
    requiredFields = {'trace','sequence','fitTrace','stateInfo','bleachFrame','lifetime'};
    for f = 1:numel(requiredFields)
        if ~isfield(state, requiredFields{f})
            error('Input state must contain field "%s".', requiredFields{f});
        end
    end

    % ---- data ---------------------------------------------------------
    y     = state.trace(:)';
    seq   = state.sequence(:)';
    fitTr = state.fitTrace(:)';
    sInfo = state.stateInfo;
    bFr   = state.bleachFrame;
    lt    = state.lifetime;

    % ---- background subtraction in plotting only ---------------------
    bgVal = 0;
    if useBgSub
        if ~isfield(state, 'background') || isempty(state.background) || isnan(state.background)
            error(['useBgSub=true, but state.background is missing, empty, ', ...
                   'or NaN.']);
        end
        bgVal = state.background;
        y     = y - bgVal;
        fitTr = fitTr - bgVal;
        yLabel = 'Background-subtracted intensity (Photons/frame)';
        titleTag = 'BG-subtracted';
    else
        yLabel = 'Intensity (Photons/frame)';
        titleTag = 'Raw';
    end

    F = numel(y);

    if numel(seq) ~= F
        error('state.sequence length does not match state.trace length.');
    end

    if numel(fitTr) ~= F
        error('state.fitTrace length does not match state.trace length.');
    end

    % ---- colour map for active states --------------------------------
    cmap = [0.00 0.45 0.74;   % S1
            0.85 0.33 0.10;   % S2
            0.47 0.67 0.19;   % S3
            0.93 0.69 0.13;   % S4
            0.49 0.18 0.56;   % S5
            0.30 0.75 0.93];  % S6
    bleachClr = [0.55 0.55 0.55];

    % ---- figure -------------------------------------------------------
    hFig = figure('Name', 'State Analysis', ...
                  'NumberTitle', 'off', ...
                  'Position', [100 150 850 300]);

    % ---- raw trace (gray) --------------------------------------------
    plot(1:F, y, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.5);
    hold on;

    % ---- fitted step trace (coloured segments) -----------------------
    segs = getSegments(seq);
    for s = 1:numel(segs)
        sg = segs(s);
        xL = sg.startF - 0.5;
        xR = sg.stopF  + 0.5;
        lv = fitTr(sg.startF);

        if sg.state == 0
            clr = bleachClr;
        else
            ci = mod(sg.state - 1, size(cmap,1)) + 1;
            clr = cmap(ci,:);
        end

        plot([xL xR], [lv lv], '-', 'Color', clr, 'LineWidth', 2);

        % vertical connector to previous segment
        if s > 1
            prevLv = fitTr(segs(s-1).startF);
            plot([xL xL], [prevLv lv], '-', ...
                 'Color', [0.35 0.35 0.35], 'LineWidth', 1);
        end
    end

    % ---- horizontal guide lines + right-side labels ------------------
    for j = 1:numel(sInfo)
        ci  = mod(sInfo(j).label - 1, size(cmap,1)) + 1;
        lvl = sInfo(j).meanIntensity - bgVal;

        plot([1 F], [lvl lvl], ':', ...
             'Color', cmap(ci,:), 'LineWidth', 0.8);

        text(F * 1.015, lvl, ...
             sprintf('S%d: %.0f \\pm %.0f', ...
                     sInfo(j).label, ...
                     sInfo(j).meanIntensity - bgVal, ...
                     sInfo(j).stdIntensity), ...
             'FontSize', 8, ...
             'FontWeight', 'bold', ...
             'Color', cmap(ci,:), ...
             'VerticalAlignment', 'middle', ...
             'Clipping', 'off');
    end

    % ---- bleach annotation -------------------------------------------
    if ~isnan(bFr)
        yl = ylim;
        xline(bFr - 0.5, '--', ...
              'Color', [0.9 0.1 0.1], 'LineWidth', 1.5);

        text(bFr + F*0.008, yl(1) + 0.92*(yl(2)-yl(1)), ...
             sprintf('Bleach\nLifetime = %d fr', lt), ...
             'FontSize', 9, ...
             'Color', [0.9 0.1 0.1], ...
             'FontWeight', 'bold', ...
             'VerticalAlignment', 'top');
    end

    % ---- axes ---------------------------------------------------------
    xlim([0.5, F * 1.12]);
    xlabel('Frame');
    ylabel(yLabel);

    if isnan(lt)
        ltStr = 'NaN';
    else
        ltStr = sprintf('%d fr', lt);
    end

    title(sprintf('%s | %d states | Lifetime = %s', ...
          titleTag, numel(sInfo), ltStr), 'FontSize', 10);

    hold off;
    box on;

    % ---- save figure --------------------------------------------------
    if saveFigFlag
        if isempty(saveDir)
            error('When ''saveFig'' is true, you must provide ''saveDir''.');
        end

        figDir = saveDir;
        if ~exist(figDir, 'dir')
            mkdir(figDir);
        end

        pngPath = fullfile(figDir, [saveName, '.png']);
        exportgraphics(hFig, pngPath, 'Resolution', 300);

    end
end


%% ======================================================================
% LOCAL HELPER — contiguous segment extraction
% ======================================================================
function segs = getSegments(seq)
    F = numel(seq);
    segs = struct('state', {}, 'startF', {}, 'stopF', {});

    if F == 0
        return;
    end

    idx = 1;
    segs(1).state  = seq(1);
    segs(1).startF = 1;

    for f = 2:F
        if seq(f) ~= seq(f-1)
            segs(idx).stopF = f - 1;
            idx = idx + 1;
            segs(idx).state  = seq(f);
            segs(idx).startF = f;
        end
    end

    segs(idx).stopF = F;
end