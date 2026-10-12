function states = analyzeStates(traceMat, method, varargin)
%ANALYZESTATES Detect discrete intensity states in traces.
%
%   states = analyzeStates(traceMat)
%   states = analyzeStates(traceMat, method, 'Name', Value, ...)
%
% INPUTS
%   traceMat   N x F numeric matrix, one trace per row
%   method     'HMM' (default) | 'changepoint'
%
% Name-Value options (common):
%   'bleachTail'  min frames at trace-end to confirm bleach   (default 10)
%   'verbose'     logical                                     (default false)
%
% HMM-specific:
%   'maxStates'   max HMM states to try                       (default 5)
%   'minStates'   min HMM states to try                       (default 2)
%   'maxIter'     max EM iterations per fit                   (default 200)
%   'nRestarts'   EM random restarts per K                    (default 3)
%   'bicPenalty'  BIC penalty multiplier                      (default 1.0)
%
% Changepoint-specific:
%   'penalty'     BIC penalty multiplier                      (default 2)
%   'minSegLen'   minimum segment length in frames            (default 5)
%   'mergeThr'    merge threshold in noise-std units          (default 2.0)
%   'refineIter'  changepoint refinement iterations           (default 5)
%
% State merging (common):
%   'bgThreshold' background intensity threshold; []=auto, 0=disable
%   'bgAutoSigma' bgLevel + bgAutoSigma * noise              (default 3)
%   'minStateSep' merge if state mean difference < this*noise (default 2.0)
%   'bgQuantile'  lowest quantile used for bg estimation      (default 0.10)
%
% OUTPUT
%   states      N x 1 struct array with fields:
%       .trace
%       .sequence
%       .transitions
%       .nStates
%       .stateInfo
%       .fitTrace
%       .bleachFrame
%       .method
%       .lifetime

% ---- input check ------------------------------------------------------
if nargin < 2 || isempty(method)
    method = 'HMM';
end

if ~isnumeric(traceMat) || ndims(traceMat) ~= 2
    error('traceMat must be a numeric N x F matrix.');
end

[N, F] = size(traceMat); %#ok<ASGLU>

% ---- options ----------------------------------------------------------
p = inputParser;

% Common
addParameter(p, 'bleachTail', 10);
addParameter(p, 'verbose', false);

% HMM
addParameter(p, 'maxStates', 5);
addParameter(p, 'minStates', 2);
addParameter(p, 'maxIter', 200);
addParameter(p, 'nRestarts', 3);
addParameter(p, 'bicPenalty', 1.0);

% Changepoint
addParameter(p, 'penalty', 2.0);
addParameter(p, 'minSegLen', 5);
addParameter(p, 'mergeThr', 2.0);
addParameter(p, 'refineIter', 5);

% State merging
addParameter(p, 'bgThreshold', []);
addParameter(p, 'bgAutoSigma', 3);
addParameter(p, 'minStateSep', 2.0);
addParameter(p, 'bgQuantile', 0.10);

parse(p, varargin{:});
opts = p.Results;

% ---- init output ------------------------------------------------------
states = repmat(struct( ...
    'trace', [], ...
    'sequence', [], ...
    'transitions', [], ...
    'nStates', [], ...
    'stateInfo', [], ...
    'fitTrace', [], ...
    'background', [], ...
    'traceBgSub', [], ...
    'fitTraceBgSub', [], ...
    'bleachFrame', [], ...
    'method', '', ...
    'lifetime', []), N, 1);

tStart = tic;

% ---- dispatch ---------------------------------------------------------
switch upper(method)
    case 'HMM'
        states = runHMM(traceMat, states, opts);
    case 'CHANGEPOINT'
        states = runChangepoint(traceMat, states, opts);
    otherwise
        error('analyzeStates:badMethod', ...
            'Method "%s" not implemented. Supported: HMM, changepoint.', method);
end

fprintf('Total time = %.3f s\n', toc(tStart));

end


%% CHANGEPOINT METHOD
function states = runChangepoint(traceMat, states, opts)

N = size(traceMat, 1);

for n = 1:N
    if opts.verbose
        fprintf('  Trace %d / %d (changepoint) ...\n', n, N);
    end

    y = double(traceMat(n, :));
    F = length(y);

    % -- precompute cumulative sums for O(1) segment RSS -----------------
    cumY  = [0, cumsum(y)];
    cumY2 = [0, cumsum(y.^2)];

    % -- Step 1: recursive binary segmentation ---------------------------
    cpInter  = binarySeg(cumY, cumY2, 1, F, opts.penalty, opts.minSegLen);
    cpBounds = sort(unique([1, cpInter(:)', F+1]));

    % -- Step 2: refine changepoint positions ----------------------------
    cpBounds = refineCP(cumY, cumY2, cpBounds, opts.minSegLen, opts.refineIter);

    % -- Step 3: segment statistics --------------------------------------
    nSeg    = numel(cpBounds) - 1;
    segMean = zeros(1, nSeg);
    segLen  = zeros(1, nSeg);
    segVar  = zeros(1, nSeg);

    for s = 1:nSeg
        a = cpBounds(s);
        b = cpBounds(s+1) - 1;
        segLen(s)  = b - a + 1;
        segMean(s) = (cumY(b+1) - cumY(a)) / segLen(s);
        segVar(s)  = segRSS(a, b, cumY, cumY2) / max(segLen(s) - 1, 1);
    end

    % -- noise estimate --------------------------------------------------
    noiseEst = sqrt(sum(segVar .* segLen) / sum(segLen));
    if noiseEst < 1e-6
        noiseEst = std(y) * 0.1;
    end

    % -- Step 4: agglomerative clustering --------------------------------
    thr = opts.mergeThr * noiseEst;
    stateLabels = agglomerativeCluster(segMean, segLen, thr);

    % -- sort states by ascending mean intensity -------------------------
    uLab  = unique(stateLabels);
    cMean = zeros(size(uLab));
    for j = 1:numel(uLab)
        m = (stateLabels == uLab(j));
        cMean(j) = sum(segMean(m) .* segLen(m)) / sum(segLen(m));
    end
    [~, sOrd] = sort(cMean);
    rankMap = zeros(1, max(uLab));
    for j = 1:numel(uLab)
        rankMap(uLab(sOrd(j))) = j;
    end
    stateLabels = rankMap(stateLabels);

    % -- build frame-by-frame state sequence -----------------------------
    seq = zeros(1, F);
    for s = 1:nSeg
        seq(cpBounds(s):cpBounds(s+1)-1) = stateLabels(s);
    end

    % -- Step 5: bleach detection ---------------------------------------
    seq = mergeCloseStates(seq, y, opts);
    [sFrame, bFrame, seq] = detectBleach(seq, y, opts.bleachTail);

    % -- Step 5b: background verification --------------------------------
    % If all remaining active states have mean intensity below background
    % threshold, the trace is essentially background (e.g. changepoint
    % failed to isolate a sharp spike, producing a single low-level state).
    if any(seq > 0)
        % Use robust (MAD-based) noise estimate for threshold calculation
        dy = diff(y(:)');
        robustNoise = median(abs(dy - median(dy))) / 0.6745 / sqrt(2);
        if robustNoise < 1e-6
            robustNoise = noiseEst;
        end

        % Compute background threshold
        if isempty(opts.bgThreshold)
            sortedY = sort(y);
            nLow = max(round(numel(y) * opts.bgQuantile), 1);
            bgLevel = mean(sortedY(1:nLow));
            bgThreshVerify = bgLevel + opts.bgAutoSigma * robustNoise;
        elseif opts.bgThreshold == 0
            bgThreshVerify = -Inf;  % disabled
        else
            bgThreshVerify = opts.bgThreshold;
        end

        % Check if all active states fall below the background threshold
        uAct = unique(seq(seq > 0));
        allBelowBg = true;
        for j = 1:numel(uAct)
            stateMean = mean(y(seq == uAct(j)));
            if stateMean > bgThreshVerify
                allBelowBg = false;
                break;
            end
        end

        if allBelowBg
            % All "active" states are actually at background level
            seq(:) = 0;
            sFrame = 1;
            bFrame = NaN;
        end
    end

    % -- Step 6: renumber active states 1..S contiguously ----------------
    uAct = unique(seq(seq > 0));
    mapping = zeros(1, max(max(seq), 1));
    for j = 1:numel(uAct)
        mapping(uAct(j)) = j;
    end
    for f = 1:F
        if seq(f) > 0
            seq(f) = mapping(seq(f));
        end
    end
    nAct = numel(uAct);

    % -- Step 7: stateInfo & fitTrace ------------------------------------
    sInfo = struct('label',{}, 'meanIntensity',{}, 'stdIntensity',{}, 'occupancy',{});
    fitTr = zeros(1, F);
    activeFr = F;
    if ~isnan(bFrame)
        activeFr = bFrame - 1;
    end

    for j = 1:nAct
        mask = (seq == j);
        sInfo(j).label         = j;
        sInfo(j).meanIntensity = mean(y(mask));
        sInfo(j).stdIntensity  = std(y(mask));
        sInfo(j).occupancy     = sum(mask) / max(activeFr, 1);
        fitTr(mask)            = sInfo(j).meanIntensity;
    end

    if any(seq == 0)
        fitTr(seq == 0) = mean(y(seq == 0));
    end

    % -- Step 8: transitions ---------------------------------------------
    tr = zeros(0, 2);
    for f = 2:F
        if seq(f) ~= seq(f-1)
            tr(end+1, :) = [f, seq(f)]; %#ok<AGROW>
        end
    end

    % -- determine fitted background ------------------------------------
    if any(seq == 0)
        bgVal = mean(y(seq == 0));
    else
        % if no explicit bleach/background segment exists,
        % use the lowest fitted state's mean as background estimate
        if ~isempty(sInfo)
            bgVal = min([sInfo.meanIntensity]);
        else
            bgVal = 0;
        end
    end

    % -- background-subtracted outputs ----------------------------------
    traceBgSub = y - bgVal;
    fitTraceBgSub = fitTr - bgVal;


    % -- Step 9: store ---------------------------------------------------
    states(n).trace         = y;
    states(n).sequence      = seq;
    states(n).transitions   = tr;
    states(n).nStates       = nAct;
    states(n).stateInfo     = sInfo;
    states(n).fitTrace      = fitTr;
    states(n).background    = bgVal;
    states(n).traceBgSub    = traceBgSub;
    states(n).fitTraceBgSub = fitTraceBgSub;
    states(n).bleachFrame   = bFrame;
    states(n).method        = 'changepoint';

    if isnan(bFrame)
        states(n).lifetime = NaN;
    else
        states(n).lifetime = bFrame - sFrame;
    end

    if mod(n, 10) == 0 || n == N
        fprintf('Progress: %.1f%%\n', (n / N) * 100);
    end
end

end


%% HMM METHOD
function states = runHMM(traceMat, states, opts)

N = size(traceMat, 1);

for n = 1:N
    if opts.verbose
        fprintf('  Trace %d / %d (HMM) ...\n', n, N);
    end

    y = double(traceMat(n, :));
    F = length(y);
    varFloor = max(std(y) * 0.01, 1e-6);

    % ---- BIC model selection ------------------------------------------
    bestBIC = Inf;
    best = [];

    for K = opts.minStates : opts.maxStates
        [mu, sig, A, pi0, LL] = baumWelch(y, K, opts.maxIter, opts.nRestarts, varFloor);
        nPar = 2*K + K*(K-1) + (K-1);
        bic  = -2*LL + opts.bicPenalty * nPar * log(F);
        if bic < bestBIC
            bestBIC = bic;
            best = struct('mu', mu, 'sig', sig, 'A', A, 'pi0', pi0, 'K', K);
        end
    end

    % ---- Viterbi decode -----------------------------------------------
    rawSeq = viterbiDecode(y, best.mu, best.sig, best.A, best.pi0);

    % ---- sort states by ascending mean intensity ----------------------
    [~, ord] = sort(best.mu);
    rank = zeros(1, best.K);
    rank(ord) = 1:best.K;
    seq = rank(rawSeq);

    % ---- bleach detection ---------------------------------------------
    if best.K == 1
        sFrame = 1;
        bFrame = NaN;
    else
        seq = mergeCloseStates(seq, y, opts);
        [sFrame, bFrame, seq] = detectBleach(seq, y, opts.bleachTail);
    end

    % ---- renumber active states 1..S contiguously ---------------------
    uLab = unique(seq(seq > 0));
    mapping = zeros(1, max(max(seq), 1));
    for j = 1:numel(uLab)
        mapping(uLab(j)) = j;
    end
    for f = 1:F
        if seq(f) > 0
            seq(f) = mapping(seq(f));
        end
    end

    % ---- stateInfo & fitTrace -----------------------------------------
    nAct = numel(uLab);
    sInfo = struct('label',{}, 'meanIntensity',{}, 'stdIntensity',{}, 'occupancy',{});
    fitTr = zeros(1, F);
    activeFr = F;
    if ~isnan(bFrame)
        activeFr = bFrame - 1;
    end

    for j = 1:nAct
        mask = (seq == j);
        sInfo(j).label         = j;
        sInfo(j).meanIntensity = mean(y(mask));
        sInfo(j).stdIntensity  = std(y(mask));
        sInfo(j).occupancy     = sum(mask) / max(activeFr, 1);
        fitTr(mask)            = sInfo(j).meanIntensity;
    end

    if ~isnan(bFrame)
        fitTr(seq == 0) = mean(y(seq == 0));
    end

    % ---- transitions --------------------------------------------------
    tr = zeros(0, 2);
    for f = 2:F
        if seq(f) ~= seq(f-1)
            tr(end+1, :) = [f, seq(f)]; %#ok<AGROW>
        end
    end
    % ---- determine fitted background ----------------------------------
    if any(seq == 0)
        bgVal = mean(y(seq == 0));
    else
        % if no explicit bleach/background segment exists,
        % use the lowest fitted state's mean as background estimate
        if ~isempty(sInfo)
            bgVal = min([sInfo.meanIntensity]);
        else
            bgVal = 0;
        end
    end

    % ---- background-subtracted outputs --------------------------------
    traceBgSub = y - bgVal;
    fitTraceBgSub = fitTr - bgVal;

    % ---- store --------------------------------------------------------
    states(n).trace         = y;
    states(n).sequence      = seq;
    states(n).transitions   = tr;
    states(n).nStates       = nAct;
    states(n).stateInfo     = sInfo;
    states(n).fitTrace      = fitTr;
    states(n).background    = bgVal;
    states(n).traceBgSub    = traceBgSub;
    states(n).fitTraceBgSub = fitTraceBgSub;
    states(n).bleachFrame   = bFrame;
    states(n).method        = 'HMM';

    if isnan(bFrame)
        states(n).lifetime = NaN;
    else
        states(n).lifetime = bFrame - sFrame;
    end

    if mod(n, 10) == 0 || n == N
        fprintf('Progress: %.1f%%\n', (n / N) * 100);
    end
end

end


%% BINARY SEGMENTATION
function cps = binarySeg(cumY, cumY2, a, b, penalty, minSeg)
cps = [];
n = b - a + 1;
if n < 2 * minSeg
    return;
end

rss0 = segRSS(a, b, cumY, cumY2);
if rss0 < 1e-20
    return;
end

bestGain = -Inf;
bestT = -1;

for t = a + minSeg - 1 : b - minSeg
    rssL = segRSS(a, t, cumY, cumY2);
    rssR = segRSS(t+1, b, cumY, cumY2);
    rss1 = max(rssL + rssR, 1e-20);
    gain = n * log(rss0 / rss1);
    if gain > bestGain
        bestGain = gain;
        bestT = t;
    end
end

if bestGain > penalty * log(n)
    cpsL = binarySeg(cumY, cumY2, a, bestT, penalty, minSeg);
    cpsR = binarySeg(cumY, cumY2, bestT+1, b, penalty, minSeg);
    cps = [cpsL, bestT+1, cpsR];
end
end


%% SEGMENT RSS
function rss = segRSS(a, b, cumY, cumY2)
n = b - a + 1;
s = cumY(b+1) - cumY(a);
s2 = cumY2(b+1) - cumY2(a);
rss = max(s2 - s*s/n, 0);
end


%% REFINE CHANGEPOINT POSITIONS
function cpBounds = refineCP(cumY, cumY2, cpBounds, minSeg, nIter)
nCP = numel(cpBounds);
if nCP <= 2
    return;
end

SEARCH_HALF = 25;

for iter = 1:nIter
    changed = false;
    for i = 2:nCP-1
        a    = cpBounds(i-1);
        bEnd = cpBounds(i+1) - 1;

        searchLo = a + minSeg;
        searchHi = bEnd - minSeg + 1;
        if searchLo > searchHi
            continue;
        end

        currT = cpBounds(i);
        lo = max(searchLo, currT - SEARCH_HALF);
        hi = min(searchHi, currT + SEARCH_HALF);

        bestRSS = Inf;
        bestT = currT;
        for t = lo:hi
            rL = segRSS(a, t-1, cumY, cumY2);
            rR = segRSS(t, bEnd, cumY, cumY2);
            if rL + rR < bestRSS
                bestRSS = rL + rR;
                bestT = t;
            end
        end

        if bestT ~= cpBounds(i)
            cpBounds(i) = bestT;
            changed = true;
        end
    end
    if ~changed
        break;
    end
end
end


%% AGGLOMERATIVE CLUSTERING
function labels = agglomerativeCluster(segMean, segLen, threshold)
nSeg = numel(segMean);
labels = 1:nSeg;

while true
    uLab   = unique(labels);
    nClust = numel(uLab);
    if nClust <= 1
        break;
    end

    cMean = zeros(1, nClust);
    for j = 1:nClust
        m = (labels == uLab(j));
        cMean(j) = sum(segMean(m) .* segLen(m)) / sum(segLen(m));
    end

    minDist = Inf;
    mergeA = -1;
    mergeB = -1;

    for j = 1:nClust
        for k = j+1:nClust
            d = abs(cMean(j) - cMean(k));
            if d < minDist
                minDist = d;
                mergeA = uLab(j);
                mergeB = uLab(k);
            end
        end
    end

    if minDist > threshold
        break;
    end

    labels(labels == mergeB) = mergeA;
end
end


%% GAUSSIAN HMM — BAUM-WELCH
function [bestMu, bestSig, bestA, bestPi, bestLL] = baumWelch(y, K, maxIter, nRest, varFloor)

F = numel(y);
bestLL = -Inf;
bestMu = [];
bestSig = [];
bestA = [];
bestPi = [];
tol = 1e-6;

for r = 1:nRest
    [mu, sig] = initParams(y, K, r, varFloor);

    if K > 1
        A = (0.02 / (K-1)) * ones(K) + (1 - 0.02) * eye(K);
    else
        A = 1;
    end
    pi0 = ones(1, K) / K;
    prevLL = -Inf;

    for it = 1:maxIter
        B = zeros(K, F);
        for k = 1:K
            B(k,:) = exp(-0.5 * ((y - mu(k)) ./ sig(k)).^2) ./ (sig(k) * sqrt(2*pi));
        end
        B = max(B, 1e-300);

        alpha = zeros(K, F);
        c = zeros(1, F);
        alpha(:,1) = pi0(:) .* B(:,1);
        c(1) = sum(alpha(:,1));
        if c(1) < 1e-300
            c(1) = 1e-300;
        end
        alpha(:,1) = alpha(:,1) / c(1);

        for t = 2:F
            alpha(:,t) = (A' * alpha(:,t-1)) .* B(:,t);
            c(t) = sum(alpha(:,t));
            if c(t) < 1e-300
                c(t) = 1e-300;
            end
            alpha(:,t) = alpha(:,t) / c(t);
        end
        LL = sum(log(c));

        beta = zeros(K, F);
        beta(:,F) = 1;
        for t = F-1:-1:1
            beta(:,t) = A * (B(:,t+1) .* beta(:,t+1));
            beta(:,t) = beta(:,t) / c(t+1);
        end

        gamma = alpha .* beta;
        gamma = gamma ./ (sum(gamma, 1) + 1e-300);

        xiSum = zeros(K);
        for t = 1:F-1
            v = B(:,t+1) .* beta(:,t+1);
            tmp = (alpha(:,t) * v') .* A;
            xiSum = xiSum + tmp / (sum(tmp(:)) + 1e-300);
        end

        pi0 = gamma(:,1)';
        pi0 = max(pi0, 1e-10);
        pi0 = pi0 / sum(pi0);

        if K > 1
            A = xiSum ./ (sum(xiSum, 2) + 1e-300);
            A = max(A, 1e-10);
            A = A ./ sum(A, 2);
        end

        for k = 1:K
            wk = gamma(k,:);
            ws = sum(wk) + 1e-300;
            mu(k) = sum(wk .* y) / ws;
            sig(k) = sqrt(sum(wk .* (y - mu(k)).^2) / ws);
            sig(k) = max(sig(k), varFloor);
        end

        if it > 1 && abs(LL - prevLL) < tol * abs(LL)
            break;
        end
        prevLL = LL;
    end

    if LL > bestLL
        bestLL = LL;
        bestMu = mu;
        bestSig = sig;
        bestA = A;
        bestPi = pi0;
    end
end
end


%% INITIALISE MEANS & STDS
function [mu, sig] = initParams(y, K, r, varFloor)

F = numel(y);
sy = sort(y);
mu = zeros(1, K);
sig = zeros(1, K);

for k = 1:K
    i1 = max(round((k-1)/K * F) + 1, 1);
    i2 = min(round(k/K * F), F);
    mu(k) = mean(sy(i1:i2));
    sig(k) = max(std(sy(i1:i2)), varFloor);
end

if r > 1
    mu = mu + randn(1, K) .* std(y) * 0.15;
    [mu, idx] = sort(mu);
    sig = sig(idx);
end
end


%% VITERBI DECODER
function seq = viterbiDecode(y, mu, sigma, A, pi0)

K = numel(mu);
F = numel(y);
logA = log(A + 1e-300);
logPi = log(pi0(:) + 1e-300);

logB = zeros(K, F);
for k = 1:K
    logB(k,:) = -0.5 * log(2*pi*sigma(k)^2) - 0.5 * ((y - mu(k)) ./ sigma(k)).^2;
end

delta = zeros(K, F);
psi = zeros(K, F);
delta(:,1) = logPi + logB(:,1);

for t = 2:F
    [mx, ix] = max(delta(:,t-1) + logA, [], 1);
    delta(:,t) = mx(:) + logB(:,t);
    psi(:,t) = ix(:);
end

seq = zeros(1, F);
[~, seq(F)] = max(delta(:,F));
for t = F-1:-1:1
    seq(t) = psi(seq(t+1), t+1);
end
end


%% BLEACH DETECTION
function [sFrame, bFrame, seq] = detectBleach(seq, y, minTail)

F = numel(seq);
uStates = unique(seq(seq > 0));

if isempty(uStates)
    sFrame = F + 1;
    bFrame = F + 1;
    seq(:) = 0;
    return;
end

sMean = arrayfun(@(s) mean(y(seq == s)), uStates);
[~, idx] = min(sMean);
bgLabel = uStates(idx);

if numel(uStates) == 1
    sFrame = 1;
    bFrame = NaN;
    return;
end

sFrame = 1;
if seq(1) == bgLabel
    firstActive = find(seq ~= bgLabel, 1, 'first');
    if isempty(firstActive)
        sFrame = F + 1;
        bFrame = F + 1;
        seq(:) = 0;
        return;
    end
    sFrame = firstActive;
    seq(1:sFrame-1) = 0;
end

bFrame = NaN;
lastActive = find(seq ~= bgLabel & seq > 0, 1, 'last');

if isempty(lastActive)
    seq(seq > 0) = 0;
    sFrame = F + 1;
    bFrame = F + 1;
    return;
end

tailLen = F - lastActive;
if tailLen >= minTail
    bFrame = lastActive + 1;
    seq(bFrame:end) = 0;
end
end


%% MERGE CLOSE STATES
function seq = mergeCloseStates(seq, y, opts)

uStates = unique(seq(seq > 0));
if isempty(uStates)
    return;
end

dy = diff(y(:)');
noiseEst = median(abs(dy - median(dy))) / 0.6745 / sqrt(2);
if noiseEst < 1e-6
    noiseEst = std(y) * 0.1;
end

if isempty(opts.bgThreshold)
    sortedY = sort(y);
    nLow = max(round(numel(y) * opts.bgQuantile), 1);
    bgLevel = mean(sortedY(1:nLow));
    bgThresh = bgLevel + opts.bgAutoSigma * noiseEst;
elseif opts.bgThreshold == 0
    bgThresh = -Inf;
else
    bgThresh = opts.bgThreshold;
end

uStates = unique(seq(seq > 0));
sMean = arrayfun(@(s) mean(y(seq == s)), uStates);

bgMask = (sMean <= bgThresh);
if sum(bgMask) > 1
    bgLabels = uStates(bgMask);
    keepLabel = bgLabels(1);
    for j = 2:numel(bgLabels)
        seq(seq == bgLabels(j)) = keepLabel;
    end
end

changed = true;
while changed
    changed = false;
    uStates = unique(seq(seq > 0));
    if numel(uStates) < 2
        break;
    end

    sMean = arrayfun(@(s) mean(y(seq == s)), uStates);
    sLen  = arrayfun(@(s) sum(seq == s), uStates);

    [sMean, sOrd] = sort(sMean);
    uStates = uStates(sOrd);
    sLen = sLen(sOrd);

    diffs = diff(sMean);
    [minD, idx] = min(diffs);

    if minD < opts.minStateSep * noiseEst
        if sLen(idx) >= sLen(idx+1)
            seq(seq == uStates(idx+1)) = uStates(idx);
        else
            seq(seq == uStates(idx)) = uStates(idx+1);
        end
        changed = true;
    end
end
end