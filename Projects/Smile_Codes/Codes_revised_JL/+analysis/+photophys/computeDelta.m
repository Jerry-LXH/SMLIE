function [stateMax, stateMin, stateMean] = computeDelta(states, minOccurrences)
%COMPUTEDELTA Compute max/min/weighted-mean state intensities for each molecule.
%
% [stateMax, stateMin, stateMean] = computeDelta(states, minOccurrences)
%
% INPUT
%   states          struct array from analyzeStates
%   minOccurrences  minimum number of contiguous runs for a state (default 3)
%
% OUTPUT
%   stateMax        N x 1, mean intensity of brightest valid state
%   stateMin        N x 1, mean intensity of darkest valid state
%   stateMean       N x 1, occupation-weighted mean intensity of valid states
%
% Notes:
%   - state 0 is excluded
%   - if fewer than 2 valid states are found, outputs are NaN

    if nargin < 2 || isempty(minOccurrences)
        minOccurrences = 3;
    end

    N = numel(states);
    stateMax  = nan(N, 1);
    stateMin  = nan(N, 1);
    stateMean = nan(N, 1);

    for n = 1:N
        seq = states(n).sequence(:)';

        if isempty(seq)
            continue;
        end

        % exclude state 0
        uStates = unique(seq(seq > 0));

        if isempty(uStates)
            continue;
        end

        % count number of contiguous appearances for each state
        occCount = zeros(size(uStates));
        for i = 1:numel(uStates)
            occCount(i) = countStateRuns(seq, uStates(i));
        end

        % keep states appearing at least minOccurrences times
        validStates = uStates(occCount >= minOccurrences);

        if numel(validStates) < 2
            continue;
        end

        % mean intensity and occupation of valid states
        validMeans = nan(size(validStates));
        validOccs  = zeros(size(validStates));
        for i = 1:numel(validStates)
            s = validStates(i);
            mask = (seq == s);
            validMeans(i) = mean(states(n).trace(mask), 'omitnan');
            validOccs(i)  = sum(mask);
        end

        stateMax(n)  = max(validMeans);
        stateMin(n)  = min(validMeans);
        stateMean(n) = sum(validMeans .* validOccs) / sum(validOccs);
    end

end


function nRuns = countStateRuns(seq, stateLabel)
%COUNTSTATERUNS Count number of contiguous runs for a given state.

    mask = (seq == stateLabel);
    d = diff([0, mask, 0]);
    nRuns = sum(d == 1);
end