function [data, exTime, gainDAC] = readSIFData(fileName, interval, convFactor)
% readSIFData  Read SIF file and convert to photons
% [interval]: read one frame per interval frames. If interval=1, read all frames.
% [convFacctor]: convert values into electrons.

    if nargin < 2 || isempty(interval)
        interval = 1;
    end
    if nargin < 3 || isempty(convFactor)
        convFactor = 5.75;
    end

    dataSet = io.readsifX(fileName,interval);

    exTime  = dataSet.exposureTime;
    gainDAC = dataSet.gainDAC;

    rawData = dataSet.imageData;

    data = rawData .* convFactor ./ gainDAC;
    
end

