%% Parameter settings
    clear all;clc;

    file_name = '/Volumes/SMILeSSD/Optics/Ti2Tests/TirfModule/Cy3_polylys/20260609/cy3poly_1.80A_0p1s_90x/20260609-225638956/TUC-001.tif';
    %file_name = '/Volumes/SMILeSSD/Optics/TE2000UTests/Cy3_polylys/20260606_cmos-1.4mW/20260606-122652746/TUC-001.tif';
    %file_name_filt = '/Volumes/SMILeSSD/Optics/Ti2Tests/TirfModule/Cy3_polylys/20260227_532_cy3_0.2s_1.82A_filter/TUC-002.tif';
    gain = 33.448; % DN/e-
    ex_time = 0.1;
    offset = 1893.68;
    interval = 0; % time interval between two adjacent exposure
    oneFrameTime = ex_time+interval;
    
%% read folder into tiff
    folderPath = '/Volumes/SMILeSSD/Optics/Ti2Tests/TirfModule/Cy3_polylys/20260206_532_cy3Test_0.2s_1.72A/2_filt';
    files = dir(fullfile(folderPath, '*.tif'));
    [~,idx] = sort({files.name});
    files = files(idx);
    % 读第一张确定尺寸
    img = imread(fullfile(folderPath,files(1).name));
    [ny,nx] = size(img);

    n = numel(files);
    raw_data = zeros(ny,nx,n,class(img));

    for i = 1:n
        raw_data(:,:,i) = imread(fullfile(folderPath,files(i).name));
    end

%% write tiff
    outFile = [file_name(1:end-4) '.tif'];
    io.writeTiffStack(raw_data, outFile);

%% read tif stack
    raw_data = io.readTiffStack(file_name);
    % raw_data_filt = io.readTiffStack(file_name_filt);

%% Convert to electrons unit
    raw_data = double(raw_data) - offset;
    raw_data = raw_data / gain;
    %raw_data_filt = double(raw_data_filt) - offset;
    %raw_data_filt = raw_data_filt / gain;

%% Load sif data (Note that data is in absolute e- unit)
    [raw_data,ex_time,gainDAC] = io.readSIFData(file_name);
    interval = 0; % time interval between two adjacent exposure. In unit of seconds.
    oneFrameTime = ex_time+interval;

%% Averaging data
    mean_raw_data = mean(raw_data(:,:,1:20),3);
    %mean_raw_data_filt = mean(raw_data_filt(:,:,:),3);

%% clear data
    clear raw_data

%% Viz
    figure;
    viz.plotImage(mean_raw_data, 1, 'hot','Mean Image (without filter)');
    %figure;
    %viz.plotImage(mean_raw_data_filt, 1, 'hot','Mean Image (with filtert)');

%% Freq Domain
    F = fftshift(abs(fft2(mean_raw_data)));
    figure;
    viz.plotImage(log(F+1), 1, 'gray','Freq Domain');

%% Lowpass Filter
    L_est = imgaussfilt(mean_raw_data, 60);
    %L_est_filt = imgaussfilt(mean_raw_data_filt, 50);
    figure;
    viz.plotImage(L_est, 1, 'hot','Lowpassed Image (without filter)');
    %figure;
    %viz.plotImage(L_est_filt, 1, 'hot','Lowpassed Image (with filtert)');
    %figure;
    %viz.plotImage(log10(L_est./L_est_filt), 1, 'hot','Lowpassed Image (ratio)');

%% Gaussian Fit
    est = optics.illumi.fitIlluminationGaussian(L_est);
    
%% Fit model
    [ny,nx] = size(L_est);
    [xx,yy] = meshgrid(1:nx,1:ny);

    A  = est.A;
    x0 = est.x0;
    y0 = est.y0;
    sx = est.sx;
    sy = est.sy;
    B  = est.B;

    fit_img = A * exp(-((xx-x0).^2)/(2*sx^2) ...
                 -((yy-y0).^2)/(2*sy^2)) + B;
    figure;
    viz.plotImage(fit_img, 1, 'hot','Fit Image');

%% Residual
    residual_map = L_est - fit_img;
    figure;
    viz.plotImage(residual_map, 1, 'hot','Residual Image');


