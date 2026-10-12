function data = readTiffStack(inFile)
% readTiffStack  Read multi-page TIFF using Tiff class
%
%   data = readTiffStack(inFile)
%
%   Output preserves original numeric type.

    t = Tiff(inFile, 'r');

    % 读取第一页信息
    width  = t.getTag('ImageWidth');
    height = t.getTag('ImageLength');
    bits   = t.getTag('BitsPerSample');
    sampleFormat = t.getTag('SampleFormat');

    % 判断数据类型
    switch sampleFormat
        case Tiff.SampleFormat.UInt
            if bits == 8
                dataType = 'uint8';
            elseif bits == 16
                dataType = 'uint16';
            else
                error('Unsupported unsigned integer bit depth: %d', bits);
            end
        case Tiff.SampleFormat.IEEEFP
            if bits == 32
                dataType = 'single';
            elseif bits == 64
                dataType = 'double';
            else
                error('Unsupported floating-point bit depth: %d', bits);
            end
        otherwise
            error('Unsupported SampleFormat.');
    end

    % 统计帧数
    nFrames = 1;
    while ~t.lastDirectory()
        t.nextDirectory();
        nFrames = nFrames + 1;
    end

    % 预分配
    data = zeros(height, width, nFrames, dataType);

    % 回到第一页
    t.setDirectory(1);

    % 逐帧读取
    for k = 1:nFrames
        data(:,:,k) = t.read();
        if ~t.lastDirectory()
            t.nextDirectory();
        end
    end

    t.close();

    fprintf('TIFF read: %s\n', inFile);
    fprintf('Size: %d x %d x %d | Type: %s\n', ...
        height, width, nFrames, dataType);
end