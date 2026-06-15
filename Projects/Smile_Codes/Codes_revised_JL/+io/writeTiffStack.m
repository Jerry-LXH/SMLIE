function writeTiffStack(data3D, outFile)
% writeTiffStack  Save 3D stack using Tiff class
%   Preserves original numeric type

    if isfile(outFile)
        warning('writeTiffStack:Overwrite', ...
            'File "%s" already exists and will be overwritten.', outFile);
    end

    [height, width, nFrames] = size(data3D);
    dataType = class(data3D);
    
    t = Tiff(outFile, 'w');

    for k = 1:nFrames

        tagStruct.ImageLength      = height;
        tagStruct.ImageWidth       = width;
        tagStruct.Photometric      = Tiff.Photometric.MinIsBlack;
        tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

        switch class(data3D)
            case 'uint16'
                tagStruct.BitsPerSample = 16;
                tagStruct.SampleFormat  = Tiff.SampleFormat.UInt;
            case 'uint8'
                tagStruct.BitsPerSample = 8;
                tagStruct.SampleFormat  = Tiff.SampleFormat.UInt;
            case 'single'
                tagStruct.BitsPerSample = 32;
                tagStruct.SampleFormat  = Tiff.SampleFormat.IEEEFP;
            case 'double'
                tagStruct.BitsPerSample = 64;
                tagStruct.SampleFormat  = Tiff.SampleFormat.IEEEFP;
            otherwise
                error('Unsupported data type: %s', class(data3D));
        end

        tagStruct.SamplesPerPixel = 1;
        tagStruct.RowsPerStrip    = height;
        tagStruct.Software        = 'MATLAB';

        t.setTag(tagStruct);
        t.write(data3D(:,:,k));

        if k < nFrames
            t.writeDirectory();
        end
    end

    t.close();

    fprintf('TIFF written: %s\n', outFile);
    fprintf('Size: %d x %d x %d | Type: %s\n', ...
        height, width, nFrames, dataType);
end