function pointsBlur2D = gen_blurred2D_gpu(xSizeUp, zSizeUp, positions, emitterNumber, PSF3d)
    % Calculate 2D image based on these points
    pointsBlur2D = zeros(xSizeUp, xSizeUp, 'single');
    
    % Normalize PSF3d once outside the loop
    PSF3d = PSF3d / sum(PSF3d(:,:,10), "all");

    % Calculate center of z slice
    zCenter = zSizeUp / 2;
    
    % Preallocate a temporary array for each iteration result
    tempBlurred = zeros(xSizeUp, xSizeUp, 'single');
    
    parfor i = 1:emitterNumber
        % Get the coordinates and intensity for this emitter
        x = positions(i, 1);
        y = positions(i, 2);
        z = positions(i, 3);
        N = positions(i, 4);

        % Calculate start and end indices for x and y
        x_start = max(1, x + 1);
        x_end   = min(xSizeUp, x + xSizeUp);
        y_start = max(1, y + 1);
        y_end   = min(xSizeUp, y + xSizeUp);
        
        % Calculate the z index, adjusting to the center slice
        z_index = max(1, min(zCenter + 1 - z, zSizeUp)); 

        % Calculate PSF selection region
        psf_x_start = max(1, x_start - x + 1);  % Ensure starting indices are within bounds
        psf_x_end   = min(xSizeUp, x_end - x + 1);
        psf_y_start = max(1, y_start - y + 1);  
        psf_y_end   = min(xSizeUp, y_end - y + 1);

        % Accumulate PSF into a temporary 2D plane (this avoids issues with parallel indexing)
        tempBlurred(x_start:x_end, y_start:y_end) = ...
            tempBlurred(x_start:x_end, y_start:y_end) + ...
            PSF3d(psf_x_start:psf_x_end, psf_y_start:psf_y_end, z_index) * N;
    end
    
    % Add the results from each parallel iteration to the final 2D image
    pointsBlur2D = pointsBlur2D + tempBlurred;
end