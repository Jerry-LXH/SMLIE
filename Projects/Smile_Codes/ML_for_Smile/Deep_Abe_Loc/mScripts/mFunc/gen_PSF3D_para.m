function PSF3d = gen_PSF3D_para(xPixelSize, xSize, zPixelSize, zSize, lambda, RI, pupilFun, pupilMask)
    % 生成 3D PSF，基于 pupil function
    
    % 计算焦外波前相位
    dConst = defocus_wavefront(xPixelSize, xSize, lambda, RI, pupilMask);
    PSF3d = zeros(xSize, xSize, zSize, 'single'); % 预分配内存
    Soz = (zSize+1)/2;
    
    % **使用 GPU 加速计算**
    useGPU = (gpuDeviceCount > 0);
    if useGPU
        pupilFun = gpuArray(pupilFun);
        dConst = gpuArray(dConst);
        PSF3d = gpuArray(PSF3d);
    end
    
    % **预计算相移矩阵**
    zPositions = ((1:zSize) - Soz) * zPixelSize;
    phaseShifts = exp(1i * dConst .* reshape(zPositions, 1, 1, []));
    
    % **并行计算 PSF**
    parfor i = 1:zSize
        pupilFun2D = pupilFun .* phaseShifts(:,:,i); % 直接使用预计算相移矩阵
        prf = fftshift(ifft2(ifftshift(pupilFun2D))); % 计算光场
        PSF3d(:,:,i) = abs(prf).^2; % 计算强度
    end
    
    % **归一化**
    PSF3d = PSF3d / max(PSF3d(:));
    
    % **从 GPU 取回数据**
    if useGPU
        PSF3d = gather(PSF3d);
    end
    
    end