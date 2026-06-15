function pointsBlur2D = gen_blurred2D(xSizeUp,zSizeUp,positions,emitterNumber,PSF3d)

% calculete * 2D * image based on these points ! 
pointsBlur2D = zeros(xSizeUp, xSizeUp, 'single');
PSF3d = PSF3d/sum(PSF3d(:,:,10), "all");

for i = 1:emitterNumber
    % 错位量表示和亮度
    x = positions(i,1);
    y = positions(i,2);
    z = positions(i,3);
    N = positions(i,4); 

    x_start = max(1, x+1);
    x_end   = min(xSizeUp, x + xSizeUp);
    y_start = max(1, y+1);
    y_end   = min(xSizeUp, y + xSizeUp);
    z_index = zSizeUp/2 + 1 - z;  % 取近似 Z 轴平面

    % 防止 z_index 超界
    z_index = max(1, min(z_index, zSizeUp));

    % 计算 PSF 选取区域
    psf_x_start = 1 + (x_start - (x+1));
    psf_x_end   = xSizeUp - ((x + xSizeUp) - x_end);
    psf_y_start = 1 + (y_start - (y+1));
    psf_y_end   = xSizeUp - ((y + xSizeUp) - y_end);

    % 累加 PSF 到 2D 平面
    pointsBlur2D(x_start:x_end, y_start:y_end) = ...
        pointsBlur2D(x_start:x_end, y_start:y_end) + ...
        PSF3d(psf_x_start:psf_x_end, psf_y_start:psf_y_end, z_index) * N;
end
