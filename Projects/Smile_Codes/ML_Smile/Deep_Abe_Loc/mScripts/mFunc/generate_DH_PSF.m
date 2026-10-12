function phase_mask = generate_DH_PSF(r, theta, idx, p_values, l_values, w, num_iter)
    % 生成 Double Helix PSF (DH-PSF) 相位掩膜
    %
    % 输入参数：
    %   r        - 极坐标中的径向坐标 (矩阵)
    %   theta    - 极坐标中的角度坐标 (矩阵)
    %   p_values - Gauss-Laguerre 模态的径向模式索引 (向量)
    %   l_values - Gauss-Laguerre 模态的角动量索引 (向量)
    %   w        - 光束腰宽 (标量)
    %   num_iter - 迭代优化次数 (标量)
    %
    % 输出：
    %   phase_mask - 生成的 DH-PSF 相位掩膜 (矩阵)
    
    % 初始化 GL 模态叠加
    GL_sum = zeros(size(r));
    
    % 遍历所有指定的 GL 模态
    for k = 1:length(p_values)
        p = p_values(k); % 分别取出第k个标签
        l = l_values(k);
        
        % 计算 Laguerre 多项式
        L_pl = laguerreL(p, abs(l), 2*r.^2/w^2);
        
        % 计算 GL 模态
        GL_mode = ((r./w).^ abs(l)) .* exp(-r.^2/w^2) .* L_pl .* exp(-1i * l * theta);
       
        % 按能量归一化
        GL_mode = GL_mode ./ sqrt(sum(abs(GL_mode(:)).^2));
        
        % 叠加所有模式
        GL_sum = GL_sum + GL_mode;
        fprintf('Mode (p=%d, l=%d): Energy(GL_mode) = %.4f, max(GL_mode) = %.4f\n', ...
        p, l, sum(abs(GL_mode(:)).^2), max(abs(GL_mode(:))));
    end
    

    % 归一化叠加后的 GL 模态
    GL_sum = GL_sum ./ max(abs(GL_sum(:)));

    % 进行频域优化
    for iter = 1:num_iter
        % 计算傅里叶变换
        GL_freq = fftshift(fft2(fftshift(GL_sum)));

        % 施加频域约束（例如仅保留某些频率范围的分量）
        GL_freq(~idx) = 0;

        % 逆傅里叶变换回到空间域
        GL_sum = ifftshift(ifft2(ifftshift(GL_freq)));

        fprintf('Iteration %d: max(GL_sum) = %.4f\n', iter, max(abs(GL_sum(:))));

        % 归一化
        GL_sum = GL_sum ./ max(abs(GL_sum(:)));

    end

    % 计算最终相位掩膜
    phase_mask = angle(GL_sum);

    % 施加 2π 取模（确保相位值在 [0, 2π] 范围内）
    phase_mask = mod(phase_mask, 2 * pi);

    % phase_mask = exp(1i*phase_mask);
    phase_mask = GL_sum;

    end