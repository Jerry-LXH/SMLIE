n1 = 1.515;            % 玻璃/物镜
n2 = 1.00;             % 样品 (水)
lambda0 = 0.532;       % um, Cy3 发射可用激发波长替代
NA = 1.49;             % 物镜 NA

n  = n2/n1;
thc = asind(n);                          % 临界角 (deg)
th_max = asind(min(NA/n1, 1));           % NA 限定的最大角
fprintf('theta_c = %.2f deg, theta_max(NA) = %.2f deg\n', thc, th_max);

th = linspace(thc+1e-3, th_max, 500);    % TIR 区间
t  = th*pi/180;

Is = 4*cos(t).^2 ./ (1 - n^2);
Ip = 4*cos(t).^2 .* (2*sin(t).^2 - n^2) ./ ...
     (n^4*cos(t).^2 + sin(t).^2 - n^2);

d  = lambda0/(4*pi) ./ sqrt(n1^2*sin(t).^2 - n2^2);   % 穿透深度(um)

figure;
subplot(2,1,1); plot(th,Ip,'-',th,Is,':','LineWidth',1.5);
xlabel('\theta (deg)'); ylabel('I(0)/I_{inc}');
legend('p','s'); xline(thc,'--'); grid on; title('界面倏逝强度');
subplot(2,1,2); plot(th,d*1000,'LineWidth',1.5);
xlabel('\theta (deg)'); ylabel('penetration depth (nm)'); grid on;


%% ===== 表观增强 vs 角度（修正版）=====
f_p     = 0.5;
sample  = 'beads';   % 'surface' | 'uniform' | 'layer' | 'shell'
h       = 0.5;         % um, 'layer'
z0      = 0.05;        % um, 'shell'
L       = 2.0;         % um, epi 等效采集深度（仅 ref='epi' 用）
R = 0.05;   % um, bead 半径 (100nm bead → 0.05)
ref     = 'epi';       % 'abs' = 绝对 I(0)/I_inc ; 'epi' = 相对落射

I0    = f_p*Ip + (1-f_p)*Is;            % p/s 加权界面强度

switch sample            % 轴向积分 ∫ρ e^{-z/d} dz  (单位 um, surface 无量纲)
     case 'surface', axial = ones(size(d)); Sepi = 1;
     case 'uniform', axial = d;             Sepi = L;
     case 'layer',   axial = d.*(1-exp(-h./d)); Sepi = h;
     case 'shell',   axial = exp(-z0./d);   Sepi = 1;   % 薄壳
     case 'sphere'                                              % 实心 bead
          % ∫(2Rz - z²)e^{-z/d}dz / (4R³/3)，逐角度数值积分
          axial = zeros(size(d));
          for k = 1:numel(d)
               zz = linspace(0, 2*R, 400);
               num = trapz(zz, (2*R*zz - zz.^2).*exp(-zz/d(k)));
               axial(k) = num / (4*R^3/3);
          end
          Sepi = 1;            % epi 下整球均匀照亮，归一化已含在分母
end

switch ref
    case 'abs',  S = I0 .* axial;          ylab = 'I(0)/I_{inc} \times 轴向因子';
    case 'epi',  S = I0 .* axial / Sepi;   ylab = 'TIRF/epi 增强倍数';
end

figure; plot(th, S, 'LineWidth', 2); hold on; xline(thc,'--');
xlabel('\theta (deg)'); ylabel(ylab);
title(sprintf('%s, f_p=%.2f, ref=%s', sample, f_p, ref)); grid on;

[Smax, im] = max(S);
fprintf('峰值 = %.2f, 峰位 = %.2f deg\n', Smax, th(im));