function [J,res] = jacobian(theta, roi, xx, yy, F_EM, use_pixel)

eps_val = 1e-6;

res0 = localize.LSQ1.residual(theta, roi, xx, yy);
m = numel(res0);
n = numel(theta);

J = zeros(m,n);

for k = 1:n
    dtheta = zeros(size(theta));
    dtheta(k) = eps_val;
    res1 = localize.LSQ1.residual(theta + dtheta, roi, xx, yy, F_EM, use_pixel);
    J(:,k) = (res1 - res0) / eps_val;
end

res = res0;

end