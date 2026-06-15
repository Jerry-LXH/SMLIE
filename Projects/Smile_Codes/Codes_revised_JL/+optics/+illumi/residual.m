function r = residual(p, L, xx, yy)

A  = p(1);
x0 = p(2);
y0 = p(3);
sx = abs(p(4));
sy = abs(p(5));
B  = p(6);

model = A * exp(-((xx-x0).^2)/(2*sx^2) ...
               -((yy-y0).^2)/(2*sy^2)) + B;

r = model - L;
r = r(:);

end