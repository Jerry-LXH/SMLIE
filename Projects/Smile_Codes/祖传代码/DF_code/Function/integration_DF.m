function ginte=integration_DF(fr_D)
n=fr_D;
% 为每个参数赋值
[u1,u2,ss1,ss2,A,th]=deal(n(1),n(2),n(3),n(4),n(5),n(6));
cc=[cos(th),sin(th);-sin(th),cos(th)];
xy=[u1;u2];
xxyy=cc^(-1)*xy;
[u1,u2]=deal(xxyy(1),xxyy(2));
f=@(x,y)A*exp(-(x-u1).^2./(2*ss1^2)-(y-u2).^2/(2*ss2^2))/(2*pi*ss1*ss2);
ginte=integral2(f,u1-10*ss1,u1+10*ss1,u2-10*ss2,u2+10*ss2);

%% 坐标不对，转换后的坐标