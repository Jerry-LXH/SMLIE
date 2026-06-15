% %%
% %function
% function x=freebody(x0,v0,t)
% %自由落体运动
% x=x0+v0.*t+1/2*9.8*t.*t;
% %.*:元素与元素直接相乘
% function [a,F]=acc(v2,v1,t2,t1,m)
% a=(v2-v1)./(t2-t1);
% F=m.*a;


function C=FtoC (F)
disp(' Temperature in F is: ');