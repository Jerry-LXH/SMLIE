function [xr,yr]=gfcenter(x,y,th)

cc=[cos(th),sin(th);-sin(th),cos(th)];
xy=[x;y];
xxyy=cc^(-1)*xy;
xr=xxyy(1);yr=xxyy(2);