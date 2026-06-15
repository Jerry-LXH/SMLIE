%%
x=[0.025,0.035,0.050,0.060,0.080];
y=[20,30,40,50,60];
fit=polyfit(x,y,1);
figure
xfit=x(1):0.0001:x(end);
yfit=fit(1)*xfit+fit(2);
plot(x,y,'ro',xfit,yfit);
grid on


%%
clear;clc
% 数值微分
h=1;
for i=1:6
x=pi/2;h=h/10;
xd=[x,x+h];y=[sin(x),sin(x+h)];
slope=diff(y)./diff(xd);
end

%%
x=0:0.01:2*pi;
y=sin(x);
slope=diff(y)./diff(x);
plot(x(1:end-1),slope,x,y)

%%
h=1;hold on;g=colormap(lines)
for i=1:3
    h=h/10;
    x=1:h:2*pi;
    y=sin((x.^2)/2).*exp(-x);
    slope=diff(y)./diff(x);
    plot(x(1:end-1),slope,'Color',g(i,:));
end
hold off;
set(gca,'Xlim',[0,2*pi]);

%% 数值积分
h=0.01;x=0:h:2;
y=4*x.^3;s=h*trapz(y)


%%
f=@(x,y) y.*sin(x)+x.*cos(y);
integral2(f,pi,2*pi,0,pi)


