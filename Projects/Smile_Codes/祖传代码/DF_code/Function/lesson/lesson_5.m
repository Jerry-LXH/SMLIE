%% advanced plot
% 2D plot/对数图
x=logspace(-1,1,100);
y=x.^2;
subplot(2,2,1);
plot(x,y);
subplot(2,2,2);
semilogx(x,y);
subplot(2,2,3);
semilogy(x,y);
subplot(2,2,4);
loglog(x,y);
set(gcf,'xgrid','on')

%% 双Y轴
x=0:0.01:20;
y1=200*exp(-0.05*x).*sin(x);
y2=0.8*exp(-0.5*x).*sin(10*x);
[AX,H1,H2]=plotyy(x,y1,x,y2);
set(get(AX(1),'ylabel'),'string','L-Ylabel');
set(get(AX(2),'ylabel'),'string','R-Ylabel');

%% 直方图histogram
y=randn(1,1000);
hist(y,50);

%% bar charts
x=[1 2 5 4 8];y=[x;1:5];
subplot(1,3,1);bar(x);title("xxx")
subplot(1,3,2);bar(y);
subplot(1,3,3);bar3(y);
% subplot(1,2,1);bar(y,'stacked');
% subplot(1,2,2);barh(y);

%% pie char
a=[10 5 20 30];
subplot(1,3,1);pie(a);
subplot(1,3,2);pie(a,[0,0,1,0]);
subplot(1,3,3);pie3(a,[0,0,0,1]);

%% polar chart
x=1:100;theta=x/10;r=log10(x);
polarplot(theta,r);

%% stairs and stem charts
x=linspace(0,4*pi,40);y=sin(x);
subplot(1,2,1);stairs(y);
subplot(1,2,2);stem(y);

%%
t1=0:0.01:10;y=sin((pi*t1.^2)/4);
t2=0:0.2:10;
hold on
plot(t1,y);
stem(t2,sin((pi*t2.^2)/4));
set(gca,'ylim',[-1,1]);
yticks([-1 -0.5 0 0.5 1]);
hold off

%% Boxplot and error bar
load carsmall
boxplot(MPG,Origin)

x=0:pi/10:pi;y=sin(x);
e=std(y)*ones(size(x));
errorbar(x,y,e);


%% fill

t=(1:2:15)*pi/8;x=sin(t);y=cos(t);
fill(x,y,'r'); 
axis square off
text(0,0,'STOP','color','w','FontSize', ...
    80,'FontWeight','bold','HorizontalAlignment','center');

%%
t=(0:3)*pi/2;x=sin(t);y=cos(t);
p=fill(x,y,'y'); 
axis square off
text(0,0,'WAIT','color','k','FontSize', ...
    80,'FontWeight','bold','HorizontalAlignment','center');
p.LineWidth=4;

%% color space
%RGB [R G B],0-255
% x=randi([0,255])   %随机数
G=[46 38 29 24 13];S=[29 27 17 26 23];B=[29 23 19 32 7];
h=bar(1:5,[G' S' B']);
a=[G' S' B'];
h(1).FaceColor='#FFD700';
h(2).FaceColor='#C0C0C0';
h(3).FaceColor='#FF9912';

%% imagesc()
clear;
[x,y]=meshgrid(-3:.2:3,-3:.2:3)
z=10000*(x.^2+x.*y+y.^2);
surf(x,y,x.^2+x.*y+y.^2);box on;
set(gca,'fontsize',16);
figure;
imagesc(z);axis square;
colorbar;
colormap("gray");
%% colormap
x=[1:10;3:12;5:14];
imagesc(x);
colorbar;
a=zeros(256,3);
for i=1:256
    a(i,:)=[0 i/256 0];
end
colormap(a);


%% 3D plot
x=0:0.1:3*pi;
z1=sin(x);z2=sin(2*x);z3=sin(3*x);
y1=zeros(size(x));y3=ones(size(x));y2=y3./2;
plot3(x,y1,z1,'r',x,y2,z2,'b',x,y3,z3,'g');grid on;
xlabel("x");ylabel("y");zlabel("z");


%% 3D surface plot
%需要先用meshgrid 创建一个网格，
x=-2:0.1:2;
y=-2:0.1:2;
[X,Y]=meshgrid(x,y);
z=X.*exp(-X.^2-Y.^2);
subplot(1,2,1);mesh(X,Y,z);
subplot(1,2,2);surf(X,Y,z);
colorbar;
figure;
% contour
%
subplot(1,3,1);contour(z,[-0.45:0.05:0.45]);axis square
subplot(1,3,2);[C,h]=contour(z);clabel(C,h);axis square;

%%
imgGrayData(rows,cols)

%%
hold on
    [C,h]=contourf(z,[-0.45:0.05:0.45]);axis square;
    clabel(C,h);
hold off

%% meshc/surfc
meshc(X,Y,z);

%% view() 看的角度
sphere(200);shading("flat");
light('position',[1 3 2]);
light('position',[-3 -1 3]);
material shiny;
axis vis3d off;
set(gcf,'color',[1 1 1]);
view(-45,20);
% light 

















