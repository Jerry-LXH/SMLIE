%%
%基础绘图
%%
%plot from data
%matlab will refresh:hold on~hold off
%调格式 linespec
%legend：命名每个曲线
x=0:0.5:4*pi;
y=sin(x);h=cos(x);
plot(x,y,'r',x,h,'g');
legend('sin(x)','cos(x)');%依次命名每个曲线
title('a')
xlabel('t=0 to 2\pi')
ylabel('e^{-x}')
%zlabel()三维图形
str='$$ \int_{0}^{2} x^2\sin(x) dx $$'
text(0.25,0.5,str,'Interpreter','latex')%text()可以用LaTeX来输入特殊符号/公式
annotation('arrow',[0,1],[0,1]);

%%
x=1:0.01:2;
y=x.^2;h=sin(2.*pi.*x);
hold on
plot(x,y,'k',x,h,'ro');
hold off
title('Mini Assigenment');
xlabel('Time (ms)');
ylabel('f(t)');
legend('t^{2}','sin(2\pit)','Location','northwest');

%% figure adjustment
% 绘图物件
%figure object

%axes object
%坐标
%line object
%线
%怎么调节图的属性
%1)找到不同元素的识别码
%gcf,gcf,allchild,ancestor,delete,findall
%2）get(),set()



%%
x=linspace(0,2*pi,1000);
y=sin(x);plot(x,y)
h=plot(x,y);
% set(gca,'XLim',[0,2*pi]);
% set(gca,'ylim',[0,2]);
% set(gca,'fontsize',25);
set(gca,'XTick',0:pi/2:2*pi);
% set(gca,'xticklabel',0:90:360);
set(gca,'fontname','latex');
set(gca,'xticklabel',{'0','p/2','p','3p/2','2p'});

%% line
set(plot(x,y),'LineStyle','-.','ColorMode','manual','Color','g','LineWidth',7);

%marker/每个点：color/face


%%
x=rand(20,1);set(gca,'fontsize',18);
plot(x,'-md','markeredgecolor','k','markerfacecolor','g','markersize',10);
xlim([1,20]);

%%
x=1:0.01:2;
y=x.^2;h=sin(2.*pi.*x);
hold on
h=plot(x,h,'ro',x,y,'k');
hold off
title('Mini Assigenment');
xlabel('Time (ms)');
ylabel('f(t)');
legend('t^{2}','sin(2\pit)','Location','northwest');
set(gca,'fontsize',15);
set(h,"LineWidth",2);
set(h,"MarkerFaceColor",[0.1,0.2,0.3]);

%% 多图
%figure
%注意gca和gcf（指向第二个刚画的）
x=-10:0.1:10;
y=x.^2-8;
h=exp(x);
figure,plot(x,y);
figure,plot(x,h);

%%
%图的大小和位置
%figure('Position',[left,buttom,width,height]);

%%
%一个界面多个图
t=0:0.1:2*pi;x=3*cos(t);y=sin(t);
subplot(2,2,1);plot(x,y);axis normal;
subplot(2,2,2);plot(x,y);axis square;
subplot(2,2,3);plot(x,y);axis equal;
subplot(2,2,4);plot(x,y);axis equal tight;
%保存图片
saveas(gcf,'exerice','png');
%bitmap
%vector矢量图
%解析度很高，print()









