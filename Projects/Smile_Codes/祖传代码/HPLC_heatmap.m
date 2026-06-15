clc;
clear all;
raw_data = readmatrix('/Users/lu_yk/Desktop/Data/HPLC/RNO-0-release/RNO-0-open.xlsx');
%%
data = raw_data(2:end,3:end-50)';
x = raw_data(2:end,1)';
y = raw_data(1,3:end-50);
clim([0 30]);
imagesc(x,y,data,clim);
set(gca,'fontname','times new roman','fontsize',20,'xtick',(7:9));
axis square;
colorbar('ticks',[]);