%选择光信号矩阵中的某一帧，选择其中某一个分子，选中其周围的9*9矩阵。计算背景亮度，即最后50帧亮度取平均。计算这个9*9矩阵在第i帧的亮度；这个亮度=（9*9矩阵亮度之和-背景亮度）/曝光时间。作图，纵轴是亮度；横轴是时间，=帧数*曝光时间
function gettrace(s,time,px,py,fra,meas,Loctext)
if nargin==5
  meas=7;%输入想要分析几*几矩阵，比如输入9代表分析9*9矩阵。必须是奇数。
end
[r,ss,t]=size(s);

mbg_br=mean(sum(sum(s((py-(meas-1)/2):(py+(meas-1)/2),(px-(meas-1)/2):(px+(meas-1)/2),(end-50):end))));

br=(sum(sum(s((py-(meas-1)/2):(py+(meas-1)/2),(px-(meas-1)/2):(px+(meas-1)/2),fra:t)))-mbg_br)/time;
bri=zeros(1,t-fra+1);
for j=1:(t-fra+1)
    bri(j)=br(j);
end

time_mat=(1:(t-fra+1))*time;

figure
set(gcf,'position',[600 600 600 500])


plot(time_mat,bri,'LineWidth',1)
if Loctext==1
    title([num2str(px),' ',num2str(py)])
end
xlabel('time(s)')
ylabel('Photon luminescence (photons/s)')
ylim([min(bri),max(bri)*1.15])
xlim([0,80])
%title([num2str(px),' ',num2str(py)])
set(gca,'FontSize',18)
end