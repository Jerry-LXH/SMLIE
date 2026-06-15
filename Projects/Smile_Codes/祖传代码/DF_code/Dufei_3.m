%%
close all
clear;
% time为模拟的时间次数，0.1为一个间隔
time=100;
% 进行100次模拟，每次模拟10S，每个数据代表0.1S
x=zeros(100,time);
% 设置初始条件，每次模拟一A状态开始
x(:,1)=1;
% pa为k(A->B,1)*t(0.1S)，同理k(B->A)=0.5
pa=0.1;pb=0.05;
% 100为100个分子
for j=1:100
% 因为第一个状态已确认为A，所以从第二个开始算
for i=2:time
    % 先判断状态，这里假设1为状态A，0为状态B
    if (x(j,i-1)==1)
        % rand(1) 在0~1取一个随机数，计算能否从状态A到B，假设随机数小于pa，则转化为B
        if (rand(1)<pa)
            x(j,i)=0;
        % rand(1) 大于pa，则无法转化为B
        else
            x(j,i)=1;
        end
    % 同上，先判断状态为B，则进一步进行概率计算
    else
        if (rand(1)<pb)
            x(j,i)=1;
        else
            x(j,i)=0;
        end
    end
end
end
% 通过sum函数进行总和计算处于A状态的数量，其中sum函数对每一列进行加和，因此y可以统计每个时间处于A状态的数量
a=sum(x);
t=0:0.1:time/10-0.1;
% 根据老师课上的公式，从而可以得到宏观状态下的分子状态随时间的变化
aa=time-(time*2/3*(1-exp(-1.5*t)));
% 之后，进行绘图，画出上述得到的两个数据
figure;plot(t,a,t,aa)
%B 状态的数目
b=100-a;
bb=100-aa;
figure;plot(t,b,t,bb)

%计算error
figure;plot(t,a-aa,t,zeros(time,1));
figure;plot(t,b-bb,t,zeros(time,1));
%%
close all
clc
% 创建一个分子的随时间的变化，在A，B态间的转变
zz=zeros(1,10);
zz(1)=1;
pa=0.1;pb=0.05;
for i=2:10000
    if (zz(i-1)==1)
        if (rand(1)<pa)
            zz(i)=0;
        else
            zz(i)=1;
        end
    else
        if (rand(1)<pb)
            zz(i)=1;
        else
            zz(i)=0;
        end
    end
end
% plot(zz);
% ylim([-0.5 1.5]);

% 计算A分子的dewll时间，即在A状态的持续时间
num=find(zz(:)==0);
num=[0;num];
dwell=diff(num)-1;
dwell(dwell(:)==0)=[];
dwell=dwell./10;
histogram(dwell)
mean(dwell)