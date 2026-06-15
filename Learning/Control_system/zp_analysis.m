numg=[6 0 1]; deng=[1 3 3 1];sysg=tf(numg,deng); % 创建传递函数对象
z=zero(sysg);
p=pole(sysg);
n1=[1 1]; n2=[1 2]; d1=[1 2*i]; d2=[1 -2*i]; d3=[1 3]; % 定义H的分子分母
numh=conv(n1, n2); denh=conv(d1, conv(d2,d3)); % 等价于多项式相乘，H展开
sysh=tf(numh, denh);
sys=sysg/sysh; % G/H 组合
pzmap(sys);