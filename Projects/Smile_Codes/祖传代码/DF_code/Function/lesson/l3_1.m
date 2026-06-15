%% first exerice
for i=1:10
    x=linspace(0,10,101);
    plot(x,sin(x+i));
    print(gcf,'-deps',strcat('plot',num2str(i),'.ps'));
end
%%
% second exerice
if rem(a,2)==0
    disp('a is even')
else
    disp('a is odd')
end

%%
n=1;
while prod(1:n)<1e100 %prod:1*2*...*n
    n=n+1;
end

%%
sum=0;
i=0;
while i<=999
    sum=sum+i;
    i=i+1;
end
sum

%%
i=1;
clear a
for n=1:2:10
    a(i)=2^n;
    a((n+1)/2);
    i=i+1;
end
disp(a)

%%
A=zeros(2000,2000);
tic
for ii=1:size(A,1)
    for jj=1:size(A,2)
        A(ii,jj)=ii+jj;
    end
end
toc

%%
A=[0 -1 4;9 -14 25;-34 49 64];
B=zeros(3,3);
C=zeros(3,3);
for ii=1:size(A,1)
    for jj=1:size(A,2)
        if A(ii,jj)<0
            C(ii,jj)=-A(ii,jj);
            B(ii,jj)=C(ii,jj);
        else
        B(ii,jj)=A(ii,jj);
        end
    end
end
A
B

%%
x=2;k=0;error=inf;
error_threshold=1e-32;
while error>error_threshold
    if k>100
        break
    end
    x=x-sin(x)/cos(x)
    error=abs(x-pi);
    k=k+1;
end

%%
% clear all
% close all
% clc
% ...：换行号
% Ctrl+c：直接停止程序

%%
% 第二单元：function



