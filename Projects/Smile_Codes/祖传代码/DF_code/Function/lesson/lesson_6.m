%% 影像处理
idlcolormapdemo
s1=seq.imageData(:,:,1);
s1_mean=mean2(s1);
m=0;
clear st;
%% 筛选符合条件的数据
for i=1:256
    for j=1:256
        if (s1(i,j)>s1_mean)
            m=m+1;
            st(m,1)=i;st(m,2)=j;st(m,3)=s1(i,j);
        end
    end
end


%%
plot3(st(:,1),st(:,2),st(:,3))
clear z;
for i=1:size(st,1)
    z(st(i,1),st(i,2))=st(i,3);
end
%%
figure;
set(gcf,'position',[200,200,700 640])
mesh(x,y,z);
[x,y]=meshgrid(1:256,1:256);
xlim([1,256]);ylim([1,256]);box on

figure;
imagesc(double(z));
axis square
colormap("gray");

figure;
contour(z);axis square

%%
mesh(x,y,s1)

%%
aa=h5read('idlcolormaps.h5','/name')






