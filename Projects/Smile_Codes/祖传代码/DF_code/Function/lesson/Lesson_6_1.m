%%
%提高对比度，等比扩大，线性方程y=kx+b，根据最大，最小值算得k，b
I=imread('rice.png');
subplot(2,2,1);imshow(I);subplot(2,2,2);imhist(I); 

m1=double(max(I(:)))
m2=double(min(I(:)))
k=255/(m1-m2)
b=255-255*m1/(m1-m2)
% for i=1:size(I,1)
%     for j=1:size(I,2)
%         I(i,j)=k*I(i,j)+b;
%     end
% end
I=I.*k+b;
subplot(2,2,3);imshow(I);subplot(2,2,4);imhist(I);


%%
imadd();


%%
% image histogram;
imhist(imread('rice.png'))
histeq();

%%
%transformation
J=imrotate(I,45,"bicubic")
imshow(J)


%%
M=zeros(600);
for i=1:size(I,1)
    for j=1:size(I,2)
        M_i=round(cos(pi/2)*i-sin(pi/2)*j)+160;
        M_j=round(sin(pi/2)*i+cos(pi/2)*j);
        M(M_i,M_j)=I(i,j);
    end
end
imagesc(M)


%%
% M=zeros(360);
clear M_i
clear M_j
clear M
for i=1:size(I,1)
    for j=1:size(I,2)
        M_i=round(cos(pi/4)*i-sin(pi/4)*j+sin(pi/4)*256+2);
        M_j=round(sin(pi/4)*i+cos(pi/4)*j+2);
            M(M_i,M_j)=I(i,j);
    end
end

subplot(1,2,1);imagesc(M);colormap('gray');axis square;

for i=2:size(M,1)-1
    for j=2:size(M,2)-1
        if(M(i,j)==0)
            q=[M(i-1,j) M(i-1,j-1) M(i-1,j+1) M(i,j-1) M(i-1,j+1) M(i+1,j-1) M(i+1,j) M(i+1,j+1)];
            M(i,j)=mean(q);
        end
    end
end

subplot(1,2,2);imagesc(M);colormap('gray');axis square;


%% 导出图片
imwrite()
