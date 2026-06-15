%% image thresholding(二值化）

idlcolormapdemo

%%
% 二值化
s=double(seq.imageData(:,:,1));
S=rescale(s);
bg=imopen(S,strel('disk',15));
s2=S-bg;
level=graythresh(s2);% find the threshold
bw=imbinarize(s2,level);
% subplot(1,2,1);imshow(bw);subplot(1,2,2);imagesc(label);colormap('gray');axis square
[label,object]=bwlabel(bw,8);
RGB_label=label2rgb(label);imshow(RGB_label);
data_1=regionprops(label,'basic');%% 确定中心点
%%
clear x y z
x=round(data_1(11).BoundingBox);
z=label(x(2):x(2)+x(4),x(1):x(1)+x(3))
%%
stats=regionprops('table',label,'Centroid',...
    'MajorAxisLength','MinorAxisLength')
centers = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2+1;
imagesc(s);colormap('gray');axis square
hold on
viscircles(centers,radii,'color','#0072BD','LineWidth',1.5,'EnhanceVisibility',false);
hold off
%% 统计每个点的大小
clear sum;
sum=zeros(1,object);
for x=1:object
    for i=1:size(label,1)
        for j=1:size(label,2)
            if(label(i,j)==x)
                sum(x)=sum(x)+1;
            end
        end
    end
end
max_s=max(sum);
histogram(sum);
%%  将数据变为红色
s3=zeros(256,256,3);
for i=1:size(label,1)
    for j=1:size(label,2)
        if(label(i,j)>0)
            s3(i,j,:)=[255,0,0];
        end
    end
end
imshow(s3);
%%
I=imread('rice.png');
level=graythresh(I);% find the threshold
bw1=im2bw_df(I,level);bw2=imbinarize(I,level);
subplot(1,3,1);imagesc(I);colormap('gray');axis square
subplot(1,3,2);imagesc(bw1);colormap('gray');axis square;
subplot(1,3,3);imagesc(bw2);colormap('gray');axis square;



%% 去除背景影响
I=imread('rice.png');
bg=imopen(I,strel('disk',15));
I2=I-bg;
level=graythresh(I2);% find the threshold
bw=imbinarize(I2,level);
subplot(1,3,1);imshow(I);
subplot(1,3,2);imshow(bg)
subplot(1,3,3);imshow(bw)
figure;surf(bg);


%%
ima1=double(seq.imageData(:,:,14));
[x,y]=meshgrid(1:512,1:512);
surf(x,y,ima1)
colormap("cool")


