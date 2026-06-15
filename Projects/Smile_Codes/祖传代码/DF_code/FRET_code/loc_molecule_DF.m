function centers=loc_molecule_DF(ima1,factor)
ima1_1=ima1';
% 二值化
S=rescale(ima1);
% 确定并去除背景值
% bg=imopen(S,strel('disk',4));
% S=S-bg;
% find the threshold
level=graythresh(S);
bw=imbinarize(S,level*factor);
% 当图像比较乱时，可以通过该方法起到一定效果
bw = bwmorph(bw,"clean");
% bw = bwmorph(bw,"open");
[label,~]=bwlabel(bw,4);
% 将regionprops函数得到的数据存入到stats
stats=regionprops('table',label,ima1,'Centroid','Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');
% 判断找到的点的数目，如果点的数目太少，减小level，重新识别
if size(stats,1)<5
    bg=imopen(S,strel('disk',4));
    S=S-bg;
    bw=imbinarize(S,level*0.8);
    % bw = bwmorph(bw,"open");
    [label,~]=bwlabel(bw,4);
    % 将regionprops函数得到的数据存入到stats
    stats=regionprops('table',label,ima1,'Centroid','Area','MajorAxisLength','MinorAxisLength','MaxIntensity','BoundingBox');
end
num_d=find(stats.Area>25);
% 找寻多个光斑连到一起的区域，重复上述运算，并将其分开
if (isempty(num_d)==0)
    for i=1:size(num_d,1)
        state_d=stats(num_d(i),:);
        x_d=state_d.BoundingBox;
        state_a=p_d(x_d,ima1_1);
        stats=[stats;state_a];
    end
end
% 删去光斑连在一起的数组
stats(num_d,:)=[];
% 当识别的点比较多时，即识别错了
if size(stats,1)>200
    stats(stats.MaxIntensity<1.5*mean(ima1(:)),:)=[];
    stats(stats.Area<5,:)=[];
end
stats(stats.MaxIntensity<2*mean(ima1(:)),:)=[];
stats(stats.Area<5,:)=[];
% 去除重叠的值
stats(stats.MaxIntensity>2.5*mean(stats.MaxIntensity),:)=[];
centers = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2+1;

% figure;
% imagesc(ima1);colormap('gray');axis square;colorbar;
% % 圈出每个数据点
% hold on
% viscircles(centers,radii,'color','#4F94CD','LineWidth',1.5,'EnhanceVisibility',false);
% % plot(fr_D(:,2),fr_D(:,1),'.','Color','red')
% hold off
