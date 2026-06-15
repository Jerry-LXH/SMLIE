function data_1=p_d(xd,ima1_1)
xd(1:2)=xd(1:2)+0.5;
row_edge=round(xd(1):xd(1)+xd(3));
col_edge=round(xd(2):xd(2)+xd(4));
col_edge(col_edge>256)=[];
col_edge(col_edge<0)=[];
row_edge(row_edge>256)=[];
row_edge(row_edge<0)=[];
imad=ima1_1(row_edge,col_edge);
imad=imad';
S=rescale(imad);
bg=imopen(S,strel('disk',15));
s2=S-bg;
% find the threshold
level=graythresh(s2)*1.5;
bw=imbinarize(s2,level);
[label,~]=bwlabel(bw,8);
% RGB_label=labelclc2rgb(label);imshow(RGB_label);
% 确定中心点
data_1=regionprops('table',label,imad,'Centroid',...
    'Area','MajorAxisLength','MinorAxisLength','MaxIntensity', ...
    'BoundingBox');
if size(data_1)>0
    data_1.Centroid(:,1)=data_1.Centroid(:,1)+xd(1)-1;
    data_1.Centroid(:,2)=data_1.Centroid(:,2)+xd(2)-1;
    data_1.BoundingBox(:,1)=data_1.BoundingBox(:,1)+xd(1)-1;
    data_1.BoundingBox(:,2)=data_1.BoundingBox(:,2)+xd(2)-1;
else
end

% centers=data_1.Centroid;
% figure;hold on;
% imagesc(imad);
% viscircles(centers,[2,2],'color','#0072BD','LineWidth',1.5,'EnhanceVisibility',false);
% hold off





