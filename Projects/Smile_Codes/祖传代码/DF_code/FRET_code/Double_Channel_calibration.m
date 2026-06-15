%% 计算双通道间的偏差，并使用三次函数进行拟合
function [fit_result,T,loc_red,loc_green]=Double_Channel_calibration(path,name)
seq=readsif([path name]);
imadata=double(seq.imageData(:,:,1));
imadata=imadata.*3.75./300;
dataedge=128;
green=imadata(dataedge:255+dataedge,1:256);
red=imadata(dataedge:255+dataedge,257:512);
gap=5;
loc_green=loc_molecule_DF(green);
loc_red=loc_molecule_DF(red);
loc_red=filter_duplicate_points(loc_red);
loc_green(loc_green(:,1)>251|loc_green(:,1)<6|loc_green(:,2)>251|loc_green(:,2)<6,:)=[];
loc_red(loc_red(:,1)>251|loc_red(:,1)<6|loc_red(:,2)>251|loc_red(:,2)<6,:)=[];

matched_points_A=[];
matched_points_B=[];
for m = 1:size(loc_red, 1)
    for n = 1:size(loc_green, 1)
        distance = norm(loc_red(m, :) - loc_green(n, :));
        if distance < gap
            matched_points_A = [matched_points_A; loc_red(m, :)];
            matched_points_B = [matched_points_B; loc_green(n, :)];
        end
    end
end
loc_red=matched_points_A;
loc_green=matched_points_B;

x = loc_red(:,1);
y = loc_red(:,2);
u = loc_green(:,1);
v = loc_green(:,2);
% 创建非线性转换矩阵
% 构造设计矩阵
X = [ones(size(x)), x, y, x.^2, y.^2, x.*y, x.^3, y.^3, x.^2.*y, x.*y.^2];
% 求解u = X * a
a = X \ u;
% 求解v = X * b
b = X \ v;
fit_result=[a b]
%运用最后一帧进行验证
loc_green1(:,1)=X*a;
loc_green1(:,2)=X*b;

writematrix(fit_result,[path 'calibration_result.xlsx'],"WriteMode","append")

writematrix(loc_red,[path 'loc_red.xlsx'],"WriteMode","append")
writematrix(loc_green,[path 'loc_green.xlsx'],"WriteMode","append")
% 计算线性转移矩阵
A = loc_red;
B = loc_green;
% 添加齐次坐标
A_homogeneous = [A, ones(size(A, 1), 1)];
B_homogeneous = [B, ones(size(B, 1), 1)];
% 使用最小二乘法求解转换矩阵
T = (A_homogeneous' * A_homogeneous) \ (A_homogeneous' * B_homogeneous);
% 将A矩阵中的坐标转换为B矩阵中的坐标
A_transformed = A_homogeneous * T;

writematrix(T,[path 'calibration_result_transformed.xlsx'],"WriteMode","append")

figure;
imagesc(green);
colormap('gray')
hold on
plot(loc_green(:,1),loc_green(:,2),'Color','green','Marker','.','LineStyle','none')
viscircles([loc_green(:,1),loc_green(:,2)],3.*ones(size(loc_green,1),1),'color','green','LineWidth',1.5,'EnhanceVisibility',false);
hold off

figure;
imagesc(red);
colormap('gray')
hold on
plot(loc_red(:,1),loc_red(:,2),'Color','red','Marker','.','LineStyle','none')
viscircles([loc_red(:,1),loc_red(:,2)],3.*ones(size(loc_red,1),1),'color','red','LineWidth',1.5,'EnhanceVisibility',false);
hold off

figure
imagesc(green);
title('transformed')
colormap('gray')
hold on
plot(A_transformed(:,1),A_transformed(:,2),'Color','green','Marker','.','LineStyle','none')
viscircles([A_transformed(:,1),A_transformed(:,2)],3.*ones(size(loc_green,1),1),'color','green','LineWidth',1.5,'EnhanceVisibility',false);
hold off


figure
imagesc(green);
title('poly')
colormap('gray')
hold on
plot(loc_green1(:,1),loc_green1(:,2),'Color','green','Marker','.','LineStyle','none')
viscircles([loc_green1(:,1),loc_green1(:,2)],3.*ones(size(loc_green,1),1),'color','green','LineWidth',1.5,'EnhanceVisibility',false);
hold off