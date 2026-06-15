clc;
clear all;

pathway = '/Users/lu_yk/Desktop/Data/UV/20240325_RNO0_DAN';
data_list = dir(fullfile(pathway,'*.csv'));

parallel_time = 3;%number of parallel tests

total = cell(parallel_time,size(data_list,1)/parallel_time);

% start_row = 11;%读取csv文件时起始行，第一行为0
% start_column = 0;

for i = 1:length(data_list)
    file_ptw = fullfile(pathway,data_list(i).name);
    data = readmatrix(file_ptw);
    x = data(:,1);
    y = data(:,2);

    range_left = 510;%需要寻找的波长
    range_right = 538;
    index_range = find(x >= range_left & x <= range_right);
    intensity_0 = max(y(index_range));

    %最后10个数据取均值扣去背景
    means = mean(y(end-10:end));

    %还原浓度荧光强度
    intensity = intensity_0-means;

    total(i) = num2cell(intensity);
end

table_0 = total';
x_label = {'0';'40';'80';'120';'160';'200'};
table = cell2table(table_0,'RowNames',x_label);
output = [pathway,'intensity.xlsx'];
writetable(table,output,'WriteRowNames',true);
clc;
disp('Completed');
