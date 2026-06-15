clc;
clear all;

pathway = '/Users/lu_yk/Desktop/Data/UV/20240202_RNO0-DAN_DMSO_NaOH';
data_list = dir(fullfile(pathway,'*.csv'));
parallel_time = 3;%number of parallel tests
total = cell(parallel_time,size(data_list,1)/parallel_time);
start_row = 11;%读取csv文件时起始行，第一行为0
start_column = 0;

for i = 1:length(data_list)
    file_ptw = fullfile(pathway,data_list(i).name);
    data = csvread(file_ptw,start_row,start_column);

    wavelength = 380;%需要寻找的波长
    intensity_0 = data(find(data(:,1) == wavelength),2);

    %还原浓度荧光强度
    intensity = 1*intensity_0;

    total(i) = num2cell(intensity);
end

table_0 = total';
x_label = {'0';'40';'80';'120';'160';'200'};
table = cell2table(table_0,'RowNames',x_label);
output = [pathway,'intensity.xlsx'];
writetable(table,output,'WriteRowNames',true);
clc;
disp('Completed');
