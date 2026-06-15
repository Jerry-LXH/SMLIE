clc;
clear all;

pathway = '/Users/lu_yk/Desktop/Data/FLS/20250407_Rh-5_Gly';
data_list = dir(fullfile(pathway,'*.txt'));


total = zeros(800-550+1,size(data_list,1)+1);%！！！！！！修改行数！！！！！！


for i = 1:length(data_list)
    file_ptw = fullfile(pathway,data_list(i).name);
    data = readmatrix(file_ptw);
    if i == 1
    % total(:,1) = data(:,1);%Abs.
    total(:,1) = data(20:end,1);%emission
    end
    %y = data(:,2);
    %means = mean(y(end-10:end));
    %total(:,i+1) = data(:,2);%-means;
    % total(:,i+1) = data(:,2)-data(end,2);%abs
    total(:,i+1) = data(20:end,2);%emission
end

table = array2table(total);
output = [pathway,'intensity.xlsx'];
writetable(table,output,'WriteRowNames',true);
clc;
disp('Completed');
