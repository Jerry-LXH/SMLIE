%% 从Abs的CSV文件读取吸收数据，并且将所有的数据整合到一个文件里
clear;clc
path='C:\Users\23342\Desktop\Data\Abs\pH Sensitive\SO-DPBF\Cy5-pH\30uM+1%HOAc\';
files = dir(fullfile(path, '*.csv*'));
% 获取该文件下的每个数据文件
datanames = {files(~[files.isdir]).name}';
% 获取数据的数量
datanumber=size(datanames,1);
% 设置保存路径
C=strsplit(path,'\');
data_save_name=strcat(path,string(C(end-1)),'_Abs.xlsx');
data=[];  
m=[];
name='Wavelength';
figure;
hold on
for x=1:datanumber
    dataname_mat=cell2mat(datanames(x,:));
    time=strsplit(dataname_mat,'-');
    time_tem=cell2mat(time(end-1));
    data_max(x,1)=str2double(time_tem(1:2));
    % 设置文件读取的路径
    pathname=strcat(path,dataname_mat);
    data_tem=readmatrix(pathname);
    m(x,1)=mean(data_tem(end-50:end,2));
    data_tem(:,2)=data_tem(:,2)-mean(data_tem(end-50:end,2));
    data(:,x+1)=data_tem(:,2);
    if x==1
        data(:,1)=data_tem(:,1);
    end
    data_max(x,2)=max(data_tem(230:250,2));
    plot(data(21:end,1),data_tem(21:end,2));
    name=[name;string(dataname_mat(1:end-4))];
end
hold off
figure;
%data_max(:,2)=data_max(:,2)./data_max(1,2);
%data_max=sortrows(data_max,1);
plot(data_max(:,1),data_max(:,2))

datat=array2table(data,"VariableNames",name);
writetable(datat,data_save_name);
writematrix(data_max,strcat(path,string(C(end-1)),'_Abs_max.xlsx'))