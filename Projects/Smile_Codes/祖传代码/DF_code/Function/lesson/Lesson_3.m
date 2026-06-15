%变量类型：逻辑、数字
%char：
% 
clear
s1=input('str:','s');
n=length(s1);
for i=1:n
    s2(i)=s1(n-i+1);
end
s2

%%
%structure:异质的数据，
student.name='a';
student.id='sdjj';
student.number=1;

%%
student(2).name='a';
student(2).id='sdjj';
student(2).number=1;

%%
%nesting structures:结构体套娃

%%
%cell array:另一种保存异质数据的方法。使用{}宣告
A(1,1)={[1 4 3; 0 5 8]};
A(1,2)={'ajdsj'};
A(2,1)={3+7i};
A(2,2)={-pi:pi:pi};
%huozhe
% A{1,1}=[1 4 3; 0 5 8];

%%
%读取Cell Array
A{1,1}

%%
num2cell()
mat2cell()

%%
%多维数组，三维：row, column, layer
% 
% cat()
% 
% reshape()%改变横纵

A=[1:3;4:6];
B =reshape(A,3,2)


%%
%file acess
clear;a=magic(4);
save mydata1.mat
save mydata2.mat -ascii

load("mydata1.mat")
load("mydata2.mat",'-ascii')

%%

xlswrite('testscore',a,1,'A1:A3')

%get text in excel


%%
%low level file input/output

s=readsif('Y-92Yb8Er-Lu_Er2017_dilu15000x_300mA_01152021_0.1s_5frames_sample01152020_test1.sif');





