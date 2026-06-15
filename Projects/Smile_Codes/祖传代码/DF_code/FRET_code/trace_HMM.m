function [state_num,likest,intensity_time]=trace_HMM(X,ex_time)
clear Hint HMMint state_num

% X是trace 
X=X';
Y=pdist(X,'euclidean');%Compute the pairwise distances
Z=linkage(Y,'average');%Generate clustering hierarchical tree according to distance information
%Calculate cophenetic correlation coefficient for the hierarchical cluster tree, a larger value indicates that the tree fits the distance well
C=cophenet(Z,Y)

%%%%%%Natural Divisions
%Set cutoff threshold according from standard deviation and cophenetic correlation coefficient
I1=find(X~=0);
cutoff=20*std(X(I1((end-50):end)));%If the number of categories
%is small,try increasing the cutoff,otherwise reducing the cutoffD
%clustering result
seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
seq=seq';

%     Determine emission value
F1=length(unique(seq));
intensity=zeros(1,F1);
for i=1:F1
    intensity(i)=mean(X(seq==i));
end

Ctrace=zeros(1,length(X));
for i=1:length(intensity)
    Ctrace(seq==i)=intensity(i);
end
R=corrcoef(Ctrace',X);
i=0.1;
% R(1,2);

while R(1,2)<0.88
    cutoff=(20-i)*std(X(end-50:end));
    i=i+0.1;
    seq = cluster(Z,'cutoff',cutoff,'Criterion','distance');
    seq=seq';
    F1=length(unique(seq));
    intensity=zeros(1,F1);
    for j=1:F1
        intensity(j)=mean(X(seq==j));
    end
    for j=1:F1
        Ctrace(seq==j)=intensity(j);
    end
    R=corrcoef(Ctrace,X);
end
clear lloc
% Baum-welch algorithm is used to estimate the transfer probability, emission probability
lloc=length(intensity);
TRANS_GUESS = eye(lloc)*1+rand(lloc,lloc)*0.0000001;
EMIS_GUESS = eye(lloc)*0.99+rand(lloc,lloc)*0.01;
[TRANS_EST2, EMIS_EST2] = hmmtrain(seq, TRANS_GUESS, EMIS_GUESS);

%%%%%% The state sequence is estimated by Viterbi algorithm
likelystates = hmmviterbi(seq, TRANS_EST2, EMIS_EST2);
likest=likelystates;
state_num=length(intensity);
for i=1:length(intensity)
    likest(likelystates==i)=intensity(i);
end
if state_num ==2
    if likest(1)==max(intensity)
        intensity_time=[max(intensity),sum(likest==max(intensity)).*ex_time];
    else
        intensity_time=[0,0];
    end

else
    intensity_time=[0,0];
end
%Draw analysis results


% if length(intensity)<=3 && length(intensity)~=1
%     good_num=good_num+1;
% end

% state_num=state_num';
% clear a
% k=1;
% for i=1:size(state_num,1)
%     if state_num(i)>2 | state_num(i)==1 | state_num(i)==0
%         a(k)=i;
%         k=k+1;
%     end
% end
