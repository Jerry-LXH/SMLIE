function location=filter_duplicate_points(loc,threshold)
% 假设你的点集存储在一个Nx2的矩阵中，X
X = loc; % 示例点集
N = size(X, 1);
keep = true(N, 1);
num=[];
n=1;
for i = 1:N
    for j = i+1:N
        if norm(X(i, :) - X(j, :)) < threshold
            keep(j) = false;
            num(n,:)=[i j];
            n=n+1;
        end
    end
end
% 判断是否进行平均
if threshold<2
    if(isempty(num)==0)
        for i=1:size(num,1)
            X(num(i,1),:)=mean([X(num(i,1),:);X(num(i,2),:)]);
        end
    end
end

location = X(keep, :);


