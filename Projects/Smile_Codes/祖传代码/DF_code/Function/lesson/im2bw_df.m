function S=im2bw_df(s,level)

S=zeros(size(s));
for i=1:size(s,1)
    for j=1:size(s,2)
        if(s(i,j)>=level*256)
            S(i,j)=1;
        else
            S(i,j)=0;
        end
    end
end

