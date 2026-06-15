function F2C
while 1>0
    show='tem in F:';
    F=input(show);
    if isempty(F)==1
        break
    end
    C=(F-32).*5./9;
    Cstr=num2str(C);
    Cshow=['→','tem in C=',Cstr];
    disp(Cshow)
end
% function handle: 