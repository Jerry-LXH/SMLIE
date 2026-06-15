function y = round_significant(x, n, method)

%Round numeric value to the first n significant digits
%
%y = round_significant(x, n, method)
%
%Input:
% x        -   Numeric input value
% n        -   Number of significant digits
% method   -   'ceil' rounds towards positive infinity, 'floor' rounds
%              towardsnegative infinity, 'round' rounds to nearest decimal
%              or integer
%
%Output:
% y        -   Numeric value round to the first n significant digits


e = floor(log10(abs(x)) - n + 1);
og = 10.^abs(e);
y = feval(method,x./og).*og;
k = find(e<0);
if ~isempty(k)
    y(k) = feval(method,x(k).*og(k))./og(k);
end
y(x==0) = 0;
end