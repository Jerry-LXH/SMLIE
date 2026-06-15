function series = time_trace(loc, data)

% integrate with 7*7 pixels
[loc_num, ~] = size(loc);
[~, ~, frames] = size(data);
series = zeros(loc_num, frames); 
for num = 1:loc_num
    for frame = 1:frames   
        series(num, frame) = sum(data(loc(num,1)-3:loc(num,1)+3, loc(num,2)-3:loc(num,2)+3, frame), 'all');
    end
end
counts = size(series,1);
disp([num2str(counts),' traces in total']);
end