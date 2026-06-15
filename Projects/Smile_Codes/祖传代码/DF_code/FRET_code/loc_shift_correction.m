function loc1=loc_shift_correction(loc1,loc2,gap)
matched_points_A=[];
matched_points_B=[];
for i = 1:size(loc1, 1)
    for j = 1:size(loc2, 1)
        distance = norm(loc1(i, :) - loc2(j, :));
        if distance < gap
            loc1(i,:) = loc2(j,:);
            matched_points_A=[matched_points_A;loc1(i,:)];
            matched_points_B=[matched_points_B;loc2(j,:)];
        end
    end
end
end

