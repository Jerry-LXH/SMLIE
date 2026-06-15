function minDist = distance_to_roi_border(P, roi)

%minDist = distance_to_roi_border(P, roi)
%
% Find the minimum distance between a point P and a polygonal region
% defined by roi.
%
%
% Input:
%   P           -   Array with two entries [xPos, yPos]
%   roi         -   2-column vector containing coordinates of a region of interest. 
%                   1st column x-coordinates (horizontal), 2nd column y-coordinates 
%                   (vertical).
%
% Output: Struct out with fields
%   minDist     -   Minimum distance between a point P and a polygonal 
%                   regiondefined by roi

minDist = inf;

for subRoiPointIdx = 1:size(roi,1)-1
    
    A = roi(subRoiPointIdx,1:2);
    B = roi(subRoiPointIdx+1,1:2);
    
    %Vector AB
    AB(1) = B(1) - A(1);
    AB(2) = B(2) - A(2);
    
    %Vector BP
    BP(1) = P(1) - B(1);
    BP(2) = P(2) - B(2);
    
    %Vector AP
    AP(1) = P(1) - A(1);
    AP(2) = P(2) - A(2);
    
    %Calculating the dot product
    AB_BP = AB(1) * BP(1) + AB(2) * BP(2);
    AB_AP = AB(1) * AP(1) + AB(2) * AP(2);
    
    
    if (AB_BP > 0)
        %B is closest
        
        %Finding the magnitude
        y = P(2) - B(2);
        x = P(1) - B(1);
        distanceToCurrentRegionPoint = sqrt(x^2 + y^2);
        
        
    elseif (AB_AP < 0)
        %A is closest
        y = P(2) - A(2);
        x = P(1) - A(1);
        distanceToCurrentRegionPoint = sqrt(x^2 + y^2);
        
        
    else
        %Finding the perpendicular distance
        
        x1 = AB(1);
        y1 = AB(2);
        x2 = AP(1);
        y2 = AP(2);
        mod = sqrt(x1 * x1 + y1 * y1);
        distanceToCurrentRegionPoint = abs(x1 * y2 - y1 * x2) / mod;
    end
    
    if distanceToCurrentRegionPoint < minDist
        minDist = distanceToCurrentRegionPoint;
    end
end

end