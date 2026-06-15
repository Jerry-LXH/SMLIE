function subRegionAssignment = assign_spots_to_regions(subRegions, spots)

%Assign a list of non-linked spots to the sub-region where they appear.
%
%
%Input: 
%   subRegions  -   Cell array containing one sub-region per cell. Each
%                   sub-region cell contains one cell per frame which again 
%                   contains one cell per region part where the coordinates
%                   are stored in a 2-column vector(one sub-region can be
%                   composed of several separate regions). If the region
%                   was drawn by hand, the region is saved in the first
%                   cell.
%                   Example: if you draw a single sub-region by hand it
%                   would be saved in subROI{1}{1}{1} because it is the
%                   first sub-region in the first frame and consists only
%                   of one part. If you draw another sub-region by hand it
%                   would be saved in subROI{2}{1}{1}. If the two
%                   sub-regions are merged so that the two drawn regions
%                   belong to the same sub-region number, they would be saved
%                   in subROI{1}{1}{1} and subROI{1}{1}{2}. If a sub-region
%                   is created using the "threshold" button, the sub-region
%                   is drawn for each frame. So subROI{1}{5}{2} would be
%                   the first sub-region in the 5th frame and the second
%                   part of multiple regions in this frame.
%   spots       -   List of spots as an array with at least 3 columns:
%                   [frame, xPos, yPos].
%   
%Output:
%   subRegionAssignment -   Array with one column that has as many entries
%                           as there are spots in the input "spots"
%                           variable. Each entry of the array represents
%                           the sub-region number in which the spot
%                           appears (0 = tracking-region).
%       

if isempty(spots)
    subRegionAssignment = [];
    return
end

nSpots = size(spots,1);
nRegions = length(subRegions);
framesWithSpots = unique(spots(:,1));
subRegionAssignment = zeros(nSpots,1);

for m = 1:nRegions
    curBoundary = subRegions{m};
    
    if length(curBoundary) > 1
        %Boundary varies with frame number because sub-region has been
        %drawn using an intensity threshold
        
        %Iterate through frames
        for frameIdx = framesWithSpots'
            %Create array of logicals to get the indices of spots in the
            %current frame
            indicesOfSpotsInCurFrame = spots(:,1) == frameIdx;
            
            %Convert array of logicals to numerical index values
            idxOfSpotsInCurFrame = find(indicesOfSpotsInCurFrame);
            
            %Iterate through all spots in the current frame
            for spotInCurFrameIdx = idxOfSpotsInCurFrame'
                
                %Get coordinates of current spot
                curSpotX = spots(spotInCurFrameIdx,2);
                curSpotY = spots(spotInCurFrameIdx,3);
                
                
                for n = 1:length(curBoundary{frameIdx})
                    %Check if spot is found in this part of the sub-region
                    isIn = inpolygon(curSpotX,curSpotY,curBoundary{frameIdx}{n}(:,1),curBoundary{frameIdx}{n}(:,2));
                    
                    if isIn
                        %Spot was found in one of the boundaries so loop can be stopped
                        
                        %Save sub-region number for current spot
                        subRegionAssignment(spotInCurFrameIdx) = m;
                        break;
                    end
                end
                
            end
            
        end
        
    else
        %Boundary does not vary with frame number because sub-region
        %was hand-drawn
        
        %Iterate through all parts of the sub-region
        for n = 1:length(curBoundary{1}) 
            %Check which spots lie within the current part of the sub-region
            isIn = inpolygon(spots(:,2),spots(:,3),curBoundary{1}{n}(:,1),curBoundary{1}{n}(:,2));
            
            %Save the sub-region number of the spots that have been found
            %in this part current sub-region
            subRegionAssignment(isIn) = m;
        end
    end
    
end

end
