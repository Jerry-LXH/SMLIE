function [subRegionAssignment, splittedTracks, deletedTracks] = assign_tracks_to_regions(subRegions, origTracks, subRoiBorderHandling)

%% Assign tracks to sub-regions.
%
%
%Input:
%   subRegions              - 	Cell array containing one sub-region per cell. Each
%                               sub-region cell contains one cell per frame which again 
%                               contains one cell per region part where the coordinates
%                               are stored in a 2-column vector(one sub-region can be
%                               composed of several separate regions). If the region
%                               was drawn by hand, the region is saved in the first
%                               cell.
%                               Example: if you draw a single sub-region by hand it
%                               would be saved in subROI{1}{1}{1} because it is the
%                               first sub-region in the first frame and consists only
%                               of one part. If you draw another sub-region by hand it
%                               would be saved in subROI{2}{1}{1}. If the two
%                               sub-regions are merged so that the two drawn regions
%                               belong to the same sub-region number, they would be saved
%                               in subROI{1}{1}{1} and subROI{1}{1}{2}. If a sub-region
%                               is created using the "threshold" button, the sub-region
%                               is drawn for each frame. So subROI{1}{5}{2} would be
%                               the first sub-region in the 5th frame and the second
%                               part of multiple regions in this frame.
%   origTracks              -	Cell array containing the tracked molecules. Has
%                               as many cells as there are tracks. If spot
%                               positions were refined with the fit_spots
%                               function each cell contains an array with 9 columns:
%                               [frame,xpos,ypos,A,BG,sigma_x,sigma_y,angle,exitflag]
%   subRoiBorderHandling    -   String that defines how to proceed with tracks that
%                               cross sub-region border. Use one of the following:
%                               'Assign by first appearance' - Tracks are assigned
%                               to the sub-region where the first detection in the
%                               track appears in.
%                               'Split tracks at border' - Tracks are split at
%                               sub-region borders and each separate track is then
%                               assigned to its sub-region.
%                               'Delete tracks crossing borders' - Delete all
%                               tracks that thouch a sub-region border.
%                               'Only use tracks crossing borders'     - Use only the
%                               tracks that touch a sub-region border.
%
%
%Output:
%   subRegionAssignment     -   Array with one column that has as many entries
%                               as there are tracks in the output variable
%                               "splittedTracks". Each entry of the array represents
%                               the sub-region number in which the track
%                               appears (0 = tracking-region).
%   splittedTracks          -   Cell array containing the tracked molecules after
%                               processing the tracks that touch the
%                               sub-region borders.
%   deletedTracks           -   Cell array containing the deleted tracks if
%                               the user chose 'Only use tracks crossing borders' or 
%                               'Delete tracks crossing borders' 
%



if isempty(origTracks)
    %There are no tracks so abort assignement to sub-regions
    subRegionAssignment = [];
    splittedTracks = [];
    return
end

%Initialize cell array containing the tracks that have to be deleted
deletedTracks = {};

%Initialize cell array of tracks after splitting
splittedTracks = origTracks;

%Get number of tracks
nTracks = length(origTracks);

%Get number sub-regions
nRegions = length(subRegions);

%Initialize matrix of logicals where for each combination of track and
%sub-region the value is true if the track lies within this sub-region
regionOfTracks = false(nRegions,nTracks);

%Initialize logical array where we save wether a track has to be deleted
if strcmp(subRoiBorderHandling, 'Only use tracks crossing borders')
    tracksToDelete = true(nTracks,1);
else
    tracksToDelete = false(nTracks,1);
end

regionsOfTracksSplitted = cell(nTracks,1);

%Iterate through all tracks and asign each track to a sub-region
parfor curTrackIdx = 1:nTracks
        
    %Get the current track
    curTrack = origTracks{curTrackIdx};
    
    %Initialize vector that contains the sub-region number for each spot of 
    %the current track
    regionsOfCurTrack = zeros(size(curTrack,1),1);
    
    %------------Iterate through all sub-regions-------------------------------------------------
    for subRoiIdx = 1:nRegions
        
        %Get sub-region with the current loop index
        curSubRegion = subRegions{subRoiIdx};
        
        if length(curSubRegion) > 1 
            %Boundary varies with frame number because sub-region has been
            %drawn using an intensity threshold
            
            %Initialize array of logicals stating which spot of the current
            %track lies within the current sub-region
            isIn = false(size(curTrack,1),1);
            
            %Iterate through each spot of the current track
            for spotIdx = 1:size(curTrack,1)
                %Get the frame number in which the spot appears
                frameOfSpotInTrack = curTrack(spotIdx,1);
                
                %Get the sub-region that belongs to the frame in which the
                %spot appears
                subROIs = curSubRegion{frameOfSpotInTrack};
                
                %Check if roi exists in this frame
                if ~isempty(subROIs)
                    %Combine boundaries (sub-region might be composed of
                    %several small regions)
                    combinedSubRois = zeros(0,2);
                    for subRoiPartIdx=1:length(subROIs)
                        combinedSubRois = [combinedSubRois; subROIs{subRoiPartIdx}; [NaN NaN]];
                    end
                    combinedSubRois(end,:) = [];
                    
                    %Check if spot is inside this region and save the result
                    isIn(spotIdx) = inpolygon(curTrack(spotIdx,2),curTrack(spotIdx,3),combinedSubRois(:,1),combinedSubRois(:,2));
                end
            end
            
        else
            %Boundary does not vary with frame
            
            %Sub-region does not vary with frame so it is saved in the
            %first cell
            subROIs = curSubRegion{1};
            
            %Combine boundaries (sub-region might be composed of
            %several small regions)
            combinedSubRois = zeros(0,2);
            for subRoiPartIdx=1:length(subROIs)
                combinedSubRois = [combinedSubRois; subROIs{subRoiPartIdx}; [NaN NaN]];
            end
            
            combinedSubRois(end,:) = [];
            
            %Check if spot is inside this region and save the result
            
            isIn = inpolygon(curTrack(:,2),curTrack(:,3),combinedSubRois(:,1),combinedSubRois(:,2));
            
        end
        
        %----------Assign tracks depending on what the user selected------------------------------
        switch subRoiBorderHandling
            case 'Assign by first appearance' 
                %Assign tracks based on the position of the first spot in
                %the track
                if isIn(1)
                    %First spot of track was found in the current
                    %sub-region so set value for this region and track to
                    %true
                    regionOfTracks(subRoiIdx,curTrackIdx) = true;
                end
            case 'Split tracks at borders' 
                %Split tracks if they cross a sub-region border
                
                if all(isIn) 
                    %Particle was found in one of the boundaries and is 
                    %completely inside boundary.
                    
                    %Set value for this track in this region to true
                    regionOfTracks(subRoiIdx,curTrackIdx) = true;
                    regionsOfCurTrack(:) = subRoiIdx;
                elseif any(isIn) 
                    %Track was found in one of the boundaries but is 
                    %only partly inside boundary so Track has to be 
                    %split.
                    
                    %Save the sub-region number of the spots that lie within
                    %the current sub-region
                    regionsOfCurTrack(isIn) = subRoiIdx;
                    
                    %Set the value of the current track in the list of
                    %tracks that have to be deleted to true
                    tracksToDelete(curTrackIdx) = true;
                end
            case 'Delete tracks crossing borders' 
                %Delete whole track if track crosses a border
                
                if all(isIn)
                    %Track was found in one of the boundaries and is 
                    %completely inside boundary.
                    
                    %Set value for this track in this region to true
                    regionOfTracks(subRoiIdx,curTrackIdx) = true;
                elseif any(isIn) 
                    %Track was found in one of the boundaries but is 
                    %only partly inside boundary so delete this track.
                    
                    %Set the value of the current track in the list of
                    %tracks that have to be deleted to true
                    tracksToDelete(curTrackIdx) = true;
                end
            case 'Only use tracks crossing borders'                 
                %Use only tracks that cross a sub-region border and delete
                %all other tracks
                
                if any(isIn) && ~all(isIn)
                    %Track was found in one of the boundaries but is only 
                    %partly inside boundary so keep this track.
                    
                    %Set value for this track in this region to true
                    regionOfTracks(subRoiIdx,curTrackIdx) = true;
                    
                    %Set the value of the current track in the list of
                    %tracks that have to be deleted to false
                    tracksToDelete(curTrackIdx) = false;
                end
        end
        
    end
    
    regionsOfTracksSplitted{curTrackIdx} = regionsOfCurTrack;
    
end


%----------Split tracks if user selected to split them------------------------------
for curTrackIdx = 1:nTracks
    
    curTrack = origTracks{curTrackIdx};
    regionsOfCurTrack = regionsOfTracksSplitted{curTrackIdx};
    
    if strcmp(subRoiBorderHandling, 'Split tracks at borders')
        %Loop to assign each spot of the current track to a sub-region is
        %finished so now we can split the track at the region border
        
        
        %Find the points in the track where the track crosses from one
        %region to another
        changepoints = find(diff(regionsOfCurTrack));
        
        %Get the amount of crossings
        nRegionsOfTrack = numel(changepoints)+1;
        
        %Check if any spot of the current track lies within a sub-region
        if nRegionsOfTrack > 1
                        
            %Iterate through all crossings to split the track
            for regionIdx = 1:nRegionsOfTrack
                
                %Get the current sub-set of the track that lies between two
                %changepoints and get the corresponding sub-region
                if regionIdx == 1
                    curSplittedTrackPart = curTrack(1:changepoints(regionIdx),:);
                    regionOfSplittedTrack = regionsOfCurTrack(1);
                elseif regionIdx == nRegionsOfTrack
                    curSplittedTrackPart = curTrack(changepoints(regionIdx-1)+1:end,:);
                    regionOfSplittedTrack = regionsOfCurTrack(end);
                else
                    curSplittedTrackPart = curTrack(changepoints(regionIdx-1)+1:changepoints(regionIdx),:);
                    regionOfSplittedTrack = regionsOfCurTrack(changepoints(regionIdx));
                end
                
                %Add the track sub-set at the end of the tracks cell-array
                splittedTracks{end+1} = curSplittedTrackPart;
                
                %Add a new column in the logical matrix defining in which
                %regions the track appears
                regionOfTracks(:,end+1) = false;
                if regionOfSplittedTrack > 0
                    %Current track sub-set is located in one of the
                    %sub-regions so set the value of the corresponding
                    %sub-region to true
                    regionOfTracks(regionOfSplittedTrack,end) = true;
                end
            end
        end
    end
end

%----------Delete obsolete tracks and initialize assignment array----------

switch subRoiBorderHandling
    case 'Assign by first appearance'
        %Assign tracks by first position so no tracks have to be deleted
        
        %Initialize array that contains the sub-region number of each track
        subRegionAssignment = zeros(length(origTracks),1);
    case 'Split tracks at borders'
        %Tracks where splitted at border so delete the originals of the
        %now splitted tracks
        
        splittedTracks(tracksToDelete) = [];
        regionOfTracks(:,tracksToDelete) = [];
        
        %Initialize array that contains the sub-region number of each track
        subRegionAssignment = zeros(length(splittedTracks),1);
    otherwise
        %Get the tracks the have to be deleted
        deletedTracks = origTracks(tracksToDelete);
        
        %Delete the tracks from the tracks cell array and list containing
        %the sub-region of each track
        splittedTracks(tracksToDelete) = [];
        regionOfTracks(:,tracksToDelete) = [];
        
        %Initialize array that contains the sub-region number of each track
        subRegionAssignment = zeros(length(splittedTracks),1);
end

%-----------Assign tracks to the "highest" region--------------------------

%Check which region lies within another region and assign tracks to the 
%"highest" level region

%Iterate through sub-regions
for curRegionIdx = 1:nRegions
    
    %Get the number of tracks in the current sub-region
    nTracksInCurRegion = numel(find(regionOfTracks(curRegionIdx,:)));
    
    %Initialize new matrix of logicals where each track is only
    %assigned to one region
    tracksInAreaIdxNew = regionOfTracks;
    
    %Iterate through all sub-regions to compare the tracks in the first
    %region loop with the tracks in all other regions
    for comparisonRegionIdx = 1:nRegions
        %Compare only different regions
        if comparisonRegionIdx ~= curRegionIdx
            
            %Find tracks that appear in multiple regions
            tracksInMultipleRegions = regionOfTracks(curRegionIdx,:) == 1 & regionOfTracks(comparisonRegionIdx,:) == 1;
            
            %Get the number of tracks in the comparison region
            nTracksInCompRegion = sum(regionOfTracks(comparisonRegionIdx,:));
            
            %Get the number of tracks that appear in multiple regions
            nTracksInMultipleRegions = sum(tracksInMultipleRegions);
            
            if nTracksInMultipleRegions < nTracksInCurRegion && nTracksInMultipleRegions > 0 && nTracksInCompRegion == nTracksInMultipleRegions 
                %Comparison Region is inside current regions so we set the
                %logical of the current region to false (the index of the
                %comparison regions stays true)
                tracksInAreaIdxNew(curRegionIdx,tracksInMultipleRegions) = false;
                
                %Other possible cases that we do not use
                %                 elseif nTracksInMultipleRegions < nTracksInCurRegion && nTracksInMultipleRegions > 0 && nTracksInCompRegion ~= nTracksInMultipleRegions  %Comparison Region is overlapping with current regions
                %                 elseif nTracksInMultipleRegions == nTracksInCurRegion  %Current Region is inside Comparison Region
                %                 elseif nTracksInMultipleRegions == 0 %No overlap between regions
            end
            
            
        end
        
    end
    
    %Save the current region number for the tracks in the current region
    %in the array that contains the sub-region number of each track
    subRegionAssignment(tracksInAreaIdxNew(curRegionIdx,:)) = curRegionIdx;
end

end