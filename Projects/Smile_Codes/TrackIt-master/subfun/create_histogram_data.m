function results = create_histogram_data(batch,moviesIdx, ui)

%
% Creates all the data needed for the TrackIt data analysis tool (data_analysis_tool.m)
%
% Input:
%   batch       -   TrackIts main variable where all movie informations,
%                   analysis parameters and results are stored 
%                   (see subfun/init_batch.m for more information)
%   moviesIdx   -   List of movies which should be included to create the results
%   ui          -   Struct array containing handles to all ui elements of the 
%                   data analysis tool. Needed to get the user input eg.
%                   "No. of jumps to consider".
%   
% Output:
%   results     -   Struct array containing all results that can be
%                   displayed in the data analyis tool (data_analysis_tool.m)

                               

%-----Get all settings from the user interface-----------------------------

%Get currently selected Region
regionNum = ui.popRegionSel.Value-1;     

%Save parameters that are essential for results calculations in a struct
para.nJumpsToConsider = str2double(ui.editJumpsToConsider.String);           %Number of jumps in on track that should be considered for mobility analysis
para.boolRemoveGaps = ui.cboxRemoveGaps.Value;                               %Wether user wants to remove jumps and angles where gap frames are involved
para.alphaThres = str2double(ui.editAlphaThres.String);                      %Get alpha value threshold for msd analysis
para.boolPixelsAndFrames = ui.btnPxFr.Value;                                 %Check if units should be displayed in pixels or frames
para.nDarkForLong = str2double(ui.editNDarkForLong.String);                  %Number of frames or dark periods to count as "long track"
para.nBrightFrames = str2double(ui.editNBrightFrames.String);                %Number of subsequent bright frames before a dark period (1 for continuous, >1 for ITM)
para.minJumpDistForAngles = str2double(ui.editAnglesMinJumpDist.String);        %Minimum jump distance of the jumps making up the jump angles
para.maxJumpDistForAngles = str2double(ui.editAnglesMaxJumpDist.String);        %Maximum jump distance of the jumps making up the jump angles


if para.boolPixelsAndFrames
    %Results are displayed in pixels and frames
    para.pixelsize = 1;
else
    %Results are displayed in microns and seconds
    para.pixelsize = str2double(ui.editPixelsize.String);
end

%-------------Initialize variables-----------------------------------------

%Get number of movies
nMovies = numel(moviesIdx);

movieNames = cell(nMovies,1);
jumpDistances = cell(nMovies,1);
startEndFrameOfTracks = cell(nMovies,1);
distToRoiBorder = cell(nMovies,1);
angles = cell(nMovies,1);
meanJumpDists = cell(nMovies,1);
trackLengths = cell(nMovies,1);
meanTrackLength = zeros(nMovies,1);
nTracks = zeros(nMovies,1);
nNonLinkedSpots = zeros(nMovies,1);
roiSize = zeros(nMovies,1);
meanTracksPerFrame = zeros(nMovies,1);
meanSpotsPerFrame = zeros(nMovies,1);
frameCycleTimes = zeros(nMovies,1);
alphaValues = cell(nMovies,1);
confRad = cell(nMovies,1);
meanJumpDistConfRad = cell(nMovies,1);
msdDiffConst = cell(nMovies,1);
trackingRadii = zeros(nMovies,1);
meanJumpDistMoviewise = zeros(nMovies,1);
nJumps = zeros(nMovies,1); 

%----------Iterate through all movies in moviesIdx-------------------------
for mIdx = 1:nMovies 
    
    %Get movie number of the current movie
    batchMovieNum = moviesIdx(mIdx);    
    
    %Get results, movie infos and tracking parameter of the current movie
    curMovieResults = batch(batchMovieNum).results;
    curMovieInfos = batch(batchMovieNum).movieInfo;
    curMovieParams = batch(batchMovieNum).params;
    
    if para.boolPixelsAndFrames
        %Results are displayed in pixels and frames
        
        %Set frame cycle time to 1 frame
        frameCycleTimes(mIdx) = 1;
        
        %Get the tracking radius in pixels
        trackingRadii(mIdx) = curMovieParams.trackingRadius;
    else
        %Unit is seconds and microns
        
        %Get the tracking radius and convert it to microns
        trackingRadii(mIdx) = curMovieParams.trackingRadius*para.pixelsize;
        
        %Get frame cycle times in seconds
        frameCycleTimes(mIdx) = curMovieInfos.frameCycleTime/1000;
    end
    
    
    %Get current movie name
    movieNames{mIdx} = curMovieInfos.fileName;
    
    %Calculate average number of tracks per frame
    nFramesAnalyzed = curMovieResults.nFramesAnalyzed;
    
    %Save first and last frame of all tracks
    startEndFrameOfTracks{mIdx} = curMovieResults.startEndFrameOfTracks;
    
    %Get tracklengths and meanTrackLengths and multiply with the
    %respective frame cycle time
    trackLengths{mIdx} = curMovieResults.trackLengths.*frameCycleTimes(mIdx);
    
    %Get cell array containing the jump distances of each track in a separate cell
    curMovieJumps = curMovieResults.jumpDistances;
    
    %Get the array containing the mean jump distance of each track
    curMovieMeanJumps = curMovieResults.meanJumpDists.*para.pixelsize;
    
    %Get cell array containing the jump angle of each track in a separate cell
    curMovieAngles = curMovieResults.angles;
    
    %Cell array containing the tracks of the current movie
    curMovieTracks = curMovieResults.tracks;
    
    %Get number of subRegions
    nRegionsInMovie = curMovieResults.nSubRegions+1;
    
    %Get distance of each track segment to closest roi border
    if isfield(curMovieResults, 'distToRoiBorder')
        %User clicked on "Calculate distance" so these data exist
        curMovieDistToRoiBorder = curMovieResults.distToRoiBorder;
    else
        %"Calculate distance" button not yet pressed so create empty array
        curMovieDistToRoiBorder = {};
    end
    
    %Check if results of all tracks can be displayed or if the results
    %of a specific sub-region must to be displayed
    if (numel(regionNum) == 1 && regionNum == 0) || (nRegionsInMovie == 1 && regionNum(1) == 1)

        %------------No subregion exists or user chose "All regions"-------
        
        %Get the number of tracks in the current movie
        nTracks(mIdx) = curMovieResults.nTracks;
        
        %Get the number of non-linked spots in the current movie
        nNonLinkedSpots(mIdx) = curMovieResults.nNonLinkedSpots;
        
        %Get the average tracklength in the current movie
        meanTrackLength(mIdx) = curMovieResults.meanTrackLength.*frameCycleTimes(mIdx);
        
        %Get the ROI size of the current movie
        roiSize(mIdx) = curMovieResults.roiSize.*para.pixelsize^2;
        
        %Get the average number of spots and tracks per frame
        meanTracksPerFrame(mIdx) =  curMovieResults.meanTracksPerFrame;
        meanSpotsPerFrame(mIdx) =  curMovieResults.meanSpotsPerFrame;
                
        %List that states which tracks belong to the current sub-region is
        %obviously a list of logical ones.
        tracksSubRoi = true(nTracks(mIdx),1);
        
        %-------Adjust results for user defined restrictions---------------
                
        %Adjust all values based on the user input that restricts our tracks
        %(eg. "nJumps to consider", "Remove jumps over gap frames etc.)      
        
        resultsIn.tracks = curMovieTracks;
        resultsIn.jumpDists = curMovieJumps;
        resultsIn.meanJumpDists = curMovieMeanJumps;
        resultsIn.angles = curMovieAngles;
        resultsIn.distToRoiBorder = curMovieDistToRoiBorder;
        
        resultsOut = adjust_jumps_and_angles(resultsIn, para);
        
        curMovieJumps = resultsOut.jumpDists;
        curMovieMeanJumps = resultsOut.meanJumpDists;
        curMovieAngles = resultsOut.angles;
        curMovieDistToRoiBorder = resultsOut.distToRoiBorder;
        meanJumpDistMoviewise(mIdx) = resultsOut.meanJumpDistMoviewise;
        
    else
        %------------User selected a subregion-----------------------------
        
        %List of logicals that indicates which tracks belong to the
        %selected sub-region(s)
        tracksSubRoi = any(curMovieResults.tracksSubRoi == regionNum-1,2);
        
        %List of logicals that indicates which non-linked spots belong to
        %the selected sub-region(s)
        nonLinkedSpotsSubRoi = any(curMovieResults.nonLinkedSpotsSubRoi == regionNum-1,2);
        
        %Get the number of tracks in the selected sub-region(s)
        nTracks(mIdx) = sum(tracksSubRoi);
        
        %Get the number of non-linked spots in the selected sub-region(s)
        nNonLinkedSpots(mIdx) = sum(nonLinkedSpotsSubRoi);
        
        %Use only tracks within selected subRoi(s)
        startEndFrameOfTracks{mIdx} = startEndFrameOfTracks{mIdx}(tracksSubRoi,:);        
        trackLengths{mIdx} = trackLengths{mIdx}(tracksSubRoi);
        meanTrackLength(mIdx) = mean(trackLengths{mIdx});
        
        if ~isempty(curMovieDistToRoiBorder)
            curMovieDistToRoiBorder = curMovieDistToRoiBorder(tracksSubRoi);
        end
        
        
        %Array containing the ROI sizes of relevant movies
        if numel(regionNum) > 1
            %If more than one sub-region is selected we don't know how to
            %handle overlapping regions so we just set the roi size to 1
            roiSize(mIdx) = 1;
        elseif nRegionsInMovie < regionNum
            %Current movie has less regions than the selected region number
            roiSize(mIdx) = 0;
        else
            %Get roi size of the selected sub-region number of the current
            %movie
            roiSize(mIdx) = curMovieResults.subRegionResults(regionNum).roiSize.*para.pixelsize^2;
        end
        
        %-------Adjust results for user restrictions-----------------------
        
        %Adjust all values based on the user input that restricts our tracks
        %(eg. "nJumps to consider", "Remove jumps over gap frames etc.)   
        
        resultsIn.tracks = curMovieTracks(tracksSubRoi);
        resultsIn.jumpDists = curMovieJumps(tracksSubRoi);
        resultsIn.meanJumpDists = curMovieMeanJumps(tracksSubRoi);
        resultsIn.angles = curMovieAngles(tracksSubRoi);
        resultsIn.distToRoiBorder = curMovieDistToRoiBorder;
        
        resultsOut = adjust_jumps_and_angles(resultsIn, para);
        
        curMovieJumps = resultsOut.jumpDists;
        curMovieMeanJumps = resultsOut.meanJumpDists;
        curMovieAngles = resultsOut.angles;
        curMovieDistToRoiBorder = resultsOut.distToRoiBorder;
        meanJumpDistMoviewise(mIdx) = resultsOut.meanJumpDistMoviewise;
                
        %Calculate average number of spots and tracks per frame
        meanTracksPerFrame(mIdx) = resultsOut.nSpotsInTracks/nFramesAnalyzed;
        meanSpotsPerFrame(mIdx) = (resultsOut.nSpotsInTracks + nNonLinkedSpots(mIdx))/nFramesAnalyzed;
        
        %------------------------------------------------------------------

    end
    
    %Catenate angles of all tracks in current movie
    curMovieAngles = vertcat(curMovieAngles{:});
    angles{mIdx} = curMovieAngles;
    
    %Catenate jump distances of all tracks in current movie
    jumpDistances{mIdx} = vertcat(curMovieJumps{:});  
    nJumps(mIdx) = numel(jumpDistances{mIdx});
    
    %Catenate distances of all tracks to region border in current movie
    if ~isempty(curMovieDistToRoiBorder)
%         distToRoiBorder{mIdx} = vertcat(curMovieDistToRoiBorder{:});
        distToRoiBorder{mIdx} = curMovieDistToRoiBorder(:);
%     else
%         distToRoiBorder{mIdx} = 0;
    end
    
    %Catenate the mean jump distances of all tracks in current movie
    meanJumpDists{mIdx} = curMovieMeanJumps(:);
    
    %Check if user clicked on "Calculate confinement radii" so that these data exist
    if isfield(curMovieResults, 'alphaValues')
        %Get the alpha valuesm confinement radii, meanjump distances and 
        %diffusion constants of the tracks in current movie
        curMovieAlphaValues = curMovieResults.alphaValues;
        curMovieConfRad = curMovieResults.confRad;
        curMovieMeanJumpDistConfRad = curMovieResults.meanJumpDistConfRad.*para.pixelsize;
        curMovieMsdDiffConst = curMovieResults.msdDiffConst./(frameCycleTimes(mIdx)/para.pixelsize^2);
        
        %Create a list of logicals indicating which alpha values are below
        %the user defined threshold and make sure that the tracks lies
        %within the sub-region(s) selected by the user
        belowAlphaThresAndInRoiIdx = (curMovieAlphaValues <= para.alphaThres) & tracksSubRoi;
        
        %Use the logical list to filter out tracks that do not match the
        %above criteria and adjust for the pixelsize
        alphaValues{mIdx} = curMovieAlphaValues(~isnan(curMovieAlphaValues) & belowAlphaThresAndInRoiIdx);
        confRad{mIdx} = curMovieConfRad(~isnan(curMovieConfRad) & belowAlphaThresAndInRoiIdx).*para.pixelsize;        
        curMovieMeanJumpDistConfRad = curMovieMeanJumpDistConfRad(belowAlphaThresAndInRoiIdx);    
        curMovieMsdDiffConst = curMovieMsdDiffConst(belowAlphaThresAndInRoiIdx);
        
        %Filter out NaNs
        meanJumpDistConfRad{mIdx} = curMovieMeanJumpDistConfRad(~isnan(curMovieMeanJumpDistConfRad));
        msdDiffConst{mIdx} = curMovieMsdDiffConst(~isnan(curMovieMsdDiffConst));
    else
        %User did not yet click on "Calculate confinement radii"
        alphaValues = {0};
        msdDiffConst = {0};
        confRad = {0};
        meanJumpDistConfRad = {0};
    end
    
end


%Calculate tracked fractions
[nLong, nShort, trackedFractions] = bf_analysis(startEndFrameOfTracks, nNonLinkedSpots, para);

results.distToRoiBorder = distToRoiBorder;
results.meanTrackLength = meanTrackLength;
results.meanSpotsPerFrame = meanSpotsPerFrame;
results.trackLengths = trackLengths;
results.angles = angles;
results.jumpDistances = jumpDistances;
results.meanJumpDists = meanJumpDists;
results.roiSize = roiSize;
results.meanTracksPerFrame = meanTracksPerFrame;
results.movieNames = movieNames;
results.trackingRadii = trackingRadii;
results.meanJumpDistMoviewise = meanJumpDistMoviewise;
results.frameCycleTimes = frameCycleTimes;
results.meanJumpDistConfRad = meanJumpDistConfRad;
results.alphaValues = alphaValues;
results.msdDiffConst = msdDiffConst;
results.confRad = confRad;
results.movieNumbers = moviesIdx';
results.nJumps = nJumps;

results.batchName = '';

results.nAllEvents = nTracks + nNonLinkedSpots;
results.nNonLinkedSpots = nNonLinkedSpots;
results.nTracks = nTracks;
results.nShort = nShort;
results.nLong =  nLong;
results.trackedFractions = trackedFractions;

end


function resultsOut = adjust_jumps_and_angles(resultsIn, para)

%Adjust all values based on the user input that restricts our tracks
        %(eg. "nJumps to consider", "Remove jumps over gap frames etc.)   


tracks = resultsIn.tracks;
jumpDists = resultsIn.jumpDists;
meanJumpDists = resultsIn.meanJumpDists;
angles = resultsIn.angles;
distToRoiBorder = resultsIn.distToRoiBorder;
if ~isempty(distToRoiBorder)    
    distToRoiBorder = resultsIn.distToRoiBorder.*para.pixelsize;
end

%Counter for the total number of spots in tracks 
%(only needed for sub-region analysis)
nSpotsInTracks = 0;

%Iterate through all tracks
for trackIdx = 1:length(tracks)
    
    %Get jump distances of the current movie and multiply by the pixelsize
    curTrackJumps = jumpDists{trackIdx}.*para.pixelsize;
    
    %Increase counter for number of spots in tracks by the amount of spots
    %in the current track
    nSpotsInTracks = nSpotsInTracks + numel(curTrackJumps)+1;
    
    %Get the angles of the current track
    curTrackAngles = angles{trackIdx};
        
    %-----Check if gaps have to be removed in current track----------------
    if para.boolRemoveGaps
        %"Remove jumps over gap frames" chechbox is enabled
        
        %Get frames where track appears
        curTrackFrames = tracks{trackIdx}(:,1);
        %Find out which spots of the tracks appear in subsequent frames
        subseqFrames = diff(curTrackFrames) == 1;
        
        %Create a list of logicals indicating where gap frames happen
        hasGapFrames = ~all(subseqFrames);
        
        %Check if there are gap frames in the track
        if hasGapFrames
            removeGapsInCurTrack = true;
        else
            removeGapsInCurTrack = false;
        end
    else
        removeGapsInCurTrack = false;
    end
    
    %------------------ Adjust Angles--------------------------------------
    
    %Get the number of jump angles in current track
    nAnglesInCurTrack = numel(curTrackAngles);
    
    %Check if there are any angles at all
    if nAnglesInCurTrack > 0
        
        %Initialize list of logicals which angles should be taken into
        %account in the results
        anglesToUse = true(nAnglesInCurTrack,1);
        
                
        %-----Remove gap frames from angles--------------------------------
        if removeGapsInCurTrack
            %Use only angles where the two track segments making up the 
            %angle do not contain a gap frame
                        
            anglesToUse(subseqFrames(1:end-1)+subseqFrames(2:end) ~= 2) = 0;
        end
        
        %------Jump distances making up the angle--------------------------
        
        %Only use angles where the jump distances making up the angles are
        %higher or lower the threshold set by the user
        if para.minJumpDistForAngles > 0 || para.maxJumpDistForAngles < inf
            %User only wants to keep angles where the two track
            %segments making up the angle have a jump distance higher or
            %lower than the user defined threshold
            
            %Jumps must be between the two thresholds
            jumpsToUseForAngles = curTrackJumps >= para.minJumpDistForAngles & curTrackJumps <= para.maxJumpDistForAngles;
            
            anglesWithAllowedJumps = (jumpsToUseForAngles(1:end-1) & jumpsToUseForAngles(2:end));
            
            %Throw out angles where the two jumps making up the angle are
            %outside the desired jump distance range
            anglesToUse(~anglesWithAllowedJumps) = 0;
        end
                
        
        %---------Adjust angles for nJumpsToConsider-----------------------
        if para.nJumpsToConsider-1 < numel(curTrackAngles)
            %Use only angles spanned up by the first nJumpsToConsider jumps
            anglesToUse(para.nJumpsToConsider:end) = 0;
        end
                
        %Get the angles that match all above criteria
        curTrackAngles = curTrackAngles(anglesToUse);
        
        
%         TODO in future: this can be used to plot the mean jump distances making up the angle vs. the fold anisotropy
%         curTrackAnglesMeanDisp = (curTrackJumps(1:end-1)+curTrackJumps(2:end))./2;
%         curTrackAnglesMeanDisp = curTrackAnglesMeanDisp(anglesToUse);
        
        %Save the angles of the current track
        angles{trackIdx} = curTrackAngles;        
    end
    
    
    %-----Adjust jump distances-----------------------------------
            
    if removeGapsInCurTrack
        %Remove jumps over gaps from the jump distances and Roi border distances
        curTrackJumps = curTrackJumps(subseqFrames);
    end
        
    %Check if track is longer than nJumpsToConsider
    adjustNJumpsInCurTrack = para.nJumpsToConsider < numel(curTrackJumps);
    
    if adjustNJumpsInCurTrack
        %Adjust jump dists for nJumpsToConside
        curTrackJumps = curTrackJumps(1:para.nJumpsToConsider);
    end
    
    %Calculate new mean jump distance
    if removeGapsInCurTrack || adjustNJumpsInCurTrack
        meanJumpDists(trackIdx) = mean(curTrackJumps);
    end
    
    %Save jump distances of the current track
    jumpDists{trackIdx} = curTrackJumps;
        
end

%Calculate average jump distance of all jumps in current movie
meanJumpDistMoviewise = mean(vertcat(jumpDists{:}));
meanJumpDistMoviewise(isnan(meanJumpDistMoviewise)) = 0;

resultsOut.meanJumpDistMoviewise = meanJumpDistMoviewise;
resultsOut.jumpDists = jumpDists;
resultsOut.meanJumpDists = meanJumpDists;
resultsOut.angles = angles;
resultsOut.distToRoiBorder = distToRoiBorder;
resultsOut.nSpotsInTracks = nSpotsInTracks;


end


function [nLong, nShort, trackedFractions] = bf_analysis(startEndFrameOfTracks, nNonLinkedSpots, para)

%Analyze tracked fractions 

%Initialize variables
nMovies = length(startEndFrameOfTracks); %Amount of movies in this Batch
nShort = zeros(nMovies,1); %Amount of tracks that last shorter than nDarkTimesThres
nLong = zeros(nMovies,1); %Amount of tracks that last longer than nDarkTimesThres
curMovieNNonLinkedSpots = zeros(nMovies,1); %Amount of non-linked spots
longVsAllEventsMoviewise = zeros(nMovies,1); %Fraction of long tracks with respect to all events: nLong/(nShort+nLong+nNonLinkedSpots)
allTracksVsAllEventsMoviewise = zeros(nMovies,1); %Fraction of all tracks with respect to all events: (nShort+nLong)/(nShort+nLong+nNonLinkedSpots)
longVsAllTracksMoviewise = zeros(nMovies,1); %Fraction of long tracks with respect to all tracks: (nLong)/(nShort+nLong)
shortVsAllEventsMoviewise = zeros(nMovies,1); %Fraction of short tracks with respect to all tracks: (nShort)/(nShort+nLong)

for movieIdx = 1:nMovies
    
    curMovieStartEndFrameOfTracks = startEndFrameOfTracks{movieIdx};
    curMovieNNonLinkedSpots(movieIdx) = nNonLinkedSpots(movieIdx);
    
    %Calculate the amount of dark times a track survives
    nDarkTimes = zeros(size(curMovieStartEndFrameOfTracks,1),1);
    for trackIdx = 1:size(curMovieStartEndFrameOfTracks,1)
        nDarkTimes(trackIdx) = floor((curMovieStartEndFrameOfTracks(trackIdx,2)-1)/para.nBrightFrames)-floor((curMovieStartEndFrameOfTracks(trackIdx,1)-1)/para.nBrightFrames);
    end
    
    %Get amount of tracks which survive less dark times than nDarkTimes
    nShort(movieIdx) = sum(nDarkTimes < para.nDarkForLong);
    %Get amount of tracks which survive at least nDarkTimes
    nLong(movieIdx) = sum(nDarkTimes >= para.nDarkForLong);
    %Get amount of all tracks
    nAllEvents = nShort(movieIdx) + nLong(movieIdx) + curMovieNNonLinkedSpots(movieIdx);
    
    %Save long immobile fractions of each movie in an array for a moviewise evaluation
    %Convert NaN to zero if divided by zero
    longVsAllEventsMoviewise(movieIdx) = max(0,nLong(movieIdx)/nAllEvents);
    allTracksVsAllEventsMoviewise(movieIdx) = max(0,(nLong(movieIdx)+nShort(movieIdx))/nAllEvents);
    longVsAllTracksMoviewise(movieIdx) = max(0,nLong(movieIdx)/(nShort(movieIdx)+nLong(movieIdx)));    
    shortVsAllEventsMoviewise(movieIdx) = max(0,nShort(movieIdx)/nAllEvents);
    
end

%Get mean values of the moviewise fractions
longVsAllEventsMean = mean(longVsAllEventsMoviewise);
allTracksVsAllEventsMean = mean(allTracksVsAllEventsMoviewise);
longVsAllTracksMean = mean(longVsAllTracksMoviewise);
shortVsAllEventsMean = mean(shortVsAllEventsMoviewise);

%Get standard error of the mean values of the moviewise fractions
longVsAllEventsStd = std(longVsAllEventsMoviewise)/sqrt(nMovies);
allTracksVsAllEventsStd = std(allTracksVsAllEventsMoviewise)/sqrt(nMovies);
longVsAllTracksStd = std(longVsAllTracksMoviewise)/sqrt(nMovies);
shortVsAllEventsStd = std(shortVsAllEventsMoviewise)/sqrt(nMovies);

%Calculate total amounts of different track classes for a pooled evaluation
nLongPooled = sum(nLong);
nShortPooled = sum(nShort);
nNonLinkedSpotsPooled = sum(curMovieNNonLinkedSpots);
ntracksPooled = nLongPooled + nShortPooled;
nMoleculesTotalSum = nNonLinkedSpotsPooled +  ntracksPooled;

%Get pooled immobile fractions: here all tracks are pooled and not split by
%movies as in the case of the moviewise factions
longVsAllEventsPooled = nLongPooled/nMoleculesTotalSum;
allTracksVsAllEventsPooled = ntracksPooled/nMoleculesTotalSum;
longVsAllTracksPooled = nLongPooled/ntracksPooled;
shortVsAllEventsPooled = nShortPooled/nMoleculesTotalSum;

%Linear error propagation using counting error as error source
errorLongVsAllEventsPooled = sqrt(nLongPooled)/nMoleculesTotalSum...
    +sqrt(nMoleculesTotalSum)*nLongPooled/(nMoleculesTotalSum)^2;

errorAllTracksVsAllEventsPooled = sqrt(ntracksPooled)/(nMoleculesTotalSum)...
    +sqrt(nMoleculesTotalSum)*ntracksPooled/(nMoleculesTotalSum)^2;


errorLongVsAllTracksPooled = sqrt(nLongPooled)/ntracksPooled...
    +sqrt(ntracksPooled)*nLongPooled/(ntracksPooled)^2;

errorShortVsAllEventsPooled = sqrt(nShortPooled)/nMoleculesTotalSum...
    +sqrt(nMoleculesTotalSum)*nShortPooled/(nMoleculesTotalSum)^2;



trackedFractions = struct(...
    'longVsAllEventsMoviewise',longVsAllEventsMoviewise,'longVsAllTracksMoviewise',longVsAllTracksMoviewise,...
    'allTracksVsAllEventsMoviewise',allTracksVsAllEventsMoviewise,'longVsAllEventsMean',longVsAllEventsMean,...
    'allTracksVsAllEventsMean',allTracksVsAllEventsMean,'longVsAllTracksMean',longVsAllTracksMean,...
    'longVsAllEventsStd',longVsAllEventsStd,'allTracksVsAllEventsStd',allTracksVsAllEventsStd,'longVsAllTracksStd',longVsAllTracksStd,...
    'longVsAllEventsPooled',longVsAllEventsPooled,'longVsAllTracksPooled',longVsAllTracksPooled,...
    'allTracksVsAllEventsPooled',allTracksVsAllEventsPooled,'errorAllTracksVsAllEventsPooled',errorAllTracksVsAllEventsPooled,...
    'errorLongVsAllEventsPooled',errorLongVsAllEventsPooled,'errorLongVsAllTracksPooled',errorLongVsAllTracksPooled,...
    'shortVsAllEventsMoviewise',shortVsAllEventsMoviewise,'shortVsAllEventsStd',shortVsAllEventsStd,...
    'shortVsAllEventsMean',shortVsAllEventsMean, 'shortVsAllEventsPooled',shortVsAllEventsPooled,...
    'errorShortVsAllEventsPooled',errorShortVsAllEventsPooled);

end