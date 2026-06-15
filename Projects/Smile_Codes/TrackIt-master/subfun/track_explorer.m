function track_explorer(stack, curMovieResults, curMovieFrameTime, trackID, trackItPixelSize)

%TrackIt tool that lets the user explore single-molecule tracks by displaying
%kymographs, intensities, jump distances, angles and mean squared
%displacement of a track.
%
%
%track_explorer(movieStack, curMovieResults, curMovieFrameTime, trackID, pixelSizeFieldValue)
%
%Input:
%   stack               -   3d-array of pixel values. 
%                                   1st dimension (column): y-coordinate of the image plane 
%                           2nd dimension (row):    x-coordinate of the image plane 
%                               3rd dimension: frame number
%   curMovieResults     -   Structure containing the TrackIt results of a single
%                           movie eg. curMovieResults = batch(1).results;
%   curMovieFrameTime   -   Frame cycle time of the movie
%   trackID             -   TrackId of the track that was selected in TrackIT
%   trackItPixelSize    -   Pixelsize in microns per pixel


%Initialize pixelsize and frame cycle time
pixelSize = 1;
frameCycleTime = 1;

%Create user Interface
ui = InitUI();

% Initialize images and plots
%Plot style
plotLineWidth = 1.5;
spotMarkerSize = 10;

%Create data for current trackID
results = createCurTrackData(stack, curMovieResults, trackID);

%Get number of tracks in movie
nTracks = curMovieResults.nTracks;

%Initialize plots and save handles in plotHandles structure
plotHandles = InitGraphs();

%Plot results into graphs
UpdateGraphs();


%Adjust Slider and frame text
nFrames = results.stackSize(3);
ui.sliderFrame.Max = nFrames;
ui.sliderFrame.SliderStep = [min(1,1/(nFrames-1)) 1/min((nFrames-1),10)];
ui.sliderFrame.Value = results.trackLength;


%% Functions
    function ui = InitUI()
        
        %% Create figure and axes
        ui.hFig = figure('Name','Track explorer',...
            'Units','normalized',...
            'Tag','track_explorer',...
            'Position',[0 .1 .7 .7]);
        
        posLeft = .02;
        sizeLeft = .45;
        posRight = .52;
        sizeRight = .28;
        
        ui.axTrackPlot     = axes(ui.hFig,'Position',[0.005 0.375 .32 0.62]);
        ui.axKymVer        = axes(ui.hFig,'Position',[posLeft 0.17 sizeLeft 0.12]);
        ui.axKymHor        = axes(ui.hFig,'Position',[posLeft 0.01 sizeLeft 0.12]);
        ui.axInt           = axes(ui.hFig,'Position',[posRight 0.7 sizeRight 0.25]);
        ui.axDis           = axes(ui.hFig,'Position',[posRight 0.38 sizeRight 0.25]);
        ui.axAngles        = polaraxes(ui.hFig,'Position',[.82 0.37 .16 0.25]);
        ui.axMSD           = axes(ui.hFig,'Position',[posRight 0.06 sizeRight 0.25]);
        
        %% Slider
        ui.textFrame       = uicontrol('Units','normalized',....
            'Position',[.01 .35 .4 .02],....
            'Style','Text',....
            'String','',....
            'HorizontalAlignment','Left');
        
        ui.sliderFrame     = uicontrol('Units','normalized',....
            'Position',[.01 .32 sizeLeft+.01 .02],....
            'Style','slider',....
            'Min',1,'Max',1,'Value',1,....
            'SliderStep',[1 1]);
        
        addlistener(ui.sliderFrame,'Value','PostSet',@(~,~)FrameSlider);
        
        %% Show frame or no image
        ui.btnGroupShowSpotImage = uibuttongroup(...
            'Units','normalized',...
            'Position', [.335  .9  .13 .06],...
            'Title','Track background',...
            'SelectionChangedFcn',@FrameSlider);
        
        uicontrol(ui.btnGroupShowSpotImage,...
            'Units','normalized',...
            'Position', [.05  .0   .5 .95],...
            'Style','radiobutton',...
            'String','Frame',...
            'Tag','Frame',...
            'HorizontalAlignment','Left');
        
        uicontrol(ui.btnGroupShowSpotImage,...
            'Units','normalized',...
            'Position', [.45  .0   .6 .95],...
            'Style','radiobutton',...
            'String','No image',...
            'Tag','Black',...
            'HorizontalAlignment','Left');
        
        %% Track coloring button group
        ui.btnGroupTrackColor = uibuttongroup(...
            'Units','normalized',...
            'Position', [.335  .78   .13 .10],...
            'Title','Track color',...
            'SelectionChangedFcn',@FrameSlider);
        
        uicontrol(ui.btnGroupTrackColor,...
            'Units','normalized',...
            'Position', [.05  .5   .5 .4],...
            'Style','radiobutton',...
            'String','Uniform',...
            'Tag','Uniform',...
            'HorizontalAlignment','Left');
        
        uicontrol(ui.btnGroupTrackColor,...
            'Units','normalized',...
            'Position', [.45  .5   .6 .4],...
            'Style','radiobutton',...
            'String','Time',...
            'Tag','Time',...
            'HorizontalAlignment','Left');
        
        uicontrol(ui.btnGroupTrackColor,...
            'Units','normalized',...
            'Position', [.05  .0   .5 .4],...
            'Style','radiobutton',...
            'String','Intensity',...
            'Tag','Intensity',...
            'HorizontalAlignment','Left');
        
        uicontrol(ui.btnGroupTrackColor,...
            'Units','normalized',...
            'Position', [.45  .0   .6 .4],...
            'Style','radiobutton',...
            'String','Jump distance',...
            'Tag','Jump distance',...
            'HorizontalAlignment','Left');
        
        %% Kymograph button group
        ui.btnGroupKymo = uibuttongroup(...
            'Units','normalized',...
            'Position', [.335  .65   .13 .10],...
            'Title','Kymograph projection window',...
            'SelectionChangedFcn',@(~,~)UpdateGraphs);
        
        uicontrol(ui.btnGroupKymo,...
            'Units','normalized',...
            'Position', [.05  .5   .95 .4],...
            'Style','radiobutton',...
            'String','Projection of whole track',...
            'Tag','TrackCentered',...
            'HorizontalAlignment','Left');
        
        uicontrol(ui.btnGroupKymo,...
            'Units','normalized',...
            'Position', [.05  .0   .95 .4],...
            'Style','radiobutton',...
            'String','Propagating with spot center',...
            'Tag','SpotCentered',...
            'HorizontalAlignment','Left');
        
        
        uicontrol(....
            'Units','normalized',....
            'Position',  [.335   .55  .09  .05],....
            'Style','Text',....
            'String','Number of frames to show after end of track',....
            'HorizontalAlignment','Left');
        
        ui.editFramesAfterBleach = uicontrol(....
            'Units','normalized',....
            'Position',  [.435  .56 .028  .035],....
            'Style','Edit',....
            'String',20,....
            'HorizontalAlignment','Right',....
            'Callback',@FramesAfterBleachChangedCB);
        
        
        %% Track ID selection
        panelTrackId                = uipanel('Position',[.82  .88 .16 .1],'Tag','panelTrackId'); %'Title','TrackID'
        
        uicontrol(panelTrackId,...
            'Units','normalized','Position',....
            [.52  .55   .45  .4],....
            'String','Next',....
            'Callback',@TrackIdCB);
        
        uicontrol(panelTrackId,....
            'Units','normalized',....
            'Position', [.025 .55   .45  .4],....
            'String','Previous',....
            'Callback',@TrackIdCB);
        
        uicontrol(panelTrackId,....
            'Units','normalized',....
            'Position',  [.1   .05  .20  .32],....
            'Style','Text',....
            'String','Track ID',....
            'HorizontalAlignment','Left');
        
        ui.editTrackId = uicontrol(panelTrackId,....
            'Units','normalized',....
            'Position',  [.52  .05 .20  .4],....
            'Style','Edit',....
            'String',num2str(trackID),....
            'HorizontalAlignment','Right',....
            'Tag','editTrackId');
        
        
        addlistener(ui.editTrackId,'String','PostSet',@TrackIdCB);
        
        uicontrol(panelTrackId,....
            'Units','normalized',....
            'Position',  [.75  .05  .20  .32],....
            'Style','Text',....
            'String',['/', num2str(curMovieResults.nTracks)],....
            'HorizontalAlignment','Left');
        
        %% Unit controls
        ui.btnGroupUnits = uibuttongroup(...
            'Units','normalized',...
            'Position', [.82  .77 .16 .1],...
            'Title','Units',...
            'SelectionChangedFcn',@UnitsChangedCB);
        
        uicontrol(ui.btnGroupUnits,...
            'Units','normalized',...
            'Position', [.05  .5  .5 .4],...
            'Style','radiobutton',...
            'String','pixels & frames',...
            'Tag','pxfr',...
            'HorizontalAlignment','Left');
        
        uicontrol(ui.btnGroupUnits,...
            'Units','normalized',...
            'Position', [.55  .5  .6 .4],...
            'Style','radiobutton',...
            'String','microns & sec',...
            'Tag', 'musec',...
            'HorizontalAlignment','Left');
        
        ui.textPixelsize = uicontrol(ui.btnGroupUnits,...
            'Units','normalized',...
            'Position', [.05  .1  .75 .25],...
            'Visible','off',...
            'Style','text',...
            'String','Pixelsize in microns per px:',...
            'HorizontalAlignment','Left');
        
        ui.editPixelsize = uicontrol(ui.btnGroupUnits,...
            'Units','normalized',...
            'Position', [.7  .1  .2 .3],...
            'Visible','off',...
            'FontSize',9,...
            'Style','Edit',...
            'String', trackItPixelSize,...
            'Callback',@UnitsChangedCB);
        
        %% Intensity button group
        ui.btnGroupIntensity = uibuttongroup(...
            'Units','normalized',...
            'Position', [.82  .7   .16 .06],...
            'Title','Spot intensity measure',...
            'SelectionChangedFcn',@(~,~)IntensityMeasureCB);
        
        uicontrol(ui.btnGroupIntensity,...
            'Units','normalized',...
            'Position', [.05  .0   .3 .95],...
            'Style','radiobutton',...
            'String','Mean',...
            'Tag','Mean',...
            'HorizontalAlignment','Left');
        
        uicontrol(ui.btnGroupIntensity,...
            'Units','normalized',...
            'Position', [.4  .0   .3 .95],...
            'Style','radiobutton',...
            'String','Max',...
            'Tag','Max',...
            'HorizontalAlignment','Left');
        
        uicontrol(ui.btnGroupIntensity,...
            'Units','normalized',...
            'Position', [.7  .01  .3 .95],...
            'Style','radiobutton',...
            'String','Fit peak',...
            'Tag','Fit',...
            'HorizontalAlignment','Left');
        
        %% MSD controls
        ui.btnGroupFitFun = uibuttongroup(...
            'Units','normalized',...
            'Position', [.82  .25 .16 .06],...
            'Title','Fitting function',...
            'SelectionChangedFcn',@MsdChangedCB);
        
        uicontrol(ui.btnGroupFitFun,...
            'Units','normalized',...
            'Position', [.05  .0   .45 .95],...
            'Style','radiobutton',...
            'String','Power law',...
            'Tag','MSD',...
            'HorizontalAlignment','Left');
        
        uicontrol(ui.btnGroupFitFun,...
            'Units','normalized',...
            'Position', [.55  .0   .45 .95],...
            'Style','radiobutton',...
            'String','Linear',...
            'Tag','Linear',...
            'HorizontalAlignment','Left');
        
        textHeight = .025;
        textWidth = .16;
        editWidth = .03;
        
        uicontrol('Units','normalized',...
            'Position',[.82  .19   textWidth textHeight],...
            'Style','Text',...
            'String','#Points to fit',...
            'HorizontalAlignment','Left');
        
        ui.editPointsToFit = uicontrol('Units','normalized',...
            'Position',[.82  .17   editWidth textHeight],...
            'Style','Edit',...
            'String','80 %',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.86  .17   editWidth textHeight],...
            'String','60 %',...
            'Tag','nPoints',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.89  .17   editWidth textHeight],...
            'String','70 %',...
            'Tag','nPoints',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.92  .17   editWidth textHeight],...
            'String','80 %',...
            'Tag','nPoints',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.95  .17   editWidth textHeight],...
            'String','90 %',...
            'Tag','nPoints',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.86  .14   editWidth textHeight],...
            'String','5',...
            'Tag','nPoints',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.89  .14   editWidth textHeight],...
            'String','10',...
            'Tag','nPoints',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.92  .14   editWidth textHeight],...
            'String','20',...
            'Tag','nPoints',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.95  .14   editWidth textHeight],...
            'String','100 %',...
            'Tag','nPoints',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.82  .09   textWidth textHeight],...
            'Style','Text',...
            'String','Maximum offset (px^2)',...
            'HorizontalAlignment','Left');
        
        ui.editOffset = uicontrol('Units','normalized',...
            'Position',[.82  .07   editWidth textHeight],...
            'Style','Edit',...
            'String','0',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.86  .07   editWidth textHeight],...
            'String','0',...
            'Tag','offset',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.89  .07   editWidth textHeight],...
            'String','0.5',...
            'Tag','offset',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.92  .07   editWidth textHeight],...
            'String','1',...
            'Tag','offset',...
            'HorizontalAlignment','Left',...
            'Callback',@MsdChangedCB);
        
        uicontrol('Units','normalized',...
            'Position',[.82  .01 .16 .04],...
            'String','Export data to Matlab workspace',...
            'Tag','offset',...
            'HorizontalAlignment','Left',...
            'Callback',@ExportCB);
        
        set(ui.hFig, 'WindowScrollWheelFcn', {@MouseWheelCB});
        
    end

    function plotHandles = InitGraphs()
        %Executed only once at the beginning to initialize all displayed plots
        
        %Plot Intensities
        plotHandles.trackIntHandle = plot(ui.axInt, NaN, NaN, '-');
        plotHandles.intTitle = title(ui.axInt,'');
        ylabel(ui.axInt,'Intensity (a.u.)')
        
        %Plot intensity after track loss
        hold(ui.axInt, 'on')
        plotHandles.afterBleachIntHandle = plot(ui.axInt,NaN, NaN,'k-');
        
        %Write mean intensity into plot
        plotHandles.intTextHandle = text(ui.axInt,.01,.95,'','Units','Normalized');
        
        %Plot marker for intensity in current frame
        plotHandles.intMarker = plot(ui.axInt, NaN, NaN,'Color','r','Marker','o','Linestyle','none');
        
        hold(ui.axInt, 'off')
        
        %Plot Distances
        plotHandles.jumpDistHandle = plot(ui.axDis, NaN, NaN,'.-');
        plotHandles.disTitle = title(ui.axDis,'Jump distance');
        
        %Write mean jump distance into plot
        hold(ui.axDis, 'on')
        plotHandles.disTextHandle = text(ui.axDis,.01,.95,'','Units','Normalized');
        
        %Plot marker for distance in current frame
        plotHandles.disMarker   = plot(ui.axDis, NaN, NaN,'Color','r','Marker','o','Linestyle','none');
        hold(ui.axDis, 'off')
        
        %Show image of spot
        plotHandles.imageHandle     = imshow(1,[],'Parent',ui.axTrackPlot);
        
        %Initialize confinement radius handle
        hold(ui.axTrackPlot, 'on')
        plotHandles.confRadiusHandle = viscircles(ui.axTrackPlot,[NaN, NaN], NaN, 'Linewidth', 1);
        hold(ui.axTrackPlot, 'off')
        
        %Plot track
        hold(ui.axTrackPlot, 'on')
        plotHandles.trackHandle     = plot(ui.axTrackPlot, NaN, NaN,'-','LineWidth',plotLineWidth,'MarkerSize',spotMarkerSize,'Color','g');
        hold(ui.axTrackPlot, 'off')
        
        %Plot Kymographs
        plotHandles.kymVerHandle = imshow(1,[], 'Parent',ui.axKymVer);
        plotHandles.kymHorHandle = imshow(1,[],'Parent',ui.axKymHor);
        
        %Plot red line of current frame and and blue spot indicatiing frame
        %and positions in vertical Kymographs
        hold(ui.axKymVer, 'on')
        plotHandles.kymVerLine = plot(ui.axKymVer, NaN, NaN,'red','linewidth',1.5);
        plotHandles.kymVerSpot = plot(ui.axKymVer, NaN, NaN,'b.','MarkerSize',10);
        colorbar(ui.axKymVer)
        title(ui.axKymVer,'kymograph vertical (yt)')
        hold(ui.axKymVer, 'off')
        
        %Plot red line of current frame and and blue spot indicatiing frame
        %and positions in horizontal Kymographs
        hold(ui.axKymHor, 'on')
        plotHandles.kymHorLine = plot(ui.axKymHor, NaN, NaN,'red','linewidth',1.5);
        plotHandles.kymHorSpot = plot(ui.axKymHor, NaN, NaN,'b.','MarkerSize',10);
        colorbar(ui.axKymHor)
        title(ui.axKymHor,'kymograph horizontal (xt)')
        hold(ui.axKymHor, 'off')
        
        %Plot MSD
        plotHandles.msd = plot(ui.axMSD,NaN,NaN,'.','Color', [0 0.4470 0.7410]);
        hold(ui.axMSD, 'on')
        plotHandles.msdFit = plot(ui.axMSD, NaN,NaN,'r');
        plotHandles.msdText = text(ui.axMSD,.01,.75,'','Units','Normalized');
        hold(ui.axMSD, 'off')
        title(ui.axMSD,'Mean squared displacement')
        
    end

    function results = createCurTrackData(imageStack, curMovieResults, trackID)
        halfKymoSize = 6;  %Half window size of movie cutout taken for plotting and analysis
        framesAfterBleach = str2double(ui.editFramesAfterBleach.String); %Amount of frames kymograph is plotted after track is bleached
        
        %Get first and last frame of track
        firstFrame = curMovieResults.tracks{trackID}(1,1);
        lastFrame = curMovieResults.tracks{trackID}(end,1);
        
        %Get borders for cutout of the movie stack
        yMin = min(round(curMovieResults.tracks{trackID}(:,2)))-halfKymoSize;
        yMax = max(round(curMovieResults.tracks{trackID}(:,2)))+halfKymoSize;
        xMin = min(round(curMovieResults.tracks{trackID}(:,3)))-halfKymoSize;
        xMax = max(round(curMovieResults.tracks{trackID}(:,3)))+halfKymoSize;
        
        %Calculate start & end positions where to cut the frame also
        %taking care if the track is close to an edge
        xStart = max(1,xMin);
        xEnd   = min(size(imageStack,1),xMax);
        yStart = max(1,yMin);
        yEnd   = min(size(imageStack,2),yMax);
        
        %Get maximum amount of frames that can be shown after track is lost
        zEnd = lastFrame + min(size(imageStack,3)-lastFrame, framesAfterBleach);
        
        %Cut out stack from original movie
        curTrackStack = imageStack(xStart:xEnd,yStart:yEnd,firstFrame:zEnd);
        
        %Save relative image position with respekt to the original movie
        relativeImagePos = [xStart, yStart, firstFrame];
        
        %% Get track data for current trackID
        track = curMovieResults.tracks{trackID};
        
        %Calculate new track coordinates and frames relative to the cutout stack
        track(:,2) = (track(:,2) - relativeImagePos(2)+1);
        track(:,3) = (track(:,3) - relativeImagePos(1)+1);
        track(:,1) = track(:,1) - relativeImagePos(3)+1;
        
        %Get jump distances and angles
        jumpDists = curMovieResults.jumpDistances{trackID}.*pixelSize;
        angles = curMovieResults.angles{trackID};
                
        %% Create projections (kymographs) of the whole track area    
        %Get stacksize and the tracklength in frames
        stackSize = size(curTrackStack);
        trackLength = track(end,1) - track(1,1) + 1;
        
        %Initialize kymographs
        kymoVerFull = zeros(stackSize(1),stackSize(3));
        kymoHorFull = zeros(stackSize(2),stackSize(3));
        
        for frameIdx = 1:stackSize(3)
            %Get current image
            curSpotImage = curTrackStack(:,:,frameIdx);
            
            %Insert projection of current image into kymograph
            kymoVerFull(:,frameIdx) = max(curSpotImage,[],2);
            kymoHorFull(:,frameIdx) = max(curSpotImage,[],1)';
        end
                
        %% Create spot centered kymographs and calculate intensities
        
        %Initialize kymographs and intensity variables
        kymoVerSpot = zeros(halfKymoSize*2+1,stackSize(3));
        kymoHorSpot = zeros(halfKymoSize*2+1,stackSize(3)); 
        
        intMean = zeros(stackSize(3),1);
        intMax = zeros(stackSize(3),1);
        intFitted = track(:,4);
        
        %Iterate trough frames
        for frameIdx = 1:stackSize(3)
            if frameIdx <= trackLength-1
                %Frame contains track
                trackFrame = find(track(:,1) <= frameIdx,1,'last');
                curYpos = round(track(trackFrame,2));
                curXpos = round(track(trackFrame,3));
            else
                %Frame after track is lost
                trackFrame = trackFrame+1;
            end
            
            %Calculate start & end positions where to cut the frame also
            %taking care if the track is close to an edge
            xStart = max(1,curXpos-halfKymoSize);
            xEnd = min(size(curTrackStack,1),curXpos+halfKymoSize);
            yStart = max(1,curYpos-halfKymoSize);
            yEnd = min(size(curTrackStack,2),curYpos+halfKymoSize);
            
            %Cut out spot image for the kymograph
            curSpotImage = curTrackStack(xStart:xEnd,yStart:yEnd,frameIdx);
            
            %Caclulate start & end position where to insert kymograph if
            %the track is close to the edge
            xStartSpotImage = abs(curXpos-xStart-halfKymoSize)+1;
            yStartSpotImage = abs(curYpos-yStart-halfKymoSize)+1;
            xEndSpotImage = abs(curXpos-xEnd-halfKymoSize)+1;
            yEndSpotImage = abs(curYpos-yEnd-halfKymoSize)+1;
            
            %Insert projection of current image into kymograph
            kymoVerSpot(xStartSpotImage:xEndSpotImage,frameIdx) = max(curSpotImage,[],2);
            kymoHorSpot(yStartSpotImage:yEndSpotImage,frameIdx) = max(curSpotImage,[],1)';
            
            %Get mean intensity of small spot area
            spotWindowHalfSize = 1; %Half window size for calculation of the mean spot intensity
            curSpotPeakImage1 = curSpotImage(halfKymoSize-spotWindowHalfSize:halfKymoSize+spotWindowHalfSize,halfKymoSize-spotWindowHalfSize:halfKymoSize+spotWindowHalfSize);
            intMean(frameIdx) = mean(curSpotPeakImage1(:));
            intMax(frameIdx) = max(curSpotPeakImage1(:));
        end
        
        %% Create Colormaps
        %Colormaps to plot the tracks according to the time, intensity or
        %displacement
        if trackLength > 3
            colMap = parula(64);
            timeCMap = interp1(linspace(1,trackLength,length(colMap)),colMap,(1:trackLength)'); % map color to y values
            timeCMap = uint8(timeCMap'*255); % need a 4xN uint8 array
            timeCMap(4,:) = 255; % last column is transparency
                        
            switch ui.btnGroupIntensity.SelectedObject.Tag
                case 'Mean'
                    intTrack = intMean;
                case 'Max'
                    intTrack = intMax;
                case 'Fit'
                    intTrack = intFitted;
            end
            
            intCMap = interp1(linspace(min(intTrack),max(intTrack),length(colMap)),colMap,intTrack); % map color to y values
            intCMap = uint8(intCMap'*255); % need a 4xN uint8 array
            intCMap(4,:) = 255; % last column is transparency
            
            dispCMap = interp1(linspace(min(jumpDists),max(jumpDists),length(colMap)),colMap,[jumpDists(1); jumpDists]); % map color to y values
            dispCMap = uint8(dispCMap'*255); % need a 4xN uint8 array
            dispCMap(4,:) = 255; % last column is transparency
            
        else
            %Track is too short to create color maps
            timeCMap = [];
            intCMap = [];
            dispCMap = [];
        end
        
        %% Save results in structure
        results = struct('curTrackStack',curTrackStack,'relativeImagePos',relativeImagePos,...
            'track',track,'jumpDists',jumpDists,'angles',angles,...
            'stackSize',stackSize,'trackLength',trackLength,...
            'kymoVerFull',kymoVerFull,'kymoHorFull',kymoHorFull,...
            'kymoVerSpot',kymoVerSpot,'kymoHorSpot',kymoHorSpot,...
            'intMean',intMean,'intMax',intMax,'intFitted',intFitted,'timeCMap',timeCMap,...
            'intCMap',intCMap,'dispCMap',dispCMap);
        
    end

    function TrackIdCB(src,~)

        %User changed track ID        
        previousTrackID = trackID;

        if ~ishandle(src)
            src = ui.editTrackId;
        end
        
        if strcmp(src.String,'Previous') && previousTrackID > 1 %User pressed 'previous' button
            trackID =  previousTrackID - 1;
        elseif strcmp(src.String,'Next') && previousTrackID < nTracks %User pressed 'next' button
            trackID = previousTrackID + 1;
        elseif str2double(src.String) <= nTracks && str2double(src.String) > 0 %User entered number into edit field
            trackID = str2double(src.String);
        end
        
        %Write new track Id into edit field
        ui.editTrackId.String = num2str(trackID);
        
        if previousTrackID ~= trackID
            %TrackID was changed so update all results and plots
            results = createCurTrackData(stack, curMovieResults, trackID);
            
            
            %Update slider and text
            ui.sliderFrame.Max = results.stackSize(3);
            ui.sliderFrame.SliderStep = [min(1,1/(results.stackSize(3)-1)) 1/min((results.stackSize(3)-1),10)];
            ui.sliderFrame.Value = results.trackLength;
            
            %Re-initialize graphs
            UpdateGraphs();
        end
        
    end

    function FramesAfterBleachChangedCB(~,~)
        %User changed number of frames after track loss so update results
        %and plots
        
        results = createCurTrackData(stack, curMovieResults, trackID);
        
        %Re-initialize graphs
        UpdateGraphs();
        
        %Update slider and text
        ui.sliderFrame.Max = results.stackSize(3);
        ui.sliderFrame.SliderStep = [min(1,1/(results.stackSize(3)-1)) 1/min((results.stackSize(3)-1),10)];
        if ui.sliderFrame.Value > ui.sliderFrame.Max
            ui.sliderFrame.Value = ui.sliderFrame.Max;
        end
        
    end

    function IntensityMeasureCB(~,~)
        %User changed number of frames after track loss so update results
        %and plots
        
        results = createCurTrackData(stack, curMovieResults, trackID);
        
        %Re-initialize graphs
        UpdateGraphs();
        
    end

    function UnitsChangedCB(~,~)
        %User changed units so update results and plots
        
        switch ui.btnGroupUnits.SelectedObject.Tag
            case 'pxfr'
                %Display results in pixels and frames
                ui.textPixelsize.Visible = 'off';
                ui.editPixelsize.Visible = 'off';
                frameCycleTime = 1;
                pixelSize = 1;
            case 'musec' %Display results in microns and seconds
                ui.textPixelsize.Visible = 'on';
                ui.editPixelsize.Visible = 'on';
                pixelSize = str2double(ui.editPixelsize.String);
                
                frameCycleTime = curMovieFrameTime/1000;
        end
        
        %Update results and plots
        results = createCurTrackData(stack, curMovieResults, trackID);
        UpdateGraphs();
        FrameSlider()
        
    end

    function FrameSlider(~,~)
        
        
        
        %Retrieve variables from results structure to speed up execution
        curTrackStack = results.curTrackStack;
        trackLength = results.trackLength;
        relativeImagePos = results.relativeImagePos;
        track = results.track;
        jumpDists = results.jumpDists;
        kymoVerSpot = results.kymoVerSpot;
        kymoHorSpot = results.kymoHorSpot;
        kymoVerFull = results.kymoVerFull;
        kymoHorFull = results.kymoHorFull;
        timeCMap = results.timeCMap;
        intCMap = results.intCMap;
        dispCMap = results.dispCMap;
        stackSize = results.stackSize;
        
        curFrame = round(ui.sliderFrame.Value);
        
        switch ui.btnGroupShowSpotImage.SelectedObject.Tag
            case 'Frame'
                %User wants to see the original movie frame
                plotHandles.imageHandle.CData = curTrackStack(:,:,curFrame);
            case 'Black'
                %User wants to have a black background
                plotHandles.imageHandle.CData = zeros(stackSize(1),stackSize(2));
        end
        
       
        
        %Find last frame that matches with a frame number of the track
        trackFrame = find(track(:,1) <= curFrame,1,'last');
        
        %Update plot of track positions
        set(plotHandles.trackHandle,'xdata',track(1:trackFrame, 2),'ydata', track(1:trackFrame, 3));

        %Get height of the currently shown kymograph
        switch ui.btnGroupKymo.SelectedObject.Tag
            case 'SpotCentered'
                kymoHeightVer = size(kymoVerSpot,1);
                kymoHeightHor = size(kymoHorSpot,1);
            case 'TrackCentered'
                kymoHeightVer = size(kymoVerFull,1);
                kymoHeightHor = size(kymoHorFull,1);
        end
        
        %Update red line showing the current frame number in the kympgraphs
        set(plotHandles.kymHorLine, 'xdata',[curFrame+0.5,curFrame+0.5],'ydata',[0.5,kymoHeightHor+.5])
        set(plotHandles.kymVerLine, 'xdata',[curFrame+0.5,curFrame+0.5],'ydata',[0.5,kymoHeightVer+.5])
        
        %Set track coloring of the plotted track
        if ~isempty(trackFrame)
            switch ui.btnGroupTrackColor.SelectedObject.Tag
                case 'Uniform' %Uniform track color (green)
                    set(plotHandles.trackHandle.Edge, 'ColorBinding','object', 'ColorData', uint8([0; 255; 0; 255]));
                    plotHandles.trackHandle.Marker = '.';
                case 'Time' %Track color changes from blue to yellow with increasing time
                    set(plotHandles.trackHandle.Edge, 'ColorBinding','interpolated', 'ColorData',timeCMap(:,1:trackFrame))
                    plotHandles.trackHandle.Marker = 'none';
                case 'Intensity' %Track color coded by the intensity in the corresponding frame
                    set(plotHandles.trackHandle.Edge, 'ColorBinding','interpolated', 'ColorData',intCMap(:,1:trackFrame))
                    plotHandles.trackHandle.Marker = 'none';
                case 'Jump distance' %Track color coded by the jump distance between the former and the current frame
                    set(plotHandles.trackHandle.Edge, 'ColorBinding','interpolated', 'ColorData',dispCMap(:,1:trackFrame))
                    plotHandles.trackHandle.Marker = 'none';
            end
        end
        
        %Check if current frame occurs in the framelist of the current track
        trackFrame = find(track(:,1) == curFrame);
        
        %Update red circle indicating the intensity of the current frame
        switch ui.btnGroupIntensity.SelectedObject.Tag
            case 'Mean'
                intMean = results.intMean;
                set(plotHandles.intMarker, 'xdata',curFrame.*frameCycleTime,'ydata',intMean(curFrame))
            case 'Max'
                intMax = results.intMax;
                set(plotHandles.intMarker, 'xdata',curFrame.*frameCycleTime,'ydata',intMax(curFrame))
            case 'Fit'
                
                if ~isempty(trackFrame)
                    intFitted = results.intFitted;
                    set(plotHandles.intMarker, 'xdata',curFrame.*frameCycleTime,'ydata',intFitted(trackFrame))
                else
                    %No spot detection in the current frame
                    set(plotHandles.intMarker, 'xdata',NaN,'ydata',NaN)
                end
        end
        
        if ~isempty(trackFrame)
            %Update red circle indicating the jump distance of the current frame
            set(plotHandles.disMarker, 'xdata',max(1,trackFrame-1).*frameCycleTime,'ydata',jumpDists(max(1,trackFrame-1)))
            
            %Update blue dot indicating the current spot position in the kymographs
            switch ui.btnGroupKymo.SelectedObject.Tag
                case 'SpotCentered'
                    set(plotHandles.kymHorSpot, 'xdata',curFrame+.5,'ydata',kymoHeightHor/2+0.5)
                    set(plotHandles.kymVerSpot, 'xdata',curFrame+.5,'ydata',kymoHeightVer/2+0.5)
                case 'TrackCentered'
                    set(plotHandles.kymHorSpot, 'xdata',curFrame+.5,'ydata',track(trackFrame, 2))
                    set(plotHandles.kymVerSpot, 'xdata',curFrame+.5,'ydata',track(trackFrame, 3))
            end
            
            ui.textFrame.String = sprintf('Frame in original movie: %d,   detection: \t %d/%d \t,   original coordinates (x / y): \t %.1f / %.1f \t ',...
                            curFrame+relativeImagePos(3)-1,curFrame, trackLength, relativeImagePos(2)-1+round(track(trackFrame,2),1), relativeImagePos(1)-1+round(track(trackFrame,3),1));
            
        else
            %No spot detection in the current frame
            set(plotHandles.kymHorSpot, 'xdata',NaN,'ydata',NaN)
            set(plotHandles.kymVerSpot, 'xdata',NaN,'ydata',NaN)
            set(plotHandles.disMarker, 'xdata',NaN,'ydata',NaN)
            set(plotHandles.intMarker, 'xdata',NaN,'ydata',NaN) 
            ui.textFrame.String = sprintf('Frame in original movie: %d', curFrame+relativeImagePos(3)-1);
        end
                
    end

    function UpdateGraphs()        
        %Retrieve variables from results structure to speed up execution
        track = results.track;
        jumpDists = results.jumpDists;
        stackSize = results.stackSize;
        trackLength = results.trackLength;
        curFrame = results.curTrackStack(:,:,trackLength);
        angles = results.angles;
                
        
        switch ui.btnGroupKymo.SelectedObject.Tag
            case 'SpotCentered'
                %Get kymograph centered around spot positions
                kymoVer = results.kymoVerSpot;
                kymoHor = results.kymoHorSpot;
            case 'TrackCentered'
                %Get kymograph showing projection of the whole track region
                kymoVer = results.kymoVerFull;
                kymoHor = results.kymoHorFull;
        end
        
        switch ui.btnGroupIntensity.SelectedObject.Tag
            case 'Mean'
                %Get mean intensity of spots
                intMean = results.intMean;
                intTrack = intMean(track(:,1));
                
                %Mean intensity at the position of the last spot after track is lost
                intAfterBleach = intMean(track(end,1)+1:end); 
                
                %Create vector with frame numbers after track loss
                frameNumsAfterBleach = track(end,1)+1:stackSize(3);
                intTitleText = 'Mean intensity (3x3 window around spot position)';
            case 'Max'
                %Get max intensity of spots
                intMax = results.intMax;
                intTrack = intMax(track(:,1));
                
                %Mean intensity at the position of the last spot after track is lost
                intAfterBleach = intMax(track(end,1)+1:end); 
                
                %Create vector with frame numbers after track loss
                frameNumsAfterBleach = track(end,1)+1:stackSize(3);
                intTitleText = 'Max intensity (3x3 window around spot position)';
            case 'Fit'
                %Get intensities from gaussin fits
                intFitted = results.intFitted;
                intTrack = intFitted(:);
                
                %After track loss intensity is not defined
                frameNumsAfterBleach = NaN;
                intAfterBleach = NaN;
                intTitleText = 'Peak of gaussian fit';
        end
        
        %Create axis labels and text field entries
        meanIntText = ['Mean intensity: ', num2str(mean(intTrack),'%.2f')];
        switch ui.btnGroupUnits.SelectedObject.Tag
            case 'pxfr'
                meanDispText = ['Mean jump distance: ', num2str(round_significant(mean(jumpDists),2,'round')), ' px'];
                ylabel(ui.axDis,'Jump distance (px)')
            case 'musec'
                meanDispText = ['Mean jump distance: ', num2str(round_significant(mean(jumpDists),2,'round')), ' µm'];
                ylabel(ui.axDis,'Jump distance (µm)')
        end
        
        %Update title of the intensity plot
        set(plotHandles.intTitle, 'String',intTitleText)
        %Update intensity plot
        set(plotHandles.trackIntHandle, 'xdata',results.track(:,1).*frameCycleTime,'ydata',intTrack)
        %Update intensity after track loss (black line)
        set(plotHandles.afterBleachIntHandle, 'xdata',frameNumsAfterBleach'.*frameCycleTime,'ydata',intAfterBleach)
        %Update text field with mean intensity
        set(plotHandles.intTextHandle, 'String',meanIntText)
        %Update jump distances plot
        set(plotHandles.jumpDistHandle, 'xdata',(1:numel(jumpDists)).*frameCycleTime,'ydata',results.jumpDists)
        %Update text field showing mean jump distance
        set(plotHandles.disTextHandle, 'String',meanDispText)
        
        
        %Update image of spot
        plotHandles.imageHandle.CData = curFrame;
        %Update axis limits
        ui.axTrackPlot.XLim = [0.5 stackSize(2)+0.5];
        ui.axTrackPlot.YLim = [0.5 stackSize(1)+0.5];
        %Update lookup table
        lowerThres = prctile(curFrame(:), 3);
        upperThres = prctile(curFrame(:), 99.8);
        ui.axTrackPlot.CLim = [lowerThres upperThres];
        
        %Update vertical kymograph
        plotHandles.kymVerHandle.CData = kymoVer;
        %Update axis limits
        kymoVerSize = size(kymoVer);
        ui.axKymVer.XLim = [0.5 kymoVerSize(2)+0.5];
        ui.axKymVer.YLim = [0.5 kymoVerSize(1)+0.5];
        %Update lookup table
        lowerThres = prctile(kymoVer(:), 3);
        upperThres = prctile(kymoVer(:), 99.8);
        ui.axKymVer.CLim = [lowerThres upperThres];
        
        %Update horizontal kymograph
        plotHandles.kymHorHandle.CData = kymoHor;
        %Update axis limits
        kymoHorSize = size(kymoHor);
        ui.axKymHor.XLim = [0.5 kymoHorSize(2)+0.5];
        ui.axKymHor.YLim = [0.5 kymoHorSize(1)+0.5];
        %Update lookup table
        lowerThres = prctile(kymoVer(:), 3);
        upperThres = prctile(kymoVer(:), 99.8);
        ui.axKymHor.CLim = [lowerThres upperThres];
        
        %Delete confinement radius circle
        delete(plotHandles.confRadiusHandle)
        
        %Plot Angles
        plotHandles.angleHandle = polarhistogram(ui.axAngles,angles,36);
        title(ui.axAngles,'Jump angles (deg)')
        %Plot MSD
        PlotMSD();
        FrameSlider()
    end

    function PlotMSD(~,~)
        
        track = results.track;
        
        %Plot MSD if track is longer than 4 frames
        if size(track,1) > 4
            
            %Get user settings: number of points or percentage to fit,
            %maximum Offset and msd or linear fit
            pointsToFit = ui.editPointsToFit.String;
            maxOffset = str2double(ui.editOffset.String);
            msdOrLinear = ui.btnGroupFitFun.SelectedObject.Tag;
            
            %Fit msd
            msdResult = msd_analysis({track},2,pointsToFit, maxOffset, msdOrLinear);
            
            %Adjust results to the unit configuration
            msdResult.msd = msdResult.msd{1}.*pixelSize^2;
            msdResult.msdDiffConst = msdResult.msdDiffConst./(frameCycleTime/pixelSize^2);
            msdResult.offset = msdResult.offset.*pixelSize^2;
            
            %Create x-Axis vektor for msd plot
            xMSD = (1:numel(msdResult.msd)).*frameCycleTime;
            
            %Create array containing the fit curve
            fittedPoints = 1:msdResult.nPointsFitted;            
            xFitCurve = fittedPoints.*frameCycleTime;
            yFitCurve = 4*frameCycleTime*msdResult.msdDiffConst*fittedPoints.^msdResult.alphaValues+msdResult.offset;            
            msdResult.msdFitCurveXY = [xFitCurve' yFitCurve'];
            
            %Save msd results in results structure
            results.msdResult = msdResult;
            
            %Creat axis labels and text field with results depending on the
            %selected units
            switch ui.btnGroupUnits.SelectedObject.Tag
                case 'pxfr' %Pixels and frames
                    msdResultsText = char(['alpha = ', num2str(round_significant(msdResult.alphaValues,2,'round'))],...
                        ['D = ', num2str(round_significant(msdResult.msdDiffConst,2,'round')), ' px^2/frame'],...
                        ['Offset = ', num2str(round_significant(msdResult.offset,2,'round')), ' px^2'],...
                        ['r_{confinement} = ', num2str(round_significant(msdResult.confRad,2,'round')), ' px']);
                    xlabel1 = 'frames';
                    ylabel1 = 'MSD (px^2)';
                case 'musec' %Microns and seconds
                    msdResultsText = char(['alpha = ', num2str(round_significant(msdResult.alphaValues,2,'round'))],...
                        ['D = ', num2str(round_significant(msdResult.msdDiffConst,2,'round')), ' µm^2/sec'],...
                        ['Offset = ', num2str(round_significant(msdResult.offset,2,'round')), ' µm^2'],...
                        ['r_{confinement} = ', num2str(round_significant(msdResult.confRad*pixelSize,2,'round')), ' µm']);
                    xlabel1 = 'time (sec)';
                    ylabel1 = 'MSD (µm^2)';
            end
            
            %Update msd plot
            set(plotHandles.msd, 'xdata',xMSD,'ydata',msdResult.msd)
            %Update red msd fit curve
            set(plotHandles.msdFit, 'xdata',xFitCurve,'ydata',yFitCurve)
            %Update text field containing fit results
            set(plotHandles.msdText, 'String',msdResultsText)
            %Update axis labels
            xlabel(ui.axMSD,xlabel1)
            ylabel(ui.axMSD,ylabel1)
            
            %Plot confinement radius into the track plot if confinement
            %radius fit was succesful
            if ~isnan(msdResult.confRad)
                hold(ui.axTrackPlot, 'on')
                plotHandles.confRadiusHandle = viscircles(ui.axTrackPlot,[mean(track(:,2)), mean(track(:,3))],msdResult.confRad, 'Linewidth', 1);
                hold(ui.axTrackPlot, 'off')
            end
        else
            set(plotHandles.msd, 'xdata',NaN,'ydata',NaN)
            set(plotHandles.msdFit, 'xdata',NaN,'ydata',NaN)
            set(plotHandles.msdText, 'String','')
        end
        
    end

    function MsdChangedCB(src,~)
        %User changed settings for msd analysis so redo the analysis
        
        delete(plotHandles.confRadiusHandle)
        
        switch src.Tag
            case 'nPoints'
                ui.editPointsToFit.String = src.String;
            case 'offset'
                ui.editOffset.String = src.String;
        end
        
        PlotMSD();
    end

    function ExportCB(~,~)
        %User clicked the export to Matlab workspace button
        
        results.msdResult.confRad = results.msdResult.confRad*pixelSize;
        trackExplorerData = struct('track', curMovieResults.tracks{trackID},...
            'jumpDistances', results.jumpDists,...
            'angles', results.angles,'msdResults', results.msdResult,...
            'kymoVerticalWholeTrack', results.kymoVerFull,...
            'kymoHorizontalWholeTrack', results.kymoHorFull,...
            'kymoVerticalSpotCentered', results.kymoVerSpot,...
            'kymoHorizontalSpotCentered', results.kymoHorSpot);
        
        assignin('base', 'trackExplorerData', trackExplorerData)
        
    end

    function MouseWheelCB(~,callbackdata)
        %Executed by mouse wheel for scrolling through frames
        curFrame = round(ui.sliderFrame.Value);
        %Take care that the frame number stays inside the frame range
        if curFrame + callbackdata.VerticalScrollCount > ui.sliderFrame.Max
            ui.sliderFrame.Value = ui.sliderFrame.Max;
        elseif curFrame + callbackdata.VerticalScrollCount < 1
            ui.sliderFrame.Value = ui.sliderFrame.Min;
        else
            ui.sliderFrame.Value = curFrame + callbackdata.VerticalScrollCount;
        end
    end


end