function data_analysis_tool(batches,currentBatchPath,curMovieIndex,pixelSize)

%
%Tool to analyze TrackIt batch files for a multitude of paramaters
%including diffusion analysis, bound fractions, track lengths etc.). Can be
%started either from within TrackIt via "Analysis" -> "Tracking data
%analysis" or by directly executing the data_analysis_tool function.
%
%
%Usage:
%data_analysis_tool(batches,currentBatchPath,curMovieIndex,pixelSize) or data_analysis_tool()
%
%
% Input: (function can also be called without inputs)
%     batches           -   initialize cell array of batches with the current
%                           batch analyzed in TrackIt
%     currentBatchPath  -   opening path for the "load batch file" dialog
%     curMovieIndex     -   Number of the movie currently visible in TrackIt main window
%     pixelSize         -   pixelsize in microns per pixel
%                   

if nargin == 4 && ~isempty(batches{1}(1).movieInfo.fileName)
    %Tracking data analysis tool was called from inside TrackIt
    
    if pixelSize == 0
       pixelSize = 1; 
    end
    
    %Get list of frame cycle times of each movie of the current batch
    frameCycleTimeMovieList = {zeros(length(batches{1}),1)};
    
    for cMovieIdx = 1:length(batches{1})
        frameCycleTimeMovieList{1}(cMovieIdx) = batches{1}(cMovieIdx).movieInfo.frameCycleTime;
    end
    
    %Get list of unique frame cycle times of the curernt batch
    frameCycleTimesList = {unique(frameCycleTimeMovieList{1})};
    
    %Create user interface
    ui = createHistogramUI();
    
    %Set name of TrackIt batch in the batch file list
    ui.popBatchSel.String = {'1: Current batch'};
           
     
else
    %Tracking data analysis tool was opended without TrackIt
    pixelSize = 1;
    
    %Add all subfolders to the matlab path
    mainFolder = fileparts(which(mfilename));
    addpath(genpath(mainFolder));
    
    %Opening path for the "load batch file" dialog
    currentBatchPath = '';  
    
    curMovieIndex = 1; %Initialize index of movie currently analyzed
    frameCycleTimeMovieList = {}; %Initialize list of frame cycle times of all movies
    frameCycleTimesList = {}; %Initialize list of unique frame cycle times
    batches = {}; %Initialize cell array of batch files
    ui = createHistogramUI(); 
end

%Initialize results structure
results = InitResults();

%Initialize structure containing values plotted in cetral graph
currentPlotValues = struct;

%Create results
BatchSelectionCB()

%Initialize plotStyle
plotStyle = 'Histogram';

%Initialize results and user interface

    function results = InitResults()
        %Initialzize results structure
        
        results.batchName = '';         %Name of the batch file
        results.movieNames = '';        %Names of movies in the current batch
        results.movieNumbers = -1;      %Movie numbers of the movies with selected frame cycle time
        results.frameCycleTimes = -1;   %List of frame cycle times of the movies in the current batch
        results.trackingRadii = -1;     %List of the tracking radii used in each movie of the current batch
         
        results.trackLengths = -1;      %Array containing the duration of all tracks in all movies of the current batch
        results.meanTrackLength = -1;   %Array containing the average track duration of all tracks in all movies of the current batch
        
        results.angles = -1;            %Array containing the angles between jumps within all tracks in all movies of the current batch
        results.jumpDistances = -1;     %Array containing the distances between jumps within all tracks in all movies of the current batch
        results.meanJumpDists = -1;     %Array containing the average jump distance of each track in all movies of the current batch
        results.meanJumpDistMoviewise = -1;%Array containing the average jump distance of each movie of the current batch
        results. nJumps = -1;           %Array containing the number of jumps in each movie of the current batch
        
        results.roiSize = -1;           %Array containing the sizes of the region of interest in all movies of the current batch
        results.meanTracksPerFrame = -1;%Array containing the average number of tracks per frame for each movie of the current batch
        results.meanSpotsPerFrame = -1; %Array containing the average number of spots per frame for each movie of the current batch
        
        results.alphaValues = -1;       %Array containing the alpha values from the msd fit of all tracks of all movies of the current batch
        results.msdDiffConst = -1;      %Array containing the diffusion coefficient calculated from the msd fit of all tracks of all movies of the current batch
        results.confRad = -1;           %Array containing the confinement radius calculated from the msd fit of all tracks of all movies of the current batch
        results.meanJumpDistConfRad = -1;%Array containing the mean jump distance of the tracks where a confinement radius was calculated. Used to plot confinement radius vs. mean jump distance
                
        results.nTracks = -1;           %Array containing the number of tracks in each movie of the current batch
        results.nShort = -1;            %Array containing the number of short tracks in each movie of the current batch (threshold is defined in ui)
        results.nLong = -1;             %Array containing the number of long tracks in each movie of the current batch (threshold is defined in ui)
        results.nNonLinkedSpots = -1;   %Array containing the number of non-linked spots each movie of the current batch
        results.nAllEvents = -1;        %Array containing the number of all events in each movie of the current batch (nTracks + nNonLinkedSpots)
        results.trackedFractions = struct;%Structure array containing the results of tracked fractions
        
        results.distToRoiBorder = -1;   %Array containing the minimum distance of the average track position from the region of interest
    end

    function ui = createHistogramUI()
        %%Figure and axes
        ui.f   = figure('Units','normalized',...
            'Position',[0.15 0.15 .77 .73],...
            'Name','Tracking data analysis',...
            'DefaultAxesFontSize',12,...
            'CloseRequestFcn',@(~,~)CloseHistogram);
        
        ui.pax  = polaraxes(ui.f,...
            'Units','normalized',...
            'visible','off',...
            'Position',[0.31 0.1 0.495 0.85]);
        
        ui.ax  = axes(ui.f,...
            'Units','normalized',...
            'Position',[0.315 0.1 0.485 0.85]);
        
 
        ui.hist = gobjects(1);
        
        %% Data selection
        
        ui.panSel = uipanel(ui.f,'Position',[0.01 0.01 0.25 0.98]);
        
        ui.btnLoadBatchFiles  = uicontrol('Parent',ui.panSel,...
            'Units','normalized',...
            'FontSize',8,...
            'Position',[0.05 0.94 .44 0.05],...
            'String','Load batch .mat file(s)',... <html><br>
            'Callback',@LoadBatchFilesCB);
        
        ui.btnRemoveBatchFiles  = uicontrol('Parent',ui.panSel,...
            'Units','normalized',...
            'FontSize',8,...
            'Position',[0.51 0.94 .44 0.05],...
            'String','Remove selected file(s)',... <html><br>
            'Callback',@RemoveBatchFilesCB);
        
        ui.popBatchSel = uicontrol('Parent',ui.panSel,...
            'Style','Listbox',...
            'Units','normalized',...
            'Position',[0.05 0.685 0.9 0.25],...
            'FontSize',9,...
            'Max',2,'Min',0,...
            'Callback',@(~,~)BatchSelectionCB);
        
        %Tl selection listbox
        ui.popTlSel = uicontrol('Parent',ui.panSel,...
            'Style','Listbox',...
            'Units','normalized',...
            'Position',[0.05 0.48 0.45 0.18],...
            'FontSize',9,...
            'String', {'Single movie', 'All movies'},...
            'Max',2,'Min',0,...
            'Callback',@(~,~)TlSelectionCB);
                
         % Movie Selection panel
        ui.panelMovieSel = uipanel(ui.panSel,...
            'Position',[0.52 0.58 0.45 0.09],'Title','Movie number'); %'BorderType','none'
        
        ui.btnNextMovie = uicontrol(ui.panelMovieSel,...
            'Units','normalized',...
            'Position', [.52  .45   .45 .5],...
            'String','Next',...
            'Callback',@MovieNumberCB);
        
        ui.btnPreviousMovie = uicontrol(ui.panelMovieSel,...
            'Units','normalized',...
            'Position', [.01  .45   .45 .5],...
            'String','Previous',...
            'Callback',@MovieNumberCB);
        ui.textMovie = uicontrol(ui.panelMovieSel,...
            'Units','normalized',...
            'Position', [.05  .06   .3 .3],...
            'Style','Text',...
            'String','Movie',...
            'HorizontalAlignment','Left');
        ui.editMovie = uicontrol(ui.panelMovieSel,...
            'Units','normalized',...
            'Position', [.4  .0  .25 .4],...
            'Style','Edit',...
            'String',curMovieIndex,...
            'HorizontalAlignment','Right',...
            'Callback',@MovieNumberCB);
        ui.textMovie2 = uicontrol(ui.panelMovieSel,...
            'Units','normalized',...
            'Position', [.69  .06  .1 .3],...
            'Style','Text',...
            'String','/',...
            'HorizontalAlignment','Left');
        ui.textMovie3 = uicontrol(ui.panelMovieSel,...
            'Units','normalized',...
            'Position', [.75  .06   .3 .3],...
            'Style','Text',...
            'String',1,...
            'HorizontalAlignment','Left');
        
        %Sub-regions (experimental)
        ui.popRegionSel = uicontrol('Parent',ui.panSel,...
            'Style','Listbox',...
            'Units','normalized',...
            'Position',[0.52 0.48 0.45 0.091],...
            'FontSize',9,...
            'String', {'All regions'},...
            'Max',2,'Min',0,...
            'Visible','on',... 
            'Callback',@(~,~)TlSelectionCB);
        
        %% Parameters selection tab group
        ui.tabGroupParam = uitabgroup('Parent',ui.panSel,...
            'Units','normalized',...
            'Position',[0.05 0.21 0.93 0.26],...
            'SelectionChangedFcn',@(~,~)PlotHistogramCB);
        
        ui.tab2 = uitab(ui.tabGroupParam,...
            'Title','Mobility');
        
        ui.popMobility = uicontrol('Parent',ui.tab2,...
            'Style','Listbox',...
            'Units','normalized',...
            'Position',[0.0 0.0 0.99 0.99],...
            'FontSize',9,...
            'String', {'Jump distances', 'Cumulative jump distances',...
            'Diffusion parameters','Mean jump distances',...
            'Cumulative mean jump distances','Mean jump distance per movie ','Jump angles',...
            'Confinement radius', 'Confinement radius vs. mean jump dist.',...
            'Alpha values from MSD fit', 'Diffusion constants from MSD fit'},...
            'Callback',@(~,~)PlotHistogramCB);
        
        ui.tab1 = uitab(ui.tabGroupParam,...
            'Title','Tracked fractions');
        
        ui.popTrackedFraction = uicontrol('Parent',ui.tab1,...
            'Style','Listbox',...
            'Units','normalized',...
            'Position',[0.0 0.0 0.99 0.99],...
            'FontSize',9,...
            'String', {'Tracks vs. all events',...
            'Long tracks vs. all events',...
            'Short tracks vs. all events',...
            'Long tracks vs. long + short tracks',...
            'No. of tracks','No. of non-linked spots'...
            'No. of all events (tracks + non-linked)',...
            'No. of long tracks', 'No. of short tracks'},...
            'Callback',@(~,~)PlotHistogramCB);
        
        ui.tab3 = uitab(ui.tabGroupParam,...
            'Title','Statistics');
        
        ui.popStatistics = uicontrol('Parent',ui.tab3,...
            'Style','Listbox',...
            'Units','normalized',...
            'Position',[0.0 0.0 0.99 0.99],...
            'FontSize',9,...
            'String', {'Track lengths','Cumulative track lengths',...
             'Avg. track length', 'Avg. no. of tracks per frame',...
             'Avg. no. of spots per frame','No. of jumps','ROI size',...
             'Distance to ROI border', 'Dist. to ROI border vs. mean jump dist.'},...
            'Callback',@(~,~)PlotHistogramCB);
               
        %% Statistics overview
        
        ui.txtStats  = uicontrol('Parent',ui.panSel,...
            'Style','text',...
            'Units','normalized',...
            'FontSize',10,...
            'Position',[0.01 0.14 1 0.06],...
            'String','Statistics overview');
        
        
        ui.tableStatistics = uitable(ui.panSel,'Units','normalized',...
            'Position',[0.01 0.01 .98 .155],...
            'ColumnName',{''},...            
            'RowName',{},...
            'Data',{'#movies';'#tracks';'#non-linked spots';'#all events'},...
            'ColumnEditable',[false,false,false,false]);
        
        %% Axis settings
        ui.panSelStat = uipanel(ui.f,'Position',[0.81 0.07 0.18 0.925]);
                
        %X-Axis
        ui.panelAxLim = uipanel(ui.panSelStat,...
            'Position',[0.02 0.76 0.96 0.23],...
            'Visible','on',...
            'BorderType','none'); %line/none
        
        ui.txtHistLimX = uicontrol('Parent',ui.panelAxLim,...
            'Style','text',...
            'Units','normalized',...
            'HorizontalAlignment','Left',...
            'FontSize',10,...
            'Position', [0.0 0.8 0.5 0.2],...
            'String', 'x-Axis Limits');
        
        ui.txtLimX = uicontrol('Parent',ui.panelAxLim,...
            'Style','text',...
            'Units','normalized',...
            'Position', [0.21 0.62 0.075 0.2],...
            'String', '-');
        ui.editLimX1  = uicontrol('Parent',ui.panelAxLim,...
            'Style','edit',...
            'String','0',...
            'Units','normalized',...
            'FontSize',9,...
            'Position',[0.0 0.68 0.22 0.16],...
            'Tag','x',...
            'Callback',@EditLimitsCB);
        ui.editLimX2  = uicontrol('Parent',ui.panelAxLim,...
            'Style','edit',...
            'String','1',...
            'FontSize',9,...
            'Units','normalized',...
            'Position',[0.3 0.68 0.22 0.16],...
            'Tag','x',...
            'Callback',@EditLimitsCB);
        
        ui.cboxLogX = uicontrol('Parent',ui.panelAxLim,...
            'Style','checkbox',...
            'Units','normalized',...
            'Position',[0.6 0.77 0.45 0.2],...
            'String','Logarithmic',...
            'Tag','logX',...
            'Callback',@EditLimitsCB);
        
        ui.cboxAutoX = uicontrol('Parent',ui.panelAxLim,...
            'Style','checkbox',...
            'Units','normalized',...
            'Value',1,...
            'Position',[0.6 0.62 0.45 0.2],...
            'String','Auto adjust',...            
            'Tag','autoX',...
            'Callback',@EditLimitsCB);
        
        %Y-Axis
       ui.txtHistLimY = uicontrol('Parent',ui.panelAxLim,...
            'Style','text',...
            'Units','normalized',...
            'FontSize',10,...  
            'HorizontalAlignment','Left',...
            'Position', [0.0 0.37 0.5 0.2],...
            'String', 'y-Axis Limits');
        
        ui.txtLimY = uicontrol('Parent',ui.panelAxLim,...
            'Style','text',...
            'Units','normalized',...
            'Position', [0.21 0.19 0.075 0.2],...
            'String', '-');
        ui.editLimY1  = uicontrol('Parent',ui.panelAxLim,...
            'Style','edit',...
            'String','0',...
            'Units','normalized',...
            'FontSize',9,...
            'Position',[0.0 0.25 0.22 0.16],...
            'Tag','y',...
            'Callback',@EditLimitsCB);
        ui.editLimY2  = uicontrol('Parent',ui.panelAxLim,...
            'Style','edit',...
            'String','1',...
            'FontSize',9,...
            'Units','normalized',...
            'Position',[0.3 0.25 0.22 0.16],...
            'Tag','y',...
            'Callback',@EditLimitsCB);
        
        ui.cboxLogY = uicontrol('Parent',ui.panelAxLim,...
            'Style','checkbox',...
            'Units','normalized',...
            'Position',[0.6 0.35 0.45 0.2],...
            'String','Logarithmic',...
            'Tag','logY',...
            'Callback',@EditLimitsCB);
        
        ui.cboxAutoY = uicontrol('Parent',ui.panelAxLim,...
            'Style','checkbox',...
            'Units','normalized',...
            'Value',1,...
            'Position',[0.6 0.20 0.45 0.2],...
            'String','Auto adjust',...            
            'Tag','AutoY',...
            'Callback',@EditLimitsCB);
        
        ui.cboxShowLegend = uicontrol('Parent',ui.panelAxLim,...
            'Style','checkbox',...
            'Units','normalized',...
            'Position',[0.00 0.0 0.45 0.2],...
            'String','Show legend',...
            'Value',1,...
            'Callback',@(~,~)PlotHistogramCB);
        
        % #bins
        ui.txtBinNum  = uicontrol('Parent',ui.panSelStat,...
            'Style','text',...
            'Units','normalized',...            
            'HorizontalAlignment','Left',...
            'FontSize',10,...
            'Position',[0.02 0.71 0.4 0.04],...
            'String','#bins');
        
        ui.editBinNum = uicontrol('Parent',ui.panSelStat,...
            'Style', 'edit',...
            'FontSize',9.5,...
            'Units','normalized',...
            'Visible','on',...
            'Position', [0.24 0.72 0.2 0.035],...
            'String','50',...
            'Callback',@(~,~)PlotHistogramCB);
        
        ui.txtLut  = uicontrol('Parent',ui.panSelStat,...
            'Style','text',...
            'Units','normalized',...            
            'HorizontalAlignment','Left',...
            'FontSize',10,...
            'Position',[0.5 0.71 0.4 0.04],...
            'String','LUT');
        
        ui.popLut = uicontrol('Parent',ui.panSelStat,...
            'Style', 'popupmenu',...
            'FontSize',9.5,...
            'Units','normalized',...
            'Visible','on',...
            'Position', [0.7 0.72 0.27 0.035],...
            'String',{'standard','winter','parula','jet','copper','gray'},...
            'Callback',@(~,~)PlotHistogramCB);
        
        %% Units panel
        
        ui.btnGroupUnits = uibuttongroup(ui.panSelStat,...
            'Units','normalized',...
            'Position', [.02  .59 .96 .12],...
            'Title','Units',...
            'SelectionChangedFcn',@btnGroupUnitsCB);
                
        ui.btnPxFr = uicontrol(ui.btnGroupUnits,...
            'Units','normalized',...
            'Position', [.05  .55  .5 .4],...
            'Style','radiobutton',...
            'String','pixels & frames',...
            'HorizontalAlignment','Left');
        
        ui.btnMiMs = uicontrol(ui.btnGroupUnits,...
            'Units','normalized',...
            'Position', [.55  .55  .6 .4],...
            'Style','radiobutton',...
            'String','microns & sec',...
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
            'Position', [.7  .1  .17 .3],...    
            'Visible','off',...
            'FontSize',9,...
            'Style','Edit',...
            'String',pixelSize,...
            'Callback',@btnGroupUnitsCB);
        
        %% Jumps to consider panel
        
        %Jumps to consider panel
        ui.panelJumpsToConsider = uipanel(ui.panSelStat,...
            'Position',[0.02 0.49 0.96 0.08],...
            'Visible','off',...
            'BorderType','none'); %line/none
        
        ui.txtJumpsToConsider = uicontrol('Parent',ui.panelJumpsToConsider,...
            'Style','text',...
            'Units','normalized',...
            'Position',[0.05 0.5 0.7 0.5],...
            'FontSize',10,...
            'HorizontalAlignment','Left',...
            'Visible','on',...
            'String','No. of jumps to consider');
        
        ui.editJumpsToConsider = uicontrol('Parent',ui.panelJumpsToConsider,...
            'Style','edit',...
            'Units','normalized',...
            'Position',[0.75 0.55 0.2 0.4],...
            'String',Inf,...
            'FontSize',9.5,...
            'Visible','on',...
            'Callback',@(~,~)CreateData);
        
        ui.cboxRemoveGaps = uicontrol('Parent',ui.panelJumpsToConsider,...
            'Style','checkbox',...
            'Units','normalized',...
            'Position',[0.05 0.0 0.95 0.5],...
            'FontSize',9,...
            'HorizontalAlignment','Left',...
            'Visible','on',...
            'Value',1,...
            'String','Remove jumps over gap frames',...
            'Callback',@(~,~)CreateData);
        
        %% Histogram normalization and style panel
        %Normalize by count or probability
        ui.panelNormalization = uipanel(ui.panSelStat,...
            'Position',[0.02 0.31 0.96 0.16],...
            'Visible','off',...
            'BorderType','none'); %line/none
                
        ui.btnGroupNormalization = uibuttongroup(ui.panelNormalization,...
            'Units','normalized',...
            'Position', [.0  .55   1 .45],...
            'Title','Normalization',...
            'SelectionChangedFcn',@(~,~)PlotHistogramCB);
                
        ui.btnCount = uicontrol(ui.btnGroupNormalization,...
            'Units','normalized',...
            'Position', [.05  .0   .5 .9],...
            'Style','radiobutton',...
            'String','Count',...
            'HorizontalAlignment','Left');
        
        ui.btnProbability = uicontrol(ui.btnGroupNormalization,...
            'Units','normalized',...
            'Position', [.55  .0   .6 .9],...
            'Style','radiobutton',...
            'String','Probability',...
            'HorizontalAlignment','Left');
        
        %Bars or Stairs
        ui.btnGroupBarsStairs = uibuttongroup(ui.panelNormalization,...
            'Units','normalized',...
            'Position', [.0  .0 1 .45],...
            'Title','Display style',...
            'Visible','off',...
            'SelectionChangedFcn',@(~,~)PlotHistogramCB);
                
        ui.btnBars = uicontrol(ui.btnGroupBarsStairs,...
            'Units','normalized',...
            'Position', [.05  .0  .5 .9],...
            'Style','radiobutton',...
            'String','Bars',...
            'HorizontalAlignment','Left');
        
        ui.btnStairs = uicontrol(ui.btnGroupBarsStairs,...
            'Units','normalized',...
            'Position', [.55  .0  .6 .9],...
            'Style','radiobutton',...
            'String','Stairs',...
            'HorizontalAlignment','Left');
         
        %% Panel show curves in cumulative plot
        ui.panelShowCurves = uipanel(ui.panSelStat,...
            'Position',[0.02 0.1 1 0.2],...
            'Visible','off',...
            'BorderType','none'); %line/none
        
        ui.cboxShowFit1 = uicontrol('Parent',ui.panelShowCurves,...
            'Style','checkbox',...
            'Units','normalized',...       
            'FontSize',9,...
            'Position',[0.05 0.9 0.8 0.1],...
            'String','Show 1-rate diffusion fit',...
            'Callback',@(~,~)PlotHistogramCB);
                
        ui.cboxShowFit2 = uicontrol('Parent',ui.panelShowCurves,...
            'Style','checkbox',...
            'Units','normalized',...   
            'FontSize',9,... 
            'Position',[0.05 0.7 0.8 0.1],...
            'String','Show 2-rate diffusion fit',...
            'Callback',@(~,~)PlotHistogramCB);
                 
        ui.cboxShowFit3 = uicontrol('Parent',ui.panelShowCurves,...
            'Style','checkbox',...
            'Units','normalized',...
            'FontSize',9,...
            'Position',[0.05 0.5 0.8 0.1],...
            'String','Show 3-rate diffusion fit',...
            'Callback',@(~,~)PlotHistogramCB);
        
        ui.cboxCumMovieWise = uicontrol('Parent',ui.panelShowCurves,...
            'Style','checkbox',...
            'Units','normalized',...
            'FontSize',9,...
            'Position',[0.05 0.3 0.8 0.1],...
            'String','Show movie-wise curves',...
            'Callback',@(~,~)PlotHistogramCB);
        
        %% Panel diffusion analysis
        ui.panelDiffParam = uipanel(ui.panSelStat,...
            'Position',[0.00 0.13 1 0.35],...
            'Visible','off',...
            'BorderType','none'); %line/none
        
        ui.btnGroupFitType = uibuttongroup(ui.panelDiffParam,...
            'Units','normalized',...
            'Position', [.02  .78   .96 .2],...
            'Title','Fit type',...
            'Visible','on',...
            'SelectionChangedFcn',@(~,~)PlotHistogramCB);
        
        ui.btnThreeRates = uicontrol(ui.btnGroupFitType,...
            'Units','normalized',...
            'Position', [.02  .0   .4 .9],...
            'Style','radiobutton',...
            'String','3 rates',...
            'HorizontalAlignment','Left');
        
        ui.btnTwoRates = uicontrol(ui.btnGroupFitType,...
            'Units','normalized',...
            'Position', [.37  .0   .4 .9],...
            'Style','radiobutton',...
            'String','2 rates',...
            'HorizontalAlignment','Left');
        
        ui.btnOneRate = uicontrol(ui.btnGroupFitType,...
            'Units','normalized',...
            'Position', [.7  .0   .4 .9],...
            'Style','radiobutton',...
            'String','1 rate',...
            'HorizontalAlignment','Left');
        
        ui.btnGroupFitVariable = uibuttongroup(ui.panelDiffParam,...
            'Units','normalized',...
            'Position', [.02  .56   .96 .2],...
            'Title','Parameter',...
            'Visible','on',...
            'SelectionChangedFcn',@(~,~)PlotHistogramCB);
        
        ui.btnShowD = uicontrol(ui.btnGroupFitVariable,...
            'Units','normalized',...
            'Position', [.02  .0   .3 .9],...
            'Style','radiobutton',...
            'String','D',...
            'HorizontalAlignment','Left');
        
        ui.btnShowA = uicontrol(ui.btnGroupFitVariable,...
            'Units','normalized',...
            'Position', [.3  .0   .3 .9],...
            'Style','radiobutton',...
            'String','A',...
            'HorizontalAlignment','Left');
        
        ui.btnShowEffectiveD = uicontrol(ui.btnGroupFitVariable,...
            'Units','normalized',...
            'Position', [.6  .0   .5 .9],...
            'Style','radiobutton',...
            'String','Effective D',...
            'HorizontalAlignment','Left');
        
        ui.textError = uicontrol('Parent',ui.panelDiffParam,...
            'Style','text',...
            'Units','normalized',...
            'Position',[0.02 0.43 0.8 0.1],...
            'HorizontalAlignment','Left',...
            'String','Displayed value and error',...
            'FontSize',9.5,...
            'Visible','on');
        
        ui.popError = uicontrol('Parent',ui.panelDiffParam,...
            'Style', 'popupmenu',...
            'FontSize',9.5,...
            'Units','normalized',...
            'Visible','on',...
            'Position', [0.02 0.34 .96 0.1],...
            'String',{'Pooled data, 95% confidence interval','Mean & stand. dev. of movie-wise values','Mean & stand. dev. of respampling values'},...
            'Callback',@(~,~)PlotHistogramCB);
        
        ui.panelResampling = uipanel(ui.panelDiffParam,...
            'Position',[0.00 0.19 1 0.12],...
            'Visible','off',...
            'BorderType','none'); %line/none
        
        ui.editNResampling = uicontrol('Parent',ui.panelResampling,...
            'Style','edit',...
            'Units','normalized',...
            'Position',[0.02 0.01 0.15 0.8],...
            'String',5,...
            'FontSize',9.5,...
            'Visible','on',...
            'Callback',@(~,~)PlotHistogramCB);
        
        ui.textNResampling1 = uicontrol('Parent',ui.panelResampling,...
            'Style','text',...
            'Units','normalized',...
            'Position',[0.2 0.0 0.4 0.75],...
            'String','resamplings with',...
            'HorizontalAlignment','Left',...
            'FontSize',9.5,...
            'Visible','on');
        
        ui.editPercResampling = uicontrol('Parent',ui.panelResampling,...
            'Style','edit',...
            'Units','normalized',...
            'Position',[0.6 0.01 0.15 0.8],...
            'String',50,...
            'FontSize',9.5,...
            'Visible','on',...
            'Callback',@(~,~)PlotHistogramCB);
               
        ui.textNResampling2 = uicontrol('Parent',ui.panelResampling,...
            'Style','text',...
            'Units','normalized',...
            'Position',[0.78 0.00 0.5 0.75],...
            'HorizontalAlignment','Left',...
            'String','% data',...
            'FontSize',9.5,...
            'Visible','on');
        
        ui.btnOverlayFitWithHist = uicontrol('Parent',ui.panelDiffParam,...
            'Units','normalized',...
            'Position',[0.02 0.11 0.96 0.12],...
            'String','Overlay fit with jump distance histogram',...
            'FontSize',9.5,...
            'Callback',@(~,~)OverlayFitWithHistCB);
        
        
        %Start values
        ui.tableStartD = uitable(ui.panSelStat,'Units','normalized',...
            'Position',[0.15 0.005 0.7 0.132],...
            'ColumnName',{'','Start value'},...         
            'Visible','off',...
            'RowName',{},...
            'Data',{'D1',0.1;'D2',1;'D3',10;},...
            'ColumnEditable',[false,true],...
            'CellEditCallback',@(~,~)PlotHistogramCB);
         
        %% Angles: Polarplot or lineplot
        
        ui.btnGroupPolarOrLine = uibuttongroup(ui.panSelStat,...
            'Units','normalized',...
            'Position', [.02  .205 .96 .08],...
            'Title','Histogram style',...
            'Visible','off',...
            'SelectionChangedFcn',@(~,~)PlotHistogramCB);
                
        ui.btnAnglesPolarplot = uicontrol(ui.btnGroupPolarOrLine,...
            'Units','normalized',...
            'Position', [.05  .0  .5 .9],...
            'Style','radiobutton',...
            'String','Polarplot',...
            'HorizontalAlignment','Left');
        
        ui.btnAnglesLineplot = uicontrol(ui.btnGroupPolarOrLine,...
            'Units','normalized',...
            'Position', [.55  .0  .6 .9],...
            'Style','radiobutton',...
            'String','Lineplot',...
            'HorizontalAlignment','Left');
        
        %% Angles: Jump distances making up the angle
        
        ui.btnGroupAnglesJumpDist = uipanel(ui.panSelStat,...
            'Units','normalized',...
            'Position', [.02  .07 .96 .12],...
            'Title','Jumps making up the angle',...
            'Visible','off');
                
        ui.btnMinJumpDist = uicontrol(ui.btnGroupAnglesJumpDist,...
            'Units','normalized',...
            'Position', [.05  .65  .5 .2],...
            'Style','text',...
            'String','Min. jump distance:',...
            'HorizontalAlignment','Left');
        
        ui.editAnglesMinJumpDist = uicontrol(ui.btnGroupAnglesJumpDist,...
            'Units','normalized',...
            'Position', [.5  .6  .18 .3],...    
            'Visible','on',...
            'FontSize',9,...
            'Style','Edit',...
            'String',0,...
            'Callback',@(~,~)CreateData);
        
        ui.txtAnglesMinJumpDist = uicontrol(ui.btnGroupAnglesJumpDist,...
            'Units','normalized',...
            'Position', [.7  .55  .2 .3],...    
            'Visible','on',...
            'FontSize',9,...
            'Style','text',...
            'String','px');
        
        ui.btnMaxJumpDist = uicontrol(ui.btnGroupAnglesJumpDist,...
            'Units','normalized',...
            'Position', [.05  .15  .6 .2],...
            'Style','text',...
            'String','Max. jump distance:',...
            'HorizontalAlignment','Left');
                                
        ui.editAnglesMaxJumpDist = uicontrol(ui.btnGroupAnglesJumpDist,...
            'Units','normalized',...
            'Position', [.5  .1  .18 .3],...    
            'Visible','on',...
            'FontSize',9,...
            'Style','Edit',...
            'String',Inf,...
            'Callback',@(~,~)CreateData);
        
        ui.txtAnglesMaxJumpDist = uicontrol(ui.btnGroupAnglesJumpDist,...
            'Units','normalized',...
            'Position', [.7  .05  .2 .3],...    
            'Visible','on',...
            'FontSize',9,...
            'Style','text',...
            'String','px');
        
       
        %% Illumination pattern panel
        ui.popParamSel.Value = 1;
        
        ui.panelITM = uipanel(ui.panSelStat,...
            'Position',[0.02 0.43 0.96 0.15],...
            'Visible','off',...
            'BorderType','none'); %'BorderType','line'
        
        ui.btnGroupITM = uibuttongroup(ui.panelITM,...
            'Units','normalized',...
            'Position', [.0  .48   1 .55],...
            'Title','Illumination pattern',...
            'SelectionChangedFcn',@bgITMselectionCB);
                
        ui.btnContinuous = uicontrol(ui.btnGroupITM,...
            'Units','normalized',...
            'Position', [.05  .0   .5 .9],...
            'Style','radiobutton',...
            'String','Continuous',...
            'HorizontalAlignment','Left');
        
        ui.btnITM = uicontrol(ui.btnGroupITM,...
            'Units','normalized',...
            'Position', [.47  .0   .6 .9],...
            'Style','radiobutton',...
            'String','Interlaced (ITM)',...
            'HorizontalAlignment','Left');
        
        ui.textNDarkForLong = uicontrol(ui.panelITM,...
            'Units','normalized',...
            'Position', [.01  .05   .65 .35],...
            'Style','text',...
            'String','Count as long track if track is longer than:',...
            'HorizontalAlignment','Left');
                
        ui.editNDarkForLong = uicontrol(ui.panelITM,...
            'Units','normalized',...
            'Position', [.8  .15  .17 .25],...         
            'FontSize',9,...
            'Style','Edit',...
            'String',3,...
            'Callback',@(~,~)bgITMselectionCB);
        
        ui.textNBrightFrames = uicontrol(ui.panelITM,...
            'Units','normalized',...
            'Position', [.01  .305   .65 .17],...
            'Visible','off',...
            'Style','text',...
            'String','#bright frames in one cycle',...
            'HorizontalAlignment','Left');
        
        ui.editNBrightFrames = uicontrol(ui.panelITM,...
            'Units','normalized',...
            'Position', [0.8 0.35 .17 .11],... 
            'Visible','off',...           
            'FontSize',9,...
            'Style','Edit',...
            'String',1,...
            'Callback',@(~,~)bgITMselectionCB);
               
        %% Group movies panel: normalize by ROI and batch or movie-wise
              
        ui.btnGroupSwarmVsMovie = uibuttongroup(ui.panSelStat,...
            'Units','normalized',...
            'Position', [.02  .34 .96 .08],...
            'Title','x-axis',...
            'Visible','off',...
            'SelectionChangedFcn',@(~,~)PlotHistogramCB);
                
        ui.btnValueVsBatchFile = uicontrol(ui.btnGroupSwarmVsMovie,...
            'Units','normalized',...
            'Position', [.05  .1  .5 .9],...
            'Style','radiobutton',...
            'String','Batch file',...
            'HorizontalAlignment','Left');
        
        ui.btnValueVsParameter = uicontrol(ui.btnGroupSwarmVsMovie,...
            'Units','normalized',...
            'Position', [.4  .1  .2 .9],...
            'Style','radiobutton',...
            'HorizontalAlignment','Left');
        
        ui.menuValueVsParameter = uicontrol(ui.btnGroupSwarmVsMovie,...
            'Units','normalized',...
            'Position', [.50  0  .45 .8],...
            'Style','popupmenu',...            
            'String',{'Movie number', 'No. of tracks', 'No. of non-linked spots', 'No. of all events',...
            'Mean jump distance', 'Avg. track length', 'Avg. no. of tracks per frame',...
            'Avg. no. of spots per frame', 'ROI size','No. of jumps'},...
            'HorizontalAlignment','Left',...
            'Callback',@(~,~)PlotHistogramCB);
                      
        %% MSD fit panel
        
        textHeight = .08;
        editWidth = 0.18;
        
        ui.panelMsdFit = uipanel(ui.panSelStat,...
            'Position',[0.02 0.02 0.95 0.45],...
            'Visible','off',...
            'BorderType','none'); %line/none
        
         ui.btnGroupFitFun = uibuttongroup(ui.panelMsdFit,...
            'Units','normalized',...
            'Position', [.0  .85 .99 .15],...
            'Title','Fit function',...
            'Tag','FitFun',...
            'SelectionChangedFcn',@PointsToFitBtnCB);
                
        ui.btnMSD = uicontrol(ui.btnGroupFitFun,...
            'Units','normalized',...
            'Position', [.05  .15  .5 .8],...
            'Style','radiobutton',...
            'String','Power law',...
            'Tag','MSD',...
            'HorizontalAlignment','Left');
        
        ui.btnLinear = uicontrol(ui.btnGroupFitFun,...
            'Units','normalized',...
            'Position', [.55  .15  .6 .8],...
            'Style','radiobutton',...
            'String','Linear',...            
            'Tag','Linear',...
            'HorizontalAlignment','Left');
        
        %Points to fit
        ui.txtPointsToFit  = uicontrol('Parent',ui.panelMsdFit,...
            'Style','text',...
            'Units','normalized',...
            'HorizontalAlignment','Left',...
            'FontSize',9,...
            'Position',[0.0 0.75 0.5 textHeight],...
            'String','#Points to fit');
                
        
        ui.editPointsToFit = uicontrol('Parent',ui.panelMsdFit,...
            'Style', 'edit',...
            'FontSize',9.5,...
            'Units','normalized',...
            'Visible','on',...
            'Position', [0.0 0.7 editWidth textHeight],...
            'String','90%');
    
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.35  .7   editWidth textHeight],...
            'String','5',...
            'Tag','PointsToFit',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);
        
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.55  .7   editWidth textHeight],...
            'String','60%',...
            'Tag','PointsToFit',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);
                
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.75  .7   editWidth textHeight],...
            'String','90%',...
            'Tag','PointsToFit',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);
        
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.35  .62   editWidth textHeight],...
            'String','10',...
            'Tag','PointsToFit',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);
        
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.55  .62   editWidth textHeight],...
            'String','75%',...
            'Tag','PointsToFit',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);
                
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.75  .62   editWidth textHeight],...
            'String','100%',...
            'Tag','PointsToFit',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);
        
        %Offset
        
        ui.txtOffset = uicontrol('Parent',ui.panelMsdFit,...
            'Style','text',...
            'Units','normalized',...
            'HorizontalAlignment','Left',...
            'FontSize',9,...
            'Position',[0.0 0.5 0.5 textHeight],...
            'String','Max. offset');
        
        ui.editOffset = uicontrol('Parent',ui.panelMsdFit,...
            'Style', 'edit',...
            'FontSize',9.5,...
            'Units','normalized',...
            'Visible','on',...
            'Position', [0.0 0.45 editWidth textHeight],...
            'String','0',...
            'Callback',@(~,~)PlotHistogramCB);
        
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.35  .45   editWidth textHeight],...
            'String','0',...
            'Tag','Offset',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);
        
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.55  .45   editWidth textHeight],...
            'String','0.5',...
            'Tag','Offset',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);
                
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.75  .45   editWidth textHeight],...
            'String','1',...
            'Tag','Offset',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);
        
       %Shortest track        
        ui.txtMsdShortestTrack  = uicontrol('Parent',ui.panelMsdFit,...
            'Style','text',...
            'Units','normalized',...            
            'HorizontalAlignment','Left',...
            'FontSize',9,...
            'Position',[0.0 0.3 0.5 textHeight],...
            'String','Shortest track');
        
        ui.editMsdShortestTrack = uicontrol('Parent',ui.panelMsdFit,...
            'Style', 'edit',...
            'FontSize',9.5,...
            'Units','normalized',...
            'Visible','on',...
            'Position', [0.5 0.3 editWidth textHeight],...
            'String','10',...
            'Callback',@(~,~)PlotHistogramCB);
        
        %Alpha Treshold
        ui.txtAlphaThres  = uicontrol('Parent',ui.panelMsdFit,...
            'Style','text',...
            'Units','normalized',...            
            'HorizontalAlignment','Left',...
            'FontSize',9,...
            'Position',[0.0 0.19 0.95 textHeight],...
            'String','Show tracks with alpha values below');
        
        ui.editAlphaThres = uicontrol('Parent',ui.panelMsdFit,...
            'Style', 'edit',...
            'FontSize',9.5,...
            'Units','normalized',...
            'Visible','on',...
            'Position', [0.0 0.13 editWidth textHeight],...
            'String','Inf',...
            'Callback',@(~,~)CreateData);
        
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.35  .13   editWidth textHeight],...
            'String','0.7',...
            'Tag','Alpha',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);
        
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.55  .13   editWidth textHeight],...
            'String','1',...
            'Tag','Alpha',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);
                
        uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'Position',[.75  .13   editWidth textHeight],...
            'String','Inf',...
            'Tag','Alpha',...
            'HorizontalAlignment','Left',...
            'Callback',@PointsToFitBtnCB);

                
        ui.btnFitMsd = uicontrol('Parent',ui.panelMsdFit,...
            'Units','normalized',...
            'FontSize',9,...
            'Visible','on',...
            'Position',[0.1 0.0 .8 0.1],...
            'String','Fit MSD',... <html><br>
            'Callback',@FitMsdCB);
              
        %% Normalize by ROI size
        
        ui.cboxNormalizeROI = uicontrol(ui.panSelStat,...
            'Style','checkbox',...
            'Units','normalized',...
            'Visible','off',...
            'Position',[0.05 0.28 0.8 0.04],...
            'String','Normalize by ROI-size',...
            'Callback',@(~,~)PlotHistogramCB);
        
        %% Distance from ROI
        
        ui.btnCalcDistance = uicontrol('Parent',ui.panSelStat,...
            'Units','normalized',...
            'FontSize',9,...
            'Visible','off',...
            'Position',[0.1 0.05 .8 0.05],...
            'String','Calculate distances',... <html><br>
            'Callback',@CalcDistanceFromRoiCB);
                        
        %% Export
        
        ui.btnCopyToWorkspace  = uicontrol('Parent',ui.f,...
            'Units','normalized',...
            'FontSize',8,...
            'Position',[0.81 0.01 .18 0.06],...
            'String','Export to Matlab workspace',... <html><br>
            'Callback',@CopyToWorkspaceCB);
        
        
        %% Other
        %Initialize datacursor
        ui.dataCoursor = datacursormode(ui.f);
        set(ui.dataCoursor,'UpdateFcn',{@dataCursorCB});
        
    end

%Data creation

    function CreateData()
        %Function to create all the results (except confinement radii and diffusion parameters)
        %Executed every time the composition of movies to show changes
        
        if isempty(batches)
            results = InitResults();
        else
            %Get amount of selected batch files
            nBatches = numel(ui.popBatchSel.Value);
            
            %Initialize results structure
            resultsInit = InitResults();
            results = repmat(resultsInit, 1, nBatches);
            
            %Initialize loop variables
            resultIdx = 1;
            nSelectedBatches = numel(ui.popBatchSel.Value);
            
            columnName = cell(nSelectedBatches,1);
            statData = cell(4,nSelectedBatches);

            %Iterate through selected batch files and create data for the
            %selected combination of batch file, frameCycleTime and
            %subRegion
            for datasetIdx = 1:nSelectedBatches
                curBatchNum = ui.popBatchSel.Value(datasetIdx);
                
                %Get batch file for this loop
                currentBatch = batches{curBatchNum};
                
                if ui.popTlSel.Value(1) == 1 
                    %User want to see results of one specific movie number
                    if length(batches{curBatchNum}) < curMovieIndex
                        moviesIdx = [];
                    else
                        moviesIdx = curMovieIndex;
                    end
                elseif ui.popTlSel.Value(1) == 2
                    %User want to see results of all movies
                    %Get movie indices of all movies in this batch
                    moviesIdx = 1:length(currentBatch);
                else
                    %User want to see movies with specific frameCycleTimes                 
                    %Find logical indices of movies matching these fcts
                    fctMovies = zeros(length(currentBatch),1);
                    for k = 1:numel(ui.popTlSel.Value)
                        curFct = str2double(ui.popTlSel.String{ui.popTlSel.Value(k)}(1:end-3));
                        fctMovies = or(fctMovies, frameCycleTimeMovieList{curBatchNum} == curFct);
                    end
                    
                    %Get indices of movies having the selected frameCycleTimes
                    moviesIdx = find(fctMovies');
                end
                
                resultsInCurBatch = [currentBatch(moviesIdx).results];
                batchName = ui.popBatchSel.String{curBatchNum}(1:end);

                if isempty(moviesIdx)
                    %There is no movie in this batch with the selected
                    %frame cycle time
                    results(resultIdx) = [];
                                        
                    nMovies = 0;
                    nTracks = 0;
                    nNonLinkedSpots = 0;
                    nAllEvents = 0;
                elseif all([resultsInCurBatch.nSpots]' == 0)
                    % none of the movies contains any spots
                    results(resultIdx) = [];
                                                   
                    nMovies = sum(moviesIdx);
                    nTracks = 0;
                    nNonLinkedSpots = 0;
                    nAllEvents = 0;
                else
                    %Create results for the current batch and the given movie indices
                    results(resultIdx) = create_histogram_data(currentBatch,moviesIdx, ui);
                    
                    %Save list of batch filenames in results structure
                    results(resultIdx).batchName = batchName;
                    
                    
                    nMovies = sum(results(resultIdx).roiSize(:) ~= 0);
                    nTracks = sum(vertcat(results(resultIdx).nTracks));
                    nNonLinkedSpots = sum(vertcat(results(resultIdx).nNonLinkedSpots));
                    nAllEvents = sum(vertcat(results(resultIdx).nAllEvents));
                    
                    
                    resultIdx = resultIdx + 1;
                    
                end
                
                
                %Create and display statistics at the lower left corner of the ui
                statData{1,datasetIdx} = nMovies;
                statData{2,datasetIdx} = nTracks;
                statData{3,datasetIdx} = nNonLinkedSpots;
                statData{4,datasetIdx} = nAllEvents;

                if curBatchNum<10
                    columnName{datasetIdx} = ['#', batchName(1)];
                else
                    columnName{datasetIdx} = ['#', batchName(1:2)];
                end
                
                
            end
            
            if ~isempty(columnName)
                ui.tableStatistics.Data = [ui.tableStatistics.Data(:,1), statData];
                
                ui.tableStatistics.ColumnName = [{''}; columnName];
            else
                ui.tableStatistics.ColumnName = {''};
                ui.tableStatistics.Data = {'#movies';'#tracks';'#non-linked spots';'#all events'};
            end
            
        end
        
        PlotHistogramCB()        
    end

    function FitMsdCB(src,~)
        %Executed when user presses "Fit MSD" button, fits the msd of all
        %tracks in all movies of all batches
        
        %Button turns red during execution
        src.BackgroundColor = 'r';
        
        %Get user settings for msd analysis
        shortestTrack = str2double(ui.editMsdShortestTrack.String);
        maxOffset = str2double(ui.editOffset.String);
        pointsToFit = ui.editPointsToFit.String;
        msdOrLinear = ui.btnGroupFitFun.SelectedObject.Tag;
        
        %Iterate through batches
        for batchIdx = 1:length(batches)
            %Get current batch
            curbatch = batches{batchIdx};
            %Get number of movies in current batch
            nMoviesInCurBatch = length(curbatch);
            %Iterate through movies of current batch
            
            for curMovieIdx = 1:nMoviesInCurBatch
                %Display fitting progress
                src.String = ['Fitting batch #',num2str(batchIdx), ', movie ', num2str(curMovieIdx), '/', num2str(nMoviesInCurBatch)];
                drawnow
                %Get tracks in current movie
                tracks = curbatch(curMovieIdx).results.tracks;
                %Fit msd and retrieve results
                msdResults = msd_analysis(tracks,shortestTrack,pointsToFit,maxOffset,msdOrLinear);
                
                curbatch(curMovieIdx).results.msdDiffConst = msdResults.msdDiffConst;
                curbatch(curMovieIdx).results.alphaValues = msdResults.alphaValues;
                curbatch(curMovieIdx).results.confRad = msdResults.confRad;
                curbatch(curMovieIdx).results.meanJumpDistConfRad = msdResults.meanJumpDist;
                
                nRegions = curbatch(curMovieIdx).results.nSubRegions+1;                
                trackSubRegionAssignment = curbatch(curMovieIdx).results.tracksSubRoi;
                
                
                for subRegionIdx = 1:nRegions
                    tracksInSubRegion  = subRegionIdx-1 == trackSubRegionAssignment;
                    curbatch(curMovieIdx).results.subRegionResults(subRegionIdx).msdDiffConst = msdResults.msdDiffConst(tracksInSubRegion);
                    curbatch(curMovieIdx).results.subRegionResults(subRegionIdx).alphaValues = msdResults.alphaValues(tracksInSubRegion);
                    curbatch(curMovieIdx).results.subRegionResults(subRegionIdx).confRad = msdResults.confRad(tracksInSubRegion);
                    curbatch(curMovieIdx).results.subRegionResults(subRegionIdx).meanJumpDistConfRad = msdResults.meanJumpDist(tracksInSubRegion);
                end
                
            end
            %Save fit results in batches structure
            batches{batchIdx} = curbatch;
        end
        
        %Reset button
        src.String = 'Fit MSD';
        src.BackgroundColor = [.94 .94 .94];
        
        %Update all results and plots
        CreateData()
    end

    function CalcDistanceFromRoiCB(src,~)

        oriString = src.String;
        src.BackgroundColor = 'r';
        drawnow
        for batchIdx =  1:length(batches)
            
            %Get current batch
            curbatch = batches{batchIdx};

            for movieIdx = 1 : length(curbatch)
                src.String = ['Batch #',num2str(batchIdx), ', Movie #', num2str(movieIdx)];
                drawnow
                
                tracksInSubRoiIdx = curbatch(movieIdx).results.tracksSubRoi;
                
                tracks = curbatch(movieIdx).results.tracks;
                
                nTracks = length(tracks);
                
%                 minDistCurMovie = cell(nTracks,1);
                minDistCurMovie = zeros(nTracks,1);
                
                for trackIdx = 1:nTracks
                    
                    %                     nDetectionsInTrack = size(tracks{trackIdx},1);
                    %
                    %                     for detectionIdx = 1:nDetectionsInTrack-1
                    %                         roiIdx = tracksInSubRoiIdx(trackIdx)+1;
                    %
                    %                         if roiIdx == 1
                    %                             subROI = curbatch(movieIdx).ROI{1};
                    %                         else
                    %                             subROI = vertcat(curbatch(movieIdx).subROI{roiIdx-1}{1}{:});
                    %                         end
                    %
                    %                         P = mean(tracks{trackIdx}(detectionIdx:detectionIdx+1,2:3));
                    %
                    %                         minDistCurMovie{trackIdx}(detectionIdx,1) = abs(p_poly_dist1(P(1), P(2), subROI(:,1), subROI(:,2)));
                    %
                    %                     end
                    
                    if isempty(tracksInSubRoiIdx)
                        roiIdx = 1;
                    else
                        
                        roiIdx = tracksInSubRoiIdx(trackIdx)+1;
                    end
                    
                    if roiIdx == 1
                        subROI = curbatch(movieIdx).ROI{1};
                    else
                        subROI = vertcat(curbatch(movieIdx).subROI{roiIdx-1}{1}{:});
                    end
                    
                    P = mean(tracks{trackIdx}(:,2:3));
                    minDistCurMovie(trackIdx) = abs(p_poly_dist1(P(1), P(2), subROI(:,1), subROI(:,2)));
                    
                end
                
                curbatch(movieIdx).results.distToRoiBorder = minDistCurMovie;
                
                nRegions = curbatch(movieIdx).results.nSubRegions+1;
                trackSubRegionAssignment = curbatch(movieIdx).results.tracksSubRoi;
                
                for subRegionIdx = 1:nRegions
                    tracksInSubRegion  = subRegionIdx-1 == trackSubRegionAssignment;
                    
                    curbatch(movieIdx).results.subRegionResults(subRegionIdx).distToRoiBorder = minDistCurMovie(tracksInSubRegion);
                    
%                     if any(tracksInSubRegion)
%                         curbatch(movieIdx).results.subRegionResults(subRegionIdx).distToRoiBorder = minDistCurMovie{tracksInSubRegion};
%                     else
%                         curbatch(movieIdx).results.subRegionResults(subRegionIdx).distToRoiBorder = {};
%                     end
                end
                
            end
            batches{batchIdx} = curbatch;
        end
        
        CreateData()
        src.String = oriString;
        src.BackgroundColor = [0.9400    0.9400    0.9400];
    end


%Update UI and data display

    function PlotHistogramCB()

        if isempty(batches) || isempty(results)
            plot(ui.ax,1)
            return
        end
        %Get number of selected batch files
        nBatches = length(results);
        
        selectedColor = ui.popLut.String{ui.popLut.Value};

        %Create colormap for plotting
        switch selectedColor
            case 'standard'
                batchColors = custom_colormap(nBatches);
            case 'jet'
                batchColors = jet(nBatches);
            case 'copper'             
                batchColors = copper(nBatches);
            case 'parula'
                batchColors = parula(nBatches);
            case 'winter'
                batchColors = winter(nBatches);
            case 'gray'
                batchColors = gray(nBatches+1);
                batchColors = batchColors(1:end-1,:);
        end
                        
        %Initialize struct array for exporting plot values to matlab
        %workspace
        currentPlotValues = repmat(struct,nBatches,1);
        
        %Hide all ui elements on the right side because we later enable the
        %ui elements corresponding to the current parameter selection
        ui.panelITM.Visible = 'off';
        ui.textNDarkForLong.Visible = 'off';
        ui.editNDarkForLong.Visible = 'off';
        ui.panelNormalization.Visible = 'off';
        ui.panelShowCurves.Visible = 'off';
        ui.cboxNormalizeROI.Visible = 'off';
        ui.panelJumpsToConsider.Visible = 'off';
        ui.tableStartD.Visible = 'off';
        ui.panelMsdFit.Visible = 'off';
        ui.panelDiffParam.Visible = 'off';
        ui.btnCalcDistance.Visible = 'off';
        ui.btnGroupAnglesJumpDist.Visible =  'off';
        ui.btnGroupSwarmVsMovie.Visible =  'off';
        ui.btnGroupPolarOrLine.Visible =  'off';
        ui.btnGroupBarsStairs.Visible =  'off';
        ui.panelResampling.Visible =  'off';
                
        if strcmp(ui.pax.Visible, 'on') 
            %Clean up polar
            %histogram
            delete(ui.polhist)
            ui.pax.Visible = 'off';
        end
        
        %Initialize x -and y labels
        xlabel1 = '';
        ylabel1 = '';
                        
        %Find currently selected analysis Parameter and adjust the ui and
        %prepare the corresponding values from the results structure
        switch ui.tabGroupParam.SelectedTab.Title            
            case 'Mobility'
                %% Mobility tab
                switch ui.popMobility.Value
                    case 1
                        % Jump Distance histogram
                        ui.panelJumpsToConsider.Visible = 'on';
                        if ui.btnPxFr.Value
                            xlabel1 = 'Jump distance (px)';
                        else
                            xlabel1 = 'Jump distance (\mum)';
                        end
                        
                        valuesY = cell(1,nBatches);
                        
                        for batchIdx = 1:nBatches   
                            valuesY{batchIdx} = vertcat(results(batchIdx).jumpDistances{:});
                        end
                                        
                        plotStyle = 'histogram';
                    case 2
                        % Cumulative jump Distance histogram
                        ui.panelJumpsToConsider.Visible = 'on';
                        
                        ui.panelShowCurves.Visible = 'on';
                        ui.tableStartD.Visible = 'on';
                        
                        
                        squareX = true;
                        
                        if ui.btnPxFr.Value
                            xlabel1 = 'd^2/(4\cdot\Deltat)  [px^2frame^{-1}]';
                        else
                            xlabel1 = 'd^2/(4\cdot\Deltat)  [\mum^{2}s^{-1}]';
                        end
                        
                        valuesY = cell(1,nBatches);
                        
                        for batchIdx = 1:nBatches                            
                            valuesY{batchIdx} = vertcat(results(batchIdx).jumpDistances{:});
                        end
                        
                        plotStyle = 'cumHistogram';
                    case 3
                        % Diffusion Parameter
                        
                        %Diffusion parameter is a special case because we
                        %need to fit the diffusion curves first. Here we
                        %adjust the visible ui elemts and do the rest
                        %later.
                        
                        ui.tableStartD.Visible = 'on';
                        ui.panelDiffParam.Visible = 'on';
                        ui.panelJumpsToConsider.Visible = 'on';
                        
                        if ui.popError.Value == 3
                            ui.panelResampling.Visible =  'on';
                        end
                        
                        plotStyle = 'diffusion parameters';
                    case 4 
                        % Mean jump distance histogram
                        ui.panelJumpsToConsider.Visible = 'on';
                        
                        if ui.btnPxFr.Value
                            xlabel1 = 'Mean jump distance (px)';
                        else
                            xlabel1 = 'Mean jump distance  (\mum)';
                        end
                        
                        valuesY = cell(1,nBatches);
                        
                        for batchIdx = 1:nBatches
                            valuesY{batchIdx} = vertcat(results(batchIdx).meanJumpDists{:});
                        end
                        
                        plotStyle = 'histogram';
                    case 5 
                        % Mean jump distance cumulative histogram
                        ui.panelJumpsToConsider.Visible = 'on';
                        
                        if ui.btnPxFr.Value
                            xlabel1 = 'Mean jump distance (px)';
                        else
                            xlabel1 = 'Mean jump distance  (\mum)';
                        end
                        
                        valuesY = cell(1,nBatches);
                        
                        squareX = false;
                        
                        for batchIdx = 1:nBatches                            
                            valuesY{batchIdx} = vertcat(results(batchIdx).meanJumpDists{:});
                        end
                        
                                                
                        plotStyle = 'cumHistogram';
                        ui.cboxShowFit1.Value = 0;
                        ui.cboxShowFit2.Value = 0;
                        ui.cboxShowFit3.Value = 0;
                    case 6
                        % Mean jump distance per movie
                        ui.panelJumpsToConsider.Visible = 'on';
                        ui.cboxNormalizeROI.Value = 0;

                        valuesY = {results(:).meanJumpDistMoviewise};
                        isTrackedFraction = 0;
                        
                        if ui.btnPxFr.Value
                            ylabel1 = 'Mean jump distance (pixel)';
                        else
                            ylabel1 = 'Mean jump distance (µm)';
                        end
                        
                        ui.btnGroupSwarmVsMovie.Visible =  'on';
                        
                        if ui.btnValueVsBatchFile.Value
                            plotStyle = 'swarmplot';
                        elseif ui.btnValueVsParameter.Value
                            plotStyle = 'valueVsParameter';
                        end
                        
                    case 7
                        % Angles
                        ui.btnGroupAnglesJumpDist.Visible =  'off';
                        ui.panelJumpsToConsider.Visible = 'on';
                        ui.btnGroupAnglesJumpDist.Visible =  'on';
                        ui.btnGroupPolarOrLine.Visible =  'on';
                        
                        valuesY = cell(1,nBatches);
                        
                        if ui.btnAnglesPolarplot.Value
                            plotStyle = 'angular histogram';
                            for batchIdx = 1:nBatches
                                valuesY{batchIdx} = vertcat(results(batchIdx).angles{:});
                            end
                        else
                            plotStyle = 'histogram';
                            for batchIdx = 1:nBatches
                                curBatchAngles = vertcat(results(batchIdx).angles{:});
                                curBatchAngles = curBatchAngles*180/pi;
                                curBatchAngles(curBatchAngles < 0) = curBatchAngles(curBatchAngles < 0)+360;
                                valuesY{batchIdx} = curBatchAngles;
                            end                            
                        end
                        
                        
                        delete(ui.hist)
                        cla(ui.ax)
                        ui.ax.Visible = 'off';
                   
                    case 8
                        % Confinement radii
                        ui.panelNormalization.Visible = 'on';
                        ui.panelMsdFit.Visible = 'on';
                        
                        if ui.btnPxFr.Value
                            xlabel1 = 'Confinement radius (px)';
                        else
                            xlabel1 = 'Confinement radius (\mum)';
                        end
                        
                        valuesY = cell(1,nBatches);
                        
                        for batchIdx = 1:nBatches                            
                            valuesY{batchIdx} = vertcat(results(batchIdx).confRad{:});
                        end
                        
                        plotStyle = 'histogram';
                    case 9
                        % Confinement radius vs. mean jump distance
                        plotStyle = 'scatterplot';
                        ui.panelMsdFit.Visible = 'on';
                        if ui.btnPxFr.Value
                            xlabel1 = 'Confinement radius (px)';
                            ylabel1 = 'Mean jump distance (px)';
                        else
                            xlabel1 = 'Confinement radius (\mum)';
                            ylabel1 = 'Mean jump distance (\mum)';
                        end
                        
                        valuesY = cell(1,nBatches);

                        for batchIdx = 1:nBatches   
                            valuesY{batchIdx} = vertcat(results(batchIdx).meanJumpDistConfRad{:});
                        end
                        
                        valuesY = vertcat(valuesY{:});
                        valuesX = cell(1,nBatches);

                        for batchIdx = 1:nBatches   
                            valuesX{batchIdx} = vertcat(results(batchIdx).confRad{:});
                        end
                        
                        valuesX = vertcat(valuesX{:});
                        
                    case 10
                        % Alpha values from MSD
                        ui.panelMsdFit.Visible = 'on';
                        
                        xlabel1 = 'Alpha value';
                        
                        valuesY = cell(1,nBatches);
                        
                        for batchIdx = 1:nBatches   
                            valuesY{batchIdx} = vertcat(results(batchIdx).alphaValues{:});
                        end
                                                
                        plotStyle = 'histogram';
                    case 11
                        % Diffusion constants from MSD
                        ui.panelMsdFit.Visible = 'on';
                        
                        if ui.btnPxFr.Value
                            xlabel1 = 'Diffusion coefficient (px^2/frame)';
                        else
                            xlabel1 = 'Diffusion coefficient (\mum^2/sec)';
                        end
                        
                        valuesY = cell(1,nBatches);
                        
                        for batchIdx = 1:nBatches   
                            valuesY{batchIdx} = vertcat(results(batchIdx).msdDiffConst{:});
                        end
                                                                       
                        plotStyle = 'histogram';
                end
            case 'Tracked fractions'
                %% Tracked fractions tab
                ui.panelITM.Visible = 'on';
                
                ui.cboxNormalizeROI.Visible = 'on';
                ui.btnGroupSwarmVsMovie.Visible =  'on';
                if ui.btnValueVsBatchFile.Value
                    plotStyle = 'swarmplot';
                elseif ui.btnValueVsParameter.Value
                    plotStyle = 'valueVsParameter';
                end
                
                switch ui.popTrackedFraction.Value
                    case 1 %Tracks vs. all events
                        trackedFractions = [results(:).trackedFractions];
                        valuesY = {trackedFractions.allTracksVsAllEventsMoviewise};
                        wholeSet = cell2mat({trackedFractions.allTracksVsAllEventsPooled});
                        wholeSetErr = cell2mat({trackedFractions.errorAllTracksVsAllEventsPooled});
                        valuesMean = cell2mat({trackedFractions.allTracksVsAllEventsMean});
                        valuesStd = cell2mat({trackedFractions.allTracksVsAllEventsStd});
                        isTrackedFraction = 1;
                        if ui.cboxNormalizeROI.Value
                            if ui.btnPxFr.Value   
                                ylabel1 = 'Tracks vs. all events (%/pixel)';
                            else
                                ylabel1 = 'Tracks vs. all events (%/µm^2)';
                            end
                        else
                            ylabel1 = 'Ttracks vs. all events';
                        end
                    case 2 %Long tracks vs. all events
                        ui.textNDarkForLong.Visible = 'on';
                        ui.editNDarkForLong.Visible = 'on';
                        trackedFractions = [results(:).trackedFractions];
                        valuesY = {trackedFractions.longVsAllEventsMoviewise};
                        wholeSet = cell2mat({trackedFractions.longVsAllEventsPooled});
                        wholeSetErr = cell2mat({trackedFractions.errorLongVsAllEventsPooled});
                        valuesMean = cell2mat({trackedFractions.longVsAllEventsMean});
                        valuesStd = cell2mat({trackedFractions.longVsAllEventsStd});
                        isTrackedFraction = 1;
                        
                        if ui.cboxNormalizeROI.Value
                            if ui.btnPxFr.Value
                                ylabel1 = 'Long tracks vs. all events (%/pixel)';
                            else
                                ylabel1 = 'Long tracks vs. all events (%/µm^2)';
                            end
                        else
                            ylabel1 = 'Long tracks vs. all events';
                        end
                    case 3 %Short tracks vs. all events
                        ui.textNDarkForLong.Visible = 'on';
                        ui.editNDarkForLong.Visible = 'on';
                        trackedFractions = [results(:).trackedFractions];
                        valuesY = {trackedFractions.shortVsAllEventsMoviewise};
                        wholeSet = cell2mat({trackedFractions.shortVsAllEventsPooled});
                        wholeSetErr = cell2mat({trackedFractions.errorShortVsAllEventsPooled});
                        valuesMean = cell2mat({trackedFractions.shortVsAllEventsMean});
                        valuesStd = cell2mat({trackedFractions.shortVsAllEventsStd});
                        isTrackedFraction = 1;
                        
                        if ui.cboxNormalizeROI.Value
                            if ui.btnPxFr.Value
                                ylabel1 = 'Short tracks vs. all events (%/pixel)';
                            else
                                ylabel1 = 'Short tracks vs. all events (%/µm^2)';
                            end
                        else
                            ylabel1 = 'Short tracks vs. all events';
                        end
                        
                        
                    case 4 %Long tracks vs. short tracks
                        ui.textNDarkForLong.Visible = 'on';
                        ui.editNDarkForLong.Visible = 'on';
                        
                        trackedFractions = [results(:).trackedFractions];
                        valuesY = {trackedFractions.longVsAllTracksMoviewise};
                        wholeSet = cell2mat({trackedFractions.longVsAllTracksPooled});
                        wholeSetErr = cell2mat({trackedFractions.errorLongVsAllTracksPooled});
                        valuesMean = cell2mat({trackedFractions.longVsAllTracksMean});
                        valuesStd = cell2mat({trackedFractions.longVsAllTracksStd});
                        isTrackedFraction = 1;
                        if ui.cboxNormalizeROI.Value
                            if ui.btnPxFr.Value   
                                ylabel1 = 'Long tracks vs. all tracks (%/pixel)';
                            else
                                ylabel1 = 'Long tracks vs. all tracks (%/µm^2)';
                            end
                        else
                            ylabel1 = 'Long tracks vs. all tracks';
                        end
                    case 5 %No. of tracks                        
                        valuesY = {results(:).nTracks};
                        isTrackedFraction = 0;
                        if ui.cboxNormalizeROI.Value
                            if ui.btnPxFr.Value   
                                ylabel1 = 'No. of tracks/pixel';
                            else
                                ylabel1 = 'No. of tracks/µm^2';
                            end
                        else
                            ylabel1 = 'No. of tracks';
                        end
                    case 6 %No. of non-linked spots                        
                        valuesY = {results(:).nNonLinkedSpots};
                        isTrackedFraction = 0;
                        if ui.cboxNormalizeROI.Value
                            if ui.btnPxFr.Value   
                                ylabel1 = 'No. of non-linked spots/pixel';
                            else
                                ylabel1 = 'No. of non-linked spots/µm^2';
                            end
                        else
                            ylabel1 = 'No. of non-linked spots';
                        end
                    case 7 %No. of all events                        
                        valuesY = {results(:).nAllEvents};
                        isTrackedFraction = 0;
                        if ui.cboxNormalizeROI.Value
                            if ui.btnPxFr.Value   
                                ylabel1 = 'No. of all events/pixel';
                            else
                                ylabel1 = 'No. of all events/µm^2';
                            end
                        else
                            ylabel1 = 'No. of all events';
                        end
                    case 8 %No. of long tracks
                        ui.textNDarkForLong.Visible = 'on';
                        ui.editNDarkForLong.Visible = 'on';
                       
                        valuesY = {results(:).nLong};
                        isTrackedFraction = 0;
                        
                        if ui.cboxNormalizeROI.Value
                            if ui.btnPxFr.Value   
                                 ylabel1 = 'No. of long tracks/pixel';
                            else
                                 ylabel1 = 'No. of long tracks/µm^2';
                            end
                        else
                             ylabel1 = 'No. of long tracks';
                        end
                    case 9 %No. of short tracks
                        ui.textNDarkForLong.Visible = 'on';
                        ui.editNDarkForLong.Visible = 'on';
                        
                        if ui.cboxNormalizeROI.Value
                            if ui.btnPxFr.Value    
                                ylabel1 = 'No. of short tracks/pixel';
                            else
                                ylabel1 = 'No. of short tracks/µm^2';
                            end
                        else
                            ylabel1 = 'No. of short tracks';
                        end

                        valuesY = {results(:).nShort};
                        isTrackedFraction = 0;
                end
            case 'Statistics'
                %% Statistics tab
                
                isTrackedFraction = 0;
                switch ui.popStatistics.Value
                    case 1% Track lengths
                        
                        plotStyle = 'histogram';
                        if ui.btnPxFr.Value
                            xlabel1 = 'Tracklength (frames)';
                        else
                            xlabel1 = 'Tracklength (sec)';
                        end
                        
                        
                        valuesY = cell(1,nBatches);
                        
                        for batchIdx = 1:nBatches
                            valuesY{batchIdx} = vertcat(results(batchIdx).trackLengths{:});
                        end
                        
                    case 2%Cumulative Track lengths
                        
                        plotStyle = 'cumHistogram';
                        
                        if ui.btnPxFr.Value
                            xlabel1 = 'Tracklength (frames)';
                        else
                            xlabel1 = 'Tracklength (sec)';
                        end
                        
                        squareX = false;
                        valuesY = cell(1,nBatches);
                        
                        for batchIdx = 1:nBatches
                            valuesY{batchIdx} = vertcat(results(batchIdx).trackLengths{:});
                        end
                        
                    case 3 %Avg. track lengths
                        
                        ui.cboxNormalizeROI.Visible = 'on';
                        
                        ui.btnGroupSwarmVsMovie.Visible =  'on';
                        if ui.btnValueVsBatchFile.Value
                            plotStyle = 'swarmplot';
                        elseif ui.btnValueVsParameter.Value
                            plotStyle = 'valueVsParameter';
                        end
                        
                        valuesY = {results(:).meanTrackLength};
                        
                        
                        if ui.btnPxFr.Value                            
                            if ui.cboxNormalizeROI.Value
                                ylabel1 = 'Average tracklength (frames/pixel)';
                            else
                                ylabel1 = 'Average tracklength (frames)';
                            end
                        else
                            
                            if ui.cboxNormalizeROI.Value
                                ylabel1 = 'Average tracklength (sec/µm^2)';
                            else
                                ylabel1 = 'Average tracklength (sec)';
                            end
                        end
                    case 4 %Avg. #tracks per frame
                        
                        ui.btnGroupSwarmVsMovie.Visible =  'on';
                        if ui.btnValueVsBatchFile.Value
                            plotStyle = 'swarmplot';
                        elseif ui.btnValueVsParameter.Value
                            plotStyle = 'valueVsParameter';
                        end
                        
                        ui.cboxNormalizeROI.Visible = 'on';
                        valuesY = {results(:).meanTracksPerFrame};
                                                
                        if ui.cboxNormalizeROI.Value
                            if ui.btnPxFr.Value    
                                ylabel1 = 'Average no. of tracks / (frame \cdot pixel)';
                            else
                                ylabel1 = 'Average no. of tracks / (frame \cdot µm^2)';
                            end
                        else
                            ylabel1 = 'Average no. of tracks / frame';
                        end
                    case 5 %Avg. #spots per frame
                        
                        ui.cboxNormalizeROI.Visible = 'on';
                        ui.btnGroupSwarmVsMovie.Visible =  'on';
                        if ui.btnValueVsBatchFile.Value
                            plotStyle = 'swarmplot';
                        elseif ui.btnValueVsParameter.Value
                            plotStyle = 'valueVsParameter';
                        end
                        
                        valuesY = {results(:).meanSpotsPerFrame};
                        if ui.cboxNormalizeROI.Value
                            if ui.btnPxFr.Value
                                ylabel1 = 'Average no. of spots / (frame \cdot pixel)';
                            else
                                ylabel1 = 'Average no. of spots / (frame \cdot µm^2)';
                            end
                        else
                            ylabel1 = 'Average no. of spots / frame';
                        end
                    case 6 %No. of jumps                     
                        ui.cboxNormalizeROI.Visible = 'on';
                        ui.btnGroupSwarmVsMovie.Visible =  'on';
                        ui.panelJumpsToConsider.Visible = 'on';
                        
                        if ui.btnValueVsBatchFile.Value
                            plotStyle = 'swarmplot';
                        elseif ui.btnValueVsParameter.Value
                            plotStyle = 'valueVsParameter';
                        end
                        
                            ylabel1 = 'No. of jumps';
                        
                        valuesY = {results(:).nJumps};
                    case 7 %ROI Size
                        
                        ui.cboxNormalizeROI.Visible = 'on';                                                
                        ui.btnGroupSwarmVsMovie.Visible =  'on';
                        if ui.btnValueVsBatchFile.Value
                            plotStyle = 'swarmplot';
                        elseif ui.btnValueVsParameter.Value
                            plotStyle = 'valueVsParameter';
                        end
                        
                        if ui.btnPxFr.Value
                            ylabel1 = 'No. of pixels';
                        else
                            ylabel1 = 'ROI size (\mum^2)';
                        end
                        
                        valuesY = {results(:).roiSize};
                    case 8 %Distance from ROI border
                        plotStyle = 'histogram';
                        
                        ui.btnCalcDistance.Visible = 'on';
                        if ui.btnPxFr.Value
                            xlabel1 = 'Distance from ROI border (px)';
                        else
                            xlabel1 = 'Distance from ROI border (\mum)';
                        end
                                                
                        valuesY = cell(1,nBatches);
                        
                        for batchIdx = 1:nBatches   
                            valuesY{batchIdx} = vertcat(results(batchIdx).distToRoiBorder{:});
                            if isempty(valuesY{batchIdx})
                                valuesY{batchIdx} = 0;
                            end
                        end

                        
                    case 9 %Mean jump dist vs. distance from border
                        plotStyle = 'scatterplot';
                        ui.btnCalcDistance.Visible = 'on';
                        if ui.btnPxFr.Value
                            ylabel1 = 'Distance from ROI border (px)';
                            xlabel1 = 'Mean jump distance (px)';
                        else
                            ylabel1 = 'Distance from ROI border (\mum)';
                            xlabel1 = 'Mean jump distance (\mum)';
                        end

                        
                        
                        valuesY = cell(1,nBatches);
                        valuesX = cell(1,nBatches);
%                         
                        for batchIdx = 1:nBatches   
                            valuesY{batchIdx} = vertcat(results(batchIdx).distToRoiBorder{:});
                            valuesX{batchIdx} = vertcat(results(batchIdx).meanJumpDists{:});
%                             valuesY{batchIdx} = vertcat(results(batchIdx).jumpDistances{:});
                        end
%                                                 
                        valuesX = vertcat(valuesX{:});
                        valuesY = vertcat(valuesY{:});

                        
%                         valuesX = valuesX(valuesX < 30);
%                         valuesY = valuesY(valuesX < 30);      
                        
                        
                end
        end
        
        switch plotStyle
            case 'histogram'
                %% Plot as histogram
                                                
                %Make count/probability panel visible
                ui.panelNormalization.Visible = 'on';
                ui.btnGroupBarsStairs.Visible =  'on';
                
                %Get maximum value in the dataset and create histogram
                %limits
                histLimits = [min(0, min(cell2mat(valuesY'))) max(cell2mat(valuesY'))];
                
                %Adjust y-label of the axis to the selected normalization
                if ui.btnCount.Value                    
                    histNorm = 'Count';
                else                    
                    histNorm = 'Probability';
                end
                                
                ylabel1 = histNorm;
                
                %Save batch file names of selected batches into
                %currentPlotValues structure
                [currentPlotValues(:).batchName] = ui.popBatchSel.String{ui.popBatchSel.Value};
                
                %Initialize legend entries variable
                legendString = cell(numel(nBatches),1);
                
                %Iterate through selected batches
                for resultsIdx = 1:nBatches
                    
                    %Create legend entry
                    if strcmp(results(resultsIdx).batchName(2), ':')
                        legendString{resultsIdx} = ['Batch #', results(resultsIdx).batchName(1)];
                    else
                        legendString{resultsIdx} = ['Batch #', results(resultsIdx).batchName(1:2)];
                    end
                    
                    if resultsIdx == 2
                        hold(ui.ax,'on')
                    end
                    
                    if ui.btnBars.Value
                        edgeColor1 = 'k';
                        edgeAlpha1 = .25;
                        displayStyle1 = 'bar';
                        faceColor1 = batchColors(resultsIdx,:);
                    elseif ui.btnStairs.Value
                        edgeColor1 = batchColors(resultsIdx,:);
                        faceColor1 = 'none';
                        edgeAlpha1 = 1;
                        displayStyle1 = 'stairs';
                    end
                    
                    %Create and plot histogram for current batch
                    ui.hist = histogram(valuesY{resultsIdx},...
                        'EdgeColor',edgeColor1,...
                        'EdgeAlpha',edgeAlpha1,...
                        'DisplayStyle',displayStyle1,...
                        'FaceColor',faceColor1,...
                        'BinLimits',histLimits,...
                        'Normalization',histNorm,...
                        'NumBins',str2double(ui.editBinNum.String),...
                        'Parent',ui.ax);
                    
                    %Save results ins currentPlotValues structure
                    currentPlotValues(resultsIdx).histogramData = [(ui.hist.BinEdges(2:end)-(ui.hist.BinEdges(2)-ui.hist.BinEdges(1))/2)', ui.hist.Values'];
                    
                end
            case 'cumHistogram'
                %% Plot as cumulative histogram
                
                %Create y-axis label entry
                ylabel1 = 'Probability';
                
                %Save batch file names of selected batches into
                %currentPlotValues structure
                [currentPlotValues(:).batchName] = ui.popBatchSel.String{ui.popBatchSel.Value};
                
                %Initialize legend entries variable
                legendString = cell(nBatches,1);
                
                %Initialize variable that is set to true if a batch file contains different
                %tracking radii
                showWarnDlg = false;
                
                % Iterate through selected batches
                for resultsIdx = 1:nBatches
                    
                    %Create legend entry
                    if strcmp(results(resultsIdx).batchName(2), ':')
                        batchNum = results(resultsIdx).batchName(1);
                    else
                        batchNum = results(resultsIdx).batchName(1:2);
                    end
                    
                    legendString{resultsIdx} = ['Batch #', batchNum];
                    
                    if resultsIdx == 2
                        hold(ui.ax,'on')
                    end
                    
                    %Create cumulative density function from valuesY for
                    %current batch
                    [y,x] = histcounts(valuesY{resultsIdx},...
                        'Normalization','cdf',...
                        'NumBins',str2double(ui.editBinNum.String));
                    
                    %Get bin centers
                    x = (x(2:end)-(x(2)-x(1))/2);                    

                    %Square x axis for diffusion fit
                    if squareX
                        
                        %Get frame cycle times as we need them to calculate the
                        %diffusion constants
                        frameCycleTimes = {results(:).frameCycleTimes};
                        
                        if ui.btnPxFr.Value
                            %User wants results in pixels and frames
                            frameCycleTime = 1;
                        else
                            %User want results in microns and seconds
                            
                            %Get frame cycle times in movies of current batch
                            curframeCycleTimes = unique(frameCycleTimes{resultsIdx});
                            %Convert frame cycle time to seconds. If more than
                            %one frame cycle time is found, use first one
                            frameCycleTime = curframeCycleTimes(1);
                        end
                        
                        x = x.^2;
                        x = x./(4*frameCycleTime);
                    end
                    
                    
                    
                    %Save in currentPlotValues structure
                    currentPlotValues(resultsIdx).histogramData = [x', y'];
                    
                    %Plot histogram
                    ui.hist = plot(ui.ax,x,y,'-','Color', batchColors(resultsIdx,:));
                    
                    if (ui.cboxShowFit1.Value || ui.cboxShowFit2.Value || ui.cboxShowFit3.Value)
                        
                        %Retreive start values for fitting
                        startD = vertcat(ui.tableStartD.Data{:,2});
                        
                        %Get tracking radius
                        allTrackingRadii = vertcat(results(resultsIdx).trackingRadii);
                        trackingRadius = max(unique(allTrackingRadii));
                        trackingRadius = trackingRadius^2/(4*frameCycleTime);
                        nTrackingRadii = numel(unique(allTrackingRadii));
                        
                        % Show fits
                        if ui.cboxShowFit1.Value
                            %Show 1-exp fit
                            
                            %Fit curve with 1-exp diffusion fit
                            outOne = dispfit_cumulative(x',y', trackingRadius, startD, 1);
                            hold(ui.ax,'on');
                            %Save results in currentPlotValues structure
                            currentPlotValues(resultsIdx).singleExpFitCurve = outOne.xy;
                            %Plot fit curve
                            plot(ui.ax, outOne.xy(:,1), outOne.xy(:,2),'g','linewidth',2)
                            
                        end
                        
                        if ui.cboxShowFit2.Value
                            %Show 2-exp fit
                            
                            %Fit curve with 2-exp diffusion fit
                            outTwo = dispfit_cumulative(x',y', trackingRadius, startD, 2);
                            hold(ui.ax,'on');
                            %Save results in currentPlotValues structure
                            currentPlotValues(resultsIdx).doubleExpFitCurve = outTwo.xy;
                            %Plot fit curve
                            plot(ui.ax, outTwo.xy(:,1), outTwo.xy(:,2),'m','linewidth',2)
                        end
                        
                        if ui.cboxShowFit3.Value
                            %Show 3-exp fit
                            
                            %Fit curve with 3-exp diffusion fit
                            outThree = dispfit_cumulative(x',y', trackingRadius, startD, 3);
                            hold(ui.ax,'on');
                            %Save results in currentPlotValues structure
                            currentPlotValues(resultsIdx).threeExpFitCurve = outThree.xy;
                            %Plot fit curve
                            plot(ui.ax, outThree.xy(:,1), outThree.xy(:,2),'r','linewidth',2)
                          
                        end
                        
                        if nTrackingRadii > 1
                            showWarnDlg = true;
                        end
                        
                        
                        if showWarnDlg
                            warndlg(['At least one batch file contains '...
                                'movies that have been analyzed with different '...
                                'tracking radii. The tracking radius is used to '...
                                'normalize the fitting function (see manual). '...
                                'Using largest tracking radius in batch for normalization.'],'Warning');
                        end
                    end
                    
                end
                
                % Show movie-wise curves
                if ui.cboxCumMovieWise.Value
                    
                    
                    %Get current amount of already plotted curves to later bring
                    %them to the front
                    nLines = numel(get(gca,'Children'));
                    hold(ui.ax,'on');
                    %Iterate through selected batches
                    for resultsIdx = 1:nBatches
                        
                        %Initialize cell array to save movie-wise curves
                        curMoviewiseCurves = cell(length(results(resultsIdx).movieNames),1);
                        
                        %Iterate through movies of current batch
                        for movieIdx = 1:length(results(resultsIdx).movieNames)
                            
                            %Retreive jump distances of all movies of
                            %current batch
                            
                            dispsCurMovie = results(resultsIdx).jumpDistances{movieIdx};
                            
                            %Create cumulative density function
                            [y,x] = histcounts(dispsCurMovie,'BinMethod','integers',...
                                'Normalization','cdf',...
                                'NumBins',str2double(ui.editBinNum.String));
                            
                            if squareX
                                %Square the x-axis
                                x = (x(2:end)-(x(2)-x(1))/2).^2;
                                x = x./(4*frameCycleTime);
                            end
                            
                            %Plot the curves
                            plot(ui.ax, x, y,'.-','linewidth',2, 'Color', [.8 .8 .8])
                            
                            %Save current curve
                            curMoviewiseCurves{movieIdx} = [x',y'];
                        end
                        %Save movie-wise curves into currentPlotValues structure
                        currentPlotValues(resultsIdx).movieWiseCurves = curMoviewiseCurves;
                    end
                    
                    
                    %Send the movie-wise curves to the background
                    h2 = get(gca,'Children');
                    set(gca,'Children', [h2(end-nLines:end); h2(1:end-nLines)])
                    
                end
                
            case 'angular histogram'
                %% Plot as angular plot
                
                                
                %Make count/probability panel visible
                ui.panelNormalization.Visible = 'on';
                ui.btnGroupBarsStairs.Visible =  'on';
                
                %Set histogram normalization according to selected value
                if ui.btnCount.Value
                    histNorm = 'Count';
                else
                    histNorm = 'Probability';
                end
                
                %Set empty legend String for cartesian axis
                legendString = {};
                pLegendString = cell(nBatches,1);
                
                for resultsIdx = 1:nBatches
                    %Create legend entry  for polar axis                  
                    if strcmp(results(resultsIdx).batchName(2), ':')
                        pLegendString{resultsIdx} = ['Batch #', results(resultsIdx).batchName(1)];
                    else
                        pLegendString{resultsIdx} = ['Batch #', results(resultsIdx).batchName(1:2)];
                    end
                    
                    %Create and plot polar histogram
                    if resultsIdx == 2
                        hold(ui.pax, 'on')
                    end
                     
                    if ui.btnBars.Value
                        edgeColor1 = 'k';
                        edgeAlpha1 = .25;
                        displayStyle1 = 'bar';
                        faceColor1 = batchColors(resultsIdx,:);
                    elseif ui.btnStairs.Value
                        edgeColor1 = batchColors(resultsIdx,:);
                        faceColor1 = 'none';
                        edgeAlpha1 = 1;
                        displayStyle1 = 'stairs';
                    end
                    
                    ui.polhist = polarhistogram(ui.pax, valuesY{resultsIdx},...
                        'BinMethod','Integer',...
                        'EdgeColor',edgeColor1,...
                        'EdgeAlpha',edgeAlpha1,...
                        'FaceColor',faceColor1,...
                        'DisplayStyle',displayStyle1,...
                        'Normalization',histNorm,...
                        'NumBins',str2double(ui.editBinNum.String));
                    
                    %Save results in currentPlotValues structure
                    currentPlotValues(resultsIdx).anglesCounts = [(ui.polhist.BinEdges(2:end)-(ui.polhist.BinEdges(2)-ui.polhist.BinEdges(1))/2)', ui.polhist.Values'];
                    currentPlotValues(resultsIdx).angles = valuesY;
                    currentPlotValues(resultsIdx).batchName = ui.popBatchSel.String{ui.popBatchSel.Value(1)}(4:end);
                end
                
                hold(ui.pax, 'off')
                if ui.cboxShowLegend.Value
                    legend(ui.pax,'boxoff')
                    legend(ui.pax,pLegendString)
                end         

            case 'swarmplot'
                %% Plot as swarmplot
                                
                %Iterate through selected batches and create labels for the
                %x-axis ticks (=batch file numbers)
                xTicks = zeros(nBatches,1);                
                for idx = 1:nBatches
                    if strcmp(results(idx).batchName(2), ':')
                        xTicks(idx) = str2double(results(idx).batchName(1));
                    else
                        xTicks(idx) = str2double(results(idx).batchName(1:2));
                    end
                end
                
                %Set x-axis label
                xlabel1 = 'Batch file #';
                
                %Retreive the roi sizes in case the user wants to normalize
                %the values by the roi size
                roiSizes = {results(:).roiSize};
                                
                %Get the number of bins entered in the ui
                nBins = str2double(ui.editBinNum.String);
                
                %Catenate all values of all batch files so we know how to
                %set the bin edges
                if ui.cboxNormalizeROI.Value                    
                    allRois = vertcat(roiSizes{:});
                    allValues = vertcat(valuesY{:})./allRois;
                else
                    allValues = vertcat(valuesY{:});
                end
                
                        
                %Calculate bin edges                
                binEdges = linspace(min(allValues),max(allValues(allValues < inf)),nBins+1);
                                
                valuesX = cell(1,nBatches);
                cellMovieNames = cell(1,nBatches);
                cellMovieNumbers = cell(1,nBatches);
                
                
                for resultsIdx = 1:nBatches
                                        
                    if ui.cboxNormalizeROI.Value
                        %Normalize values by roi size
                        currentValuesY = valuesY{resultsIdx}./roiSizes{resultsIdx};
                    else
                        currentValuesY = valuesY{resultsIdx};
                    end
                    
                    %Sort values for nice visualization
                    [currentValuesY, I] = sort(currentValuesY);
                    
                    
                    %Save movie names and movie numbers. This is needed for
                    %showing to which movie a certain value belongs if the
                    %user selects a datapoint with the datatip tool.
                    currentMovieNames = results(resultsIdx).movieNames;
                    currentMovieNumbers = results(resultsIdx).movieNumbers;
                    cellMovieNames{resultsIdx} = currentMovieNames(I);
                    cellMovieNumbers{resultsIdx} = currentMovieNumbers(I);
                             
                    %Initialize x-axis values to the number of the
                    %respective batch file number
                    currentValuesX = ones(numel(currentValuesY),1)*resultsIdx;

                    %Save current x values in cell array
                    valuesX{resultsIdx} = currentValuesX;
                    valuesY{resultsIdx} = currentValuesY;
                    
                    %Create swarm plot using the given y values and their
                    %created x values. Save the values and the respective movie
                    %names and number in user data. This is for later used in
                    %the DataCursorCB function for knowing which datapoint has
                    %been clicked.
                    
                    if resultsIdx == 2
                        hold(ui.ax,'on')
                    end
                    
                    %TODO for future versions: implement matlab swarmchart including XJitterWidth
%                     if verLessThan('matlab','9.9')
                        
                        %Sort values into bins
                        [binNum] = discretize(currentValuesY,binEdges);
                        
                        %Iterate through bins an distribute x-axis values depending
                        %on the amount of values in a given bin
                        for binIdx = 1:nBins
                            %Find which values belong to current bin
                            valuesIdx = find(binNum == binIdx);
                            %Get amount of values in current bin
                            nValues = numel(valuesIdx);
                            
                            if nValues > 1
                                %Calculate distribution width
                                curWidth = 0.7*(1-exp(-0.1*nValues));
                                %Get width of one element
                                widthElement = curWidth / (nValues-1);
                                
                                %Initialize offset depending if current bin has
                                %even or odd amount of values
                                if mod(nValues,2) == 0
                                    offset = widthElement / 2;
                                else
                                    offset = eps;
                                end
                                
                                %Iterate though all values in this bin and
                                %calculate their position on the x-axis
                                for valueIdxCurBin = 1:nValues
                                    %Calculate x-value
                                    currentValuesX(valuesIdx(valueIdxCurBin),1) = resultsIdx + offset;
                                    %Increase offset
                                    offset = offset - sign(offset) * widthElement * valueIdxCurBin;
                                end
                            end
                            
                        end
                        

                        ui.hist = scatter(ui.ax,currentValuesX,currentValuesY,[],batchColors(resultsIdx,:),'.','SizeData',200,'UserData',{cellMovieNames{resultsIdx}; cellMovieNumbers{resultsIdx}});
%                     else
%                         ui.hist = swarmchart(ui.ax,currentValuesX,currentValuesY,[],batchColors(resultsIdx,:),'.','SizeData',200,'XJitterWidth',.4,'UserData',{cellMovieNames{resultsIdx}; cellMovieNumbers{resultsIdx}});
%                     end
                end
                                
                hold(ui.ax,'on')
                
                %Save all created vales in the currentPlotValues structure
                [currentPlotValues(:).batchName] = ui.popBatchSel.String{ui.popBatchSel.Value};
                [currentPlotValues(:).movieNames] = cellMovieNames{:};
                [currentPlotValues(:).movieNumbers] = cellMovieNumbers{:};
                [currentPlotValues(:).movieWiseValues] = valuesY{:};
                
                
                %Set labels and ticks for the x-axis                
                xticks(ui.ax, 1:nBatches)
                xticklabels(xTicks)                
                
                if isTrackedFraction == 1 && ~ui.cboxNormalizeROI.Value
                    %User is looking at tracked fractions
                    
                    %Create legend entries
                    legendString = {'Movie-wise mean + standard error','Pooled fraction'};
                    
                    %Add movie-wise mean value plus standard error as error
                    %bar
                    err(1) = errorbar(ui.ax,(1:nBatches)-0.3,valuesMean,valuesStd,...
                        '.','MarkerSize',20,'Color','r','LineWidth',1.5,'CapSize',8,'Userdata','meanMoviewise');
                    %Add pooled tracked fraction plus error bar
                    err(2) = errorbar(ui.ax,(1:nBatches)+0.3,wholeSet',wholeSetErr',...;
                        '.','MarkerSize',20,'Color','k','LineWidth',1.5,'CapSize',8,'Userdata','wholeSet');
                    
                    %Save all data to currentPlotValues structure
                    cellMean = num2cell(valuesMean);
                    [currentPlotValues(:).mean] = cellMean{:};
                    cellStd = num2cell(valuesStd);
                    [currentPlotValues(:).stdError] = cellStd{:};
                    cellPooled = num2cell(wholeSet);
                    [currentPlotValues(:).pooledTrackedFraction] = cellPooled{:};
                    cellPooledError = num2cell(wholeSetErr);
                    [currentPlotValues(:).pooledTrackedFractionError] = cellPooledError{:};
                else
                    %User is not looking at tracked fractions
                    
                    %Create legend entries
                    legendString = {'Mean + std. dev.'};
                    
                    %Calculate mean values and standard deviations of the values in each batch file                    
                    valuesMean = cellfun(@(x) mean(x, 'omitnan'),valuesY);                      
                    valuesStd = cellfun(@(x) std(x, 'omitnan'),valuesY);
                    
                    %Plot movie-wise mean values with standard deviation as
                    %error bar
                    err(1) = errorbar(ui.ax,(1:nBatches)-0.3,valuesMean,valuesStd,...
                        '.','MarkerSize',20,'Color','r','LineWidth',1.5,'CapSize',8,'Userdata','meanMoviewise');
                    
                    %Save all data to currentPlotValues structure
                    cellMean = num2cell(valuesMean);
                    [currentPlotValues(:).mean] = cellMean{:};
                    cellStd = num2cell(valuesStd);
                    [currentPlotValues(:).stdDev] = cellStd{:};
                end
            case 'valueVsParameter'
                %% Plot values vs. movie number
                
                
                selectedParameter = ui.menuValueVsParameter.String{ui.menuValueVsParameter.Value};

                switch selectedParameter
                    case 'Movie number'
                        
                        valuesX = {results(:).movieNumbers};
                        
                        xlabel1 = 'Movie number';
                        
                    case 'No. of tracks'
                        valuesX = {results(:).nTracks};
                        
                        xlabel1 = 'No. of tracks';
                    case 'No. of non-linked spots'
                        valuesX = {results(:).nNonLinkedSpots};
                        
                        xlabel1 = 'No. of non-linked spots';
                    case 'No. of all events'
                        valuesX = {results(:).nAllEvents};
                        
                        xlabel1 = 'No. of all events';
                    case 'Mean jump distance'
                        valuesX = {results(:).meanJumpDistMoviewise};
                        xlabel1 = 'Average jump distance';
                    case 'Avg. track length'
                        valuesX = {results(:).meanTrackLength};
                        
                        xlabel1 = 'Average track length';
                    case 'Avg. no. of tracks per frame'
                        valuesX = {results(:).meanTracksPerFrame};
                        
                        xlabel1 = 'Average no. of tracks per frame';
                    case 'Avg. no. of spots per frame'
                        valuesX = {results(:).meanSpotsPerFrame};
                        
                        xlabel1 = 'Average no. of spots per frame';
                    case 'ROI size'
                        valuesX = {results(:).roiSize};
                        
                        xlabel1 = 'ROI size';
                    case 'No. of jumps'
                        valuesX = {results(:).nJumps};
                        
                        xlabel1 = 'No. of jumps';
                end

                %Initialize legend entries variable
                legendString = cell(nBatches,1);
                
                %Iterate through selected batches
                for resultsIdx = 1:nBatches
                    
                    %Create legend entry
                    if strcmp(results(resultsIdx).batchName(2), ':')
                        batchNum = results(resultsIdx).batchName(1);
                    else
                        batchNum = results(resultsIdx).batchName(1:2);
                    end
                    
                    currentMovieNames = results(resultsIdx).movieNames;
                    currentMovieNumbers = results(resultsIdx).movieNumbers;
                    
                    legendString{resultsIdx} = ['Batch #', batchNum];
                    
                    curMovieValuesX = valuesX{resultsIdx};
                    curMovieValuesY = valuesY{resultsIdx};
                                     
                    if ui.cboxNormalizeROI.Value    
                        curMovieRoiSize = results(:).roiSize;
                        curMovieValuesY = curMovieValuesY./curMovieRoiSize;
                    end
                    
                                                  
                    if resultsIdx == 2
                        hold(ui.ax,'on')
                    end
                    
                    %Plot histogram
                    ui.hist = scatter(ui.ax,curMovieValuesX,curMovieValuesY,[],batchColors(resultsIdx,:),'.','SizeData',100,'UserData',{currentMovieNames; currentMovieNumbers'});
                    
                    currentPlotValues(resultsIdx).batchName = results(resultsIdx).batchName;
                    currentPlotValues(resultsIdx).movieNames = results(resultsIdx).movieNames;
                    currentPlotValues(resultsIdx).values = curMovieValuesY;
                    
                end                    
            case 'scatterplot'
                %% Plot as scatterplot
                
                if strcmp(selectedColor, 'standard')
                    colormap('default')
                else
                    colormap(selectedColor)
                end
                
                if ~isempty(valuesY) && any(valuesY)
                    %Scatterplot is plotted directly in the AdjustAxisLimits()
                    %function so that density calculations are renewed when
                    %the user changes the axis limits.
                    
                    ui.hist = [valuesX, valuesY];
                else
                    ui.hist = [];
                end
                
                colorbar('off')                
            case 'diffusion parameters'
                %% Plot as barchart
                
                %Make sure axis is not set to logarithmic
                ui.cboxLogX.Value = 0;
                                
                %Get frame cycle times as we need them to calculate the
                %diffusion constants
                frameCycleTimes = {results(:).frameCycleTimes};
                
                %Initialize lengend entries
                legendString = cell(numel(nBatches),1);
                
                %Initialize cell arrays
                movieOrResamplingValues = cell(nBatches,1);
                resamplingValues = cell(nBatches,1);
                
                %Initialize variable that is set to true if a batch file contains different
                %tracking radii
                showWarnDlg = false;
                
                %Retreive start values for fitting
                startD = vertcat(ui.tableStartD.Data{:,2});
                
                %Iterate through all selected batches and fit diffusion
                %model
                for resultsIdx = 1:nBatches
                    
                    %Create legend entry
                    if strcmp(results(resultsIdx).batchName(2), ':')
                        batchNum = results(resultsIdx).batchName(1);
                    else
                        batchNum = results(resultsIdx).batchName(1:2);
                    end
                    
                    legendString{resultsIdx} = ['Batch #', batchNum];
                    
                    if ui.btnPxFr.Value
                        %User wants results in pixels and frames
                        frameCycleTime = 1;
                    else                      
                        %User want results in microns and seconds
                        
                        %Get frame cycle times in movies of current batch
                        curframeCycleTimes = unique(frameCycleTimes{resultsIdx});
                        %Convert frame cycle time to seconds. If more than
                        %one frame cycle time is found, use first one
                        frameCycleTime = curframeCycleTimes(1);
                    end
                    
                    %Get tracking radius
                    allTrackingRadii = vertcat(results(resultsIdx).trackingRadii);
                    trackingRadius = max(unique(allTrackingRadii));
                    trackingRadius = trackingRadius^2/(4*frameCycleTime);
                    
                                        
                    %Create pooled jump histogram including all movies of current batch
                    jumpDistsPooled = vertcat(results(resultsIdx).jumpDistances{:});
                    
                    %Check if 1, 2 or 3 exponential rates are selected
                    if ui.btnOneRate.Value
                        nRates = 1;
                    elseif ui.btnTwoRates.Value
                        nRates = 2;
                    elseif ui.btnThreeRates.Value
                        nRates = 3;
                    end
                                        
                    if ui.popError.Value == 1
                        [pooledY,pooledX] = histcounts(jumpDistsPooled,'BinMethod','integers',...
                            'Normalization','cdf',...
                            'NumBins',str2double(ui.editBinNum.String));
                        
                        %Get bin centers and square the x-axis to calculate the diffusion coefficient
                        pooledX = (pooledX(2:end)-(pooledX(2)-pooledX(1))/2).^2;
                        pooledX = pooledX./(4*frameCycleTime);
                        
                        %Fit the curve
                        fitResults(resultsIdx) = dispfit_cumulative(pooledX',pooledY', trackingRadius, startD, nRates);
                        
                        %Save results
                        results(resultsIdx).diffParams = fitResults(resultsIdx);
                        
                        nTrackingRadii = numel(unique(allTrackingRadii));
                        
                        if nTrackingRadii > 1
                            showWarnDlg = true;
                        end
                    elseif ui.popError.Value == 2 %Show movie-wise values
                        
                        %User wants to see movie-wise fitted results
                        
                        %Iterate through each movie of current batch
                        for movieIdx = 1:length(results(resultsIdx).movieNames)
                            %Get jump distances of current movie
                            dispsCurMovie = results(resultsIdx).jumpDistances{movieIdx};
                            if ~isempty(dispsCurMovie)
                                
                                %Create cumulative density function
                                [y,x] = histcounts(dispsCurMovie,'BinMethod','integers',...
                                    'Normalization','cdf',...
                                    'NumBins',str2double(ui.editBinNum.String));
                                
                                
                                %Get bin centers and square the x-axis to calculate the diffusion coefficient
                                x = (x(2:end)-(x(2)-x(1))/2).^2;
                                x = x./(4*frameCycleTime);
                                
                                %Create fit results
                                fitResults{resultsIdx}(movieIdx) = dispfit_cumulative(x',y', trackingRadius, startD, nRates);
                                
                                %Adjust fit results with respect to the given
                                %frame cycle time
                            else
                                fitResults{resultsIdx}(movieIdx)=struct('D',[0 0 0],'Derr',[0 0 0], 'A', [0 0 0], 'Aerr',[0 0 0],'EffectiveD',0,'Ajd_R_square',0,'Message',0,'SSE',0,'xy',[0 0]);
                            end
                        end
                        
                        %Save results
                        results(resultsIdx).diffParams = fitResults{resultsIdx};
                    elseif ui.popError.Value == 3 %Show resampling
                        ui.editNResampling.BackgroundColor = 'r';
                        drawnow
                        rng('default')
                        %User also wants to see movie-wise fitted results
                        nResampling = str2double(ui.editNResampling.String);
                        percResampling = str2double(ui.editPercResampling.String);
                        
                        %Iterate through each movie of current batch                        
                        nJumps = numel(jumpDistsPooled);
                                                
                        for resamplingIdx = 1:nResampling
                            %Get jump distances of current movie
                            randInd = randperm(nJumps, round(nJumps*percResampling/100));
                            
                            jumpDistRand = jumpDistsPooled(randInd);
                            
                            %Create cumulative density function
                            [y,x] = histcounts(jumpDistRand,'BinMethod','integers',...
                                'Normalization','cdf',...
                                'NumBins',str2double(ui.editBinNum.String));
                            
                            %Square x-axis
                            x = (x(2:end)-(x(2)-x(1))/2).^2;
                            x = x./(4*frameCycleTime);
                            
                            %Create fit results
                            fitResults{resultsIdx}(resamplingIdx) = dispfit_cumulative(x',y', trackingRadius, startD, nRates);
                        end
                        
                        %Save results
                        results(resultsIdx).diffParams = fitResults{resultsIdx};
                        
                        nTrackingRadii = numel(unique(allTrackingRadii));
                        
                        if nTrackingRadii > 1
                            showWarnDlg = true;
                        end
                    end
                    
                end
                
                if ui.btnShowD.Value      
                    %User wants to see diffusion constants
                    
                    %Create label for y-axis
                    if ui.btnPxFr.Value
                        ylabel1 = 'Diffusion coefficient (px^2/frame)';
                    else
                        ylabel1 = 'Diffusion coefficient (\mum^2/sec)';
                    end
                                        
                    %Create x-tick labels
                    if ui.btnOneRate.Value                        
                        xTickLabel1 = {'D_1'};
                    elseif ui.btnTwoRates.Value
                        xTickLabel1 = {'D_1','D_2'};
                    elseif ui.btnThreeRates.Value
                        xTickLabel1 = {'D_1','D_2','D_3'};
                    end
                    
                    if ui.popError.Value == 1
                        %Prepare pooled confidence intervall and diffusion
                        %constants for display
                        error = vertcat(fitResults(:).Derr);
                        yValues = vertcat(fitResults(:).D);
                    else
                        %Prepare movie-wise or resampling values for display
                        for resultsIdx = 1:length(fitResults)
                            movieOrResamplingValues{resultsIdx} = vertcat(fitResults{resultsIdx}(:).D);
                        end
                        
                        yValues = cellfun(@mean, movieOrResamplingValues, 'UniformOutput', false);
                        yValues = vertcat(yValues{:});
                        error = cellfun(@std, movieOrResamplingValues, 'UniformOutput', false);
                        error = vertcat(error{:});
                    end
                               
                elseif ui.btnShowA.Value   
                    %User wants to see amplitudes
                    
                    %Create label for y-axis
                    ylabel1 = 'Fraction';
                     
                    %Create x-ticj labels
                    if ui.btnOneRate.Value
                        xTickLabel1 = {'A_1','A_2'};
                    elseif ui.btnTwoRates.Value
                        xTickLabel1 = {'A_1','A_2'};
                    elseif ui.btnThreeRates.Value
                        xTickLabel1 = {'A_1','A_2','A_3'};
                    end
                    
                    if ui.popError.Value == 1
                        %Prepare pooled confidence intervall and diffusion
                        %constants for display
                        error = vertcat(fitResults(:).Aerr);
                        yValues = vertcat(fitResults(:).A);
                    else
                        %Prepare movie-wise or resampling values for display
                        for resultsIdx = 1:length(fitResults)
                            movieOrResamplingValues{resultsIdx} = vertcat(fitResults{resultsIdx}(:).A);
                        end
                        
                        yValues = cellfun(@mean, movieOrResamplingValues, 'UniformOutput', false);
                        yValues = vertcat(yValues{:});
                        error = cellfun(@std, movieOrResamplingValues, 'UniformOutput', false);
                        error = vertcat(error{:});
                    end
                    
                elseif ui.btnShowEffectiveD.Value   
                    %User wants to see Deff
                    
                    %Create label for y-axis
                    if ui.btnPxFr.Value
                        ylabel1 = 'Effective diffusion coefficient (px^2/frame)';
                    else
                        ylabel1 = 'Effective diffusion coefficient (\mum^2/sec)';
                    end
                    
                    %Create x-axis label
                    xTickLabel1 = {''};
               
                    if ui.popError.Value == 1
                        %Prepare pooled confidence intervall and diffusion
                        %constants for display
                        error = 0;
                        yValues = vertcat(fitResults(:).EffectiveD);
                    else
                        %Prepare movie-wise or resampling values for display
                        for resultsIdx = 1:length(fitResults)
                            movieOrResamplingValues{resultsIdx} = vertcat(fitResults{resultsIdx}(:).EffectiveD);
                        end
                        
                        yValues = cellfun(@mean, movieOrResamplingValues, 'UniformOutput', false);
                        yValues = vertcat(yValues{:});
                        error = cellfun(@std, movieOrResamplingValues, 'UniformOutput', false);
                        error = vertcat(error{:});
                    end
                                           
                end
                
                if size(yValues,2) == 1 && size(yValues,1) > 1
                    %Dirty workaround for displaying only one diffusion
                    %coefficient in a bar chart for Matlab versions < 2020a
                    yPlotValues = yValues;
                    yPlotValues(1,2) = 0;
                    
                    %Create bar chart
                    hBar = bar(ui.ax, 1:2, yPlotValues', 0.8 ,'FaceColor','flat');
                    
                    for idx = 1:numel(hBar)
                        hBar(idx).XData = hBar(idx).XData(1);
                        hBar(idx).YData = hBar(idx).YData(1);
                    end
                else
                    %Create bar chart
                    hBar = bar(ui.ax, 1:size(yValues,2), yValues', 0.8 ,'FaceColor','flat');
                end
                                
                hold(ui.ax, 'on')
                
                %Set bar colors and display moviewise or resampling values
                for k1 = 1:nBatches                    
                    hBar(k1).CData = batchColors(k1,:);
                    
                    %Get centers of bar chart 
                    center(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');      % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
                               
                    if ui.popError.Value > 1
                        %Add movie-wise values to bar chart
                        movieOrResamplingX = repmat(center(k1,:),size(movieOrResamplingValues{k1},1),1);
                        movieOrResamplingY = movieOrResamplingValues{k1};
                        plot(ui.ax, movieOrResamplingX,movieOrResamplingY, '.', 'Color','k');
                    end
                end
                
                
                if error ~= 0
                    errorbar(ui.ax, center, yValues, error, '.k', 'Capsize', 20);
                end
                
                %Save results in currentPlotValues variable
                switch ui.popError.Value
                    case 1
                        currentPlotValues = struct('batchName',{ui.popBatchSel.String(ui.popBatchSel.Value)},'pooledValues',yValues,'confidenceIntervall', error);
                    case 2
                        currentPlotValues = struct('batchName',{ui.popBatchSel.String(ui.popBatchSel.Value)},'meanValue',yValues,'stdMoviewise', error,'moviewiseValues',{movieOrResamplingValues});
                    case 3
                        currentPlotValues = struct('batchName',{ui.popBatchSel.String(ui.popBatchSel.Value)},'meanValue',yValues,'stdResampling', error,'resamplingValues',{movieOrResamplingValues});
                end
                                  
                hold(ui.ax, 'off')
                
                if showWarnDlg
                    warndlg(['At least one batch file contains '...
                        'movies that have been analyzed with different '...
                        'tracking radii. The tracking radius is used to '...
                        'normalize the fitting function (see manual). '...
                        'Using largest tracking radius in batch for normalization.'],'Warning');
                end
                
                %Set x-tick labels 
                set(ui.ax,'XTickLabel',xTickLabel1)
                
                %Set empty x axis labels
                xlabel1 = {};
                ui.editNResampling.BackgroundColor = [1 1 1];
        end
        
        if ui.cboxShowLegend.Value && ~strcmp(plotStyle, 'scatterplot')
            %Show legend
            if strcmp(plotStyle, 'swarmplot')
                legend(err,legendString)
            else
                legend(ui.ax,legendString)
            end
            legend('boxoff')
        else
            %Hide legend
            legend(ui.ax,'off')
            legend(ui.pax,'off')
        end
        
        if ~isempty(ylabel1) || ~isempty(xlabel1)
            %Set x and y-axis labels
            xlabel(ui.ax,xlabel1)
            ylabel(ui.ax,ylabel1)
        end
        hold(ui.ax,'off')
        
        
        AdjustAxisLimits()
    end

    function AdjustAxisLimits()
        %Executed either when the axis limits have changed or after the
        %plot was updated by PlotHistogramCB function
        
        
        
        %Check if y-axis has to be plotted logarithmic or linear
        if ui.cboxLogY.Value
            set(ui.ax,'YScale','log')
        else
            set(ui.ax,'YScale','linear')
        end
        
        %Check if x-axis has to be plotted logarithmic or linear
        if ui.cboxLogX.Value
            set(ui.ax,'XScale','log')
        else
            set(ui.ax,'XScale','linear')
        end
        
        nBatches = numel(ui.popBatchSel.Value);
            
        if ui.cboxAutoX.Value
            %Auto adjust x-axis selected so calculate limits
            if strcmp(plotStyle, 'swarmplot')
                %Plot style is swarmplot so give a little more space to the
                %left and right
                ui.editLimX1.String = 0.5;
                ui.editLimX2.String = nBatches+0.5;
            elseif strcmp(plotStyle, 'scatterplot')
                if ~isempty(ui.hist)
                    %Plot style is scatterplot
                    ui.editLimX1.String = min(ui.hist(:,1));
                    ui.editLimX2.String = max(ui.hist(:,1));
                end
            else
                axis(ui.ax,'tight')
                xLimits = xlim(ui.ax);
                if isnumeric(xLimits)
                    ui.editLimX1.String = round_significant(xLimits(1),2,'floor');
                    ui.editLimX2.String = round_significant(xLimits(2),2,'ceil');
                else
                    ui.editLimX1.String = xLimits(1);
                    ui.editLimX2.String = xLimits(2);
                end
            end
            
        end
        
        if ui.cboxAutoY.Value
            if strcmp(plotStyle, 'scatterplot')
                %Plot style is scatterplot
                if ~isempty(ui.hist)
                    ui.editLimY1.String = min(ui.hist(:,2));
                    ui.editLimY2.String = max(ui.hist(:,2));
                end
            else
                %Auto adjust y-axis selected so calculate limits
                axis(ui.ax,'tight')
                yLimits = ylim(ui.ax);
                
                ui.editLimY1.String = round_significant(max(0,yLimits(1)),12,'floor');
                ui.editLimY2.String = round_significant(yLimits(2),2,'ceil');
            end
        end
        
        
        if strcmp(plotStyle, 'scatterplot')
            %Plot scatterplot with density heat-map using the user defined axis limits
            
            if ~isempty(ui.hist)
                valuesX = ui.hist(:,1);
                valuesY = ui.hist(:,2);
                
                
                xIdx = valuesX < str2double(ui.editLimX2.String);
                yIdx = valuesY < str2double(ui.editLimY2.String);
                
                allIdx = xIdx & yIdx;
                valuesX = valuesX(allIdx);
                valuesY = valuesY(allIdx);
                scatplot(valuesX,valuesY);
                colorbar('off')
            else
                ui.hist = scatter(ui.ax,NaN,NaN);
            end
        end
                
        %Set x -and y-axis limits
        xlim(ui.ax,[str2double(ui.editLimX1.String) str2double(ui.editLimX2.String)]);
        ylim(ui.ax,[str2double(ui.editLimY1.String) str2double(ui.editLimY2.String)]);
        
    end

    function cursorText = dataCursorCB(~,eventHandle)
        %Executed when the cursor is hovered over a data entry or when a data
        %point is selected
        
        graphObjHandle = get(eventHandle,'Target');
        pos = get(eventHandle,'Position');
               
        if strcmp(plotStyle, 'histogram') || strcmp(plotStyle, 'angular histogram')
            %Plotstyle is histogram or angular histogram so display the
            %value and the bin edges
            
            upperEdgeIdx = find(pos(1) <= graphObjHandle.BinEdges,1,'first');
            binEdges = graphObjHandle.BinEdges(upperEdgeIdx-1:upperEdgeIdx);
            cursorText = {['Value: ',num2str(pos(2))],...
                    ['Bin edges: [', num2str(binEdges(1)), ' ',num2str(binEdges(2)), ']']};
        elseif strcmp(plotStyle, 'swarmplot') || strcmp(plotStyle, 'valueVsParameter')
            %Plotstyle is swarmplot
            
            %Check which datapoint has been selected
            if strcmp(graphObjHandle.UserData,'meanMoviewise')
                if find(ui.popParamSel.Value == [8,9,10])
                cursorText = {'Mean of moviewise values',...
                    ['Mean: ', num2str(graphObjHandle.YData(ceil(pos(1))))],...
                    ['Standard error: ', num2str(graphObjHandle.YPositiveDelta(ceil(pos(1))))]};
                else
                    cursorText = {'Mean of moviewise values',...
                    ['Mean: ', num2str(graphObjHandle.YData(ceil(pos(1))))],...
                    ['Standard deviation: ', num2str(graphObjHandle.YPositiveDelta(ceil(pos(1))))]};
                end
            elseif strcmp(graphObjHandle.UserData,'wholeSet')
                cursorText = {'Tracked fraction all movies pooled',...
                    ['Tracked fraction: ', num2str(graphObjHandle.YData(floor(pos(1))))],...
                    ['Error: \pm', num2str(graphObjHandle.YPositiveDelta(floor(pos(1))))]};
            else
                
                movieNames = graphObjHandle.UserData{1};
                movieNumbers = graphObjHandle.UserData{2};
                
                index = graphObjHandle.XData == pos(1) & graphObjHandle.YData == pos(2);

                movieNumber = movieNumbers(index);
                movieName = char(movieNames{index});
                movieName(strfind(movieName, '_')) = ' ';
                
                cursorText = {['Value: ',num2str(pos(2))],...
                    ['Movienumber: ',num2str(movieNumber)],...
                    ['Filename: ', movieName]};
            end
        else
            %Plotstyle is either diffusion parameter or scatterplot so just
            %show the x and y value of the selected datapoint
            cursorText = {['X: ',num2str(pos(1))],...
                ['Y: ',num2str(pos(2))]};
        end
    end

%Load/remove batch file, select batch file, tl condition or movie number

    function LoadBatchFilesCB(src,~)
        %Executed when "Load .mat batch file(s)" is pressed
        
        %Open file dialog box
        [fileNameList,pathName] = uigetfile('*.mat','Select .mat batch file(s)','MultiSelect','on',currentBatchPath);
        
        if isequal(fileNameList,0)  
            %User didn't choose a file
            return
        elseif ~iscell(fileNameList) 
            %Check if only one file has been chosen
            fileNameList = {fileNameList};
        end
        
        %Save selected path for using it as the starting path in the next
        %file dialog box
        currentBatchPath = pathName;
        
        %Number of batches to load
        nNewBatches = length(fileNameList);
        
        %Number of batches currently loaded in the track analyser
        nOldBatches = length(batches);
        
        %Iterate through files and load each batch file
        newBatchFiles = cell(nNewBatches,1);
        
        for fileIdx = 1:nNewBatches
            %Monitor progress
            src.String = ['Loading ', num2str(fileIdx), ' of ', num2str(nNewBatches)];
            drawnow
            curBatchFile = fullfile(pathName,fileNameList{fileIdx});
            loadedBatch = load(curBatchFile);
            loadedBatch = loadedBatch.batch; 
            
            newBatchFiles{fileIdx} = loadedBatch;
            %Add number to the list of filenames to later display the filenames in the ui
            fileNameList{fileIdx} = [num2str(nOldBatches+fileIdx),': ', fileNameList{fileIdx}];
        end
        
        %Catenate the old and new batches
        batches = vertcat(batches, newBatchFiles);
        
        %List containing unique frame cycle times of each batch
        frameCycleTimesList = cell(length(batches),1);
        %List containing all the frame cycle times of all movies of each batch
        frameCycleTimeMovieList = cell(length(batches),1);
        
        %Iterate through all batches and create lists of frame cycle times
        for batchIdx = 1:length(batches)
           currentBatch = batches{batchIdx};
            
            frameCycleTimeMovieList{batchIdx} = zeros(length(currentBatch),1);
            for movieIdx = 1:length(currentBatch)
                frameCycleTimeMovieList{batchIdx}(movieIdx) = currentBatch(movieIdx).movieInfo.frameCycleTime;
            end
            
            frameCycleTimesList{batchIdx} = unique(frameCycleTimeMovieList{batchIdx}); 
        end
        
         
        %Reset String of the load button
        src.String = 'Load batch .mat file(s)';
        %Set selected batch file to first entry
        ui.popBatchSel.Value = 1;
        %Update list of batch files
        ui.popBatchSel.String = [ui.popBatchSel.String; fileNameList'];
        BatchSelectionCB()
        
    end

    function RemoveBatchFilesCB(~,~)
        %Executed when "remove selected file(s)" button is pressed
        
        %Remove selected batch file names from ui list
        ui.popBatchSel.String(ui.popBatchSel.Value) = [];
        
        %Create new list of batch files
        for idx = 1:length(ui.popBatchSel.String)
            ui.popBatchSel.String{idx} = [num2str(idx),': ', ui.popBatchSel.String{idx}(4:end)];
        end
        
        %Remove selected batches from the frame cycle time lists and from
        %the batches structure
        frameCycleTimesList(ui.popBatchSel.Value) = [];
        frameCycleTimeMovieList(ui.popBatchSel.Value) = [];        
        batches(ui.popBatchSel.Value) = [];
        
        %Set seleted batch file to the first entry
        ui.popBatchSel.Value = 1;  
        
        %Selection of batch files changed so call BatchSelectionCB function
        BatchSelectionCB()
    end

    function BatchSelectionCB()
        %Executed whenever the selection of batch files changes
        
        if isempty(batches)
            %No batch files in the list
            ui.popTlSel.String = {};
            ui.editMovie.String = 1;
            ui.textMovie3.String = 0;
            ui.tableStatistics.ColumnName = {''};
            ui.tableStatistics.Data = {'#movies';'#tracks';'#non-linked spots';'#all events'};            
        else
            %Get selected batch files
            batchIdx = ui.popBatchSel.Value;
            
            %Create unique list of frame cycle times accouring in the
            %current selection of batch files
            allTl = [];
            nFilesPerBatch = [];
            for fileIdx = batchIdx
                allTl = [allTl; frameCycleTimesList{fileIdx}];
                nFilesPerBatch = [nFilesPerBatch; length(batches{fileIdx})];
            end
            
            allTl = unique(allTl);
            
            %Create cell array containing all frame cycle times plus
            %"single" and "all movies" in the list of frame cycle times
            
            tlList = {'Single movie', 'All movies'};
            
            for k = 1:numel(allTl)
                tlList{k+2} = [num2str(allTl(k)), ' ms'];
            end
            
            %Update list of frame cycle times in the ui
            ui.popTlSel.String = tlList;
            
            %Set maximum amount of movies in the single movie selection
            %panel to the maximum number of movies in one batch
            ui.textMovie3.String = max(nFilesPerBatch);
            
            %Take care that the current movie number is not higher than the
            %maximum amount of movies in a batch
            if str2double(ui.editMovie.String) > max(nFilesPerBatch)
                curMovieIndex = max(nFilesPerBatch);
                ui.editMovie.String = max(nFilesPerBatch);
            end
            
            %Take care that the selected value in the frame cycle time list 
            %is not higher than the amount of enries in this list
            if ui.popTlSel.Value(end)-2 > numel(allTl)
                ui.popTlSel.Value = length(tlList);
            end
        end
        
        %Update all results and plots
        UpdateRegionsList()
        CreateData()
    end

    function TlSelectionCB()
        %Executed when the user interacts with the frame cycle time list or
        %region list
                
        if ui.popTlSel.Value ==1
            %"Single movie" selected so show movie number ui elements
            ui.panelMovieSel.Visible = 'on';
        else
            %Hide movie number ui elements
            ui.panelMovieSel.Visible = 'off';
        end
        
        nRegionsSelected = numel(ui.popRegionSel.Value);
        
        if ui.popRegionSel.Value(1) == 1 && nRegionsSelected > 1
            ui.popRegionSel.Value = 1;            
        end
        
        if nRegionsSelected > 1
            ui.cboxNormalizeROI.Value = 0;
        end
        
        UpdateRegionsList()
        CreateData()
        
        if nRegionsSelected > 1
            ui.cboxNormalizeROI.Visible = 'off';
        end
    end

    function MovieNumberCB(src,~)
        %Executed when the movie number is changed by the user
        
        previousMovieNumber = curMovieIndex;
        
        if strcmp(src.String,'Previous') && previousMovieNumber > 1
            curMovieIndex =  previousMovieNumber - 1;
        elseif strcmp(src.String,'Next') && previousMovieNumber < str2double(ui.textMovie3.String)
            curMovieIndex = previousMovieNumber + 1;
        elseif str2double(src.String) <= str2double(ui.textMovie3.String) && str2double(src.String) > 0
            curMovieIndex = str2double(src.String);
        end
                
        ui.editMovie.String = curMovieIndex;
        
        UpdateRegionsList()
        CreateData()
    end

    function UpdateRegionsList()
        %Sub-region analysis will be available in following versions and is
        %currently under developement
        
        %Updates the list of displayed sub-regions
        
        if isempty(batches)
            ui.popRegionSel.String = {};
            return
        end
                
        nMaxRegions = 0;
        
        %Iterate through selected batch files and find the maximum amount
        %of sub-regions a movie contains
        for batchSelVal = ui.popBatchSel.Value
            
            currentBatch = batches{batchSelVal};
            
            %Current movie, all movies or selection of movies
            if ui.popTlSel.Value(1) == 1 %Current movie
                if length(batches{batchSelVal}) < curMovieIndex
                    continue
                else
                    tlMoviesIdx = curMovieIndex;
                end
            elseif ui.popTlSel.Value(1) == 2 %All movies
                tlMoviesIdx = 1:length(currentBatch);
            else %Specific TL
                tlMovies = zeros(length(currentBatch),1);
                for k = 1:numel(ui.popTlSel.Value)
                    curTlCond = str2double(ui.popTlSel.String{ui.popTlSel.Value(k)}(1:end-3));
                    tlMovies = or(tlMovies, frameCycleTimeMovieList{batchSelVal} == curTlCond);
                end
                tlMoviesIdx = find(tlMovies');
            end

            for movieID = tlMoviesIdx
                nMaxRegions = max(nMaxRegions, batches{batchSelVal}(movieID).results.nSubRegions+1);
            end
        end
        
        if ui.popRegionSel.Value > nMaxRegions+1
            ui.popRegionSel.Value = 1;
        end
        
        %Create cell array containing the displayed region strings

        if nMaxRegions == 1
            ui.popRegionSel.Visible = 'off';
            regionList{1} = 'All regions';
        else            
            ui.popRegionSel.Visible = 'on';
            
            regionList = cell(nMaxRegions,1);
            for regionIdx = 0:nMaxRegions
                if regionIdx == 0
                    regionList{regionIdx+1} = 'All regions';
                elseif regionIdx == 1
                    regionList{regionIdx+1} = ['Region ',num2str(regionIdx), ' (tracking-region)'];
                else
                    regionList{regionIdx+1} = ['Region ',num2str(regionIdx)];
                end
            end
        end
        
        %Update regions list
        ui.popRegionSel.String = regionList;
    end

%Other small stuff

    function OverlayFitWithHistCB(~,~)
        
        answer = inputdlg('Enter number of bins displayed in the jump distance histogram (this is just for visualization and not for fitting!).','Input');
        
        nBins = str2double(answer);

        
        %Get number of selected batch files
        nBatches = length(results);
        
        selectedColor = ui.popLut.String{ui.popLut.Value};
        %Create colormap for plotting
        switch selectedColor
            case 'standard'
                batchColors = custom_colormap(nBatches);
            case 'jet'
                batchColors = jet(nBatches);
            case 'copper'
                batchColors = copper(nBatches);
            case 'parula'
                batchColors = parula(nBatches);
            case 'winter'
                batchColors = winter(nBatches);
            case 'gray'
                batchColors = gray(nBatches+1);
                batchColors = batchColors(1:end-1,:);
        end
        
        %Iterate through selected batches
        for resultsIdx = 1:nBatches
            
            %Create Titel
            if strcmp(results(resultsIdx).batchName(2), ':')
                figTitle = ['Batch #', results(resultsIdx).batchName(1)];
            else
                figTitle = ['Batch #', results(resultsIdx).batchName(1:2)];
            end
            
            figure('Name',figTitle)
            
            curBatchJumps = vertcat(results(resultsIdx).jumpDistances{:});
            
            %Get maximum value in the dataset and create histogram
            %limits
            histLimits = [min(0, min(curBatchJumps)) max(curBatchJumps)];
            
            %Create and plot histogram for current batch
            hist = histogram(curBatchJumps,...
                'EdgeColor','k',...
                'EdgeAlpha',.25,...
                'DisplayStyle','bar',...
                'FaceColor',batchColors(resultsIdx,:),...
                'BinLimits',histLimits,...
                'Normalization','Probability',...
                'NumBins',nBins);
            
            
            %Create cumulative density function from valuesY for
            %current batch
            
            [y,edges] = histcounts(curBatchJumps,...
                'Normalization','cdf',...
                'NumBins',str2double(ui.editBinNum.String));
            
            %Get bin centers
            x = (edges(2:end)-(edges(2)-edges(1))/2);
            
            
            %Get frame cycle times as we need them to calculate the
            %diffusion constants
            frameCycleTimes = {results(:).frameCycleTimes};
            
            if ui.btnPxFr.Value
                %User wants results in pixels and frames
                frameCycleTime = 1;
            else
                %User want results in microns and seconds
                
                %Get frame cycle times in movies of current batch
                curframeCycleTimes = unique(frameCycleTimes{resultsIdx});
                %Convert frame cycle time to seconds. If more than
                %one frame cycle time is found, use first one
                frameCycleTime = curframeCycleTimes(1);
            end
            
            %Square x axis for diffusion fit
            xSq = x.^2;
            xSq = xSq./(4*frameCycleTime);
            
            %Retreive start values for fitting
            startD = vertcat(ui.tableStartD.Data{:,2});
            
            %Get tracking radius
            allTrackingRadii = vertcat(results(resultsIdx).trackingRadii);
            trackingRadius = max(unique(allTrackingRadii));
            trackingRadius = trackingRadius^2/(4*frameCycleTime);
            
            %Get number of diffusive species
            if ui.btnOneRate.Value
                nRates = 1;
            elseif ui.btnTwoRates.Value
                nRates = 2;
            elseif ui.btnThreeRates.Value
                nRates = 3;
            end
            
            %Fit curve with n-exp diffusion fit
            outDiff = dispfit_cumulative(xSq',y', trackingRadius, startD, nRates);
            xEdges = hist.BinEdges;
            deltaX = xEdges(2)-xEdges(1);
            
            D = outDiff.D;
            A = outDiff.A;
            
            %Create y-values
            switch nRates
                case 1
                    y = (1/(2*frameCycleTime))*deltaX*x.*(A(1)/D(1)*exp(-x.^2/(4*frameCycleTime*D(1))));
                case 2
                    y = (1/(2*frameCycleTime))*deltaX*x.*(A(1)/D(1)*exp(-x.^2/(4*frameCycleTime*D(1)))+A(2)/D(2)*exp(-x.^2/(4*frameCycleTime*D(2)))/(1-exp(-trackingRadius/(D(2)))));
                case 3
                    y = (1/(2*frameCycleTime))*deltaX*x.*(A(1)/D(1)*exp(-x.^2/(4*frameCycleTime*D(1)))+A(2)/D(2)*exp(-x.^2/(4*frameCycleTime*D(2)))+A(3)/D(3)*exp(-x.^2/(4*frameCycleTime*D(3)))/(1-exp(-trackingRadius/(D(3)))));
            end
            
            %Plot complete fit curve
            hold on
            plot(x, y,'Color','k','linewidth',2)
            
            %Plot first component
            y = (1/(2*frameCycleTime))*deltaX*x.*(A(1)/D(1)*exp(-x.^2/(4*frameCycleTime*D(1))));
            plot(x, y,'r','linewidth',2)
            
            if nRates > 1
                %Plot second component
                y = (1/(2*frameCycleTime))*deltaX*x.*(A(2)/D(2)*exp(-x.^2/(4*frameCycleTime*D(2)))/(1-exp(-trackingRadius/(D(2)))));
                plot(x, y,'g','linewidth',2)
            end
            if nRates > 2
                %Plot thrid component
                y = (1/(2*frameCycleTime))*deltaX*x.*(A(3)/D(3)*exp(-x.^2/(4*frameCycleTime*D(3)))/(1-exp(-trackingRadius/(D(3)))));
                plot(x, y,'y','linewidth',2)
            end
            
            if ui.btnPxFr.Value
                ylabel('Probability');
                xlabel('Jump distance (px)');
            else
                ylabel('Probability');
                xlabel('Jump distance (\mum)');
            end
            
        end
    end

    function bgITMselectionCB(~,~)
        %Executed when user switches between "ITM" and "continuous" in the track fractions tab
        
        %editNBrightFrames states the amount of bright frames separated by a long dark time.
        %Currently not used but can be changed if an itm scheme contains
        %more than 2 frames in a row
        if ui.btnITM.Value            
            ui.editNBrightFrames.String = 2;
            ui.textNDarkForLong.String = '#survived dark periods to count as long track';
        elseif ui.btnContinuous.Value
            %ui.editNBrightFrames must be one for continuous illumination
            ui.editNBrightFrames.String = 1;
            ui.textNDarkForLong.String = 'Count as long track if number of survived frames is greater than:';
        end
        
        CreateData()
    end

    function btnGroupUnitsCB(~,~)
        %Executed when the units are changed or the pixelsize has been
        %changed
        
        
        if ui.btnPxFr.Value
            %User wants to display results in pixels and frames
            ui.textPixelsize.Visible = 'off';
            ui.editPixelsize.Visible = 'off';
            ui.txtAnglesMinJumpDist.String = 'px';
            ui.txtAnglesMaxJumpDist.String = 'px';
        elseif ui.btnMiMs.Value            
            %User wants to display results in microns and seconds so show
            %field where pixelsize can be entered
            ui.textPixelsize.Visible = 'on';
            ui.editPixelsize.Visible = 'on';
            ui.txtAnglesMinJumpDist.String = 'µm';
            ui.txtAnglesMaxJumpDist.String = 'µm';
        end
        
        %Save pixelsize in variable for later returning it to TrackIt when
        %figure is closed
        pixelSize = str2double(ui.editPixelsize.String);
        
        CreateData()
    end

    function EditLimitsCB(src,~)
        %Executed whenever axis limits are changed or the "Auto adjust" or
        %"Logarithmic" checkboxes are pressed
        
        if strcmp(src.Tag,'x')
            %User entered an axis limit so uncheck the "Auto adjust"
            %checkbox for the x-axis
            ui.cboxAutoX.Value = 0;
        elseif strcmp(src.Tag,'y')
            %User entered an axis limit so uncheck the "Auto adjust"
            %checkbox for the y-axis
            ui.cboxAutoY.Value = 0;
        end
        
        AdjustAxisLimits()
    end

    function CopyToWorkspaceCB(src,~)
        %Executed when "Export to Matlab workspace" is pressed
        
        src.BackgroundColor = 'r';
        drawnow
        
        assignin('base','allAnalysisResults',results);
        assignin('base','valuesInCurrentPlot',currentPlotValues);
        
        
        src.BackgroundColor = [.94 .94 .94];
        
    end

    function CloseHistogram()
        delete(gcf)
    end

    function PointsToFitBtnCB(src,~)
        %Executed when user presses buttons in the msd analysis panel
        
        switch src.Tag
            case 'PointsToFit'
                ui.editPointsToFit.String = src.String;
            case 'Offset'
                ui.editOffset.String = src.String;
            case 'FitFun'
            case 'Alpha'
                ui.editAlphaThres.String = src.String;
                TlSelectionCB()
        end

    end


end

