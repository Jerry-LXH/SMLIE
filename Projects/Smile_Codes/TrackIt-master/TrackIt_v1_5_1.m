%% Copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TrackIt

% Copyright (C) 2020 Timo Kuhn, Johannes Hettich and J. Christof M. Gebhardt
% timo.kuhn@uni-ulm.de, johannes.hettich@uni-ulm.de, christof.gebhardt@uni-ulm.de

% https://gitlab.com/GebhardtLab/TrackIt

% Publication:
% Timo Kuhn, Johannes Hettich, Rubina Davtyan, J. Christof M. Gebhardt
% Single molecule tracking and analysis framework including theory-predicted parameter settings
% Scientific Reports 11, 9465 (2021). doi: https://doi.org/10.1038/s41598-021-88802-7

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% TrackIt Code
function TrackIt_v1_5_1()

%Add all subfolders to the matlab path
mainFolder = fileparts(which(mfilename));
addpath(genpath(mainFolder));

%Position and size of the figure (Hor. pos, vertical pos, hor. size, vert. size)
figurePosition = [0 .06 .75 .85];

%Initialize pixel size (µm/px)´
pixelSize = 0;

%Set a starting path for movies here e.g. searchPath{1} = 'D:\Data\CDX2'; 
searchPath{1} = pwd;

%Set a starting path for batch files e.g. searchPath{2} = 'D:\Data\CDX2';
searchPath{2} = searchPath{1};

%Filepath of the current batch file
searchPath{3} = '';

%Initialize plot settings
subROIColors = distinguishable_colors(10,{'k'});
plotSettings = struct(... 
    'colMap',               'gray',...  %LuT
    'invertColMap',         false,...   %Invert pixel colormap
    'bgColor',              'black',... %LuT pixel color for 0 value
    'trackLinewidth',       1.4,...     %Plotted track linewidth
    'trackMarkerSize',      8,...       %Size of the plotted spot points in a track
    'initialPosMarker',     '.',...     %Marker for initial position
    'roiColor',             'w',...     %Color of the tracking region
    'roiLinewidth',         1.1,...     %Tracking region linewidth
    'spotColor',            'b',...     %Color of the plotted spots
    'spotMarkerSize',       14,...      %Size of the plotted spots
    'spotMarker',           'o',...     %Marker of the plotted spots
    'singleColor',          'y',...     %Color in which non-linked spots are plotted
    'singleMarkerSize',     7,...       %Size in which events < min. track length are plotted
    'singleMarker',         '.',...     %Marker of events < min. track length
    'frameRate',            20,...      %Rate at which frames are shown if the "Play" button is pressed
    'sbDim',                [0 0],...   %Scale bar dimension in pixels
    'sbPos',                [10 10],... %Position of the lower left corner of the scale bar in pixels
    'sbTextString',         '',...      %Scale bar text
    'sbTextPos',            [10 11],... %Position of the scale bar text
    'sbTextSize',           15,...      %Size of the scale bar text
    'scalingFactor',        2,...       %Scaling factor at which "Detection mapping" or "Jump distance mapping" is plotted
    'nFramesTrackIsVisible',0,...       %Number of frames a track is still visible after it disappeared
    'ITM',                  false,...   %Wether to use an interlaced timelapse scheme to distinguish long and short bound molecules
    'subRoiColors',         subROIColors(2:end,:),...%Color order of the sub-regions
    'curTrackSel',          0,...       %Track Id of the currently selected track
    'tracksInFrame',        {[]},...    %Matrix of logicals indicating which track has to be shown in each frame (will be set later)
    'trackColors',          [],...      %List of the colors for each track (will be set later));
    'rectZoomWinSize',      50);        %Window size that is used to zoom into a rectangular region 

%-------------Do not change anything below here----------------------------------------------

%Index of currently shown movie
curMovieIndex = 1;

%Timelapse dependent tracking parameters
tlDependentTrackingParams = zeros(0,5);

%Initialize batch structure and table of movie filenames
[batch,filesTable] = init_batch();

%Stacks for: 1: tracking, 2: original second channel, 3: filtered second
%channel ,4: Z-Projection of tracking stack or TALM/Jump distance mapping
movieStacks = cell(4,1);

%Timer for movie playback
timerObj = timer('ExecutionMode', 'FixedRate', ...
    'Period', 1/20, ...
    'TimerFcn', {@Playing});

%Initialize user interface
ui = InitUI();

%Initialize plot

imageHandle = imshow(1,'Parent',ui.axes1);
hold on;
ROIHandle = plot(ui.axes1,NaN, NaN,'LineStyle','--');
subROIHandle = gobjects(0);
trackHandles = gobjects(0);
trackInitialPosHandle = scatter(ui.axes1,NaN,NaN, 'MarkerFaceAlpha',1,'MarkerEdgeColor','flat','MarkerFaceColor','none');
spotsHandle = plot(ui.axes1,NaN,NaN,'LineStyle','none');
spotsNonLinkedHandle = plot(ui.axes1,NaN,NaN,'LineStyle','none');
scalebarHandle = fill(ui.axes1,NaN, NaN, 'w');
scalebarTextHandle = text(ui.axes1,0,0,'','Color','w','FontSize',15,'VerticalAlignment','top');
hold off

%Initialize colormap to the one specified in plotSettings.colMap
initColMap = find(strcmp(plotSettings.colMap, ui.popLut.String));
if ~isempty(initColMap)
    ui.popLut.Value = initColMap;
else
    ui.popLut.Value = 1;
end
PlotSettingsChangedCB(ui.popLut)

%Initialize user interface

    function ui = InitUI()
        
        %% Set menu and toolbar
        moviePanelWidth = .125;
        ui.hFig = figure('Units','normalized','Position',figurePosition,'MenuBar','None','toolBar','figure','Name',mfilename,'NumberTitle','off');
        ui.axes1 = axes('Units','normalized','Position',[.135 .1 .73 .89],'XTickMode','manual','YTickMode','manual');
        
        
        ui.menuFile = uimenu(ui.hFig,'Label','File');
        uimenu(ui.menuFile,'Label','Load batch file','Callback',@LoadBatchFileCB);
        uimenu(ui.menuFile,'Label','Save batch file as...','Tag','saveAs','Callback',@SaveBatchCB);
        uimenu(ui.menuFile,'Label','Merge multiple batch files','Callback',@MergeBatchFilesCB);
        uimenu(ui.menuFile,'Label','Subdivide batch file','Callback',@SplitBatchCB);
        uimenu(ui.menuFile,'Label','Re-analyze multiple batch files','Callback',@MultiBatchAnalyzerCB);
        uimenu(ui.menuFile,'Label','Export all data to Matlab workspace','Separator','on','Callback',@CopyToWorkspaceCB);
        uimenu(ui.menuFile,'Label','Export tracks to .mat or .csv','Callback',@(~,~)ExportTracksCB);
        uimenu(ui.menuFile,'Label','Create .avi movie','Callback',@CreateMovieCB);
        
        ui.menuTools = uimenu(ui.hFig,'Label','Tools');
        
        uimenu(ui.menuTools,'Label','Track explorer','Tag','menucall','Callback',@TrackExplorerCB);
        uimenu(ui.menuTools,'Label','Kymograph','Callback',@KymographToolCB);
        uimenu(ui.menuTools,'Label','Movie splitter','Callback',@SplitMoviesCB);
        uimenu(ui.menuTools,'Label','Predict tracking radii (GRID)','Callback',@PredictParamsCB);      
        uimenu(ui.menuTools,'Label','Find stack 2 for all movies','Separator','on','Callback',@FindAllStack2CB);        
        uimenu(ui.menuTools,'Label','Remove ROI from all movies','Callback',@RemoveRoiFromAllMoviesCB);  
        uimenu(ui.menuTools,'Label','Reload all ROIs from ROI Folder','Callback',@ReloadRoisCB);
        uimenu(ui.menuTools,'Label','Change filenames in current batch','Callback',@ChangeFilenamesCB);
        
        ui.menuAnalysis = uimenu(ui.hFig,'Label','Analysis');
        
        uimenu(ui.menuAnalysis,'Label','Tracking data analysis','Callback',@TrackingDataAnalysisCB);
        uimenu(ui.menuAnalysis,'Label','Spots statistics','Callback',@SpotStatisticsCB);
        uimenu(ui.menuAnalysis,'Label','Analyse dissociation rates (GRID)','Tag','GRID','Callback',@StartGridCB);
        
        
        ui.menuHelp = uimenu(ui.hFig,'Label','Help');
        uimenu(ui.menuHelp,'Label','Open manual','Callback',@OpenManualCB);
        uimenu(ui.menuHelp,'Label','About','Callback',@AboutCB);
        
                
        % Delete unnecessary toolbar icons
        if ~verLessThan('matlab','9.4')
            addToolbarExplorationButtons(ui.hFig)
        end
        
        hTools1 = findall(ui.hFig,'Type','uipushtool');
        hTools2 = findall(ui.hFig,'Type','uitogglesplittool');
        hTools3 = findall(ui.hFig,'Type','uitoggletool');
        
        
        delete(hTools1)
        delete(hTools2)
        
        for c = 1:length(hTools3)
            if ~strcmp(hTools3(c).Tag, 'Exploration.DataCursor') && ~strcmp(hTools3(c).Tag, 'Exploration.ZoomIn') &&...
                    ~strcmp(hTools3(c).Tag, 'Exploration.ZoomOut') && ~strcmp(hTools3(c).Tag, 'Exploration.Pan')
                delete(hTools3(c))
            end
        end
        
        % Read ellipse image
        [imgEll,mapEll] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','tool_ellipse.gif'));
        
        [imgRect,mapRect] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','tool_rectangle.gif'));
        
        %Convert image from indexed to truecolor
        roiIcon = ind2rgb(imgEll,mapEll);
        subRoiIcon = ind2rgb(imgEll,mapEll);
        squareIcon = ind2rgb(imgRect,mapRect);
        
        %Create cross in rectangle
        imgRect(4:13,8:9) = 0;
        imgRect(8:9,4:13) = 0;        
        squarePlusIcon = ind2rgb(imgRect,mapRect);
        
        
        %Convert to blue
        firstPlane = subRoiIcon(:,:,3);
        firstPlane(firstPlane < 1) = 1;
        subRoiIcon(:,:,3) = firstPlane;
        
        %Create toolbar icons
        defaultToolbar = findall(ui.hFig,'Type','uitoolbar');
        ui.scaleToRoiTool = uitoggletool(defaultToolbar,'CData',roiIcon,'TooltipString','Scale to ROI','Tag','ScaleToRoi','ClickedCallback',@ScaleToRoiCB);
        ui.scaleToSubRoiTool = uitoggletool(defaultToolbar,'CData',subRoiIcon,'TooltipString','Scale to subROI','Tag','ScaleToSubRoi','ClickedCallback',@ScaleToRoiCB);
        ui.scaleToSquare = uitoggletool(defaultToolbar,'CData',squareIcon,'TooltipString','Scale to square region','Tag','ScaleToSquare','ClickedCallback',@ScaleToRoiCB);
        ui.zoomToSquare = uipushtool(defaultToolbar,'CData',squarePlusIcon,'TooltipString','Zoom to square region with specific window size','Tag','ZoomToSquare','ClickedCallback',@ScaleToRoiCB);
        
        % set(groot,'DefaultUIControlFontSize',9)
        
        %% Movie panel
        
        MoviePanelHeight = 0.2;
        btnSize  = [.94  .19];
        btnSize2 = [.45  .19];
        txtSize = [.2   .15];
        panelMovie                      = uipanel('Title','Movie','Position',[.005 .8 moviePanelWidth MoviePanelHeight]);
        ui.btnSelectMovies              = uicontrol(panelMovie,'Units','normalized','Position', [.025 .8  btnSize],'String','Movie selector','Callback',@MovieSelectorCB);
        ui.btnNextMovie                 = uicontrol(panelMovie,'Units','normalized','Position', [.52  .6   btnSize2],'String','Next','Callback',@MovieNumberCB);
        ui.btnPreviousMovie           	= uicontrol(panelMovie,'Units','normalized','Position', [.025 .6   btnSize2],'String','Previous','Callback',@MovieNumberCB);
        ui.textMovie                    = uicontrol(panelMovie,'Units','normalized','Position', [.1   .38   txtSize],'Style','Text','String','Movie','HorizontalAlignment','Left');
        ui.editMovie                    = uicontrol(panelMovie,'Units','normalized','Position', [.52  .41 txtSize],'Style','Edit','String','0','HorizontalAlignment','Right','Callback',@MovieNumberCB);
        ui.textMovie2                   = uicontrol(panelMovie,'Units','normalized','Position', [.75  .38   txtSize],'Style','Text','String','/0','HorizontalAlignment','Left');
        ui.btnRemoveMovie               = uicontrol(panelMovie,'Units','normalized','Position', [.025 .03  btnSize],'String','Remove current movie','Callback',@RemoveMovieCB);
        ui.popMovieSelection            = uicontrol(panelMovie,'Units','normalized','Position', [.025 .24  btnSize(1) .15],'Style','popupmenu','String',{''},'Tag','popMenu','Callback',@MovieNumberCB);
        
        %% Second Stack panel
        ui.panelSecondStack = uipanel(ui.hFig,'Title','Second movie / image','Position',[.005 .55 moviePanelWidth .25],'Visible','on');
        
        cboxsize = [0.8 0.10];
        btnsize = [0.90 0.15];
        editsize = [0.4 0.08];
        
        ui.cboxShowStack2              = uicontrol(ui.panelSecondStack,'Units','normalized','Position',[.05 .87 cboxsize],'Style','checkbox','String','Show stack 2','Callback',@SelectStack2CB);
        ui.textReplaceString1           = uicontrol(ui.panelSecondStack,'Units','normalized','Position',[.05 .77 editsize],'Style','text','String','Find');
        ui.textReplaceString2           = uicontrol(ui.panelSecondStack,'Units','normalized','Position',[.55 .77 editsize],'Style','text','String','Replace with');
        ui.editReplaceString1           = uicontrol(ui.panelSecondStack,'Units','normalized','Position',[.05 .7 editsize],'Style','edit','String','');
        ui.editReplaceString2           = uicontrol(ui.panelSecondStack,'Units','normalized','Position',[.55 .7 editsize],'Style','edit','String','');
        ui.btnLoadSecondStack           = uicontrol(ui.panelSecondStack,'Units','normalized','Position',[.05 .52 btnsize],'String','Load 2nd movie / image','Callback',@SelectStack2CB);
        
        ui.tableFilterOptions = uitable(ui.panelSecondStack,'Units','normalized',...
            'Position',[0.01 0.01 .985 .45],...
            'ColumnName',{},...
            'ColumnWidth',{20,125,25},...
            'RowName',{},...
            'Data',{false,'Flip horizontally','';false,'Average frames',2;false,'Moving frame average',2;false,'Gaussian filter (px)', 2;false, 'Remove background (px)', 20},...
            'ColumnEditable',[true,false,true],...
            'CellEditCallback',@ProcessStack2CB);
        
        %% Roi panel
        ui.txtROI                       = uicontrol(ui.hFig,'Units','normalized','Position',[.005 .49 moviePanelWidth .05],'Style','text','String','Regions of interest','FontSize',9);

        ui.popDrawMode                  = uicontrol(ui.hFig,'Units','normalized','Position',[.005 .47 moviePanelWidth .05],'Style','popupmenu','String',{'Draw freehand','Assisted freehand','Draw polygon','Draw ellipse'});
        ui.btnLoadRoiFromTxt            = uicontrol(ui.hFig,'Units','normalized','Position',[.005 .46 moviePanelWidth .03],'String','Load regions from .roi file','Tag','MainROI','Callback',@ImportRoiFomFileCB);

        if verLessThan('matlab','9.5')
            ui.popDrawMode.Visible = 'off';
        end
        
        btnsize = [0.95 0.8];
        ui.panelRoi = uipanel(ui.hFig,'Position',[.005 .395 moviePanelWidth .06],'Title', 'Tracking-region');%,'BorderType','none'
        ui.btnDrawMainROI               = uicontrol(ui.panelRoi,'Units','normalized','Position',[.025 .1 btnsize(1)/2-.01 btnsize(2)],'String','Draw ROI','Tag','MainROI','Callback',@AddMainROICB);
        ui.btnDeleteMainROI             = uicontrol(ui.panelRoi,'Units','normalized','Position',[.52 .1 btnsize(1)/2-.01 btnsize(2)],'String','Delete ROI','Tag','MainROI','Callback',@RemoveMainROICB);
        
        %% subRoi panel
        btnsize = [0.95 0.15];
        
        ui.panelSubROI = uipanel(ui.hFig,'Position',[.005 .165 moviePanelWidth .23],'Title', 'Sub-regions','Visible','on');%,'BorderType','none'
        ui.popBorderHandling              = uicontrol(ui.panelSubROI,'Units','normalized','Position',[.025 .81 .95 btnsize(2)],'Style','popupmenu','String',{'Assign by first appearance', 'Split tracks at borders', 'Delete tracks crossing borders', 'Only use tracks crossing borders'});
        ui.btnDrawSubROI                = uicontrol(ui.panelSubROI,'Units','normalized','Position',[.025 .64 btnsize(1)/2-.01 btnsize(2)],'String','Hand-drawn','Tag','SubROI','Callback',@DrawSubROICB);
        ui.btnDrawSubROI                = uicontrol(ui.panelSubROI,'Units','normalized','Position',[.52 .64 btnsize(1)/2-.01 btnsize(2)],'String','Threshold','Callback',@DrawBrightnessROICB);
        ui.btnDeleteSubROI              = uicontrol(ui.panelSubROI,'Units','normalized','Position',[.025 .48 btnsize(1)/2-.01 btnsize(2)],'String','Delete selection','Tag','SubROI','Callback',@DeleteSubROICB);
        ui.btnMergeSubROI               = uicontrol(ui.panelSubROI,'Units','normalized','Position',[.52 .48 btnsize(1)/2-.01 btnsize(2)],'String','Merge selection','Tag','SubROI','Callback',@MergeSubROICB);
        ui.listSubRoi                   = uicontrol(ui.panelSubROI,'Units','normalized','Position',[.025 .025 .95 .4],'Max',2,'Style','Listbox');
        
        %% Tracking panel
        ui.tabgpTracking = uipanel(ui.hFig,'Title','Tracking','Position',[.874 .55 moviePanelWidth .45],'Visible','on');
        
        panelFindSpots                  = uipanel(ui.tabgpTracking,'Title','Spot detection parameters','Position',[.028 .73 .96 0.26]);
        editsize = [.25 .23];
        txtsize = [.21 .35];
        ui.cboxFindSpots                = uicontrol(panelFindSpots,'Units','normalized','Position',[.05 .7 .9 0.25],'Style','Checkbox','String','Find spots','Value',1,'Callback',@cboxFindSpotsCB);
        ui.textThresFactor              = uicontrol(panelFindSpots,'Units','normalized','Position',[.05 .25 .6 txtsize(2)],'Style','Text','String','Threshold factor','HorizontalAlignment','Left');
        ui.editThresFactor              = uicontrol(panelFindSpots,'Units','normalized','Position',[.70 .40 editsize],'Style','edit', 'String',2);
        ui.textFramerange1              = uicontrol(panelFindSpots,'Units','normalized','Position',[.05 .0 txtsize],'Style','Text','String','Framerange','HorizontalAlignment','Left');
        ui.editFramerange1              = uicontrol(panelFindSpots,'Units','normalized','Position',[.35 .05 editsize],'Style','Edit','String',1);
        ui.textFramerange2              = uicontrol(panelFindSpots,'Units','normalized','Position',[.63 -.08 txtsize],'Style','Text','String','-','HorizontalAlignment','Left');
        ui.editFramerange2              = uicontrol(panelFindSpots,'Units','normalized','Position',[.70 .05 editsize],'Style','Edit','String',inf);
        
        editsize = [0.25 0.11];
        txtsize = [0.6 0.16];
        panelTrackingParameters         = uipanel(ui.tabgpTracking,'Title','Tracking parameters','Position',[.023 .2 .96 .52]);
        ui.popTrackingMethod            = uicontrol(panelTrackingParameters,'Units','normalized','Position', [.05 .95 .92 .02],'Style','popupmenu','String',{'Nearest neighbour','u-track random motion','u-track linear+random motion'});
        ui.cboxTlDependentTr            = uicontrol(panelTrackingParameters,'Units','normalized','Position', [.05 .68 0.99 0.15],'Style','checkbox','String','Timelapse specific','Callback',@(~,~)TrackingParamsCB);
        ui.btnTlDependentTr             = uicontrol(panelTrackingParameters,'Units','normalized','Position', [.78 .69 0.18 0.11],'String','edit','Callback',@EditTrackingParamsCB);
        ui.textTrackingRadius           = uicontrol(panelTrackingParameters,'Units','normalized','Position', [.05 .47 txtsize],'Style','Text','String','Tracking radius (px)','HorizontalAlignment','Left');
        
        ui.editTrackingRadius           = uicontrol(panelTrackingParameters,'Units','normalized','Position', [.70 .53 editsize],'Style','Edit','String',1);
        cmenu=uicontextmenu;
        
        uimenu(cmenu,'label','Convert from µm', 'Callback', @ConvertTrackingRadiusCB);
        set(ui.editTrackingRadius,'uicontextmenu',cmenu);
        ui.textMinTrackLength           = uicontrol(panelTrackingParameters,'Units','normalized','Position', [.05 .31 txtsize],'Style','Text','String','Min. track length','HorizontalAlignment','Left');
        ui.editMinTrackLength           = uicontrol(panelTrackingParameters,'Units','normalized','Position', [.70 .37 editsize],'Style','Edit','String',2);
        ui.textGapFrames                = uicontrol(panelTrackingParameters,'Units','normalized','Position', [.05 .15 txtsize],'Style','Text','String','Gap frames','HorizontalAlignment','Left');
        ui.editGapFrames                = uicontrol(panelTrackingParameters,'Units','normalized','Position', [.70 .21 editsize],'Style','Edit','String',1);
        ui.textMinLengthBeforeGap       = uicontrol(panelTrackingParameters,'Units','normalized','Position', [.05 .02 txtsize],'Style','Text','String','Min. track length before gap frame','HorizontalAlignment','Left');
        ui.editMinLengthBeforeGap       = uicontrol(panelTrackingParameters,'Units','normalized','Position', [.70 .05 editsize],'Style','Edit','String',0);
        
        ui.btnAnalyzeThisMovie          = uicontrol(ui.tabgpTracking,'Units','normalized','Position', [.05 .1 .90 .075],'String','Analyze current movie','Tag','Current','Callback',@StartTrackingCB);
        ui.btnStartBatch                = uicontrol(ui.tabgpTracking,'Units','normalized','Position', [.05 .01 .90 .075],'String','Analyze all movies','Tag','All','Callback',@StartTrackingCB);
        
        %% Plotting panel
        
        ui.tabgpPlot = uipanel(ui.hFig,'Title','Plot properties','Position',[.874 .17 moviePanelWidth .37],'Visible','on');
        ui.cboxSpots                    = uicontrol(ui.tabgpPlot,'Units','normalized','Position',[.05 .9 .9 .07],'Style','checkbox','Value',1,'String','Show Spots','Callback',@(~,~)FrameSlider);
        ui.cboxTracks                   = uicontrol(ui.tabgpPlot,'Units','normalized','Position',[.05 .815 .1 .07],'Style','checkbox','Value',1,'String','','Callback',@(~,~)FrameSlider);
        ui.popTracks                    = uicontrol(ui.tabgpPlot,'Units','normalized','Position',[.15 .8 .8 .08],'Style','popupmenu','String',{'Show tracks in Range','Show all Tracks','Show initial positions','Show all tracked positions'},'Tag','nFramesTrackIsVisible','Callback',@PlotSettingsChangedCB);
        ui.textTrackLength              = uicontrol(ui.tabgpPlot,'Units','normalized','Position',[.05 .665 .7 .1],'Style','Text','String','#frames track is visible','HorizontalAlignment','Left');
        ui.editTrackLength              = uicontrol(ui.tabgpPlot,'Units','normalized','Position',[.75 .71 .2 .065],'Style','Edit','String',num2str(plotSettings.nFramesTrackIsVisible),'Tag','nFramesTrackIsVisible','Callback',@PlotSettingsChangedCB);
        ui.cboxShowSingle               = uicontrol(ui.tabgpPlot,'Units','normalized','Position',[.05 .6 .9 .07],'Style','checkbox','Value',0,'String','Show non-linked spots','Callback',@PlotSettingsChangedCB);
        ui.popColoredTrackLengths       = uicontrol(ui.tabgpPlot,'Units','normalized','Position',[.05 .5 .9 .07],'Style','popupmenu','String',{'Random colored tracks','Colored by frame of appearance','Colored by track length','Colored by track length regime','Colored by mean jump distance','Colored by mean jump distance regime','Colored by sub-region'},'Tag','coloredTrackLengths','Callback',@PlotSettingsChangedCB);
        ui.textNTrackColors             = uicontrol(ui.tabgpPlot,'Units','normalized','Position',[.05 .4 .7 .065],'Style','Text','String','#track-colors','HorizontalAlignment','Left');
        ui.editNTrackColors             = uicontrol(ui.tabgpPlot,'Units','normalized','Position',[.75 .4 .2 .065],'Style','Edit','String',7,'Tag','nTrackColors','Callback',@PlotSettingsChangedCB);
        
        editsize = [0.25 0.2];
        txtsize = [0.55 0.3];
        cboxSize = [0.9 0.15];
        
        ui.panelTrackLengthRegimes         = uipanel(ui.tabgpPlot,'Title','','Position',[0 .18 1 .3],'BorderType','none','Visible','off');
        ui.cboxShowShort                = uicontrol(ui.panelTrackLengthRegimes,'Units','normalized','Position',[.05 .85 cboxSize],'Style','checkbox','Value',1,'String','Show short tracks','Tag','coloredTrackLengths','Callback',@PlotSettingsChangedCB);
        ui.cboxShowLong                 = uicontrol(ui.panelTrackLengthRegimes,'Units','normalized','Position',[.05 .58 cboxSize],'Style','checkbox','Value',1,'String','Show long tracks','Tag','coloredTrackLengths','Callback',@PlotSettingsChangedCB);
        ui.textNDarkForLong             = uicontrol(ui.panelTrackLengthRegimes,'Units','normalized','Position', [.05 .15 txtsize],'Style','Text','String','Min. length to count as long track','HorizontalAlignment','Left');
        ui.editNDarkForLong             = uicontrol(ui.panelTrackLengthRegimes,'Units','normalized','Position', [.70 .2 editsize],'Style','Edit','String',3,'Tag','coloredTrackLengths','Callback',@PlotSettingsChangedCB);
        
        
        ui.panelJumpDistRegimes = uipanel(ui.tabgpPlot,'Title','','Position',[0 .18 1 .3],'BorderType','none','Visible','off');
        ui.cboxShowLowDist              = uicontrol(ui.panelJumpDistRegimes,'Units','normalized','Position',[.05 .85 cboxSize],'Style','checkbox','Value',1,'String','Show immobile tracks','Tag','coloredTrackLengths','Callback',@PlotSettingsChangedCB);
        ui.cboxShowHighDist             = uicontrol(ui.panelJumpDistRegimes,'Units','normalized','Position',[.05 .58 cboxSize],'Style','checkbox','Value',1,'String','Show mobile tracks','Tag','coloredTrackLengths','Callback',@PlotSettingsChangedCB);
        ui.textDistForHigh              = uicontrol(ui.panelJumpDistRegimes,'Units','normalized','Position', [.05 .15 txtsize],'Style','Text','String','Distance cut-off for mobile regime','HorizontalAlignment','Left');
        ui.editDistForHigh              = uicontrol(ui.panelJumpDistRegimes,'Units','normalized','Position', [.70 .2 editsize],'Style','Edit','String',1.5,'Tag','coloredTrackLengths','Callback',@PlotSettingsChangedCB);
        
        
        ui.btnPlotProperties            = uicontrol(ui.tabgpPlot,'Units','normalized','Position', [.05 .01 .40 .08],'String','Advanced','Tag','advanced','Callback',@PlotSettingsChangedCB);
        ui.btnScalebar            = uicontrol(ui.tabgpPlot,'Units','normalized','Position', [.55 .01 .40 .08],'String','Scale bar','Tag','scalebar','Callback',@PlotSettingsChangedCB);
        ui.btnPlotSelectedTrack         = uicontrol(ui.tabgpPlot,'Units','normalized','Position', [.05 .1 .90 .08],'String','Plot selected track','Tag','trackselected','Callback',@TrackExplorerCB);
        
        %% Lower part of the GUI
        ui.editFeedbackWin              = uicontrol('Units','normalized','Position',[.005 .01 moviePanelWidth .15],'Style','Text','BackgroundColor',[1 1 1],'HorizontalAlignment','left','Callback',@UpdateSettingsCB);
                
        ui.btnPlay                      = uicontrol('Units','normalized','Position',[.20 .01 .035 .02],'String','Play','Callback',@PlayCB);
        ui.sliderFrame                  = uicontrol('Units','normalized','Position',[.25 .01 .56 .02],'Style','slider','Min',1,'Max',1,'Value',1,'SliderStep',[1 1]);
        addlistener(ui.sliderFrame,'Value','PostSet',@(~,~)FrameSlider);
        
        % %Checbox for showing spots and tracks
        
        ui.textFPS                      = uicontrol('Units','normalized','Position',[.135  .01 .015 .02],'Style','text','String','FPS','HorizontalAlignment','left');
        ui.editFPS                      = uicontrol('Units','normalized','Position',[.155  .01 .03 .02],'Style','edit','String',num2str(plotSettings.frameRate));
        
        ui.textLut                      = uicontrol('Units','normalized','Position',[.874 .136 .08 .02],'Style','text','String','LUT','Visible','on','HorizontalAlignment','left');
        ui.popLut                       = uicontrol('Units','normalized','Position',[.892 .142 .04 .02],'Style','popupmenu','String',{'gray','jet','parula','hot','inferno','magma','plasma','viridis'},'Tag','lut','Callback',@PlotSettingsChangedCB);
        
        ui.textScalingFactor            = uicontrol('Units','normalized','Position',[.94 .136 .04 .028],'Style','text','String','Scaling factor','Visible','off','HorizontalAlignment','left');
        ui.editScalingFactor            = uicontrol('Units','normalized','Position',[.973  .138 .025 .023],'Style','edit','String',num2str(plotSettings.scalingFactor),'Visible','off','Callback',@CreateZProjection);
        
        ui.cboxShowZProjection          = uicontrol('Units','normalized','Position',[.874 .11 .04 .018],'Style','checkbox','String','','Callback',@CreateZProjection);
        ui.popShowImage                 = uicontrol('Units','normalized','Position',[.888 .11 .11 .02],'Style','popupmenu','String',{'Standard deviation','Average intensity','Max. intensity','Wavelet filtered image',...
            'No image','Detection map incl. non-linked','Detection map w/o non-linked','Jump distance map'},'Callback',@CreateZProjection);
        
        ui.textFileName                 = uicontrol('Units','normalized','Position',[.135 .04 .70 .02],'Style','Text','Enable','Inactive','HorizontalAlignment','Left','String','Filename: ','ButtonDownFcn',@CopyFilenameCB);
        ui.textFrame                    = uicontrol('Units','normalized','Position',[.66 .04 .15 .02],'Style','Text','String','0/0','HorizontalAlignment','Right');
        
        ui.btnBrightness                = uicontrol('Units','normalized','Position',[.874 .07 .08 .03],'String','Autoadjust brightness','Callback',@(~,~) AdjustBrightnessCB('autoAdjust')); %<html>Auto<br>adjust
        ui.cboxContAutoAdj              = uicontrol('Units','normalized','Position',[.955 .07 .085 .03],'Style','checkbox','String','Always','Callback',@(~,~) AdjustBrightnessCB('autoAdjust')); %<html>Auto<br>adjust
        ui.sliderBlack                  = uicontrol('Units','normalized','Position',[.874 .04 .123 .02],'Style','slider','Min',0,'Max',1,'Value',0,'SliderStep',[1/100 1/12],'Tag','sliderBlack');
        ui.sliderWhite                  = uicontrol('Units','normalized','Position',[.874 .01 .123 .02],'Style','slider','Min',0,'Max',1,'Value',1,'SliderStep',[1/100 1/12],'Tag','sliderWhite');
        addlistener(ui.sliderBlack,'Value','PostSet',@(~,~)AdjustBrightnessCB('slider'));
        addlistener(ui.sliderWhite,'Value','PostSet',@(~,~)AdjustBrightnessCB('slider'));
        
        ui.editBlack                    = uicontrol('Units','normalized','Position',[.84 .04 .03 .018],'Style','edit','String','0','Tag','textInput','HorizontalAlignment','Right','Callback',@(~,~) AdjustBrightnessCB('textInput'));
        ui.editWhite                    = uicontrol('Units','normalized','Position',[.84 .01 .03 .018],'Style','edit','String','0','Tag','textInput','HorizontalAlignment','Right','Callback',@(~,~) AdjustBrightnessCB('textInput'));
        
        %% Other
        %Set mouse wheel handle for scrolling through frames
        set(ui.hFig, 'WindowScrollWheelFcn', {@MouseWheelCB});
        
        
        %Set callback for zoom
        ui.zoomeHandle = zoom;
        ui.zoomeHandle.ActionPostCallback = @ZoomCB;
        
        %Initialize datacursor
        ui.dataCoursor = datacursormode(ui.hFig);
        set(ui.dataCoursor,'UpdateFcn',{@DataCursorCB});
    end

%Movie Tab

    function MovieSelectorCB(src,~)
        
        %Open movie selector window
        [newFilesTable, canceled] = movie_selector(filesTable,searchPath{1});
        
        if canceled
            return
        end
        
        %Get number of files in movie selector
        nFiles = height(newFilesTable);
        
        if nFiles == 0
            %No files have been selected
            [batch,~] = init_batch();
            ui.popMovieSelection.String = {''};
        elseif nFiles
            %Initialize batch structure
            [batchInit,~] = init_batch();
            newBatch = repmat(batchInit, nFiles, 1);
            
            %Counter for movies that appear more than once
            nMultipleFiles = 0;
            
            %Go through all movies and either pass existing results to new
            %batch structure or treat as new movie
            for fileIdx = 1:nFiles
                %Check if movie was there before
                existingIdx = strcmp(newFilesTable.FileName{fileIdx}, filesTable.FileName(:));
                
                nExisting = sum(existingIdx);
                
                if nExisting > 1
                    nMultipleFiles = nMultipleFiles + 1;
                end
                
                if nExisting == 1
                    %Movie was there before so pass old results to new batch
                    newBatch(fileIdx) = batch(existingIdx);
                    %Update tl condition in case user changed it
                    newBatch(fileIdx).movieInfo.frameCycleTime = newFilesTable.frameCycleTime(fileIdx);
                else
                    %Movie is new or exists more than once
                    %Save filename, pathname and timelapse condition
                    newBatch(fileIdx).movieInfo.fileName = newFilesTable.FileName{fileIdx};
                    newBatch(fileIdx).movieInfo.pathName = newFilesTable.PathName{fileIdx};
                    newBatch(fileIdx).movieInfo.frameCycleTime = newFilesTable.frameCycleTime(fileIdx);
                    
                    % Check if ROI file exists
                    [~,name,~] = fileparts(newFilesTable.FileName{fileIdx});
                    roiFileName = fullfile(newFilesTable.PathName{fileIdx},'ROI',strcat(name, '.roi'));
                    txtFileName = fullfile(newFilesTable.PathName{fileIdx},'ROI',strcat(name, '_ROI.txt'));
                    
                    if exist(roiFileName, 'file') == 2 %Load ROI from .roi file
                        loadedFile = load(roiFileName,'-mat');
                        newBatch(fileIdx).ROI = loadedFile.ROI;
                        newBatch(fileIdx).subROI = loadedFile.subROI;
                    elseif exist(txtFileName, 'file') == 2 %Load ROI from .txt file
                        newBatch(fileIdx).ROI{1} = load(txtFileName);
                    end
                end
            end
            
            if nMultipleFiles > 0
                warndlg(strcat(num2str(nMultipleFiles), ' filenames have been found more than once.' ))
            end
            
            batch = newBatch;
            
            ui.popMovieSelection.String = newFilesTable.FileName;
        end
        
        filesTable = newFilesTable;
        
        %Set movienumber of current movie
        curMovieIndex = 1;
        ui.popMovieSelection.Value = 1;
        ui.editMovie.String = 1;
        ui.textMovie2.String = ['/' num2str(length(batch))];
        
        %Reset timelapse dependent tracking parameters
        ui.cboxTlDependentTr.Value = 0;
        TrackingParamsCB()
        
        src.BackgroundColor = 'r'; %Button turns red during computation
        drawnow
        
        %If current filename of movie has changed, adjust ui to new movie
        if ~strcmp(ui.textFileName.String(10:end), batch(curMovieIndex).movieInfo.fileName)
            AdjustUiToNewMovie()
        end
        
        src.BackgroundColor = [.94 .94 .94];  %Button turns red during computation
    end

    function MovieNumberCB(src,~)
        
        try
            %User either selected a movie in the popup menu, clicked "Previous"
            %or "Next" button or inserted a movie number in the edit text-field
            
            persistent movieIsLoading % Shared with all calls of MovieNumberCB. Makes sure that movie is loaded completely before loading the next movie
            
            % If there are no recent clicks, execute the single
            % click action. If there is a recent click, ignore new clicks.
            if isempty(movieIsLoading)
                movieIsLoading = 1;
                
                %Save previous movie number
                previousMovieNumber = curMovieIndex;
                
                %Get new movie number and make sure that it is not out of range
                if strcmp(src.Tag,'popMenu')
                    %User used the popup menu to select a file
                    curMovieIndex = ui.popMovieSelection.Value;
                else
                    if strcmp(src.String,'Previous') && previousMovieNumber > 1
                        %User pressed 'previous' button
                        curMovieIndex =  previousMovieNumber - 1;
                    elseif strcmp(src.String,'Next') && previousMovieNumber < length(batch)
                        %User pressed 'next' button
                        curMovieIndex = previousMovieNumber + 1;
                    elseif str2double(src.String) <= length(batch) && str2double(src.String) > 0
                        %User entered number into edit field
                        curMovieIndex = str2double(src.String);
                    end
                    %Update the selected movie in the popup menu list
                    ui.popMovieSelection.Value = curMovieIndex;
                    
                end
                
                %Update new movie number in the edit field
                ui.editMovie.String = num2str(curMovieIndex);
                
                %Turn selected ui boject red as long as movie is beeing loaded
                src.BackgroundColor = 'r';
                drawnow
                
                %Check if movie number was changed, if yes load movie and update UI
                if previousMovieNumber ~= curMovieIndex
                    AdjustUiToNewMovie()
                end
                
                %Reset color of the ui object
                src.BackgroundColor = [.94 .94 .94];
                
                %Reset persistent variable so next movie can be loaded as
                %soon as the user wants it.
                movieIsLoading = [];
                
            else
                %A new movie is currently beeing loaded and the user
                %selected a different movie before loading was finished. So
                %reset the selected movie and the movie number.
                ui.popMovieSelection.Value = curMovieIndex;
                ui.editMovie.String = num2str(curMovieIndex);
            end
            
        catch ex
            %Some error occured so reset persistent variable so next movie
            %can be loaded as soon as the user wants it.
            movieIsLoading = [];
        end
        
    end

    function RemoveMovieCB(~,~)
        %User clicked "Remove current movie" button
        
        if length(batch) > 1
            %Remove current movie from batch structure and table of files;
            batch(curMovieIndex) = [];
            filesTable(curMovieIndex,:) = [];
            
            %Update list of movies in the popup menu
            ui.popMovieSelection.String = filesTable.FileName;
        else
            %User wants to remove the last remaining movie in batch so
            %re-initialize batch structure and the movie popup menu list
            [batch,filesTable] = init_batch();
            ui.popMovieSelection.String = {''};
        end
        
        if length(batch) < curMovieIndex
            %User wants to remove the last element in the batch structure array
            %so current  movie index has to be adjusted
            curMovieIndex = length(batch);
            
            %Update selected movies in the popup menu
            ui.popMovieSelection.Value = curMovieIndex;
        end
        
        %Reset timelapse dependent tracking parameters because the amount
        %of different time-lapse conditions might have changed
        ui.cboxTlDependentTr.Value = 0;
        TrackingParamsCB()
        
        %Set movienumber of current movie
        ui.editMovie.String = num2str(curMovieIndex);
        
        %Update amount of movies
        ui.textMovie2.String = ['/' num2str(length(batch))];
        
        %New movie has to be shown so update ui
        AdjustUiToNewMovie()
    end

    function AdjustUiToNewMovie()
        %The movie has changed so this function loads the new movie and
        %updates all ui settings accordingly
        
        %Get filename of current movie
        fileName = batch(curMovieIndex).movieInfo.fileName;
        pathName = batch(curMovieIndex).movieInfo.pathName;
        fullFilePath = fullfile(pathName,fileName);
        
        %Delete all pixel data
        movieStacks = cell(4,1);
        
        %--------Display empty image and abort if necessary----------------
        
        if ~(exist(fullFilePath, 'file') == 2) && batch(curMovieIndex).movieInfo.height == 0
            %The file is not found and no information about the movie size
            %is present so clean up plot and ui
            ui.sliderFrame.Max = 1;
            ui.sliderFrame.SliderStep = [1/1 1/1];
            ui.textFileName.String = 'Filename:';
            imageHandle.CData = 1;
            ui.axes1.XLim = [0.5 1];
            ui.axes1.YLim = [0.5 1];
            ui.axes1.CLim = [0 1];
            return
        end
        
        %-----------Check if movie has to be loaded------------------------
        
        %If the movie dimensions are known and the user selected "hide
        %image", "detection mapping" or "jump distance mapping", then the pixel
        %information is not needed and the movie will not be loaded
        if ~(ui.cboxShowZProjection.Value && ui.popShowImage.Value > 4) || batch(curMovieIndex).movieInfo.height == 0
            
            searchPath{1} = batch(curMovieIndex).movieInfo.pathName;
            movieStacks{1} = load_stack(searchPath{1}, fileName, ui);
            
            if ~isempty(movieStacks{1})
                %Movie has been loaded
                [height, width, nFrames] = size(movieStacks{1});
                
                %Save movie & properties
                batch(curMovieIndex).movieInfo.height = height;
                batch(curMovieIndex).movieInfo.width = width;
                batch(curMovieIndex).movieInfo.frames = nFrames;
            else
                %Movie has not been loaded because file was not found but
                %has been loaded/analyzed before -> Hide image
                ui.cboxShowZProjection.Value = 1;
                
                if ui.popShowImage.Value < 5
                    ui.popShowImage.Value = 5;
                end
            end
        end
        
        %Get new number of frames from batch structure
        nFrames = batch(curMovieIndex).movieInfo.frames;
        
        %--------Display tracking parameters in ui-------------------------
        
        %Check if movie has already been analyzed. If yes, display the analysis
        %parameters in the user interface
        
        if ~isempty(batch(curMovieIndex).results.spotsAll)
            ui.editThresFactor.String = batch(curMovieIndex).params.thresholdFactor;
            ui.editFramerange1.String = batch(curMovieIndex).params.frameRange(1);
            
            if batch(curMovieIndex).params.frameRange(2) == nFrames
                ui.editFramerange2.String = Inf;
            else
                ui.editFramerange2.String = batch(curMovieIndex).params.frameRange(2);
            end
            
            ui.editTrackingRadius.String = batch(curMovieIndex).params.trackingRadius;
            ui.editMinTrackLength.String = batch(curMovieIndex).params.minTrackLength;
            ui.editGapFrames.String = batch(curMovieIndex).params.gapFrames;
            ui.editMinLengthBeforeGap.String = batch(curMovieIndex).params.minLengthBeforeGap;
            
%             %Sub-rois where introduced in TrackIt v1.3 so we first have
%             %to check wether this parameter exists
%             if isfield(batch(curMovieIndex).params,'subRoiBorderHandling')
%                 subRoiBorderHandling = batch(curMovieIndex).params.subRoiBorderHandling;
%                 ui.popBorderHandling.Value = find(strcmp(subRoiBorderHandling, ui.popBorderHandling.String));
%             end
        end
        
        %----------Load second stack---------------------------------------
        %Make sure second stack is shown/loaded if the checkbox is checked
        %or uncheck the checkbox if no information on a second movie exists
        LoadStack2()
        
        %----------Update FeedbackWin--------------------------------------
        feedbackWin = char(...
            ['Movie ', num2str(curMovieIndex), ' :'],...
            [num2str(batch(curMovieIndex).results.nSpots),' detected spots'],...
            [num2str(batch(curMovieIndex).results.nTracks),' tracks'],...
            [num2str(batch(curMovieIndex).results.nNonLinkedSpots),' non-linked spots'],...
            [num2str(round(batch(curMovieIndex).results.meanTrackLength,1)),' frames mean tracklength'],...
            [num2str(round(batch(curMovieIndex).results.trackedFraction*100,1)),'% of spots linked to tracks']);
        
        ui.editFeedbackWin.String = char('Plotting', feedbackWin);
        drawnow
        
        %--------ROI and sub-ROI related adjustments-----------------------
        
        %Set whole movie as ROI if none was drawn
        if isempty(batch(curMovieIndex).ROI )
            batch(curMovieIndex).ROI{1} = [.5 .5; .5 batch(curMovieIndex).movieInfo.height+.5;...
                batch(curMovieIndex).movieInfo.width+.5 batch(curMovieIndex).movieInfo.height+.5;...
                batch(curMovieIndex).movieInfo.width+.5 .5; .5 .5];
        end
        
        
        %--------Set frame slider limits and update plot-------------------
        %Set framenumber to 1 and call FrameSlider to update plot
        %Update Frame slider
        ui.sliderFrame.Max = nFrames;
        ui.sliderFrame.SliderStep = [min(1,1/(nFrames-1)) 1/min((nFrames-1),10)];
        
        %Create a logical matrix to know which track has to be displayed in each frame
        GetTracksInFrame();
        
        %Make sure a z-projection is created if the z projection checkbox is checked
        CreateZProjection('newMovie')
                
        if ui.sliderFrame.Value > 1
            %Set frame slider to first frame and call FrameSlider() through the listener callback
            ui.sliderFrame.Value = 1;
        else
            %Slider is already at first frame so just call FrameSlider()
            FrameSlider()
        end
        
        AdjustBrightnessCB('autoAdjust')
        
        
        %Update list of sub-regions and adjust axis to new dimensions
        RoiChanged(false)
        
        %Plotting is finished so remove "plotting" message from feedback window
        ui.editFeedbackWin.String = feedbackWin;
                
    end

%Tracking Tab

    function TrackingParamsCB()
        %User clicked the "Timelapse specific" checkbox"
        if ui.cboxTlDependentTr.Value
            ui.editTrackingRadius.Enable = 'off';
            ui.editMinTrackLength.Enable = 'off';
            ui.editGapFrames.Enable = 'off';
            ui.editMinLengthBeforeGap.Enable = 'off';
            EditTrackingParamsCB()
        else
            ui.editTrackingRadius.Enable = 'on';
            ui.editMinTrackLength.Enable = 'on';
            ui.editGapFrames.Enable = 'on';
            ui.editMinLengthBeforeGap.Enable = 'on';
        end
        
    end

    function EditTrackingParamsCB(~,~)
        %User clicked on the "edit" button for timelapse (tl) specific tracking parameters
        
        %Create a list of the timelapse conditions in the current batch
        tlUni = unique(filesTable.frameCycleTime);
        
        %Create a list where the user can enter the tl parameters
        newParamsList = zeros(numel(tlUni),5);
        
        %Fill first row with tl conditions
        newParamsList(:,1) = tlUni;
        
        %Iterate through tl conditions
        for tlIdx = 1:numel(tlUni)
            %Get current tl
            curTL = newParamsList(tlIdx,1);
            
            %Get previously set tl specific parameters
            prevParamsForCurTl = tlDependentTrackingParams(tlDependentTrackingParams(:,1) == curTL,:);
            
            %Check if parameters have been set before by the user
            if isempty(prevParamsForCurTl)
                %No parameters have been set befor so fill the list with the parameters in the main ui
                newParamsList(tlIdx,2) = str2double(ui.editTrackingRadius.String);
                newParamsList(tlIdx,3) = str2double(ui.editMinTrackLength.String);
                newParamsList(tlIdx,4) = str2double(ui.editGapFrames.String);
                newParamsList(tlIdx,5) = str2double(ui.editMinLengthBeforeGap.String);
            else
                %Fill the list with the previously set parameters
                newParamsList(tlIdx,:) = prevParamsForCurTl;
            end
        end
        
        %Create figure where user can manipulated the list of tl dependent
        %tracking parameter settings
        f = figure('Units','normalized',...
            'Position',[0.1 0.5 .33 .25],...
            'MenuBar','None',...
            'Name','Timelapse dependent tracking parameters',...
            'NumberTitle','off',...
            'CloseRequestFcn', @ClosedCB);
        tb = uitable('Units','normalized',...
            'Position',[0.05 0.05 .9 .9],...
            'ColumnName',{'TL condition (ms)';'Tracking radius';'Min. track length';'Gap Frames';'Min. length before gap frame'},...
            'Data', newParamsList,...
            'ColumnEditable',[false,true,true,true,true]);
        
        %Wait for user to close figure
        uiwait(f)
        
        function ClosedCB(~,~)
            %User has closed the figure so save the new parameters
            tlDependentTrackingParams = tb.Data;
            delete(f)
        end
        
    end

    function StartTrackingCB(src,~)

        %User either clicked "Analyze current movie" or "Analyze all movies" button
        
        
        %-------UI preparations-------------------------------------------
        
        %Delete current feedback window entry
        ui.editFeedbackWin.String = '';
        
        %Change button appearance while tracking
        src.BackgroundColor = 'r';
        src.String = 'Press Ctrl + X to stop';
        drawnow;
        
        %-------Get tracking parameters------------------------------------
        
        %Get settings for spot detection
        para.boolFindSpots = ui.cboxFindSpots.Value;
        para.thresholdFactor = str2double(ui.editThresFactor.String);
        para.frameRange = [str2double(ui.editFramerange1.String),...
            str2double(ui.editFramerange2.String)];
        
        %Get settings for tracking
        para.trackingMethod = ui.popTrackingMethod.String{ui.popTrackingMethod.Value};
        
        if ui.cboxTlDependentTr.Value
            %"Timelapse specific" checkbox is checked so we will use a
            %unique set of tracking parameters for each tl condition
            para.tlConditions = tlDependentTrackingParams(:,1);
            para.trackingRadius = tlDependentTrackingParams(:,2);
            para.minTrackLength = tlDependentTrackingParams(:,3);
            para.gapFrames = tlDependentTrackingParams(:,4);
            para.minLengthBeforeGap = tlDependentTrackingParams(:,5);
        else
            %Use same tracking parameters for all tl conditions defined in
            %the main ui. The tlConditions variable can be left empty.
            para.tlConditions = [];
            para.trackingRadius = str2double(ui.editTrackingRadius.String);
            para.minTrackLength = str2double(ui.editMinTrackLength.String);
            para.gapFrames = str2double(ui.editGapFrames.String);
            para.minLengthBeforeGap = str2double(ui.editMinLengthBeforeGap.String);
        end
        
        %Get sub-region handling option
        para.subRoiBorderHandling = ui.popBorderHandling.String{ui.popBorderHandling.Value};
        
        %Get trackIt version
        para.trackItVersion = mfilename;
        
        
        %---------Do tracking----------------------------------------------
        
        switch src.Tag
            %Check if user pressed "Analyze current movie" or "Analyze all movies"
            case 'Current'
                if ~isempty(movieStacks{1})
                    %Stack is already loaded in the main gui so use that one
                    movieStack = movieStacks{1};
                else
                    %Stack was not loaded so pass an empty variable
                    movieStack = [];
                end
                
                origButtonString = 'Analyze current movie';
                
                %Start main tracking routine and analyze only current movie
                [out, boolCancelled] = tracking_routine(batch(curMovieIndex), para, movieStack, ui);
                                
                if ~boolCancelled
                    batch(curMovieIndex) = out;
                end
                
            case 'All'
                origButtonString = 'Analyze all movies';
                
                %Start main tracking routine
                [out, boolCancelled] = tracking_routine(batch, para, [], ui);
                
                if ~boolCancelled
                    batch = out;
                end
        end
        
        %--------Tracking finished -> update ui----------------------------
        
        
        ui.editFeedbackWin.String = char('Plotting', ui.editFeedbackWin.String);
        drawnow
        
        %Make sure a z-projection is created if the z projection checkbox is checked
        CreateZProjection('startBatch')
        
        %Create a logical matrix to know which track has to be displayed in each frame
        GetTracksInFrame()
        
        %Update plot
        FrameSlider()
        
        %Plotting is finished so remove "plotting" notice from feedback window
        ui.editFeedbackWin.String = ui.editFeedbackWin.String(2:end,:);
        
        %Reset button color to gray
        src.BackgroundColor = [.94 .94 .94];
        
        %Reset button string
        src.String = origButtonString;
    end

%ROI tab

    function ImportRoiFomFileCB(~,~)
        %User pressed "Load ROI from .roi file" button
        
        [filename,pathName] = uigetfile({'*.roi;*.txt;',...
            '*.roi,*.txt'; '*.*',  'All Files (*.*)'},'Select ROI.txt or.roi file',searchPath{1});
        
        
        if filename == 0
            %User didn't choose a file
            return
        end
        
        searchPath{1} = pathName;
        
        %Get file extension
        [~,~,ext] = fileparts(filename);
        
        if strcmp(ext,'.roi')
            %Load ROI from .roi file
            loadedFile = load(fullfile(pathName,filename),'-mat');
            batch(curMovieIndex).ROI = loadedFile.ROI;
            batch(curMovieIndex).subROI = loadedFile.subROI;
        elseif strcmp(ext,'.txt')
            %Load ROI from .txt file
            batch(curMovieIndex).ROI{1} = load(fullfile(pathName,filename));
        end
        
        RoiChanged(false)
        FrameSlider()
    end

    function AddMainROICB(src,~)
        try
            %User pressed "Draw ROI" button
            
            %Make sure current character is not esc key
            set(gcf,'currentch',char(1))
            
            if src.BackgroundColor(1) == 1
                %User is already drawing and pressed the same button again
                return
            end
            
            %Initialize ROI
            curROI = [];
            
            %---------Create ROI--------------
            
            if verLessThan('matlab','9.5')
                %Use imfreehand for old Matlab versions
                ROIFreehandHandle = imfreehand();
                
                if ~isempty(ROIFreehandHandle)
                    %Drawing is finished so save ROI in List
                    curROI = getPosition(ROIFreehandHandle);
                end
            else
                src.BackgroundColor = 'r';
                src.String = '<html>Press enter<br> to Finish';
                drawnow
                if ui.popDrawMode.Value == 1
                    %Draw freehand
                    ROIFreehandHandle = drawfreehand('Color','w');
                elseif ui.popDrawMode.Value == 2
                    %Use drawassisted
                    ROIFreehandHandle = drawassisted('Color','w');
                elseif ui.popDrawMode.Value == 3
                    %Draw polygon
                    ROIFreehandHandle = drawpolygon('Color','w');
                elseif ui.popDrawMode.Value == 4
                    %Draw ellipse
                    ROIFreehandHandle = drawellipse('Color','w');
                end
                
                if ~isempty(ROIFreehandHandle.Position)
                    % Wait for the user to press enter
                    while true
                        w = waitforbuttonpress;
                        switch w
                            case 1 %keyboard press
                                key = get(gcf,'currentcharacter');
                                if key == 13  %User pressed enter key
                                    break
                                end
                        end
                    end
                    
                    %Get ROI coordinates
                    if ui.popDrawMode.Value == 4
                        binaryMask = createMask(ROIFreehandHandle);
                        boundariesBinaryMask = bwboundaries(binaryMask);
                        curROI(:,1) = boundariesBinaryMask{1}(:,2);
                        curROI(:,2) = boundariesBinaryMask{1}(:,1);
                    else
                        curROI = ROIFreehandHandle.Position;
                    end
                end
            end
            
            if ~isempty(curROI)
                %Delete handle to ROI
                delete(ROIFreehandHandle)
                
                %Connect end with starting point
                ROI = cat(1, curROI, curROI(1,:));
                
                %Restrict ROI to image size if user drew outside and
                %account for scaling factor
                height = batch(curMovieIndex).movieInfo.height;
                width = batch(curMovieIndex).movieInfo.width;
                scalingFactor = plotSettings.scalingFactor;
                
                ROI(ROI(:) < 0) = 0.5;
                ROI(ROI(:,2) > height*scalingFactor,2) = height*scalingFactor+0.5;
                ROI(ROI(:,1) > width*scalingFactor,1) =width*scalingFactor+0.5;
                ROI = ROI./plotSettings.scalingFactor;
                ROI = {ROI};
                
                %Save ROI in batch structure
                batch(curMovieIndex).ROI = ROI;
                
                %Save ROI to file
                RoiChanged(true)
                
                %Reset results
                [batchInit,~] = init_batch();
                batch(curMovieIndex).results = batchInit.results;
            end
            
            FrameSlider()
            
            src.BackgroundColor = [.94 .94 .94];
            src.String = 'Draw ROI';
            
        catch ex
            disp(ex)
            %Reset button appearance and update sub-region list in UI
            src.BackgroundColor = [.94 .94 .94];
            src.String = 'Hand-drawn';
        end
        
    end

    function RemoveMainROICB(~,~)
        %User pressed "Detlete ROI" button
        
        %Get filename of current movie
        [~,fileName,~] = fileparts(batch(curMovieIndex).movieInfo.fileName);
        
        if ~isempty(fileName)
            %Reset ROI
            batch(curMovieIndex).ROI = {[.5 .5; .5 batch(curMovieIndex).movieInfo.height+.5;...
                batch(curMovieIndex).movieInfo.width+.5 batch(curMovieIndex).movieInfo.height+.5;...
                batch(curMovieIndex).movieInfo.width+.5 .5; .5 .5]};
            
            
            %Save ROI to file
            RoiChanged(true)
            
            %Reset results
            [batchInit,~] = init_batch();
            batch(curMovieIndex).results = batchInit.results;
            
            FrameSlider()
        end
    end

    function SelectStack2CB(src,~)
        %User pressed "Show Stack 2" or "Load 2nd movie / image" button
        
        %Check if user has to select a file: User pressed "load second movie" button or
        %pressed the "Show stack 2" checkbox and no file was chosen yet
        if  strcmp(src.Style, 'pushbutton') || (ui.cboxShowStack2.Value && isempty(batch(curMovieIndex).movieInfo.fileName2))
            
            %Check if user entered an insert string
            if ~isempty(ui.editReplaceString2.String)
                
                %Current filename and pathname of tracking movie
                fileName = batch(curMovieIndex).movieInfo.fileName;
                pathName = batch(curMovieIndex).movieInfo.pathName;
                
                %Get insert string
                insert = ui.editReplaceString2.String;
                newFilename = '';
                
                if isempty(ui.editReplaceString1.String)
                    %No replacement string was entered so just search for
                    %the insert string in the current folder
                    
                    %First check for direct match
                    
                    %Create search pattern
                    searchPattern = [insert,'.tif*'];
                    
                    %Create list of files containing the insert string
                    fileList = dir(fullfile(pathName,searchPattern));
                    
                    if isempty(fileList)
                        %No direct match found check if string is part of
                        %and filename
                        
                        %Create search pattern
                        searchPattern = ['*',insert,'*'];
                        
                        %Create list of files containing the insert string
                        fileList = dir(fullfile(pathName,searchPattern));
                    end
                    
                    %Get resulting file names
                    fileNameList = {fileList.name};
                    
                    if ~isempty(fileNameList)
                        %Use first file that has been found
                        newFilename = fileNameList{1};
                    end
                else %Find string and replacement string entered
                    %String which should be replaced in tracking movie filename
                    replaceString1 = ui.editReplaceString1.String;
                    %Split filename at the position of replacement string
                    splittedFilename = strsplit(fileName, replaceString1);
                    
                    %Create new filename with the fileparts of the tracking
                    %movie and the insert string from the user interface
                    if length(splittedFilename)>1
                        newFilename = strcat(splittedFilename{1},insert,splittedFilename{2});
                    end
                    
                end
                
                if ~isempty(newFilename)
                    %Create full directory
                    openingPath = fullfile(pathName,newFilename);
                else
                    openingPath = searchPath{1};
                end
            else
                %No replacement string was entered so use the path of the
                %current tracking movie for opening the dialog
                openingPath = searchPath{1};
            end
            
            if isfile(openingPath) && ui.cboxShowStack2.Value
                %User pressed show stack 2 and a file was found by
                %replacing the "find" string with "replace" string
                fileName = newFilename;
            else
                %Open file dialog box
                [fileName,pathName] = uigetfile({'*.tiff;*.tif;',...
                    '*.tif,*.tiff'; '*.*',  'All Files (*.*)'},'Select second movie',openingPath);
                
                if fileName == 0
                    %User didn't choose a file
                    ui.cboxShowStack2.Value = 0;
                    return
                end
            end
            %Save file and path of second movie in batch structure
            batch(curMovieIndex).movieInfo.fileName2 = fileName;
            batch(curMovieIndex).movieInfo.pathName2 = pathName;
            
            
            %Reinitialize variable where the original second movie and the filtered second movie is stored
            movieStacks{2} = {};
            movieStacks{3} = {};
            
            %Check the "Show stack 2" checkbox
            ui.cboxShowStack2.Value = 1;
        end
        
        if ui.cboxShowStack2.Value
            LoadStack2()
        else
            %Stack 2 is not shown (anymore) so display the filename of the tracking movie
            ui.textFileName.String = strcat('Filename:', batch(curMovieIndex).movieInfo.fileName);
        end
        
        CreateZProjection('stack2')
        FrameSlider()
        AdjustBrightnessCB('autoAdjust')
    end

    function LoadStack2()
        
        %Get file name of the second movie
        fileName = batch(curMovieIndex).movieInfo.fileName2;
        
        %Check if a fileName for a second movie has been specified
        if ~isempty(fileName) && ui.cboxShowStack2.Value
            %Check if movie has not been loaded already
            if isempty(movieStacks{2})
                movieStacks{2} = load_stack(batch(curMovieIndex).movieInfo.pathName2, fileName,ui);
                %Apply image processing with user defined filters
                ProcessStack2CB('StackLoaded')
            end
            %Update UI
            ui.textFileName.String = strcat('Filename:', fileName);
        else
            ui.cboxShowStack2.Value = 0;
            %Display filename of first stack
            fileName = batch(curMovieIndex).movieInfo.fileName;
            ui.textFileName.String = strcat('Filename:', fileName);
        end
        
    end

    function ProcessStack2CB(src,~)
        %User changed movie or changed settings in the table with stack2
        %filter settings
        
        filtered = movieStacks{2};
        
        filterOptions = ui.tableFilterOptions.Data;
        
        %Time average
        if filterOptions{2,1}
            nAvgFrames = filterOptions{2,3};
            divisor = floor(size(filtered,3)/nAvgFrames);
            boolRemaining = logical(mod(size(filtered,3),nAvgFrames));
            nFramesNewStack = divisor + boolRemaining;
            dummyStack = zeros(size(filtered,1),size(filtered,2),nFramesNewStack);
            
            for i=1:nFramesNewStack
                if i == nFramesNewStack  %Last sequence might have less frames
                    dummyStack(:,:,i) = mean(filtered(:,:,(i-1)*nAvgFrames+1:end),3);
                else
                    dummyStack(:,:,i) = mean(filtered(:,:,(i-1)*nAvgFrames+1:(i-1)*nAvgFrames+nAvgFrames),3);
                end
            end
            
            
            filtered = dummyStack;
        end
        
        %Moving time average
        if filterOptions{3,1}
            nAvgFrames = filterOptions{3,3};
            filtered = movmean(filtered,nAvgFrames,3);
        end
        
        %Rotate
        if filterOptions{1,1}
            filtered = fliplr(filtered);
        end
        
        %Gaussian Filter
        if filterOptions{4,1}
            nPixels = filterOptions{4,3};
            for i=1:size(filtered,3)
                filtered(:,:,i) = imgaussfilt(filtered(:,:,i),nPixels);
            end
        end
        
        %Background subtraction
        if filterOptions{5,1}
            nPixels = filterOptions{5,3};
            for i=1:size(filtered,3)
                se = strel('disk',nPixels);
                I = filtered(:,:,i);
                background = imopen(I,se);
                filtered(:,:,i) = I - background;
            end
        end
        
        movieStacks{3} = filtered;
        
        if ishandle(src) %User changed filter settings
            FrameSlider()
            AdjustBrightnessCB('autoAdjust')
        end
        
    end

    function DrawSubROICB(src,~)
        %User pressed "Draw" button
        try
            set(gcf,'currentch',char(1)) %Make sure current character is not esc key
            if src.BackgroundColor(1) == 1
                return
            end
            
            %Initialize ROI
            curROI = [];
            
            %---------Create ROI--------------
            nROIS = size(batch(curMovieIndex).subROI,2);
            curROIList = {};
            
            if verLessThan('matlab','9.5')
                ROIFreehandHandle = imfreehand();
                
                if ~isempty(ROIFreehandHandle)
                    curROI = getPosition(ROIFreehandHandle); %save ROI in List
                end
            else
                src.BackgroundColor = 'r';
                src.String = 'Press enter to Finish';
                drawnow
                if ui.popDrawMode.Value == 1
                    ROIFreehandHandle = drawfreehand();
                elseif ui.popDrawMode.Value == 2
                    %Use drawassisted
                    ROIFreehandHandle = drawassisted();
                elseif ui.popDrawMode.Value == 3
                    %Draw polygon
                    ROIFreehandHandle = drawpolygon();
                elseif ui.popDrawMode.Value == 4
                    %Draw ellipse
                    ROIFreehandHandle = drawellipse();
                end
                
                % Wait for the user to press escape
                if ~isempty(ROIFreehandHandle.Position)
                    
                    % Wait for the user to press enter
                    while true
                        w = waitforbuttonpress;
                        switch w
                            case 1 %keyboard press
                                key = get(gcf,'currentcharacter');
                                if key == 13  %User pressed enter key
                                    break
                                end
                        end
                    end
                    
                    %Get ROI coordinates
                    if ui.popDrawMode.Value == 4
                        binaryMask = createMask(ROIFreehandHandle);
                        boundariesBinaryMask = bwboundaries(binaryMask);
                        curROI(:,1) = boundariesBinaryMask{1}(:,2);
                        curROI(:,2) = boundariesBinaryMask{1}(:,1);
                    else
                        curROI = ROIFreehandHandle.Position;
                    end
                end
            end
            if ~isempty(curROI)
                %Delete handle to ROI
                delete(ROIFreehandHandle)
                
                %Connect start and end of ROI
                curROI = cat(1, curROI, curROI(1,:));
                
                %Restrict ROI to image size if user drew outside and
                %account for scaling factor
                height = batch(curMovieIndex).movieInfo.height;
                width = batch(curMovieIndex).movieInfo.width;
                scalingFactor = plotSettings.scalingFactor;
                
                curROI(curROI(:) < 0) = 0.5;
                curROI(curROI(:,2) > height*scalingFactor,2) = height*scalingFactor+0.5;
                curROI(curROI(:,1) > width*scalingFactor,1) =width*scalingFactor+0.5;
                
                %Append to current region list
                curROIList{end+1} = curROI./plotSettings.scalingFactor;
                
                if ~isempty(curROIList)
                    %Append to list of all sub-regions
                    batch(curMovieIndex).subROI{nROIS+1}{1} = curROIList;
                    RoiChanged(true);
                    FrameSlider()
                end
            end
            
            %Reset button appearance and update sub-region list in UI
            src.BackgroundColor = [.94 .94 .94];
            src.String = 'Hand-drawn';
        catch ex
            disp(ex)
            %Reset button appearance and update sub-region list in UI
            src.BackgroundColor = [.94 .94 .94];
            src.String = 'Hand-drawn';
        end
    end

    function DrawBrightnessROICB(~,~)
        %User pressed "Threshold" Button
        
        curFrameStack1 = round(ui.sliderFrame.Value);
        subROI = batch(curMovieIndex).subROI;
        filteredStack = double(movieStacks{3});
        
        if isempty(filteredStack)
            return
        end
        
        sizeStack2 = size(filteredStack);
        
        divisor = max([1, size(movieStacks{1},3)/size(filteredStack,3)]);
        curFrameStack2 = ceil(curFrameStack1/divisor);
        
        if ~isempty(batch(curMovieIndex).ROI)
            ROI = batch(curMovieIndex).ROI{1};
            ROImask = poly2mask(ROI(:,1),ROI(:,2),sizeStack2(1),sizeStack2(2));
        else
            sizeStack1 = size(filteredStack);
            ROImask = ones(sizeStack1(1),sizeStack1(2));
        end
        
        %Set pixels outside ROI to zero
        imageMask = filteredStack(:,:,curFrameStack2).*ROImask;
        
        %Get lowest pixel value inside ROI (start value of slider)
        minThreshold = floor(min(nonzeros(imageMask)));
        max1 = max(imageMask(:));
        otsu = max1*otsuthresh(imageMask(:));
        
        %Initialize
        minAreaSize = 0;
        boolFillHoles = false;
        
        canceled = 0;
        %-----------Create interface where user can drag two sliders for adjusting the intensity threshold
        d = dialog('Position',[300 300 500 200],'Name','Choose a threshold','WindowKeyPressFcn',@KeyPressFcnCB,'CloseRequestFcn',@CloseRequestCB);
        textTh1 = uicontrol('Units','normalized','Parent',d,'Style','text','Position',[.05 .84 .2 .10], 'String','Threshold: ','HorizontalAlignment','Left');
        editTh1 = uicontrol('Units','normalized','Parent',d,'Style','edit','Position',[.18 .87 .15 .08],'Tag','th1', 'String',num2str(otsu),'Callback',@EditThresholdCB);
        sliderFrameTh1 = uicontrol('Units','normalized','Enable','on','Style','slider','Position',[.05 .75 .9 .1],'Min',minThreshold,'Max',max(imageMask(:)),'Value',max(minThreshold,otsu),'SliderStep',[1/20 1/10]);
        addlistener(sliderFrameTh1,'Value','PostSet',@(~,~) SliderThListenerCB);
        
        textTh2 = uicontrol('Units','normalized','Parent',d,'Style','text','Position',[.05 .59 .2 .10],'String','Min. area size:','HorizontalAlignment','Left');
        editTh2 = uicontrol('Units','normalized','Parent',d,'Style','edit','Position',[.22 .62 .15 .08],'Tag','th2','String','0','Callback',@EditThresholdCB);
        sliderFrameTh2 = uicontrol('Units','normalized','Enable','on','Style','slider','Position',[.05 .5 .9 .1],'Min',0,'Max',10000,'Value',0,'SliderStep',[1/100000 1/10]);
        addlistener(sliderFrameTh2,'Value','PostSet',@(~,~) SliderThListenerCB);
        
        cboxFillHoles = uicontrol('Units','normalized','Style','checkbox','Position',[.05 .25 .9 .1],'String','Fill holes','Callback',@(~,~)SliderThListenerCB);
        
        uicontrol('Units','normalized','Parent',d,'Position',[.05 .1 .2 .1],'String','OK', 'Callback','delete(gcf)');
        uicontrol('Units','normalized','Parent',d,'Position',[.3 .1 .2 .1],'String','Cancel', 'Callback',@CloseRequestCB);
        SliderThListenerCB()
        
        %------------Interactive interface to choose the boundary. Is executed until user presses "OK" button.
        function SliderThListenerCB()
            minThreshold = sliderFrameTh1.Value; %Get threshold values from sliders
            editTh1.String = num2str(round(minThreshold)); %Display threshold value in interface
            
            minAreaSize = round(sliderFrameTh2.Value);
            editTh2.String = num2str(minAreaSize); %Display threshold value in interface
            
            binary = imageMask >= minThreshold; %Calculate bindary image according to the intensity threshold
            
            binary=bwareaopen(binary,minAreaSize);
            
            boolFillHoles = cboxFillHoles.Value;
            if boolFillHoles
                binary = imfill(binary,'holes');
            end
            
            curBoundaries = cellfun(@fliplr,bwboundaries(binary),'UniformOutput',false); %Create boundary with respect to the binary image
            hold (ui.axes1, 'on');
            
            %Plot existing ROIs
            nROIhandles = 1;
            nROIs = length(batch(curMovieIndex).subROI);
            
            for roiIdx = 1:nROIs
                if length(batch(curMovieIndex).subROI{roiIdx}) > 1
                    boundaryFrame = curFrameStack1;
                else
                    boundaryFrame = 1;
                end
                
                for f=1:length(batch(curMovieIndex).subROI{roiIdx}{boundaryFrame})
                    curROI = batch(curMovieIndex).subROI{roiIdx}{boundaryFrame}{f};
                    if nROIhandles>numel(subROIHandle)
                        subROIHandle(nROIhandles) = plot(ui.axes1,curROI(:,1), curROI(:,2), '-','Color',plotSettings.subRoiColors(mod(roiIdx-1,7)+1,:),'Linewidth',1);
                    else
                        set(subROIHandle(nROIhandles),'xdata',curROI(:,1),'ydata', curROI(:,2),'Color',plotSettings.subRoiColors(mod(roiIdx-1,7)+1,:));
                    end
                    nROIhandles = nROIhandles +1 ;
                end
            end
            
            %Plot new ROI
            
            for f=1:length(curBoundaries)
                b1 = curBoundaries{f};
                if nROIhandles>numel(subROIHandle)
                    subROIHandle(nROIhandles) = plot(b1(:,1), b1(:,2), '-','Color',plotSettings.subRoiColors(mod(nROIs,7)+1,:),'Linewidth',1,'parent',ui.axes1);
                else
                    set(subROIHandle(nROIhandles),'xdata',b1(:,1),'ydata', b1(:,2),'Color',plotSettings.subRoiColors(mod(nROIs,7)+1,:));
                end
                nROIhandles = nROIhandles + 1;
            end %Plot boundaries into the image
            
            %Delete unecessary line handles
            for q = nROIhandles:numel(subROIHandle)
                delete(subROIHandle(nROIhandles));
                subROIHandle(nROIhandles) = [];
            end
            
            hold (ui.axes1, 'off');
        end
        
        function EditThresholdCB(src,~)
            if strcmp(src.Tag, 'th1')
                sliderFrameTh1.Value = max(sliderFrameTh1.Min,min(sliderFrameTh1.Max,str2double(src.String)));
            else
                sliderFrameTh2.Value = max(sliderFrameTh2.Min,min(sliderFrameTh2.Max,str2double(src.String)));
            end
        end
        
        function CloseRequestCB(~,~)            
            canceled = 1;
            delete(gcf)
        end
        
        function KeyPressFcnCB(~,event)
            if strcmp(event.Key, 'escape')
                delete(gcf)
                canceled = 1;
            end
        end
        
        uiwait(d); % Wait for d to close before running to completion
        
        if ~canceled
            %---------------Create boundaries for each frame taking into account the user defined thresholds
            
            nROIS = size(subROI,2);
            
            nFramesFirstStack = size(movieStacks{1},3);
            nFramesSecondStack = size(filteredStack,3);
            divisor = max([1, size(movieStacks{1},3)/nFramesSecondStack]);
            
            lastFrame = 0;
            
            for i=1:nFramesFirstStack %Iterate through frames and create boundaries for each frame
                curMovieFrame = min(nFramesSecondStack,ceil(i/divisor));
                if lastFrame ~= curMovieFrame
                    curFrame = filteredStack(:,:,curMovieFrame).*ROImask;
                    
                    binary1 = curFrame >= minThreshold;
                    
                    binary1=bwareaopen(binary1,minAreaSize);
                    
                    if boolFillHoles
                        binary1 = imfill(binary1,'holes');
                    end
                    
                    boundary1 = bwboundaries(binary1);
                    
                    lastFrame = curMovieFrame;
                end
                
                subROI{nROIS+1}{i} = cellfun(@fliplr,boundary1,'UniformOutput',false);
            end
            
            batch(curMovieIndex).subROI = subROI;
            RoiChanged(true);
            FrameSlider()
        end
        
        FrameSlider()
        
    end

    function DeleteSubROICB(~,~)
        %User pressed "Delete selected region" button
        subROI = batch(curMovieIndex).subROI;
        nSubRois = length(subROI);
        
        selection = ui.listSubRoi.Value;
        
        selection = selection(selection ~=1) - 1;
        nSelected = numel(selection);
        
        selectionIdx = true(nSubRois,1);
        selectionIdx(selection) = false;
                
        if ~isempty(subROI) && nSelected > 0
            %Remove selected sub-region and save in batch struct
            batch(curMovieIndex).subROI = subROI(selectionIdx);
            RoiChanged(true);
            FrameSlider()
        end
    end

    function MergeSubROICB(~,~)
        %User pressed "Merge selection" button
        
        %Get number of subRois
        nSubRois = length(ui.listSubRoi.String)-1;
        
        %Get selected regions
        selection = ui.listSubRoi.Value-1;
        
        %Ignore region 1 (tracking-region)
        selection = selection(selection ~=0);
        
        %Get number of selected regions
        nSelected = numel(selection);
        
        %Create logical array of selected sub-rois
        selectionIdx = false(nSubRois,1);
        selectionIdx(selection) = true;
        
        %Get subrois from batch
        subROI = batch(curMovieIndex).subROI;
        newSubROI = {};
        
        for idx = selection
            %Merge selected sub rois
            if length(subROI{idx}) > 1
                msgbox('Cannot merge regions drawn with threshold.')
                return
            end
            newSubROI =  [newSubROI, subROI{idx}{1}];
        end
        
        %Get subrois that habe not been merged
        subROI = subROI(~selectionIdx);
        
        %Get number of subrois after merge
        nSubRoisAfterMerge = nSubRois - nSelected + 1;
        
        %Save new subroi at the end of the sub roi array
        subROI{nSubRoisAfterMerge}{1} = newSubROI;
        
        %Save in batch structure
        batch(curMovieIndex).subROI = subROI;
        
        RoiChanged(true);
        FrameSlider()
        
    end

    function RoiChanged(boolSaveToFile)
        
        %Get sub regions from batch structure
        subROI = batch(curMovieIndex).subROI;
        
        %Create sub region list entries
        newList = cell(length(subROI),1);
        if ~isempty(subROI)
            newList{1} = 'Region 1 (tracking-region)';
            for roiIdx = 2:length(subROI)+1
                newList{roiIdx} = ['Region ',num2str(roiIdx)];
            end
        end
        
        %Set selected sub-region to first entry
        ui.listSubRoi.Value = 1;
        
        %Save new sub region list
        ui.listSubRoi.String = newList;
        
        if boolSaveToFile
            %Save to .roi file
            
            %Create ROI folder if it doesn't exists
            roiPath = fullfile(batch(curMovieIndex).movieInfo.pathName,'ROI');
            if ~exist(roiPath, 'dir')
                mkdir(roiPath)
            end
            
            %Save ROI and sub-rois in .roi file
            [~,fileName,~] = fileparts(batch(curMovieIndex).movieInfo.fileName);
            fullROIname = fullfile(roiPath,strcat(fileName,'.roi'));
            ROI = batch(curMovieIndex).ROI;
            save(fullROIname,'ROI','subROI');
        end
        
        
        ScaleToRoiCB('')
    end

%Plotting functions

    function cursorText = DataCursorCB(~,eventHandle)
        %User used the data tip tool and clicked into the plot window
        
        %Reset current track selection
        plotSettings.curTrackSel = 0;
        if strcmp(ui.dataCoursor.Enable,'on')
            
            %Get current frame number
            curFrame = round(ui.sliderFrame.Value);
            
            %Get handle to the selected object
            graphObjHandle = eventHandle.Target;
            
            %Get position value of the selected object
            pos = eventHandle.Position;
            
            
            
            if isgraphics(graphObjHandle,'image')
                %User clicked on a pixel so show coordinates and pixel intensity
                cursorText = {['X: ',num2str(pos(1)),'   Y: ',num2str(pos(2))],...
                    ['Intensity: ',num2str(imageHandle.CData(pos(2),pos(1)))]};
            elseif graphObjHandle == spotsHandle
                %User clicked on a spot
                
                %Scale the position with the scaling factor
                pos = pos./plotSettings.scalingFactor;
                
                %Get spots in current frame
                spotsInFrame = batch(curMovieIndex).results.spotsAll{curFrame};
                
                %Find the index of the selected spot by comparing the selected position
                %with all spot positions in the current frame
                spotIdx = spotsInFrame(:,1) == pos(1) & spotsInFrame(:,2) == pos(2);
                
                %Get data on the selected spot
                selectedSpot = spotsInFrame(spotIdx,:);
                
                %Calculate the SNR
                [SNRspots, intPeak] = calc_snr(spotsInFrame,movieStacks{1}(:,:,curFrame));
                
                %Display data tip text
                cursorText = {['X: ',num2str(round(pos(1),1)),'   Y: ',num2str(round(pos(2),1))],...
                    ['Fitted int.: ',num2str(round(selectedSpot(3)))],...
                    ['Max. int.: ',num2str(intPeak(spotIdx))],...
                    ['SNR: ',num2str(SNRspots(spotIdx))]};
                
            elseif any(graphObjHandle == trackHandles)
                %User clicked on a track
                
                %Scale the position with the scaling factor
                pos = pos./plotSettings.scalingFactor;
                
                %Get the trackId of the selected track
                trackId = eventHandle.Target.UserData;
                
                %Save the trackID in the plotSettings in case we need it
                %for the track explorer
                plotSettings.curTrackSel = trackId;
                
                %Get the selected track
                curTrack = batch(curMovieIndex).results.tracks{trackId};
                
                %Find the index of the selected spot of the track
                spotIdx = curTrack(:,2) == pos(1) & curTrack(:,3) == pos(2);
                
                %Get the detection number of the selected spot of current track
                spotNum = find(spotIdx);
                
                %Get the data of the selected spot of current track
                curSpotInTrack = curTrack(spotIdx,:);
                
                %Display the data tip text
                cursorText = {['X: ',num2str(round(pos(1),1)),'   Y: ',num2str(round(pos(2),1))],...Â´
                    ['Track ID: ', num2str(trackId)],...
                    ['Detection: ', num2str(spotNum), '/',num2str(numel(spotIdx))],...
                    ['Fitted int: ',num2str(round(curSpotInTrack(4)))],...
                    ['Mean jump dist.: ', num2str(round(batch(curMovieIndex).results.meanJumpDists(trackId),2))]};
            else
                %User clicked anything else eg. the ROI so just display the
                %coordinates
                cursorText = {['X: ',num2str(pos(1)),'   Y: ',num2str(pos(2))]};
            end
        end
    end

    function AdjustBrightnessCB(src,~)
        %Executed whenever the lookup table of the displayed image has to be adjusted
        
        %Get the currently displayed image and save it as single datatype
        I = single(imageHandle.CData(:));
        
        if ~isempty(I)
            %Get minimum and maximum pixel value
            [frameMin, frameMax] = bounds(I);
            
            if strcmp(src, 'autoAdjust')
                %Either user clicked the "Auto adjust brightness" button or
                %function was called from within another function
                
                %Get pixels with values above 0
                pixAboveZero = I(I > 0);
                
                if ui.popShowImage.Value > 4 && ui.cboxShowZProjection.Value
                    %Detection or jump distance map is currently displayed
                    %so set dark value to zero

                    sliderBlack = 0;
                else
                    %Calculate the black threshold slider value so that the
                    %lowest 3% of pixels are displayed black
                    sliderBlack = max(0,(prctile(pixAboveZero, 3)-frameMin)/frameMax);
                    
                end
                
                %Calculate the white threshold slider value so that the
                %highest 0.3% of pixels are displayed white
                sliderWhite = min(1,prctile(pixAboveZero, 99.7)/frameMax);
                
                %Calculate pixel threshold values from the slider values
                lowerThres = frameMax * sliderBlack+frameMin;
                upperThres = frameMax * sliderWhite;
                
            elseif strcmp(src,'textInput')
                %User entered a manual value in the edit fields to set the
                %brightnes
                
                %Get the treshold values from the edit fields of the ui
                lowerThres = str2double(ui.editBlack.String);
                upperThres = str2double(ui.editWhite.String);
                
                
                %Calculate the corresponding slider values
                sliderBlack = max(0,(lowerThres-frameMin)/frameMax);
                sliderWhite = min(1,upperThres/frameMax);
                
            elseif strcmp(src,'slider')
                %Slider value was changed
                
                %Calculate pixel threshold values from the slider values
                lowerThres = frameMax * ui.sliderBlack.Value+frameMin;
                upperThres = frameMax * ui.sliderWhite.Value;
                
                %Get the treshold values from the edit fields of the ui
                editBlack = str2double(ui.editBlack.String);
                editWhite = str2double(ui.editWhite.String);
                
                %Make sure that values are not overwritte if the user entered
                %a lower value in the edit field than the lowest possible
                %slider value
                if ui.sliderBlack.Value == 0 && editBlack < lowerThres
                    lowerThres = editBlack;
                end
                
                %Make sure that values are not overwritte if the user entered
                %a higher value in the edit field than the highest possible
                %slider value
                if ui.sliderWhite.Value == 1 && editWhite > upperThres
                    upperThres = editWhite;
                end
                
                sliderBlack = [];
                sliderWhite = [];
            end
            
            %Make sure lower threshold is lower than upper threshold
            if lowerThres >= upperThres
                upperThres = lowerThres+0.001;
            end
            
            %Apply threshold values to axis
            ui.axes1.CLim = [lowerThres upperThres];
            
            %Write new pixel threshold values into the edit fields of the ui
            ui.editBlack.String = num2str(lowerThres,'%.1f');
            ui.editWhite.String = num2str(upperThres,'%.1f');
            
            %Update slider value if function was not called by a slider movement
            if ~isempty(sliderBlack)
                ui.sliderBlack.Value = sliderBlack;
                ui.sliderWhite.Value = sliderWhite;
            end
        end
    end

    function GetTracksInFrame()
        %Creates a matrix of logicals indicating which track has to be
        %shown in each frame
        
        nFrames = batch(curMovieIndex).movieInfo.frames;
        tracks = batch(curMovieIndex).results.tracks;
        nTracks = batch(curMovieIndex).results.nTracks;
        nFramesTrackIsVisible = plotSettings.nFramesTrackIsVisible;
        tracksInFrame = cell(nFrames,1);
        frameTrackMatrix = false(nTracks,nFrames);
        
        for trackIdx = 1:nTracks
            frameTrackMatrix(trackIdx, tracks{trackIdx}(1,1):min(tracks{trackIdx}(end,1)+nFramesTrackIsVisible, nFrames)) = true;
        end
        
        for frameIdx = 1:nFrames
            tracksInFrame{frameIdx} = find(frameTrackMatrix(:,frameIdx)).';
        end
        
        plotSettings.tracksInFrame = tracksInFrame;
        plotSettings.tracksChanged = 1;
        SetTrackColors()
    end

    function PlayCB(~,~)
        %User pressed "play" button
        
        if strcmp(ui.btnPlay.String,'Play')
            %User pressed "play" so start timer
            tic
            ui.btnPlay.String = 'Pause';
            timerObj.Period = round(1/str2double(ui.editFPS.String),3);
            timerObj.UserData = ui.sliderFrame.Value; %Start frame
            start(timerObj);
        else
            %User pressed "Pause" so stop timer
            ui.btnPlay.String = 'Play';
            stop(timerObj);
        end
    end

    function Playing(~,~)
        %Called by the timer object in the frequency specified by "FPS"
        
        if timerObj.Period ~= round(1/str2double(ui.editFPS.String),3)
            %User changed the "FPS" field
            stop(timerObj)
            timerObj.Period = round(1/str2double(ui.editFPS.String),3);
            timerObj.UserData = ui.sliderFrame.Value; %Start frame
            tic
            start(timerObj);
        end
        
        newFrame = floor(toc/timerObj.Period) + timerObj.Userdata;
        
        if newFrame >= size(movieStacks{1},3)
            stop(timerObj);
            ui.sliderFrame.Value = 1;
            timerObj.UserData = 1;
            tic
            start(timerObj);
        else
            ui.sliderFrame.Value = newFrame;
            drawnow %Make sure drawing has finished
        end
    end

    function ScaleToRoiCB(src,~)
        %Executed when the user presses "Scale to ROI" or "Scale to subROI"
        %toolbar icon, a new movie is loaded, a region or subregion was changed
        %or the scaling factor changed
        
        %Get the currently set scaling factor
        scalingFactor = plotSettings.scalingFactor;
        
        %Get movie dimensions
        width = batch(curMovieIndex).movieInfo.width;
        height = batch(curMovieIndex).movieInfo.height;
        
        %Calculate minimum and maximum axis limits
        minX = 0.5;
        maxX = width*scalingFactor+0.5;
        minY = 0.5;
        maxY = height*scalingFactor+0.5;
        
        if ~isempty(src) && strcmp(src.Tag, 'ZoomToSquare')
            
            prompt = {'Zoom into a rectangular region with specific window size. Specify the window size in pixels below, press ok and point to a region where you want to zoom in.'};
            dlgtitle = 'Zoom into rectangular region';
            answer = inputdlg(prompt,dlgtitle,1,{num2str(plotSettings.rectZoomWinSize)});
            
            if isempty(answer)
                %Canceled
                return
            end
            
            %Save selected pixelSize in plotSettings struct
            plotSettings.rectZoomWinSize = str2double(answer);   
            
            halfWindowSize = str2double(answer)/2;    
            %Allow user to select point where to zoom in
            pointHandle = drawpoint(ui.axes1);
            
            %Get selected pixel position
            selectedPixel = round(pointHandle.Position);         
            
            %Delete handle of point selection tool
            delete(pointHandle)
            
            %Get window positions
            minX = round(selectedPixel(1)-halfWindowSize);
            maxX = round(selectedPixel(1)+halfWindowSize);
            minY = round(selectedPixel(2)-halfWindowSize);
            maxY = round(selectedPixel(2)+halfWindowSize);
            
            ui.scaleToRoiTool.State = 'off';
            ui.scaleToSubRoiTool.State = 'off';
            ui.scaleToSquare.State = 'off';
            
        elseif (isempty(src) || strcmp(src.Tag, 'ScaleToRoi')) && strcmp(ui.scaleToRoiTool.State, 'on')
            %Scale to roi toolbar icon is enabled
            
            %Make sure not both "scale to roi" and "scale to subroi"
            %are enabled simultaneously
            ui.scaleToSubRoiTool.State = 'off';
            ui.scaleToSquare.State = 'off';
            
            %Check if ROI exists
            if ~isempty(batch(curMovieIndex).ROI)
                %Get ROI
                ROI = batch(curMovieIndex).ROI{1};
                
                %Shrink minimum and maximum axis limits to the ROI
                minX = (min(ROI(:,1))-3)*scalingFactor;
                maxX = (max(ROI(:,1))+3)*scalingFactor;
                minY = (min(ROI(:,2))-3)*scalingFactor;
                maxY = (max(ROI(:,2))+3)*scalingFactor;
            end
            
            
        elseif (isempty(src) || strcmp(src.Tag, 'ScaleToSubRoi')) && strcmp(ui.scaleToSubRoiTool.State, 'on')
            %Scale to subroi toolbar icon is enabled
            
            %Make sure not both "scale to roi" and "scale to subroi"
            %are enabled simultaneously
            ui.scaleToRoiTool.State = 'off';
            ui.scaleToSquare.State = 'off';
            
            %Get subROI
            subROI = batch(curMovieIndex).subROI;
            
            %Check if subROI exists
            if ~isempty(subROI)
                
                %Initialize the min and max values
                minX = width*scalingFactor+0.5;
                maxX = 0.5;
                minY = height*scalingFactor+0.5;
                maxY = 0.5;
                
                %Iterate through all parts of all subregions and find
                %biggest rectangle that contains all sub-regions
                for idxSubRoiNum = 1:length(subROI)
                    %Get current subRegion
                    curSubRoi = subROI{idxSubRoiNum};
                    
                    %Get number of frames of the current sub-region
                    nFrames = length(curSubRoi);
                    
                    %Iterate through all frames of the current subRegion
                    %(either nFrames=1 in case of a hand-drawn subregion or
                    %nFrames = number of frames in the movie for
                    %threshold-drawn subregion.
                    for idxFrame = 1:nFrames
                        %Get the subRegion of the current frame
                        curFrameSubRoi = curSubRoi{idxFrame};
                        
                        %Get the number of separated sub-regions in the
                        %current frame
                        nParts = length(curFrameSubRoi);
                        
                        %Iterate through all parts of the current sub-region of the current frame
                        for idxPart = 1:nParts
                            %Get the current sub-region part
                            curPartSubRoi = curFrameSubRoi{idxPart};
                            
                            %Adjust the axis limits
                            minX = min(minX, (min(curPartSubRoi(:,1))-3)*scalingFactor);
                            maxX = max(maxX, (max(curPartSubRoi(:,1))+3)*scalingFactor);
                            minY = min(minY, (min(curPartSubRoi(:,2))-3)*scalingFactor);
                            maxY = max(maxY, (max(curPartSubRoi(:,2))+3)*scalingFactor);
                        end
                    end
                end
                
            end
        elseif (isempty(src) || strcmp(src.Tag, 'ScaleToSquare')) && strcmp(ui.scaleToSquare.State, 'on')
            %Scale to subroi toolbar icon is enabled
            
            %Make sure not both "scale to roi" and "scale to subroi"
            %are enabled simultaneously
            ui.scaleToRoiTool.State = 'off';
            ui.scaleToSubRoiTool.State = 'off';
            
            %Get smaller side
            
            nPixels = min(maxY, maxX);
            
            maxX = nPixels;
            maxY = nPixels;
        end
        
        %Save new axis limits
        ui.axes1.XLim = [minX maxX];
        ui.axes1.YLim = [minY maxY];
        
    end

    function ZoomCB(~,~)
        %Make sure axis limits are correct so user doesn't zoom out too far
        h = zoom;
        
        [height, width, ~] = size(imageHandle.CData);
        
        if (ui.axes1.XLim(1)<0 || ui.axes1.YLim(1) < 0 || ui.axes1.XLim(2)>width || ui.axes1.YLim(2)>height) %strcmp(h.Direction, 'out')
            ui.axes1.XLim = [0.5 width+0.5]; %Set Zoom
            ui.axes1.YLim = [0.5 height+0.5];
        end
        
        FrameSlider()
    end

    function CreateZProjection(src,~)
        %This function creates an image of whatever is selected in the
        %pop-up menu in the bottom right of the main gui (eg. Standard
        %deviation, detection map, wavelet filtered image etc.). The image
        %is saved in the variable movieStacks{4} and displayed when the
        %checkbox besides the popup menu is enabled.
        
        
        if ui.popShowImage.Value < 5 || ~ui.cboxShowZProjection.Value
            %Hide scaling factor ui elements
            ui.editScalingFactor.Visible = 'off';
            ui.textScalingFactor.Visible = 'off';
            scalingFactor = 1;
            
            %Check if movie has to be loaded
            if isempty(movieStacks{1})
                %Load movie
                movieStacks{1} = load_stack(batch(curMovieIndex).movieInfo.pathName, batch(curMovieIndex).movieInfo.fileName, ui);
            end
            
            if ui.cboxShowStack2.Value
                %"Show stack 2" checkbox is checked so use stack 2
                curStack = movieStacks{3};
            else
                %"Show stack 2" checkbox is not checked so use stack 1
                curStack = movieStacks{1};
            end
        else
            %Show scaling factor ui elements
            ui.editScalingFactor.Visible = 'on';
            ui.textScalingFactor.Visible = 'on';
            
            %Factor by which the amount of pixels in each direction will increase
            scalingFactor = str2double(ui.editScalingFactor.String);
        end
        
        if ui.cboxShowZProjection.Value
            switch ui.popShowImage.Value
                case 1
                    %%Show standard deviation
                    movieStacks{4} = std(single(curStack),0,3);
                case 2
                    %%Show average intensity
                    movieStacks{4} = mean(curStack,3);
                case 3
                    %%Show maximum intensity
                    movieStacks{4} = max(curStack,[],3);
                case 4
                    %%Show wavelet filtered image (calculated in FrameSlider)
                case 5
                    %%No Image
                    movieStacks{4} = zeros(batch(curMovieIndex).movieInfo.height*scalingFactor,batch(curMovieIndex).movieInfo.width*scalingFactor);
                case 6
                    %% Show detection map inlcuding non-linked spots
                    
                    %Get movie dimensions
                    movieInfo = batch(curMovieIndex).movieInfo;
                    height = round(movieInfo.height*scalingFactor);
                    width = round(movieInfo.width*scalingFactor);
                    
                    %Initialize talm-map with zeros
                    talmMap = zeros(height,width);
                    
                    %Get spots from batch structure
                    spotsAll = batch(curMovieIndex).results.spotsAll;
                    
                    %Get number of frames
                    nFrames = length(spotsAll);
                    
                    %Iterate through all frames
                    for frameIdx = 1:nFrames
                        %Get spots in current frame
                        curSpots = spotsAll{frameIdx};
                        
                        %Iterate through all spots in current frame
                        for spotIdx = 1:size(curSpots,1)
                            %Get spot coordinates
                            curSpotY = round(curSpots(spotIdx,1)*scalingFactor);
                            curSpotX = round(curSpots(spotIdx,2)*scalingFactor);
                            %Increase pixel value of the spot coordinates by one
                            talmMap(curSpotX,curSpotY) = talmMap(curSpotX,curSpotY) + 1;
                        end
                    end
                    
                    %Save talm map in the movieStacks variable
                    movieStacks{4} = talmMap;
                case 7
                    %% Show detection map without non-linked spots
                    
                    %Get movie dimensions
                    movieInfo = batch(curMovieIndex).movieInfo;
                    height = round(movieInfo.height*scalingFactor);
                    width = round(movieInfo.width*scalingFactor);
                    
                    %Initialize talmmap with zeros
                    talmMap = zeros(height,width);
                    
                    %Get tracks from batch structure because we only want
                    %tracked spots in our localization map
                    tracks = batch(curMovieIndex).results.tracks;
                    
                    %Iterate through tracks
                    for trackIdx = 1:length(tracks)
                        %Get current track
                        curTrack = tracks{trackIdx};
                        %Iterate through all spots in current track
                        for spotIdx = 1:size(curTrack,1)
                            %Get spot coordinates
                            posY = round(curTrack(spotIdx,2)*scalingFactor);
                            posX = round(curTrack(spotIdx,3)*scalingFactor);
                            %Increase pixel value of the spot coordinates by one
                            talmMap(posX,posY) = talmMap(posX,posY) + 1;
                        end
                    end
                    
                    %Save talm map in the movieStacks variable
                    movieStacks{4} = talmMap;
                case 8
                    %% Show jump distance map
                    
                    %Get movie dimensions
                    movieInfo = batch(curMovieIndex).movieInfo;
                    height = round(movieInfo.height*scalingFactor);
                    width = round(movieInfo.width*scalingFactor);
                    
                    %Initialize jump distance map with zeros
                    jumpDistMap = zeros(height,width);
                    
                    %Initialize normalization map with zeros
                    normMap = zeros(height,width);
                    
                    tracks = batch(curMovieIndex).results.tracks;
                    jumpDistances = batch(curMovieIndex).results.jumpDistances;
                    
                    %Between all the pixels inbetween to jumping points
                    %insert the corresponding jump distance and normalize
                    %it by the amount of events in each pixel
                    for trackIdx = 1:length(tracks)
                        curTrack = tracks{trackIdx};
                        curTrackJumpDistances = jumpDistances{trackIdx};
                        
                        for spotIdx = 1:size(curTrack,1)-1
                            
                            %Make sure to not display gap jump
                            if curTrack(spotIdx+1,1)-curTrack(spotIdx,1) == 1
                                
                                %Spot position 1 in previous frame
                                curSpotY1 = curTrack(spotIdx,2)*scalingFactor;
                                curSpotX1 = curTrack(spotIdx,3)*scalingFactor;
                                %Spot position 2 (in next frame)
                                curSpotY2 = curTrack(spotIdx+1,2)*scalingFactor;
                                curSpotX2 = curTrack(spotIdx+1,3)*scalingFactor;
                                
                                
                                %Create a line from point 1 to point 2
                                spacing = .4;
                                numSamples = max(1,ceil(sqrt((curSpotX2-curSpotX1)^2+(curSpotY2-curSpotY1)^2) / spacing));
                                x = linspace(curSpotX1, curSpotX2, numSamples);
                                y = linspace(curSpotY1, curSpotY2, numSamples);
                                
                                %Round positions on the line to the next pixel
                                xy = round([x',y']);
                                
                                %Find line positions which where rounded to the same pixel
                                dxy = abs(diff(xy, 1));
                                duplicateRows = [0; sum(dxy, 2) == 0];
                                
                                if size(xy,1) == 1
                                    %Spot stayed in same pixel
                                    duplicateRows = 0;
                                end
                                
                                %Get final pixel positions on the line between
                                %point 1 and point 2
                                finalxy = xy(~duplicateRows,:);
                                finalx = finalxy(:, 1);
                                finaly = finalxy(:, 2);
                                
                                %Distance between spots
                                dist = curTrackJumpDistances(spotIdx);
                                
                                %For every pixel on a line between the two
                                %points add the jump distance between the two
                                %points
                                for pixelIdx = 1:numel(finalx)
                                    curSpotX = finalx(pixelIdx);
                                    curSpotY = finaly(pixelIdx);
                                    jumpDistMap(curSpotX,curSpotY) = jumpDistMap(curSpotX,curSpotY) + dist;
                                    %Create map containing the amount of events
                                    %in each pixel which will be used for
                                    %normalization
                                    normMap(curSpotX,curSpotY) = normMap(curSpotX,curSpotY) + 1;
                                end
                            end
                        end
                    end
                    %Replace zeros with -1 to prevent division by 0
                    
                    normMap(normMap == 0) = -1;
                    
                    %Normalize the accumulated jump distances in each pixel with the normalization map
                    jumpDistMap = jumpDistMap./(normMap);
                    
                    movieStacks{4} = jumpDistMap;
                    
            end
        end
        
        %Set tracksChanged value to 1 so that plotted tracks are renewed
        %afterwards in the FrameSlider function
        plotSettings.tracksChanged = 1;
        plotSettings.scalingFactor = scalingFactor;
        
        if ishandle(src)
            %Function was called by a user interaction
            
            %Adjust the axis limits to the currently used scaling factor
            ScaleToRoiCB('')
            %Update plot
            FrameSlider()
            %Auto adjust brightness
            AdjustBrightnessCB('autoAdjust')
        end
        
    end

    function PlotSettingsChangedCB(src,~)
        %Called whenever the user changed anything inside the "Plot
        %properties" panel
        
        switch src.Tag
            case 'nFramesTrackIsVisible'
                switch ui.popTracks.Value
                    case 1
                        %"Show tracks in range"
                        ui.textTrackLength.Visible = 'on';
                        ui.editTrackLength.Visible = 'on';
                        plotSettings.nFramesTrackIsVisible = str2double(ui.editTrackLength.String);
                        GetTracksInFrame()
                    case {2,3,4}
                        %"Show all tracks" or "Show initial positions"
                        ui.textTrackLength.Visible = 'off';
                        ui.editTrackLength.Visible = 'off';
                end
                
            case 'advanced'
                %User pressed "advanced plot properties" button
                newCurPlot = plot_properties(plotSettings);
                
                if ~isempty(newCurPlot)
                    plotSettings = newCurPlot;
                    if plotSettings.ITM
                        ui.textNDarkForLong.String = 'Min. #darktimes for long regime';
                    else
                        ui.textNDarkForLong.String = 'Min. tracklength for long regime';
                    end                    
                end
                                
            case 'scalebar'
                %User pressed button
                newCurPlot = set_scalebar(ui, scalebarHandle, scalebarTextHandle, plotSettings);
                if ~isempty(newCurPlot)
                    plotSettings = newCurPlot;
                end
                
                %Get scalebar coordinates
                xPos = plotSettings.sbPos(1);
                yPos = plotSettings.sbPos(2);
                xDim = plotSettings.sbDim(1);
                yDim = plotSettings.sbDim(2);
                
                %Write scalebar coordinates into the handle for the scalebar
                %graphics object
                set(scalebarHandle,'xdata',[xPos, xPos, xPos+xDim, xPos+xDim, xPos],...
                    'ydata',[yPos, yPos+yDim, yPos+yDim, yPos, yPos]);
                
                %Get settings for the scalebar text
                scalebarTextHandle.String = plotSettings.sbTextString;
                scalebarTextHandle.Position = [xPos, yPos+yDim+1];
                scalebarTextHandle.FontSize = plotSettings.sbTextSize;
                scalebarTextHandle.Position = plotSettings.sbTextPos;
                
            case 'lut'
        end
        
        
        colMapName = ui.popLut.String{ui.popLut.Value};
        
        %Create colormap for plotting
        switch colMapName
            case {'gray'; 'jet'; 'parula'; 'hot'}
                colMap = colormap(ui.axes1,colMapName);
            case 'inferno'
                colMap = inferno;
            case 'magma'             
                colMap = magma;
            case 'plasma'
                colMap = plasma;
            case 'viridis'
                colMap = viridis;
        end
        
        
        if plotSettings.invertColMap
            colMap = flipud(colMap);
        end
        
        switch plotSettings.bgColor
            case 'black'
                colMap(1,:) = [0 0 0];
            case 'white'
                colMap(1,:) = [1 1 1];
        end
        
        colormap(ui.axes1,colMap);
        
        plotSettings.colMap = colMap;        
        plotSettings.tracksChanged = 1;
        drawnow
        SetTrackColors()
        FrameSlider()
    end

    function SetTrackColors()
        %Creates a vector with containing the colors for each track
        
        %User changed colors of tracks
        
        ui.textNTrackColors.Visible = 'off';
        ui.editNTrackColors.Visible = 'off';
        ui.panelTrackLengthRegimes.Visible = 'off';
        ui.panelJumpDistRegimes.Visible = 'off';
        
        
        colorMap = parula(64);
        
        switch ui.popColoredTrackLengths.Value
            case 1
                %% Random colored tracks
                ui.textNTrackColors.Visible = 'on';
                ui.editNTrackColors.Visible = 'on';
                
                nColors = str2double(ui.editNTrackColors.String);
                nTracks = batch(curMovieIndex).results.nTracks;
                trackColors = repmat(distinguishable_colors(nColors, {'k'}), ceil(nTracks/nColors) ,1);
                plotSettings.trackColors = trackColors(1:nTracks,:);
                
            case 2
                %% Colored by frame of appearance
                
                values = batch(curMovieIndex).results.startEndFrameOfTracks;
                
                if isempty(values)
                    plotSettings.trackColors = 0;
                    return
                end
                
                nFrames = batch(curMovieIndex).movieInfo.frames;
                values = values(:,1);
                
                trackColors = interp1(linspace(1,nFrames,length(colorMap)),colorMap,values); % map color to y values
                trackColors = uint8(trackColors*255); % need a 4xN uint8 array
                plotSettings.trackColors = trackColors;
                
            case 3
                %% Color-coded track lengths/mean jump distance
                values = batch(curMovieIndex).results.trackLengths;
                
                if isempty(values)
                    plotSettings.trackColors = 0;
                    return
                end
                
                
                trackColors = interp1(linspace(min(values),max(values),length(colorMap)),colorMap,values); % map color to y values
                trackColors = uint8(trackColors*255); % need a 4xN uint8 array
                plotSettings.trackColors = trackColors;
            case 4
                %% Color-coded track length regimes
                
                ui.panelTrackLengthRegimes.Visible = 'on';
                nDarkForLong = str2double(ui.editNDarkForLong.String);
                trackLengths = batch(curMovieIndex).results.trackLengths;
                nTracks = batch(curMovieIndex).results.nTracks;
                if plotSettings.ITM
                    tracks = batch(curMovieIndex).results.tracks;
                    longIdx = false(length(tracks),1);
                    nBrightFrames = 2;
                    for trackID = 1:length(tracks)
                        curTrack = tracks{trackID};
                        longIdx(trackID) = floor((curTrack(end,1)-1)/nBrightFrames)-floor((curTrack(1,1)-1)/nBrightFrames)>= nDarkForLong;
                    end
                else
                    longIdx = trackLengths >= nDarkForLong;
                end
                shortIdx = ~longIdx;
                nLong = sum(longIdx);
                nShort = sum(shortIdx);
                
                colorLow = [1 0 0];
                colorHigh = [0 1 0];
                
                if ~ui.cboxShowLong.Value
                    colorHigh = colorHigh.*-1;
                end
                
                if ~ui.cboxShowShort.Value
                    colorLow = colorLow.*-1;
                end
                
                trackColors = ones(nTracks,3);
                trackColors(longIdx,:) = repmat(colorHigh,nLong,1);
                trackColors(shortIdx,:) = repmat(colorLow,nShort,1);
                plotSettings.trackColors = trackColors;
            case 5
                %% Color-coded mean jump distance
                values = batch(curMovieIndex).results.meanJumpDists;
                
                if isempty(values)
                    plotSettings.trackColors = 0;
                    return
                end
                
                
                trackColors = interp1(linspace(min(values),max(values),length(colorMap)),colorMap,values); % map color to y values
                trackColors = uint8(trackColors*255); % need a 4xN uint8 array
                plotSettings.trackColors = trackColors;
            case 6
                %% Colored by mean jump distance regimes
                
                ui.panelJumpDistRegimes.Visible = 'on';
                
                distForHigh = str2double(ui.editDistForHigh.String);
                meanJumpDist = batch(curMovieIndex).results.meanJumpDists;
                nTracks = batch(curMovieIndex).results.nTracks;
                
                highIdx = meanJumpDist >= distForHigh;
                
                lowIdx = ~highIdx;
                nHigh = sum(highIdx);
                nLow = sum(lowIdx);
                
                colorLow = [1 0 0];
                colorHigh = [0 1 0];
                
                if ~ui.cboxShowHighDist.Value
                    colorHigh = colorHigh.*-1;
                end
                
                if ~ui.cboxShowLowDist.Value
                    colorLow = colorLow.*-1;
                end
                
                
                
                trackColors = ones(nTracks,3);
                trackColors(highIdx,:) = repmat(colorHigh,nHigh,1);
                trackColors(lowIdx,:) = repmat(colorLow,nLow,1);
                plotSettings.trackColors = trackColors;
                
                
            case 7
                %% Color-coded subregions
                subROI = batch(curMovieIndex).subROI;
                nSubROIs = length(subROI);
                colors = distinguishable_colors(nSubROIs+1, {'k'});
                
                if size(colors,1) <= 4 %Easy workaround to brighten or darken the most used colors
                    %                     colors = min(1, colors+0.5)
                    colors = max(0, colors-0.2);
                else
                    colors(1:4,:) = min(1, colors(1:4,:)+0.5);
                end
                
                nTracks = batch(curMovieIndex).results.nTracks;
                trackColors = repmat(colors(1,:),nTracks,1);
                subRegionAssignment = batch(curMovieIndex).results.tracksSubRoi;
                
                for regionIdx = 1:nSubROIs
                    trackIdx = subRegionAssignment == regionIdx;
                    nTracksLongerMinLengthInCurRegion = sum(trackIdx);
                    trackColors(subRegionAssignment == regionIdx,:) = repmat(colors(regionIdx+1,:),nTracksLongerMinLengthInCurRegion,1);
                end
                
                plotSettings.trackColors = trackColors;
        end
        
    end

    function FrameSlider()
        
                
        %% Save field variables in local variables to speed up plotting
        scalingFactor = plotSettings.scalingFactor;
        trackMarkerSize = plotSettings.trackMarkerSize;
        trackLinewidth = plotSettings.trackLinewidth;
        trackColors= plotSettings.trackColors;
        
        %% ----------Show Image--------------------------------------------
        
        %Get current frame from the slider
        curFrame = round(ui.sliderFrame.Value);
        
        %Check wether "Show stack 2" checkbox is enabled
        if ~ui.cboxShowStack2.Value
            %Show stack 1
            curMovieFrame = curFrame;
            %Get current image from image stack
            I = movieStacks{1}(:,:,min(size(movieStacks{1},3),curMovieFrame));
            %Update the current frame number
            ui.textFrame.String = [num2str(curFrame) '/' num2str(batch(curMovieIndex).movieInfo.frames)];
        else
            %Show stack 2
            
            %Frames of stack 2 are distributed equally among frames of stack 1
            divisor = max([1, size(movieStacks{1},3)/size(movieStacks{3},3)]);
            
            %Calculate frame in stack 2 corresponding to the selected frame in stack 1
            curMovieFrame = ceil(curFrame/divisor);
            
            %Get frame from image stack
            I = movieStacks{3}(:,:,min(size(movieStacks{3},3),curMovieFrame));
            
            %Update the current frame number
            ui.textFrame.String = sprintf('Stack 1: \t %d/%d \t Stack 2: %d/%d',...
                curFrame,batch(curMovieIndex).movieInfo.frames,curMovieFrame,size(movieStacks{3},3));
        end
        
        %Check if the checkbox to show a z-projection (eg. standard
        %deviation, detection map etc.) is checked
        if ui.cboxShowZProjection.Value
            if ui.popShowImage.Value ~= 5
                %Show colorbar if any field but "no image" is selected
                colorbar
            else
                colorbar('off')
            end
            
            if ui.popShowImage.Value == 4
                %Show wavelet filtered image
                I = wavelet_filter(I);
            else
                %Show z-projection created in the CreateZProjection function
                I = movieStacks{4};
            end
        else
            colorbar('off')
        end
        
        %Write image to the image handle
        imageHandle.CData = I;
        
        %% ----------Plot tracking-ROI and sub-regions---------------------
        
        hold(ui.axes1,'on')
        
        %Get ROI from batch structure
        trackingROI = batch(curMovieIndex).ROI;
        
        if ~isempty(trackingROI)
            %Multiply ROI coordinates with the currently set scaling factor
            %and write it into the roi handle
            set(ROIHandle,'xdata',trackingROI{1}(:,1).*scalingFactor,'ydata', trackingROI{1}(:,2).*scalingFactor,'Color',plotSettings.roiColor,'Linewidth',plotSettings.roiLinewidth);
        else
            %No Roi exists
            set(ROIHandle,'xdata',NaN,'ydata', NaN);
        end
        
        %Get sub-regions from batch structure
        subROIs = batch(curMovieIndex).subROI;
        
        subRoiHandleNr = 1;
        
        %Get number of existing handles to plotted sub-regions
        nSubRoiHandles = numel(subROIHandle);
        
        %Get sub-region colors
        subRoiColors = plotSettings.subRoiColors;
        
        %Iterate through all subregions
        for nROIs = 1:size(subROIs,2)
            
            if length(subROIs{nROIs}) > 1
                %Sub-region was drawn via thrshold so get the sub-region
                %of the current frame
                subRegionFrame = curFrame;
            else
                %Sub-region was hand-drawn so only one "frame" exists
                subRegionFrame = 1;
            end
            
            %Iterate through all parts of the sub-region in the current frame
            for subRegionPartIdx=1:length(subROIs{nROIs}{subRegionFrame})
                %Get the current part and multiply coordinates with the
                %scaling factor
                curPart = subROIs{nROIs}{subRegionFrame}{subRegionPartIdx}.*scalingFactor;
                if subRoiHandleNr>nSubRoiHandles
                    %Create new plot element
                    subROIHandle(subRoiHandleNr) = plot(ui.axes1,curPart(:,1), curPart(:,2), '-','Color',subRoiColors(mod(nROIs-1,7)+1,:),'Linewidth',4);
                else
                    %User handle to existing plot element and change the
                    %coordinates (faster)
                    set(subROIHandle(subRoiHandleNr),'xdata',curPart(:,1),'ydata', curPart(:,2),'Color',subRoiColors(mod(nROIs-1,7)+1,:),'Linewidth',4);
                end
                subRoiHandleNr = subRoiHandleNr +1 ;
            end
        end
        
        %Delete unused handles
        for q = subRoiHandleNr:numel(subROIHandle)
            delete(subROIHandle(subRoiHandleNr));
            subROIHandle(subRoiHandleNr) = [];
        end
        
        %% ----------Plot tracks-------------------------------------------
        
        %Get tracks from the batch structure
        tracks = batch(curMovieIndex).results.tracks;
        handleNr = 1;
        
        %Get number of existing handles to plotted tracks
        nTrackHandles = numel(trackHandles);
        
        %Plot tracks
        if ~isempty(tracks) && ui.cboxTracks.Value
            
            switch ui.popTracks.Value
                case 1
                    %% Show tracks in range
                    
                    %Matrix of logicals indicating which track has to be shown in each frame
                    tracksInCurFrame = plotSettings.tracksInFrame{curFrame};
                    
                    %Preallocate the array of graphic object handles
                    trackHandles(end+1:numel(tracksInCurFrame)) = gobjects(numel(tracksInCurFrame)-numel(trackHandles),1);
                    
                    %Get number of different track colors
                    uniqueTrackColors = unique(trackColors,'rows');
                    nColors = size(uniqueTrackColors,1);
                    
                    %Iterate through track colors
                    for curSpotColor = 1:nColors
                        
                        if any(uniqueTrackColors(curSpotColor,:) < 0)
                            %Color has to be skipped (if user chose "color coded
                            %by track length regime" or "color coded by mean jump distance regime"
                            %and wants to see eg. only long tracks)
                            continue
                        end
                        
                        %Find tracks that match the current color
                        curColorTracks = find(all(trackColors == uniqueTrackColors(curSpotColor,:),2));
                        
                        for trackID = tracksInCurFrame(ismember(tracksInCurFrame,curColorTracks))
                            %Plot tracks only until current frame
                            maskToPlot = tracks{trackID}(:,1)<=curFrame;
                            
                            if handleNr>nTrackHandles
                                %Create new track
                                trackHandles(handleNr) = plot(ui.axes1,tracks{trackID}(maskToPlot, 2).*scalingFactor, tracks{trackID}(maskToPlot, 3).*scalingFactor, '.-','Color',trackColors(trackID,:),'Linewidth',trackLinewidth,'MarkerSize',trackMarkerSize,'UserData',trackID);
                            else
                                %Use existing track and just change its data/properties
                                set(trackHandles(handleNr),'xdata',tracks{trackID}(maskToPlot, 2).*scalingFactor,'ydata', tracks{trackID}(maskToPlot, 3).*scalingFactor,'Color',trackColors(trackID,:),'Linewidth',trackLinewidth,'MarkerSize',trackMarkerSize,'UserData',trackID);
                            end
                            handleNr=handleNr+1;
                        end
                    end
                case 2
                    %% Show all Tracks
                    if plotSettings.tracksChanged
                        %Tracks are plotted only if track
                        %settings have been changed (No need to plot tracks
                        %in each frame)
                        
                        %Preallocate the array of graphic object handles
                        trackHandles(end+1:length(tracks)) = gobjects(length(tracks)-numel(trackHandles),1);
                        
                        %Get number of different track colors
                        uniqueTrackColors = unique(trackColors,'rows');
                        nColors = size(uniqueTrackColors,1);
                        
                        for curSpotColor = 1:nColors
                            
                            if any(uniqueTrackColors(curSpotColor,:) < 0)
                                %Color has to be skipped (if user chose "color coded
                                %by track length regime" or "color coded by mean jump distance regime"
                                %and wants to see eg. only long tracks)
                                continue
                            end
                            
                            %Create logical array indicating the tracks
                            %that have to be plotted in the current color
                            tracksWithCurColor = all(trackColors == uniqueTrackColors(curSpotColor,:),2);
                            
                            for trackID = find(tracksWithCurColor)'
                                if handleNr>nTrackHandles
                                    %Create new track
                                    trackHandles(handleNr) = plot(ui.axes1,tracks{trackID}(:, 2).*scalingFactor, tracks{trackID}(:, 3).*scalingFactor, '.-','Color',trackColors(trackID,:),'Linewidth',trackLinewidth,'MarkerSize',trackMarkerSize,'UserData',trackID);
                                else
                                    %Use existing track and just change its data/properties
                                    set(trackHandles(handleNr),'xdata',tracks{trackID}(:, 2).*scalingFactor,'ydata', tracks{trackID}(:, 3).*scalingFactor,'Color',trackColors(trackID,:),'Linewidth',trackLinewidth,'MarkerSize',trackMarkerSize,'UserData',trackID);
                                end
                                
                                handleNr=handleNr+1;
                            end
                        end
                        
                        plotSettings.tracksChanged = 0;
                    else
                        handleNr = numel(trackHandles)+1;
                    end
                case 3
                    %% Show initial positions
                    
                    %Iterate through all tracks and get the first position
                    %of each track
                    curTrackPos = zeros(length(tracks),2);
                    
                    for trackId =1 :length(tracks)
                        curTrackPos(trackId,:) = tracks{trackId}(1, 2:3).*scalingFactor;
                    end
                    
                    
                    set(trackInitialPosHandle,'xdata',curTrackPos(:, 1),'ydata',...
                        curTrackPos(:, 2),'Marker',plotSettings.initialPosMarker,'LineWidth',trackLinewidth,'SizeData',trackMarkerSize^2,...
                        'MarkerFaceColor','none','MarkerEdgeColor','flat','CData', trackColors);
                    case 4
                    %% Show all tracked positions
                    
                    trackPos = vertcat(tracks{:});
                    trackPos = trackPos(:,2:3).*scalingFactor;
                    
                    trackCol = cell(length(tracks),1);
                    
                    for trackId =1:length(tracks)
                        nSpotsCurTrack = size(tracks{trackId},1);
                        trackCol{trackId} = repmat(trackColors(trackId, :),nSpotsCurTrack,1);
                    end
                    
                    trackCol = vertcat(trackCol{:});

                    set(trackInitialPosHandle,'xdata',trackPos(:, 1),'ydata',...
                        trackPos(:, 2),'Marker',plotSettings.initialPosMarker,'LineWidth',trackLinewidth,'SizeData',trackMarkerSize^2,'CData',trackCol);
            end
        else
            plotSettings.tracksChanged = 1;
        end
        
        %Delete initial positions if they are not needed
        if ui.popTracks.Value < 3 || ~ui.cboxTracks.Value || isempty(tracks)
            trackInitialPosHandle.XData = NaN;
            trackInitialPosHandle.YData = NaN;
            trackInitialPosHandle.CData = NaN;
        end
        
        % Delete unsused track handles
        for trackSpotIdx = handleNr:numel(trackHandles)
            delete(trackHandles(handleNr));
            trackHandles(handleNr) = [];
        end
        
        %% ----------Plot non-linked spots---------------------------------
        
        if ~isempty(batch(curMovieIndex).results.nonLinkedSpots) && ui.cboxShowSingle.Value
            %Get nonLinkedSpots from batch structure and multiply the
            %coordinates with the scaling factor
            nonLinkedSpots = batch(curMovieIndex).results.nonLinkedSpots.*scalingFactor;
            
            %Save coordinates and display options in the existing handle
            spotsNonLinkedHandle.XData = nonLinkedSpots(:,2);
            spotsNonLinkedHandle.YData = nonLinkedSpots(:,3);
            spotsNonLinkedHandle.MarkerSize = plotSettings.singleMarkerSize;
            spotsNonLinkedHandle.MarkerEdgeColor = plotSettings.singleColor;
            spotsNonLinkedHandle.Marker = plotSettings.singleMarker;
        else
            %Hide nonLinkedSpots
            spotsNonLinkedHandle.XData = NaN;
            spotsNonLinkedHandle.YData = NaN;
        end
        
        %% ----------Plot Spots--------------------------------------------
        
        if ~isempty(batch(curMovieIndex).results.spotsAll) && ui.cboxSpots.Value
            spotsHandle.XData = batch(curMovieIndex).results.spotsAll{curFrame}(:,1).*scalingFactor;
            spotsHandle.YData = batch(curMovieIndex).results.spotsAll{curFrame}(:,2).*scalingFactor;
            spotsHandle.MarkerSize = plotSettings.spotMarkerSize;
            spotsHandle.MarkerEdgeColor = plotSettings.spotColor;
            spotsHandle.Marker = plotSettings.spotMarker;
            uistack(spotsHandle,'top')
        else
            spotsHandle.XData = NaN;
            spotsHandle.YData = NaN;
        end
        
        %Put scalebar on top if scalebar dimension is > 0
        if sum(plotSettings.sbDim > 0)
            uistack(scalebarHandle,'top')
        end
        
        hold(ui.axes1,'off')
        
        %Call function to automatically adjust brightness if "Always"
        %checkbox is checked
        if ui.cboxContAutoAdj.Value
            AdjustBrightnessCB('autoAdjust')
        end
    end

%-------------------Menubar callbacks--------------------------------------

%Callbacks from file menu

    function LoadBatchFileCB(~,~)
        %User pressed File -> Load batch file
        
        %% Load batch file
        
        %Open file selection dialog
        [fileName,pathName] = uigetfile('*.mat','Select .mat batch file',searchPath{2},'MultiSelect','off');
        
        if isequal(fileName,0)
            %User didn't choose a file
            return
        end
        
        ui.editFeedbackWin.String = char('Loading batch file');
        drawnow
        
        %Load batch form .mat file
        fullBatchFilePath = fullfile(pathName,fileName);
        loadedBatchFile = load(fullBatchFilePath);
        
        %Check if .mat file contains batch variable
        if ~isfield(loadedBatchFile,'batch')
            ui.editFeedbackWin.String = char('Please choose a valid batch file');
            return
        end
        
        %Get batch and filesTable from the loaded mat file
        loadedBatch = loadedBatchFile.batch;
        loadedFilesTable = loadedBatchFile.filesTable;
        
        %Get number of movies in batch
        nFiles = length(loadedBatch);
        
        %% Check if all movies are found
        missingMoviesIndices = [];
        
        %Go through all files and check if they exist
        for n = 1:nFiles
            %In earlier TackIt versions file and pathnames have been stored as strings which
            %sometimes makes trouble so make sure that they are character arrays
            if ~ischar(loadedBatch(n).movieInfo.fileName) || ~ischar(loadedBatch(n).movieInfo.pathName)
                loadedBatch(n).movieInfo.fileName = char(loadedBatch(n).movieInfo.fileName);
                loadedBatch(n).movieInfo.pathName = char(loadedBatch(n).movieInfo.pathName);
            end
            
            if ~ischar(loadedBatch(n).movieInfo.fileName2) || ~ischar(loadedBatch(n).movieInfo.pathName2)
                loadedBatch(n).movieInfo.fileName2 = char(loadedBatch(n).movieInfo.fileName2);
                loadedBatch(n).movieInfo.pathName2 = char(loadedBatch(n).movieInfo.pathName2);
            end
            
            %Check if current movie exists
            curMovie = fullfile(loadedBatch(n).movieInfo.pathName,loadedBatch(n).movieInfo.fileName);
            if exist(curMovie, 'file') ~= 2
                %File not found so save index of current file in the
                %missingMoviesIndices array
                missingMoviesIndices = [missingMoviesIndices n];
            end
        end
        
        %Get number of missing movies
        nMissing = numel(missingMoviesIndices);
        
        %Ask the user to search for missing movies until either all movies
        %have been found or the user aborts
        while nMissing
            %Create message for the user listing the first 10 movies which have not
            %been found
            message = [num2str(nMissing) ' movies have not been found:' newline];
            for n = missingMoviesIndices
                curFileName = loadedBatch(n).movieInfo.fileName;
                message = [message curFileName newline];
                if n > 10
                    message = [message '...' newline];
                    break
                end
            end
            
            message = [message newline 'Do you want to choose a folder where the files should be searched?'];
            
            %Opden question dialog to ask if the user wants to search for
            %the missing files
            answer = questdlg(message,'Movies not found','Yes','No','Cancel','Yes');
            
            if isempty(answer) || strcmp(answer,'Cancel')
                %Usre pressed "cancel"
                return
            elseif strcmp(answer,'No')
                %User pressed "No'"
                break
            end
            
            %Open folder dialog box folder to ask user where movies should be searched
            newPathName = uigetdir(searchPath{2},'Choose a parent folder containing the movies');
            
            if newPathName == 0
                %User did not choose a folder
                return
            end
            
            curIndex = 1;
            m = 1;
            
            %Create list of all files with .tif or .tiff extension
            fileList = dir(fullfile(newPathName,'**\*.tif*'));
            fileNameList = {fileList.name};
            pathNameList = {fileList.folder};
            
            %Go through all missing movies and compare each with the list of files in the folder specified by the user
            while m <= numel(missingMoviesIndices)
                
                %Get movie index from the array containing indices of missing movies
                curMovieNum = missingMoviesIndices(m);
                
                %Display search progress in feedback window
                ui.editFeedbackWin.String = char(...
                    ['Searching Movie ' num2str(curIndex) ' of ' num2str(nMissing), newline,...
                    'Movies found: ' num2str(nMissing-numel(missingMoviesIndices)), newline]);
                
                curIndex = curIndex + 1;
                drawnow
                
                
                %Compare the the list of .tif files with the current filename
                Index = find(cellfun(@(s) strcmp(loadedBatch(curMovieNum).movieInfo.fileName, s),fileNameList));
                
                if Index
                    %Movie was found in the folder specified by the user
                    if numel(Index) > 1
                        %More than one file was found so warn the user
                        ui.editFeedbackWin.String = char(...
                            ['Warning: more than one file with the name  ' loadedBatch(curMovieNum).movieInfo.fileName ' was found. Using ' fullfile(pathNameList{Index(1)},fileNameList{Index(1)})]);
                    end
                    
                    %Save the new filepath into the batch structure and the files table
                    loadedBatch(curMovieNum).movieInfo.fileName = fileNameList{Index(1)};
                    loadedBatch(curMovieNum).movieInfo.pathName = pathNameList{Index(1)};
                    loadedBatch(curMovieNum).movieInfo.pathName2 = pathNameList{Index(1)};
                    loadedFilesTable.FileName(curMovieNum) = fileNameList(Index(1));
                    loadedFilesTable.PathName(curMovieNum) = pathNameList(Index(1));
                    
                    %Movie was found in new path so delete it from the list
                    %of missing movies
                    missingMoviesIndices(m) = [];
                else
                    %Movie was not found, so go to the next movie
                    m = m+1;
                end
            end
            
            ui.editFeedbackWin.String = char(['Search finished: ', num2str(nMissing-numel(missingMoviesIndices)),' of ',num2str(nMissing), ' Movies have been found']);
            %Get number of remaining movies that are missing
            nMissing = numel(missingMoviesIndices);
        end
        
        %Save loaded batch structure and loaded filesTable
        batch = loadedBatch;
        filesTable = loadedFilesTable;
        
        %% Update UI
        
        %Amount of timelapse conditions might have changed so uncheck
        %"timelpase specific" checkbox and update tl specific tracking
        %parameter list
        ui.cboxTlDependentTr.Value = 0;
        TrackingParamsCB()
        
        %Update list of movies in the popup menu
        ui.popMovieSelection.String = filesTable.FileName;
        
        %Set movienumber of current movie to first movie and update ui
        %accordingly
        curMovieIndex = 1;
        ui.popMovieSelection.Value = 1;
        ui.editMovie.String = 1;
        %Display amount of movies in ui
        ui.textMovie2.String = ['/' num2str(length(batch))];
        
        %Set movie search path to first movie in batch
        searchPath{1} = batch(1).movieInfo.pathName;
        %Set batch search path to the folder of the selected batch file
        searchPath{2} = pathName;
        %Set batch file patch to the selected filepath
        searchPath{3} = fullBatchFilePath;
        
        %Adjust ui to new movie
        AdjustUiToNewMovie()
        
    end

    function SaveBatchCB(~,~)
        %User pressed File -> Save batch file as...
        
        if isempty(searchPath{3})
            %No batch file has been loaded before so suggest a new filename
            startingPath = fullfile(searchPath{2},[datestr(now,'yymmdd'), '_trackit_batch', '.mat']);
        else
            %Batch file has been loaded before so suggest last opened filename
            startingPath = searchPath{3};
        end
        
        %Open file selection dialog
        [fileName, pathName] = uiputfile('*.mat','Choose filename for saving batch job in .mat file',startingPath);
        
        if fileName == 0
            %User did not choose a file
            return
        end
        
        %Set search path for batch files to current path
        searchPath{2} = pathName;
        %Save filepath of current patch file
        searchPath{3} = fullfile(pathName, fileName);
        
        ui.editFeedbackWin.String = 'Saving batch file, please wait...';
        drawnow
        
        %Create structure containing filesTable and batch struct
        S.filesTable = filesTable;
        S.batch = batch;
        
        %Save .mat batch file
        save(fullfile(pathName,fileName),'-struct','S')
        ui.editFeedbackWin.String = 'Saving finished';
    end

    function MergeBatchFilesCB(~,~)
        %User pressed File -> Merge multiple batch files...
        
        %Open batch file merger
        batch_file_merger(searchPath{2})
    end

    function MultiBatchAnalyzerCB(~,~)
        %User pressed File -> Re-analyze multiple batch files
        multiple_batches_analyzer(searchPath{1});
    end

    function SplitBatchCB(~,~)
        %User pressed File -> change filenames in current batch
               
        %Display user dialog
        prompt = {['You can use this tool to divide your batch file into '...
            'smaller batches where each batch file contains a subset of your movies. '...
            'Enter a regular expression to search for a '...
            'pattern in your movie filenames. Example: \d{6} represents a '...
            'sequence of 6 numeric digits and you can use it if you want '...
            'one batch file per measurement day. Or use _c\d+ if you have '...
            'filenames containing _c01, c_02 etc. '...
            'See https://de.mathworks.com/help/matlab/matlab_prog/regular-expressions.html.']};
        dlgtitle = 'Divide batch file into smaller batches';
        answer = inputdlg(prompt,dlgtitle);
        
        if isempty(answer)
            %Canceled
            return
        end
        
        %Open folder selection dialog
        parentFolder = uigetdir(searchPath{2},'Choose a folder where the batch files should be saved');
        
        if parentFolder == 0
            return
        end
        
        searchString = answer{1};
        
        %Get cell array of filenames
        fileNames = filesTable.FileName;
        
        %Convert fileNames into a  string array
        if iscell(fileNames)
            newFileNameList = "";
            for idx = 1:length(fileNames)
                
                curName = fileNames{idx};
                
                if ischar(curName)
                    curName = string(curName);
                end
                
                newFileNameList(idx,1) = curName;
                
            end  
            fileNames = newFileNameList;
        end
        
        %Find filenames which contain the string pattern and return the
        %matching part of the filenames
        [matchResults, ~] = regexp(fileNames, searchString,'match','split');
        
        %Make sure that no movie has more than one match
        for idx = 1:length(matchResults)
           if numel(matchResults{idx}) > 1               
              msgbox(['The name of movie ', num2str(idx), ' was matched more than once to the provided regular expression pattern: ', sprintf('%s, ',matchResults{idx})]) 
              return              
           end            
        end
        
        %Catenate the cell array to create an array of strings
        matchResults = vertcat(matchResults{:});
                
        %Create a unique list of matches
        uniqueMatches = unique(matchResults);
        
        %Iterate through the list of unique matches
        for stageIdx = 1:size(uniqueMatches,1)
            %Get current matched string pattern
            curMatch = uniqueMatches(stageIdx,1);
            
            %Get list of logicals matching the current string
            movieIndices = strcmp(matchResults, curMatch);
            
            %Create batch file and files table with subset of movies
            curBatch = batch(movieIndices);
            curFilesTable = filesTable(movieIndices,:);
            
            %Save new batch to file
            S.filesTable = curFilesTable;
            S.batch = curBatch;            
            save(fullfile(parentFolder,strcat(curMatch,'.mat')),'-struct','S')
        end
        
        msgbox('Batch file splitting finished.')
        
    end

    function ChangeFilenamesCB(~,~)
        %User pressed Stuff -> Change filenames in current batch
        
        %Display user dialog
        prompt = {'Search for:','Replace with:'};
        dlgtitle = 'Replace filenames';
        dims = [1 40];
        answer = inputdlg(prompt,dlgtitle,dims);
        
        if isempty(answer)
            return
        end
        
        %Get number of movies
        nFiles = length(batch);
        
        %Iterate through movies
        for fileIdx = 1:nFiles
            %Get filename of current movie
            curFilename = batch(fileIdx).movieInfo.fileName;
            
            %Check if a part of the filename matches the regular expression pattern
            [matchedString, splitResult] = regexp(curFilename, answer{1},'match','split');
                        
            if ~isempty(matchedString)
                
                %Create new filename
                newFilename = '';
                
                for idx = 1:length(splitResult)-1                    
                    newFilename = strcat(newFilename, splitResult{idx}, answer{2});                    
                end
                
                newFilename = char([newFilename , splitResult{end}]);
                
                %Save new filename
                batch(fileIdx).movieInfo.fileName = newFilename;
                
                %Make sure that filename is saved as cell (compatibility
                %with earlier TrackIt versions)
                if iscell(filesTable.FileName(fileIdx))
                    filesTable.FileName{fileIdx} = newFilename;
                else
                    filesTable.FileName(fileIdx) = newFilename;
                end
                
            end
            
        end
        
        %Question dialog to save batch file
        answer = questdlg('Would you like to save the batch file?', ...
            'Save batch file?', ...
            'Yes','No thank you','Yes');
        
        % Handle response
        switch answer
            case 'Yes'
                SaveBatchCB()
        end
        
        
    end

    function CopyToWorkspaceCB(~,~)
        %User pressed File -> Export all data to Matlab workspace
        
        assignin('base','trackitBatch',batch);
        assignin('base','trackitMoviestacks',movieStacks);
        assignin('base','trackitFilesTable',filesTable);
        
    end

    function ExportTracksCB()
        %User pressed File -> Export tracks to .mat or to .csv
        
        %Open export dialog window
        out = export_dialog(pixelSize,searchPath{2});
        
        if isempty(out.destination)
            %User canceled
            return
        end
        
        ui.editFeedbackWin.String = 'Exporting tracks...';
        drawnow
        
        %% Get timelapse conditions and amount of movies in each timelapse
        
        if out.filesPerMovie == 2
            %Export one file per cycle time
            
            %Get frame cycle time of each movie
            frameCycleTimeMovieList = zeros(1,length(batch));
            for movieIdx = 1:length(batch)
                frameCycleTimeMovieList(movieIdx) = batch(movieIdx).movieInfo.frameCycleTime;
            end
            
            %Create unique list of frame cycle times
            frameCycleTimesList = unique(frameCycleTimeMovieList);
            
            %Initialize an array where each entry represents a counter for
            %the amount of movies in each tl condition. The movie counter of the
            %corresponging tl condition is increased when iterating through
            %the movies.
            movieCounter = ones(1,numel(frameCycleTimesList));
        end
        
        
        %% Create export data
        
        %Get pixelsize
        pixelSize = out.pixelSize;
        
        %Get information on what parameters user wants to export
        additionalExport = cell2mat(out.additionalExport(:,1));
        
        %Initialize
        movieInfoAllMovies = struct;
        trackedParAllMovies = cell(1,length(batch));
        subRegionAssignment = cell(1,length(batch));
        
        %Iterate through all movies of the current batch
        for movieIdx = 1:length(batch)
            if additionalExport(1)
                %Export settings (movie info and tracking parameters)
                movieInfoAllMovies(movieIdx).movieInfo = batch(movieIdx).movieInfo;
                movieInfoAllMovies(movieIdx).movieInfo.pixelSize = pixelSize;
                movieInfoAllMovies(movieIdx).trackingParams = batch(movieIdx).params;
            end
            if additionalExport(2)
                %Export ROI
                movieInfoAllMovies(movieIdx).roi = batch(movieIdx).ROI;
            end
            if additionalExport(3)
                %Export ROI size
                movieInfoAllMovies(movieIdx).roiSize = batch(movieIdx).results.roiSize;
            end
            if additionalExport(4)
                %Export track lengths
                movieInfoAllMovies(movieIdx).tracklengths = batch(movieIdx).results.trackLengths;
            end
            if additionalExport(5)
                %Export jump distances
                movieInfoAllMovies(movieIdx).jumpDistances = batch(movieIdx).results.jumpDistances;
            end
            if additionalExport(6)
                %Export jump angles
                movieInfoAllMovies(movieIdx).angles = batch(movieIdx).results.angles;
            end
            if additionalExport(7)
                %Export number of detections
                movieInfoAllMovies(movieIdx).nSpots = batch(movieIdx).results.nSpots;
            end
            if additionalExport(8)
                %Export number of tracks
                movieInfoAllMovies(movieIdx).nTracks = batch(movieIdx).results.nTracks;
            end
            if additionalExport(9)
                %Export number of non-linked spots
                movieInfoAllMovies(movieIdx).nNonLinkedSpots = batch(movieIdx).results.nNonLinkedSpots;
            end
            if numel(additionalExport) > 9 && additionalExport(10)
                %Export sub-region results
                movieInfoAllMovies(movieIdx).subRegionResults = batch(movieIdx).results.subRegionResults;
            end
            
            %Get tracks of current movie
            curTracks = batch(movieIdx).results.tracks;
            %Initialize trackedPar struct array
            trackedPar = repmat(struct('xy',[],'Frame',[],'TimeStamp',[],'Movie',[]),1,length(curTracks));
            
            %Iterate through all tracks
            for trackIdx = 1:length(curTracks)
                %Get coordinates of current track
                trackedPar(trackIdx).xy = curTracks{trackIdx}(:,2:3).*pixelSize;
                
                %Get frames of current track
                curTrackFrames = curTracks{trackIdx}(:,1);
                
                %Save frames of current track in trackedPar struct
                trackedPar(trackIdx).Frame = curTrackFrames;
                
                if out.periodicIllumination
                    %Illumination pattern is set to "Continuous"
                    
                    %Get frame cycle time of current movie
                    cycleTime = batch(movieIdx).movieInfo.frameCycleTime;
                    
                    %Convert to seconds and save in trackedPar struct
                    trackedPar(trackIdx).TimeStamp = curTrackFrames*cycleTime/1000;
                else
                    %Illumination pattern is set to "ITM"
                    
                    %Get number of subsequent bright frames
                    nBrightFrames = out.nBrightFrames;
                    %Get exposure time
                    expTime = out.itmExpTime;
                    %Get dark time
                    darkTime = out.itmDarkTime;
                    
                    %Calculate frame times of current track in seconds and
                    %save in trackedPar struct
                    trackedPar(trackIdx).TimeStamp = curTrackFrames*expTime/1000+...
                        floor((curTrackFrames-1)/nBrightFrames)*darkTime/1000;
                end
                
                
            end
            
            if out.filesPerMovie == 1 || out.filesPerMovie == 3 || out.filesPerMovie == 4
                %One file per movie, all tracks in one file or one file per sub-region
                
                %Save movie nummber of all tracks in this movie in the
                %trackedPar struct.
                movieNum = repmat({movieIdx},1,length(curTracks));
                [trackedPar.Movie] = movieNum{:};
            elseif out.filesPerMovie == 2
                %One file per cycle time
                
                %Get tl condition of current movie
                curTlCond = batch(movieIdx).movieInfo.frameCycleTime;
                
                %Save movie nummber of all tracks in this movie in the
                %trackedPar struct.
                movieNum = repmat({movieCounter(curTlCond == frameCycleTimesList)},1,length(curTracks));
                [trackedPar.Movie] = movieNum{:};
                
                %Increase count of the number of movies with this tl
                %conditions by one
                movieCounter(curTlCond == frameCycleTimesList) = movieCounter(curTlCond == frameCycleTimesList) + 1;
            end
            
            %Save the array containing information to which sub-region the
            %tracks in the current movie belong
            subRegionAssignment{movieIdx} = batch(movieIdx).results.tracksSubRoi';
            %Save trackPar struct of this movie in a cell array for all
            %movies
            trackedParAllMovies{movieIdx} = trackedPar;
        end
        
        
        %% Write to .mat or csv
        
        if out.fileFormatMat
            %% Export as .mat
            
            if out.filesPerMovie == 1
                %% One file per movie
                
                %Iterate through all movies of batch
                for movieIdx = 1:length(batch)
                    %Get trackedPar struct for current movie
                    trackedPar = trackedParAllMovies{movieIdx};
                    
                    if out.combatibility == 2
                        %Convert for vbtSPT combatibility
                        trackedPar = struct2cell(trackedPar);
                        trackedPar = permute(trackedPar,[3 1 2]);
                        trackedPar = trackedPar(:,1);
                    end
                    
                    %Get filename of current movie and create saving filename
                    [~,fileName,~] = fileparts(batch(movieIdx).movieInfo.fileName);
                    fullFileName = fullfile(out.destination,strcat(fileName, '_trackedPar.mat'));
                    if sum(additionalExport)
                        %Also save additional information requested by the user
                        %such as ROI or jump distances
                        moviewiseData = movieInfoAllMovies(movieIdx);
                        save(fullFileName,'-mat','trackedPar','moviewiseData')
                    else
                        %Save only tracks
                        save(fullFileName,'-mat','trackedPar')
                    end
                end
            elseif out.filesPerMovie == 2
                %% One file per cycle time
                
                %Iterate through all timelapse conditions
                for curTlCond = frameCycleTimesList
                    %Get all tracks of movies with current tl condition
                    trackedPar = [trackedParAllMovies{frameCycleTimeMovieList == curTlCond}];
                    if out.combatibility == 2
                        %Covnert for vbtSPT combatibility
                        trackedPar = struct2cell(trackedPar);
                        trackedPar = permute(trackedPar,[3 1 2]);
                        trackedPar = trackedPar(:,1);
                    end
                    
                    %Create saving filename
                    fullFileName = fullfile(out.destination,strcat('trackedPar_',num2str(curTlCond),'ms.mat'));
                    
                    if sum(additionalExport)
                        %Also save additional information requested by the user
                        %such as ROI or jump distances
                        moviewiseData = movieInfoAllMovies(frameCycleTimeMovieList == curTlCond);
                        save(fullFileName,'-mat','trackedPar','moviewiseData')
                    else
                        %Save only tracks
                        save(fullFileName,'-mat','trackedPar')
                    end
                end
            elseif out.filesPerMovie == 3
                %% All tracks in one file
                
                %Catenate all tracks of all movies
                trackedPar = [trackedParAllMovies{:}];
                if out.combatibility == 2
                    %Convert for vbtSPT combatibility
                    trackedPar = struct2cell(trackedPar);
                    trackedPar = permute(trackedPar,[3 1 2]);
                    trackedPar = trackedPar(:,1);
                end
                if sum(additionalExport)
                    %Also save additional information requested by the user
                    %such as ROI or jump distances
                    moviewiseData = movieInfoAllMovies;
                    save(out.destination,'-mat','trackedPar','moviewiseData')
                else
                    %Save only tracks
                    save(out.destination,'-mat','trackedPar')
                end
            elseif out.filesPerMovie == 4
                %% One file per sub-region
                
                %Catenate the sub-region numbers of all tracks of all movies
                subRegionAssignment = [subRegionAssignment{:}];
                
                %Catenate all tracks of all movies
                trackedParAllMovies = [trackedParAllMovies{:}];
                
                %Get number of sub-regions
                nRegions = max(subRegionAssignment);
                
                %Iterate through all regions
                for curSubRegion = 0:nRegions
                    %Get tracks in current sub-region
                    trackedPar = trackedParAllMovies(curSubRegion == subRegionAssignment);
                    if out.combatibility == 2
                        %Convert for vbtSPT combatibility
                        trackedPar = struct2cell(trackedPar);
                        trackedPar = permute(trackedPar,[3 1 2]);
                        trackedPar = trackedPar(:,1);
                    end
                    
                    if curSubRegion == 0
                        %Save tracks in region 0 which corresponds to the tracking region
                        fullFileName = fullfile(out.destination,strcat('trackedPar_trackingRegion.mat'));
                    else
                        %Save tracks in current sub-region
                        fullFileName = fullfile(out.destination,strcat('trackedPar_subRegion_',num2str(curSubRegion),'.mat'));
                    end
                    
                    if sum(additionalExport)
                        %Also save additional information requested by the user
                        %such as ROI or jump distances
                        moviewiseData = movieInfoAllMovies;
                        save(fullFileName,'-mat','trackedPar','moviewiseData')
                    else
                        %Save only tracks
                        save(fullFileName,'-mat','trackedPar')
                    end
                end
            end
            
        else
            %% Export as .csv
            if out.filesPerMovie == 1
                %% One file per movie
                
                %Iterate though all movies of current batch
                for movieIdx = 1:length(batch)
                    %Get fileName of current movie
                    [~,fileName,~] = fileparts(batch(movieIdx).movieInfo.fileName);
                    
                    %Create saving filename
                    fullFileName = fullfile(out.destination,strcat(fileName, '_trackedPar', '.csv'));
                    
                    %Open .csv file
                    fid = fopen(fullFileName,'wt');
                    
                    %Write header into .csv file
                    fprintf(fid, ',frame,t,trajectory,x,y\n');
                    
                    %Get tracks in current movie
                    curMovieTracks = trackedParAllMovies{movieIdx};
                    
                    %Initialize spot counter
                    spotNum = 0;
                    
                    %Iterate through all tracks of the current movie
                    for trackIdx = 1:length(curMovieTracks)
                        %Iterate through all detections of the current track
                        for spotIdx = 1:length(curMovieTracks(trackIdx).Frame)
                            %Print current detection to .csv file
                            fprintf(fid, '%d,%d,%.4f,%d,%f,%f\n',...
                                spotNum, ...
                                curMovieTracks(trackIdx).Frame(spotIdx),...
                                curMovieTracks(trackIdx).TimeStamp(spotIdx),...
                                trackIdx,...
                                curMovieTracks(trackIdx).xy(spotIdx,1),...
                                curMovieTracks(trackIdx).xy(spotIdx,2));
                            
                            %Increase spot count
                            spotNum = spotNum+1;
                        end
                    end
                    
                    %Close .csv file
                    fclose(fid);
                    
                end
            elseif out.filesPerMovie == 2
                %% One file per cycle time
                
                %Iterate through all tl conditions
                for curTlCond = frameCycleTimesList
                    
                    %Create saving filename
                    fullFileName = fullfile(out.destination,strcat('trackedPar_',num2str(curTlCond),'ms.csv'));
                    
                    %Open .cvs file
                    fid = fopen(fullFileName,'W');
                    
                    %Write header to .csv file
                    fprintf(fid, ',frame,t,trajectory,x,y\n');
                    
                    %Initialize track and spot counter
                    trackNum = 1;
                    spotNum = 0;
                    
                    %Iterate through all movies of current tl condition
                    for movieIdx = find(frameCycleTimeMovieList == curTlCond)
                        %Get all tracks in current movie
                        curMovieTracks = trackedParAllMovies{movieIdx};
                        
                        %Iterate through all tracks
                        for trackIdx = 1:length(curMovieTracks)
                            %Iterate through all detections of the current track
                            for spotIdx = 1:length(curMovieTracks(trackIdx).Frame)
                                %Print current detection to .csv file
                                fprintf(fid, '%d,%d,%.4f,%d,%f,%f\n', ...
                                    spotNum, ...
                                    curMovieTracks(trackIdx).Frame(spotIdx),...
                                    curMovieTracks(trackIdx).TimeStamp(spotIdx),...
                                    trackNum,...
                                    curMovieTracks(trackIdx).xy(spotIdx,1),...
                                    curMovieTracks(trackIdx).xy(spotIdx,2));
                                
                                %Increase spot count
                                spotNum = spotNum+1;
                            end
                            %Increase track count
                            trackNum = trackNum + 1;
                        end
                        
                    end
                    
                    %Close .csv file
                    fclose(fid);
                    
                end
            elseif out.filesPerMovie == 3
                %% All tracks in one file
                
                %Open .cvs file
                fid = fopen(out.destination,'W');
                
                %Write header to .csv file
                fprintf(fid, ',frame,t,trajectory,x,y\n');
                
                %Initialize track and spot counter
                trackNum = 0;
                spotNum = 0;
                
                %Iterate through all movies
                for movieIdx = 1:length(batch)
                    %Get all tracks in current movie
                    curMovieTracks = trackedParAllMovies{movieIdx};
                    
                    %Iterate through all tracks
                    for trackIdx = 1:length(curMovieTracks)
                        %Iterate through all detections of the current track
                        for spotIdx = 1:length(curMovieTracks(trackIdx).Frame)
                            %Print current detection to .csv file
                            fprintf(fid, '%d,%d,%.4f,%d,%f,%f\n',...
                                spotNum, ...
                                curMovieTracks(trackIdx).Frame(spotIdx),...
                                curMovieTracks(trackIdx).TimeStamp(spotIdx),...
                                trackNum,...
                                curMovieTracks(trackIdx).xy(spotIdx,1),...
                                curMovieTracks(trackIdx).xy(spotIdx,2));
                            
                            %Increase spot count
                            spotNum = spotNum+1;
                        end
                        %Increase track count
                        trackNum = trackNum + 1;
                    end
                end
                
                %Close .csv file
                fclose(fid);
            elseif out.filesPerMovie == 4
                %% One file per sub-region
                
                %Catenate the sub-region numbers of all tracks of all movies
                subRegionAssignment = [subRegionAssignment{:}];
                
                %Catenate all tracks of all movies
                trackedParAllMovies = [trackedParAllMovies{:}];
                
                %Get number of sub-regions
                nRegions = max(subRegionAssignment);
                
                %Iterate through all regions
                for curSubRegion = 0:nRegions
                    if curSubRegion == 0
                        %Create filename for tracks in region 0 which corresponds to the tracking region
                        fullFileName = fullfile(out.destination,strcat('trackedPar_trackingRegion.csv'));
                    else
                        %Create filename for tracks in current sub-region
                        fullFileName = fullfile(out.destination,strcat('trackedPar_subRegion_',num2str(curSubRegion),'.csv'));
                    end
                    
                    %Open .cvs file
                    fid = fopen(fullFileName,'W');
                    
                    %Write header to .csv file
                    fprintf(fid, ',frame,t,trajectory,x,y\n');
                    
                    %Initialize track and spot counter
                    spotNum = 0;
                    trackNum = 1;
                    
                    %Get all tracks with current sub-region number
                    curSubRegionTracks = trackedParAllMovies(curSubRegion == subRegionAssignment);
                    
                    %Iterate through tracks
                    for trackIdx = 1:length(curSubRegionTracks)
                        %Iterate through all detections of the current track
                        for spotIdx = 1:length(curSubRegionTracks(trackIdx).Frame)
                            %Print current detection to .csv file
                            fprintf(fid, '%d,%d,%.4f,%d,%f,%f\n',...
                                spotNum, ...
                                curSubRegionTracks(trackIdx).Frame(spotIdx),...
                                curSubRegionTracks(trackIdx).TimeStamp(spotIdx),...
                                trackNum,...
                                curSubRegionTracks(trackIdx).xy(spotIdx,1),...
                                curSubRegionTracks(trackIdx).xy(spotIdx,2));
                            %Increase spot count
                            spotNum = spotNum+1;
                        end
                        %Increase track count
                        trackNum = trackNum + 1;
                    end
                    
                    %Close .csv file
                    fclose(fid);
                end
                
            end
        end
        
        ui.editFeedbackWin.String = 'Track export finished';
        
    end

    function CreateMovieCB(~,~)
        %User pressed File -> Create .avi movie
        
        %% ------Open framerange dialog------------------------------------
        nFramesOriginal = size(movieStacks{1},3);
        prompt = {'First frame:','Last frame:'};
        dlgtitle = 'Framerange';
        dims = [1 40];
        definput = {'1',num2str(nFramesOriginal)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        frames = str2double(answer{1}):str2double(answer{2});
        nFramesInAvi = frames(end) - frames(1);
        
        %% ------Open save dialog------------------------------------------
        [~,name,~] = fileparts(batch(curMovieIndex).movieInfo.fileName);
        [newFileName, newPathName] = uiputfile(fullfile(searchPath{2},strcat(name, '.avi')),'Choose a location to save Video');
        searchPath{2} = newPathName;
        fullVideoPath = fullfile(newPathName,newFileName);
        
        %% ------Create movie----------------------------------------------
        if ~isequal(newFileName,0)
            try
                
                %Create write object
                writerObj = VideoWriter(fullVideoPath);
                
                %Specify framerate as defined in the ui as "FPS"
                writerObj.FrameRate = str2double(ui.editFPS.String);
                
                %Open file
                open(writerObj)
                
                % Iterate through frames
                for curFrame = frames
                    
                    %Calculate progress and display in the feedback window
                    percentDone = round((curFrame-frames(1)) / nFramesInAvi*100);
                    if mod(percentDone,1) == 0
                        ui.editFeedbackWin.String = sprintf('Please wait while Movie is created: %3.0f %%', percentDone);
                    end
                    
                    %Change the value of the frame slider. This calls the
                    %FramSlider() function to update the plot to the
                    %current frame.
                    ui.sliderFrame.Value = curFrame;
                    
                    
                    %Get frame from the trackIt figure
                    frame = getframe(ui.axes1);
                    
                    %Cut off white borders
                    frame.cdata = frame.cdata(2:end-1,2:end-1,:);
                    
                    %Write frame to writer object
                    writeVideo(writerObj,frame);
                end
                
                ui.editFeedbackWin.String = char('Creating Movie: Complete', ui.editFeedbackWin.String(2:end,:));
                
                %Close file
                close(writerObj);
                
            catch ex
                close(writerObj);
                errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                    ex.stack(1).name, ex.stack(1).line, ex.message);
                fprintf(1, '%s\n', errorMessage)
            end
        end
    end

%Callbacks from tools menu

    function TrackExplorerCB(src,~)
        %User pressed Tools -> Track explorer or "Plot selected track" button
        %in the main figure

        
        %Check if movie has to be loaded
        if isempty(movieStacks{1})
            
            %Load movie
            movieStacks{1} = load_stack(batch(curMovieIndex).movieInfo.pathName, batch(curMovieIndex).movieInfo.fileName, ui);
            
            if isempty(movieStacks{1})
                warndlg('Movie not found. The track explorer requires the original movie data.','Movie not found')
                return
            end
        end
        
        
        switch src.Tag
            case 'trackselected'
                %User pressed "Plot selected track"
                
                %Get ID of selected track
                trackID = plotSettings.curTrackSel;
                
                %Get data cursor
                hDataTip = datacursormode(ui.hFig);
                
                if strcmp(hDataTip.Enable , 'off') && ~trackID
                    %Data cursor is not enabled so display a warning dialog
                    datacursormode on
                    msgbox(char('Please select a track with the Data Tips tool (see toolbar).', 'Data Tips tool has been enabled automatically for you.'), 'Error','error');
                    return
                elseif ~trackID
                    %Data cursor is enabled but no track is selected
                    msgbox('Please select a track.', 'Error','error');
                    return
                end
                
                
            case 'menucall'
                %User pressed "Plot selected track" button
                trackID = 1;
        end
        
        %Open track explorer 
        %TODO for future versions: Allow only one instance of the track explorer        
        
%         figs = findobj(allchild(groot), 'flat', 'Tag', 'track_explorer');
%         if ~isempty(figs)
%             test = findobj(allchild(figs), 'flat', 'Tag', 'panelTrackId');
%             test = findobj(allchild(test), 'flat', 'Tag', 'editTrackId');
%             test.String = trackID;
%             figure(figs)
%         else
            track_explorer(movieStacks{1}, batch(curMovieIndex).results, batch(curMovieIndex).movieInfo.frameCycleTime, trackID, pixelSize);
%         end
        
    end

    function KymographToolCB(~,~)
        %User pressed Tools -> Kymograph
        
        %Show note in the feedback window to draw a rectangular region
        feedbackWin = ui.editFeedbackWin.String;
        ui.editFeedbackWin.BackgroundColor = [.5 1 .5];
        ui.editFeedbackWin.String = char('Please draw a rectangular region into the movie');
        
        %Start rectangular drawing tool
        ROIFreehandHandle = drawrectangle('Color','w');
        
        %Reset feedback window
        ui.editFeedbackWin.String = feedbackWin;
        ui.editFeedbackWin.BackgroundColor = 'w';
        
        %Get coordinates of rectangle
        curROI = round(ROIFreehandHandle.Position);
        
        %Delete handle to drawn region
        delete(ROIFreehandHandle)
        
        %Get original movie stack
        curStack = movieStacks{1};
        
        %Define borders where to cut original movie
        yMin = curROI(2);
        yMax = curROI(2)+curROI(4)-1;
        xMin = curROI(1);
        xMax = curROI(1)+curROI(3)-1;
        
        %Cut out rectangular region from original movie
        imageStack = curStack(yMin:yMax,xMin:xMax,:);
        
        %Start kymograph tool
        kymograph_tool(imageStack)
        
    end

    function PredictParamsCB(~,~)
        %User pressed Tools -> predict tracking radii
        
        %Open plop gui and initialize plop object with current batch and
        %empty subregion number variable
        gui_plop(plop(batch, []))
    end

    function SplitMoviesCB(~,~)
        %User pressed Tools -> Movie splitter
        movie_splitter(searchPath{1});
    end

    function FindAllStack2CB(~,~)
        %User pressed Stuff -> Find stack 2 for all movies
        
        
        if isempty(ui.editReplaceString2.String)            
           msgbox(['Please enter a string in the fields "Find" and "Replace with". ',...
               'Eg. if your tracking movie name is CDX2_561_n001.tif and '...
               'the corresponding second movie name is CDX2_488_n001.tif, enter "561" in '...
               'the "Find" field and "488" in the "Replace with" field. '...
               'TrackIt will set all the second movie filenames by replacing "561" with "488".']) 
           
           return
        end
        
               
        %Iterate through all movies of the current batch
        for idx = 1:length(batch)
            
            %Current filename and pathname of tracking movie
            filename = batch(idx).movieInfo.fileName;
            pathname = batch(idx).movieInfo.pathName;
            
            %Get insert string
            insert = ui.editReplaceString2.String;
            newFilename = '';
            
            if isempty(ui.editReplaceString1.String)
                %No replacement string was entered so just search for
                %the insert string in the current folder
                
                %First check for direct match
                
                %Create search pattern
                searchPattern = [insert,'.tif*'];
                
                %Create list of files containing the insert string
                fileList = dir(fullfile(pathname,searchPattern));
                
                if isempty(fileList)
                    %No direct match found check if string is part of
                    %and filename
                    
                    %Create search pattern
                    searchPattern = ['*',insert,'*'];
                    
                    %Create list of files containing the insert string
                    fileList = dir(fullfile(pathname,searchPattern));
                end
                
                %Get resulting file names
                fileNameList = {fileList.name};
                
                if ~isempty(fileNameList)
                    %Use first file that has been found
                    newFilename = fileNameList{1};
                end
            else
                %String which should be replaced in tracking movie filename
                replaceString1 = ui.editReplaceString1.String;
                %Split filename at the position of replacement string
                splittedFilename = strsplit(filename, replaceString1);
                
                %Create new filename with the fileparts of the tracking
                %movie and the insert string from the user interface
                if length(splittedFilename)>1
                    newFilename = strcat(splittedFilename{1},insert,splittedFilename{2});
                end
            end
            
            %Save file and path of second movie in batch structure
            batch(idx).movieInfo.fileName2 = newFilename;
            batch(idx).movieInfo.pathName2 = pathname;
        end
        
    end

    function RemoveRoiFromAllMoviesCB(~,~)
        %User pressed Stuff -> Remove ROI from all movies
        
        %Show user dialog
        prompt = {'Enter the number of the region you want to have removed from all movies in this batch (1 = tracking region)'};
        dlgtitle = 'Remove regions';
        answer = inputdlg(prompt,dlgtitle);
        
        if isempty(answer)
            %Canceled
            return
        end
        
        %Get the entered roi number
        roiNum = str2double(answer{1})-1;
        
        %Get number of movies
        nMovies = length(batch);
        %Iterate through all movie of the current batch
        for movieIdx = 1:nMovies
            
            if roiNum == 0
                %Remove tracking ROI
                batch(movieIdx).ROI = {};
            elseif length(batch(movieIdx).subROI) >= roiNum
                %Remove sub-region
                batch(movieIdx).subROI(roiNum) = [];
            end
        end
        
        msgbox(['Region ', num2str(roiNum+1),' deleted from all movies'])
        
        %Update ui
        RoiChanged(false);
        FrameSlider()
        
    end

    function ReloadRoisCB(~,~)
        %User pressed Stuff -> Reload ROIs for all movies from ROI folder
        
        %Get number of movies
        nFiles = length(batch);
        
        for fileIdx = 1:nFiles
            
            % Check if ROI file exists
            [~,name,~] = fileparts(filesTable.FileName{fileIdx});
            roiFileName = fullfile(filesTable.PathName{fileIdx},'ROI',strcat(name, '.roi'));
            
            if exist(roiFileName, 'file') == 2 
                %Load ROI from .roi file
                loadedFile = load(roiFileName,'-mat');
                
                %Save roi and subroi in batch struct
                batch(fileIdx).ROI = loadedFile.ROI;
                batch(fileIdx).subROI = loadedFile.subROI;
            end

        end
        
        %Update region list and plot
        RoiChanged(false);
        FrameSlider()
        msgbox('Done!')

    end

%Callbacks from Analysis menu

    function TrackingDataAnalysisCB(~,~)
        %User pressed Analysis -> Tracking data analysis
        
        %Start data analysis tool
        data_analysis_tool({batch}, searchPath{2}, curMovieIndex, pixelSize);
    end

    function SpotStatisticsCB(~,~)
        %User pressed Analysis -> Spot statistics
        
        %Start spot analysis tool
        histograms_spots(batch,curMovieIndex,movieStacks{1})
    end

    function StartGridCB(~,~)
        %User pressed Analysis -> analyse dissociation rates
        
        %Get number of movies in current batch
        nMovies = length(batch);
        regionChoice = [];
        
        %Iterate through all movies and find maximum number of sub-regions
        %in a movie
        nSubRegions = 0;
        for movieIdx = 1:nMovies
            nSubRegions = max(nSubRegions, batch(movieIdx).results.nSubRegions);
        end
        
        if nSubRegions > 0
            %At least one movie has a sub-region so create a dialog window
            %where user can select the sub-region in which dissociation
            %rates should be analyzed
            
            
            canceled = false;
            
            %Create list of sub-regions
            popList = cell(nSubRegions+2,1);
            popList{1} = 'All regions';
            popList{2} = 'Region 1 (tracking-region)';
            for roiIdx = 2:nSubRegions+1
                popList{roiIdx+1} = ['Region ',num2str(roiIdx)];
            end
            
            %Create dialog
            d = dialog('Position',[300 300 200 200],'Name','Select a region','WindowKeyPressFcn',@KeyPressFcnCB,'CloseRequestFcn',@CloseRequestCB);
            uicontrol('Units','normalized','Parent',d,'Style','text','Position',[.05 .7 .8 .2], 'String','In which region do you want to analyze dissociation rates?','HorizontalAlignment','Left');
            uicontrol('Units','normalized','Parent',d,'Style','listBox','Position',[.05 .3 .9 .35], 'String',popList,'Callback',@PopUpCB);
            uicontrol('Units','normalized','Parent',d,'Position',[.05 .1 .4 .15],'String','OK', 'Callback',@OkCB);
            uicontrol('Units','normalized','Parent',d,'Position',[.5 .1 .4 .15],'String','Cancel', 'Callback',@CloseRequestCB);
            
            % Wait for d to close before running to completion
            uiwait(d);
            
            if canceled
                %User aborted
                return
            end
            
        end
        
        function PopUpCB(src,~)
            %User changed the selected region
            regionChoice = src.Value - 1;
            
            if regionChoice == 0
                %User selected "All regions"
                regionChoice = [];
            end
        end
        
        function CloseRequestCB(~,~)
            
            canceled = true;
            
            delete(gcf)
        end
        
        function OkCB(~,~)            
            delete(gcf)
        end
        
        function KeyPressFcnCB(~,event)
            if strcmp(event.Key, 'escape')
                %User pressed escape
                delete(gcf)
                canceled = true;
            end
        end


        if isempty(batch(1).results.spotsAll)
            %Batch is empty so start GRID without passing data
            GRID()
        else
            %Create histograms of track lengths for the selected region
            data = create_grid_data_from_batch_file(batch, regionChoice);
            %Pass track length histograms to GRID
            GRID(data)
        end
    end


%Callbacks from Help menu

    function AboutCB(~,~)
        %User pressed Help -> about
        
        msgbox(["TrackIt v1.5";...
            "";...
            "Publication:";...
            "Single molecule tracking and analysis framework including theory-predicted parameter settings. Sci Rep 11, 9465 (2021).";...
            "";...
            "Copyright (C) 2020 Timo Kuhn, Johannes Hettich, Rubina Davtyan and J. Christof M. Gebhardt";...
            "";...
            "This program is free software: you can redistribute it and/or modify";...
            "it under the terms of the GNU General Public License as published by";...
            "the Free Software Foundation, either version 3 of the License, or";...
            "(at your option) any later version.";...
            "";...
            "This program is distributed in the hope that it will be useful,";...
            "but WITHOUT ANY WARRANTY; without even the implied warranty of";...
            "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the";...
            "GNU General Public License for more details."])
        
    end

    function OpenManualCB(~,~)
        %User pressed Help -> Open manual
        open("TrackIt_manual.pdf")        
    end


%Other small stuff

    function CopyFilenameCB(src,~)
        %User clicked on the movie filename so copy the filename to the
        %clipboard
        
        %Save current feedback window entry
        feedbackWin = ui.editFeedbackWin.String;
        
        %Display current filename in command window and copy to clipboard
        filename = src.String(10:end);
        disp(filename)        
        clipboard('copy',filename)
        
        %Notify in feedback window and change background color to green
        ui.editFeedbackWin.String = 'Filename copied to clipboard and Matlab command window';
        ui.editFeedbackWin.BackgroundColor = [.5 1 .5];  
        
        %Update ui
        drawnow
        
        %Wait 1 second
        pause(1)
        
        %Change feedback window to original state
        ui.editFeedbackWin.BackgroundColor = 'w';
        ui.editFeedbackWin.String = feedbackWin;
        
        
    end

    function MouseWheelCB(~,callbackdata)
        %Executed by mouse wheel for scrolling through frames. Changing the
        %value of the frame slider automatically calls the FrameSlider()
        %function
        curFrame = round((get(ui.sliderFrame,'Value')));
        %Take care that the frame number stays inside the frame range
        if curFrame + callbackdata.VerticalScrollCount > ui.sliderFrame.Max
            ui.sliderFrame.Value = ui.sliderFrame.Max;
        elseif curFrame + callbackdata.VerticalScrollCount < 1
            ui.sliderFrame.Value = ui.sliderFrame.Min;
        else
            ui.sliderFrame.Value = curFrame + callbackdata.VerticalScrollCount;
        end
    end

    function cboxFindSpotsCB(src,~)
        %Enables or disables the ui elements associated with spot finding
        %when user enables or disables the "Find spots" checkbox
        if src.Value
            ui.editThresFactor.Enable = 'on';
            ui.editFramerange1.Enable = 'on';
            ui.editFramerange2.Enable = 'on';
        else
            ui.editThresFactor.Enable = 'off';
            ui.editFramerange1.Enable = 'off';
            ui.editFramerange2.Enable = 'off';
        end
        
    end

    function ConvertTrackingRadiusCB(~,~)
        %Executed when the user makes a right click on the tracking radius
        %field and presses "convert from Âµm"
        
        %Get tracking radius from ui
        trackingRadiusPx = ui.editTrackingRadius.String;
        
        %Open tracking radius conversion window
        [trackingRadiusPx, newPixelSize] = tracking_radius_conversion(trackingRadiusPx, pixelSize);
        
        
        if ~isempty(newPixelSize) && ~isempty(trackingRadiusPx)
            %Pixelsize and new tracking radius
            pixelSize = newPixelSize;
            ui.editTrackingRadius.String = trackingRadiusPx;
        end
    end

end
 