function kymograph_tool(imageStack)


%GUI based tool to display kymographs (see TrackIt manual)
%
%
%kymograph_tool(imageStack)
%
%Input:
%
%   imageStack   -   3d-array of pixel values.
%                       1st dimension (column): y-coordinate of the image plane
%                       2nd dimension (row):    x-coordinate of the image plane
%                       3rd dimension: frame number

stackSize = size(imageStack);

%% Create figure, axes and ui controls
hFig = figure('Units','normalized',...
    'Position',[0.2 .2 .6 .6]);

posLeft = .02;
sizeLeft = .45;
posRight = .5;
sizeRight = .45;

axFrameImage     = axes(hFig,'Position',[posLeft 0.38 sizeLeft 0.6]);
axKymVer        = axes(hFig,'Position',[posLeft 0.17 .95 0.12]);
axKymHor        = axes(hFig,'Position',[posLeft 0.01 .95 0.12]);
axInt           = axes(hFig,'Position',[posRight 0.5 sizeRight 0.45]);

%Slider and text box with frame number
textFrame       = uicontrol('Units','normalized','Position',[.01 .354 .3 .025],'Style','Text','String','Frame in original movie: 1','HorizontalAlignment','Left');
sliderFrame     = uicontrol('Units','normalized','Position',[.01 .32 .95 .025],'Style','slider','Min',1,'Max',stackSize(3),'Value',1,'SliderStep',[min(1,1/(stackSize(3)-1)) 1/min((stackSize(3)-1),10)]);
addlistener(sliderFrame,'Value','PostSet',@(~,~)FrameSlider);

%Field where number of frames shown in kymograph can be entered
uicontrol('Units','normalized',...
    'Position',[posRight  .38   .15 .05],...
    'Style','Text',...
    'FontSize',10,...
    'String','No. of frames shown in kymograph and intensity plot',...
    'HorizontalAlignment','Left');

editFramesInKym = uicontrol('Units','normalized',...
    'Position',[posRight + .18  .38   .05 .04],...
    'Style','Edit',...
    'String',stackSize(3),...
    'HorizontalAlignment','Left',...
    'Callback',@FrameSlider);


uicontrol('Units','normalized',...
    'Position',[posRight+.28  .38   .15 .04],...
    'FontSize',10,...
    'String','Export to Matlab workspace',...
    'HorizontalAlignment','Left',...
    'Callback',@(~,~)ExportCB);


set(hFig, 'WindowScrollWheelFcn', {@MouseWheelCB});

%% Initialize images and plots
imageHandle     = imshow(imageStack(:,:,1),[],'Parent',axFrameImage);

%Create kymograph projections
[kymoVer, kymoHor, intMax] = createKymographs();

%Plot intensity
hold(axInt, 'on')
hIntPlot = plot(axInt, intMax,'k.-');
%Create red circle indicating intensity in current frame
intMarker     = plot(axInt, NaN, NaN,'Color','r','Marker','o','Linestyle','none');
hold(axInt, 'off')
%Add title to intensity plot
title(axInt,'Max. intensity in frame')

%Show image of vertical kymograph
hKymoVerImage = imshow(kymoVer,[], 'Parent',axKymVer);
hold(axKymVer, 'on')
%Plot red line indicating current frame in kymograph
kymVerLine = plot(axKymVer, NaN,NaN,'red','linewidth',1.5);
%Add colorbar and title
colorbar(axKymVer)
title(axKymVer,'kymograph vertical (yt)')
hold(axKymVer, 'off')

%Show image of horizontal kymograph
hKymoHorImage = imshow(kymoHor,[],'Parent',axKymHor);
hold(axKymHor, 'on')
%Plot red line indicating current frame in kymograph
kymHorLine = plot(axKymHor, NaN,NaN,'red','linewidth',1.5);
%Add colorbar and title
colorbar(axKymHor)
title(axKymHor,'kymograph horizontal (xt)')
hold(axKymHor, 'off')

FrameSlider()

%% Functions

    function FrameSlider(~,~)
        %Get current frame
        curFrame = round(sliderFrame.Value);
        
        %Get number of frames to show
        nFramesInKym = str2double(editFramesInKym.String);
        
        %Caculate first and last frame of kymographs to be shown
        kymoFirstFrame = min(max(1,curFrame-nFramesInKym/2), stackSize(3)-nFramesInKym+1);
        kymoLastFrame = kymoFirstFrame+nFramesInKym-1;
        
        %Update kymograph images
        hKymoVerImage.CData = kymoVer(:, kymoFirstFrame:kymoLastFrame);
        hKymoHorImage.CData = kymoHor(:, kymoFirstFrame:kymoLastFrame);
        
        %Update x axis limits and lookup table for kymographs
        axKymVer.XLim = [0.5 nFramesInKym];
        axKymVer.CLim = [min(hKymoVerImage.CData,[],'all') max(hKymoVerImage.CData,[],'all')];
        axKymHor.XLim = [0.5 nFramesInKym];
        axKymHor.CLim = [min(hKymoHorImage.CData,[],'all') max(hKymoHorImage.CData,[],'all')];
        
        %Update intensity plot 
        hIntPlot.XData = kymoFirstFrame:kymoLastFrame;
        hIntPlot.YData = intMax(kymoFirstFrame:kymoLastFrame);
        axInt.XLim = [kymoFirstFrame kymoLastFrame];
        
        %Update current frame
        imageHandle.CData = imageStack(:,:,curFrame);
        axFrameImage.CLim = [min(imageHandle.CData,[],'all') max(imageHandle.CData,[],'all')];
        
        %Update text indicating current frame
        textFrame.String = sprintf('Frame in original movie: %d', curFrame);
                
        %Update red lines and circle indicating current frames in
        %kymographs nad intensity plot
        kymoHeightVer = size(kymoVer,1);
        kymoHeightHor = size(kymoHor,1);
        
        set(kymHorLine, 'xdata',[curFrame+1.5-kymoFirstFrame,curFrame+1.5-kymoFirstFrame],'ydata',[0.5,kymoHeightHor+.5])
        set(kymVerLine, 'xdata',[curFrame+1.5-kymoFirstFrame,curFrame+1.5-kymoFirstFrame],'ydata',[0.5,kymoHeightVer+.5])
        set(intMarker, 'xdata',curFrame,'ydata',intMax(curFrame))
                
    end


    function [kymoVerFull,kymoHorFull, intMean] = createKymographs()
        %Creates a maximum projection in both horizontal and vertical
        %direction over time
        
        kymoVerFull = zeros(stackSize(1),stackSize(3));
        kymoHorFull = zeros(stackSize(2),stackSize(3));
        intMean = zeros(stackSize(3),1);
        
        %Iterate through frames and create projection of each frame and
        %instert it into the kymograph
        for frameIdx = 1:stackSize(3)
            %Get current frame
            curSpotImage = imageStack(:,:,frameIdx);
            
            %Insert projection of current image into kymograph
            kymoVerFull(:,frameIdx) = max(curSpotImage,[],2);
            kymoHorFull(:,frameIdx) = max(curSpotImage,[],1)';
            intMean(frameIdx) = max(curSpotImage(:));            
        end
               
    end

    function ExportCB()
        %User pressed "Export to matlab workspace"
        
        kymographExport.horizontalKymograph = hKymoHorImage.CData;
        kymographExport.verticalKymograph = hKymoVerImage.CData;
        kymographExport.spotImage = imageHandle.CData;
        kymographExport.intensityPlot = [hIntPlot.XData', hIntPlot.YData'];
        
        assignin('base', 'kymographExport', kymographExport)
        
    end


    function MouseWheelCB(~,callbackdata)
        %Executed by mouse wheel for scrolling through frames
        curFrame = round((get(sliderFrame,'Value')));
        %Take care that the frame number stays inside the frame range
        if curFrame + callbackdata.VerticalScrollCount > sliderFrame.Max
            sliderFrame.Value = sliderFrame.Max;
        elseif curFrame + callbackdata.VerticalScrollCount < 1
            sliderFrame.Value = sliderFrame.Min;
        else
            sliderFrame.Value = curFrame + callbackdata.VerticalScrollCount;
        end
    end
end