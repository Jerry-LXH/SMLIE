function histograms_spots(batch,curMovieIndex,curStack)


%Spot statistics tool which can be opened in the TrackIt main GUI via the
%menubar by selecting "Tools" -> "Spot statistics". Can be used to display
%statistics on the detected molecules in the movies of the batch file.
%
%histograms_spots(batch,curMovieIndex,curStack)
%
%Input:
%   batch           -   TrackIt batch struct array. See init_batch.m file
%                       for a description.
%   curMovieIndex   -   Index of the currently selected movie in TrackIt
%   curStackstack   -   3d-array of pixel values. 
%                       	1st dimension (column): y-coordinate of the image plane 
%                           2nd dimension (row):    x-coordinate of the image plane 
%                           3rd dimension: frame number

%Get results of current movie in TrackIT
curMovieResults = batch(curMovieIndex).results;

%Create user interface
S = createHistogramUI();

%Initialize variable for export of current histogram data
histData = struct;

%Create all results for current movie
CreateData()

%Update user interface and histogram
PopParamSelectionCB()


    function S = createHistogramUI()
        S.f   = figure('Units','normalized',...
            'Position',[0 0.3 .5 .55],...
            'Name','Spot statistics tool',...
            'CloseRequestFcn',@(~,~)CloseHistogram);
        
        
        S.ax  = axes(S.f,...
            'Units','normalized',...
            'Position',[0.45 0.1 0.5 0.85]);
        
        
        S.pan = uipanel(S.f,'Position',[0.01 0.02 0.35 0.95]);        

        S.btnNextMovie                 = uicontrol(S.pan,'Units','normalized','Position', [.52  .9   .3 .06],'String','Next movie','Callback',@MovieNumberCB);
        S.btnPreviousMovie           	= uicontrol(S.pan,'Units','normalized','Position', [.2  .9   .3 .06],'String','Previous movie','Callback',@MovieNumberCB);
        S.textMovie                    = uicontrol(S.pan,'Units','normalized','Position', [.4   .85   .3 .035],'Style','Text','String','Movie','HorizontalAlignment','Left');
        S.editMovie                    = uicontrol(S.pan,'Units','normalized','Position', [.52  .85  .15 .04],'Style','Edit','String',num2str(curMovieIndex),'HorizontalAlignment','Right','Callback',@MovieNumberCB);
        S.textMovie2                   = uicontrol(S.pan,'Units','normalized','Position', [.68  .85   .1 .035],'Style','Text','String','/','HorizontalAlignment','Left');
        S.textMovie3                   = uicontrol(S.pan,'Units','normalized','Position', [.70  .85   .3 .035],'Style','Text','String',num2str(length(batch)),'HorizontalAlignment','Left');
        S.textMovieName                 = uicontrol(S.pan,'Units','normalized','Position', [.05  .8   .9 .04],'Style','Text','String',batch(curMovieIndex).movieInfo.fileName);

        
        S.popParamSel = uicontrol('Parent',S.pan,...
            'Style','Listbox',...
            'Units','normalized',...
            'Position',[0.1 0.5 0.8 0.25],...
            'FontSize',11,...
            'String', {'Peak spot intensity','Fitted spot intensity','Spot SNR', 'Spot width', '#Spots in frames'},...
            'Callback',@(~,~)PopParamSelectionCB);
        
        S.txtHistLim = uicontrol('Parent',S.pan,...
            'Style','text',...
            'Units','normalized',...
            'FontSize',10,...
            'Position', [0.1 0.41 0.425 0.075],...
            'String', 'Histogram Limits');
        S.txtLim = uicontrol('Parent',S.pan,...
            'Style','text',...
            'Units','normalized',...
            'Position', [0.27 0.38 0.075 0.0625],...
            'String', '-');
        S.editLim1  = uicontrol('Parent',S.pan,...
            'Style','edit',...
            'Units','normalized',...
            'String','0',...
            'FontSize',9.5,...
            'Position',[0.1 0.4 0.15 0.05],...
            'Callback',@(~,~)UpdateHistogramCB);
        S.editLim2  = uicontrol('Parent',S.pan,...
            'Style','edit',...
            'String','1',...
            'FontSize',9.5,...
            'Units','normalized',...
            'Position',[0.35 0.4 0.15 0.05],...
            'Callback',@(~,~)UpdateHistogramCB);
        
        S.txtBinNum  = uicontrol('Parent',S.pan,...
            'Style','text',...
            'Units','normalized',...
            'FontSize',10,...
            'Position',[0.58 0.41 0.325 0.075],...
            'String','# Bins');
        S.editBinNum = uicontrol('Parent',S.pan,...
            'Style', 'edit',...
            'FontSize',9.5,...
            'Units','normalized',...
            'Visible','on',...
            'Position', [0.6625 0.4 0.15 0.05],...
            'String','15',...
            'Callback',@(~,~)UpdateHistogramCB);
        
        
        S.btnCopyToWorkspace  = uicontrol('Parent',S.pan,...
            'Units','normalized',...
            'FontSize',8,...
            'Position',[0.1 0.3 .8 0.06],...
            'String','Export to Matlab workspace',...
            'Callback',@CopyToWorkspaceCB);
        
        
    end

    function CloseHistogram()
        delete(gcf)
    end

    function CreateData()
        
        stdAll = [];
        intFittedAll = [];
        
        spots = curMovieResults.spotsAll;
        nFrames = length(spots);
        
        %% Get width and fitted intensities of spots
        for k = 1:nFrames
            if ~isempty(spots{k})
                stdAll  = [stdAll; spots{k}(:,5)];
                intFittedAll  = [intFittedAll; spots{k}(:,3)];
            end
        end
        
        %% Get spots per frame
        [spotsPerFrame,~] = cellfun(@size, curMovieResults.spotsAll, 'UniformOutput', false);
        spotsPerFrame = cell2mat(spotsPerFrame);
        
        %% Calculate SNR and brightest pixel
        
        intPeak = cell(1,nFrames);
        SNRspots = cell(1,nFrames);
        
        spotsAll = curMovieResults.spotsAll;
        
        parfor k = 1:nFrames     
            [SNRspots{k}, intPeak{k}] = calc_snr(spotsAll{k} ,curStack(:,:,k));
        end
        
        histData.intFittedAll = intFittedAll;
        histData.stdAll  = stdAll;
        histData.spotsPerFrame = spotsPerFrame;
        histData.SNRall = [SNRspots{:}];
        histData.intPeakAll= [intPeak{:}];
    end

    function UpdateHistogramCB()
        
        %Switch between selected parameter in listbox and update histogram
        switch S.popParamSel.Value
            case 1%Peak intensity
                S.hist = histogram(histData.intPeakAll,...
                    'BinLimits',[str2double(S.editLim1.String) str2double(S.editLim2.String)],...
                    'NumBins',str2double( S.editBinNum.String),...
                    'Parent',S.ax);
                
                histData.currentHistogramData = [(S.hist.BinEdges(2:end)-(S.hist.BinEdges(2)-S.hist.BinEdges(1))/2)', S.hist.Values'];
                
                xlabel(S.ax,'Intensity [a.u.]')
                ylabel(S.ax,'Spot Count')
            case 2%Fitted intensity
                S.hist = histogram(histData.intFittedAll,...
                    'BinLimits',[str2double(S.editLim1.String) str2double(S.editLim2.String)],...
                    'NumBins',str2double( S.editBinNum.String),...
                    'Parent',S.ax);
                
                histData.currentHistogramData = [(S.hist.BinEdges(2:end)-(S.hist.BinEdges(2)-S.hist.BinEdges(1))/2)', S.hist.Values'];
                
                xlabel(S.ax,'Intensity [a.u.]')
                ylabel(S.ax,'Spot Count')
            case 3%SNR
                S.hist = histogram(histData.SNRall,str2double( S.editBinNum.String),...
                    'BinLimits',[str2double(S.editLim1.String), str2double(S.editLim2.String)],...
                    'Parent',S.ax);
                
                histData.currentHistogramData = [(S.hist.BinEdges(2:end)-(S.hist.BinEdges(2)-S.hist.BinEdges(1))/2)', S.hist.Values'];
                
                xlabel(S.ax,'SNR')
                ylabel(S.ax,'Spot Count')
            case 4%Width
                S.hist = histogram(histData.stdAll,...
                    'BinLimits',[str2double(S.editLim1.String) str2double(S.editLim2.String)],...
                    'NumBins',str2double(S.editBinNum.String),...
                    'Parent',S.ax);
                
                histData.currentHistogramData = [(S.hist.BinEdges(2:end)-(S.hist.BinEdges(2)-S.hist.BinEdges(1))/2)', S.hist.Values'];
                
                xlabel(S.ax,'Spot Width [px]')
                ylabel(S.ax,'Spot Count')
            case 5 %Spots in frames
                %                 S.hist = plot(str2double(S.editLim1.String):str2double(S.editLim2.String),histData.spotsPerFrame(str2double(S.editLim1.String):str2double(S.editLim2.String)),...
                %                     'Parent',S.ax);
                
                S.hist = stairs(str2double(S.editLim1.String):str2double(S.editLim2.String),histData.spotsPerFrame(str2double(S.editLim1.String):str2double(S.editLim2.String)),...
                    'Parent',S.ax);
                
                histData.currentHistogramData = histData.spotsPerFrame;
                xlabel(S.ax,'Frame number')
                ylabel(S.ax,'Number of spots')
        end
        
    end

    function PopParamSelectionCB()
        %Executed when user selects a statistics parameter in the listbox
        
        S.editBinNum.Enable = 'on';
        
        if isempty(S.editLim1.String) || isempty(S.editLim2.String)
            S.editLim1.String = '0';
            S.editLim2.String = '1';
        end
        
        %Switch between selected parameter in listbox and update histogram
        switch S.popParamSel.Value
            
            case 1 %Peak intensity
                if histData.SNRall == -1
                    CreateData()
                end                
                S.editLim1.String = num2str(round(min(histData.intPeakAll),3));
                S.editLim2.String = num2str(round(max(histData.intPeakAll),3));
                S.editBinNum.String = '100';
            case 2 %Fitted intensity
                if histData.SNRall == -1
                    CreateData()
                end
                
                S.editLim1.String = num2str(round(min(histData.intFittedAll),3));
                S.editLim2.String = num2str(round(max(histData.intFittedAll),3));
                S.editBinNum.String = '100';
            case 3 %SNR
                if histData.SNRall == -1
                    CreateData()
                end
                S.editLim1.String = '0';
                S.editLim2.String = num2str(ceil(max(histData.SNRall)));
                S.editBinNum.String = num2str((ceil(max(histData.SNRall))-floor(min(histData.SNRall)))*2);
            case 4 %Width
                S.editLim1.String = '0';
                S.editLim2.String = num2str(round(max(histData.stdAll),1));
                S.editBinNum.String = '100';
            case 5 %Spots in frames
                S.editLim1.String = num2str(1);
                S.editLim2.String = num2str(size(curStack,3));
                S.editBinNum.String = length(str2double(S.editLim1.String):str2double(S.editLim2.String));
                S.editBinNum.Enable = 'off';
        end
        
        UpdateHistogramCB()
    end

    function CopyToWorkspaceCB(src,~)
        %Executed when user presses "Export to Matlab workspace" button
        
        src.BackgroundColor = 'r';
        drawnow
        
        if histData.SNRall == -1
            calc_snr()
        end
        assignin('base','spotHistogramResults',histData);
        
        src.BackgroundColor = [.94 .94 .94];
        
    end

    function MovieNumberCB(src,~)
        %Executed if the movie number ist changed by the user
        
        previousMovieNumber = curMovieIndex;
        
        %Check which is the new movie number taking care that it is not
        %below 1 and not higher than the amount of movies in the batch
        if strcmp(src.String,'Previous movie') && previousMovieNumber > 1
            %User clicked prevuious movie button
            curMovieIndex =  previousMovieNumber - 1;
        elseif strcmp(src.String,'Next movie') && previousMovieNumber < length(batch)
            %User clicked next movie button
            curMovieIndex = previousMovieNumber + 1;
        elseif str2double(src.String) <= length(batch) && str2double(src.String) > 0
            %User entered a movie number in the edit field
            curMovieIndex = str2double(src.String);
        end
        
        %Switch background color of ui control until movie has been loaded
        src.BackgroundColor = 'r';
        drawnow
        
        uiDummy.editFeedbackWin.String = '';
        
        %Load movie
        curStack = load_stack(batch(curMovieIndex).movieInfo.pathName, batch(curMovieIndex).movieInfo.fileName, uiDummy);   
        
        %Retreive results of current movie from batch structure
        curMovieResults = batch(curMovieIndex).results;
        
        %Create statistics results
        CreateData()
        
        %Update ui and histogram
        PopParamSelectionCB()
        
        %Reset ui control color
        src.BackgroundColor = [.94 .94 .94];
        
        %Update movie number and movie name in ui
        S.editMovie.String = num2str(curMovieIndex);
        S.textMovieName.String = batch(curMovieIndex).movieInfo.fileName;
        
    end

end











