function multiple_batches_analyzer(varargin)


if nargin == 1
    startingPath = varargin{1};
else
    startingPath = pwd;
end


%Initialize table with filenames and paths
filesTable = cell2table(cell(0,2));
filesTable.Properties.VariableNames = {'FileName','PathName'};

%Create user interface
ui = CreateMovieSelectorUI();

    function ui = CreateMovieSelectorUI()
        
        ui.f = figure('Units','normalized',...
            'Position',[0.4 0.4 .42 .45],...
            'MenuBar','None',...
            'Name','Multiple batch file analyzer',...
            'NumberTitle','off',...
            'WindowKeyPressFcn',@KeyPressFcnCB,...
            'CloseRequestFcn',@CloseRequestCB);
        ui.lb = uicontrol('Style','Listbox'...
            ,'Units','normalized',...
            'Value', [],...
            'Min', 0, 'Max', 2,...
            'Position',[0.05 0.48 .7 .5]);
        ui.tb = uitable('Units','normalized',...
            'Position',[0.05 0.01 .7 .45],...
            'ColumnName',{'Threshold factor';'Tracking radius';'Min. track length';'Gap frames';'Min. track length before gap frame'},...
            'Data',{'3','2','2','1','0';'3','2','2','1','0'},...
            'ColumnEditable',[true,true,true,true,true]);
        
        buttonPosHor = .76;
        buttonWidth = .22;
        buttonHeight = .06;
        
        filesButton = uicontrol('String','Select batch files',...
            'Units','normalized',...
            'Position',[buttonPosHor .92 buttonWidth buttonHeight],...
            'Callback',@AddFilesCB);
        
        removeButton = uicontrol('String','Remove selected files',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.85 buttonWidth buttonHeight],...
            'Callback',@RemoveFilesCB);
        
        ui.findSpotsCheckbox = uicontrol('String','Find spots',...
            'Units','normalized','Style','checkbox',...
            'Value',1,...
            'Position',[buttonPosHor .78 buttonWidth buttonHeight],...
            'Callback',@FindSpotsCB);

         frameRangeText1 = uicontrol('String','Framerange',...
            'Units','normalized','Style','Text',...
            'HorizontalAlignment','Left',...
            'Position',[buttonPosHor .69 buttonWidth-.1 buttonHeight]);
        
        ui.frameRangeEdit1 = uicontrol('String','1',...
            'Units','normalized','Style','edit',...
            'Position',[buttonPosHor+.09 .71 buttonWidth-.17 buttonHeight-.01]);
        
        frameRangeText2 = uicontrol('String','-',...
            'Units','normalized','Style','text',...
            'Position',[buttonPosHor+.14 .69 buttonWidth-.2 buttonHeight]);
                
        ui.frameRangeEdit2 = uicontrol('String','Inf',...
            'Units','normalized','Style','edit',...
            'Position',[buttonPosHor+.17 .71 buttonWidth-.17 buttonHeight-.01]);

        trackingMethodText = uicontrol('String','Tracking algorithm',...
            'Units','normalized','Style','text',...
            'Position',[buttonPosHor .62 buttonWidth buttonHeight]);
        
        ui.trackingMethodMenu = uicontrol('String',{'Nearest neighbour','u-track random motion','u-track linear+random motion'},...
            'Units','normalized','Style','popupmenu',...
            'Position',[buttonPosHor .58 buttonWidth buttonHeight]);
                
        trackingMethodText = uicontrol('String','Sub-region assignment',...
            'Units','normalized','Style','text',...
            'Position',[buttonPosHor .50 buttonWidth buttonHeight]);
        
        ui.subRoiBorderHandlingMenu = uicontrol('String',...
            {'Assign by first appearance', 'Split tracks at borders', 'Delete tracks crossing borders', 'Only use tracks crossing borders'},...
            'Units','normalized','Style','popupmenu',...
            'Position',[buttonPosHor .46 buttonWidth buttonHeight]);

        nParametersText = uicontrol('String','Amount of tracking parameter sets',...
            'Units','normalized','Style','Text',...
            'HorizontalAlignment','Left',...
            'Position',[buttonPosHor+.01 .38 buttonWidth-.05 buttonHeight],...
            'Callback',@AddFilesCB);
        
        nParametersEdit = uicontrol('String','2',...
            'Units','normalized','Style','Edit',...
            'Position',[.92 .385 .05 buttonHeight-.01],...
            'Callback',@NParametersCB);
        
        ui.editProgressFeedbackWin = uicontrol('String','',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.29 buttonWidth .07],...
            'Style','Text','BackgroundColor',[.9 .9 .9]);
        
        ui.editFeedbackWin = uicontrol('String','',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.08 buttonWidth .20],...
            'Style','Text','BackgroundColor',[.9 .9 .9]);
        
        startButton = uicontrol('String','Start',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.01 buttonWidth buttonHeight],...
            'Callback',@StartCB);
        
    end

    function KeyPressFcnCB(~,event)
        %Close figure if esc is pressed
        if strcmp(event.Key, 'escape')
            delete(gcf)
        end
    end

    function CloseRequestCB(~,~)
        delete(gcf)
    end

    function FindSpotsCB(src,~)
        
        data = ui.tb.Data;
        
        if src.Value
            ui.tb.ColumnEditable = [true,true,true,true,true];
            ui.tb.Data(:,1) = repmat({'3'},size(data,1),1);
            ui.frameRangeEdit1.Enable = 'on';
            ui.frameRangeEdit2.Enable = 'on';
        else
            ui.tb.ColumnEditable = [false,true,true,true,true];
            ui.tb.Data(:,1) = repmat({''},size(data,1),1);
            ui.frameRangeEdit1.Enable = 'off';
            ui.frameRangeEdit2.Enable = 'off';
        end
    end

    function NParametersCB(src,~)
        %Executed when "Amount of tracking parameters" is changed
        
        %Get current table data
        data = ui.tb.Data;
        
        %Get amount of desired parameters sets
        nSequences = str2double(src.String);
        
        %Calculate difference between current amount of parameters sets and new
        %amount of parameters sets
        nAdditional = nSequences - size(data,1);
        
        if nAdditional > 0
            %More parameters sets required
            
            ui.tb.Data = [data; repmat(data(end,:),nAdditional,1)];
            
        elseif nAdditional < 0
            %Less parameters sets required so just delete amount of excess rows
            ui.tb.Data = data(1:end+nAdditional,:);
        end
    end

    function AddFilesCB(~,~)
        %User pressed "Add files" button
        
        %Open file dialog box
        [fileNameListNew,pathName] = uigetfile({'*.mat*'},'Select files you want to split', 'MultiSelect', 'on',startingPath);
        
        if isequal(fileNameListNew,0) %User didn't choose a file
            return
        elseif ~iscell(fileNameListNew) %Check if only one file has been chosen
            fileNameListNew = {fileNameListNew};
        end
        
        %Save path for next usgage of file dialog box
        startingPath = pathName;
        
        %Create table containing new filenames and paths
        pathNameListNew = repmat({pathName},1,length(fileNameListNew));
        filesTableNew = table(fileNameListNew',pathNameListNew','VariableNames',{'FileName','PathName'});
        
        %Add new filetable to existing one
        filesTable = [filesTable;filesTableNew];
        
         %Show filenames in ui
        ui.lb.String = filesTable.FileName;
        
        if size(ui.lb.String,1) > 1
            %Set selected value to first entry
            ui.lb.Value = 1;
            %Make sure more than one file can be selected in the ui list
            ui.lb.Max = 2;
        end
    end

    function RemoveFilesCB(~,~)
        %User pressed "Remove selected files" button
        
        %Remove selection from filestable        
        filesTable(ui.lb.Value,:) = [];
        
        %Update list of files in ui
        ui.lb.String = filesTable.FileName;
        
        if size(ui.lb.String,1) <= 1
            %Less than 2 files are less so disable multiselection and set
            %selected value to 1
            ui.lb.Value = 1;
            ui.lb.Max = 0;
        elseif size(ui.lb.String,1) < ui.lb.Value(end)
            %Make sure that selected value is not higher than the amount of
            %files in list
            ui.lb.Value = numel(ui.lb.String);
        end
    end

    function StartCB(src,~)
        %Make sure current character is not str+s so user can press strg + s to cancel
        set(gcf,'currentch',char(1))
        
        %User pressed "Start" button
        src.BackgroundColor = 'r';
        src.String = 'Press Ctrl + X to stop';
        ui.editProgressFeedbackWin.String = '';
        drawnow;

        %Get number of batch files to analyze
        nBatchFiles = height(filesTable);
        
        %Get tracking settings
        trackingParams = cellfun(@str2double,ui.tb.Data);
        boolFindSpots = ui.findSpotsCheckbox.Value;
        trackingMethod = ui.trackingMethodMenu.String{ui.trackingMethodMenu.Value};
        subRoiBorderHandling = ui.subRoiBorderHandlingMenu.String{ui.subRoiBorderHandlingMenu.Value};
        
        %Get amount of tracking parameter sets
        nParameterSets = size(trackingParams,1);
        
        if ~isempty(filesTable) && ~isempty(trackingParams)
            
            %Iterate through files
            tic
            for batchFileIdx = 1:nBatchFiles
                %Monitor progress in ui
                ui.editProgressFeedbackWin.String = ['Tracking batch ' num2str(batchFileIdx) ' of ' num2str(nBatchFiles)];
                drawnow
                
                %Get current filename and pathname
                curFileName = filesTable.FileName{batchFileIdx};
                curPathName = filesTable.PathName{batchFileIdx};
                
                loadedFile = load(fullfile(curPathName,curFileName));
                
                if ~isfield(loadedFile,'batch')
                    ui.editFeedbackWin.String = char('Please choose a valid batch file');
                    return
                end
                                
                loadedBatch = loadedFile.batch;
                loadedFilesTable = loadedFile.filesTable;
                
                [~,curFileName,~] = fileparts(curFileName);
                
                for parameterSetIdx = 1:nParameterSets
                    ui.editProgressFeedbackWin.String = char(ui.editProgressFeedbackWin.String(1,:), ['Using parameter set ' num2str(parameterSetIdx) ' of ' num2str(nParameterSets)]);                    
                    %-------------------------------
                    %Get settings for spot detection
                    para.boolFindSpots = boolFindSpots;
                    para.thresholdFactor = trackingParams(parameterSetIdx,1);
                    para.frameRange = [str2double(ui.frameRangeEdit1.String),...
                        str2double(ui.frameRangeEdit2.String)];
                    
                    %Get settings for tracking
                    para.trackingMethod = trackingMethod;
                    
                    para.tlConditions = [];
                    para.trackingRadius = trackingParams(parameterSetIdx,2);
                    para.minTrackLength = trackingParams(parameterSetIdx,3);
                    para.gapFrames = trackingParams(parameterSetIdx,4);
                    para.minLengthBeforeGap = trackingParams(parameterSetIdx,5);
                    
                    %Get sub-region handling option
                    para.subRoiBorderHandling = subRoiBorderHandling;
                    
                    %Get trackIt version and time stamp
                    para.trackItVersion = '';
                    
                    %--------------------------------------------------
                    
                    [newBatch, boolCancelled]  = tracking_routine(loadedBatch, para, [], ui);
                    
                    if boolCancelled
                        src.BackgroundColor = [.94 .94 .94]; %Reset button color to gray
                        src.String = 'Start'; %Reset button string
                        return
                    end
                    
                    %Save batch file
                    ui.editProgressFeedbackWin.String = char(ui.editProgressFeedbackWin.String(1,:), ['Saving Batch file...' num2str(parameterSetIdx) ' of ' num2str(nParameterSets)]);
                    drawnow
                    
                    
                    
                    settingsName = ['_tf_', num2str(para.thresholdFactor,'%.2f'),...
                        '_tr_',num2str(para.trackingRadius,'%.2f'),...
                        '_mtl_',num2str(para.minTrackLength),...
                        '_gf_',num2str(para.gapFrames),...
                        '_befGf_',num2str(para.minLengthBeforeGap)];
                    S.filesTable = loadedFilesTable;
                    S.batch = newBatch;
                    
                    fileName = [curFileName, settingsName, '.mat'];
                    save(fullfile(curPathName,fileName),'-struct','S')
                end
            end
            
            if  ~isempty(get(gcf,'CurrentCharacter')) &&  double(get(gcf,'CurrentCharacter')) == 24
                %Analysis aborted by user 
                ui.editProgressFeedbackWin.String = char('Analysis stopped by user');
            else
                %Analysis finished
                ui.editProgressFeedbackWin.String = ['Analysis finished in ', num2str(floor(toc/60)),' min, ', num2str(ceil(toc - floor(toc/60) * 60)), ' sec'];
                ui.editFeedbackWin.String = '';
            end
            
        end
        
        
        src.BackgroundColor = [.94 .94 .94]; %Reset button color to gray
        src.String = 'Start'; %Reset button string
    end

end
