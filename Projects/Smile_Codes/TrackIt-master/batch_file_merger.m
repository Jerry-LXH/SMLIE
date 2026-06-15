function batch_file_merger(startingPath)
%%Merge multiple batch files 

if nargin == 0
    %Starting path for opening the file dialog was not passed in function
    %call so intialize it to the current Matlab path
    startingPath = pwd;
end

%Initialize table with filenames and paths
filesTable = cell2table(cell(0,2));
filesTable.Properties.VariableNames = {'FileName','PathName'};

%Create user interface
S = CreateMovieSelectorUI();

    function S = CreateMovieSelectorUI()
        
        S.f = figure('Units','normalized',...
            'Position',[0.1 0.5 .3 .2],...
            'MenuBar','None',...
            'Name','Batch file merger',...
            'NumberTitle','off',...
            'WindowKeyPressFcn',@KeyPressFcnCB,...
            'CloseRequestFcn',@CloseRequestCB);
        
        S.editDestination = uicontrol('Style','Edit'...
            ,'Units','normalized',...
            'HorizontalAlignment','Left',...
            'String', [],...
            'Position',[0.05 0.8 .7 .1]);
        
        S.lb = uicontrol('Style','Listbox'...
            ,'Units','normalized',...
            'Value', [],...
            'Min', 0, 'Max', 2,...
            'Position',[0.05 0.05 .7 .7]);
        
        buttonPosHor = .78;
        buttonWidth = .2;
        buttonHeight = .15;
        
        filesButton = uicontrol('String','Add files',...
            'Units','normalized',...
            'Position',[buttonPosHor .8 buttonWidth buttonHeight],...
            'Callback',@AddFilesCB);
        
        destinationButton = uicontrol('String','Destination',...
            'Units','normalized',...
            'Position',[buttonPosHor .6 buttonWidth buttonHeight],...
            'Callback',@DestinationCB);
        
        removeButton = uicontrol('String','Remove selected files',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.42 buttonWidth buttonHeight],...
            'Callback',@RemoveFilesCB);
        
        S.feedbackWin = uicontrol('String','',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.22 buttonWidth .17],...
            'Style','Text','BackgroundColor',[.9 .9 .9]);%
        
        okButton = uicontrol('String','Start merge',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.05 buttonWidth buttonHeight],...
            'Callback',@StartMergeCB);
        
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

    function AddFilesCB(~,~)
        %User pressed "Add files" button
        
        %Open file dialog box
        [fileNameListNew,pathName] = uigetfile({'*.mat*'},'Select files you want to merge', 'MultiSelect', 'on',startingPath);
        
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
        S.lb.String = filesTable.FileName;
        
        
        if size(S.lb.String,1) > 1
            %Set selected value to first entry
            S.lb.Value = 1;
            %Make sure more than one file can be selected in the ui list
            S.lb.Max = 2;
        end
        
        if isempty(S.editDestination.String)
            %Set destination field
            S.editDestination.String = [startingPath,'merged_batch.mat'];
        end
        
    end

    function DestinationCB(~,~)
        %User pressed destiantion button
        
        %Open file dialog box
        [fileName,pathName] = uiputfile('*.mat','Choose filename for saving merged batch job in .mat file',startingPath);
        
        %Save path for next usgage of file dialog box  
        startingPath = pathName;
        
        %Write filename and path into ui text field
        S.editDestination.String = fullfile(pathName,fileName);
    end

    function RemoveFilesCB(~,~)
        %User pressed "Remove selected files" button
        
        %Remove selection from filestable
        filesTable(S.lb.Value,:) = [];
        
        %Update list of files in ui
        S.lb.String = filesTable.FileName;
        
        if size(S.lb.String,1) <= 1
            %Less than 2 files are less so disable multiselection and set
            %selected value to 1
            S.lb.Value = 1;
            S.lb.Max = 0;
        elseif size(S.lb.String,1) < S.lb.Value(end)
            %Make sure that selected value is not higher than the amount of
            %files in list
            S.lb.Value = numel(S.lb.String);
        end
    end

    function StartMergeCB(~,~)
        %User pressed "Start merge" button
        
        newFilePath = S.editDestination.String;
        
        if ~isempty(filesTable)
            %Get number of files in table
            nBatches = height(filesTable);
            
            %Initialize variables for containing all batches and files
            %tables
            batches = cell(nBatches,1);
            filesTables = cell(nBatches,1);
            
            %Iterate through all files and concatenate batch structure arrays
            for fileIdx = 1:nBatches
                %Monitor progress in ui
                S.feedbackWin.String = char('Loading batch file ', [num2str(fileIdx), ' of ' , num2str(nBatches)]);
                drawnow
                
                %Get current filename and load current batch file
                curBatchFile = fullfile(filesTable.PathName{fileIdx},filesTable.FileName{fileIdx});
                loadedBatch = load(curBatchFile);
                
                %Check if file contains a batch
                if ~isfield(loadedBatch,'batch')
                    S.feedbackWin.String = 'Please choose a valid batch file';
                    return
                end
                
                %Write batch and filestable into cell arrays
                batches{fileIdx} = loadedBatch.batch;
                filesTables{fileIdx} = loadedBatch.filesTable;
            end
            
            S.feedbackWin.String = 'Merging batch files';            
            drawnow
            
            %Catenate batches cell array into a structure array containing
            %data of all movies of all files
            newBatch = vertcat(batches{:});
            
            %Catenate cell array into a table containing all filenames,
            %filepaths and frame cycle times of all movies
            newFilesTable = vertcat(filesTables{:});
            
            S.editFeedbackWin.String = 'Saving batch file';
            drawnow
            
            %Prepare structure variables for saving
            newBatchStruct.filesTable = newFilesTable;
            newBatchStruct.batch = newBatch;
            
            %Save new batch file
            save(newFilePath,'-struct','newBatchStruct')
            
            % Update UI
            S.feedbackWin.String = 'Merge finished';
            drawnow
            
        end
        
    end


end
