function [filesTable, canceled] = movie_selector(filesTable,startingPath)


%Ui based tool where the user can select multiple tiff files. The list of
%files is used to create TrackIts batch structure (see init_batch.m).
%
%[filesTable, canceled] = movie_selector(filesTable,startingPath)
%
%Input:
%   filesTable      -   Existing filesTable passed by TrackIt main gui.
%                       Table with 3 columns:
%                           FileName: Char array containing the filename of the tiff
%                           PathName: Char array containing the path of of the file
%                           frameCycleTime: the frame cycle time of the
%                           movie in miliseconds.
%
%   startingPath	-   Path that is displayed when the file or folder
%                       dialog is opened
%   
%Output:
%
%
%   filesTable      -   New filesTable passed by TrackIt main gui.
%                       Table with 3 columns:
%                           FileName: Char array containing the filename of the tiff
%                           PathName: Char array containing the path of of the file
%                           frameCycleTime: the frame cycle time of the
%                           movie in miliseconds.
%
%
%   canceled        -   Whether user closed the window without pressing "OK"
%
%


%Create ui
S = CreateMovieSelectorUI();

%Get frame cycle times and number of movies in each frame cycle time
S.tb.Data = getTimelapseTable(filesTable.frameCycleTime);


uiwait(S.f);

    function S = CreateMovieSelectorUI()
        
        buttonPosHor = .75;
        buttonWidth = .22;
        buttonHeight = .075;
        
        S.f = figure('Units','normalized',...
            'Position',[0.1 0.5 .3 .35],...
            'MenuBar','None',...
            'Name','Movie selector',...
            'NumberTitle','off',...
            'WindowKeyPressFcn',@KeyPressFcnCB,...
            'CloseRequestFcn',@CloseRequestCB);
        S.tb = uitable('Units','normalized',...
            'Position',[0.02 0.63 .7 .35],...
            'ColumnName',{'#Files';'Frame cycle time (tl condition) [ms]'},...
            'ColumnEditable',[false,true],...
            'CellEditCallback',@TableEditCB,...
            'CellSelectionCallback',@TableSelectionCB);
        S.lb = uicontrol('Style','Listbox'...
            ,'Units','normalized',...
            'Position',[0.02 0.0 .7 .6],...
            'Max',2,'Min',0);
        
        filesButton = uicontrol('String','Add movies',...
            'Units','normalized',...
            'Position',[buttonPosHor .9 buttonWidth buttonHeight],...
            'Callback',@AddFilesCB);
        
        
        textSearchString = uicontrol('String','Search string for folder scan (optional):',...
            'Style','Text',...
            'Units','normalized',...
            'Position',[buttonPosHor  .77  buttonWidth  .08]); 
        
        S.editSearchString = uicontrol('String','',...
            'Style','Edit',...
            'Units','normalized',...
            'Position',[buttonPosHor .72 buttonWidth buttonHeight-.025]);               
        
        folderButton = uicontrol('String','Scan folder for movies',...
            'Units','normalized',...
            'Position',[buttonPosHor  .63  buttonWidth buttonHeight],...
            'Callback',@AddFolderCB);
        
        
        removeAllButton = uicontrol('String','Remove all files',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.515 buttonWidth buttonHeight],...
            'Callback',@RemoveAllFilesCB);
        removeButton = uicontrol('String','Remove selected files',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.42 buttonWidth buttonHeight],...
            'Callback',@RemoveFilesCB);
                
        S.bgSortMovies = uibuttongroup('Title','Sort movies by',...
            'Units','normalized',...
            'Position',[buttonPosHor  .1  buttonWidth  .3]); 

        S.rbFileName	= uicontrol(S.bgSortMovies,'Units','normalized',...
            'Position',[.05 .75 .99 .25],...
            'Style','radiobutton',...
            'String','Filename');
        S.rbTlCond = uicontrol(S.bgSortMovies,'Units','normalized',...
            'Position',[.05 .525 .99 .25],...
            'Style','radiobutton',...
            'String','Frame cycle time');
        
         S.rbStrPattern = uicontrol(S.bgSortMovies,'Units','normalized',...
            'Position',[.05 .3 .99 .25],...
            'Style','radiobutton',...
            'String','String pattern:');
        
        S.editStrPattern = uicontrol(S.bgSortMovies,'Style','Edit',...
            'Units','normalized',...
            'Position',[.2 .1 .7 .2],...
            'String','');        
                
        okButton = uicontrol('String','OK',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.01 buttonWidth buttonHeight],...
            'Callback',@OkCB);
        
    end

    function CloseRequestCB(~,~)
        canceled = 1;
        delete(gcf)
    end

    function KeyPressFcnCB(~,event)
        %Close figure if user pressed escape
        if strcmp(event.Key, 'escape')
            canceled = 1;
            delete(gcf)
        end
    end

    function OkCB(~,~)
        %User pressed "Ok" button
        
        canceled = 0;
        if S.rbFileName.Value
            %Sort movie by their filename
            filesTable = sortrows(filesTable,1);
        elseif S.rbTlCond.Value
            %Sort movies by their frame cycle time
            filesTable = sortrows(filesTable,3);
        else
            %Sort movies according to a string entered by the user
            
            %Get cell array of filenames
            fileNames = filesTable.FileName;
            
            %Find filenames which contain the string pattern and return the
            %matching part of the filenames
            [matchResults, ~] = regexp(fileNames, [S.editStrPattern.String, '.*'],'match','split');
            
            %Catenate only files where a match was found
            matchedFileParts = [matchResults{:}];
            
            %Sort found findparts
            [~,sortingIdx] = sort(matchedFileParts);
            
            %Get indices of files which match the string pattern
            matchedIdx = ~cellfun(@isempty,matchResults);
     
            %Create files table containing the sorted matched filenames
            matchedFilesTable = filesTable(matchedIdx, :);
            
            %Create files table containing the filenames where no match was
            %found
            nonMatchedFilesTable = filesTable(~matchedIdx, :);

            %Catenate matched and non-matched filenames
            filesTable = [matchedFilesTable(sortingIdx,:); nonMatchedFilesTable];
        end

        delete(gcf)
    end

    function AddFilesCB(~,~)
        %User pressed "Add files" button

        %Open file dialog box
        [fileNameListNew,pathName] = uigetfile({'*.tif*'},'Select files to track spots', 'MultiSelect', 'on',startingPath);
        
        if isequal(fileNameListNew,0) %User didn't choose a file
            return
        elseif ~iscell(fileNameListNew) %Check if only one file has been chosen
            fileNameListNew = {fileNameListNew};
        end
        
        %Save path for next usgage of file dialog box
        startingPath = pathName;
        
        %Create list of filepaths
        pathNameListNew = repmat({pathName},1,length(fileNameListNew));
        
        %Get list of frame cycle times found in the filenames
        tlListNew = getTimelapseList(fileNameListNew);
        
        %Create table containing new filenames, paths and frame cycle times
        filesTableNew = table(fileNameListNew',pathNameListNew',tlListNew,'VariableNames',{'FileName','PathName','frameCycleTime'});
        
        %Find movies that are already part of fileTable
        doubleIdx = false(height(filesTableNew),1);
        
        for fileIdx = 1:height(filesTableNew)
            doubleIdx(fileIdx) = any(strcmp(filesTableNew.FileName{fileIdx}, filesTable.FileName(:)));
        end
        
        if any(doubleIdx)
            %Movies which are already in the list where found, so inform
            %the user which movies will be ignored
            msg = ['The following movies were already added and will be ignored'; filesTableNew.FileName(doubleIdx)];
            
            nFiles = sum(doubleIdx);
            if nFiles > 10
                msg = [msg(1:11,:); ['Plus ', num2str(nFiles-10), ' additional movies']];
            end
            msgbox(msg);
            
            %delete the relevant files from the files table
            filesTableNew(doubleIdx,:) = [];    
        end
        
        %Catenate new filesTable with existing filesTable
        filesTable = [filesTable;filesTableNew];
        
        %Update timelapse list
        S.tb.Data = getTimelapseTable(filesTable.frameCycleTime);
    end

    function AddFolderCB(~,~)
        %User pressed "Search folder" button
        
        %Open folder selection dialog
        parentFolder = uigetdir(startingPath,'Choose a Folder');
        
        if parentFolder == 0
            return
        end
        
        %Check if user entered a specific string to search for in the folder
        if isempty(S.editSearchString.String)
            %No string entered so search for all .tif or .tiff files
            searchPattern = '**\*.tif*';
        else
            %String pattern was entered so define a specific search pattern
            searchPattern = ['**\*',S.editSearchString.String,'*.tif*'];
        end
        
        %Save path for next usgage of a path dialog box
        startingPath = parentFolder;
        
        %Get list of files matching the search pattern
        fileList = dir(fullfile(parentFolder,searchPattern));
        
        %Create list of filenames
        fileNameList = {fileList.name};
        
        
        %Get list of frame cycle times found in the filenames
        tlListNew = getTimelapseList(fileNameList);
        
        %Create new files table from the new files
        filesTableNew = table(fileNameList',{fileList.folder}',tlListNew,'VariableNames',{'FileName','PathName','frameCycleTime'});
        
         %Find movies that are already part of fileTable
        doubleIdx = false(height(filesTableNew),1);
        
        for fileIdx = 1:height(filesTableNew)
            doubleIdx(fileIdx) = any(strcmp(filesTableNew.FileName{fileIdx}, filesTable.FileName(:)));
        end
                
        if any(doubleIdx)
            %Movies which are already in the list where found, so inform
            %the user which movies will be ignored
            msg = ['The following movies were already added and will be ignored'; filesTableNew.FileName(doubleIdx)];
            
            nFiles = sum(doubleIdx);
            if nFiles > 10
                msg = [msg(1:11,:); ['Plus ', num2str(nFiles-10), ' additional movies']];
            end
            msgbox(msg);
            
            %delete the relevant files from the files table
            filesTableNew(doubleIdx,:) = [];    
        end
        
        %Catenate new filesTable with existing filesTable
        filesTable = [filesTable;filesTableNew];
        
        %Update list of frame cycle times
        S.tb.Data = getTimelapseTable(filesTable.frameCycleTime);
    end

    function RemoveAllFilesCB(~,~)
        %User pressed "Remove all files" button so delete all entries from
        %the files table and update the ui
        
        filesTable(:,:) = [];
        S.tb.Data = getTimelapseTable(filesTable.frameCycleTime);
        S.lb.String = filesTable.FileName(filesTable.frameCycleTime == S.tb.UserData);
    end

    function RemoveFilesCB(~,~)
        %User pressed "remove selected files" button
        
        %Check if any files are selected
        if any(S.tb.UserData)
            %Find indices of filenames corresponding to currently selected 
            %frame cycle time. S.tb.UserData contains the currently
            %selected frame cycle time and was set in TableSelectionCB
            %function.
            currentListboxList = filesTable.frameCycleTime == S.tb.UserData;
            
            %Convert logical array to indices
            currentListBoxInTable = find(currentListboxList);
            
            %Find selected files in currentListBoxInTable and delete these
            %from the filesTable
            filesTable(currentListBoxInTable(S.lb.Value),:) = [];
            
            %Update frame cycle time list
            S.tb.Data = getTimelapseTable(filesTable.frameCycleTime);
            
            %Update list of filenames
            S.lb.String = filesTable.FileName(filesTable.frameCycleTime == S.tb.UserData);

            if size(S.lb.String,1) == 1
                %If only one filename left, set selected field to first
                %entry and disable multiple selections
                S.lb.Value = 1;
                S.lb.Max = 0;
            elseif size(S.lb.String,1) < S.lb.Value(end)
                %Make sure selected field number is not higher than entries
                %in the listbox so select first field in list
                S.lb.Value = 1;
            end
        end
    end

    function TableSelectionCB(tableHandle,selectedCell)
        %User selected a specific frame cycle time in the table so display
        %the filenames corresponing to this frame cycle time in the listbox
        
        %Set seleted file to first entry
        S.lb.Value = 1;
        if any(selectedCell.Indices)
            %Display the filenames corresponding to the selected frame
            %cycle time in the listbox
            S.lb.String = filesTable.FileName(filesTable.frameCycleTime == tableHandle.Data(selectedCell.Indices(1), 2));
            
            %Save currently selected frame cycle time in the user data property so we
            %know which files to delete if the user presses "delete selected files"
            S.tb.UserData = tableHandle.Data(selectedCell.Indices(1), 2);             
            if size(S.lb.String,1) > 1
                %More than one filename exisits so allow multiple selection
                S.lb.Max = 2;
            else
                %Less than two filenames exisits so disable multiple selection
                S.lb.Max = 0;
            end
        end
    end

    function TableEditCB(~,editedData)
        %User changed the frame cycle time so update the frame cycle time
        %in filesTable
        filesTable.frameCycleTime(filesTable.frameCycleTime == editedData.PreviousData) = str2double(editedData.EditData);
        
        %Update the frame cycle times and number of movies in each frame cycle time
        S.tb.Data = getTimelapseTable(filesTable.frameCycleTime);
    end

    function tlList = getTimelapseList(fileNameList)
        %Check list of filenames wether they contain either "s" (seconds), 
        %"ms" (miliseconds), "Hz" (Hertz) or "_t" (tens of miliseconds) and
        %return a list of frame cycle times
        
        tlList = ones(length(fileNameList),1).*-1;
        
        %Iterate through filenames
        for i = 1:length(fileNameList)
            curTl1 = regexp(fileNameList{i}, '_\d+ms','match'); %Search for _+Number+ms pattern
            curTl2 = regexp(fileNameList{i}, '_\d+s','match');  %Search for _+Number+s pattern
            curTl3 = regexp(fileNameList{i}, '_t\d+','match');  %Search for _t+Number pattern (This convention is in tens of miliseconds!)
            curTl4 = regexp(fileNameList{i}, '_\d+Hz','match');  %Search for _+Number+Hz pattern
            
            if ~isempty(curTl1)
                %Filename contains "ms"
                curTl1 = regexp(curTl1{1}, '\d+','match');
                tlList(i) = str2double(curTl1{1}(1:end));
            elseif ~isempty(curTl2)
                %Filename contains "s"
                curTl2 = regexp(curTl2{1}, '\d+','match');
                tlList(i) = str2double(curTl2{1}(1:end))*1000;
            elseif ~isempty(curTl4)
                %Filename contains "Hz"
                curTl4 = regexp(curTl4{1}, '\d+','match');
                tlList(i) = 1000/str2double(curTl4{1}(1:end));
            elseif ~isempty(curTl3)
                %Filename contains "_t"
                curTl3 = regexp(curTl3{1}, '\d+','match');
                tlList(i) = str2double(curTl3{1}(1:end))*10;
            end
        end
        
    end

    function tlTable = getTimelapseTable(tlList)
        %Get frame cycle times and number of movies in each frame cycle time
        
        %Get unique set of frame cycle times
        tluni = unique(tlList);
        
        %Iterate through entries and count occurences of each frame cycle time
        N = zeros(length(tluni),1);        
        for i = 1:length(tluni)
            N(i) = sum(tlList == tluni(i));
        end
        
        %Save amount and frame cycle time to the return variable
        tlTable = [N,tluni];
        
    end

end
