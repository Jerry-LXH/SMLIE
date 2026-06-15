function movie_splitter(varargin)

if nargin == 1
    startingPath = varargin{1};
else
    startingPath = pwd;
end

%Initialize table with filenames and paths
filesTable = cell2table(cell(0,2));
filesTable.Properties.VariableNames = {'FileName','PathName'};

%Create user interface
S = CreateMovieSelectorUI();

    function S = CreateMovieSelectorUI()
        
        S.f = figure('Units','normalized',...
            'Position',[0.5 0.5 .35 .3],...
            'MenuBar','None',...
            'Name','Movie splitter',...
            'NumberTitle','off',...
            'WindowKeyPressFcn',@KeyPressFcnCB,...
            'CloseRequestFcn',@CloseRequestCB);
        S.lb = uicontrol('Style','Listbox'...
            ,'Units','normalized',...
            'Value', [],...
            'Min', 0, 'Max', 2,...
            'Position',[0.05 0.38 .7 .6]);
        S.tb = uitable('Units','normalized',...
            'Position',[0.05 0.01 .7 .35],...
            'ColumnName',{'#frames in sequence';'Name';'Create .tiff file?'},...
            'Data',{'1','seq1',true;'1','seq2',true},...
            'ColumnEditable',[true,true,true]);
        
        buttonPosHor = .76;
        buttonWidth = .22;
        buttonHeight = .1;
        
        filesButton = uicontrol('String','Select files',...
            'Units','normalized',...
            'Position',[buttonPosHor .88 buttonWidth buttonHeight],...
            'Callback',@AddFilesCB);
        removeButton = uicontrol('String','Remove selected files',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.75 buttonWidth buttonHeight],...
            'Callback',@RemoveFilesCB);
        
        S.cboxCreateMetadata = uicontrol('String','<html>Create .txt file containing metadata of original movie',...
            'Units','normalized','Style','checkbox',...
            'Value',1,...
            'Position',[buttonPosHor .55  buttonWidth .2]);
        nSequencesText = uicontrol('String','Amount of splits',...
            'Units','normalized','Style','Text',...
            'Position',[buttonPosHor .45 buttonWidth buttonHeight],...
            'Callback',@AddFilesCB);
        nSequencesEdit = uicontrol('String','2',...
            'Units','normalized','Style','Edit',...
            'Position',[.82 .4 .1 .08],...
            'Callback',@NSplitsCB);
        
        S.editFeedbackWin = uicontrol('String','',...
            'Units','normalized',...
            'Position',[buttonPosHor 0.15 buttonWidth .2],...
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

    function NSplitsCB(src,~)
        %Executed when "Amount of splits" is changed
        
        %Get current table data
        data = S.tb.Data;
        
        %Get amount of desired splits
        nSequences = str2double(src.String);
        
        %Calculate difference between current amount of splits and new
        %amount of splits
        nAdditional = nSequences - size(data,1);
        
        if nAdditional > 0
            %More splits required
            
            %Create vector containing numbers of new rows
            newRowsVec = size(data,1)+1:size(data,1)+nAdditional;
            
            %Create cell array containing sequence strings for display
            strArray = cell(numel(newRowsVec),1);
            for m = 1:numel(newRowsVec)
                strArray{m} = ['seq',num2str(newRowsVec(m))];
            end
            
            %Create additional table entries
            S.tb.Data = [data; repmat({'1','',true},nAdditional,1)];
            %Write sequence strings into "Name" column
            S.tb.Data(size(data,1)+1:end,2) = strArray;
        elseif nAdditional < 0
            %Less splits required so just delete amount of excess rows
            S.tb.Data = data(1:end+nAdditional,:);
        end
    end

    function AddFilesCB(~,~)
        %User pressed "Add files" button
        
        %Open file dialog box
        [fileNameListNew,pathName] = uigetfile({'*.tif*'},'Select files you want to split', 'MultiSelect', 'on',startingPath);
        
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

    function StartCB(~,~)
        %User pressed "Start" button
        
        %Get table containing information on how to split
        splitInfo = S.tb.Data;
        
        %Get information if .txt files with metadata should be created
        createTxt = S.cboxCreateMetadata.Value;
        
        if ~isempty(filesTable) && ~isempty(splitInfo)
            
            %Create array containing the sequences
            sequences = cellfun(@str2double,splitInfo(:,1));
            
            %Get amount of splits
            nSequences = numel(sequences);
            
            %Create cell array containing the add-ons to original filenames
            nameArray = splitInfo(:,2);
            
            %Additionally we create a cell array with unique sequence name
            %entries so we now which sequences belong together because they
            %share the same name.
            [uniqueNameArray,ia,~] = unique(nameArray,'stable');
            
            %Create array containing information if a .tiff file should be
            %saved
            createTiff = splitInfo(ia,3);
            
            %Get amount of unique sequences
            nUniqueSequences = length(uniqueNameArray);
            
            %Get number of files to split
            nFiles = height(filesTable);
            
            %Iterate through files
            for n = 1:nFiles
                %Monitor progress in ui
                S.editFeedbackWin.String = ['Loading Movie ' num2str(n) ' of ' num2str(nFiles)];
                drawnow
                
                %Get current filename and pathname
                curFileName = filesTable.FileName{n};
                curPathName = filesTable.PathName{n};
                
                %----------Write text file with metadata-----------------------
                if createTxt
                    %Retreive metadata using bioformats
                    metaData = read_tiff_metadata(fullfile(curPathName,curFileName));
                    
                    %Split metadata at commas and sort them
                    sortedMetaData = sort(strsplit(char(metaData{1}),','));
                    
                    %Create filename for text file
                    [~,curFileNamePart] = fileparts(curFileName);
                    txtFileName = fullfile(curPathName,strcat(curFileNamePart,'.txt'));
                    
                    %Open a text file 
                    fid = fopen(txtFileName,'w');
                    %Write metadata to text file
                    for r=1:size(sortedMetaData,1)
                        fprintf(fid,'%s\r\n',sortedMetaData{r,:});
                    end
                    %Close text file
                    fclose(fid);
                end
                
                %----------Load Stack--------------------------------------
                original = load_stack(curPathName, curFileName, S);
                
                %--------Divide movie into substacks---------------------
                S.editFeedbackWin.String = ['Splitting Movie ' num2str(n) ' of ' num2str(nFiles)];
                drawnow
                
                %Initialize variable indicating current substack number
                curSubStackNum = 1;
                
                %Get satck size
                oriStackSize = size(original);
                
                %Get bit depth of original movie
                oriClass = class(original);
                
                %Initialize cell array to contain one substack pero cell
                substacks = repmat({zeros(oriStackSize(1),oriStackSize(2),0,oriClass)},nUniqueSequences,1);
                
                %Initialize frame counter indicating amount of frames in
                %current sequence
                curSequenceFrameCounter = 1;
                
                %Create variable containing information on how many frames
                %the current sequence should be containing. Start with
                %first sequence.                
                nFramesInCurSeq = sequences(1);
                
                %Get amount of frames in original stack
                nFrames = oriStackSize(3);
                
                feedbackWin = S.editFeedbackWin.String;
                
                curSequenceStack = zeros(oriStackSize(1),oriStackSize(2),0); 
                
                %Iterate through frames
                for m = 1:nFrames
                    %Write current frame to corresponding substack
                    curSequenceStack(:,:,curSequenceFrameCounter) = original(:,:,m); 
                    
                    %Check if current frame should be in the next sequence
                    %or if the end of the original movie is reached
                    if curSequenceFrameCounter == nFramesInCurSeq || m == nFrames
                        %Get name of current Sequence
                        curName = nameArray{curSubStackNum};
                        %Find index of corresponding substack
                        index = find(strcmp(curName,uniqueNameArray));
                        %Catenate current sequence into substack
                        substacks{index} = cat(3,substacks{index},curSequenceStack(:,:,1:curSequenceFrameCounter));
                        
                        %Go to next sequence
                        curSubStackNum = mod(curSubStackNum, nSequences)+1;
                        %Get number of frames in next sequence
                        nFramesInCurSeq = sequences(curSubStackNum);
                        %Reinitialize current sequence
                        curSequenceStack = zeros(oriStackSize(1),oriStackSize(2),nFramesInCurSeq);
                        %Reinitialize frame counter for current Sequence
                        curSequenceFrameCounter = 0;
                    end
                    
                    
                    
                    %Increase counter indicating amount of frames in
                    %current stack
                    curSequenceFrameCounter = curSequenceFrameCounter + 1;
                    
                    %Monitor progress in ui
                    percentDone = round(m * 100/ nFrames);                    
                    if mod(percentDone,5) == 0
                        S.editFeedbackWin.String = char(sprintf('Splitting progress: %3.0f %%', percentDone), feedbackWin);
                        drawnow
                        if double(get(gcf,'CurrentCharacter')) == 27
                            t.close();
                            break
                        end
                    end
                    
                end
                
                %----------Get tiff tags from original movie-------------
                fullFileOriginal = char(fullfile(curPathName,curFileName));
                
                warning('off'); %Supress warnings for unrecognized tif tags
                TifLink = Tiff(fullFileOriginal, 'r');
                
                %Required .tiff tags
                tagstruct.ImageWidth = getTag(TifLink,'ImageWidth');
                tagstruct.ImageLength = getTag(TifLink,'ImageLength');
                tagstruct.BitsPerSample = getTag(TifLink,'BitsPerSample');
                tagstruct.SamplesPerPixel = getTag(TifLink,'SamplesPerPixel');
                tagstruct.Compression = getTag(TifLink,'Compression');
                tagstruct.PlanarConfiguration = getTag(TifLink,'PlanarConfiguration');
                tagstruct.Photometric = getTag(TifLink,'Photometric');
                
                %Additional .tiff tags
                tagstruct.RowsPerStrip = getTag(TifLink,'RowsPerStrip');
                tagstruct.Orientation = getTag(TifLink,'Orientation');
                tagstruct.SampleFormat = getTag(TifLink,'SampleFormat');
                close(TifLink);
                
                warning('on');

                %-------------Write .tiff files------------------------
                for m = 1:nUniqueSequences
                    if createTiff{m} %User chose 'save .tiff file' for this sequence
                        
                        %Get original filename
                        [~,fileWithoutExt,~] = fileparts(curFileName);
                        
                        %Get add-on to original filename
                        curFileAddon = splitInfo{m,2};
                        
                        %Create new filename
                        newFullFilename = char(fullfile(curPathName, strcat(fileWithoutExt,'_', curFileAddon, '.tiff')));
                        
                        %Open tiff file for writing
                        t = Tiff(newFullFilename,'w');
                        
                        %Set tiff tiags
                        t.setTag(tagstruct);
                        
                        %Write first frame
                        t.write(substacks{m}(:,:,1));
                        
                        %Get amount of frames
                        nFrames = size(substacks{m},3);
                        
                        %Iterate through the rest of frames
                        for k=2:nFrames
                            %Write current frame and tags
                            t.writeDirectory();
                            t.setTag(tagstruct);
                            t.write(substacks{m}(:,:,k));
                            
                            %Monitor progress
                            percentDone = round(k * 100/ nFrames);
                            if mod(percentDone,5) == 0
                                S.editFeedbackWin.String = char(sprintf('Writing Sequence %3.0f: %3.0f %%',m, percentDone), feedbackWin);
                                drawnow
                                if double(get(gcf,'CurrentCharacter')) == 27
                                    t.close();
                                    break
                                end
                            end
                            
                        end
                        t.close();
                    end
                end
                
            end
            S.editFeedbackWin.String = char('Splitting Finished');
        end
    end

end
