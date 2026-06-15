function [NPC_ShiftX,NPC_ShiftY] = NP_Cloud_NoSort(ArRefXY, ArCmpXY, TolInPixel)  
%Call function to calculate the relative shift between two set of SMLM data through NP-Cloud

%---IMPORTANT: Both RefXY and CmpXY need to be already sorted by Y! Otherwise uncomment the below two lines to sort
%    ArRefXY = sortrows(ArRefXY,2);
%    ArCmpXY = sortrows(ArCmpXY,2);
%----------------------------------

%---Inputs------
%ArRefXY: Array storing the X,Y positions of the reference SMLM data. See note above: Already sorted by Y
%ArCmpXY: Array storing the X,Y positions of the comparing SMLM data. See note above: Already sorted by Y
%TolInPixel: Search radius
%The above parameters are all in the unit of pixel
%---------------

%---Outputs------
%NPC_ShiftX,NPC_ShiftY: The X,Y shift values between the comparing and reference SMLM data through the NP-Cloud calculation
%---------------

%Note: In the notations below, "localizations" are sometimes referred to as "molecules". We apologize for this inconsistency.

    PositionCloudSizeLimit=uint32(90000000);  %Preset an upper limit of cloud size    
    
    TolInit=TolInPixel*1.25; %Initial processing uses a larger initial search radius for later refinement
    TolInitSq=TolInit*TolInit; %Square of initial search radius for easier comparison
    TolSq=TolInPixel*TolInPixel; %Square of search radius for easier comparison
    
    DistanceCloudArrayX=zeros(PositionCloudSizeLimit,1,'single'); %DistanceCloudArrayX and DistanceCloudArrayY store the cloud of displacements
    DistanceCloudArrayY=zeros(PositionCloudSizeLimit,1,'single');
    
    PointerSearchYLow=uint32(1); %Pointers in the reference data: with presorted Y, we only need to scan from top to bottom in increasing Y.
    PointerSearchYHi=uint32(1);
    CountInCloud=uint32(0);
    NumMoleculeWithMatches=uint32(0);

    ArRefX=ArRefXY(:,1);
    ArRefY=ArRefXY(:,2);
    ArRefY(end+1)= inf; %Y is patched at the end to cap the final value

    ArCmpX=ArCmpXY(:,1);
    ArCmpY=ArCmpXY(:,2);
    MoleNumCmp=uint32(length(ArCmpX));
        
    MoleculeMatchPositionArrary=zeros(PositionCloudSizeLimit/2,1,'uint32');
    
    for CurrentMoleInCmp=uint32(1):MoleNumCmp %Go through all the localizations in the comparing SMLM data to pair with the localizations the reference data

        CurrentCmpY=ArCmpY(CurrentMoleInCmp);
        CmpYLow=CurrentCmpY-TolInit;
        CmpYHi=CurrentCmpY+TolInit;

        while(CmpYLow>ArRefY(PointerSearchYLow)) %Lower pointer for a y range within the initial search radius 
            PointerSearchYLow=PointerSearchYLow+1;
        end

        while(CmpYHi>ArRefY(PointerSearchYHi)) %Upper pointer for a y range within the initial search radius
            PointerSearchYHi=PointerSearchYHi+1;            
        end            

        FoundY=PointerSearchYLow:PointerSearchYHi-1; %-1 as the search has just passed

        if ~isempty(FoundY) %If localizaitons are found within the y range

            CurrentCmpX=ArCmpX(CurrentMoleInCmp);

            ArMatchedX=ArRefX(FoundY);
            ArMatchedY=ArRefY(FoundY);
            ArDistanceX=CurrentCmpX-ArMatchedX;
            ArDistanceY=CurrentCmpY-ArMatchedY;

            IndexMatchingX=find(abs(ArDistanceX) < TolInit); %Find positions with X range also within the initial search radius

            if ~isempty(IndexMatchingX)

                UsefulDeltaX=ArDistanceX(IndexMatchingX);
                UsefulDeltaY=ArDistanceY(IndexMatchingX);

                ArDistanceSq=UsefulDeltaX.*UsefulDeltaX+UsefulDeltaY.*UsefulDeltaY; %Calculate the distance square
                IndexWithinDistance=find(ArDistanceSq<TolInitSq); %Find those within the square of initial search radius

                if ~isempty(IndexWithinDistance)

                    NumMatched=length(IndexWithinDistance);

                    StartRecordPos=CountInCloud+1;                        
                    CountInCloud=CountInCloud+NumMatched;

                    NumMoleculeWithMatches=NumMoleculeWithMatches+1;
                    MoleculeMatchPositionArrary(NumMoleculeWithMatches)=CountInCloud; %Mark the end of matched point of each molecule in the array

                    DistanceCloudArrayX(StartRecordPos:CountInCloud)= UsefulDeltaX(IndexWithinDistance); %Add all matched displacements to the cloud
                    DistanceCloudArrayY(StartRecordPos:CountInCloud)= UsefulDeltaY(IndexWithinDistance);

                end
            end
        end 
    end


    CloudX=DistanceCloudArrayX(1:CountInCloud);
    CloudY=DistanceCloudArrayY(1:CountInCloud);

    Current_centerX=mean(CloudX);
    Current_centerY=mean(CloudY);
    CurrentCenterShiftSq=Current_centerX*Current_centerX+Current_centerY*Current_centerY;
    PreviousShiftSq=100;

    Current_ShiftX=single(0);
    Current_ShiftY=single(0);


    while (CurrentCenterShiftSq<PreviousShiftSq) %Continue to shift the center until converges, i.e., the new mean position of the cloud is no longer closer to the origin than the previous round

        PreviousShiftSq=CurrentCenterShiftSq;
        Current_centerX=Current_centerX+Current_ShiftX; %Accumulative shift of the averaged center of the cloud
        Current_centerY=Current_centerY+Current_ShiftY;

        CurrentPassX=CloudX-Current_centerX;
        CurrentPassY=CloudY-Current_centerY;

        CurrentPassRSq=CurrentPassX.*CurrentPassX + CurrentPassY.*CurrentPassY;

        CurrentSel=(CurrentPassRSq<TolSq); %Find the displacements that are within the precise search radius

        CurrentPassXSel=CurrentPassX(CurrentSel);
        CurrentPassYSel=CurrentPassY(CurrentSel);

        Current_ShiftX=mean(CurrentPassXSel);
        Current_ShiftY=mean(CurrentPassYSel);


        CurrentCenterShiftSq = Current_ShiftX*Current_ShiftX + Current_ShiftY*Current_ShiftY;

    end

    %--After converged based on all in range, now shift to that rough center and search again using only the nearest neighbors

    CurrentSelPoints=find(CurrentSel); %These initial cloud points are used.
    CurrentSelPoints(end+1)=inf; %Patch at the end for scanning.

    ScanMatchMolecule=uint32(1); %For scanning in the initial cloud points with multiple matches.
    NumNearestMole=uint32(0);
    CurrentPassSelSq=CurrentPassXSel.*CurrentPassXSel + CurrentPassYSel.*CurrentPassYSel;

    for CurrentMatchMoleucle=1:NumMoleculeWithMatches %Find nearest neighbor for each molecule

        NearestDistanceSq=inf;
        CurrentCloudRangeUpper=MoleculeMatchPositionArrary(CurrentMatchMoleucle);

        while(CurrentSelPoints(ScanMatchMolecule)<=CurrentCloudRangeUpper) %For all matched molecules of the current molecule
            CmpDistanceSq=CurrentPassSelSq(ScanMatchMolecule);

            if (CmpDistanceSq<NearestDistanceSq) %Record the nearest pair
                RecordMole=ScanMatchMolecule;
                NearestDistanceSq=CmpDistanceSq;
            end

            ScanMatchMolecule=ScanMatchMolecule+1;

        end

        if isfinite(NearestDistanceSq) %If the current molecule has found a nearest pair, since started each molecule as Inf.
            NumNearestMole=NumNearestMole+1;
            CurrentPassXSel(NumNearestMole)=CurrentPassXSel(RecordMole); %Reusing CurrentPassXSel to store each molecule - will not overwrite as NumNearestMole <= RecordMole;
            CurrentPassYSel(NumNearestMole)=CurrentPassYSel(RecordMole);
        end    				

    end

    NPC_ShiftX=mean(CurrentPassXSel(1:NumNearestMole))+Current_centerX; %Total shift between the two channels is the shift in the final refining step plus the accumulated shifts above in the above iterations
    NPC_ShiftY=mean(CurrentPassYSel(1:NumNearestMole))+Current_centerY;