function [NPC_ShiftX,NPC_ShiftY,NPC_ShiftZ] = NP_Cloud_NoSort_3D_Pass1(ArRefXYZ, ArCmpXYZ, TolXYInPixel, TolZnm)  
%Call function to calculate the relative shift between two set of 3D-SMLM data through NP-Cloud-3D
%For Pass 1, relative z shift is calculated as the Z position of the cloud center by simply averaging all Z values of the cloud

%---IMPORTANT: Both ArRefXYZ and ArCmpXYZ need to be already sorted by Y! Otherwise uncomment the below two lines to sort
%    ArRefXYZ = sortrows(ArRefXYZ,2);
%    ArCmpXYZ = sortrows(ArCmpXYZ,2);
%----------------------------------

%---Inputs------
%ArRefXYZ: Array storing the X,Y,Z positions of the reference SMLM data. See note above: Already sorted by Y
%ArCmpXYZ: Array storing the X,Y,Z positions of the comparing SMLM data. See note above: Already sorted by Y
%TolXYInPixel: In-plane (XY) search radius. Unit is pixel
%TolZnm: Z search radius in the unit of nm.
%---------------

%---Outputs------
%NPC_ShiftX,NPC_ShiftY,NPC_ShiftZ: The X,Y,Z shift values between the comparing and reference SMLM data through the NP-Cloud-3D calculation
%---------------

%Note: In the notation below, "localizations" are sometimes referred to as "molecules". We apologize for this inconsistency.

    PositionCloudSizeLimit=uint32(90000000);  %Preset an upper limit of cloud size    
    
    TolXYInit=TolXYInPixel*1.25; %Initial processing uses a larger initial search radius for later refinement
    TolXYInitSq=TolXYInit*TolXYInit; %Square of initial search radius for easier comparison
    TolXYSq=TolXYInPixel*TolXYInPixel; %Square of search radius for easier comparison
    
    DistanceCloudArrayX=zeros(PositionCloudSizeLimit,1,'single'); %DistanceCloudArrayX, DistanceCloudArrayY, DistanceCloudArrayZ store the cloud of displacements
    DistanceCloudArrayY=zeros(PositionCloudSizeLimit,1,'single');
    DistanceCloudArrayZ=zeros(PositionCloudSizeLimit,1,'single');
    
    PointerSearchYLow=uint32(1); %Pointers in the reference data: with presorted Y, we only need to scan from top to bottom in increasing Y.
    PointerSearchYHi=uint32(1);
    CountInCloud=uint32(0);
    NumMoleculeWithMatches=uint32(0);

    ArRefX=ArRefXYZ(:,1);
    ArRefY=ArRefXYZ(:,2);
    ArRefZ=ArRefXYZ(:,3);
    ArRefY(end+1)= inf; %Y is patched at the end to cap the final value

    ArCmpX=ArCmpXYZ(:,1);
    ArCmpY=ArCmpXYZ(:,2);
    ArCmpZ=ArCmpXYZ(:,3);
    MoleNumCmp=uint32(length(ArCmpX));
        
    MoleculeMatchPositionArrary=zeros(PositionCloudSizeLimit/2,1,'uint32');
    
    for CurrentMoleInCmp=uint32(1):MoleNumCmp %Go through all the localizations in the comparing SMLM data to pair with the localizations the reference data

        CurrentCmpY=ArCmpY(CurrentMoleInCmp);
        CmpYLow=CurrentCmpY-TolXYInit;
        CmpYHi=CurrentCmpY+TolXYInit;

        while(CmpYLow>ArRefY(PointerSearchYLow)) %Lower pointer for a y range within the initial search radius 
            PointerSearchYLow=PointerSearchYLow+1;
        end

        while(CmpYHi>ArRefY(PointerSearchYHi)) %Upper pointer for a y range within the initial search radius
            PointerSearchYHi=PointerSearchYHi+1;            
        end            

        FoundY=PointerSearchYLow:PointerSearchYHi-1; %-1 as the search has just passed

        if ~isempty(FoundY) %If localizaitons are found within the y range

            CurrentCmpX=ArCmpX(CurrentMoleInCmp);
            CurrentCmpZ=ArCmpZ(CurrentMoleInCmp);

            ArMatchedX=ArRefX(FoundY);
            ArMatchedY=ArRefY(FoundY);
            ArMatchedZ=ArRefZ(FoundY);
            ArDistanceX=CurrentCmpX-ArMatchedX;
            ArDistanceY=CurrentCmpY-ArMatchedY;
            ArDistanceZ=CurrentCmpZ-ArMatchedZ;

            IndexMatchingX=find(abs(ArDistanceX) < TolXYInit & abs(ArDistanceZ) < TolZnm); %Find positions with X & Z ranges also within the initial search radius

            if ~isempty(IndexMatchingX)

                UsefulDeltaX=ArDistanceX(IndexMatchingX);
                UsefulDeltaY=ArDistanceY(IndexMatchingX);
                UsefulDeltaZ=ArDistanceZ(IndexMatchingX);

                ArDistanceSq=UsefulDeltaX.*UsefulDeltaX+UsefulDeltaY.*UsefulDeltaY; %Calculate the distance square
                IndexWithinDistance=find(ArDistanceSq<TolXYInitSq); %Find those within the square of initial search radius; note Z range is already filtered above

                if ~isempty(IndexWithinDistance)

                    NumMatched=length(IndexWithinDistance);

                    StartRecordPos=CountInCloud+1;                        
                    CountInCloud=CountInCloud+NumMatched;

                    NumMoleculeWithMatches=NumMoleculeWithMatches+1;
                    MoleculeMatchPositionArrary(NumMoleculeWithMatches)=CountInCloud; %Mark the end of matched point of each molecule in the array

                    DistanceCloudArrayX(StartRecordPos:CountInCloud)= UsefulDeltaX(IndexWithinDistance); %Add all matched displacements to the cloud
                    DistanceCloudArrayY(StartRecordPos:CountInCloud)= UsefulDeltaY(IndexWithinDistance);
                    DistanceCloudArrayZ(StartRecordPos:CountInCloud)= UsefulDeltaZ(IndexWithinDistance);
                    
                end
            end
        end 
    end


    CloudX=DistanceCloudArrayX(1:CountInCloud);
    CloudY=DistanceCloudArrayY(1:CountInCloud);
    CloudZ=DistanceCloudArrayZ(1:CountInCloud);

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

        CurrentSel=(CurrentPassRSq<TolXYSq); %Find the displacements that are within the precise search radius

        CurrentPassXSel=CurrentPassX(CurrentSel);
        CurrentPassYSel=CurrentPassY(CurrentSel);

        Current_ShiftX=mean(CurrentPassXSel);
        Current_ShiftY=mean(CurrentPassYSel);


        CurrentCenterShiftSq = Current_ShiftX*Current_ShiftX + Current_ShiftY*Current_ShiftY;

    end
    

    %--After converged based on all in range, now shift to that rough center and search again using only the nearest neighbors

    CurrentSelPoints=find(CurrentSel); %These initial cloud points are used.
    CurrentSelPoints(end+1)=inf; %Patch at the end for scanning.

    CurrentPassZSel=CloudZ(CurrentSel);
    
    
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
            CurrentPassZSel(NumNearestMole)=CurrentPassZSel(RecordMole);
        end    				

    end

    NPC_ShiftX=mean(CurrentPassXSel(1:NumNearestMole))+Current_centerX; %Total shift between the two channels is the shift in the final refining step plus the accumulated shifts above in the above iterations
    NPC_ShiftY=mean(CurrentPassYSel(1:NumNearestMole))+Current_centerY;
    NPC_ShiftZ=mean(CurrentPassZSel(1:NumNearestMole)); %Relative z shift is calculated as the cloud Z center from simple Z averaging
    