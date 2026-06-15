function [Xc,Yc,Xd1,Yd1,Xd2,Yd2, NPC_Time, RR_NPC_TotalTime]=NPC_RRNPC_CallFunction(ArrayX,ArrayY,ArrayFrames,DC_SegmentSize,MaxSearchRadius_Pix1,MaxSearchRadius_Pix2,ReSampleFold)

%Call function for NP-Cloud (NPC) and RR-NP Cloud (RR-NPC)
    
%---Inputs------
%ArrayX,ArrayY,ArrayFrames: 1D arrays storing the X position, Y position, and frame number of each localized molecule. All X,Y values are given in pixel unit. We assume frame number always goes up in the data.
%DC_SegmentSize: Drift-correction segment size: Frames in each segment for drift correction. Default is 70.
%MaxSearchRadius_Pix1: unit in pixel; max search radius between segments in Pass 1. Default is 0.3 (48 nm). Increase for larger drifts or localization uncertainties.
%MaxSearchRadius_Pix2: unit in pixel; max search radius between segments in Pass 2 of RR-NPC, which is set a bit smaller than TolInPixel1 (0.75*TolInPixel1 in the main code) as the image is already drift-corrected once and we are refining.
%ReSampleFold: Resample fold when compared to the initial count of localizations in each segment. Default is 12.

%---Outputs-------
%Xc,Yc: 1D arrays storing the X position, Y position of each localization in the final drift-corrected data.
%Xd1,Yd1: Drift curves from the NPC calculation
%Xd2,Yd2: Drift curves from the second pass of RR-NPC
%NPC_Time: Time spent by NPC.
%RR_NPC_TotalTime: Total time spent by RR-NPC 

%-----------------	
	tic
    ArrayFrames=int32(ArrayFrames);
    FrForSearch = ArrayFrames;
    StartFrame=FrForSearch(1);
    EndFrame=FrForSearch(end);
    FrForSearch(end+1)=intmax; %Patch the end to facilitate pointer scanning

    NumSegment=int32((EndFrame-StartFrame)/DC_SegmentSize); 

    SegmentStarts=zeros(NumSegment+1,1, 'uint32'); %Stores the starting position of each segment in the data arrays
    SegmentStarts(1)=1;
  
    CurrentFrame=int32(0);
    CurrentMole=uint32(1);
    
    ArrayXY=[ArrayX ArrayY];
    
    for CurrentSegment=1:NumSegment %Define segment boundaries in the localization list
        AfterCurrentEndFrame=StartFrame+DC_SegmentSize*CurrentSegment;
        while(CurrentFrame<AfterCurrentEndFrame)
            CurrentMole=CurrentMole+1;
            CurrentFrame=FrForSearch(CurrentMole);
        end
        SegmentStarts(CurrentSegment+1)=CurrentMole;
    end

	DriftX=zeros(NumSegment,1); %Stores the calculated drift
	DriftY=zeros(NumSegment,1);

    %----1st round of single-referenced NPC
    
    ArRefXY=ArrayXY(1:SegmentStarts(2)-1,:); %Reference of first round: Segment 1
    ArRefXY=sortrows(ArRefXY,2); %Sort by y for later use of "NP_Cloud_NoSort". 
   
       
    for CurrentSegment=2:NumSegment
        if ~mod(CurrentSegment,200)
            fprintf(1,'Done 1st: %d of %d segments\n', CurrentSegment, NumSegment);
        end
        
        ArCmpXY = ArrayXY(SegmentStarts(CurrentSegment):(SegmentStarts(CurrentSegment+1)-1),:); %All the XY positions in the current segment 
        ArCmpXY = sortrows(ArCmpXY,2); %Sort by y for later use of "NP_Cloud_NoSort"
        ArCmpXY = [ArCmpXY(:,1)-DriftX(CurrentSegment-1)  ArCmpXY(:,2)-DriftY(CurrentSegment-1)]; %First shift by the drift of the previous segment to bring close to the reference segment

        [NPC_ShiftX,NPC_ShiftY] = NP_Cloud_NoSort(ArRefXY, ArCmpXY, MaxSearchRadius_Pix1); %Call NP_Cloud_NoSort to find out the relative shift between the pre-shifted positions in the current segment vs. the reference
        
        DriftX(CurrentSegment)=NPC_ShiftX+DriftX(CurrentSegment-1); %Goes back to the original shift related to the reference segment
        DriftY(CurrentSegment)=NPC_ShiftY+DriftY(CurrentSegment-1);
    end
       
    SegmentCenterFrames = ((double(StartFrame) + double(DC_SegmentSize-1)/2): double(DC_SegmentSize): double(EndFrame)).';
    fi = double((StartFrame:EndFrame).');
    Xd1 = interp1(SegmentCenterFrames,DriftX,fi,'linear','extrap'); %Linear interpolation of the drift in each frame based on the calculated drift of each segment.
    Yd1 = interp1(SegmentCenterFrames,DriftY,fi,'linear','extrap'); 

    NPC_Time = toc
    
    fiSeg = floor((fi-double(StartFrame))/DC_SegmentSize)+1;
    DriftX1=[DriftX; DriftX(end)];
    DriftY1=[DriftY; DriftY(end)];
    Xd1NoExtr = DriftX1(fiSeg); 
    Yd1NoExtr = DriftY1(fiSeg); 

    RelativeFrames=ArrayFrames-StartFrame+1;
    CorrectXY1NoExtr=[Xd1NoExtr(RelativeFrames) Yd1NoExtr(RelativeFrames)];
    
    ArrayXc1Yc1NoExtr=ArrayXY-CorrectXY1NoExtr; %These are without interpolation/extrapolation, so all XY values in each segment are drift-corrected by the same amount
    
    ArrayXc1Yc1Extr=[ArrayX-Xd1(RelativeFrames) ArrayY-Yd1(RelativeFrames)]; 

    
    %----2nd round: resample localizations from the drift-corrected localization list, so that the reference image samples the entire movie, and has a size that is ReSampleFoldx of the 1st round
%     tic

    SampleEvery=round(NumSegment/ReSampleFold); %In the previous 1st step, data was divided according to frames by NumSegment, so the localization count in the reference data is the counts in the first segment, =TotalLocalizationCount/NumSegment. For RR-NPC here, the whole dataset is sampled for every (NumSegment/ReSampleFold) localizations, so the reference samples ReSampleFoldx more counts.
   
    ArrayXYSampledRef=ArrayXc1Yc1Extr(1:SampleEvery:end, :); %These localizations of the image (after 1st drift correction, interpolated) are sampled as the reference for 2nd round
    ArrayXYSampledRef=sortrows(ArrayXYSampledRef,2); %Sort by y for later use of "NP_Cloud_NoSort"

	DriftX2=zeros(NumSegment,1);
	DriftY2=zeros(NumSegment,1);    
    
    for CurrentSegment=2:NumSegment
        if ~mod(CurrentSegment,100)
            fprintf(1,'Done 2nd: %d of %d segments\n', CurrentSegment, NumSegment);
        end
        
        ArCmpXY = ArrayXc1Yc1NoExtr(SegmentStarts(CurrentSegment):(SegmentStarts(CurrentSegment+1)-1),:); %Select current segment (after 1st drift correction. not interpolated). All segments should be relatively close to no drift, so no need for pre-shift here.
        ArCmpXY = sortrows(ArCmpXY,2); %Sort by y for "NP_Cloud_NoSort"

        [NPC_ShiftX,NPC_ShiftY] = NP_Cloud_NoSort(ArrayXYSampledRef, ArCmpXY, MaxSearchRadius_Pix2); %Call NP_Cloud_NoSort to find out the relative shift between current segment (after 1st round drift correction) and the new reference.
        DriftX2(CurrentSegment)=NPC_ShiftX; %New adjustments with reference to the sampled global image
        DriftY2(CurrentSegment)=NPC_ShiftY;
    end    
        
    DriftX=DriftX+DriftX2; %The total drift is the drift from the 1st pass + New adjustments from the 2nd pass
    DriftY=DriftY+DriftY2;

    Xd2 = interp1(SegmentCenterFrames,DriftX,fi,'linear','extrap'); %Linear interpolation of the drift in each frame based on the calculated drift of each segment.
    Yd2 = interp1(SegmentCenterFrames,DriftY,fi,'linear','extrap'); 
    
    RR_NPC_TotalTime = toc
   
    Xc=ArrayX-Xd2(RelativeFrames); %Correct each localization by Fr to give the drift-corrected localizations
    Yc=ArrayY-Yd2(RelativeFrames); 