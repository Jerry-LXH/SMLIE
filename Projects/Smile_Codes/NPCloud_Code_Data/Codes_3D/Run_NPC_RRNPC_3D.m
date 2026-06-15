% Example main function calling NP-Cloud-3D (NPC-3D) and RR-NP Cloud-3D (RR-NPC-3D) with example data and plotting the results.

% Load a mat file. Uncomment a line below to load the corresponding file.

load(fullfile('..', 'SimulatedData', '3D-SimulatedSMLM_HighDensity100.mat')); SubDirName='3D-Simulated_HighDensity100';
% load(fullfile('..', 'ExperimentalData', '3D-TOM20-WithFocalLock.mat')); SubDirName='3D-Exp-TOM20-FocalLock';
% load(fullfile('..', 'ExperimentalData', '3D-TOM20-NoFocalLock.mat')); SubDirName='3D-Exp-TOM20-NoFocalLock';

% Or: paste the data into an array named "SMLM_Data"
% SMLM_Data=0; SubDirName='Custom';
% Format: 4 columns single-precision: Frame, X, Y, Z. 
% All X,Y values are given in the pixel unit. The provided test data are based on a pixel size of 160 nm.
% Z values are in the unit of nm.
% Note: Our codes assume frame number always goes up in the data. Make sure to sort your data by frame before running.

%--------------Adjusted parameter---------------------
DC_SegmentSize = 100; %Drift-correction segment size: Frames in each segment for drift correction. Default 100 works well for most data. Smaller values allow finer following of drift but potentially less pairing of molecules and hence increased uncertainties.

%--------------Typically fixed parameters-------------
MaxSearchRadius_Pix1 = 0.35; %pixel / max search radius between segments in NPC Pass 1. Default is 0.35 for typical STORM sigma of 0.0625 (10 nm). Increase this for larger drifts or localization uncertainties.
MaxSearchRadius_Pix2 = 0.28; %pixel / Max search radius between segments in RR-NPC Pass 2, which is set a bit smaller, as the image is already drift-corrected once and we are refining.
TolZnm = 110; %nm / Max tolerance in Z in nm in NPC Pass 1. Default is 110 nm for typical 3D-STORM sigma of 22 nm in Z.
TolZnm2 = 90; %nm / Max tolerance of Z in RR-NPC Pass 2, which is set a bit smaller (90 nm by default), as the image is already drift-corrected once and we are refining.
ReSampleFold=12; %Resample fold when compared to the initial count of molecules in each segment. Default is 12.
%-------------------------------------------

close all;

ArraryFrames=int32(SMLM_Data(:,1)); %1D array storing the frame numbers
ArrayX=single(SMLM_Data(:,2)); %1D array storing the X positions
ArrayY=single(SMLM_Data(:,3)); %1D array storing the Y positions
ArrayZ=single(SMLM_Data(:,4)); %1D array storing the Z positions

[Xc,Yc,Zc,Xd1,Yd1,Zd1,Xd2,Yd2,Zd2, NPC_Time, RR_NPC_TotalTime]=NPC_RRNPC_3D_CallFunction(ArrayX,ArrayY,ArrayZ,ArraryFrames,DC_SegmentSize,MaxSearchRadius_Pix1,MaxSearchRadius_Pix2,TolZnm,TolZnm2,ReSampleFold); %Perform NP-Cloud (NPC) and RR-NP Cloud (RR-NPC)

MaxFrameNum=ArraryFrames(end);
DriftFrames=single((1:MaxFrameNum)');

DriftValues1=zeros(MaxFrameNum, 3, 'single'); %Here and below: 1 is for the first pass, i.e., single-referenced NPC
DriftValues2=zeros(MaxFrameNum, 3, 'single'); %Here and below: 2 is for the second pass of RR-NPC

DriftValues1(ArraryFrames(1):MaxFrameNum,:)=[Xd1 Yd1 Zd1];
DriftValues2(ArraryFrames(1):MaxFrameNum,:)=[Xd2 Yd2 Zd2];

DriftComb1=[DriftFrames DriftValues1]; %Combining frame number and drift in each frame for writing out later
DriftComb2=[DriftFrames DriftValues2];

OutDirName=fullfile('DriftCurves',SubDirName); %Resulting drift-correction curves are txt files saved in this new folder
mkdir(OutDirName) 
OutFileNameHead1=sprintf('%s/NPC_DriftCurve_prd=%d,DriftRange=%g',OutDirName,DC_SegmentSize,MaxSearchRadius_Pix1);
OutFileNameHead2=sprintf('%s_RR%d',OutFileNameHead1,ReSampleFold);

figure();
plot(DriftFrames,DriftComb1(:,2),'g',DriftFrames,DriftComb1(:,3),'b',DriftFrames,DriftComb2(:,2),'m--',DriftFrames,DriftComb2(:,3),'r--')
xlabel('Frame number')
ylabel('Drift in X/Y')
legend('NPC X','NPC Y','RR-NPC X','RR-NPC Y');
saveas(gcf,[OutFileNameHead2 '_driftXY.png']);


figure();
plot(DriftFrames,DriftComb1(:,4),'m',DriftFrames,DriftComb2(:,4),'g--')
xlabel('Frame number')
ylabel('Drift in Z')
legend('NPC Z','RR-NPC Z');
saveas(gcf,[OutFileNameHead2 '_driftZ.png']);


dlmwrite([OutFileNameHead1,'.txt'], DriftComb1, 'delimiter', '\t', 'precision', 7); %Text files with calculated drift curves
dlmwrite([OutFileNameHead2,'.txt'], DriftComb2, 'delimiter', '\t', 'precision', 7);


