% Example main function calling NP-Cloud (NPC) and RR-NP Cloud (RR-NPC) with example data and plotting the results.

% Load a mat file. Uncomment a line below to load the corresponding file.

% load(fullfile('..', 'SimulatedData', 'SimulatedSMLM_LowDensity10.mat')); SubDirName='Simulated_LowDensity10';
load(fullfile('..', 'SimulatedData', 'SimulatedSMLM_HighDensity100.mat')); SubDirName='Simulated_HighDensity100';
% load(fullfile('..', 'ExperimentalData', 'TMOD.mat')); SubDirName='Exp-TMOD';
% load(fullfile('..', 'ExperimentalData', 'Spectrin_CTerm.mat')); SubDirName='Exp-Spectrin_CTerm';
% load(fullfile('..', 'ExperimentalData', 'Actin.mat')); SubDirName='Exp-Actin';
% load(fullfile('..', 'ExperimentalData', 'Microtubules.mat')); SubDirName='Exp-Microtubules';
% load(fullfile('..', 'ExperimentalData', 'ER.mat')); SubDirName='Exp-ER';
% load(fullfile('..', 'ExperimentalData', '3D-TOM20-WithFocalLock.mat')); SubDirName='Exp-TOM20-FocalLock-2D';

% Or: paste the data into an array named "SMLM_Data"
% SMLM_Data=0; SubDirName='Custom';
% Format: 3 columns single-precision: Frame, X, Y. All X,Y values are given in pixel unit. The provided test data are all based on a pixel size of 160 nm.
% Note: Our codes assume frame number always goes up in the data. Make sure to sort your data by frame before running.

%--------------Adjusted parameter---------------------
DC_SegmentSize = 70; %Drift-correction segment size: Frames in each segment for drift correction. Default 70 works well for most data. Smaller values allow finer following of drift but potentially less pairing of molecules and hence increased uncertainties.

%--------------Typically fixed parameters-------------
MaxSearchRadius_Pix1 = 0.3; %pixel / max search radius between segments in NPC Pass 1. Default is 0.3 (48 nm). Increase this for larger drifts or localization uncertainties.
MaxSearchRadius_Pix2 = MaxSearchRadius_Pix1 *0.75; %pixel / Max search radius between segments in Pass 2, which is set a bit smaller, as the image is already drift-corrected once and we are refining.
ReSampleFold=12; %Resample fold when compared to the initial count of molecules in each segment. Default is 12.
%-------------------------------------------

close all;

ArraryFrames=int32(SMLM_Data(:,1)); %1D array storing the frame numbers
ArrayX=single(SMLM_Data(:,2)); %1D array storing the X positions
ArrayY=single(SMLM_Data(:,3)); %1D array storing the Y positions

[Xc,Yc,Xd1,Yd1,Xd2,Yd2, NPC_Time, RR_NPC_TotalTime]=NPC_RRNPC_CallFunction(ArrayX,ArrayY,ArraryFrames,DC_SegmentSize,MaxSearchRadius_Pix1,MaxSearchRadius_Pix2,ReSampleFold); %Perform NP-Cloud (NPC) and RR-NP Cloud (RR-NPC)

MaxFrameNum=ArraryFrames(end);
DriftFrames=single((1:MaxFrameNum)');

DriftValues1=zeros(MaxFrameNum, 2, 'single'); %Here and below: 1 is for the first pass, i.e., single-referenced NPC
DriftValues2=zeros(MaxFrameNum, 2, 'single'); %Here and below: 2 is for the second pass of RR-NPC

DriftValues1(ArraryFrames(1):MaxFrameNum,1:2)=[Xd1 Yd1];
DriftValues2(ArraryFrames(1):MaxFrameNum,1:2)=[Xd2 Yd2];

DriftComb1=[DriftFrames DriftValues1]; %Combining frame number and drift in each frame for writing out later
DriftComb2=[DriftFrames DriftValues2];

plot(DriftFrames,DriftComb1(:,2),'g',DriftFrames,DriftComb1(:,3),'b',DriftFrames,DriftComb2(:,2),'m--',DriftFrames,DriftComb2(:,3),'r--')
xlabel('Frame number')
ylabel('Drift in X/Y')

legend('NPC X','NPC Y','RR-NPC X','RR-NPC Y');

OutDirName=fullfile('DriftCurves',SubDirName); %Resulting drift-correction curves are txt files are saved in this new folder
mkdir(OutDirName) 
OutFileNameHead1=sprintf('%s/NPC_DriftCurve_prd=%d,DriftRange=%g',OutDirName,DC_SegmentSize,MaxSearchRadius_Pix1);
OutFileNameHead2=sprintf('%s_RR%d',OutFileNameHead1,ReSampleFold);

saveas(gcf,[OutFileNameHead2 '_drift.png']);
dlmwrite([OutFileNameHead1,'.txt'], DriftComb1, 'delimiter', '\t', 'precision', 7); %Text files with calculated drift curves
dlmwrite([OutFileNameHead2,'.txt'], DriftComb2, 'delimiter', '\t', 'precision', 7);


