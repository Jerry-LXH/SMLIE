Codes and Data for nearest paired cloud (NP-Cloud)

"Fast and robust drift correction for single-molecule localization microscopy"

-------------------------------------------------------

Code and data are provided in MATLAB format, workable with MATLAB 2015 and higher versions. 

Installation: Extract all folders and files.

-------------------------------------------------------
For 2D drift correction: Please work with codes in the “Codes_2D” folder:

Open “Run_NPC_RRNPC.m”, and uncomment one of the “load” commands at the beginning to load example datasets. Run in MATLAB to perform both NP-Cloud and RR-NP Cloud.

By default, the starting code loads ..\SimulatedData\SimulatedSMLM_HighDensity100.mat, which is the simulated SMLM dataset of 66,000 frames with features of 100 clusters/µm2. The typical running time is ~20 s to complete both NP-Cloud and RR-NP Cloud.

There is only one “Adjusted parameter” at the beginning of the code after the “load” section, namely the number of frames in each segment for drift correction. This is followed by a section of other “Typically fixed parameters”. See notes in comments.

Output: A new folder “DriftCurves” will be generated under the “Codes” folder, with another subfolder for each dataset, defined by “SubDirName”. This folder will contain two txt files corresponding to the NP-Cloud and RR-NP-Cloud calculated drift curves, respectively. Both files list 3 columns: Frame, Drift in X, and Drift in Y. The folder also saves a PNG plot comparing the drift curves calculated by the two methods.

Alternatively, open the variables “DriftComb1” “DriftComb2” in MATLAB to copy out the values.

To run on your own SMLM data:

Paste the data into an array named "SMLM_Data".

Format: A list of all localized molecules as 3 columns single-precision: Frame, X, Y. 
-We assume the Frame number always goes up. If not, please sort it by frame first before running the code. 
-Note also that the XY coordinates are in the unit of pixels. Our experiment data were based on a pixel size of 160 nm.  

-------------------------------------------------------
For 3D drift correction: Please work with codes in the “Codes_3D” folder:

Open “Run_NPC_RRNPC_3D.m”, and uncomment one of the “load” commands at the beginning to load example 3D datasets. The codes work similarly to the 2D codes above, except that they load 3D-SMLM data and perform NP-Cloud-3D and RR-NP Cloud-3D drift corrections. Output txt files list 4 columns: Frame, Drift in X, Drift in Y, and Drift in Z. PNG plots are saves separately for in-plane and z drifts.

To run on your own 3D-SMLM data: Paste the data into an array named "SMLM_Data". 

Format: A list of all localized molecules as 4 columns single-precision: Frame, X, Y, Z. Note that Z is in the unit of nm. The XY coordinates are still in the unit of pixels.
-We assume the Frame number always goes up. If not, please sort it by frame first before running the code. 