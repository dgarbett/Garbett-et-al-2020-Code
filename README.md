# Garbett-et-al-2020-Code
Garbett, D. et al. “T-Plastin reinforces membrane protrusions to bridge matrix gaps during cell migration”
Software code guide and brief instructions

General Instructions
All code is run via Matlab and requires working knowledge of its use. The example code provided here has been run on both Mac and PC with only minor changes to the scripts required to use on your data of interest. Code has been run on multiple Matlab versions ranging between 2014b-2019a. Matlab install times vary by computer and network connection, installation at Stanford can take upwards of one hour when installing all toolkit packages. Note that some scripts call functions either built into Matlab toolkits which occasionally vary by version. For custom sub-functions we have provided them in the general functions folder or with their corresponding specialized scripts and these folders should be added to the Matlab path in order to be called correctly. 

Scratch Assay Tracking
Figures: S4E, S4G, S5B
This code uses cell trace data generated using tracking code described previously(1,2) and is also available at https://github.com/scappell/Cell_tracking
An example set of trace data is provided (tracedata_2_2_1.mat) and can be analyzed using the provided script (Scratch Assay Tracking.m). It tracks average cell velocities based on a nuclear marker and can be configured to only analyze the outer leader cells near a scratch or all cells further back from the scratch.

Ratiometric Analysis
Figures: 1H, 2F, 3A, 3F, 3G, 4A
This code generates a parula colormap of two channels as a ratio and is similar to that previously described(1). An alternate method of masking is provided if imaging variables change. Sample data is provided in the folder ‘Example Ratio and Edge Velocity RawData’.

Edge Velocity and Cross Correlation
Figures: 3B-E
This code generates a parula colormap of the outer edge of the cell mask divided into windows and has previously been described3. The main script used is called ‘Edge_Velocity_main.m’ and calls on many subfunctions which are also provided. It is also used to generate kymographic panels showing edge regions velocities (protrusion vs retraction) and ratiometric enrichment of inputted sensors. ‘XCorr.m’ is then used to generate cross correlations of sensors with protrusion events. Sample data is provided in the folder ‘Example Ratio and Edge Velocity RawData’.

Ladder Assay Analysis
Figures: 1I, 4C-F, S5E
This code is used to measure cell area changes once cells make contact with fibronectin ladder patterns. The main script used is called ‘PRIMO_tracking.m’ and required sub-functions are provided. The script ‘getAreaChange.m’ was used to extract cell area over time data to plot (as shown in Fig 4D).

ICQ Colocalization Analysis
Figure: 6D
This code is used to measure the colocalization in two-color STED images and its methodology has been previously described(4). ‘coloc.m’ is the master program that loads data and calls the other functions required. Sample images are provided in the directory for analysis. 

References
1.	Hayer, A. et al. Engulfed cadherin fingers are polarized junctional structures between collectively migrating endothelial cells. Nat. Cell Biol. 18, 1311–1323 (2016).
2.	Cappell, S. D. et al. EMI1 switches from being a substrate to an inhibitor of APC/CCDH1 to start the cell cycle. Nature 558, 313–317 (2018).
3.	Yang, H. W., Collins, S. R. & Meyer, T. Locally excitable Cdc42 signals steer cells during chemotaxis. Nat. Cell Biol. 18, 191–201 (2016).
4.	Li, Q. et al. A Syntaxin 1, Gαo, and N-Type Calcium Channel Complex at a Presynaptic Nerve Terminal: Analysis by Quantitative Immunocolocalization. J. Neurosci. 24, 4070–4081 (2004).

