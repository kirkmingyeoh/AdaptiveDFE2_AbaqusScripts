# AdaptiveDFE2_AbaqusScripts

Adaptive Direct FE2 is a newly proposed adaptive multiscale modelling method which seeks to reduce the cost of multiscale computational homogenisation. It does so by introducing RVEs only in regions of nonlinearity while homogenised macroscale elastic properties are used in the remaining regions. Through such an approach, the framework is able to reduce the computational costs required to perform a multiscale analysis significantly. Moreover, Python scripts are used to drive the analysis automatically. No further user effort is needed once the inputs are provided, simplifying the user task significantly. Further information on the Adaptive Direct FE2 method is detailed in the paper below. 

Adaptive multiscale modelling with Direct FE2 
DOI: 

This repository contains the following files: 

(1) Python scripts - These are the Python scripts used to run the Adaptive Direct FE2 analysis, which are further detailed below. 

(2) Instructions.docx - Instructions on how to prepare the user-provided input files as well as the central control script 1) **DFE2_Adaptive_Main** for running an Adaptive Direct FE2 analysis. 

--------------
Python Scripts
--------------
The following are the Python scripts to set up and run an Adaptive Direct FE2 analysis. 
Their functions in the adaptive analysis are briefly described below. 
Scripts 1) to 6) are those referred to in Appendix A of the paper mentioned and linked above. 
Script 7) is a backend Python script that is called multiple times to write the Abaqus .inp file for the next part of the analysis.
Different versions of the same script are used for different types of adaptive analysis such as 2D/3D, different element types at the macroscale and microscale, additional DOFs etc. 
Labels for these different versions are appended the end of the script name. 

1) **DFE2_Adaptive_Main** - This is the central command script which stores user provided information (detailed in (2) Instructions.docx above), reads and stores .inp file information for the macroscale and microscale, and calls the different scripts as necessary to run the Adaptive DFE2 analysis. 

2) **DFE2_Adaptive_RVEAnalysis** - This script reads and processes the microscale data to perform RVE unit strain linear analyses. It then further processes the results to obtain the homogenised macroscale elastic properties as well as the RVE stresses (which may be processed into the nonlinearity transition criterion) to be used for nonlinearity checks, which will be passed back to Script 1).  

3) **DFE2_Adaptive_MacroDataExtract** - This script reads and processes the macroscale data to initiate the adaptive multiscale analysis. It sets up the first analysis which contains only the macroscale mesh with the homogenised macroscale elastic properties provided by Script 2). 

4) **DFE2_Adaptive_OdbDataCheck** - This script performs the nonlinearity checks for the adaptive analysis using the RVE stresses or nonlinearity transition criterion. When a new increment is solved, the script is called to perform the check. If there are sufficient macroscale elements which have turned nonlinear, the adaptive process is triggered and the current analysis is terminated. 

5) **DFE2_Adaptive_AdaptPreLoad** - This script performs the first part of the adaptive process to obtain the deformed configurations of the RVEs. The script first converts the macroscale elements to be adapted into Direct FE2 elements in the new analysis .inp file. In the analysis, the macroscale elements will be deformed along with their RVEs to obtain the deformations of the RVEs at the last elastic increment. The deformed configurations of all other Instances are directly imported from the previous main multiscale analysis.  

8) **DFE2_Adaptive_AdaptMainLoad** - This script performs the second part of the adaptive process to continue on with the main multiscale analysis. The script sets up a fresh analysis for the main multiscale analysis. It then imports the deformed configurations of all Instances and reapplies the boundary conditions to allow the main multiscale analysis to continue.

11) **DFE2_Adaptive_InpFileGen** - This script is used in the backend to generate the .inp file for each Abaqus job required for the adaptive analysis. 
