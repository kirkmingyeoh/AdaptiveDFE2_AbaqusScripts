# AdaptiveDFE2_AbaqusScripts

This is a collection of Python scripts to run an Adaptive Direct FE2 analysis on Abaqus. Summary of Adaptive DFE2. Mention about paper and link/DOI

It contains versions of the following files. 

(1) Instructions.pdf - 

(2) Python scripts - 

(3) Example.zip - 

--------------
Python Scripts
--------------
The following are the Python scripts to set up and run an Adaptive Direct FE2 analysis. 
Their functions in the adaptive analysis are briefly described below. 
Scripts 1) to 6) are those referred to in Appendix A of the paper mentioned and linked above. 
Script 7) is a backend Python script that is called multiple times to write the Abaqus .inp file for the next part of the analysis.
Different versions of the same script are used for different types of adaptive analysis such as 2D/3D, different element types at the macroscale and microscale, additional DOFs etc. 
Indicators for these different versions are appended the end of the script name. 

1) DFE2_Adaptive_Main - This is the central command script which stores user provided information (detailed in (1) Instructions.pdf above), reads and stores .inp file information for the macroscale and microscale, and calls the different scripts as necessary to run the Adaptive DFE2 analysis. 

2) DFE2_Adaptive_RVEAnalysis - This script reads and processes the microscale data to perform RVE unit strain linear analyses. It then further processes the results to obtain the homogenised macroscale elastic properties 

3) DFE2_Adaptive_MacroDataExtract - 

4) DFE2_Adaptive_OdbDataCheck - 

5) DFE2_Adaptive_AdaptPreLoad - 

6) DFE2_Adaptive_AdaptMainLoad - 

7) DFE2_Adaptive_InpFileGen - 
