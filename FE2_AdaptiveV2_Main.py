# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 09:51:50 2023

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

**Note
- This is the central command script to drive the adaptive Direct FE2
- Sub-processes which are not general (changes with dimensionality or element type etc) are executed with other Python scripts
- Plug and play approach allows the central command script to be flexible and stay focused on process flow
- ***Run this on Abaqus CAE nogui to minimise the runtime 

**Instructions
Macroscale
- Create and mesh the macroscale part without assigning any material properties
- Define any sets or surfaces at the part level
- Add an instance of the macroscale into assembly
- Create the Step and impose the boundary conditions on to the macroscale using Sets
- Include any necessary macroscale interactions
- Create a macroscale .inp file

Microscale
- Create and mesh the RVE, assign the material properties using Sets and give it the unscaled thickness (1 by default)
- Define any sets or surfaces at the part level
- Add an instance of the RVE into assembly
- Include any necessary internal interactions for the RVE
- Create a microscale .inp file

- Update the parameters below

"""

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from abaqusConstants import *
from odbAccess import *

import os
import time
import shutil
import math
import numpy as np
import copy

start_time0 = time.time()

### Update the following parameters
# Adaptive process parameters
BaseDir = 'E:\\Kirk Ming Abaqus\\AdaptiveDFE2-2' # Main working directory
ModelName = 'Ex6_NLGeom' # This will also be used as the base to name subsequent files
NL_PreCrit = 0.9 # NL factor for pre-emptive adaptivity
Frac_Ele = 0.4 # 0.4 % of elements to turn NL before activating an adaptive step
Frac_Time = 0.3 # 0.3 % of time step that has NL elements before activating an adaptive step
N_DataPoints = 25 # Number of data points for each macroscale strain when determining the nonlinearity criterion
NLGeom = 'YES'

# Macroscale parameters
MacroInpName = 'Ex6_Macro_NLGeom_BC1'
MacroPart = 'Macro'
Macro_Set = ['Load','Fixed'] #'Bottom','Contact'
Macro_Surface = []
Macro_Contact_Set = [] #'Contact'
#Macro_Load = ['BC-1','BC-2','BC-3'] 

# Microscale parameters
RVEInpName = 'Ex2_Micro_Ref' #'Ex2_Micro2_Ref_Mat2'
RVEPart = 'RVE'
RVE_Mat_Set = ['All'] # For assigning RVE materials
RVE_Set = []
RVE_Surface = ['Pore']
RVE_Mat = ['HDPE']
RVE_Yield = [29.4] # If assume elastic, put 'NA'

# Python script modules
RVEAnalysis = "FE2_AdaptiveV2_RVEAnalysisV2_Solid_2D.py" # See below
"""
Read and process RVE data
Setting up the .inp files for unit strain analysis
Read the .odb files and extract the corresponding data 
"""
MacroDataExtract = "FE2_AdaptiveV2_MacroDataExtract_Solid_2D.py" # Read and process macroscale part data, then create the new .inp file data
InpFileGen = "FE2_AdaptiveV2_InpFileGen_Solid_2D.py" # Generates the next .inp file based on inp array (not sure if need further refinement)
OdbDataCheck = "FE2_AdaptiveV2_OdbDataCheckV2_Solid_2D.py" # Checks the running ODB to see if nonlinearity criterion is met
AdaptPreLoad = "FE2_AdaptiveV2_AdaptPreLoadV4_Solid_2D.py" # Writes the new .inp file after adaptivity to solve for the pre-load state with RVE 
AdaptMainLoad = "FE2_AdaptiveV2_AdaptMainLoadV2_Solid_2D.py" # Writes the new .inp file after adaptivity to continue solving

'''
Beginning of adaptive code, no further user input required
'''

### Global parameters
N_Model = 1 # Version model of the number, R2 is after the first adaptive process
N_load = 0.0 # Load increment status, ends at 1.0
prev_N_load = 0.0 # Previous load increment, for checking NL

os.chdir(BaseDir)

Time_log = open(ModelName+'_TimeLog.dat','w')

### Define functions
# Search inp array for a particular part's elements and nodes; currently assumes no generate in Abaqus .inp file provided
def Search1(inp,key1,key2,out1,type1): # If dealing with nodes and coordinates, set type1 to 1 
    a = inp.index(key1)
    b = inp[a:].index(key2)
    
    for temp_line in inp[a+b+1:]:
        if (temp_line == '') or (temp_line.count("*") != 0):
            break
        temp_line = (temp_line.replace(',',' ')).split()
        temp_line.pop(0) # Removes element number in connectivity and node number in nodal coordinates
        for i in range(len(temp_line)):
            if type1 == 1:
                temp_line[i] = float(temp_line[i])
            else:
                temp_line[i] = int(temp_line[i])-1 # All node/element labels are -1 for use in Python as indices
        out1.append(temp_line)   

# Search inp array for a particular part's sets         
def Search2(inp,key1,key2,out2,type1): # If extracted output are not node/element labels, set type1 to 1
    a = inp.index(key1)
    b = inp[a:].index(key2)
    
    for temp_line in inp[a+b+1:]:
        if (temp_line == '') or (temp_line.count("*") != 0):
            break
        temp_line = ((temp_line.strip('[]')).replace(',',' ')).split()
        for i in temp_line:
            if type1 == 1:
                out2.append(float(i))
            else:
                out2.append(int(i)-1) # All node/element labels are -1 for use in Python as indices

### Extracting information from old input files
# Macroscale input file
MacroInp = []
f1 = open('%s.inp'%(MacroInpName),'r+')

while 1:
    line = f1.readline()
    if not line:
        break
    line = line.strip() # Removes additional white spaces on left and right
    MacroInp.append(line)
    
f1.close()

for i in reversed(range(len(MacroInp))): # Removing 'generate' for easier processing
    if (MacroInp[i].count('generate')!=0):
        Temp = (MacroInp[i+1].replace(',',' ')).split()
        k = 0 # Term counter
        m = 0 # Extra line counter
        for j in range(int(Temp[0]),int(Temp[1])+1,int(Temp[2])):
            if k==0:
                Temp2 = str(j)
                k = k+1
            elif k==16:
                MacroInp.insert(i+2+m,Temp2)
                m = m+1
                Temp2 = str(j)
                k = 1
            else:
                Temp2 = Temp2+', '+str(j)
                k = k+1
        MacroInp.insert(i+2+m,Temp2)
        MacroInp[i] = MacroInp[i][0:len(MacroInp[i])-10]
        del MacroInp[i+1]

# RVE input file
RVEInp = []
f1 = open('%s.inp'%(RVEInpName),'r+')

while 1:
    line = f1.readline()
    if not line:
        break
    line = line.strip() # Removes additional white spaces on left and right
    RVEInp.append(line)
    
f1.close()

for i in reversed(range(len(RVEInp))):
    if (RVEInp[i].count('generate')!=0):
        Temp = (RVEInp[i+1].replace(',',' ')).split()
        k = 0 # term counter
        m = 0 # extra line counter
        for j in range(int(Temp[0]),int(Temp[1])+1,int(Temp[2])):
            if k==0:
                Temp2 = str(j)
                k = k+1
            elif k==16:
                RVEInp.insert(i+2+m,Temp2)
                m = m+1
                Temp2 = str(j)
                k = 1
            else:
                Temp2 = Temp2+', '+str(j)
                k = k+1
        RVEInp.insert(i+2+m,Temp2)
        RVEInp[i] = RVEInp[i][0:len(RVEInp[i])-10]
        del RVEInp[i+1]
        
print('Original .inp data extraction completed')
Time = time.time() - start_time0
print('Script time: %ss'%str(Time))
print('Total time: %ss'%str(Time))

print>>Time_log,'Original .inp data extraction completed'
print>>Time_log,'Script time: %ss'%str(Time)
print>>Time_log,'Total time: %ss'%str(Time)

### Setting up RVE unit load analysis
RVEData = [[] for x in range(len(RVE_Mat_Set))]
# Store the RVE unit load output data, be it stress or strain or other parameters
# Output will be RVEData[Material Set][Element][Integration point][Unit Load][Stress component]
MacroProps = [] 
# Store macro properties from RVE analysis to be used later when setting up unadapted macro elements
execfile(RVEAnalysis)

### Setting up the macroscale elements and new .inp file
JobName = ModelName+'_R'+str(N_Model)
Heading = []
Macro_Parts = []
RVE_Parts = []
Insts = []
RVE_Sets = []
Macro_Sets = []
Surfs = []
Ties = []
Eqns = []
Mats = []
IntProps = []
Ints = []
RVEPreInts = []
Steps = []
StepsPre = []
StepsMain = []
StepEnd = []
execfile(MacroDataExtract) # Read and process macroscale part data, then create the new .inp file data
""" 
Extract nodal connect and nodal coord
Then use that information to generate the individual parts with nodes and connectivities
Put the information into the inp array, which will be updated as we go along
"""
execfile(InpFileGen) # Generates the next .inp file based on inp data arrays
os.system("abaqus job="+str(JobName))

### Start the solution and adaptivity loop
Ele_Macro = [(x+1) for x in range(len(MacroNodalConnectOld))] # MacroNodalConnectOld array should come from MacroDataExtract python
Ele_DFE2 = []
Ele_NewDFE2 = []
# RVE Parts info set separately such that the Parts can be reused by other Instances (not deleted after adapting)
RVE_SF = []
while N_load<1:
    """
    Check data
    Adapt and create new input file
    Submit job and restart loop
    
    """
    prev_inc = 0 # Job increment and N_load are different as they would mismatch once we adapt and restart 
    Res_inc = 0 # Increment for the restart file
    N_load_NL = 1 # Used to track the starting N_load where NL appears
    Criteria = 0 # Marker to indicate which adaptive process criteria was satisfied
    Adapt = 0 # Marker to activate adaptive process
    
    # Checks the .sta file for a new converged increment before checking the .odb for nonlinearity
    while 1: 
        NLInc1 = 0 # For use in case NL occurs at first increment
        if (not(os.path.isfile(str(JobName)+'.sta'))):
            continue

        sta_lines = []        
        f1 = open(JobName+'.sta','r')
        while 1:
            line = f1.readline()
            if not line:
                break
            line = line.strip() # Removes additional white spaces on left and right
            sta_lines.append(line)
        f1.close()
        
        if len(sta_lines)<5: # Job just started
            continue
        elif ((len(sta_lines[-1].split())==9)and((sta_lines[-1].split()[2]).count('U')!=0)): # Increment did not converge
            continue
        elif ((len(sta_lines[-1].split())==9)and((sta_lines[-1].split()[1])==prev_inc)): # Next increment has not been solved
            continue
        elif (sta_lines[-1]=='THE ANALYSIS HAS NOT BEEN COMPLETED'):
            break
        else: # Job finished or next valid increment solved
            current_inc = prev_inc + 1 # Primary used to check if NL occurs at increment 1 and force ODBDataCheck to check increment by increment
            execfile(OdbDataCheck) # Checks the running ODB to see if nonlinearity criterion is met
            
            ### What the analysis should do going forward
            if Adapt == 1: # Need to do adaptivity
                prev_N_load = N_load_NL
                N_load = N_load_NL # This resets the N_load to where it was before any NL set in, for use in AdaptPreLoad and AdaptMainLoad
                # Also tells the NL checking loop to not stop above at (N_load < 1) criteria even when adapt kicked in at N_load = 1
                break
            elif NLInc1 == 1: # To restart the job with a smaller starting increment due to NL at first increment
                continue
            elif N_load == 1: # No need to do adaptivity and job completed
                break
            else: # No need to do adaptivity but job not yet complete, continue checking loop
                prev_inc = current_inc
                continue
            
    if N_load == 1 and Adapt == 0:
        print('The analysis has completed')
        print>>Time_log,'The analysis has completed'
        break
    elif (sta_lines[-1]=='THE ANALYSIS HAS NOT BEEN COMPLETED'):
        print('The adaptive job did not complete due to non-convergence of '+JobName)
        print>>Time_log,'The adaptive job did not complete due to non-convergence of '+JobName
        break
    
    os.system("abaqus terminate job="+str(JobName))
    
    execfile(AdaptPreLoad) # Writes the new .inp file after adaptivity to solve for the pre-load state with RVE
    # N_model + a
    # JobName should contain the updated N_model
    execfile(InpFileGen)
    os.system("abaqus job="+str(JobName))
    
    # Check for pre-load step completion
    while 1:
        if (os.path.isfile(str(JobName)+'.lck')):
            continue
        
        elif (not(os.path.isfile(str(JobName)+'.sta'))):
            continue
        
        else:
            break
        #*** Check sta file, update pre_load to 1 when it is fully solved
        # When fully solved (last line has a step time of 1), break the while loop
    execfile(AdaptMainLoad) # Writes the new .inp file after adaptivity to continue solving
    # N_model + 1b
    # JobName should contain the updated N_model
    execfile(InpFileGen)
    os.system("abaqus job="+str(JobName))
    
Time_log.close() 
















































