# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 15:26:25 2023

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

**Note
- This script prepares the .inp file data for the adaptive main load continuation step
- It will import all the deformed state from the pre-load step 
- All other parts of the .inp file have the same informtaion as the original base .inp file
- The script will reset the stored .inp data to the original, except for the Parts and Insts

"""

start_time = time.time()

### Updating job name
JobName = ModelName+'_R'+str(N_Model)+'b'
Heading[1] = '** Job name: '+JobName+' Model name: '+JobName

### Reset all the .inp data
Macro_Parts = []
Macro_Parts.append('** PARTS')
Macro_Parts.append('**')
# RVE_Parts left unchanged
# Insts will be updated later to add in the newly adapted elements and RVE
# RVE_Sets left unchanged
Macro_Sets = copy.deepcopy(Macro_Sets_Full)
Macro_Sets_Full = []
Ties = copy.deepcopy(Ties_Full)
Ties_Full = []
NodeTieGroups = copy.deepcopy(NodeTieGroups_Full)
NodeTieGroups_Full = []
# Eqns left unchanged
# Mats left unchanged 
Steps = copy.deepcopy(StepsMain)
if N_load < 0.9:
    Start_inc = 0.1/(1-N_load)
else:
    Start_inc = 1.0
Steps[6] = '%s, 1., 1e-05, 1.' %(str(Start_inc))

# Modify the load increment for the rigid body
for i in range(7,len(Steps)):
    if Steps[i].count('Name'):
        Start = i
        break
    
for i in range(Start,len(Steps)):
    if Steps[i].count('Name') != 0:
        RigidMark = 'NA'
        LName = Steps[i][9:len(Steps[i])-28]    
        for j in range(len(RigidLoad)):
            if LName == RigidLoad[j][0]:
                RigidMark = j
                break
    if RigidMark == 'NA':
        continue
    if Steps[i].count('*') == 0:
        Line = Steps[i].split(',')
        if (Line[-2].strip() == '1') or (Line[-2].strip() == '2'):
            DOF = int(Line[-2].strip())
            Line[-1] = str(float(Line[-1])-Disp_Rigid[RigidMark][DOF-1])
        Line2 = Line[0]
        for j in range(1,len(Line)):
            Line2 = Line2 + ',' + str(Line[j])
        Steps[i] = Line2

for i in range(len(StepEnd)):
    Steps.append(StepEnd[i])

### Import deformed state of all instances from the pre-load stage odb
Insts = []
Insts.append('**')
Insts.append('** ASSEMBLY')
Insts.append('**')
Insts.append('*Assembly, name=Assembly')
Insts.append('**')
for i in range(len(Ele_Macro)):
    Insts.append('*Instance, library='+ModelName+'_R'+str(N_Model)+'a, instance=Macro'+str(Ele_Macro[i]))
    Insts.append('**')
    Insts.append('** PREDEFINED FIELD')
    Insts.append('**')
    Insts.append('** Name: Predefined Field-1   Type: Initial State')
    Insts.append('*Import, state=yes, update=no')
    Insts.append('*End Instance')
    Insts.append('**')
for i in range(len(Ele_DFE2)):
    Insts.append('*Instance, library='+ModelName+'_R'+str(N_Model)+'a, instance=Macro'+str(Ele_DFE2[i]))
    Insts.append('**')
    Insts.append('** PREDEFINED FIELD')
    Insts.append('**')
    Insts.append('** Name: Predefined Field-1   Type: Initial State')
    Insts.append('*Import, state=yes, update=no')
    Insts.append('*End Instance')
    Insts.append('**')
    for j in range(4):
        Insts.append('*Instance, library='+ModelName+'_R'+str(N_Model)+'a, instance=Ele'+str(Ele_DFE2[i])+'-RVE'+str(j+1))
        Insts.append('**')
        Insts.append('** PREDEFINED FIELD')
        Insts.append('**')
        Insts.append('** Name: Predefined Field-1   Type: Initial State')
        Insts.append('*Import, state=yes, update=no')
        Insts.append('*End Instance')
        Insts.append('**')
for i in range(len(Ele_NewDFE2)):
    Insts.append('*Instance, library='+ModelName+'_R'+str(N_Model)+'a, instance=Macro'+str(Ele_NewDFE2[i]))
    Insts.append('**')
    Insts.append('** PREDEFINED FIELD')
    Insts.append('**')
    Insts.append('** Name: Predefined Field-1   Type: Initial State')
    Insts.append('*Import, state=yes, update=no')
    Insts.append('*End Instance')
    Insts.append('**')
    for j in range(4):
        Insts.append('*Instance, library='+ModelName+'_R'+str(N_Model)+'a, instance=Ele'+str(Ele_NewDFE2[i])+'-RVE'+str(j+1))
        Insts.append('**')
        Insts.append('** PREDEFINED FIELD')
        Insts.append('**')
        Insts.append('** Name: Predefined Field-1   Type: Initial State')
        Insts.append('*Import, state=yes, update=no')
        Insts.append('*End Instance')
        Insts.append('**')
    # Transfer the newly adapted element into the previously adapted list
    Ele_DFE2.append(Ele_NewDFE2[i])
Ele_NewDFE2 = []

### Extracting other information from the macroscale input file
for i in range(len(MacroInp)):
    if (MacroInp[i].count('*Part') != 0) and (MacroInp[i] != '*Part, name='+str(MacroPart)):
        for j in range(i,len(MacroInp)):
            Macro_Parts.append(MacroInp[j])
            if MacroInp[j] == '*End Part':
                Macro_Parts.append('**')
                break        
    
    # Other instances
    if (MacroInp[i].count('*Instance') != 0) and (MacroInp[i] != '*Instance, name='+str(MacroPart)+'-1, part='+str(MacroPart)):
        
        # Look for name of rigid body instance
        Start = MacroInp[i].index('=')
        End = MacroInp[i][Start:].index(',')
        Name = MacroInp[i][Start+1:Start+End]
        
        # Match to RP set label using RigidLoad grouping
        for j in range(len(RigidLoad)):
            if RigidLoad[j][2] == Name:
                RigidMark1 = j
                break            
        
        # Editing the new start location of the rigid body
        for j in range(i,len(MacroInp)):
            line = MacroInp[j]
            if j == i+1:
                line = line.split(',')
                for k in range(2):
                    line[k] = str(float(line[k]) + Disp_Rigid[RigidMark1][k])
                del RigidMark1
                line2 = line[0]
                for k in range(1,len(line)):
                    line2 = line2 + ',' + line[k]
                line = line2
            Insts.append(line)
            if MacroInp[j] == '*End Instance':
                Insts.append('**')
                break
            
    # Assembly reference points
    N_Assem = MacroInp.index('** ASSEMBLY')
    if (MacroInp[i] == '*Node') and (i > N_Assem):
        
        # Match to RP set label using RigidLoad grouping
        for j in range(len(RigidLoad)):
            if RigidLoad[j][2] == Name:
                RigidMark2 = j
                break
        
        # Editing the new start location of the RP
        for j in range(i,len(MacroInp)):
            line = MacroInp[j]
            if j == i+1:
                line = line.split(',')
                for k in range(1,3):
                    line[k] = str(float(line[k]) + Disp_Rigid[RigidMark2][k-1])
                del RigidMark2
                line2 = line[0]
                for k in range(1,len(line)):
                    line2 = line2 + ',' + line[k]
                line = line2
            Insts.append(line)
            if (MacroInp[j+1].count('*') != 0): 
                break

print('Adaptive main load step .inp data preparation completed')
Time1 = time.time() - start_time
Time2 = time.time() - start_time0
print('Script time: %ss'%str(Time1))
print('Total time: %ss'%str(Time2))

print>>Time_log,'Adaptive main load step .inp data preparation completed'
print>>Time_log,'Script time: %ss'%str(Time1)
print>>Time_log,'Total time: %ss'%str(Time2)











