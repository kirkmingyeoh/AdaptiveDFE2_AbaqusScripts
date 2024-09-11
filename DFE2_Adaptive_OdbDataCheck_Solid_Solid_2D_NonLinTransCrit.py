# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 11:36:12 2023

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

**Note
- This script checkes each of the non-DFE2 macroscale elements to see if they have reached nonlinearity
- When called, it extracts the macroscale elements' strains, and use it to calculate the stresses at the RVEs by superposition
- If any of the RVE element integration point hits nonlinearity, it skips to the next increment
- The N_load required for it to just be under the nonlinearity criterion is marked and checked against the other elements
- The smallest N_load is taken as the starting point for the adapted model 

"""

start_time = time.time()

### Copy the odb file to avoid clashing with Abaqus, open the copied version
original = r''+BaseDir+'\\'+JobName+'.odb'
target = r''+BaseDir+'\\'+JobName+'2.odb'
shutil.copyfile(original,target)

# Check if the copied file has errors due to copying when Abaqus is writing to it
# Proceeds when it can be succesfully opened
while 1:
    try:
        ODBFile = BaseDir+'\\'+JobName+'2.odb'
        odb = openOdb(path=ODBFile)
    except:
        os.remove(BaseDir+'\\'+JobName+'2.odb')
        original = r''+BaseDir+'\\'+JobName+'.odb'
        target = r''+BaseDir+'\\'+JobName+'2.odb'
        shutil.copyfile(original,target)
    else: 
        break

print('ODB file opened')
Time1 = time.time() - start_time
Time2 = time.time() - start_time0
print('Script time: %ss'%str(Time1))
print('Total time: %ss'%str(Time2))
print>>Time_log,'ODB file opened'
print>>Time_log,'Script time: %ss'%str(Time1)
print>>Time_log,'Total time: %ss'%str(Time2)

N_inc_ODBMax = len(odb.steps['Step-1'].frames)

while current_inc <= (N_inc_ODBMax-1):
    start_time1 = time.time()
    
    # Looping to make sure that the current N_load is updated correctly
    # Seems like the .sta file is updated more slowly than the .odb file in the final increment
    sta_mark = 0
    while 1:
        sta_lines = []        
        f1 = open(JobName+'.sta','r')
        while 1:
            line = f1.readline()
            if not line:
                break
            line = line.strip() # Removes additional white spaces on left and right
            sta_lines.append(line)
        f1.close()
        
        for i in range(5,len(sta_lines)):
            if (int(sta_lines[i].split()[1]) == current_inc) and (sta_lines[i].split()[2].count('U') == 0):
                N_load_inc = float(sta_lines[i].split()[7])
                sta_mark = 1
                break
        
        if sta_mark == 1:
            break
    
    N_load_current = (1.0 - prev_N_load)*N_load_inc + prev_N_load
    N_load = N_load_current
    
    ### Extracting the macroscale strain values and check the RVE element stresses 
    Ele_NLPreCrit = [] # Element array for those to be pre-emptively adapted
    for i in reversed(range(len(Ele_Macro))): # Macroscale elements that have not been adapted to DFE2
        # We used reverse so removing the current would not affect the subsequent element indices 
        NL = 0 # Used as a marker for nonlinear macroscale element 
        # Checking is skipped for the rest of the macroscale element GP if one is nonlinear
        # Given now that we check every increment
        N_Ele = Ele_Macro[i]
        Ele = odb.rootAssembly.instances['MACRO'+str(Ele_Macro[i])].elements[0]
        
        # Which strain measure to take
        if NLGeom == 'YES':
            Strain = 'LE'
        elif NLGeom == 'NO':
            Strain = 'E'
        
        field = odb.steps['Step-1'].frames[current_inc].fieldOutputs[Strain]
        for j in range(4): # Integration points of the macroscale element being checked
            Val = field.getSubset(region=Ele).values[j]
            
            if NLGeom == 'YES':
                E11 = (math.e)**(Val.data[0]) - 1
                E22 = (math.e)**(Val.data[1]) - 1
                E12 = (math.e)**(Val.data[3]) - 1
            elif NLGeom == 'NO':            
                E11 = Val.data[0]
                E22 = Val.data[1]
                E12 = Val.data[3]
            
            #NLCrit_Coeff
            NLCrit = NLCrit_Coeff[0]*E11+NLCrit_Coeff[1]*E22+NLCrit_Coeff[2]*E12+NLCrit_Coeff[3]*(E11**2)+NLCrit_Coeff[4]*(E22**2)+NLCrit_Coeff[5]*(E12**2)+NLCrit_Coeff[6]*E11*E22+NLCrit_Coeff[7]*E11*E12+NLCrit_Coeff[8]*E22*E12+NLCrit_Coeff[9]*E11*E22*E12
            if NLCrit >= 1.0:
                if current_inc == 1:
                    NLInc = 1
                    break # Stop check for current macroscale element
                
                else:
                    NL = 1
                    if N_Ele not in Ele_NewDFE2: 
                        Ele_NewDFE2.append(Ele_Macro[i])
                        
                        # Removes the element from the pre-emptive adapt group to avoid repetition
                        if N_Ele in Ele_NLPreCrit:
                            Ele_Ind = Ele_NLPreCrit.index(N_Ele)
                            del Ele_NLPreCrit[Ele_Ind]
                        
                        del Ele_Macro[i]
                    break
            
            elif NLCrit >= (1.0*NL_PreCrit):
                if N_Ele not in Ele_NLPreCrit:
                    Ele_NLPreCrit.append(N_Ele)
                    
        if NLInc1 == 1:
            break # Stop check for all macro elements                    

    ### Calculate the appropriate N_load for the restart at previous increment
    # If this is the first time NL has appeared for this adaptive stage
    # Unless NLInc1 is satisfied
    if (len(Ele_NewDFE2)>0) and (N_load_NL == 1) and (NLInc1 == 0):
        Res_inc = current_inc-1
        # Obtain the restart increment and N_load accounting for the unconverged increments
        for i in range(5,len(sta_lines)):
            if ((int(sta_lines[i].split()[1]) == Res_inc) and (sta_lines[i].split()[2].count('U') == 0)):
                 N_load_res_inc = float(sta_lines[i].split()[7])
                 break      
        N_load_NL = (1.0 - prev_N_load)*N_load_res_inc + prev_N_load # Marks that NL first occured after this N_load for this adaptive stage
        # This N_load_NL will also be used as prev_N_load in case restart is needed 

    ### Check if adaptive criteria is satisfied and end check for this increment    
    if ((len(Ele_NewDFE2)>0) and ((len(Ele_NewDFE2)+len(Ele_NLPreCrit))>=(Frac_Ele*len(MacroNodalConnectOld)))):
        Criteria = 1 # Frac_ele criteria
        Adapt = 1
    elif ((len(Ele_NewDFE2)>0) and ((N_load-N_load_NL)>=Frac_Time)):
        Criteria = 2 # Frac_time criteria
        Adapt = 1
    elif ((len(Ele_NewDFE2)>0) and (N_load==1)):
        Criteria = 3 # N_load = 1 criteria
        Adapt = 1
        
    print('ODB data check for increment %s has been completed'%(str(current_inc)))
    if NLInc1 == 1:
        print('Job to be resubmitted due to nonlinearity at increment 1')
        Time1 = time.time() - start_time1
        Time2 = time.time() - start_time0
        print('Script time: %ss'%str(Time1))
        print('Total time: %ss'%str(Time2))
    elif Adapt == 1:
        # Transferring the pre-emptive elements into the list of element to be adapted, if the latter is not 0
        # Also removes it from the Ele_Macro list
        for i in range(len(Ele_NLPreCrit)):
            if Ele_NLPreCrit[i] not in Ele_NewDFE2: # THis is just a double check, the algorithm above should be prevented duplicates between Ele_NLPreCrit and Ele_NewDFE2
                Ele_NewDFE2.append(Ele_NLPreCrit[i])
            Ele_Ind = Ele_Macro.index(Ele_NLPreCrit[i])
            del Ele_Macro[Ele_Ind]
        Ele_NLPreCrit = []
        
        print('Adaptive criteria satisfied')
        print('Nonlinear elements')
        print(Ele_NewDFE2)
        print('%s nonlinear elements')%(str(len(Ele_NewDFE2)))
        print('N_load %s'%str(N_load)) # This prints where the analysis is right now
        print('Res_inc %s'%str(Res_inc))
        print('Criteria %s '%str(Criteria))
        Time1 = time.time() - start_time1
        Time2 = time.time() - start_time0
        print('Script time: %ss'%str(Time1))
        print('Total time: %ss'%str(Time2))
    else:
        print('Nonlinear elements')
        print(Ele_NewDFE2)
        print('%s nonlinear elements')%(str(len(Ele_NewDFE2)))
        print('N_load %s'%str(N_load)) # This prints where the analysis is right now
        Time1 = time.time() - start_time1
        Time2 = time.time() - start_time0
        print('Script time: %ss'%str(Time1))
        print('Total time: %ss'%str(Time2))
    
    print>>Time_log,'ODB data check for increment %s has been completed'%(str(current_inc))
    if NLInc1 == 1:
        print>>Time_log,'Job to be resubmitted due to nonlinearity at increment 1'
        print>>Time_log,'Script time: %ss'%str(Time1)
        print>>Time_log,'Total time: %ss'%str(Time2)
        break
    elif Adapt == 1:
        print>>Time_log,'Adaptive criteria satified'
        print>>Time_log,'Nonlinear elements'
        print>>Time_log,Ele_NewDFE2
        print>>Time_log,'%s nonlinear elements'%(str(len(Ele_NewDFE2)))
        print>>Time_log,'N_load %s'%str(N_load) # This prints where the analysis is right now
        print>>Time_log,'Res_inc %s'%str(Res_inc)
        print>>Time_log,'Criteria %s'%str(Criteria)
        print>>Time_log,'Script time: %ss'%str(Time1)
        print>>Time_log,'Total time: %ss'%str(Time2)
        current_inc = current_inc + 1
        break
    else:
        print>>Time_log,'Nonlinear elements'
        print>>Time_log,Ele_NewDFE2
        print>>Time_log,'%s nonlinear elements'%(str(len(Ele_NewDFE2)))
        print>>Time_log,'N_load %s'%str(N_load) # This prints where the analysis is right now
        print>>Time_log,'Script time: %ss'%str(Time1)
        print>>Time_log,'Total time: %ss'%str(Time2)
        current_inc = current_inc + 1

# To account for the additional +1 when the loop is broken
# This does not come into play when NLInc1 or NL is satisfied, only when the next ODB needs to be opened
current_inc = current_inc - 1

### Collects the necessary data when adaptive criteria is satisfied    
if Adapt == 1:    
    # Extracting nodal displacement data for adaptive pre-load step
    # Will be skipped if resubmission due to NLInc1 == 1 as Ele_NewDFE2 would be empty
    Disp_NewDFE2 = []
    # Disp_NewDFE2[which newly adapted element in Ele_NewDFE2 set][which node of that element][dof]
    for i in range(len(Ele_NewDFE2)):
        Ele = odb.rootAssembly.instances['MACRO'+str(Ele_NewDFE2[i])]
        field = odb.steps['Step-1'].frames[Res_inc].fieldOutputs['U']
        Ele_nodes = []
        for j in range(4):
            Val = field.getSubset(region=Ele).values[j]
            Ele_nodes.append([Val.data[0],Val.data[1]])
        Disp_NewDFE2.append(Ele_nodes)
        
    Disp_Contact = []
    for i in range(len(MacroContactNodes)):
        Temp = []
        for j in range(len(MacroContactNodes[i])):
            Ele = odb.rootAssembly.instances['MACRO'+str(NodeTieGroups[MacroContactNodes[i][j]][0][0])]
            field = odb.steps['Step-1'].frames[Res_inc].fieldOutputs['U']
            Val = field.getSubset(region=Ele).values[NodeTieGroups[MacroContactNodes[i][j]][0][1]-1]
            Temp.append([Val.data[0],Val.data[1]])
        Disp_Contact.append(Temp)
        
    if N_Model == 1:  
        Disp_Rigid = []
        for i in range(len(RigidLoad)):
            Temp = []
            RP = odb.rootAssembly.nodeSets[RigidLoad[i][1]]
            field = odb.steps['Step-1'].frames[Res_inc].fieldOutputs['U']
            Val = field.getSubset(region=RP).values[0]
            Temp = [Val.data[0],Val.data[1]]
            
            Disp_Rigid.append(Temp)
    else:
        for i in range(len(RigidLoad)):
            RP = odb.rootAssembly.nodeSets[RigidLoad[i][1]]
            field = odb.steps['Step-1'].frames[Res_inc].fieldOutputs['U']
            Val = field.getSubset(region=RP).values[0]
            Disp_Rigid[i][0] = Disp_Rigid[i][0] + Val.data[0]
            Disp_Rigid[i][1] = Disp_Rigid[i][1] + Val.data[1]

odb.close()            
os.remove(BaseDir+'\\'+JobName+'2.odb')

print('ODB file closed')
Time1 = time.time() - start_time
Time2 = time.time() - start_time0
print('Script time: %ss'%str(Time1))
print('Total time: %ss'%str(Time2))
print>>Time_log,'ODB file closed'
print>>Time_log,'Script time: %ss'%str(Time1)
print>>Time_log,'Total time: %ss'%str(Time2)

### Follow up due to nonlinearity at increment 1, by resubmitting the current job with a smaller starting increment
if NLInc1 == 1:
    # Stop the job that is currently running
    os.system("abaqus terminate job="+str(JobName))
    
    # Remove all the relevant job files when the job has stopped completely
    while 1:    
        if (not(os.path.isfile(str(JobName)+'.lck'))):
            while 1:
                try:
                    os.remove(str(JobName)+'.com')
                except:
                    continue
                else:
                    break
            while 1:
                try:
                    os.remove(str(JobName)+'.dat')
                except:
                    continue
                else:
                    break
            while 1:
                try:
                    os.remove(str(JobName)+'.inp')
                except:
                    continue
                else:
                    break
            while 1:
                try:
                    os.remove(str(JobName)+'.log')
                except:
                    continue
                else:
                    break
            while 1:
                try:
                    os.remove(str(JobName)+'.mdl')
                except:
                    continue
                else:
                    break
            while 1:
                try:
                    os.remove(str(JobName)+'.msg')
                except:
                    continue
                else:
                    break
            while 1:
                try:
                    os.remove(str(JobName)+'.odb')
                except:
                    continue
                else:
                    break
            while 1:
                try:
                    os.remove(str(JobName)+'.prt')
                except:
                    continue
                else:
                    break
            while 1:
                try:
                    os.remove(str(JobName)+'.res')
                except:
                    continue
                else:
                    break
            while 1:
                try:
                    os.remove(str(JobName)+'.sim')
                except:
                    continue
                else:
                    break
            while 1:
                try:
                    os.remove(str(JobName)+'.sta')
                except:
                    continue
                else:
                    break
            while 1:
                try:
                    os.remove(str(JobName)+'.stt')
                except:
                    continue
                else:
                    break
            break
        else: 
            continue
    
    # Set the new starting increment at half the previous, write the input file and submit the job    
    Start_inc = Start_inc*0.5
    Steps[6] = '%s, 1., 1e-05, 1.' %(str(Start_inc))
    execfile(InpFileGen)
    os.system("abaqus job="+str(JobName))
    print('Job resubmitted')
    Time1 = time.time() - start_time
    Time2 = time.time() - start_time0
    print('Script time: %ss'%str(Time1))
    print('Total time: %ss'%str(Time2))
    print>>Time_log,'Job resubmitted'
    print>>Time_log,'Script time: %ss'%str(Time1)
    print>>Time_log,'Total time: %ss'%str(Time2)











