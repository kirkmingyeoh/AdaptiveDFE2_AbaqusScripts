# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 15:26:25 2023

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

**Note
- This script prepares the .inp file data for the adaptive pre-load step
- For each newly adapted Direct FE2 element, RVEs will be added and macro properties changed to null
- These elements will be unTied and pre-loaded using nodal displacements 
- All other elements will have the Ties maintained, import deformed state from previous run
- At the unTied boundaries, also apply pre-load using nodal displacements 
- These data will be directly written as the next .inp file for stage a 

"""

start_time = time.time()

### Define functions
def Bilin_Interpolation(tsi,eta):
    N1=float(0.25*(1-tsi)*(1-eta))
    N2=float(0.25*(1+tsi)*(1-eta))
    N3=float(0.25*(1+tsi)*(1+eta))
    N4=float(0.25*(1-tsi)*(1+eta))
    return N1,N2,N3,N4

### Updating job name
N_Model += 1
JobName = ModelName+'_R'+str(N_Model)+'a'
Heading[1] = '** Job name: '+JobName+' Model name: '+JobName

### Establish the Instance information for macroscale elements which are not adapted at this stage    
Insts = []
Insts.append('**')
Insts.append('** ASSEMBLY')
Insts.append('**')
Insts.append('*Assembly, name=Assembly')
Insts.append('**')

if N_Model==2:
    for i in range(len(Ele_Macro)):
        Insts.append('*Instance, library='+ModelName+'_R'+str(N_Model-1)+', instance=Macro'+str(Ele_Macro[i]))
        Insts.append('**')
        Insts.append('** PREDEFINED FIELD')
        Insts.append('**')
        Insts.append('** Name: Predefined Field-1   Type: Initial State')
        Insts.append('*Import, increment='+str(Res_inc)+', state=yes, update=no')
        Insts.append('*End Instance')
        Insts.append('**')
else:
    for i in range(len(Ele_Macro)):
        Insts.append('*Instance, library='+ModelName+'_R'+str(N_Model-1)+'b, instance=Macro'+str(Ele_Macro[i]))
        Insts.append('**')
        Insts.append('** PREDEFINED FIELD')
        Insts.append('**')
        Insts.append('** Name: Predefined Field-1   Type: Initial State')
        Insts.append('*Import, increment='+str(Res_inc)+', state=yes, update=no')
        Insts.append('*End Instance')
        Insts.append('**')
    for i in range(len(Ele_DFE2)):
        Insts.append('*Instance, library='+ModelName+'_R'+str(N_Model-1)+'b, instance=Macro'+str(Ele_DFE2[i]))
        Insts.append('**')
        Insts.append('** PREDEFINED FIELD')
        Insts.append('**')
        Insts.append('** Name: Predefined Field-1   Type: Initial State')
        Insts.append('*Import, increment='+str(Res_inc)+', state=yes, update=no')
        Insts.append('*End Instance')
        Insts.append('**')
        for j in range(4):
            Insts.append('*Instance, library='+ModelName+'_R'+str(N_Model-1)+'b, instance=Ele'+str(Ele_DFE2[i])+'-RVE'+str(j+1))
            Insts.append('**')
            Insts.append('** PREDEFINED FIELD')
            Insts.append('**')
            Insts.append('** Name: Predefined Field-1   Type: Initial State')
            Insts.append('*Import, increment='+str(Res_inc)+', state=yes, update=no')
            Insts.append('*End Instance')
            Insts.append('**')

### Preserving old lists for use in next adaptive step
Macro_Sets_Full = copy.deepcopy(Macro_Sets)
Ties_Full = copy.deepcopy(Ties)
NodeTieGroups_Full = copy.deepcopy(NodeTieGroups)
MacroLoadNodes_Full = copy.deepcopy(MacroLoadNodes)
StepsPre = copy.deepcopy(StepsMain)

### Edit the loading steps for main load's pre-load value, as well as the increment size and limits
# Works for boundary conditions, not sure about other types of loads yet
StepsPre[6] = '1., 1., 1e-05, 1.'

for i in range(7,len(StepsPre)):
    if StepsPre[i].count('Name') != 0:
        Factor = N_load
        LName = StepsPre[i][9:len(StepsPre[i])-28]
#        print(StepsPre[i])
        for j in range(len(RigidLoad)):
            if LName == RigidLoad[j][0]:
                Factor = 0.0
                break    
    if StepsPre[i].count('*') == 0:
        Line = StepsPre[i].split(',')
        Line[-1] = str(float(Line[-1])*Factor)
        Line2 = Line[0]
        for j in range(1,len(Line)):
            Line2 = Line2 + ',' + str(Line[j])
        StepsPre[i] = Line2

### Transfer original loads to Steps for separate handling later in case sets are removed by PreLoad
Steps = copy.deepcopy(StepsPre)
StepsPre = []
        
### Extracting RVE Interaction, Properties and related Sets/Surfaces
# Seems to not be necessary as Interaction and Properties extraced in RVEAnalysis
# RVE sets and surfaces should have been included in the part definition as it copies everything from start to section

### Including any RVE interaction properties for the first adapt
if N_Model==2:
    for i in range(len(RVEIntProp)):
        for j in range(len(RVEIntProp[i])):
            line = RVEIntProp[i][j]
            if line.count('name') != 0:
                line = line.replace('name=','name=RVE-')            
            IntProps.append(line)

### Looping across all the elements that need to be adapted
GP = [[-3**-0.5,-3**-0.5],[3**-0.5,-3**-0.5],[3**-0.5,3**-0.5],[-3**-0.5,3**-0.5]]
Macro_Parts = []
Macro_Parts.append('** PARTS')
Macro_Parts.append('**')

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

### Setting up the new DFE2 elements, impose the PBCs and handle the NodeTieGroups
for i in range(len(Ele_NewDFE2)):
    N = Ele_NewDFE2[i]-1
    
    # Calculating element nodal connectivity based on DFE2 convention
    x1 = MacroNodalCoordOld[MacroNodalConnectOld[N][0]][0]
    x2 = MacroNodalCoordOld[MacroNodalConnectOld[N][1]][0]
    x3 = MacroNodalCoordOld[MacroNodalConnectOld[N][2]][0]
    x4 = MacroNodalCoordOld[MacroNodalConnectOld[N][3]][0]
    y1 = MacroNodalCoordOld[MacroNodalConnectOld[N][0]][1]
    y2 = MacroNodalCoordOld[MacroNodalConnectOld[N][1]][1]
    y3 = MacroNodalCoordOld[MacroNodalConnectOld[N][2]][1]
    y4 = MacroNodalCoordOld[MacroNodalConnectOld[N][3]][1]
    X0 = [x1,x2,x3,x4]
    Y0 = [y1,y2,y3,y4]
    X = []
    Y = []
    Nodes = []
    Ang = []
    
    X_ave = sum(X0)/4.0
    Y_ave = sum(Y0)/4.0
    for j in range(4):
        X1 = X0[j]-X_ave
        Y1 = Y0[j]-Y_ave
        if X1==0:
            if Y1 > 0:
                theta = math.pi*0.5
            else:
                theta = math.pi*1.5
        else:       
            theta = math.atan(Y1/X1)
        if X1<0:
            theta = theta+math.pi
        if theta<0:
            theta = theta+2*(math.pi)
        Ang.append(theta*360/(2*(math.pi)))
    SAng = sorted(Ang)
    for j in range(4):
        Order1 = [2,3,0,1]
        Ind_Nodes = Ang.index(SAng[Order1[j]])
        Nodes.append(Ind_Nodes)
        X.append(MacroNodalCoordOld[MacroNodalConnectOld[N][Ind_Nodes]][0])
        Y.append(MacroNodalCoordOld[MacroNodalConnectOld[N][Ind_Nodes]][1])
    
    # Establish the Part information
    Macro_Parts.append('*Part, name=Macro'+str(N+1))
    Macro_Parts.append('*Node')
    Macro_Parts.append('1, '+str(x1)+', '+str(y1))
    Macro_Parts.append('2, '+str(x2)+', '+str(y2))
    Macro_Parts.append('3, '+str(x3)+', '+str(y3))
    Macro_Parts.append('4, '+str(x4)+', '+str(y4))
    Macro_Parts.append('*Element, type=CPS4')
    Macro_Parts.append('1, 1, 2, 3, 4')
    Macro_Parts.append('*Elset, elset=Mat')
    Macro_Parts.append(' 1,')
    Macro_Parts.append('** Section: Solid_ElasticM'+str(N+1))
    Macro_Parts.append('*Solid Section, elset=Mat, material=Elastic_M'+str(N+1))
    Macro_Parts.append(str(Thickness)+',')
    Macro_Parts.append('*End Part')
    Macro_Parts.append('**')
    
    # Establish the Instance information
    Insts.append('*Instance, name=Macro'+str(N+1)+', part=Macro'+str(N+1))
    Insts.append('*End Instance')
    Insts.append('**')
    
    # Edit the Material info
    Ind_Mat = Mats.index('*Material, name=Elastic_M'+str(N+1))
    Mats[Ind_Mat+1] = '*Elastic'
    Mats[Ind_Mat+2] = '1e-10, 1e-10'
    del Mats[Ind_Mat+3]
    
    # Edit Ties and pre-loads
    for j in range(len(NodeTieGroups)):
        MacroLoad_Mark = 0 # Marker for whether it is a macro load node
        MacroContact_Mark = 0 # Marker for whether it is a macro contact node
        
        # Checking for node groups involved in contact
        # Handle them separately
        for k in range(len(MacroContactNodes)):
            for l in range(len(MacroContactNodes[k])):
                if j in MacroContactNodes[k]:
                    MacroContact_Mark = 1 
                    break
                    # The newly adapted macro element node is involved in a macro contact
                    # To be handled separately
            if MacroContact_Mark == 1:
                break
        if MacroContact_Mark == 1:
            continue
        
        # Only 1 node in Tie group
        if (len(NodeTieGroups[j]) == 1): 
            if NodeTieGroups[j][0][0] == N+1: # This newly adapted macro element is involved
                # Removing nodes of adapted elements from loading set 
                for k in range(len(MacroLoadNodes)):
                    if j in MacroLoadNodes[k]:
                        NInd = MacroLoadNodes[k].index(j)
                        del MacroLoadNodes[k][NInd]
                        for l in reversed(range(len(Macro_Sets))):
                            if (Macro_Sets[l] == ('*Nset, nset='+str(MacroLoadSet[k])+', instance=Macro'+str(NodeTieGroups[j][0][0]))):
                                if Macro_Sets[l+1] == str(NodeTieGroups[j][0][1]):
                                    del Macro_Sets[l+1]
                                    del Macro_Sets[l]
                                    break
                Macro_Sets.append('*Nset, nset=PreLoad'+str(j+1)+', instance=Macro'+str(N+1))
                Macro_Sets.append(str(NodeTieGroups[j][0][1]))
                StepsPre.append('** Name: PreLoad'+str(j+1)+' Type: Displacement/Rotation')
                StepsPre.append('*Boundary')
                StepsPre.append('PreLoad'+str(j+1)+',1,1,'+str(Disp_NewDFE2[i][NodeTieGroups[j][0][1]-1][0]))
                StepsPre.append('PreLoad'+str(j+1)+',2,2,'+str(Disp_NewDFE2[i][NodeTieGroups[j][0][1]-1][1]))
                # Disp_NewDFE2[which newly adapted element in Ele_NewDFE2 set][which node of that element][dof] 
                del NodeTieGroups[j][0]
        
        # More than 1 node in Tie group
        else:
            # Should not have to reverse loop here, as k loop will break after one of the process below is done
            # Break k loop as each macro element should only appear once in each NTG[j], so at most only 1 k value is relevant
            for k in range(len(NodeTieGroups[j])):
                if NodeTieGroups[j][k][0] == N+1: # This newly adapted macro element is involved
                    if len(NodeTieGroups[j]) == 2: # Only 2 nodes in tie                                                
                        for l in range(len(Ties)): # Find and remove the Ties, we ignore the Sets and Surfaces for convenience
                            if (('** Constraint: Tie-'+str(j+1)) in Ties[l]):
                                del Ties[l+2]
                                del Ties[l+1]
                                del Ties[l]
                                break  
                        # Removing nodes of adapted elements from loading set
                        for l in range(len(MacroLoadNodes)):
                            if j in MacroLoadNodes[l]:
                                NInd = MacroLoadNodes[l].index(j)
                                del MacroLoadNodes[l][NInd]
                                for m in reversed(range(len(Macro_Sets))):
                                    if (Macro_Sets[m] == ('*Nset, nset='+str(MacroLoadSet[l])+', instance=Macro'+str(NodeTieGroups[j][0][0]))):
                                        if Macro_Sets[m+1] == str(NodeTieGroups[j][0][1]):
                                            del Macro_Sets[m+1]
                                            del Macro_Sets[m]
                                            break
                        # Add the nodes to a pre-load set
                        Macro_Sets.append('*Nset, nset=PreLoad'+str(j+1)+', instance=Macro'+str(NodeTieGroups[j][0][0]))
                        Macro_Sets.append(str(NodeTieGroups[j][0][1]))
                        Macro_Sets.append('*Nset, nset=PreLoad'+str(j+1)+', instance=Macro'+str(NodeTieGroups[j][1][0]))
                        Macro_Sets.append(str(NodeTieGroups[j][1][1]))
                        # If the pre-load step does not exist, create the pre-load
                        if ('** Name: PreLoad'+str(j+1)+' Type: Displacement/Rotation') not in StepsPre:
                            StepsPre.append('** Name: PreLoad'+str(j+1)+' Type: Displacement/Rotation')
                            StepsPre.append('*Boundary')
                            StepsPre.append('PreLoad'+str(j+1)+',1,1,'+str(Disp_NewDFE2[i][NodeTieGroups[j][k][1]-1][0]))
                            StepsPre.append('PreLoad'+str(j+1)+',2,2,'+str(Disp_NewDFE2[i][NodeTieGroups[j][k][1]-1][1]))
                        del NodeTieGroups[j][1]
                        del NodeTieGroups[j][0]                        
                        # Remove set, surface and tie
                        # Check if this group has a pre-load
                        # If yes, add into pre-load set
                        # If no, call both into a load set and apply pre-load
                        # Pop the [j] from NTG, if so loop NTG in reverse???
                    elif (len(NodeTieGroups[j]) > 2) and (k == 0): # Originally master node
                        for l in range(len(Macro_Sets)):
                            # Edit the Master to be the 2nd node
                            # Remove the Slave set for the 2nd node
                            if (('*Nset, nset=TieGroupM'+str(j+1)) in Macro_Sets[l]):
                                # Change 2nd node in Ties group to Master
                                Macro_Sets[l] = '*Nset, nset=TieGroupM'+str(j+1)+', instance=Macro'+str(NodeTieGroups[j][1][0])
                                Macro_Sets[l+1] = str(NodeTieGroups[j][1][1])
                                # Remove 2nd node from original Slave Set
                                del Macro_Sets[l+3]
                                del Macro_Sets[l+2]
                                break 
                        # Removing nodes of adapted elements from loading set
                        for l in range(len(MacroLoadNodes)):
                            if j in MacroLoadNodes[l]:
                                NInd = MacroLoadNodes[l].index(j)
                                del MacroLoadNodes[l][NInd]
                                for m in reversed(range(len(Macro_Sets))):
                                    if (Macro_Sets[m] == ('*Nset, nset='+str(MacroLoadSet[l])+', instance=Macro'+str(NodeTieGroups[j][0][0]))):
                                        if Macro_Sets[m+1] == str(NodeTieGroups[j][0][1]):
                                            del Macro_Sets[m+1]
                                            del Macro_Sets[m]
                                            break
                        # Creating the pre-load set for the new Master
                        Macro_Sets.append('*Nset, nset=PreLoad'+str(j+1)+', instance=Macro'+str(NodeTieGroups[j][1][0]))
                        Macro_Sets.append(str(NodeTieGroups[j][1][1]))
                        if ('** Name: PreLoad'+str(j+1)+' Type: Displacement/Rotation') not in StepsPre:
                            # If the group had no pre-load, create the set for the old Master as well as load
                            Macro_Sets.append('*Nset, nset=PreLoad'+str(j+1)+', instance=Macro'+str(NodeTieGroups[j][0][0]))
                            Macro_Sets.append(str(NodeTieGroups[j][0][1]))
                            StepsPre.append('** Name: PreLoad'+str(j+1)+' Type: Displacement/Rotation')
                            StepsPre.append('*Boundary')
                            StepsPre.append('PreLoad'+str(j+1)+',1,1,'+str(Disp_NewDFE2[i][NodeTieGroups[j][k][1]-1][0]))
                            StepsPre.append('PreLoad'+str(j+1)+',2,2,'+str(Disp_NewDFE2[i][NodeTieGroups[j][k][1]-1][1]))
                        del NodeTieGroups[j][0]
                        # Change master and slave node sets
                        # Check if this group has a pre-load before
                        # If yes, add into the pre-load set
                        # If no, call the K node and new master into the pre-laod set, then apply pre-load
                        # Del the [0] from NTG[j]
                    elif (len(NodeTieGroups[j]) > 2) and (k != 0): # Originally slave node
                        for l in range(len(Macro_Sets)):
                            # Finding the Tie group node sets
                            if ('*Nset, nset=TieGroupM'+str(j+1)) in Macro_Sets[l]:
                                # Deleting the corresponding Slave node set
                                del Macro_Sets[l+2*k+1]
                                del Macro_Sets[l+2*k]
                                break
                        # Removing nodes of adapted elements from loading set
                        for l in range(len(MacroLoadNodes)):
                            if j in MacroLoadNodes[l]:
                                NInd = MacroLoadNodes[l].index(j)
                                del MacroLoadNodes[l][NInd]
                                for m in reversed(range(len(Macro_Sets))):
                                    if (Macro_Sets[m] == ('*Nset, nset='+str(MacroLoadSet[l])+', instance=Macro'+str(NodeTieGroups[j][0][0]))):
                                        if Macro_Sets[m+1] == str(NodeTieGroups[j][0][1]):
                                            del Macro_Sets[m+1]
                                            del Macro_Sets[m]
                                            break
                        # Create the pre-load set for this newly independent node
                        Macro_Sets.append('*Nset, nset=PreLoad'+str(j+1)+', instance=Macro'+str(NodeTieGroups[j][k][0]))
                        Macro_Sets.append(str(NodeTieGroups[j][k][1]))
                        if ('** Name: PreLoad'+str(j+1)+' Type: Displacement/Rotation') not in StepsPre:
                            # If there was no pre-load, create the set for the Master node as well as load
                            Macro_Sets.append('*Nset, nset=PreLoad'+str(j+1)+', instance=Macro'+str(NodeTieGroups[j][0][0]))
                            Macro_Sets.append(str(NodeTieGroups[j][0][1]))
                            StepsPre.append('** Name: PreLoad'+str(j+1)+' Type: Displacement/Rotation')
                            StepsPre.append('*Boundary')
                            StepsPre.append('PreLoad'+str(j+1)+',1,1,'+str(Disp_NewDFE2[i][NodeTieGroups[j][k][1]-1][0]))
                            StepsPre.append('PreLoad'+str(j+1)+',2,2,'+str(Disp_NewDFE2[i][NodeTieGroups[j][k][1]-1][1]))
                        del NodeTieGroups[j][k]
                        # Remove from slave node sets
                        # Check if this group has a pre-load before
                        # If yes, add into the pre-load set
                        # If no, call the K node and new master into the pre-laod set, then apply pre-load
                        
                    break # Break from the k loop within NTG[j], each NTG will at most involve the adapted macro element once
                    
    # Create node sets for RVE PBCs
    for j in range(4):
        RVE_Sets.append('*Nset, nset=Ele'+str(N+1)+'-N'+str(j+1)+', instance=Macro'+str(N+1))
        RVE_Sets.append(str(Nodes[j]+1)+',')
    
    # Calculating RVE information
    C = np.array([[1,-1,-1,1],[1,1,-1,-1],[1,1,1,1],[1,-1,1,-1]])
    C_inv = np.linalg.inv(C)
    [a0,a1,a2,a3] = np.dot(C_inv,[X[0],X[1],X[2],X[3]])
    [b0,b1,b2,b3] = np.dot(C_inv,[Y[0],Y[1],Y[2],Y[3]])
    
    for j in range(4):
        [tsi,eta] = GP[j]
        J = np.array([[a1+a3*eta,b1+b3*eta],[a2+a3*tsi,b2+b3*tsi]])
        J_RVE = abs(np.linalg.det(J)/(B_RVE*H_RVE))
        
        # Create a new RVE part if the current scaling factor has not been used before
        if round(J_RVE,5) in RVE_SF:
            Ind_RVE = RVE_SF.index(round(J_RVE,5))
        else:
            Ind_RVE = len(RVE_SF)
            RVE_SF.append(round(J_RVE,5))
            RVE_Parts.append('*Part, name='+str(RVEPart)+str(Ind_RVE))
            RVE_Parts.append('*Node')
            for k in range(len(RVENodalCoord)):
                RVE_Parts.append(str(k+1)+', '+str(RVENodalCoord[k][0])+','+str(RVENodalCoord[k][1]))
            Start0 = RVEInp.index('*Part, name='+str(RVEPart))
            Start1 = Start0 + RVEInp[Start0:].index('*Element, type=CPS4')
            for k in range(Start0,len(RVEInp)):
                if RVEInp[k].count('** Section:') != 0:
                    End = k
                    break
            for k in range(Start1,End):
                RVE_Parts.append(RVEInp[k])
            for k in range(len(RVE_Mat)):
                RVE_Parts.append('** Section: Solid_%s' %(str(RVE_Mat[k])))
                RVE_Parts.append('*Solid Section, elset=%s, material=%s' %(str(RVE_Mat_Set[k]),str(RVE_Mat[k])))
                RVE_Parts.append(str(Thickness*J_RVE)+',')
            RVE_Parts.append('*End Part')
            RVE_Parts.append('**')
        
        # Calculating RVE Instance position and adding the Instance
        N1,N2,N3,N4 = Bilin_Interpolation(tsi,eta)
                
        RVE_X = N1*X[0] + N2*X[1] + N3*X[2] + N4*X[3]
        RVE_Y = N1*Y[0] + N2*Y[1] + N3*Y[2] + N4*Y[3]
        
        Insts.append('*Instance, name=Ele'+str(N+1)+'-RVE'+str(j+1)+', part=RVE'+str(Ind_RVE))
        Insts.append(str(RVE_X)+', '+str(RVE_Y)+', 0.')
        Insts.append('*End Instance')
        Insts.append('**')
        
        # Creating sets and constraints for the RVE PBC
        J_inv = np.linalg.inv(J)
        dN1 = [-0.25*(1-eta),-0.25*(1-tsi)]
        dN2 = [0.25*(1-eta),-0.25*(1+tsi)]
        dN3 = [0.25*(1+eta),0.25*(1+tsi)]
        dN4 = [-0.25*(1+eta),0.25*(1-tsi)]
        N_NatDeriv = [dN1,dN2,dN3,dN4]
        N_GloDeriv = [[],[],[],[]]
        
        dx = B_RVE
        dy = H_RVE
        
        for k in range(4):
            N_GloDeriv[k] = np.dot(J_inv,np.transpose(np.array(N_NatDeriv[k]))) # Shape function gradients
            
        for k in range(len(RVE_Set_List)):
            if RVE_Set_Type[k] == 'N':
               RVE_Sets.append('*Nset, nset=Ele'+str(N+1)+'-RVE'+str(j+1)+'-'+str(RVE_Set[k])+', instance=Ele'+str(N+1)+'-RVE'+str(j+1))
               for l in range(len(RVE_Set_List[k])):
                   RVE_Sets.append(RVE_Set_List[k][l])
                   
            if RVE_Set_Type[k] == 'El':
               RVE_Sets.append('*Elset, elset=Ele'+str(N+1)+'-RVE'+str(j+1)+'-'+str(RVE_Set[k])+', instance=Ele'+str(N+1)+'-RVE'+str(j+1))
               for l in range(len(RVE_Set_List[k])):
                   RVE_Sets.append(RVE_Set_List[k][l])
                   
        for k in range(len(RVE_Surface_List)):
            for l in range(len(RVE_Surface_List[k])):
                line = RVE_Surface_List[k][l]
                if line.count('name=') != 0:
                    line = line.replace('name=','name=Ele'+str(N+1)+'-RVE'+str(j+1)+'-')
                else:
                    line = 'Ele'+str(N+1)+'-RVE'+str(j+1)+'-'+line
                Surfs.append(line)
        
        for k in range(len(ParingFacesLR)):
            RVE_Sets.append('*Nset, nset=Ele'+str(N+1)+'-RVE'+str(j+1)+'-NodeL'+str(k+1)+', instance=Ele'+str(N+1)+'-RVE'+str(j+1))
            RVE_Sets.append(str(ParingFacesLR[k][0]+1))
            RVE_Sets.append('*Nset, nset=Ele'+str(N+1)+'-RVE'+str(j+1)+'-NodeR'+str(k+1)+', instance=Ele'+str(N+1)+'-RVE'+str(j+1))
            RVE_Sets.append(str(ParingFacesLR[k][1]+1))
            
            for l in range(2):
                Eqns.append('** Constraint: Ele'+str(N+1)+'-RVE'+str(j+1)+'-LR'+str(k+1)+'-dof'+str(l+1))
                Eqns.append('*Equation')
                Eqns.append('6')
                Eqns.append('Ele'+str(N+1)+'-RVE'+str(j+1)+'-NodeR'+str(k+1)+','+str(l+1)+',-1.0')
                Eqns.append('Ele'+str(N+1)+'-RVE'+str(j+1)+'-NodeL'+str(k+1)+','+str(l+1)+',1.0')
                Eqns.append('Ele'+str(N+1)+'-N1,'+str(l+1)+','+str(N_GloDeriv[0][0]*dx))
                Eqns.append('Ele'+str(N+1)+'-N2,'+str(l+1)+','+str(N_GloDeriv[1][0]*dx))
                Eqns.append('Ele'+str(N+1)+'-N3,'+str(l+1)+','+str(N_GloDeriv[2][0]*dx))
                Eqns.append('Ele'+str(N+1)+'-N4,'+str(l+1)+','+str(N_GloDeriv[3][0]*dx))
                
        for k in range(len(ParingFacesBT)):
            RVE_Sets.append('*Nset, nset=Ele'+str(N+1)+'-RVE'+str(j+1)+'-NodeB'+str(k+1)+', instance=Ele'+str(N+1)+'-RVE'+str(j+1))
            RVE_Sets.append(str(ParingFacesBT[k][0]+1))
            RVE_Sets.append('*Nset, nset=Ele'+str(N+1)+'-RVE'+str(j+1)+'-NodeT'+str(k+1)+', instance=Ele'+str(N+1)+'-RVE'+str(j+1))
            RVE_Sets.append(str(ParingFacesBT[k][1]+1))
            
            for l in range(2):
                Eqns.append('** Constraint: Ele'+str(N+1)+'-RVE'+str(j+1)+'-BT'+str(k+1)+'-dof'+str(l+1))
                Eqns.append('*Equation')
                Eqns.append('6')
                Eqns.append('Ele'+str(N+1)+'-RVE'+str(j+1)+'-NodeT'+str(k+1)+','+str(l+1)+',-1.0')
                Eqns.append('Ele'+str(N+1)+'-RVE'+str(j+1)+'-NodeB'+str(k+1)+','+str(l+1)+',1.0')
                Eqns.append('Ele'+str(N+1)+'-N1,'+str(l+1)+','+str(N_GloDeriv[0][1]*dy))
                Eqns.append('Ele'+str(N+1)+'-N2,'+str(l+1)+','+str(N_GloDeriv[1][1]*dy))
                Eqns.append('Ele'+str(N+1)+'-N3,'+str(l+1)+','+str(N_GloDeriv[2][1]*dy))
                Eqns.append('Ele'+str(N+1)+'-N4,'+str(l+1)+','+str(N_GloDeriv[3][1]*dy))
                
        for k in range(len(RVE_V)):
            RVE_Sets.append('*Nset, nset=Ele'+str(N+1)+'-RVE'+str(j+1)+'-V'+str(k+1)+', instance=Ele'+str(N+1)+'-RVE'+str(j+1))
            RVE_Sets.append(str(RVE_V[k]+1))
        
        for k in range(2):
            Eqns.append('** Constraint: Ele'+str(N+1)+'-RVE'+str(j+1)+'-V12'+'-dof'+str(k+1))
            Eqns.append('*Equation')
            Eqns.append('6')
            Eqns.append('Ele'+str(N+1)+'-RVE'+str(j+1)+'-V2,'+str(k+1)+',-1.0')
            Eqns.append('Ele'+str(N+1)+'-RVE'+str(j+1)+'-V1,'+str(k+1)+',1.0')
            Eqns.append('Ele'+str(N+1)+'-N1,'+str(k+1)+','+str(N_GloDeriv[0][0]*dx))
            Eqns.append('Ele'+str(N+1)+'-N2,'+str(k+1)+','+str(N_GloDeriv[1][0]*dx))
            Eqns.append('Ele'+str(N+1)+'-N3,'+str(k+1)+','+str(N_GloDeriv[2][0]*dx))
            Eqns.append('Ele'+str(N+1)+'-N4,'+str(k+1)+','+str(N_GloDeriv[3][0]*dx))
        
        for k in range(2):
            Eqns.append('** Constraint: Ele'+str(N+1)+'-RVE'+str(j+1)+'-V23'+'-dof'+str(k+1))
            Eqns.append('*Equation')
            Eqns.append('6')
            Eqns.append('Ele'+str(N+1)+'-RVE'+str(j+1)+'-V3,'+str(k+1)+',-1.0')
            Eqns.append('Ele'+str(N+1)+'-RVE'+str(j+1)+'-V2,'+str(k+1)+',1.0')
            Eqns.append('Ele'+str(N+1)+'-N1,'+str(k+1)+','+str(N_GloDeriv[0][1]*dy))
            Eqns.append('Ele'+str(N+1)+'-N2,'+str(k+1)+','+str(N_GloDeriv[1][1]*dy))
            Eqns.append('Ele'+str(N+1)+'-N3,'+str(k+1)+','+str(N_GloDeriv[2][1]*dy))
            Eqns.append('Ele'+str(N+1)+'-N4,'+str(k+1)+','+str(N_GloDeriv[3][1]*dy))
            
        for k in range(2):
            Eqns.append('** Constraint: Ele'+str(N+1)+'-RVE'+str(j+1)+'-V14'+'-dof'+str(k+1))
            Eqns.append('*Equation')
            Eqns.append('6')
            Eqns.append('Ele'+str(N+1)+'-RVE'+str(j+1)+'-V4,'+str(k+1)+',-1.0')
            Eqns.append('Ele'+str(N+1)+'-RVE'+str(j+1)+'-V1,'+str(k+1)+',1.0')
            Eqns.append('Ele'+str(N+1)+'-N1,'+str(k+1)+','+str(N_GloDeriv[0][1]*dy))
            Eqns.append('Ele'+str(N+1)+'-N2,'+str(k+1)+','+str(N_GloDeriv[1][1]*dy))
            Eqns.append('Ele'+str(N+1)+'-N3,'+str(k+1)+','+str(N_GloDeriv[2][1]*dy))
            Eqns.append('Ele'+str(N+1)+'-N4,'+str(k+1)+','+str(N_GloDeriv[3][1]*dy))
            
        for k in range(2):
            Eqns.append('** Constraint: Ele'+str(N+1)+'-RVE'+str(j+1)+'-V1'+'-dof'+str(k+1))
            Eqns.append('*Equation')
            Eqns.append('5')
            Eqns.append('Ele'+str(N+1)+'-RVE'+str(j+1)+'-V1,'+str(k+1)+',-1.0')
            Eqns.append('Ele'+str(N+1)+'-N1,'+str(k+1)+','+str(N1-0.5*N_GloDeriv[0][0]*dx-0.5*N_GloDeriv[0][1]*dy))
            Eqns.append('Ele'+str(N+1)+'-N2,'+str(k+1)+','+str(N2-0.5*N_GloDeriv[1][0]*dx-0.5*N_GloDeriv[1][1]*dy))
            Eqns.append('Ele'+str(N+1)+'-N3,'+str(k+1)+','+str(N3-0.5*N_GloDeriv[2][0]*dx-0.5*N_GloDeriv[2][1]*dy))
            Eqns.append('Ele'+str(N+1)+'-N4,'+str(k+1)+','+str(N4-0.5*N_GloDeriv[3][0]*dx-0.5*N_GloDeriv[3][1]*dy))
            
        # Including interactions
        for k in range(len(RVEInt)):
            for m in range(len(RVEInt[k])):
                line = RVEInt[k][m]
                if line.count(str(RVEPart)+'.') != 0:
                    line = line.replace(str(RVEPart)+'.','Ele'+str(N+1)+'-RVE'+str(j+1)+'-')  
                if line.count('Interaction:') != 0:
                    line = line.replace('Interaction: ','Interaction: Ele'+str(N+1)+'-RVE'+str(j+1)+'-')
                if line.count('interaction=') != 0:
                    line = line.replace('interaction=','interaction=RVE-')                                 
                RVEPreInts.append(line)

### Process and handle the nodes involved in contact
for i in range(len(MacroContactNodes)):
    for j in range(len(MacroContactNodes[i])):
        NTG = NodeTieGroups[MacroContactNodes[i][j]] # NodeTieGroup label, for convenience
        
        if len(NTG) == 1: 
            # Similar process to Tie removal for 1 node in Tie group
            Macro_Sets.append('*Nset, nset=PreLoadC'+str(i+1)+'-'+str(j+1)+', instance=Macro'+str(NTG[0][0]))
            Macro_Sets.append(str(NTG[0][1]))
            StepsPre.append('** Name: PreLoadC'+str(i+1)+'-'+str(j+1)+' Type: Displacement/Rotation')
            StepsPre.append('*Boundary')
            StepsPre.append('PreLoadC'+str(i+1)+'-'+str(j+1)+',1,1,'+str(Disp_Contact[i][j][0]))
            StepsPre.append('PreLoadC'+str(i+1)+'-'+str(j+1)+',2,2,'+str(Disp_Contact[i][j][1]))
            
        else:
            # Check if any elements involved are adapted
            Adapt_Mark = 0
            for k in range(len(NTG)):
                if NTG[k][0] in Ele_NewDFE2:
                    Adapt_Mark = 1
                    break
            
            # None of the elements involved in the group are adapted
            if Adapt_Mark == 0:
                Macro_Sets.append('*Nset, nset=PreLoadC'+str(i+1)+'-'+str(j+1)+', instance=Macro'+str(NTG[0][0]))
                Macro_Sets.append(str(NTG[0][1]))
                StepsPre.append('** Name: PreLoadC'+str(i+1)+'-'+str(j+1)+' Type: Displacement/Rotation')
                StepsPre.append('*Boundary')
                StepsPre.append('PreLoadC'+str(i+1)+'-'+str(j+1)+',1,1,'+str(Disp_Contact[i][j][0]))
                StepsPre.append('PreLoadC'+str(i+1)+'-'+str(j+1)+',2,2,'+str(Disp_Contact[i][j][1]))
                
            else:
                # Find and delete the Tie constraint
                for k in range(len(Ties)):
                    if ('** Constraint: Tie-'+str(MacroContactNodes[i][j]+1)) in Ties[k]:
                        del Ties[k+2]
                        del Ties[k+1]
                        del Ties[k]
                        break  
                    
                # Call all nodes into a set to preload
                for k in range(len(NTG)):
                    Macro_Sets.append('*Nset, nset=PreLoadC'+str(i+1)+'-'+str(j+1)+', instance=Macro'+str(NTG[k][0]))
                    Macro_Sets.append(str(NTG[k][1]))
                StepsPre.append('** Name: PreLoadC'+str(i+1)+'-'+str(j+1)+' Type: Displacement/Rotation')
                StepsPre.append('*Boundary')
                StepsPre.append('PreLoadC'+str(i+1)+'-'+str(j+1)+',1,1,'+str(Disp_Contact[i][j][0]))
                StepsPre.append('PreLoadC'+str(i+1)+'-'+str(j+1)+',2,2,'+str(Disp_Contact[i][j][1]))

### Remove original BCs in case the sets are removed by PreLoads
for i in range(len(MacroLoadNodes)):
    if len(MacroLoadNodes[i]) == 0:
        # set name is MacroLoadSet[i]
        Start = 0
        End = 0
        for j in range(len(Steps)):
            if ((Steps[j].count(MacroLoadSet[i]+',') != 0) and (Start == 0)):
                Start = j
            if ((Steps[j].count('**') != 0) and (Start != 0)):
                End = j
                break
        for j in reversed(range(Start-2,End)):
            del Steps[j]
                
Steps.extend(StepsPre)

for i in range(len(StepEnd)):
    Steps.append(StepEnd[i])

print('Adaptive pre-load step .inp data preparation completed')
Time1 = time.time() - start_time
Time2 = time.time() - start_time0
print('Script time: %ss'%str(Time1))
print('Total time: %ss'%str(Time2))

print>>Time_log,'Adaptive pre-load step .inp data preparation completed'
print>>Time_log,'Script time: %ss'%str(Time1)
print>>Time_log,'Total time: %ss'%str(Time2)











