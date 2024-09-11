# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 10:36:52 2023

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

**Note
- This script extracts and processes the RVE data from the original .inp file
- It then sets up the .inp files for the RVE unit strain analysis
- Once the jobs are submitted, it reads the .odb files and extracts the RVE element strain data as well as effective macroscale properties

"""

start_time = time.time()

Unit_Factor = 100.0 # Scaling factor for the unit strain analyses to avoid excessive deformation on the elements

### Define functions
def SortListofNodes1D(faceN,coordinate): # 0 for x; 1 for y; 2 for z; gives node numbers ready to be called with python (already -1)
    newlist = []
    oldlist = []
    for i in range(len(faceN)):
        oldlist.append(RVENodalCoord[faceN[i]][coordinate])
    
    orderedlist = sorted(oldlist)
    for j in range(len(orderedlist)):
        ind = oldlist.index(orderedlist[j])
        newlist.append(faceN[ind])
    
    return newlist

def TakeVertexOut(face):
    face.pop(0)
    face.pop(-1)
    return face 

### Extracting RVE info from RVE .inp file
# Nodes and elements
RVENodalConnect,RVENodalCoord = [],[]
Search1(RVEInp,'*Part, name='+RVEPart,'*Element, type=CPS4',RVENodalConnect,0)
Search1(RVEInp,'*Part, name='+RVEPart,'*Node',RVENodalCoord,1)

# Direct FE2 node sets
V1,V2,V3,V4 = [],[],[],[]
Search2(RVEInp,'*Part, name='+RVEPart,'*Nset, nset=V1',V1,0)
Search2(RVEInp,'*Part, name='+RVEPart,'*Nset, nset=V2',V2,0)
Search2(RVEInp,'*Part, name='+RVEPart,'*Nset, nset=V3',V3,0)
Search2(RVEInp,'*Part, name='+RVEPart,'*Nset, nset=V4',V4,0)
RVE_V = [int(V1[0]),int(V2[0]),int(V3[0]),int(V4[0])]

FaceLNodes,FaceRNodes,FaceTNodes,FaceBNodes = [],[],[],[] # On FaceBa
Search2(RVEInp,'*Part, name='+RVEPart,'*Nset, nset=FaceL',FaceLNodes,0)
Search2(RVEInp,'*Part, name='+RVEPart,'*Nset, nset=FaceR',FaceRNodes,0)
Search2(RVEInp,'*Part, name='+RVEPart,'*Nset, nset=FaceT',FaceTNodes,0)
Search2(RVEInp,'*Part, name='+RVEPart,'*Nset, nset=FaceB',FaceBNodes,0)

# Unscaled section thickness
Ind = RVEInp.index('*Solid Section, elset='+str(RVE_Mat_Set[0])+', material='+str(RVE_Mat[0]))
if RVEInp[Ind+1] == ',':
    Thickness = 1.0
else:
    Thickness = float(RVEInp[Ind+1].strip(','))

# Start and end of RVE part details, except name and sections which will keep changing
RVEStart = RVEInp.index('*Node')     
for j in range(RVEStart,len(RVEInp)):
    if RVEInp[j].count('** Section:') != 0:
        RVEEnd = j
        break
    
# RVE part surfaces
RVE_Surface_List = []
for i in range(len(RVE_Surface)):
    Mark = 0
    for j in range(len(RVEInp)):
        if (RVEInp[j].count(RVE_Surface[i]) != 0) and (RVEInp[j].count('*Surface') != 0):
            Mark = j
            break # Skipping when reach to next surface is handled by the break in the next loop

    Temp = []
    Temp.append(RVEInp[Mark])
    for j in range(Mark+1,len(RVEInp)):
        if (RVEInp[j] == '') or (RVEInp[j].count("*") != 0):
            break
        Line = ((RVEInp[j].strip('[]')).replace(',',' ')).split()
        RVE_Set.append(Line[0])
        Temp.append(RVEInp[j])
    RVE_Surface_List.append(Temp)
    
# RVE part sets
RVE_Set_List = []
RVE_Set_Type = []
for i in range(len(RVE_Set)):
    for j in range(len(RVEInp)):
        if (RVEInp[j].count(RVE_Set[i]) != 0) and (RVEInp[j].count('set') != 0):
            Mark = j
            if RVEInp[j].count('Nset') != 0:
                Type = 'N'
                break
            elif RVEInp[j].count('Elset') != 0:
                Type = 'El'
                break
            
    Temp = []
    for j in range(Mark+1,len(RVEInp)):
        if (RVEInp[j] == '') or (RVEInp[j].count("*") != 0):
            break
        Temp.append(RVEInp[j])
    RVE_Set_List.append(Temp)
    RVE_Set_Type.append(Type)

# RVE interactions and interaction properties
RVEIntProp = []
if '** INTERACTION PROPERTIES' in RVEInp:
    Start = RVEInp.index('** INTERACTION PROPERTIES')
    for j in range(Start+2,len(RVEInp)):
        if RVEInp[j].count('name') != 0:
            if (j!=(Start+2)):
                RVEIntProp.append(Prop)
            Prop = []
        Prop.append(RVEInp[j])
        if RVEInp[j+1].count('** INTERACTIONS') != 0:
            RVEIntProp.append(Prop)
            break

RVEInt = []
if '** INTERACTIONS' in RVEInp:
    Start = RVEInp.index('** INTERACTIONS')
    for j in range(Start+2,len(RVEInp)):
        if RVEInp[j].count('** Interaction:') != 0:
            if (j!=(Start+2)):
                RVEInt.append(Int)
            Int = []
        line = RVEInp[j]
        if line.count(RVEPart+'-1') != 0:
            line = line.replace(RVEPart+'-1',RVEPart)
        Int.append(line)
        if RVEInp[j+1].count('** ----------------------------------------------------------------') != 0:
            RVEInt.append(Int)
            break

# RVE materials 
RVEMat = []
Start = RVEInp.index('** MATERIALS')
for j in range(Start+2,len(RVEInp)):
    if RVEInp[j].count('*Material') != 0:
        if (j!=(Start+2)):
            RVEMat.append(Mat)
        Mat = []
    Mat.append(RVEInp[j])
    if (RVEInp[j].count('**') != 0) or (j==(len(RVEInp)-1)):
        RVEMat.append(Mat)
        break

### Calculating RVE offsets and dimensions
OffsetX = (RVENodalCoord[RVE_V[0]][0]+RVENodalCoord[RVE_V[1]][0]+RVENodalCoord[RVE_V[2]][0]+RVENodalCoord[RVE_V[3]][0])/4
OffsetY = (RVENodalCoord[RVE_V[0]][1]+RVENodalCoord[RVE_V[1]][1]+RVENodalCoord[RVE_V[2]][1]+RVENodalCoord[RVE_V[3]][1])/4
Offset = [OffsetX,OffsetY]
B_RVE = (-RVENodalCoord[RVE_V[0]][0]+RVENodalCoord[RVE_V[1]][0]+RVENodalCoord[RVE_V[2]][0]-RVENodalCoord[RVE_V[3]][0])/2
H_RVE = (-RVENodalCoord[RVE_V[0]][1]-RVENodalCoord[RVE_V[1]][1]+RVENodalCoord[RVE_V[2]][1]+RVENodalCoord[RVE_V[3]][1])/2

for i in RVENodalCoord:
    for j in range(len(i)):
        i[j] = i[j]-Offset[j]

### Pairing nodes for PBCs later (to be used in subsequent scripts as well)
FaceLNodes = TakeVertexOut(SortListofNodes1D(FaceLNodes,1))
FaceRNodes = TakeVertexOut(SortListofNodes1D(FaceRNodes,1))
ParingFacesLR = []
for i in range(len(FaceLNodes)):
    Temp = []
    Temp.append(FaceLNodes[i])
    Temp.append(FaceRNodes[i])
    ParingFacesLR.append(Temp)

FaceBNodes = TakeVertexOut(SortListofNodes1D(FaceBNodes,0))
FaceTNodes = TakeVertexOut(SortListofNodes1D(FaceTNodes,0))
ParingFacesBT = []
for i in range(len(FaceBNodes)):
    Temp = []
    Temp.append(FaceBNodes[i])
    Temp.append(FaceTNodes[i])
    ParingFacesBT.append(Temp)

### Loads in order: epsilon_x, epsilon_y, gamma_xy
Loads = [[[1/Unit_Factor,0],[0,0]],[[0,1/Unit_Factor],[0,0]],[[0,0],[0.5/Unit_Factor,0.5/Unit_Factor]]]
Load_type = ['e_x','e_y','g_xy']

### Creating and submitting the unit load .inp files
for i in range(len(Loads)):
    JobName = ModelName+'_RVEUnitLoad_'+str(Load_type[i])
    
    f1 = open(JobName+'.inp','w')
    
    a = """*Heading
** Job name: %s Model name: %s
** Generated by: Abaqus/CAE 2017
*Preprint, echo=NO, model=NO, history=NO, contact=NO
**
** PARTS""" %(JobName,JobName)
    print>>f1,a
    print>>f1,'**'
    print>>f1,'*Part, name='+str(RVEPart)
    for j in range(RVEStart,RVEEnd): # Prints the RVE part nodes, elements, sets and surfaces
        print>>f1,RVEInp[j]
    for j in range(len(RVE_Mat_Set)):
        print>>f1,'** Section: Solid_%s' %(str(RVE_Mat[j]))
        print>>f1,'*Solid Section, elset=%s, material=%s' %(str(RVE_Mat_Set[j]),str(RVE_Mat[j]))
        print>>f1,str(Thickness)+','    
    print>>f1,'*End Part'
    a = """**  
**
** ASSEMBLY
**
*Assembly, name=Assembly
**  
*Instance, name=RVE, part=%s
*End Instance
**
*Node
1, 0, 0, 0
*Node
2, 0, 0, 0
*Nset, nset=RP-1
1,
*Nset, nset=RP-2
2,""" %(str(RVEPart))
    print>>f1,a
    
    ### Create existing sets and surfaces
    for j in range(len(RVE_Set_List)):
        if RVE_Set_Type[j] == 'N':
           print>>f1,'*Nset, nset='+str(RVE_Set[j])+', instance='+str(RVEPart)
           for k in range(len(RVE_Set_List[j])):
               print>>f1,RVE_Set_List[j][k]
               
        if RVE_Set_Type[j] == 'El':
           print>>f1,'*Elset, elset='+str(RVE_Set[j])+', instance='+str(RVEPart)
           for k in range(len(RVE_Set_List[j])):
               print>>f1,RVE_Set_List[j][k]
               
    for j in range(len(RVE_Surface_List)):
        for k in range(len(RVE_Surface_List[j])):
            print>>f1,RVE_Surface_List[j][k]
    
    ### Create sets and constraints for PBC
    for j in range(len(ParingFacesLR)):
        print>>f1,'*Nset, nset=RVE-NodeL'+str(j+1)+', instance=RVE'
        print>>f1,str(ParingFacesLR[j][0]+1)
        print>>f1,'*Nset, nset=RVE-NodeR'+str(j+1)+', instance=RVE'
        print>>f1,str(ParingFacesLR[j][1]+1)
        
    for j in range(len(ParingFacesBT)):
        print>>f1,'*Nset, nset=RVE-NodeB'+str(j+1)+', instance=RVE'
        print>>f1,str(ParingFacesBT[j][0]+1)
        print>>f1,'*Nset, nset=RVE-NodeT'+str(j+1)+', instance=RVE'
        print>>f1,str(ParingFacesBT[j][1]+1)
        
    for j in range(len(RVE_V)):
        print>>f1,'*Nset, nset=RVE-V'+str(j+1)+', instance=RVE'
        print>>f1,str(RVE_V[j]+1)
        
    for j in range(len(ParingFacesLR)):        
        print>>f1,'** Constraint: RVE-LR'+str(j+1)+'-dof1'
        print>>f1,'*Equation'
        print>>f1,'3'
        print>>f1,'RVE-NodeR'+str(j+1)+',1,-1.0'
        print>>f1,'RVE-NodeL'+str(j+1)+',1,1.0'
        print>>f1,'RP-1,1,'+str(B_RVE)
        
        print>>f1,'** Constraint: RVE-LR'+str(j+1)+'-dof2'
        print>>f1,'*Equation'
        print>>f1,'3'
        print>>f1,'RVE-NodeR'+str(j+1)+',2,-1.0'
        print>>f1,'RVE-NodeL'+str(j+1)+',2,1.0'
        print>>f1,'RP-2,1,'+str(B_RVE)
        
    for j in range(len(ParingFacesBT)):
        print>>f1,'** Constraint: RVE-BT'+str(j+1)+'-dof1'
        print>>f1,'*Equation'
        print>>f1,'3'
        print>>f1,'RVE-NodeT'+str(j+1)+',1,-1.0'
        print>>f1,'RVE-NodeB'+str(j+1)+',1,1.0'
        print>>f1,'RP-2,2,'+str(H_RVE)
        
        print>>f1,'** Constraint: RVE-BT'+str(j+1)+'-dof2'
        print>>f1,'*Equation'
        print>>f1,'3'
        print>>f1,'RVE-NodeT'+str(j+1)+',2,-1.0'
        print>>f1,'RVE-NodeB'+str(j+1)+',2,1.0'
        print>>f1,'RP-1,2,'+str(H_RVE)
    
    print>>f1,'** Constraint: RVE-V12-dof1'
    print>>f1,'*Equation'
    print>>f1,'3'
    print>>f1,'RVE-V2,1,-1.0'
    print>>f1,'RVE-V1,1,1.0'
    print>>f1,'RP-1,1,'+str(B_RVE)
    
    print>>f1,'** Constraint: RVE-V12-dof2'
    print>>f1,'*Equation'
    print>>f1,'3'
    print>>f1,'RVE-V2,2,-1.0'
    print>>f1,'RVE-V1,2,1.0'
    print>>f1,'RP-2,1,'+str(B_RVE)
    
    print>>f1,'** Constraint: RVE-V14-dof1'
    print>>f1,'*Equation'
    print>>f1,'3'
    print>>f1,'RVE-V4,1,-1.0'
    print>>f1,'RVE-V1,1,1.0'
    print>>f1,'RP-2,2,'+str(H_RVE)
    
    print>>f1,'** Constraint: RVE-V14-dof2'
    print>>f1,'*Equation'
    print>>f1,'3'
    print>>f1,'RVE-V4,2,-1.0'
    print>>f1,'RVE-V1,2,1.0'
    print>>f1,'RP-1,2,'+str(H_RVE)
    
    print>>f1,'** Constraint: RVE-V23-dof1'
    print>>f1,'*Equation'
    print>>f1,'3'
    print>>f1,'RVE-V3,1,-1.0'
    print>>f1,'RVE-V2,1,1.0'
    print>>f1,'RP-2,2,'+str(H_RVE)
    
    print>>f1,'** Constraint: RVE-V23-dof2'
    print>>f1,'*Equation'
    print>>f1,'3'
    print>>f1,'RVE-V3,2,-1.0'
    print>>f1,'RVE-V2,2,1.0'
    print>>f1,'RP-1,2,'+str(H_RVE)
    
    a = """*End Assembly
**
** MATERIALS
**""" 
    print>>f1,a
    for j in range(len(RVEMat)):
        print>>f1,RVEMat[j][0]
        Start = RVEMat[j].index('*Elastic')
        for k in range(Start,len(RVEMat[j])):
            if ((RVEMat[j][k].count('*') != 0) and (k != Start)):
                break
            print>>f1,RVEMat[j][k]
            
    # Interaction properties
    if len(RVEIntProp) != 0:
        print>>f1,'**'
        print>>f1,'** INTERACTION PROPERTIES'
        print>>f1,'**'
        for j in range(len(RVEIntProp)):
            for k in range(len(RVEIntProp[j])):
                print>>f1,str(RVEIntProp[j][k])
    
    # Interactions
    if len(RVEInt) != 0:
        print>>f1,'** INTERACTIONS'
        print>>f1,'**'
        for j in range(len(RVEInt)):
            for k in range(len(RVEInt[j])):
                print>>f1,str(RVEInt[j][k])
    
    a = """** ----------------------------------------------------------------
** 
** STEP: Step-1
** 
*Step, name=Step-1, nlgeom=%s
*Static
1., 1., 1e-05, 1.
** 
** BOUNDARY CONDITIONS
** """%(NLGeom)
    print>>f1,a
    
    # Loads
    print>>f1,'** Name: RP1_Load Type: Displacement/Rotation'
    print>>f1,'*Boundary'
    print>>f1,'RP-1, 1, 1, '+str(Loads[i][0][0])
    print>>f1,'RP-1, 2, 2, '+str(Loads[i][0][1])
    
    print>>f1,'** Name: RP2_Load Type: Displacement/Rotation'
    print>>f1,'*Boundary'
    print>>f1,'RP-2, 1, 1, '+str(Loads[i][1][0])
    print>>f1,'RP-2, 2, 2, '+str(Loads[i][1][1])    
    
    print>>f1,'** Name: V1_RigidBody Type: Displacement/Rotation'
    print>>f1,'*Boundary'
    print>>f1,'RVE-V1, 1, 1, 0'
    print>>f1,'RVE-V1, 2, 2, 0'
    
    a = """** 
** OUTPUT REQUESTS
** 
*Restart, write, frequency=0    
** 
** FIELD OUTPUT: F-Output-1
** 
*Output, field, variable=PRESELECT
** 
** HISTORY OUTPUT: H-Output-1
** 
*Output, history, variable=PRESELECT
*End Step"""
    print>>f1,a

    f1.close()
    
    os.system("abaqus job="+JobName)
    print(JobName+' has been submitted')
    print>>Time_log,JobName+' has been submitted'
        
### Reading the unit load .odb files
RFData = []
RVEVol = H_RVE*B_RVE*Thickness
for i in range(len(Loads)):
    JobName = ModelName+'_RVEUnitLoad_'+str(Load_type[i])
    
    # Check that the .lck file is gone and the .sta file is also present to make sure the job is completed        
    while 1:
        if (os.path.isfile(str(JobName)+'.lck')):
            continue
        
        elif (not(os.path.isfile(str(JobName)+'.sta'))):
            continue
        
        else:
            break 
            
    odb = openOdb(path=JobName+'.odb')
    
    field = odb.steps['Step-1'].frames[-1].fieldOutputs['S']
    
    # Read the RVE element integration point stress data by material Set (easier to sort for comparison later)
    # Output will be RVEData[Material Set][Element][Integration point][Unit Load][Stress val]
    for j in range(len(RVE_Mat_Set)):
        Set = odb.rootAssembly.instances['RVE'].elementSets[RVE_Mat_Set[j].upper()]
        if i==0:
            RVEData[j] = [[[[0,0,0] for z in range(3)] for y in range(4)] for x in range(len(Set.elements))]
        
        if RVE_Yield[j] == 'NA': # Ignore elastic & no yield sections
            continue
        
        for k in range(len(Set.elements)):
            Ele = Set.elements[k]
            
            for m in range(4):
                Val = field.getSubset(region=Ele).values[m]
                RVEData[j][k][m][i][0] = Val.data[0]*Unit_Factor
                RVEData[j][k][m][i][1] = Val.data[1]*Unit_Factor
                RVEData[j][k][m][i][2] = Val.data[3]*Unit_Factor
       
    # Extract macro properties
    TempRF = []
    field = odb.steps['Step-1'].frames[-1].fieldOutputs['RF']
    Set1 = odb.rootAssembly.instances['ASSEMBLY'].nodes[0]
    Set2 = odb.rootAssembly.instances['ASSEMBLY'].nodes[1]
    
    TempRF.append(field.getSubset(region=Set1).values[0].data[0]*Unit_Factor/float(RVEVol))
    TempRF.append(field.getSubset(region=Set1).values[0].data[1]*Unit_Factor/float(RVEVol))
    TempRF.append(field.getSubset(region=Set2).values[0].data[0]*Unit_Factor/float(RVEVol))

    RFData.append(TempRF)    
        
    odb.close()

### Calculating the coefficients for the nonlinearity criterion
E11_Arr = []
E22_Arr = []
E12_Arr = []

# Working out the max and min values for all 3 strains
Max_Vals = []
for i in range(3):
    for j in range(2):
        E = (-1)**(j+1)
        
        RVE_Yield_Index = []
        for k in range(len(RVE_Mat_Set)): # Material group within an RVE
            RVE_Mat_Set_Yield_Index = []
            if RVE_Yield[k] == 'NA': # Ignore elastic, no yield sections
                continue

            for l in range(len(RVEData[k])): # Elements within this material group
                for m in range(4): # Element's integration point
                    S11 = E*RVEData[k][l][m][i][0]
                    S22 = E*RVEData[k][l][m][i][1]
                    S12 = E*RVEData[k][l][m][i][2]
                    
                    # Need to do NL condition check here, and back calculate appropriate N_load
                    S_Mises = math.sqrt(S11**2 - S11*S22 + S22**2 + 3*(S12**2))
                    RVE_Mat_Set_Yield_Index.append(S_Mises/RVE_Yield[k])
                    
            RVE_Yield_Index.append(max(RVE_Mat_Set_Yield_Index))
            
        Yield_Index = max(RVE_Yield_Index)
        Max_Vals.append(E/Yield_Index)
        
# Working out full range of the yield surface        
E11_Arr2 = np.linspace(Max_Vals[0],Max_Vals[1],N_DataPoints)
E22_Arr2 = np.linspace(Max_Vals[2],Max_Vals[3],N_DataPoints)
E12_Arr2 = np.linspace(Max_Vals[4],Max_Vals[5],N_DataPoints)

# E11 and E12
for i in range(len(E11_Arr2)):
    for j in range(len(E12_Arr2)):
        E11 = E11_Arr2[i]
        E12 = E12_Arr2[j]
        E = [E11,0,E12]
        
        RVE_Yield_Index = []
        RVE_Yield_Ele = []
        RVE_E22Pos = []
        RVE_E22Neg = []
        for k in range(len(RVE_Mat_Set)): # Material group within an RVE  len(RVE_Mat_Set)
            RVE_Mat_Set_Yield_Index = []
            RVE_Mat_Set_E22Pos = []
            RVE_Mat_Set_E22Neg = []
            if RVE_Yield[k] == 'NA': # Ignore elastic, no yield sections
                continue
            
            for l in range(len(RVEData[k])): # Elements within this material group
                for m in range(4): # Element's integration point
                    S11a = 0.0
                    S22a = 0.0
                    S12a = 0.0
                    for n in range(3):
                        if n == 1:
                            S11b = RVEData[k][l][m][n][0]
                            S22b = RVEData[k][l][m][n][1]
                            S12b = RVEData[k][l][m][n][2]
                        else: 
                            S11a += E[n]*RVEData[k][l][m][n][0]
                            S22a += E[n]*RVEData[k][l][m][n][1]
                            S12a += E[n]*RVEData[k][l][m][n][2]
                    
                    # For calculating E11-E22 yield surface when E12 = 0
                    S_Mises = math.sqrt(S11a**2 - S11a*S22a + S22a**2 + 3*(S12a**2))
                    RVE_Mat_Set_Yield_Index.append(S_Mises/RVE_Yield[k])
                    
                    # For calculating all other points on yield surface when E12 != 0
                    if (S_Mises/RVE_Yield[k])<0.995:
                        a = S11b**2 - S11b*S22b + S22b**2 + 3*(S12b**2)
                        b = 2*S11a*S11b - S11a*S22b - S11b*S22a + 2*S22a*S22b + 6*S12a*S12b
                        c = S11a**2 - S11a*S22a + S22a**2 + 3*(S12a**2) - (RVE_Yield[k]**2)
                        RVE_Mat_Set_E22Pos.append((-b + math.sqrt(b**2-4*a*c))/(2*a))
                        RVE_Mat_Set_E22Neg.append((-b - math.sqrt(b**2-4*a*c))/(2*a))                    
                    
            RVE_Yield_Index.append(max(RVE_Mat_Set_Yield_Index))
            RVE_Yield_Ele.append(RVE_Mat_Set_Yield_Index.index(max(RVE_Mat_Set_Yield_Index)))
            RVE_E22Pos.append(min(RVE_Mat_Set_E22Pos))
            RVE_E22Neg.append(max(RVE_Mat_Set_E22Neg))
            
        Yield_Index = max(RVE_Yield_Index)
        if Yield_Index >= 0.995:
            E11_Arr.append(E11/Yield_Index)
            E12_Arr.append(E12/Yield_Index)
            E22_Arr.append(0)
        else:
            E11_Arr.append(E11)
            E12_Arr.append(E12)
            E22_Arr.append(min(RVE_E22Pos))
            E11_Arr.append(E11)
            E12_Arr.append(E12)
            E22_Arr.append(max(RVE_E22Neg))
            
# E22 and E12
for i in range(len(E22_Arr2)):
    for j in range(len(E12_Arr2)):
        E22 = E22_Arr2[i]
        E12 = E12_Arr2[j]
        E = [0,E22,E12]
        
        RVE_Yield_Index = []
        RVE_Yield_Ele = []
        RVE_E11Pos = []
        RVE_E11Neg = []
        for k in range(len(RVE_Mat_Set)): # Material group within an RVE  len(RVE_Mat_Set)
            RVE_Mat_Set_Yield_Index = []
            RVE_Mat_Set_E11Pos = []
            RVE_Mat_Set_E11Neg = []
            if RVE_Yield[k] == 'NA': # Ignore elastic, no yield sections
                continue
            
            for l in range(len(RVEData[k])): # Elements within this material group
                for m in range(4): # Element's integration point
                    S11a = 0.0
                    S22a = 0.0
                    S12a = 0.0
                    for n in range(3):
                        if n == 0:
                            S11b = RVEData[k][l][m][n][0]
                            S22b = RVEData[k][l][m][n][1]
                            S12b = RVEData[k][l][m][n][2]
                        else: 
                            S11a += E[n]*RVEData[k][l][m][n][0]
                            S22a += E[n]*RVEData[k][l][m][n][1]
                            S12a += E[n]*RVEData[k][l][m][n][2]
                    
                    # For calculating E11-E22 yield surface when E12 = 0
                    S_Mises = math.sqrt(S11a**2 - S11a*S22a + S22a**2 + 3*(S12a**2))
                    RVE_Mat_Set_Yield_Index.append(S_Mises/RVE_Yield[k])
                    
                    # For calculating all other points on yield surface when E12 != 0
                    if (S_Mises/RVE_Yield[k])<0.995:
                        a = S11b**2 - S11b*S22b + S22b**2 + 3*(S12b**2)
                        b = 2*S11a*S11b - S11a*S22b - S11b*S22a + 2*S22a*S22b + 6*S12a*S12b
                        c = S11a**2 - S11a*S22a + S22a**2 + 3*(S12a**2) - (RVE_Yield[k]**2)
                        RVE_Mat_Set_E11Pos.append((-b + math.sqrt(b**2-4*a*c))/(2*a))
                        RVE_Mat_Set_E11Neg.append((-b - math.sqrt(b**2-4*a*c))/(2*a))                    
                    
            RVE_Yield_Index.append(max(RVE_Mat_Set_Yield_Index))
            RVE_Yield_Ele.append(RVE_Mat_Set_Yield_Index.index(max(RVE_Mat_Set_Yield_Index)))
            RVE_E11Pos.append(min(RVE_Mat_Set_E11Pos))
            RVE_E11Neg.append(max(RVE_Mat_Set_E11Neg))
            
        Yield_Index = max(RVE_Yield_Index)
        if Yield_Index >= 0.995:
            E22_Arr.append(E22/Yield_Index)
            E12_Arr.append(E12/Yield_Index)
            E11_Arr.append(0)
        else:
            E22_Arr.append(E22)
            E12_Arr.append(E12)
            E11_Arr.append(min(RVE_E11Pos))
            E22_Arr.append(E22)
            E12_Arr.append(E12)
            E11_Arr.append(max(RVE_E11Neg))
            
# Fitting the coefficients
A = []
b = []
for i in range(len(E11_Arr)):
    E11 = E11_Arr[i]
    E22 = E22_Arr[i]
    E12 = E12_Arr[i]
    
    A.append([E11,E22,E12,E11**2,E22**2,E12**2,E11*E22,E11*E12,E22*E12,E11*E12*E22])
    b.append([1])

A = np.array(A)
b = np.array(b)
AT = np.transpose(A)
ATA_inv = np.linalg.inv(np.dot(AT,A))
A_inv = np.dot(ATA_inv,AT)
x = np.dot(A_inv,b)
x = np.transpose(x)
NLCrit_Coeff = list(x[0])

### Calculate effective macro elastic properties from the RF and unit strain data        
# RFData[UnitLoadCase][Sig1, Sig2, Sig3]    
G12 = RFData[2][2]    

S12 = 1.0/(RFData[0][1]-((RFData[1][1]*RFData[0][0])/RFData[1][0]))
S11 = -RFData[1][1]/RFData[1][0]*S12
S22 = -RFData[0][0]/RFData[0][1]*S12

E1 = 1.0/S11
E2 = 1.0/S22
v12 = -S12*E1

MacroProps = [E1,E2,v12,G12]
    
print('RVE pre-analysis completed')
Time1 = time.time() - start_time
Time2 = time.time() - start_time0
print('Script time: %ss'%str(Time1))
print('Total time: %ss'%str(Time2))

print>>Time_log,'RVE pre-analysis completed'
print>>Time_log,'Script time: %ss'%str(Time1)
print>>Time_log,'Total time: %ss'%str(Time2)
      




















