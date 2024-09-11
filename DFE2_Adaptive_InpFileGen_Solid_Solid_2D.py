# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 16:06:59 2023

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

**Note
- This script writes the .inp file for the next analysis based on the information from earlier analyses

"""

start_time = time.time()

f1 = open(JobName+'.inp','w')

for i in range(len(Heading)):
    print>>f1,Heading[i]

for i in range(len(Macro_Parts)):
    print>>f1,Macro_Parts[i]
    
for i in range(len(RVE_Parts)):
    print>>f1,RVE_Parts[i]
    
for i in range(len(Insts)):
    print>>f1,Insts[i]
    
for i in range(len(RVE_Sets)):
    print>>f1,RVE_Sets[i]
    
for i in range(len(Macro_Sets)):
    print>>f1,Macro_Sets[i]
    
for i in range(len(Surfs)):
    print>>f1,Surfs[i]

for i in range(len(Ties)):
    print>>f1,Ties[i]
    
for i in range(len(Eqns)):
    print>>f1,Eqns[i]
    
for i in range(len(Mats)):
    print>>f1,Mats[i]
    
for i in range(len(IntProps)):
    print>>f1,IntProps[i]
    
for i in range(len(Ints)):
    print>>f1,Ints[i]
    
for i in range(len(Steps)):
    print>>f1,Steps[i]

f1.close()

print('.inp file generation for %s has been completed'%str(JobName))
Time1 = time.time() - start_time
Time2 = time.time() - start_time0
print('Script time: %ss'%str(Time1))
print('Total time: %ss'%str(Time2))

print>>Time_log,'.inp file generation for %s has been completed'%str(JobName)
print>>Time_log,'Script time: %ss'%str(Time1)
print>>Time_log,'Total time: %ss'%str(Time2)
