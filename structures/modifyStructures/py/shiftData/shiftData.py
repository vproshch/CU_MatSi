#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 10:13:06 2018

@author: ty

mave atoms in Si/Ge interfaces that have drifted to the other side during 
relaxation back to thier original positions

read in the atomic positions. find the ge atoms in the 'leftmost' part and
cut them from the structures. shift the structure to xmin = 0.

find several intermonolayer distances on the right side, i.e. 50% of Ge region.
average the distances and append the ge atoms cut from the si side to the ge
region with the calculated monolayer distance. Write to file.

'data.relaxed' has atom positions followed some other numbers that 
dont matter, and then all atomic velocities. Only read in positions. Velocities
are 're-seeded' during NEMD.

"""

import numpy as np
import copy as cp

infile = 'data.relaxed'
yz = 8 #no of unit cells in y,z direction

no = yz*yz*2 #no. atoms per monolayer

### READ DATA ###
with open(infile, 'r') as fid:
    fid.readline() #skip comments
    fid.readline()
    
    size = int(fid.readline().strip().split()[0]) #Number of atoms
    fid.readline()
    fid.readline()
    
    tmp = fid.readline().strip().split()
    xlo = float(tmp[0])
    xhi = float(tmp[1])
    tmp = fid.readline().strip().split()
    ylo = float(tmp[0])
    yhi = float(tmp[1])
    tmp = fid.readline().strip().split()
    zlo = float(tmp[0])
    zhi = float(tmp[1])
    
    for i in range(8):
        fid.readline() #skip the rest of the comments
    
    data = np.zeros((size,5)) #original positions
    for i in range(size):
        tmp = fid.readline().strip().split()
        data[i,0] = int(tmp[0]) #id
        data[i,1] = int(tmp[1]) #type
        data[i,2] = float(tmp[2]) #x
        data[i,3] = float(tmp[3]) #y
        data[i,4] = float(tmp[4]) #z
        
### SHIFT GE ATOMS ###
data = data[data[:,2].argsort(),:] #sort by x coord

ids = np.argwhere(data[0:(size/3),1] == 2) #Ge atoms in first 3rd of structure
cut = data[ids[:,0],:] #Ge atoms
data = np.delete(data,ids[:,0],axis=0) #cut Ge atoms from structure

dmin = data[:,2].min() #minimum x coord of structure
data[:,2] = data[:,2]-dmin
nsize = len(data) #size of new array
dmax = data[nsize-no:nsize,2].mean() #average coord of last monolayer

dist = 0 #avg interlayer distance
for i in range(nsize/8/no):
    lo1 = nsize-(i+2)*128
    hi1 = nsize-(i+1)*128
    lo2 = nsize-(i+1)*128
    hi2 = nsize-i*128
    tmp = data[lo2:hi2,2].mean()-data[lo1:hi1,2].mean() 
    dist = dist+tmp
    
dist = dist/(nsize/8/no) #avg interlayer distance

del hi1, hi2, lo1, lo2, tmp #clean up variables

if len(cut) > yz*yz/2: #avg of lo coord of ge atoms
    cmin = cut[0:yz*yz/2,2].mean()
else:
    cmin = cut[:,2].mean()
    
if cmin < 0:
    cut[:,2] = cut[:,2]+abs(cmin)+dmax+dist #shift ge atoms
else:
    cut[:,2] = cut[:,2]-abs(cmin)+dmax+dist #shift ge atoms
    
data = np.append(data,cut,axis=0)
data = data[np.lexsort((data[:,4],data[:,3],data[:,2]))] #sort data

### WRITE TO FILE ###
xlo = data[:,2].min()
xhi = data[:,2].max()
vac = 20 #len of vaccuum to add to each end

with open('EXAMPLEshifted.xyz', 'w') as fid: #xyz file to view structure
    fid.write(str(size)+'\n')
    fid.write('comments\n')
    for i in range(size-1):
        if data[i,1] == 1:
            fid.write('Si\t'+str(data[i,2])+'\t'+str(data[i,3])
            +'\t'+str(data[i,4])+'\n')
        if data[i,1] == 2:
            fid.write('Ge\t'+str(data[i,2])+'\t'+str(data[i,3])
            +'\t'+str(data[i,4])+'\n')
    if data[i,1] == 1:
            fid.write('Si\t'+str(data[i,2])+'\t'+str(data[i,3])
            +'\t'+str(data[i,4]))
    if data[i,1] == 2:
            fid.write('Ge\t'+str(data[i,2])+'\t'+str(data[i,3])
            +'\t'+str(data[i,4]))

with open('data.shifted', 'w') as fid: #lammps data file
    fid.write(str('LAMMPS DATA FILE\n\n'))
    fid.write(str(size) + ' atoms\n')
    fid.write('2 atom types\n\n')
    fid.write(str(xlo-vac)+' '+str(xhi+vac)+' xlo'+' xhi\n')
    fid.write(str(ylo)+' '+str(yhi)+' ylo'+' yhi\n')
    fid.write(str(zlo)+' '+str(zhi)+' zlo'+' zhi\n\n')
    fid.write('Masses\n\n')
    fid.write('1 28.0855\n2 72.64\n\n')
    fid.write('Atoms\n\n')
    for i in range(size-1):
        fid.write(str(int(i+1))+' '+str(int(data[i,1]))+' '+str(data[i,2])
        +' '+str(data[i,3])+' '+str(data[i,4])+'\n')
    fid.write(str(size)+' '+str(int(data[-1,1]))+' '+str(data[-1,2])+' '+
              str(data[-1,3])+' '+str(data[-1,4]))






    









 
