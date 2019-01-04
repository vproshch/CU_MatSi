#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 11:51:46 2018

@author: ty
"""
import copy as cp
import numpy as np

num_period=8 #number of periods in SL

##### multilayer #####
infile = ('Si.10.Ge.10x8x8' + '.xyz')

with open(infile, 'r') as fid:
    
    size = int(fid.readline()) #get system size
    total = cp.deepcopy(size)
    fid.readline() #skip comments
    
    data = np.zeros((size,5))
    for i in range(size):
        tmp = fid.readline().strip().split() #read the coordinates
        
        data[i,0] = (i+1) #id
        if tmp[0] == 'Si':
            data[i,1] = 1 #type
        if tmp[0] == 'Ge':
            data[i,1] = 2 #type
        data[i,2] = float(tmp[1]) #x coord
        data[i,3] = float(tmp[2]) #y coord
        data[i,4] = float(tmp[3]) #z coord
#####
     
##### Make Sl #####
data[:,2] = data[:,2]*4/5.431
data[:,3] = data[:,3]*4/5.431
data[:,4] = data[:,4]*4/5.431

data2 = cp.deepcopy(data)

for i in range(num_period-1):
    
    tmp=cp.deepcopy(data2)
    maxcoord = data.max(axis=0)
    tmp[:,2]=tmp[:,2]+(maxcoord[2]+1)
    data=np.append(data,tmp,axis=0)
    i=i+1
#####
    
##### Add baths to each end #####
    infile = ('Si.12x8x8' + '.xyz')

with open(infile, 'r') as fid:
    
    size = int(fid.readline()) #get system size
    total = cp.deepcopy(size)
    fid.readline() #skip comments
    
    si = np.zeros((size,5))
    for i in range(size):
        tmp = fid.readline().strip().split() #read the coordinates
        
        si[i,0] = (i+1) #id
        if tmp[0] == 'Si':
            si[i,1] = 1 #type
        if tmp[0] == 'Ge':
            si[i,1] = 2 #type
        si[i,2] = float(tmp[1]) #x coord
        si[i,3] = float(tmp[2]) #y coord
        si[i,4] = float(tmp[3]) #z coord

si[:,2] = si[:,2]*4/5.431
si[:,3] = si[:,3]*4/5.431
si[:,4] = si[:,4]*4/5.431

infile = ('Ge.12x8x8' + '.xyz')

with open(infile, 'r') as fid:
    
    size = int(fid.readline()) #get system size
    total = cp.deepcopy(size)
    fid.readline() #skip comments
    
    ge = np.zeros((size,5))
    for i in range(size):
        tmp = fid.readline().strip().split() #read the coordinates
        
        ge[i,0] = (i+1) #id
        if tmp[0] == 'Si':
            ge[i,1] = 1 #type
        if tmp[0] == 'Ge':
            ge[i,1] = 2 #type
        ge[i,2] = float(tmp[1]) #x coord
        ge[i,3] = float(tmp[2]) #y coord
        ge[i,4] = float(tmp[3]) #z coord
        
ge[:,2] = ge[:,2]*4/5.431
ge[:,3] = ge[:,3]*4/5.431
ge[:,4] = ge[:,4]*4/5.431

maxcoord = si.max(axis=0)
data[:,2]=data[:,2]+(maxcoord[2]+1)
data=np.append(si,data,axis=0)

maxcoord = data.max(axis=0)
ge[:,2]=ge[:,2]+(maxcoord[2]+1)
data=np.append(data,ge,axis=0)
#####
        
    
##### Print structure #####

data[:,2] = data[:,2]*5.431/4
data[:,3] = data[:,3]*5.431/4
data[:,4] = data[:,4]*5.431/4
                
with open('data.sl', 'w') as f:
    
    maxcoord = data.max(axis=0)
    mincoord = data.min(axis=0)
    xmax = maxcoord[2]
    ymax = maxcoord[3]
    zmax = maxcoord[4]
    xmin = mincoord[2]
    ymin = mincoord[3]
    zmin = mincoord[4]
    
    buff=0.09
    masses = np.zeros((2,1))
    masses[0] = 28.0855 #mass of si
    masses[1] = 72.6400 #mass of ge
     
    f.write(str('LAMMPS DATA FILE\n'))
    f.write('\n' + str(len(data)) + ' atoms\n')
    f.write('\n' + str(len(masses)) + ' atom types\n')
    f.write('\n' + str(xmin-buff)+' '+str(xmax+buff)+' xlo'+' xhi\n')
    f.write(str(ymin-buff)+' '+str(ymax+buff)+' ylo'+' yhi\n')
    f.write(str(zmin-buff)+' '+str(zmax+buff)+' zlo'+' zhi\n')
    f.write('\nMasses\n')
    for i in range(len(masses)):
        f.write('\n' + str(i+1) + ' ' + str(float(masses[i])))
    f.write('\n\nAtoms\n\n')
    for i in range(len(data)-1):
        f.write(str(int(i+1)) + ' ' + str(int(data[i,1])) + ' ' + str(data[i,2]) + ' ' +
                str(data[i,3]) + ' ' + str(data[i,4]) + '\n')
    f.write(str(len(data)) +  ' ' + str(int(data[-1,1])) + ' ' + str(data[-1,2]) + ' ' + 
            str(data[-1,3]) + ' ' + str(data[-1,4]))

data[0:1024,1]=7
data[93184:94208,1]=7    
data[1024:6144,1]=8
data[88064:93184,1]=9
    
example = str('SL_EXAMPLE.xyz')
with open(example, 'w') as fid:
    fid.write(str(len(data))+'\n')
    fid.write('xmax'+'\t0'+'\t0'+'\t0\t'+'ymax'+'\t0'+'\t0'+'\t0\t'+
              'zmax'+'\n')
    for i in range(len(data)):
        if data[i,1] == 1:
            fid.write('Si\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+
                      str(data[i][4])+'\n')
            if i == (len(data)-1):
                fid.write('Si\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+
                          str(data[-1][4]))
        if data[i,1] == 2:
            fid.write('Ge\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+
                      str(data[i][4])+'\n')
            if i == (len(data)-1):
                fid.write('Ge\t'+str(data[-1][2])+'\t'+str(data[-1][3])+
                          '\t'+str(data[-1][4]))
        if data[i,1] == 7:
            fid.write('O\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+
                      str(data[i][4])+'\n')
            if i == (len(data)-1):
                fid.write('O\t'+str(data[-1][2])+'\t'+str(data[-1][3])+
                          '\t'+str(data[-1][4]))
        if data[i,1] == 8:
            fid.write('B\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+
                      str(data[i][4])+'\n')
            if i == (len(data)-1):
                fid.write('B\t'+str(data[-1][2])+'\t'+str(data[-1][3])+
                          '\t'+str(data[-1][4]))
        if data[i,1] == 9:
            fid.write('C\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+
                      str(data[i][4])+'\n')
            if i == (len(data)-1):
                fid.write('C\t'+str(data[-1][2])+'\t'+str(data[-1][3])+
                          '\t'+str(data[-1][4]))
#####