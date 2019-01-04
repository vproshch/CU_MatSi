#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 10:39:39 2018

@author: ty
"""
import copy as cp
import numpy as np

infile = 'config.xyz'
etc = 'strained'
buff = 0.6
a = 5.5445

with open(infile, 'r') as fid:
    
    size = int(fid.readline()) #get system size
    total = cp.deepcopy(size)
    fid.readline() #skip comments
    
    data = np.zeros((size,5))
    for i in range(size):
        tmp = fid.readline().strip().split() #read the coordinates
        
        data[i,0] = (i+1) #id
        if i == 0:
            atom_type = tmp[0]
        if tmp[0] == 'Si':
            data[i,1] = 1 #type
        else:
            data[i,1] = 2
        data[i,2] = float(tmp[1]) #x coord
        data[i,3] = float(tmp[2]) #y coord
        data[i,4] = float(tmp[3]) #z coord
###########
size = len(data)

vsi = 5.431**3
vge = 5.658**3

asi = vsi/a**2
age = vge/a**2

si = data[:size/2,:]
ge = data[size/2:,:]



si[:,2] = si[:,2]*asi/5.431
si[:,3] = si[:,3]*a/5.431
si[:,4] = si[:,4]*a/5.431

ge[:,2] = ge[:,2]*age/5.431
ge[:,3] = ge[:,3]*a/5.431
ge[:,4] = ge[:,4]*a/5.431


dxsi = np.unique(si[:,2])[1]-np.unique(si[:,2])[0]
dxge = np.unique(ge[:,2])[1]-np.unique(ge[:,2])[0]

dx = (dxsi+dxge)/2

shift = ge[0,2]-si[-1,2]

ge[:,2] = ge[:,2]-shift+dx


###find max and min coords          
maxcoords = data.max(axis=0)
mincoords = data.min(axis=0)
xmax = maxcoords[2]
ymax = maxcoords[3]
zmax = maxcoords[4]
xmin = mincoords[2]
ymin = mincoords[3]
zmin = mincoords[4]
del maxcoords, mincoords


###########  WRITE TO .xyz FILE    ##################
filename = etc + '.xyz'
with open(filename, 'w') as f:
    f.write(str(len(data))+'\n')
    f.write(str(xmax)+'\t0'+'\t0'+'\t0\t'+str(ymax)+'\t0'+'\t0'+'\t0\t'+str(zmax)+'\n')
    for i in range(len(data)-1):
        if data[i][1] == 1:
            f.write('Si\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+str(data[i][4])+'\n')
        if data[i][1] == 2:
            f.write('Ge\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+str(data[i][4])+'\n')
    if data[-1][1] == 1:
        f.write('Si\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+str(data[-1][4]))
    if data[-1][1] == 2:
        f.write('Ge\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+str(data[-1][4]))
####################################################



###########  WRITE LAMMPS DATA  ###################
    
masses = ['28.08556', '72.6400']
 
datafile = 'data.' + etc
with open(datafile, 'w') as f:
    
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
##################################################
