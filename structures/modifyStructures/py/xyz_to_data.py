#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 22:51:26 2018

@author: ty
"""
import copy as cp
import numpy as np
buff = 0.6

infile = raw_input('\nEnter a the name of an XYZ file minus the ".xyz"'
                   ' extension\n>')

outfile = str('data.' + infile)
infile = str(infile + '.xyz')

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
        data[i,1] = 1 #type
        data[i,2] = float(tmp[1]) #x coord
        data[i,3] = float(tmp[2]) #y coord
        data[i,4] = float(tmp[3]) #z coord
###########


data = data[data[:,2].argsort()] #sort by x
maxcoords = data.max(axis=0)
mincoords = data.min(axis=0)
xmax = maxcoords[2]
ymax = maxcoords[3]
zmax = maxcoords[4]
xmin = mincoords[2]
ymin = mincoords[3]
zmin = mincoords[4]
del maxcoords, mincoords

with open(outfile, 'w') as f:
    
    f.write(str('LAMMPS DATA FILE\n'))
    
    f.write('\n' + str(len(data)) + ' atoms\n')
    f.write('\n1 atom types\n')
    f.write('\n' + str(xmin-buff)+' '+str(xmax+buff)+' xlo'+' xhi\n')
    f.write(str(ymin-buff)+' '+str(ymax+buff)+' ylo'+' yhi\n')
    f.write(str(zmin-buff)+' '+str(zmax+buff)+' zlo'+' zhi\n')
    f.write('\nMasses\n' + '\n1 12.0107\n' + '\nAtoms\n\n')
    for i in range(len(data)-1):
        f.write(str(int(i+1)) + ' ' + str(data[i,1]) + ' ' + str(data[i,2]) + ' ' +
                str(data[i,3]) + ' ' + str(data[i,4]) + '\n')
    f.write(str(len(data)) +  ' ' + str(data[-1,1]) + ' ' + str(data[-1,2]) + ' ' + 
            str(data[-1,3]) + ' ' + str(data[-1,4]))
#####################################
