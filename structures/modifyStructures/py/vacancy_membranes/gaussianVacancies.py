#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 19:34:34 2018

@author: ty
"""
import numpy as np
import copy as cp
import scipy.stats as st
import sys
import random
import math
pi = 3.1415927


#####  GET INFO AND READ DATA  ###### 
############
filename = raw_input('\nPlease enter a filename: (USAGE: a filename string'
                    ' minus the ".xyz" file extension)\n>')
infile = (filename + '.xyz')

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
#############

#############
xran = raw_input('\nWould you like to select a range in the x-direction to'
                 ' include defects?\n(USAGE: yes or no)\n>')
if (xran != 'yes') and (xran != 'no'):
    sys.exit('\n\tUSAGE ERROR: answer must be yes or no\n')
#############

#############    
if xran == 'yes':
    data = data[data[:,2].argsort()]
    maxcoords = data.max(axis=0)
    mincoords = data.min(axis=0)
    xmax = maxcoords[2]
    xmin = mincoords[2]
    
    dx = abs(xmin) + xmax
    nx = int(round(dx/3.57))
    
    try:
        xlo = float(raw_input('\nThere are ' + str(nx) + ' unit cells in the '
                             'x-direcion.\n\nChoose how many unit cells\nto'
                             ' exclude on the "left" side. (USAGE: an integer'
                             ' from 1 to ' + str(nx)+ ')\n>'))
        xhi = float(raw_input('\nChoose how many unit cells\nto exclude on the'
                            ' "right" side. (USAGE: an integer from 1 to '
                            + str(nx-xlo)+ ')\n>'))
    except ValueError:
        sys.exit('\n\tUSAGE ERROR: enter integer values from 1 to ' + str(nx)
                 + '.')
        
    xlo = 0 + 0.36 + xlo*3.57
    xhi = (xmax - 0.35) - xhi*3.57
    
    len1 = (data[:,2] <= xlo+0.1).sum()
    len2 = (data[:,2] >= xhi-0.1).sum()
    ids = [0]*(len1+len2)
    ids[:len1] = range(len1)
    ids[-len2:] = range((size-len2),size)
    exclude = data[ids,:]
    data = np.delete(data, ids, axis=0)
    size = int(len(data))
else:
    print ('\n\tNo range selected. System size is ' + str(size) + ' atoms.')
#############

#############    
try:
    v = int(raw_input('\nHow many vacancies would you like? (USAGE: an integer'
                      ' between 1 and ' + str(size) + ' or 0 for none)\n>'))
    v2 = int(raw_input('\nHow many vacancie pairs would you like? NOTE: total'
                       ' number of vacancies will be twice the number of'
                       ' pairs.\n(USAGE: an integer between 1 and ' +
                       str(size-v) + ' or 0 for none)\n>'))
except ValueError:
    sys.exit('\n\tUSAGE ERROR: answer must be a positive integer\n')
if (v < 0) or (v2 < 0) or ((v + 2*v2) >= size):
    sys.exit('\n\tUSAGE ERROR: Total number of vacancy centers must'
             '\n\tbe between 0 and ' + str(size) + '\n')
#############

############
outfile = (filename + '.' + 'v.' + str(v) + '.v2.' + str(v2) + '.xyz')
dataout = ('data.' + filename + '.' + 'v.' + str(v) + '.v2.' + str(v2))
data = data[data[:,4].argsort()] #sort by z

###collect atoms in each z layer
layerid = np.zeros((len(data),1)) #number of atoms in each layer

end = np.zeros((1,5))
end[0,4] = -1

data = np.append(data,end,axis=0)

j = 0
k = 1
while True: #divide structure into half unit cells, z-dir
    z_min = data[j,4] #find lowest z coord
    if data[j,4] == -1: #if at end of data, exit
        break
    layerid[j] = k #index of layer
    if data[j+1,4] != z_min: #if at next layer, next
        k = k+1
    j = j+1 #index over all atoms

num = int(layerid[-1])/2 #number of layers
layers = np.zeros((num,1)) 
data = np.delete(data, (-1), axis=0)

for i in range(num): #count number of atoms in each layer, half unit cells
    tmp = ((layerid[:] == i+1).sum() + (layerid[:] == i+2).sum())
    layers[i] = tmp
    i = i+1
del end, layerid, z_min
####################################





#####   GENERATE GAUSSIAN PROFILE  ######
vs = v+v2
num_v = np.zeros((num/2,1))
std = vs/2.0 #std dev of standard normal dist.

cdf = 0
dx = (4.0/((num/2))) #size of dev from mean

for i in range(num/2): 
    x = dx+i*dx
    tmp = st.norm.cdf(x)-0.5
    d_cdf = tmp-cdf
    cdf = tmp
    num_v[i] = int(round(d_cdf*vs)) #num of vacancies per layer
del cdf, d_cdf, dx, std, x
####################################





#####   REMOVE VACANCIES FROM DATA   ####
inds = [0]*(int(num_v.sum()*2))
j = 0 #atom index on bottom of bottom layer
k = int(layers[-1]) #index on top of bottom layer
tmp = 0
for i in range(len(num_v)):
    l = size-k #index of bottom of top layer
    m = size-j #index of top of top layer
    rand1 = np.array(random.sample(range(int(j),int(k)),int(num_v[i])))
    rand2 = np.array(random.sample(range(int(l),int(m)),int(num_v[i])))
    rand1 = np.append(rand1,rand2,axis=0)
    for n in range(len(rand1)):
        inds[tmp+n] = rand1[n]
    j = j + int(layers[-1]) #iterate to next layer
    k = k + int(layers[-1])
    tmp = tmp+len(rand1)

vacs = np.zeros((len(inds),5))
vacs = data[inds,:]
data = np.delete(data, inds, axis=0)
#del rand1, rand2, layers


inds = random.sample(range(len(vacs)),v2)

vpairs = np.zeros((v2,5)) #remove from vacs to form pairs
vpairs = vacs[inds,:]
vacs = np.delete(vacs, inds, axis=0)
#########################################





#####   FIND NN's AND FORM PAIRS   #####
if v2 != 0:
    nn_dist = np.zeros((len(data),2))
    #inds = np.zeros((len(vpairs),1))
    inds = [0]*len(vpairs)
    for i in range(len(vpairs)):
        for j in range(len(data)):
            nn_dist[j,0] = j
            nn_dist[j,1] = (math.sqrt((vpairs[i,2]-data[j,2])**2 + 
                   (vpairs[i,3]-data[j,3])**2 + (vpairs[i,4]-data[j,4])**2))
        nn_dist = nn_dist[nn_dist[:,1].argsort()]
        nn = (nn_dist[:,1] <= 1.6).sum()
        if nn == 0:
            nn = (nn_dist[:,1] <= 1.65).sum()
        rand = random.sample(range(nn),1)
        inds[i] = int(nn_dist[rand,0])
    del rand, nn_dist, nn
    
    tmp = np.zeros((len(inds),5))
    tmp = data[inds,:]
    vpairs = np.append(vpairs, tmp, axis=0)
    data = np.delete(data, inds, axis=0)
########################################





#####    WRITE DATA TO FILE    #########
if xran == 'yes':
    data = np.append(data,exclude,axis=0)

## SET VACS AS TYPE 2 TO SHOW THEM AS 'Si' ##
vacs[:,1] = 2
vpairs[:,1] = 3
data_ex = np.append(data,vacs,axis=0)
data_ex = np.append(data_ex,vpairs,axis=0)

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

example = str('EXAMPLE.' + outfile)
with open(example, 'w') as fid:
    fid.write(str(len(data_ex))+'\n')
    fid.write(str(xmax)+'\t0'+'\t0'+'\t0\t'+str(ymax)+'\t0'+'\t0'+'\t0\t'+
              str(zmax)+'\n')
    for i in range(len(data_ex)):
        if data_ex[i,1] == 1:
            fid.write('C\t'+str(data_ex[i][2])+'\t'+str(data_ex[i][3])+'\t'+
                      str(data_ex[i][4])+'\n')
            if i == (len(data_ex)-1):
                fid.write('C\t'+str(data_ex[-1][2])+'\t'+str(data_ex[-1][3])+'\t'+
                          str(data_ex[-1][4]))
        if data_ex[i,1] == 2:
            fid.write('Si\t'+str(data_ex[i][2])+'\t'+str(data_ex[i][3])+'\t'+
                      str(data_ex[i][4])+'\n')
            if i == (len(data_ex)-1):
                fid.write('Si\t'+str(data_ex[-1][2])+'\t'+str(data_ex[-1][3])+
                          '\t'+str(data_ex[-1][4]))
        if data_ex[i,1] == 3:
            fid.write('Ge\t'+str(data_ex[i][2])+'\t'+str(data_ex[i][3])+
                      '\t'+str(data_ex[i][4])+'\n')
            if i == (len(data_ex)-1):
                fid.write('Ge\t'+str(data_ex[-1][2])+'\t'+
                          str(data_ex[-1][3])+'\t'+str(data_ex[-1][4]))

with open(outfile, 'w') as fid:
    fid.write(str(len(data))+'\n')
    fid.write(str(xmax)+'\t0'+'\t0'+'\t0\t'+str(ymax)+'\t0'+'\t0'+'\t0\t'+
              str(zmax)+'\n')
    for i in range(len(data)):
        fid.write('C\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+
                  str(data[i][4])+'\n')
        if i == (len(data)-1):
            fid.write('C\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+
                      str(data[-1][4]))
            
with open(dataout, 'w') as f:
    buff = 0.09
    
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

del xmax, xmin, ymax, ymin, zmax, zmin
print ('\n\n\t----------------------\n\tFinished with no errors!\n\n\tFind "'
       + outfile + '", "' + example + '",\n\t and "' + dataout + 
       '" in your current directory.\n')
######################################
   
