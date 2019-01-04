#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 14:58:06 2018

@author: ty
"""
import numpy as np

velsfile = 'vels.dat'

steps = 2000000 #run time
dt = 0.5e-15
dn = 10 #10 #print frequency
prints = steps/dn #times data is printed
split = 8 #times to split data for averaging
tn = prints/split
win = 0.3

om = np.arange(0,tn)*2*np.pi/(tn*dt*dn) #angular frequency
thz = om/2/np.pi*1e-12
dom = om[1]-om[0]

dosSi = np.zeros((tn,1)).astype(complex)
dosGe = np.zeros((tn,1)).astype(complex)
xdosSi = np.zeros((tn,1)).astype(complex)
xdosGe = np.zeros((tn,1)).astype(complex)
ydosSi = np.zeros((tn,1)).astype(complex)
ydosGe = np.zeros((tn,1)).astype(complex)
zdosSi = np.zeros((tn,1)).astype(complex)
zdosGe = np.zeros((tn,1)).astype(complex)

with open(velsfile, 'r') as fid:
    idsSi = np.zeros((1,1))
    idsGe = np.zeros((1,1))
    for j in range(3):
        fid.readline()
    N = int(fid.readline().strip().split()[0]) #number of atoms
    for j in range(5):
        fid.readline()
    for j in range(N):
        tmp = fid.readline().strip().split()
        if int(tmp[1]) == 1:
            idsSi = np.append(idsSi,j)
        else:
            idsGe = np.append(idsGe,j)
    idsSi = np.delete(idsSi,0)
    idsGe = np.delete(idsGe,0)
    idsSi = idsSi.astype(int)
    idsGe = idsGe.astype(int)
    fid.seek(0)
  
    for i in range(split):
        vels = np.zeros((tn,N*3))
        for j in range(tn):
            data = np.zeros((N,3))
            for k in range(9):
                fid.readline()
            for k in range(N):
                tmp = fid.readline().strip().split()
                data[k,0] = tmp[2]
                data[k,1] = tmp[3]
                data[k,2] = tmp[4]
            vels[j,:] = np.reshape(data,N*3)
            
        velsfft = np.fft.fft(vels,axis=0)*dt*dn
        Dos = (np.multiply(abs(velsfft),abs(velsfft))/
               np.tile(np.multiply(vels,vels).mean(axis=0),(tn,1))/(tn*dt*dn))
        DosGe = np.zeros((tn,len(idsGe)*3)).astype(complex)
        DosSi = np.zeros((tn,len(idsSi)*3)).astype(complex)
        
        DosGe[:,0::3] = Dos[:,idsGe*3] #vx of right atoms
        DosGe[:,1::3] = Dos[:,idsGe*3+1] #vy
        DosGe[:,2::3] = Dos[:,idsGe*3+2] #vz
        xDosGe = DosGe[:,0::3]
        yDosGe = DosGe[:,1::3]
        zDosGe = DosGe[:,2::3]
        
        DosSi[:,0::3] = Dos[:,idsSi*3] #vx of right atoms
        DosSi[:,1::3] = Dos[:,idsSi*3+1] #vy
        DosSi[:,2::3] = Dos[:,idsSi*3+2] #vz
        xDosSi = DosSi[:,0::3]
        yDosSi = DosSi[:,1::3]
        zDosSi = DosSi[:,2::3]
        
        DosGe = DosGe.mean(axis=1)
        xDosGe = xDosGe.mean(axis=1)
        yDosGe = yDosGe.mean(axis=1)
        zDosGe = zDosGe.mean(axis=1)
                
        DosSi = DosSi.mean(axis=1)
        xDosSi = xDosSi.mean(axis=1)
        yDosSi = yDosSi.mean(axis=1)
        zDosSi = zDosSi.mean(axis=1)    
        
        dosSi[:,0] = dosSi[:,0]+DosSi[:]
        dosGe[:,0] = dosGe[:,0]+DosGe[:]
        xdosSi[:,0] = xdosSi[:,0]+xDosSi[:]
        xdosGe[:,0] = xdosGe[:,0]+xDosGe[:]
        ydosSi[:,0] = ydosSi[:,0]+yDosSi[:]
        ydosGe[:,0] = ydosGe[:,0]+yDosGe[:]
        zdosSi[:,0] = zdosSi[:,0]+zDosSi[:]
        zdosGe[:,0] = zdosGe[:,0]+zDosGe[:]
        
dosSi = dosSi/split
dosGe = dosGe/split
xdosSi = xdosSi/split
xdosGe = xdosGe/split
ydosSi = ydosSi/split
ydosGe = ydosGe/split
zdosSi = zdosSi/split
zdosGe = zdosGe/split

## GAUSSIAN SMOOTHING ##
win = round(win*1e12*2*np.pi/dom) #I'm not certain how this works ...
gwin = win-1
n = np.arange(0,gwin+1)-(gwin/2) 
n = (2.5*n/(gwin/2.0))
gauss = np.exp(-(1.0/2.0)*np.multiply(n,n)) #gaussian fx to convolve with
gauss = gauss/np.sum(gauss) #normalized gaussian

sdosGe = np.convolve(dosGe[:,0],gauss,mode='same') #smoothened dos
sdosSi = np.convolve(dosSi[:,0],gauss,mode='same') #smoothened dos
xsdosGe = np.convolve(xdosGe[:,0],gauss,mode='same') #smoothened dos
xsdosSi = np.convolve(xdosSi[:,0],gauss,mode='same') #smoothened dos
ysdosGe = np.convolve(ydosGe[:,0],gauss,mode='same') #smoothened dos
ysdosSi = np.convolve(ydosSi[:,0],gauss,mode='same') #smoothened dos
zsdosGe = np.convolve(zdosGe[:,0],gauss,mode='same') #smoothened dos
zsdosSi = np.convolve(zdosSi[:,0],gauss,mode='same') #smoothened dos


        
        
        
            
        
    
            
            
            
            
            
            
            
            
            
            
            
            
            