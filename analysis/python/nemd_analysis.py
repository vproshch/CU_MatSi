#!/usr/bin/env python
# -*- coding: utf-8 -*-

#PAULS ATTEMPT AT CONSOLIDATING THE SCRIPTS


##IMPORTING MODULES
import datetime
import os
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
plt.switch_backend('TKagg')
import scipy.stats

##DETERMINE WHICH SCRIPT TO EXECUTE (kappa or KR)
script = raw_input("Which script to execute (0 for kappa, 1 for KR): ")



##PREPARING NEEDED INDICES
print "Check if restart timestep value needs to be altered" 
restart_timestep = 0

# Check if all the input is provided 
if len(sys.argv ) <1 :
	sys.exit ( "Usage: executable timestart_arrayindex(default:1) heatfluxcolumnO(default:11) heatfluxcolumnI(default:12)")

if len(sys.argv) > 2:
	timestart = int(sys.argv[1])
	indexO = int(sys.argv[2])
	indexI = int(sys.argv[3])
else:
	timestart = 1
	indexO = 1
	indexI = 2

	LammpsHeatO_index = 11
	LammpsHeatI_index = 12

if len(sys.argv) > 4:
   oldThickn = float(sys.argv[4])
   newThickn = float(sys.argv[5])
else:
   oldThickn = 1.00
   newThickn = 1.00



##GETTING HEATFLUX DATA
if os.path.isfile('./heatflux'):
	print 'heatflux file exists'
else:
	readHeatFlux = 0;
	print 'heatflux file does not exist'
	heatfluxFile = open("heatflux","a+")
	log_lammps = open('log.lammps')
	count = 0;
	for line in log_lammps:
		#print line;
		if (line.isspace()==0):
			
			log_lines = line.split()
			#print str(log_lines[0])
			if (log_lines[0] == "Loop"):
				readHeatFlux = 0;
			
			if (readHeatFlux == 1) and (line.isspace()==0):
			
				if (count > 0):
					heatfluxFile.write('\n')
				
				#print line
				line = log_lines[0] + " " + log_lines[LammpsHeatO_index] + " " + log_lines[LammpsHeatI_index]
				heatfluxFile.write(line.rstrip(' \t\n\r'));
				count = count + 1;


			if (log_lines[0] == "Step"):
					try:
						if (log_lines[11] == "HFluxO" or log_lines[11] =="v_HFluxO" ): 
							print "found data"
							readHeatFlux = 1;
					except IndexError:
						continue;
					
		
	
	heatfluxFile.close()


##PUTTING HEATFLUX AND TIMESTEPS INTO VECTORS
tmp = np.loadtxt('heatflux')
#tmp = np.load('heatflux')

np.save('heatflux.npy',tmp)
tmp = np.load('heatflux.npy')

time = tmp[:,0]-restart_timestep
firstprintJ = time[timestart]
printintervalJ = time[timestart+1]-time[timestart]
heatfluxO = tmp[:,indexO]*(tmp[:,0]/time)
heatfluxI = tmp[:,indexI]*(tmp[:,0]/time)

print "Test for reading data; check for errors!"
print firstprintJ, printintervalJ, heatfluxO[4], heatfluxI[4]

# Done

# Selection of stationary region and block averaging

print "Note the beginning and end of data region suitable for block averaging"

plt.plot(time,heatfluxO)
plt.show()

begintime = float(raw_input("Beginning of stationary current? "))
endtime = float(raw_input("End of data? "))
numblock = int(raw_input("Number of blocks for averaging? "))

intervalJ = (endtime-begintime)/numblock
linein = ((begintime-firstprintJ)/printintervalJ)+1
print intervalJ, linein

with open('heatflux.blockavg','w+') as f:
	for j in range(0, numblock):
		lineprev = linein+(j*intervalJ)/printintervalJ
		linerun = linein+((j+1)*intervalJ)/printintervalJ
		lineprev=int(lineprev)
		linerun=int(linerun)
		meanJO = np.mean(heatfluxO[lineprev:linerun], axis=None)
		stdJO = np.std(heatfluxO[lineprev:linerun], axis=None, dtype=np.float64,ddof=1)
		meanJI = np.mean(heatfluxI[lineprev:linerun], axis=None)
		stdJI = np.std(heatfluxI[lineprev:linerun], axis=None, dtype=np.float64,ddof=1)
		print lineprev, linerun
		print 'j', j, begintime+j*intervalJ, begintime+(j+1)*intervalJ, 'meanJO',meanJO, 'stdO', stdJO, 'meanJI',meanJI, 'stdI', stdJI
		print >> f, j, begintime+j*intervalJ, begintime+(j+1)*intervalJ, meanJO, stdJO, meanJI, stdJI



#------------- Temperature block averaging -----------------------------

#numslice = int(raw_input("Number of bins in temperature file? "))
tmp_profile = str('tmp.profile.'+raw_input("suffix of temperature profile file(e.g. 0, 1): "))
# input file
tmp = open(tmp_profile,'r')
nslice = -1
numdatablock = -1
numlines = -1

bintmp = []
temperature = []


#added by Paul to automatically get number of bins
gotBins = 0;


for line in tmp:

	entries = line.split()
	if (entries[0][0] == '#'): #this makes us skip lines that start with # (comments)
		continue
	
	
	numlines = numlines + 1
	
#	print entries, len(entries)

	if (gotBins == 0):
		#print entries
		numslice = int(entries[1])
		gotBins = 1
		temperatureSlice = range(1,numslice+1) #Counter array for temperature bins
		print 'number of bins from file', numslice

	if len(entries) == 3:
#		print line,
		nslice = -1
		numdatablock = numdatablock + 1
		if numdatablock == 0:
			firstprint = int(entries[0])
		if numdatablock == 1:
			printinterval = int(entries[0])-firstprint
	else:
		nslice = nslice + 1
#		 print 'nslice',nslice, 'numdatablock',numdatablock, entries[3]
		bintmp.append(entries[1])
		temperature.append(entries[3])
print 'nslice',nslice, 'numdatablock',numdatablock, 'numlines', numlines
a = np.zeros(len(bintmp))
b = np.zeros(len(bintmp))
for j in range(0, len(bintmp)-1):
	a[j] = bintmp[j]
	b[j] = temperature[j]
a = np.reshape(a,(numdatablock+1, nslice+1),order='C')
b = np.reshape(b,(numdatablock+1, nslice+1),order='C')

print firstprint, printinterval
blockin = ((begintime-firstprint)/printinterval)
blockout = ((endtime-firstprint)/printinterval)
interval = (blockout-blockin)/numblock

print blockin, blockout, interval

for l in range(0, numblock):
	blockprev = blockin+(l*interval)
	blockrun = blockin+((l+1)*interval)
#	print blockprev, blockrun
	with open(str(str(blockprev)+'_'+str(blockrun)),'w+') as f:	
		for j in range(0,numslice):
			tmptemp = []
			for k in range(int(blockprev),int(blockrun)):
#				print j, k,b[k,j]
				tmptemp.append(b[k,j])
#			print len(tmptemp)
#			print a[5,j],np.mean(tmptemp,dtype=np.float64,axis=None), np.std(tmptemp, axis=None, dtype=np.float64,ddof=1)
			print >> f, a[5,j],np.mean(tmptemp,dtype=np.float64,axis=None), np.std(tmptemp, axis=None, dtype=np.float64,ddof=1)
# 5 is any block, it should be same for all blocks







if script == "0": #KAPPA SCRIPT
    print "kappa script"
    #--------Fitting of temperature profile ------------------------

    slope = []
    intercept = []

    #celllength = float(raw_input("Length of cell in real units? "))
    celllength = 1.0 #Since we asked the output of tmp.profile to be printed in units box it is in real units and not scaled

    with open(str('temperature_profile_fit'),'w+') as f:	
	    for j in range(0,numblock):
		    tempkbegin = j*interval+blockin
		    tempkend = (j+1)*interval+blockin
		    print tempkbegin, tempkend
		    filename = str(str(tempkbegin)+'_'+str(tempkend))
		    temperature = np.loadtxt(filename)
		    if j == 0:
			    print "Note the beginning and end of flat region for fitting from Figure"
    #			plt.errorbar(numslice*temperature[:,0],temperature[:,1],temperature[:,2],color='black',marker='o')
			    plt.errorbar(temperatureSlice[:],temperature[:,1],temperature[:,2],color='black',marker='o')
			    plt.grid(b=None, which='major',axis='both')
			    plt.show()
			    slicebegin = int(raw_input("Beginning of temperature slice? "))
			    sliceend = int(raw_input("End of temperature slice? "))		
		    x = []
		    y = []
		    for l in range(slicebegin, sliceend):
			    x.append(celllength*temperature[l,0])
			    y.append(temperature[l,1])
		    p = np.polyfit(x,y,1)
		    plt.errorbar(celllength*temperature[:,0],temperature[:,1],temperature[:,2],color='black',marker='o')
		    plt.grid(b=None, which='major',axis='both')
		    plt.xlabel('length(Angstrom)')
		    plt.ylabel('Temperature(K)')
		    plt.plot(x,np.polyval(p,x),color='red')
		    plt.show()
		    slope.append(p[0])
		    intercept.append(p[1])
		    print slope[j], intercept[j]
		    print >> f, slope[j], intercept[j]
	    print >> f, np.mean(slope,dtype=np.float64,axis=None), np.mean(intercept,dtype=np.float64,axis=None)
	    print np.mean(slope,dtype=np.float64,axis=None), np.mean(intercept,dtype=np.float64,axis=None)

    #----------------Calculation of thermal conductivity--------------------------------------------

    meangradient = np.mean(slope,dtype=np.float64,axis=None)
    stdgradient = np.std(slope, axis=None, dtype=np.float64,ddof=1) 

    dataJ = np.loadtxt('heatflux.blockavg')
    Jblocks1 = dataJ[:,3]
    Jblocks2 = -dataJ[:,5]
    meanJ1 = np.mean(Jblocks1,dtype=np.float64,axis=None)
    stdJ1 = np.std(Jblocks1, axis=None, dtype=np.float64,ddof=1)
    print meanJ1, stdJ1
    meanJ2 = np.mean(Jblocks2,dtype=np.float64,axis=None)
    stdJ2 = np.std(Jblocks2, axis=None, dtype=np.float64,ddof=1)
    print meanJ2, stdJ2
    kappa1 = (meanJ1/meangradient)*(1e-10)
    kappastd1 = kappa1*np.sqrt((stdJ1/meanJ1)**2+(stdgradient/meangradient)**2)
    print kappa1, kappastd1
    kappa2 = (meanJ2/meangradient)*(1e-10)
    kappastd2 = kappa2*np.sqrt((stdJ2/meanJ2)**2+(stdgradient/meangradient)**2)
    print kappa2, kappastd2

    #----------------Combine two measurements to get composite kappa value--------------
    ## Source: http://www.burtonsys.com/climate/composite_standard_deviations.html

    kappa = (kappa1+kappa2)/2
    sum_of_error_of_sum_squares = (kappastd1**2)*(numblock-1)+(kappastd2**2)*(numblock-1)
    overall_group_sum_of_squares = ((kappa-kappa1)**2)*numblock + ((kappa-kappa2)**2)*numblock
    kappastd = np.sqrt((sum_of_error_of_sum_squares+overall_group_sum_of_squares)/(2*numblock-1))

    print 'kappa: ', kappa, '+/-',kappastd

    date = datetime.datetime.today()


    kappaData = open("kappaData","a+")

    kappaData.write(str(kappa) +" +/- " + str(kappastd) +  "	begin timestep: " +str(begintime) +"   "+ "end timestep: "+ str(endtime) + "   date: " + str(date)+'\n')

    kappaData.close();





if script == "1":  ##KR SCRIPT

    #-------------Interface positions --------------------------------------

    #numInterfaces = (int(raw_input("How many interfaces? ")))
    numInterfaces = 1
    numInterfaceMarkers = numInterfaces + 2 #Boundaries
    interfacePositions = np.zeros(numInterfaceMarkers)

    meansurfaceT = np.zeros((numblock, numInterfaceMarkers))
    deltaT = np.zeros((numblock, numInterfaces))
    temp = np.zeros(numblock)
    temp1 = np.zeros(numblock)


    #-------------Interface temperature data--------------------------------

    #slope = []
    #intercept = []

    #celllength = float(raw_input("Length of cell in real units? "))

    with open(str('temperature_profile_fit'),'w+') as f:	
	    for j in range(0,numblock):
		    slope = np.zeros(numInterfaceMarkers-1)
		    intercept = np.zeros(numInterfaceMarkers-1)
		    tempkbegin = j*interval+blockin
		    tempkend = (j+1)*interval+blockin
    #		print tempkbegin, tempkend
		    filename = str(str(tempkbegin)+'_'+str(tempkend))
		    temperature = np.loadtxt(filename)
		    if j == 0:
    #			plt.errorbar(numslice*temperature[:,0],temperature[:,1],temperature[:,2],color='black',marker='o')
			    plt.errorbar(temperatureSlice[:],temperature[:,1],temperature[:,2],color='black',marker='o')
			    plt.grid(b=None, which='major',axis='both')
			    plt.show()
						
			
			    intermediate = []
			    interfacePositions=[]
			    for n in range(0,numInterfaces+4):								#collects the numbers input by the user
				    print "Marker1: Left edge; Marker 2: Left Intermediate; Marker3: Middle of Interface; Marker 4: Right intermediate; Marker5: Right edge"
				    if (n == 1 or n == 3): #if intermediate points
					    intermediate.append (int(raw_input("interface_bin_marker ")) )
				    else:
					    interfacePositions.append( int(raw_input("interface_bin_marker ")) )
				    #print n, interfacePositions[n]

			    xpoints = np.zeros(numInterfaceMarkers)
			    xpoints[0] = (temperature[interfacePositions[0],0]+temperature[interfacePositions[0]+1,0])/2 #only if boundary is at the ends of simulation cell
			    for n in range(1,numInterfaceMarkers-1):
				    xpoints[n] = (temperature[interfacePositions[n]-2,0]+temperature[interfacePositions[n]+1,0])/2
			    xpoints[numInterfaceMarkers-1] = temperature[interfacePositions[numInterfaceMarkers-1]-1,0] #same as before
			    print "xpoints", xpoints
			
			    # In case the interface layers needs to be adjusted
			    check = int(raw_input("Need to change interface marker positions? (y/n->1/0)"))
			    shiftl = np.zeros(numInterfaceMarkers)
			    shiftr = np.zeros(numInterfaceMarkers)
			    if check == 1:
				    shiftl = np.zeros(numInterfaceMarkers)
				    shiftr = np.zeros(numInterfaceMarkers)
				    for n in range(0,numInterfaceMarkers-1):
					    shiftl[n] = int(raw_input("interface_bin_shift_left "))
					    shiftr[n] = int(raw_input("interface_bin_shift_right "))
					    print n, interfacePositions[n]
			    else:
				    shiftr = [0,int(intermediate[1] - interfacePositions[1]),0]
				    print shiftr
				    shiftl = [0, int(interfacePositions[1] - intermediate[0]), 0]
				    print shiftl

    #			Optional choice of region for temperature fitting
    #			slicebegin = int(raw_input("Beginning of temperature slice? "))
    #			sliceend = int(raw_input("End of temperature slice? "))		

		    for m in range(0,numInterfaceMarkers-1):	
			    x = []
			    y = []
    #			for l in range(int(interfacePositions[m]+shiftr[m]-1), int(interfacePositions[m+1]-shiftl[m+1]-1)):
			    for l in range(int(interfacePositions[m]+shiftr[m]-1), int(interfacePositions[m+1]-shiftl[m+1])):
				    x.append(temperature[l,0])
				    y.append(temperature[l,1])
			    p = np.polyfit(x,y,1)	
			    slope[m] = p[0]
			    intercept[m] = p[1]
			    plt.errorbar(temperature[:,0],temperature[:,1],temperature[:,2],color='black',marker='o')
			    plt.grid(b=None, which='major',axis='both')
			    plt.xlabel('length(scaled to box length)')
			    plt.ylabel('Temperature(K)')
			    plt.plot(x,np.polyval(p,x),color='red')
			    plt.show()

		    for n in range(0,numInterfaceMarkers-2):
			    deltaT[j,n] = abs((slope[n]-slope[n+1])*xpoints[n+1]+(intercept[n]-intercept[n+1]))
			    meansurfaceT[j,n] = 0.5*abs((slope[n]+slope[n+1])*xpoints[n+1]+(intercept[n]+intercept[n+1]))
		    temp[j] = np.mean(deltaT[j,:],dtype=np.float64,axis=None) 
    #		temp1[j] = np.std(deltaT[j,:], axis=None, dtype=np.float64,ddof=1) 
    #		Multiple interfaces
    #		print tempkbegin, tempkend, deltaT[j,0], deltaT[j,1], deltaT[j,3], deltaT[j,4], temp[j], temp1[j]
    #		print >> f, tempkbegin, tempkend, deltaT[j,0], deltaT[j,1], deltaT[j,3], deltaT[j,4], temp[j], temp1[j]
		    print tempkbegin, tempkend, deltaT[j,0], temp[j]
		    print >> f, tempkbegin, tempkend, deltaT[j,0], temp[j]

		    plt.errorbar(temperature[:,0],temperature[:,1],temperature[:,2],color='black',marker='o')
		    plt.grid(b=None, which='major',axis='both')
		    plt.xlabel('length(Angstrom)')
		    plt.ylabel('Temperature(K)')
		    plt.plot(x,np.polyval(p,x),color='red')
		    plt.show()

    #		plt.errorbar(xpoints[1:numInterfaces], meansurfaceT[j,:], deltaT[j,:], color='red',marker='o')
    #		plt.show()

    #meansurfaceTcold = np.mean(np.concatenate([meansurfaceT[:,0],meansurfaceT[:,4]]),dtype=np.float64,axis=None)
    #meansurfaceThot = np.mean(np.concatenate([meansurfaceT[:,1],meansurfaceT[:,3]]),dtype=np.float64,axis=None)
    #print meansurfaceTcold, meansurfaceThot


    #--------------Calculation of thermal conductance-----------------------


    #kappa1 = (meanJ1/meangradient)*(1e-10)
    #kappastd1 = kappa1*np.sqrt((stdJ1/meanJ1)**2+(stdgradient/meangradient)**2)
    #print kappa1, kappastd1
    #kappa2 = (meanJ2/meangradient)*(1e-10)
    #kappastd2 = kappa2*np.sqrt((stdJ2/meanJ2)**2+(stdgradient/meangradient)**2)
    #print kappa2, kappastd2


    dataJ = np.loadtxt(str('heatflux.blockavg'))
    #JblocksO = dataJ[:,3]/1e6/2
    #JblocksI = -dataJ[:,5]/1e6/2 #Need to lookup for factor 2 in reference
    JblocksO = dataJ[:,3]/1e6
    JblocksI = -dataJ[:,5]/1e6
    meanJO = np.mean(JblocksO,dtype=np.float64,axis=None)
    stdJO = np.std(JblocksO, axis=None, dtype=np.float64,ddof=1)
    meanJI = np.mean(JblocksI,dtype=np.float64,axis=None)
    stdJI = np.std(JblocksI, axis=None, dtype=np.float64,ddof=1)
    #print meanJO, stdJO, meanJI, stdJI

    temps = np.concatenate([deltaT[:,0]])
    #temps = np.concatenate([deltaT[:,0], deltaT[:,1], deltaT[:,3], deltaT[:,4]]) For multiple interfaces
    meanT = np.mean(temps,dtype=np.float64,axis=None)
    stdT = np.std(temps, axis=None, dtype=np.float64,ddof=1)
    print "Interface temperature:", meanT, stdT

    KRO = (meanJO/meanT)
    KROstd = KRO*np.sqrt((stdJO/meanJO)**2+(stdT/meanT)**2)
    print 'KRO', KRO, '+/-', KROstd

    KRI = (meanJI/meanT)
    KRIstd = KRO*np.sqrt((stdJI/meanJI)**2+(stdT/meanT)**2)
    print 'KRI', KRI, '+/-', KRIstd

    inverseJO = 1/meanJO
    inverseJOstd = inverseJO*stdJO/meanJO

    inverseJI = 1/meanJI
    inverseJIstd = inverseJI*stdJI/meanJI

    #tempcold = np.concatenate([deltaT[:,0],deltaT[:,4]])
    #temphot = np.concatenate([deltaT[:,1],deltaT[:,3]])
    #meantempcold = np.mean(tempcold,dtype=np.float64,axis=None)
    #stdtempcold = np.std(tempcold, axis=None, dtype=np.float64,ddof=1)
    #meantemphot = np.mean(temphot,dtype=np.float64,axis=None)
    #stdtemphot = np.std(temphot, axis=None, dtype=np.float64,ddof=1)

    #KRcold = meanJ/meantempcold
    #KRcoldstd = KRcold*np.sqrt((stdJ/meanJ)**2+(stdtempcold/meantempcold)**2)

    #KRhot = meanJ/meantemphot
    #KRhotstd = KRhot*np.sqrt((stdJ/meanJ)**2+(stdtemphot/meantemphot)**2)

    #print 'cold', KRcold, '+/-', KRcoldstd, 'hot', KRhot, '+/-', KRhotstd
    print meanT, stdT, meanJO, stdJO, inverseJO, inverseJOstd, KRO, KROstd, meanJI, stdJO, inverseJI, inverseJIstd, KRI, KRIstd

    #----------------Combine two measurements to get composite kappa value--------------
    ## Source: http://www.burtonsys.com/climate/composite_standard_deviations.html

    KR = (KRO+KRI)/2
    sum_of_error_of_sum_squares = (KROstd**2)*(numblock-1)+(KRIstd**2)*(numblock-1)
    overall_group_sum_of_squares = ((KR-KRO)**2)*numblock + ((KR-KRI)**2)*numblock
    KRstd = np.sqrt((sum_of_error_of_sum_squares+overall_group_sum_of_squares)/(2*numblock-1))

    print 'KR: ', KR, '+/-',KRstd, 'in MW/m^2/K' #Check for proper unit conversion

    with open(str('deltaT_std_JO_std_invJO_std_KRO_std_JI_std_invJI_std_KRI_std_KR_std'),'w+') as f:
	    print >> f, meanT, stdT, meanJO, stdJO, inverseJO, inverseJOstd, KRO, KROstd, meanJI, stdJI, inverseJI, inverseJIstd, KRI, KRIstd, KR, KRstd
	    f.close()

    date = datetime.datetime.today()


    kappaDataInterface = open("kappaDataInterface","a+")

    kappaDataInterface.write(str(KR) + " +/- " + str(KRstd) + "   begin timestep: " +str(begintime) +"	 "+ "end timestep: "+ str(endtime) + "	 date: " + str(date)+'\n')

    kappaDataInterface.close();


