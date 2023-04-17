#!/usr/bin/env python3

import math
from sre_constants import IN
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import os
from ase import io

'''
USER SET PARAMETERS
'''

#add a translation factor variable for the reference state - This is to adjust the reference state energy to the same average as the clean MD simulation 
traspositionFactor=0.212195576


#detection window settings for low and high energy windows
detectionWindowSize=1 #minimum size of high/low window in ps
proportionExceeding=0.6 #proportion of the window that must be above/below the threshold
sdevProportion=0.6 #standard deveiation threshold for high/low windows

equilTime=10 #time to be excluded from analysis for equilibration the start of the MD run in ps
selectEvery=80 #Resolution of sampling in timesteps - important for getting the timescales correct

timestep=1 #timestep of the AIMD simulation in fs

'''
MAIN script
'''


#add a warning when the translation factor is turned on
if traspositionFactor!=0:
    print("WARNING: transposition factor is set: ",traspositionFactor)


matplotlib.use("Agg")

#Load in the data from the command line arguments
mainData=sys.argv[1] #this is the adsorption energy data from the main MD simulation
shadowData=sys.argv[2] #this is the adsorption energy data from the shadow MD simulation - i.e. the reference state

#convert the timescales to steps
equilSteps=math.ceil(equilTime*1000/(selectEvery*timestep))
detectionWindowSteps=math.ceil(detectionWindowSize*1000/selectEvery)
windowThreshold=math.ceil(proportionExceeding*detectionWindowSteps)

#try to load the datafile in the first argument 
try:
    mainData = np.loadtxt(mainData)
except:
    print("Could not load main MD data from first argument - stopping")
    exit()

#try to load the second datafile
try:
    shadowData = np.loadtxt(shadowData)
except:
    print("Could not shadow MD data from second argument - stopping")
    exit()

#translate the reference state by the transposition factor (if set)
if traspositionFactor!=0:
    shadowData[:,2]=shadowData[:,2]-traspositionFactor

#check two arrays are same length and if not, shorten one:
if len(mainData[:,0]) != len(shadowData[:,0]):
    print("Two input arrays are not the same length - exiting")
    exit()

#set the font size to 16
plt.rc('font', size=16)
width=20

#Calculate the difference between the two data sets
differenceData=(mainData[:,1]-shadowData[:,2])

"""
Processing for calculating the maxima and minima of the adsorption energy. These quantities are sometimes relevant, but not always
used. 
"""
steps=mainData[:,0]
diffsteps=steps.copy()
differenceDataMod=differenceData.copy()
differenceDataMod, diffsteps =(list(t) for t in zip(*sorted(zip(differenceDataMod[equilSteps:], diffsteps[equilSteps:]))))

k=10 #number of lowest and highest values to print
maximumVals=diffsteps[len(diffsteps)-k:]
minimumVals=diffsteps[:k]

def checkListUpper(data, upperThreshold,windowThreshold, count):
    upperOverages = [x for x in data if (x > upperThreshold)]
    if len(upperOverages) >= windowThreshold:
        return np.arange(count, count+len(upperOverages))
    else: 
        return None

def checkListLower(data, lowerThreshold,windowThreshold,count):
    lowerOverages = [x for x in data if (x < lowerThreshold)]
    if len(lowerOverages) >= windowThreshold:
        return np.arange(count, count+len(lowerOverages))
    else: 
        return None

detectHighLow=True
if detectHighLow:
    differenceDataEquilibrated=differenceData[equilSteps:]
    plt.hist(differenceDataEquilibrated, bins="auto",color="grey")
    standardDev=np.std(differenceDataEquilibrated)
    mean=np.mean(differenceDataEquilibrated)
    print("Average adsorption energy (after equilibration): ", mean, "eV")
    highThreshold=mean+sdevProportion*standardDev
    lowThreshold=mean-sdevProportion*standardDev



    plt.xlabel("Adsorption energy / eV")
    plt.ylabel("Frequency")
    plt.tight_layout()  
    plt.savefig("hist.png")
    plt.clf()

    highWindows=[]
    lowWindows=[]

    count=0
    for i in range(len(differenceDataEquilibrated)):
        upperIndices=checkListUpper(differenceDataEquilibrated[i:i+detectionWindowSteps],highThreshold,windowThreshold, i)
        lowerIndices=checkListLower(differenceDataEquilibrated[i:i+detectionWindowSteps],lowThreshold,windowThreshold, i)
        if upperIndices is not None:
            highWindows.extend(upperIndices)
        if lowerIndices is not None:
            lowWindows.extend(lowerIndices)


    highWindows=np.unique(highWindows)+equilSteps
    lowWindows=np.unique(lowWindows)+equilSteps

    #separate the low windows where their values are non-consequtive
    separatedLow=np.split(lowWindows, np.where(np.diff(lowWindows) != 1)[0]+1)
    #do the same for the high windows
    separatedHigh=np.split(highWindows, np.where(np.diff(highWindows) != 1)[0]+1)
    
    #find the lowest energy low window
    prevMean=999999
    for window in separatedLow:
        meanEnergy=np.mean([differenceDataEquilibrated[j-equilSteps] for j in window])
        if meanEnergy<prevMean:
            lowWindows=window.copy()
            prevMean=meanEnergy

    #do the same for high
    prevMean=-999999
    for window in separatedHigh:
        meanEnergy=np.mean([differenceDataEquilibrated[j-equilSteps] for j in window])
        if meanEnergy>prevMean:
            highWindows=window.copy()
            prevMean=meanEnergy


print("High: ", highWindows)

print("Low: ", lowWindows)



printAdsE=True


a=io.read("SelectedImg.traj", ":")


if len(highWindows) != 0:
    highAv=np.mean([differenceDataEquilibrated[i-equilSteps] for i in highWindows])
    print("High avereage:", highAv)
else:
    print("No high windows detected with current settings - adjust and try again")

if len(lowWindows) != 0:
    lowAv=np.mean([differenceDataEquilibrated[i-equilSteps] for i in lowWindows])
    print("low avereage:", lowAv)
else:
    print("No low windows detected with current settings - adjust and try again")

'''
dirname="HighRegionStructs"
try:
    os.mkdir(dirname)
except:
    print("High directory exits: ", dirname)
os.chdir(dirname)

for i, j in enumerate(highWindows):
    io.write("highRegionStruct"+str(j)+".POSCAR",a[int(j)])

os.chdir("../")
dirname="LowRegionStructs"
try:
    os.mkdir(dirname)
except:
    print("Low directory exits: ", dirname)
os.chdir(dirname)

for i, j in enumerate(lowWindows):
    io.write("lowRegionStruct"+str(j)+".POSCAR",a[int(j)])
os.chdir("../")



print(k, "steps with the lowest energy (lowest first): ", minimumVals)
if printAdsE:
    print(differenceDataMod[:k])
    for i, j in enumerate(minimumVals):
        io.write(str(i)+"lowestE.POSCAR",a[int(j)])
print(k, "steps with the highest energy (highest last): ", maximumVals)
if printAdsE:
    print(differenceDataMod[len(diffsteps)-k:])
    for i,j in enumerate(maximumVals):
        io.write(str(i)+"highestE.POSCAR",a[int(j)])
'''



time=np.arange(0,len(differenceData))*selectEvery*timestep/(1000)
plt.plot(time[:],differenceData[:],color="black", linewidth=1.5)

plt.axhline(mean, color="darkred", xmin=equilTime/time[-1], linewidth=1.5, linestyle="--")
if detectHighLow:
    for val in highWindows:
        plt.axvline(x = val*selectEvery/1000, color = 'r', label = 'High energy', alpha=0.05)

    for val in lowWindows:
        plt.axvline(x = val*selectEvery/1000, color = 'b', label = 'Low energy', alpha=0.05)



#y axis between 0 and 10
plt.xlabel("Time / ps")

plt.ylabel("Adsorption energy / $eV$")
plt.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.show()
plt.savefig("adsorptionEnergy.svg",dpi=1200)

print("Done")