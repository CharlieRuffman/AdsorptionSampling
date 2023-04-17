#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys

traspositionFactor=0
#traspositionFactor=-30.205455 #Methanol
#traspositionFactor=-21.97165359+1/2*-6.4040815 #Formate
#-30.205455 - -3.38362605*2 # methanol without H2


#steps to omit - first 10 ps 
omitSteps=125

#print a wanrning if the transposition factor is set to anything but zero
if traspositionFactor!=0:
    print("WARNING: transposition factor is set: ",traspositionFactor)


#ensure there are two arguments
if len(sys.argv)!=3:
    print("Usage: python3 EnergyHistogram.py <energy file> <Type, one of: clean/reference/adsorbed/adsenergy>")
    exit()

#read in an input file and make a histogram plot out of the energies
file=sys.argv[1]
type=sys.argv[2].lower()
if type not in ["clean","reference","adsorbed","adsenergy"]:
    print("Type must be clean, reference, adsorbed, or adsenergy")
    exit()

data=np.loadtxt(file)[omitSteps:,-1]

print("Number of data points: ",len(data))

#transpose the data if the transposition factor is set
if traspositionFactor!=0:
    data=data-traspositionFactor

#calcuate the mean and standard deviation of the data
mean=np.mean(data)
std=np.std(data)
#print these
print("Mean: ",mean)
print("Standard Deviation: ",std)

fig, ax1 = plt.subplots()
ax1.set_ylabel("Count")
if type=="clean":
    ax1.hist(data, bins=50, color="black")
    ax1.set_xlabel("Energy / eV")
elif type=="reference":
    ax1.hist(data, bins=50, color="dimgrey")
elif type=="adsorbed":
    ax1.hist(data, bins=50, color="darkgrey")
    ax1.set_xlabel("Energy / eV")
elif type=="adsenergy":
    ax1.hist(data, bins=50, color="lightgrey")
    ax1.set_xlabel("Adsorption Energy / eV")
h, bins = np.histogram(data, bins=50)

#set the x-axis limits to a consistent range - this should be adjusted by the user 
ax1.set_xlim(-93,-88.5)

#save the data to textfiles
np.savetxt(type+"rawHistData.txt",np.column_stack((h,bins[:-1])),fmt="%i %f")

 #save the figure
plt.savefig(type+"_hist.svg",dpi=1200)



