#!/usr/bin/env python3

from select import select
import struct
import ase
from ase import io
import numpy as np
import logging
import multiprocessing as mp 
import os



class ShadowMD():
    def __init__(self,
        trajectory=None,
        parallel=False,
        startPoint=0,
        endPoint=None):
        
        print("Loading trajectory")
        try: 
            self.trajectory=io.read(trajectory, ":")
        except:
            self.logfile.error("ERROR: Could not load trajectory file with name: ", trajectory)
            self.logfile.error("Stopping")
            exit()
        print("Loaded trajectory file succesfully.")
        #prepare the trajectory
        #self.trajectory=self.trajectory[::selectEvery]

        #set the default startpoint and endpoint - these can be changed manually
        if startPoint is None:
            self.startPoint=0
        else:
            self.startPoint=startPoint
        if endPoint is None:
            self.endPoint=len(self.trajectory)
        else:
            self.endPoint=endPoint
        print("Starting at:", self.startPoint, "Ending at:", endPoint, " For a total of", endPoint-self.startPoint, "images")
        self.trajectory=self.trajectory[self.startPoint:endPoint]

        
        print("Trajectory has "+str(len(self.trajectory))+" images.")

    '''
        self.parallel=parallel
        if self.parallel:
            print("Running in parallel")
            self.cpus=int(os.environ['SLURM_NTASKS'])
            print("Number of CPUs: ", self.cpus)
            self.return_dict = mp.Manager().dict()
            self.blockCount=1
        else:
            print("Running in series")

    '''
    
    def moveAdsorbate(self, structure):
        #find the carbon in formate
        atomsInAdsorbate=[i.index for i in structure if i.symbol=="C" or i.symbol=="O" or i.symbol=="H"] #note, this takes common organics, but needs to be adjusted if you have adsorbates with other elements
        translateArray=np.zeros((len(structure),3))
        for atomNum in atomsInAdsorbate:
            translateArray[atomNum,2]=7
        structure.translate(translateArray)

    def removeAdsorbate(self, structure):
        atomsInAdsorbate=[i.index for i in structure if i.symbol=="C" or i.symbol=="O" or i.symbol=="H"] #note, this takes common organics, but needs to be adjusted if you have adsorbates with other elements
        del structure[[atom.index for atom in structure if atom.index in atomsInAdsorbate]]

    '''
    def parallelEvalE(self, structures, block_num):
        stepEnergies=None
        computedTrajectory=ase.io.Trajectory("OutputTrajectory"+str(block_num)+".traj", mode="a")
        for i, structure in enumerate(structures):
            structure.calc=self.calc
            self.removeFormate(structure)
            print("Running vasp calc on image "+ str(i) +" in block "+str(block_num))
            energy=structure.get_potential_energy()
            print("Done")
            computedTrajectory.write(structure)
            if i==0:
                stepEnergies=np.array([[i, energy, energy+formateBoxE]])
            else:
                stepEnergies=np.append(stepEnergies,[[i, energy, energy+formateBoxE]], axis=0)
            np.savetxt("StepEnergy"+str(block_num)+".txt",stepEnergies)
            self.return_dict[block_num]=stepEnergies.copy()
        computedTrajectory.close()
    '''

    def run(self):
        #takes a SINGLE structure and performs a single-point energy evaluation
        for i, structure in enumerate(self.trajectory[:]):
            self.removeAdsorbate(structure)
            name="image"+str(i+self.startPoint)
            os.system("cp -r template "+name)
            os.chdir(name)
            io.write("POSCAR", structure)
            os.chdir("../")




if __name__ == '__main__':
    import sys
    try: 
        loadfile=sys.argv[1]
    except:
        print("No file to load given as input argument")
        exit()
    #check if sys.argv[2] and sys.argv[3] exist, and if not set them to None
    try:
        sys.argv[2]
    except:
        startPoint=None
    try:
        sys.argv[3]
    except:
        endPoint=None
    #check if sys.argv[2] and sys.argv[3] are integers, and confrim that sys.argv[2] is less than sys.argv[3]
    if startPoint is not None and endPoint is not None:
        try:
            startPoint=int(startPoint)
            endPoint=int(endPoint)
        except:
            print("Start and end points must be integers")
            exit()
        if startPoint>endPoint:
            print("Start point must be less than end point")
            exit()
    shadowMD=ShadowMD(trajectory=loadfile, startPoint=startPoint, endPoint=endPoint)
    shadowMD.run()



        
        
        


