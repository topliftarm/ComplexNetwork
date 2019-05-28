import datetime
import os, shutil
import ode_solver
import numpy as np
from time  import time
from network import make_graph
from modules import *
from numpy import pi
import matplotlib.pyplot as plt
import pandas as pd

#DirPath = "/storage/users/fbaharifard/ComplexNetworks/initData"
DirPath = "~/Documents/Complex network/c/src"

NumberOfNodes = 100
couplingStrength = 0.27
NumberOfIterations = 1000
NumbertOfSteps = 400
tfinal = 100.0
tinitial = 0.0
dt = 0.1
times = np.arange(0,tfinal, dt)

graph  = make_graph()

#motherDirPath = "/storage/users/fbaharifard/ComplexNetworks/CompleteData"
motherDirPath = "~/Documents/Complex network/c/CompleteData"


start = time()
runNumber = 1
print("runNumber = ", str(runNumber))
dirName = "Run"+str(runNumber)
os.mkdir(motherDirPath+"/"+dirName)
currentPath = motherDirPath+"/"+dirName

fileName_adj = 'FinalAdj.txt'
Adj = np.loadtxt(DirPath+'/'+fileName_adj).reshape(-1,NumberOfNodes*NumberOfNodes).ravel().tolist()

fileName_omega = 'Omega.txt'
Omega = pd.read_csv(DirPath+'/'+fileName_omega, header=None).stack().tolist()

fileName_Y = 'finalY.txt'
InitialCondition = pd.read_csv(DirPath+'/'+fileName_Y, header=None).stack().tolist()

print("start simulation...")

obj = ode_solver.ODE(NumberOfNodes, tfinal, dt,	couplingStrength, InitialCondition, Omega, NumbertOfSteps, NumberOfIterations)
#sol = obj.integrate(Adj, rewire=True, selfish=True)
sol = obj.integrate(Adj, rewire=True, selfish=True, NumberOfSelfishNodes=0)


sol = np.asarray(sol)
r_glob, psi = obj.get_order_parameters()
#print('r_glob = ',len(r_glob))
MeanRinEachIteration = obj.getMeanRinEachIteration()
acceptanceRateRewiring = obj.getAcceptanceRewiring()
#MeanYPrime = obj.getMeanYPrime();

finalAdj = obj.getCij()
finalY = obj.getFinalY()
final_adj_mat = np.asarray(finalAdj).reshape((NumberOfNodes, NumberOfNodes))
sumKin = np.array([sum(final_adj_mat.T[i]) for i in range(NumberOfNodes)])
numberOfBins = 20
_sumOfHist, _bins = np.histogram(Omega, bins=numberOfBins)
digitized = np.digitize(Omega, _bins)
sumKin_bin_means = [sumKin[digitized == i].mean() for i in range(1, len(_bins)+1)]
sliceOfOmega = _bins
print("simulation...done")
del obj, sol
# ------------------cal time----------------------- #
display_time(time()-start)
# ------------------------------------------------- #
currentTime = datetime.datetime.now()
hour = currentTime.hour;
runNumber += 1
#-------------------- save --------------------
print("saving Results ...")
#if (os.path.isfile("/home/fatemeh/Documents/ComplexNetwork/src/output.txt")):
#    shutil.move("/home/fatemeh/Documents/ComplexNetwork/src/output.txt",
#            currentPath+"/"+"output.txt")

with open(currentPath+'/'+'FinalAdj.txt', 'w+') as FinalAdjFile:
    np.savetxt(FinalAdjFile, final_adj_mat)

with open(currentPath+'/'+'acceptanceRateRewiring.txt', 'w+') as acceptanceRateRewiringFile:
    np.savetxt(acceptanceRateRewiringFile, np.array(acceptanceRateRewiring).T)

with open(currentPath+'/'+'MeanRinEachIteration.txt', 'w+') as MeanRinEachIterationFile:
    np.savetxt(MeanRinEachIterationFile, np.array(MeanRinEachIteration).T)

with open(currentPath+'/'+'finalY.txt', 'w+') as finalYFile:
    np.savetxt(finalYFile, np.array(finalY).T)

with open(currentPath+'/'+'r_glob.txt', 'w+') as r_globFile:
    np.savetxt(r_globFile, np.array(r_glob).T)

with open(currentPath+'/'+'sumKin.txt', 'w+') as sumKinFile:
    np.savetxt(sumKinFile, np.array(sumKin).T)

with open(currentPath+'/'+'sumKin_bin_means.txt', 'w+') as sumKin_bin_meansFile:
    np.savetxt(sumKin_bin_meansFile, np.array(sumKin_bin_means).T)

with open(currentPath+'/'+'sliceOfOmega.txt', 'w+') as sliceOfOmegaFile:
    np.savetxt(sliceOfOmegaFile, np.array(sliceOfOmega).T)

print("saving Results ... done")
#-------------------- save --------------------

#os.system("mpg123 "+" finish.mp3")
