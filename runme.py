import datetime
import os, shutil
import ode_solver

import numpy as np
from time  import time
from network import make_graph
from modules import *
from numpy import pi
import matplotlib.pyplot as plt

print("Make clean...")
os.system("make clean ")
print("Make...")
os.system("make")

print("configuration...\n")
seed = 12358
np.random.seed(seed)

NumberOfNodes = 100
couplingStrength = 0.27
NumberOfIterations = 500
NumbertOfSteps = 200
tfinal = 100.0
tinitial = 0.0
dt = 0.1

times = np.arange(0,tfinal, dt)

graph  = make_graph()

currentTime = datetime.datetime.now()
hour = currentTime.hour;

#motherDirPath = "/storage/users/fbaharifard/ComplexNetworks/CompleteData"
motherDirPath = "/home/vahid/Documents/Complex network/c/CompleteData"
runNumber = 1
start = time()


#Adj = graph.barabasi(NumberOfNodes, m=2)
#Adj = graph.erdos_renyi_graph(NumberOfNodes, 0.04)
Adj = graph.random_k_out_graph(NumberOfNodes, 6, seed=np.random.randint(100000))
adj_mat = np.asarray(Adj).reshape((NumberOfNodes, NumberOfNodes))
InitialSumKin = np.array([sum(adj_mat.T[i]) for i in range(NumberOfNodes)])
#degree = np.sum(adj_mat, axis=1)

# Omega = np.random.normal(0, 1, NumberOfNodes).tolist()
Omega = np.random.uniform(-1, 1, size=NumberOfNodes).tolist()
#InitialCondition = np.random.uniform(0, 2*pi, NumberOfNodes).tolist()
InitialCondition = np.random.normal(0, 2*pi, NumberOfNodes).tolist()


while(runNumber < 3): #hour < 14):

    print("runNumber = ", str(runNumber))
    dirName = "Run"+str(runNumber)
    os.mkdir(motherDirPath+"/"+dirName)
    currentPath = motherDirPath+"/"+dirName

    #-------------------- save --------------------
    print("saving Initialization data ...")

    with open(currentPath+'/'+'Config.txt', 'w+') as ConfigFile:
        ConfigFile.write("NumberOfNodes= %d\r\n" % NumberOfNodes)
        ConfigFile.write("couplingStrength= %f\r\n" % couplingStrength)
        ConfigFile.write("NumberOfIterations= %d\r\n" % NumberOfIterations)
        ConfigFile.write("NumbertOfSteps= %d\r\n" % NumbertOfSteps)

    with open(currentPath+'/'+'InitAdj.txt', 'w+') as InitAdjFile:
        np.savetxt(InitAdjFile, adj_mat)

    with open(currentPath+'/'+'Omega.txt', 'w+') as OmegaFile:
        np.savetxt(OmegaFile, np.array(Omega).T)

    with open(currentPath+'/'+'Y0.txt', 'w+') as Y0File:
        np.savetxt(Y0File, np.array(InitialCondition).T)

    with open(currentPath+'/'+'InitialSumKin.txt', 'w+') as InitialSumKinFile:
        np.savetxt(InitialSumKinFile, np.array(InitialSumKin).T)


    print("saving Initialization data ... done")
    #-------------------- save --------------------

    print("start simulation...")
    #obj = ode_solver.ODE(NumberOfNodes, tfinal, dt, couplingStrength, InitialCondition, Omega, NumbertOfSteps, NumberOfIterations)
    #sol = obj.integrate(Adj, rewire=True, selfish=True)
    if (runNumber == 1):
        obj = ode_solver.ODE(NumberOfNodes, tfinal, dt,	couplingStrength, InitialCondition, Omega, NumbertOfSteps, NumberOfIterations)
        sol = obj.integrate(Adj, rewire=True, selfish=True, NumberOfSelfishNodes=0)
    else:
        obj = ode_solver.ODE(NumberOfNodes, tfinal, dt,	couplingStrength, InitialCondition, Omega, NumbertOfSteps, NumberOfIterations)
        sol = obj.integrate(Adj, rewire=True, selfish=False, NumberOfSelfishNodes=2)
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
    #if (os.path.isfile("/home/vahid/Documents/Complex network/c/src/output.txt")):
    #    shutil.move("/home/vahid/Documents/Complex network/c/src/output.txt",
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




# -----------------------  Plotes -----------------
#sizeFig = (9,7)
#_dpi = 80
#plt.figure(figsize=sizeFig, dpi=_dpi)
#plt.plot(range(1, 1+len(acceptanceRateRewiring)), acceptanceRateRewiring)
#plt.title('acceptanceRateRewiring_global')
#plt.savefig('acceptanceRateRewiring_global.png')

#plt.figure(figsize=sizeFig, dpi=_dpi)
#plt.plot(MeanYPrime)
#plt.title('MeanYPrime_selfish')
#plt.savefig('MeanYPrime_selfish.png')

#plt.figure(figsize=sizeFig, dpi=_dpi)
#plt.plot(r_glob)

#plt.title('r_glob_global')
#plt.savefig('r_glob_global.png')

#plt.figure(figsize=sizeFig, dpi=_dpi)
#plt.plot(MeanRinEachIteration)
#plt.title('MeanRinEachIteration_global')
#plt.savefig('MeanRinEachIteration_global.png')

#plt.figure(figsize=sizeFig, dpi=_dpi)
#plt.plot(sliceOfOmega, sumKin_bin_means, '*')
#plt.title('omega_kin')
#plt.savefig('omega_kin.png')

#plt.figure(figsize=sizeFig, dpi=_dpi)
#plt.plot(sumKin, '*')
#plt.title('final kin of all nodes')
#plt.savefig('finalKinOfAllNodes.png')

#plt.figure(figsize=sizeFig, dpi=_dpi)
#plt.plot(InitialSumKin, '*')
#plt.title('Initial kin of all nodes')
#plt.savefig('InitialKinOfAllNodes.png')

#plt.show()
# ----------------------------------------- #
#fig,axs = pl.subplots(1,figsize=(8,6))
#axs.plot(g, R,   marker=">", markersize=3, lw=2, c='r', label=r'$R_{global}$')
#axs.set_ylabel("R")
#axs.legend()
#axs.set_xlabel(r"$K$",  fontsize=13)
#pl.savefig(adrsf+"R.png")
# pl.close()
# pl.show()

#os.system("mpg123 "+" finish.mp3")
