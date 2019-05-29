import datetime
import os, shutil
import ode_solver
import numpy as np
from time  import time
from network import make_graph
from modules import *
from numpy import pi
import matplotlib.pyplot as plt

# ===================================== functions ========================
def make():
    print("Make clean...")
    os.system("make clean ")
    print("Make...")
    os.system("make")

def configGraph():
    #Adj = graph.barabasi(NumberOfNodes, m=2)
    #Adj = graph.erdos_renyi_graph(NumberOfNodes, 0.04)
    Adj = graph.random_k_out_graph(NumberOfNodes, NumberOfEdges, seed=np.random.randint(100000))
    adj_mat = np.asarray(Adj).reshape((NumberOfNodes, NumberOfNodes))
    InitialSumKin = np.array([sum(adj_mat.T[i]) for i in range(NumberOfNodes)])
    #degree = np.sum(adj_mat, axis=1)
    Omega = np.random.uniform(-1, 1, size=NumberOfNodes).tolist()
    InitialCondition = np.random.normal(0, 2*pi, NumberOfNodes).tolist()

    return Adj, adj_mat, InitialSumKin, Omega, InitialCondition;

def configParameters():
    print("configuration...\n")
    seed = 12358
    np.random.seed(seed)
    NumberOfNodes = 100
    NumberOfEdges = 6
    couplingStrength = 0.27
    NumberOfIterations = 500
    NumbertOfSteps = 200
    tfinal = 100.0
    tinitial = 0.0
    dt = 0.1
    times = np.arange(0,tfinal, dt)
    graph  = make_graph()
    numberOfBins = 20
    return (NumberOfNodes, NumberOfEdges, couplingStrength, NumberOfIterations,
            NumbertOfSteps, tfinal, tinitial, dt, times, graph, numberOfBins)

def createDir(runNumber, motherDirPath):
    print("runNumber = ", str(runNumber))
    dirName = "Run"+str(runNumber)
    os.mkdir(motherDirPath+"/"+dirName)
    currentPath = motherDirPath+"/"+dirName
    return currentPath

def saveInitialData(currentPath):
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

def getParametersFromC():
    r_glob, psi = obj.get_order_parameters()
    MeanRinEachIteration = obj.getMeanRinEachIteration()
    acceptanceRateRewiring = obj.getAcceptanceRewiring()
    #MeanYPrime = obj.getMeanYPrime();
    finalAdj = obj.getCij()
    finalY = obj.getFinalY()
    final_adj_mat = np.asarray(finalAdj).reshape((NumberOfNodes, NumberOfNodes))
    sumKin = np.array([sum(final_adj_mat.T[i]) for i in range(NumberOfNodes)])
    _sumOfHist, _bins = np.histogram(Omega, bins=numberOfBins)
    digitized = np.digitize(Omega, _bins)
    sumKin_bin_means = [sumKin[digitized == i].mean() for i in range(1, len(_bins)+1)]
    sliceOfOmega = _bins
    return (r_glob, psi, MeanRinEachIteration, acceptanceRateRewiring, finalAdj,
            finalY, final_adj_mat, sumKin, sumKin_bin_means, sliceOfOmega)

def saveFinalData(currentPath):
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

# ===================================== --------- ========================


make()
NumberOfNodes, NumberOfEdges, couplingStrength, \
NumberOfIterations, NumbertOfSteps, tfinal, tinitial, \
dt, times, graph, numberOfBins = configParameters()

Adj, adj_mat, InitialSumKin, \
Omega, InitialCondition = configGraph()

currentTime = datetime.datetime.now()
hour = currentTime.hour;
start = time()
#motherDirPath = "/storage/users/fbaharifard/ComplexNetworks/CompleteData"
motherDirPath = "/home/vahid/Documents/Complex network/c/CompleteData"
runNumber = 1
while(runNumber < 3): #hour < 14):
    currentPath = createDir(runNumber, motherDirPath)
    saveInitialData(currentPath)
    print("start simulation...")
    if (runNumber == 1):
        obj = ode_solver.ODE(NumberOfNodes, tfinal, dt,	couplingStrength, \
                   InitialCondition, Omega, NumbertOfSteps, NumberOfIterations)
        sol = obj.integrate(Adj, rewire=True, selfish=True, NumberOfSelfishNodes=0)
    else:
        obj = ode_solver.ODE(NumberOfNodes, tfinal, dt,	couplingStrength, \
                   InitialCondition, Omega, NumbertOfSteps, NumberOfIterations)
        sol = obj.integrate(Adj, rewire=True, selfish=False, NumberOfSelfishNodes=2)
    sol = np.asarray(sol)
    r_glob, psi, MeanRinEachIteration, acceptanceRateRewiring, \
    finalAdj,finalY, final_adj_mat, sumKin, sumKin_bin_means, \
    sliceOfOmega = getParametersFromC()
    print("simulation...done")
    del obj, sol
    saveFinalData(currentPath)
    display_time(time()-start)
    currentTime = datetime.datetime.now()
    hour = currentTime.hour;
    runNumber += 1



os.system("mpg123 "+" finish.mp3")
