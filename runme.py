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

def saveData(saveList):
    print("saving data ...")
    for _FileName, _Content in saveList.items():
        with open(currentPath+'/'+_FileName, 'w+') as File:
            if (np.array(_Content).ndim > 0):
                np.savetxt(File, np.array(_Content))
            else:
                File.write(str(_Content))
    print("saving data ... done")

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
    initialList = {
        'NumberOfNodes.txt':NumberOfNodes,
        'couplingStrength.txt':couplingStrength,
        'NumberOfIterations.txt':NumberOfIterations,
        'NumbertOfSteps.txt':NumbertOfSteps,
        'InitAdj.txt':adj_mat,
        'Omega.txt':np.array(Omega).T,
        'Y0.txt':np.array(InitialCondition).T,
        'InitialSumKin.txt':np.array(InitialSumKin).T
    }
    saveData(initialList)
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
    finalList = {
        'FinalAdj.txt':final_adj_mat,
        'acceptanceRateRewiring.txt':np.array(acceptanceRateRewiring).T,
        'MeanRinEachIteration.txt':np.array(MeanRinEachIteration).T,
        'finalY.txt':np.array(finalY).T,
        'r_glob.txt':np.array(r_glob).T,
        'sumKin.txt':np.array(sumKin).T,
        'sumKin_bin_means.txt':np.array(sumKin_bin_means).T,
        'sliceOfOmega.txt':np.array(sliceOfOmega).T
    }
    saveData(finalList)
    display_time(time()-start)
    currentTime = datetime.datetime.now()
    hour = currentTime.hour;
    runNumber += 1

os.system("mpg123 "+" finish.mp3")
