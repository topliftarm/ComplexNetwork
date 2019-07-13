import datetime
import os, shutil
import ode_solver
import numpy as np
from time  import time
from network import make_graph
from modules import *
from numpy import pi
import matplotlib.pyplot as plt
import math

# ===================================== functions ========================
def make():
    print("Make clean...")
    os.system("make clean ")
    print("Make...")
    os.system("make")

def deleteParallelEdges(_adj_mat, NumberOfNodes):
    for nodeNumber in range(NumberOfNodes):
        index = np.where(_adj_mat[nodeNumber]>1)
        index = np.concatenate(index).tolist()
        #print('node=',nodeNumber, 'index=',index)
        for element in index:
            diff = int(_adj_mat[nodeNumber][element]-1)
            _adj_mat[nodeNumber][element] = 1
            rand_element_with_self = np.random.permutation(NumberOfNodes)
            rand_element_without_self = np.delete(rand_element_with_self.tolist(), \
                            np.where(rand_element_with_self == [nodeNumber]))
            l = 0
            while(diff>0):
                if(_adj_mat[nodeNumber][rand_element_without_self[l]] == 0):
                    _adj_mat[nodeNumber][rand_element_without_self[l]] = 1
                    l = 0
                    diff = diff-1
                else:
                    l = l+1
    #print('np.where(_adj_mat>1))=',np.where(_adj_mat>1))
    return _adj_mat

def deleteBiEdges(adj_mat_with_Bi_edges, NumberOfNodes, NumberOfEdges):
    for i in range(NumberOfNodes):
        for j in range(i+1,NumberOfNodes):
            if (adj_mat_with_Bi_edges[i][j]):
                if(adj_mat_with_Bi_edges[i][j]==adj_mat_with_Bi_edges[j][i]):
                    adj_mat_with_Bi_edges[j][i] = 0
    return adj_mat_with_Bi_edges

def FixDegree(adj_mat_LowDegree, NumberOfNodes, NumberOfEdges):
    listNodes = list()
    for i in range(NumberOfNodes):
        if (sum(adj_mat_LowDegree[i])<NumberOfEdges):
            [listNodes.append(i) for _i in range(int(NumberOfEdges-sum(adj_mat_LowDegree[i])))]

    for i in listNodes:
        index = np.where(adj_mat_LowDegree[i]<1)
        index = np.concatenate(index).tolist()
        if i in index: index.remove(i)
        for ii in index:
            if(adj_mat_LowDegree[ii][i]):index.remove(ii)
        index = np.random.permutation(index)
        adj_mat_LowDegree[i][index[0]] = 1

    return adj_mat_LowDegree

def addBiEdges(adj_mat_with_Bi_edges, NumberOfNodes, NumberOfBiEdges):
    for i in range(NumberOfBiEdges):
        index_ = np.random.randint(1, NumberOfNodes, size=2)
        x, y = index_[0], index_[1]
        if(adj_mat_with_Bi_edges[x][y] == 0):
           if(adj_mat_with_Bi_edges[y][x] == 0):
              index = np.where(adj_mat_with_Bi_edges[x]>0)
              index = np.concatenate(index).tolist()
              index = np.random.permutation(index)
              adj_mat_with_Bi_edges[x][index[0]] = 0
              index = np.where(adj_mat_with_Bi_edges[y]>0)
              index = np.concatenate(index).tolist()
              index = np.random.permutation(index)
              adj_mat_with_Bi_edges[y][index[0]] = 0

              adj_mat_with_Bi_edges[x][y] = 1
              adj_mat_with_Bi_edges[y][x] = 1
    return adj_mat_with_Bi_edges

def configGraph():
    #Adj = graph.barabasi(NumberOfNodes, m=2)
    #Adj = graph.erdos_renyi_graph(NumberOfNodes, 0.04)
    _Adj = graph.random_k_out_graph(NumberOfNodes, NumberOfEdges, seed=np.random.randint(100000))
    adj_mat_with_parallel_edges = np.asarray(_Adj).reshape((NumberOfNodes, NumberOfNodes))
    adj_mat_with_Bi_edges = deleteParallelEdges(adj_mat_with_parallel_edges, NumberOfNodes)
    # adj_mat_with_Bi_edges = addBiEdges(adj_mat_with_Bi_edges, NumberOfNodes, 50)
    #adj_mat_LowDegree = deleteBiEdges(adj_mat_with_Bi_edges, NumberOfNodes, NumberOfEdges)
    #adj_mat = FixDegree(adj_mat_LowDegree, NumberOfNodes, NumberOfEdges)
    adj_mat = FixDegree(adj_mat_with_Bi_edges, NumberOfNodes, NumberOfEdges)
    Adj = adj_mat.reshape(-1)
    InitialSumKin = np.array([sum(adj_mat.T[i]) for i in range(NumberOfNodes)])
    #degree = np.sum(adj_mat, axis=1)
    #l = np.linspace(-1,+1,21)
    #Omega = l.repeat(5,0)
    #np.random.shuffle(Omega)
    #Omega = np.random.uniform(-1.6, 1.6, size=NumberOfNodes).tolist()
    Omega = np.random.normal(0, 0, size=NumberOfNodes).tolist()
    InitialCondition = np.random.normal(0, pi/2.0, NumberOfNodes).tolist()

    return Adj, adj_mat, InitialSumKin, Omega, InitialCondition;

def configParameters():
    print("configuration...\n")
    seed = 56872
    np.random.seed(seed)
    NumberOfNodes = 100
    NumberOfEdges = 6
    couplingStrength = 0.3
    NumberOfIterations = 1000
    NumbertOfSteps = 1000
    rewire = True
    tfinal = 100.0
    tinitial = 0.0
    dt = 0.1
    times = np.arange(0,tfinal, dt)
    graph  = make_graph()
    numberOfBins = 20
    return (rewire, NumberOfNodes, NumberOfEdges, couplingStrength, NumberOfIterations,
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
    return (r_glob, MeanRinEachIteration, acceptanceRateRewiring, finalAdj,
            finalY, final_adj_mat, sumKin, sumKin_bin_means, sliceOfOmega)

# ===================================== --------- ========================

make()
rewire, NumberOfNodes, NumberOfEdges, couplingStrength, \
NumberOfIterations, NumbertOfSteps, tfinal, tinitial, \
dt, times, graph, numberOfBins = configParameters()

Adj, adj_mat, InitialSumKin, \
Omega, InitialCondition = configGraph()

currentTime = datetime.datetime.now()
hour = currentTime.hour;
start = time()
#motherDirPath = "/storage/users/fbaharifard/ComplexNetworks/CompleteData"
motherDirPath = "/home/vahid/Documents/Complex network/c/CompleteData"

NumberOfSelfishNodes = 70
runNumber = 0
while(runNumber < 1): #hour < 14):
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
    obj = ode_solver.ODE(NumberOfNodes, tfinal, dt,	couplingStrength, \
                   InitialCondition, Omega, NumbertOfSteps, NumberOfIterations)
    sol = obj.integrate(Adj, rewire=rewire, currentPath=currentPath, NumberOfSelfishNodes=NumberOfSelfishNodes)
    sol = np.asarray(sol)
    r_glob, MeanRinEachIteration, acceptanceRateRewiring, \
    finalAdj,finalY, final_adj_mat, sumKin, sumKin_bin_means, \
    sliceOfOmega = getParametersFromC()
    print("simulation...done")
    del obj, sol
    finalList = {
        'FinalAdj.txt':final_adj_mat,
        #'acceptanceRateRewiring.txt':np.array(acceptanceRateRewiring).T,
        #'MeanRinEachIteration.txt':np.array(MeanRinEachIteration).T,
        'finalY.txt':np.array(finalY).T,
        #'r_glob.txt':np.array(r_glob).T,
        'sumKin.txt':np.array(sumKin).T,
        'sumKin_bin_means.txt':np.array(sumKin_bin_means).T,
        'sliceOfOmega.txt':np.array(sliceOfOmega).T
    }
    saveData(finalList)
    display_time(time()-start)
    currentTime = datetime.datetime.now()
    hour = currentTime.hour;
    runNumber += 1

if (motherDirPath == "/home/vahid/Documents/Complex network/c/CompleteData"):
    os.system("mpg123 "+" finish.mp3")
