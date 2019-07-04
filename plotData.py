import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys

mpl.rcParams['agg.path.chunksize'] = 10000

# ========================================================
def plot_rglob():
    Name = 'r_glob'
    fileName1 = Name + '.txt'
    #DirPath1 = "/home/vahid/Documents/Complex network/c/CompleteData/Run50"
    DirPath1 = dirPath+'/'+dirName
    data1 = pd.read_csv(DirPath1+'/'+fileName1, header=None)
    plt.figure(figsize=(15, 15), dpi=120)
    plt.plot(data1)
    plt.title(Name)
    print(DirPath1+'/'+Name+'.png')
    plt.savefig(DirPath1+'/'+Name+'.png')
    #plt.show()

# --------------------------------------------------------
def plot_MeanRinEachIteration():
    Name = 'MeanRinEachIteration'
    fileName1 = Name + '.txt'
    DirPath1 = dirPath+'/'+dirName
    data1 = pd.read_csv(DirPath1+'/'+fileName1, header=None)
    plt.figure(figsize=(15, 15), dpi=120)
    plt.plot(data1, '-r')
    plt.title(Name)
    plt.savefig(DirPath1+'/'+Name+'.png')
    # plt.show()

# --------------------------------------------------------
def plot_acceptanceRateRewiring():
    Name = 'acceptanceRateRewiring'
    fileName1 = Name + '.txt'
    DirPath1 = dirPath+'/'+dirName
    data1 = pd.read_csv(DirPath1+'/'+fileName1, header=None)
    plt.figure(figsize=(15, 15), dpi=120)
    plt.plot(data1, '-r')
    plt.title(Name)
    plt.savefig(DirPath1+'/'+Name+'.png')
    #plt.show()

# --------------------------------------------------------
def plot_sumKinBinMeans():
    Name = 'sumKin_bin_means'
    fileName1 = Name + '.txt'
    fileName2 = 'sliceOfOmega.txt'
    DirPath1 = dirPath+'/'+dirName
    data1_mean_degree = pd.read_csv(DirPath1+'/'+fileName1, header=None)
    data_omega = pd.read_csv(DirPath1+'/'+fileName2, header=None)
    plt.figure(figsize=(15, 15), dpi=120)
    plt.plot(data_omega, data1_mean_degree, 'r*')
    plt.title(Name)
    plt.savefig(DirPath1+'/'+Name+'.png')
    # plt.show()

# --------------------------------------------------------
def plot_sumKin():
    Name = 'sumKin'
    fileName1 = Name + '.txt'
    fileName2 = 'Omega.txt'
    fileName3 = 'selfishNodes.txt'
    DirPath1 = dirPath+'/'+dirName
    data1_degree = pd.read_csv(DirPath1+'/'+fileName1, header=None)
    data_omega = pd.read_csv(DirPath1+'/'+fileName2, header=None)
    data_selfishName = np.array(pd.read_csv(DirPath1+'/'+fileName3, header=None).values).reshape(-1)
    plt.figure(figsize=(15, 15), dpi=120)
    plt.plot(data_omega, data1_degree, 'r*', data_omega[0][data_selfishName], data1_degree[0][data_selfishName], 'g*',)
    plt.legend(('Non Selfish Nodes', 'Selfish Nodes'), loc='upper right')
    plt.title(Name)
    plt.savefig(DirPath1+'/'+Name+'.png')
    # plt.show()

# ========================================================
dirName = sys.argv[1]

if (len(sys.argv) == 1):
    print("Please Insert Dir. Name")
    sys.exit(1)
elif (len(sys.argv)==2):
    print("argv==2")
    dirPath = "/home/vahid/Documents/Complex network/c/CompleteData"
elif (len(sys.argv)==3):
    print("argv==3")
    dirPath = sys.argv[2]


print(dirPath+'/'+dirName)
plot_rglob()
plot_MeanRinEachIteration()
plot_acceptanceRateRewiring()
plot_sumKinBinMeans()
plot_sumKin()

plt.show()
# ========================================================








# fileName1 = 'r_glob.txt'
# DirPath = "/home/vahid/Documents/Complex network/c/CompleteData/Run0"
# data1 = pd.read_csv(DirPath+'/'+fileName1, header=None)
# ##plt.figure()
# ##plt.plot(data1)
# DirPath_2 = "/home/vahid/Documents/Complex network/c/CompleteData/Run4"
# data2 = pd.read_csv(DirPath_2+'/'+fileName1, header=None)
# plt.figure()
# plt.plot(data1, '-r', data2, '-g')

#---------------------------------------------------

# DirPath_2 = "/home/vahid/Documents/Complex network/c/CompleteData/Run4"
# data2 = pd.read_csv(DirPath_2+'/'+fileName1, header=None)
# plt.figure()
# plt.plot(data1, '-r', data2, '-g')

#--------------------------------------------

# DirPath2 = "/home/vahid/Documents/Complex network/c/CompleteData/Run4"
# data2 = pd.read_csv(DirPath2+'/'+fileName1, header=None)
# plt.figure()
# plt.plot(data1, '-r', data2, '-g')

#---------------------------------------------------


# fileName1 = 'sumKin.txt'
# fileName2 = 'Omega.txt'
# DirPath = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/Run1"
# data1_degree = pd.read_csv(DirPath+'/'+fileName1, header=None)
# data_omega = pd.read_csv(DirPath+'/'+fileName2, header=None)
# DirPath_2 = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/Run2"
# data2_degree = pd.read_csv(DirPath_2+'/'+fileName1, header=None)
# plt.figure()
# plt.plot(data1_degree, abs(data_omega), 'r*', data2_degree, abs(data_omega), 'g*')

#---------------------------------------------------

#fileName1 = 'sumKin.txt'
#fileName2 = 'Omega.txt'
#fileName3 = 'InitialSumKin.txt'
#DirPath = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/Run2"
#data1_degree = pd.read_csv(DirPath+'/'+fileName1, header=None)
#data_omega = pd.read_csv(DirPath+'/'+fileName2, header=None)
#DirPath_2 = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/Run2"
#data2_degree = pd.read_csv(DirPath_2+'/'+fileName3, header=None)
#plt.figure()
#plt.plot(data1_degree, abs(data_omega), 'r*', data2_degree, abs(data_omega), 'g*')


#---------------------------------------------------

#fileName1 = 'acceptanceRateRewiring.txt'
#DirPath = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/Run1"
#data = pd.read_csv(DirPath+'/'+fileName1, header=None)
#plt.figure()
#plt.plot(data, '+-')


#---------------------------------------------------


# fileName1_2 = "/home/vahid/Documents/Complex network/c/CompleteData/result8/rewiring/r_glob.txt"
# data1_2 = pd.read_csv(fileName1_2, header=None)
# plt.figure()
# plt.plot(data1,'r', data1_2, 'g')

#fileName1_2 = '/home/vahid/Documents/Complex network/c/CompleteData/global-m(2)-uniform-1000-k(0.3)/r_glob.txt'

# fileName2 = 'MeanRinEachIteration.txt'

# fileName4 = 'sumKin_bin_means.txt'
#
#DirPath = "/home/vahid/Documents/Complex network/c/CompleteData/result9/rewiring"
#fileName5 = 'sumKin.txt'
#fileName6 = 'Omega.txt'
# data5 = pd.read_csv(DirPath+'/'+fileName5, header=None)
# data6 = pd.read_csv(DirPath+'/'+fileName6, header=None)
# plt.figure()
# plt.plot(data5, abs(data6), 'r*')

#data5 = pd.read_csv(DirPath+'/'+fileName5, header=None)
#data6 = pd.read_csv(DirPath+'/'+fileName6, header=None)
#data5_2 = pd.read_csv(DirPath_2+'/'+fileName5, header=None)
#data6_2 = pd.read_csv(DirPath_2+'/'+fileName6, header=None)
#plt.figure()
#plt.plot(data5, abs(data6), 'r*', data5_2, abs(data6_2), 'g*')


#data2 = pd.read_csv(DirPath+'/'+fileName2, header=None)

# DirPath = "/home/vahid/Documents/Complex network/c/CompleteData/result6/rewiring"
# fileName3 = 'acceptanceRateRewiring.txt'
# data3 = pd.read_csv(DirPath+'/'+fileName3, header=None)
# plt.plot(data3, '-*')

# DirPath_2 = "/home/vahid/Documents/Complex network/c/CompleteData/result5/rewiring"
# data3_2 = pd.read_csv(DirPath_2+'/'+fileName3, header=None)
# plt.figure()
# plt.plot(data3, '-*r', data3_2, '-*g')

# data4 = pd.read_csv(DirPath+'/'+fileName4, header=None)

#
#

# plt.figure()
# plt.plot(data2)
# plt.figure()
# plt.plot(data4, '*')

# plt.figure()
# plt.plot(data6, '*')


#plt.show()
