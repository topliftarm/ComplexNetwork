import pandas as pd
import matplotlib.pyplot as plt
import os

# fileName1 = 'r_glob.txt'
# #DirPath = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/Run0"
# DirPath1 = "/home/vahid/Documents/Complex network/c/CompleteData/Run0"
# data1 = pd.read_csv(DirPath1+'/'+fileName1, header=None)
# plt.figure()
# plt.plot(data1)

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

fileName1 = 'MeanRinEachIteration.txt'
DirPath1 = "/home/vahid/Documents/Complex network/c/CompleteData/Run0"
data1 = pd.read_csv(DirPath1+'/'+fileName1, header=None)
plt.figure()
plt.plot(data1, '-r')

# DirPath_2 = "/home/vahid/Documents/Complex network/c/CompleteData/Run4"
# data2 = pd.read_csv(DirPath_2+'/'+fileName1, header=None)
# plt.figure()
# plt.plot(data1, '-r', data2, '-g')

#--------------------------------------------
fileName1 = 'acceptanceRateRewiring.txt'
DirPath1 = "/home/vahid/Documents/Complex network/c/CompleteData/Run0"
data1 = pd.read_csv(DirPath1+'/'+fileName1, header=None)
plt.figure()
plt.plot(data1, '-r')
# DirPath2 = "/home/vahid/Documents/Complex network/c/CompleteData/Run4"
# data2 = pd.read_csv(DirPath2+'/'+fileName1, header=None)
# plt.figure()
# plt.plot(data1, '-r', data2, '-g')


#--------------------------------------------

#fileName1 = 'sumKin.txt'
#fileName2 = 'Omega.txt'
#DirPath = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/Run1"
#data1_degree = pd.read_csv(DirPath+'/'+fileName1, header=None)
#data_omega = pd.read_csv(DirPath+'/'+fileName2, header=None)
#DirPath_2 = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/Run2"
#data2_degree = pd.read_csv(DirPath_2+'/'+fileName1, header=None)
#plt.figure()
#plt.plot(data1_degree, abs(data_omega), 'r*', data2_degree, abs(data_omega), 'g*')

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

fileName1 = 'sumKin_bin_means.txt'
fileName2 = 'sliceOfOmega.txt'
#DirPath = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/Run1"
DirPath1 = "/home/vahid/Documents/Complex network/c/CompleteData/Run0"
data1_mean_degree = pd.read_csv(DirPath1+'/'+fileName1, header=None)
data_omega = pd.read_csv(DirPath1+'/'+fileName2, header=None)
plt.figure()
plt.plot(data_omega, data1_mean_degree, 'r*')

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


plt.show()
