import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fileName1 = 'MeanRinEachIteration.txt' 
DirPath1 = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/LongResult/1/1"
DirPath2 = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/LongResult/1/2"
DirPath3 = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/LongResult/1/3"
DirPath4 = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/LongResult/1/4"
data1 = pd.read_csv(DirPath1+'/'+fileName1, header=None).to_numpy()
data2 = pd.read_csv(DirPath2+'/'+fileName1, header=None).to_numpy()
data3 = pd.read_csv(DirPath2+'/'+fileName1, header=None).to_numpy()
data4 = pd.read_csv(DirPath2+'/'+fileName1, header=None).to_numpy()
temp = np.append(data1, data2)
temp2 = np.append(data3, data4)
plt.figure()
plt.plot(np.append(temp, temp2))



#fileName1 = 'acceptanceRateRewiring.txt'
#DirPath1 = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/Run1_"
#DirPath2 = "/home/fatemeh/Documents/ComplexNetwork/CompleteData/Run1"  
#data1 = pd.read_csv(DirPath1+'/'+fileName1, header=None).to_numpy()
#data2 = pd.read_csv(DirPath2+'/'+fileName1, header=None).to_numpy()
#plt.figure()
#plt.plot(np.append(data1, data2[2:len(data2)]))




plt.show()
