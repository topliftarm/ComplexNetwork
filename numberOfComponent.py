import numpy as np 
import networkx as nx 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys

mpl.rcParams['agg.path.chunksize'] = 10000


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


	
DirPath1 = dirPath+'/'+dirName
fileName_adj = 'FinalAdj.txt' 

Adj = np.loadtxt(DirPath1+'/'+fileName_adj).reshape(-1, 100*100).ravel().tolist() 
init_adj_mat = np.asarray(Adj).reshape((100,100))  
print('degree=', np.sum(init_adj_mat, axis=1) )
G=nx.from_numpy_matrix(init_adj_mat) 
L = nx.normalized_laplacian_matrix(G).todense() 
a = np.linalg.eig(L) 

print(a)