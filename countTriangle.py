import numpy as np 
import networkx as nx 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys


def countTriangle(g, isDirected): 
    nodes = len(g) 
    count_Triangle = 0 #Initialize result 
    # Consider every possible triplet of edges in graph 
    for i in range(nodes): 
        for j in range(nodes): 
            for k in range(nodes): 
                # check the triplet if it satisfies the condition 
                if( i!=j and i !=k and j !=k and 
                        g[i][j] and g[j][k] and g[k][i]): 
                    count_Triangle += 1
    # if graph is directed , division is done by 3 
    # else division by 6 is done 
    return count_Triangle/3 if isDirected else count_Triangle/6


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
adj_mat = np.asarray(Adj).reshape((100,100))  

print(countTriangle(adj_mat, True))