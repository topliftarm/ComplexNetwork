import matplotlib.pyplot as plt;
import numpy as np
import pandas as pd
import networkx as nx
import os, sys

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


fileName3 = 'selfishNodes.txt'
#data_selfishName = np.array(pd.read_csv(fileName3, header=None).values).reshape(-1).tolist()
data_selfishName = np.loadtxt(DirPath1+'/'+fileName3).reshape(-1).tolist()

fileName_adj = 'InitAdj.txt'
Adj = np.loadtxt(DirPath1+'/'+fileName_adj).reshape(-1, 105*105).ravel().tolist()
adj_mat = np.asarray(Adj).reshape((105,105))

G1=nx.from_numpy_matrix(adj_mat, create_using=nx.MultiDiGraph() )

pos = nx.spring_layout(G1)


fileName_adj = 'FinalAdj.txt'
Adj = np.loadtxt(DirPath1+'/'+fileName_adj).reshape(-1, 105*105).ravel().tolist()
adj_mat = np.asarray(Adj).reshape((105,105))

G2=nx.from_numpy_matrix(adj_mat, create_using=nx.MultiDiGraph())

plt.figure(figsize=(22, 22), dpi=120)
fig1 = nx.draw_networkx(G1, pos = pos, node_size=140, font_size=9, arrows=True)
nx.draw_networkx_nodes(G1,pos,
                       nodelist=data_selfishName,
                       node_color='g',
                       node_size=140,
                       alpha=0.8)
plt.axis('off')
plt.savefig(DirPath1+'/'+"init.png", format = "png", dpi = 450)

plt.figure(figsize=(22, 22), dpi=120)
fig2 = nx.draw_networkx(G2, pos = pos, node_size=140, font_size=9, arrows=True)
nx.draw_networkx_nodes(G2,pos,
                       nodelist=data_selfishName,
                       node_color='g',
                       node_size=140,
                       alpha=0.8)
plt.axis('off')
plt.savefig(DirPath1+'/'+"final.png", format = "png", dpi = 450)


plt.show()
