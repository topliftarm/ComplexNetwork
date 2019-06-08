import matplotlib.pyplot as plt; 
import numpy as np
import networkx as nx

fileName_adj = 'InitAdj.txt' 
Adj = np.loadtxt(fileName_adj).reshape(-1, 100*100).ravel().tolist() 
adj_mat = np.asarray(Adj).reshape((100,100))      

G1=nx.from_numpy_matrix(adj_mat, create_using=nx.MultiDiGraph() )

pos = nx.spring_layout(G1)


fileName_adj = 'FinalAdj.txt' 
Adj = np.loadtxt(fileName_adj).reshape(-1, 100*100).ravel().tolist() 
adj_mat = np.asarray(Adj).reshape((100,100))      

G2=nx.from_numpy_matrix(adj_mat, create_using=nx.MultiDiGraph())

plt.figure(figsize=(22, 22), dpi=120)
fig1 = nx.draw_networkx(G1, pos = pos, node_size=140, font_size=9, arrows=True)
plt.savefig("init.png", format = "png", dpi = 450)

plt.figure(figsize=(22, 22), dpi=120)
fig2 = nx.draw_networkx(G2, pos = pos, node_size=140, font_size=9, arrows=True)
plt.savefig("final.png", format = "png", dpi = 450)


plt.show() 