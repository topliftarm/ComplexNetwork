    
import ode_solver
import pylab as pl 
import numpy as np 
from time  import time
from network import make_graph
from modules import *
from numpy import pi
import matplotlib.pyplot as plt

print("configuration...\n")
seed = 12358
np.random.seed(seed)
start = time()

NumberOfNodes = 100
graph  = make_graph()
couplingStrength = 0.23

tfinal = 100.0
tinitial = 0.0
dt = 0.1
times = np.arange(0,tfinal, dt)

print("configuration...done\n")
print("create graph...\n")
#Adj = graph.erdos_renyi_graph(NumberOfNodes, 0.1)
Adj = graph.random_k_out_graph(NumberOfNodes, 6, seed=np.random.randint(100000))
adj_mat = np.asarray(Adj).reshape((NumberOfNodes, NumberOfNodes))
#degree = np.sum(adj_mat, axis=1)

Omega = np.random.uniform(0, 2*pi, size=NumberOfNodes).tolist()
InitialCondition = np.random.uniform(0, 2*pi, NumberOfNodes).tolist()
print("create graph...done")
print("start simulation...")
obj = ode_solver.ODE(NumberOfNodes, tfinal, dt,	couplingStrength, InitialCondition, Omega)
#obj.set_matrices(Adj)
sol = obj.integrate(Adj)
sol = np.asarray(sol)
r_glob, psi = obj.get_order_parameters()
acceptanceRateRewiring = obj.getAcceptanceRewiring()
MeanYPrime = obj.getMeanYPrime();


# -----------------------  Plotes -----------------
sizeFig = (9,7)
_dpi = 80
plt.figure(figsize=sizeFig, dpi=_dpi)
plt.plot(acceptanceRateRewiring)
plt.title('acceptanceRateRewiring_selfish')
#plt.savefig('acceptanceRateRewiring_selfish.png')

plt.figure(figsize=sizeFig, dpi=_dpi)
plt.plot(MeanYPrime)
plt.title('MeanYPrime_selfish')
#plt.savefig('MeanYPrime_selfish.png')

plt.figure(figsize=sizeFig, dpi=_dpi)
plt.plot(r_glob)
plt.title('r_glob_selfish')
#plt.savefig('r_glob_selfish.png')
# plt.show()

print("simulation...done")
del obj, sol
# ----------------------------------------- #
display_time(time()-start)
# ----------------------------------------- #
#fig,axs = pl.subplots(1,figsize=(8,6))
#axs.plot(g, R,   marker=">", markersize=3, lw=2, c='r', label=r'$R_{global}$')
#axs.set_ylabel("R")
#axs.legend()
#axs.set_xlabel(r"$K$",  fontsize=13)
#pl.savefig(adrsf+"R.png")
# pl.close()
# pl.show()
