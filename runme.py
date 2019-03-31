
import ode_solver
import pylab as pl 
import numpy as np 
from time  import time
from network import make_graph
from modules import *
from numpy import pi


print("configuration...\n")
seed = 12358
np.random.seed(seed)
start = time()

NumberOfNodes = 100
graph  = make_graph()
couplingStrength = 0.23

tfinal = 10000.0
tinitial = 0.0
dt = 0.1
times = np.arange(0,tfinal, dt)

#R = np.zeros(times)

PLOT_RT = True # False
if PLOT_RT:
	figrt, axrt = pl.subplots(1, figsize=(7, 7))
	rt = []
print("configuration...done\n")
print("create graph...\n")
Adj = graph.random_k_out_graph(NumberOfNodes, 6, seed=np.random.randint(100000))
adj_mat = np.asarray(Adj).reshape((NumberOfNodes, NumberOfNodes))
#degree = np.sum(adj_mat, axis=1)

Omega = np.random.uniform(-1, 1, size=NumberOfNodes).tolist()
InitialCondition = np.random.uniform(-pi/2.0, pi/2.0, NumberOfNodes).tolist()
print("create graph...done")
print("start simulation...")
obj = ode_solver.ODE(NumberOfNodes, tfinal, dt,	couplingStrength, InitialCondition, Omega)
#obj.set_matrices(Adj)
sol = obj.integrate(Adj)
sol = np.asarray(sol)
r_glob, psi = obj.get_order_parameters()
print("simulation...done")
if PLOT_RT:
     axrt.plot(times, r_glob, lw=2)
     rt.append(r_glob)
 
#R[i] += np.mean(r_glob[id_tcut:])
del obj, sol


#R *= inv_num_sim

if PLOT_RT:
	axrt.set_ylabel("R(t)")
	axrt.set_xlabel("time")
	pl.savefig(adrsf+"rt.png")
	rt = np.asarray(rt)
	rt = rt.T
	rtfile = open(adrsf+"rt.txt", "w")
	for i in range(len(times)):
		rtfile.write("%15.4f" % times[i])
		rtfile.write("\n")
	rtfile.close()
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
