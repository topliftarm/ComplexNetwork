import ode_solver
import pylab as pl 
import numpy as np 
from time  import time
from network import make_graph as mg
from modules import *
from numpy import pi
from sys import exit
from copy import copy
# ----- Object creation -----
seed = 12358
np.random.seed(seed)
start = time()

# ----- Network Options -------------------- #
N = 100
gr  = mg()

couplingStrength = 0.23
# ----- Simulation Options ----------------- #
tfinal   = 500.0
tcut     = 100.0
tinitial = 0.0
dt       = 0.05
alpha = [0]
id_tcut  = int (tcut/dt)
g        = np.arange(5, 100, 5)
#num_sim  = 2
num_sim  = 1
inv_num_sim = 1.0/num_sim
times = np.arange(0,tfinal, dt)
ng = len(g)
links = 15
PLOT_RT = False # True


R = np.zeros((ng, num_sim))
R_k = np.zeros((ng, num_sim))
R1_std = np.zeros(ng)
R2_std = np.zeros(ng)

if PLOT_RT:
    fig0, ax = pl.subplots(1, figsize=(10, 6))
    # rt = []
    # rt_k = []

# Omega = np.random.normal(loc=0, scale=0.5, size=N).tolist()
# IC = np.random.uniform(-pi/2.0, pi/2.0, N).tolist()
for i in range(ng):
    for j in range(len(alpha)):
        for ens in range(num_sim):
            
            Adj = gr.barabasi(N, links, seed=np.random.randint(100000))
            # Adj = gr.erdos_renyi_graph(N, 0.03, seed=np.random.randint(100000))            
            adj_mat = np.asarray(Adj).reshape((N,N))
            degree = np.sum(adj_mat, axis=1)

            print("g = %8.3f, sim = %d" %(g[i], ens))

            #Omega = np.random.normal(loc=0, scale=0.5, size=N).tolist()
            Omega = copy(degree)
            # Omega = [0.0] * N 


            IC = np.random.uniform(-pi/2.0, pi/2.0, N).tolist()
            c = ode_solver.ODE( 
                            N, tfinal, 
                            tcut, dt, 
                            couplingStrength, IC, 
                            Omega, alpha[j])

            c.set_matrices(Adj)

            c.integrate()
            r_glob, r_glob_k = c.get_order_parameters()

            R[i, ens] = np.mean(r_glob[id_tcut:])
            R_k[i, ens] = np.mean(r_glob_k[id_tcut:])

            # if it is needed to plot r(t)
            if PLOT_RT:
                ax.plot(times, r_glob, lw=2)
                # rt.append(r_glob)
                # rt_k.append(r_glob_k)

            del c


if num_sim > 1:
    R1 =  np.mean(R, axis=1)
    R2 = np.mean(R_k, axis=1)
    R1_std = np.std(R, axis=1)
    R2_std = np.std(R_k, axis=1)
else:
    R1 = R[:,0]
    R2 = R_k[:,0]

if PLOT_RT:
	ax.set_ylabel("R(t)")
	ax.set_xlabel("time")
	pl.savefig(adrsf+"rt.png")
	# rt = np.asarray(rt)
	# rt = rt.T
	# rtfile = open(adrsf+"rt.txt", "w")
	# for i in range(len(times)):
	# 	rtfile.write("%15.4f" % times[i])
	# 	for j in range(num_sim):
	# 			rtfile.write("%15.6f " % rt[i, j])
	# 	rtfile.write("\n")
	# rtfile.close()
# ----------------------------------------- #
np.savetxt(adrsf+"R.txt", (g, R1, R2, R1_std, R2_std), fmt="%15.6f")
display_time(time()-start)
# ----------------------------------------- #
fig,axs = pl.subplots(1,figsize=(8,6))
if num_sim==1:
    axs.plot(g, R1,   marker=">", markersize=3, lw=2, c='r', label=r'$R$')
    axs.plot(g, R2,   marker=">", markersize=3, lw=2, c='k', label=r'$R_k$')
else:
    axs.errorbar(g, R1, yerr=R1_std, fmt='--o', markersize=4, lw=2, c='r', label=r'$R$')
    axs.errorbar(g, R2, yerr=R2_std, fmt='--o', markersize=4, lw=2, c='k', label=r'$R_k$')

axs.set_ylabel("Order parameter")
axs.legend(loc='best', fontsize=12)
axs.set_xlabel(r"$K$",  fontsize=13)
pl.savefig(adrsf+"R.png")
# pl.close()
# pl.show()
