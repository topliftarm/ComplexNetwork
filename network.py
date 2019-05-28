import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import jgraph as igraph
# import community
import networkx as nx
from copy import copy
from sys import exit
import ode_solver
from modules import adrsf
# np.set_printoptions(threshold=np.nan)



class make_graph:
    ''' make different graphs ans return their adjacency matrices
    as a 1 dimensional double vector in stl library'''

    def __init__(self):
        self.G = 0

    #---------------------------------------------------------------#

    def convert_to_1D_vector(self, G):
        ''' ----- '''
        N = self.N
        A = nx.to_numpy_matrix(G)
        A = A.reshape(N*N).tolist()[0]
        M = ode_solver.DoubleVector(N*N)
        for i in range(N*N):
            M[i] = A[i]

        return M

    def print_adjacency_matrix(self, G):
        '''print adjacency matrix of graph on the screen '''
        M = nx.to_numpy_matrix(G)
        M = np.array(M)
        for row in M:
            for val in row:
                print('%.0f'%val)
            print()

    #---------------------------------------------------------------#
    def complete_graph(self, N, plot_adj=False):
        ''' returns all to all graph '''

        self.N = N
        self.G = nx.complete_graph(N)
        M = nx.to_numpy_matrix(self.G)

        if plot_adj:
            self.imshow_plot(M, "con")

        return tuple(np.asarray(M).reshape(-1))
    #---------------------------------------------------------------#
    def complete_weighted(self, N, clusters, weights, d, plot_adj=False):
        '''
        return a complete weighted graph, weights distribute in clusters
        with the same size.
        weights: list with length 2, shows the weights of edges in and between clusters, respectivly
        clusters : list of cluster lengths
        d : delay for all nodes
        '''

        self.N = N
        # I use delay as weight
        self.modular_graph(N, 1, 1, clusters, weights[0],weights[1])
        M = self.D
        if plot_adj:
            self.imshow_plot(np.asarray(M).reshape(N,N), "con")

        self.D   = (d,)*N*N
        return   M
    #---------------------------------------------------------------#
    def complete_hmn_weighted(self, n_M0, level, weights, delay,
                              plot_adj=False, seed=124):
        '''
        return a complete weighted graph, weights distribute in a hierarchical
        modular form.
        Single delay for every node.
        '''
        self.hierarchical_modular_graph(n_M0, level, 1, 1, 1, weights,
                                        plot_adj=False, seed=124)
        N = self.N
        M = nx.to_numpy_matrix(self.G, weight='weight')

        if plot_adj:
            self.imshow_plot(M, "con")
        self.D = (delay,)*N*N

        return tuple(np.asarray(M).reshape(-1))

    #---------------------------------------------------------------#
    def erdos_renyi_graph(self, N, p, seed=123, directed=False, plot_adj=False):
        ''' returns Erdos Renyi graph '''

        self.N = N
        self.G = nx.erdos_renyi_graph(N,p,seed=seed, directed=False)
        M = nx.to_numpy_matrix(self.G)

        if plot_adj:
            self.imshow_plot(M, "con")

        return tuple(np.asarray(M).reshape(-1))

    #---------------------------------------------------------------#
    def random_k_out_graph(self, N, k, seed=123, plot_adj=False):
        self.N = N
        _alpha = 1
        self.G = nx.random_k_out_graph(N, k, _alpha, seed=seed)
        M = nx.to_numpy_matrix(self.G)
        if plot_adj:
            self.imshow_plot(M, "con")

        return tuple(np.asarray(M).reshape(-1))
    #---------------------------------------------------------------#

    def modular_graph(self, N, pIn, pOut, lengths, dIn, dOut, plot_adj=False):
        ''' returns a modular networks with :
        N : number of nodes
        pIn : conectivity of nodes insede clusters
        pOut: conectivity of nodes between clusters
        n_cluster :  number of clusters in graph
        dIn  : delay between nodes inside the clusters
        dOut : delay between nodes outside the clusters
        '''

        self.N = N
        M = np.zeros((N,N))
        D = np.zeros((N,N))
        n_cluster = len(lengths)

        for i in range(N):
            for j in range(i+1,N):
                r = np.random.rand()
                if r < pOut :
                    M[i,j] = M[j,i] = 1.0;
                    D[i,j] = D[j,i] = dOut

        # empty inside the clusters
        s = 0;
        for k in range(n_cluster):
            if k > 0:
                 s += lengths[k-1]
            for i in range(s,(lengths[k]+s)):
                for j in range(i+1,(lengths[k]+s)):
                    M[i,j] = M[j,i] = 0.0
                    D[i,j] = D[j,i] = 0.0

        # fill inside the clusters
        s = 0;
        for k in range(n_cluster):
            if k > 0:
                 s += lengths[k-1]
            for i in range(s,(lengths[k]+s)):
                for j in range(i+1,(lengths[k]+s)):
                    r = np.random.rand()
                    if r < pIn:
                        M[i,j] = M[j,i] = 1.0
                        D[i,j] = D[j,i] = dIn

        # print delay matrix
        def print_delay_matrix():
            ofi = open("delay.txt","w")
            for i in range(N):
                for j in range(N):
                    ofi.write("%2.0f"%D[i,j])
                ofi.write("\n")
            ofi.close()

        self.G = nx.from_numpy_matrix(M)
        self.D = tuple(D.reshape(-1))

        if plot_adj:
            self.imshow_plot(M, "con")
            # self.print_adjacency_matrix(self.G)

        return tuple(M.reshape(-1))
    #---------------------------------------------------------------#
    def hierarchical_modular_graph(self, n_M0, level, prob0, prob, alpha,
                                   delays, plot_adj=False, seed=125 ):
        '''
        n_M0 : size of module at level 1
        s    : number of levels
        n_modules : number of modules
        N : number of nodes
        ps : probability of conection in each level
             level one is 1 and the others determine with prob function
        delays: delay in each level as a list
        '''
        def probability(l,a=1,p=0.25):
            if l==0:
                from sys import exit
                print("level shold be integer and > 0", l)
                exit(0)
            else:
                return a* p**l

        s = level
        n_modules = int(2**(s-1))    #number of modules
        N = int(n_modules*n_M0)      #number of nodes
        self.N = N

        # M0 = nx.complete_graph(n_M0)
        M0 = nx.erdos_renyi_graph(n_M0, prob0, seed=seed)
        for e in nx.edges(M0):
            M0.add_edge(*e, weight=delays[0]) #delays act in weight attribute

        ps = [prob0]+[probability(i,alpha,prob) for i in range(1,s)]

        for l in range(1,s):
            if l==1:
                M_pre = [M0] * n_modules
            else:
                M_pre  = copy(M_next)
            M_next = []
            k = 0
            for ii in range(n_modules/(2**l)):
                step = 2**(l-1)
                tG   = nx.convert_node_labels_to_integers(M_pre[k+1], step*n_M0)
                tG1  = nx.compose(M_pre[k],tG)
                edge = 0
                effort = 1
                ''' make sure that connected modules are not isolated '''
                while edge<1:
                    # print "effort ", effort
                    effort +=1
                    for i in range(len(tG1)):
                        for j in range(i+1,len(tG1)):
                            if (i<step*n_M0) & (j>step*n_M0-1) & (np.random.rand()<ps[l]):
                                tG1.add_edge(i,j,weight=delays[l])
                                edge += 1

                M_next.append(tG1)
                k += 2
        self.G = M_next[0]

        if plot_adj:
            M = nx.to_numpy_matrix(self.G,weight=None)
            self.imshow_plot(M, "con")
            # self.print_adjacency_matrix(self.G)
        D = nx.to_numpy_matrix(self.G, weight='weight')
        self.D = tuple(np.asarray(D).reshape(-1))

        M = nx.to_numpy_matrix(self.G,weight=None)

        return tuple(np.asarray(M).reshape(-1))
    #---------------------------------------------------------------#
    def from_adjacency_matrix_graph(self,filename, plot_adj=False):
        ''' makes a graph from adjacency matrix in filename
        and return 1D double vector in stl '''

        A = np.genfromtxt(filename)
        self.N = len(A)

        if plot_adj:
            self.imshow_plot(A,"con")

        return tuple(A.reshape(-1))
    #---------------------------------------------------------------#
    def imshow_plot(self,data,name):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        fig = pl.figure(140,figsize=(6,6))
        pl.clf()
        ax = pl.subplot(111)
        im = ax.imshow(data,interpolation='nearest', cmap='afmhot') # , cmap=pl.cm.ocean
        ax.invert_yaxis()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        pl.colorbar(im, cax=cax)
        pl.savefig(adrsf+name+".png")
        # pl.close()
    #---------------------------------------------------------------#
    def multilevel(self, data):
        conn_indices = np.where(data)
        weights = data[conn_indices]
        edges = zip(*conn_indices)
        G = igraph.Graph(edges=edges, directed=False)
        G.es['weight'] = weights
        comm = G.community_multilevel(weights=weights, return_levels=False)
        return comm
    #---------------------------------------------------------------#
    def read_ctx(self, filename, delays, plot_adj=False):
        def check_symmetric(a, tol=1e-8):
            return np.allclose(a, a.T, atol=tol)
        data = np.genfromtxt(filename)
        data = data/float(np.max(data))
        if ~check_symmetric(data):
            data = np.maximum(data, data.transpose())
        print("Matrix is symmetric: %r"% check_symmetric(data))
        N = data.shape[0]
        self.N = N
        comm = self.multilevel(data)
        n_comm = max(comm.membership) + 1
        lengths = [len(comm[i]) for i in range(n_comm)]
        self.cluster_lengths = lengths
        print('There are %d communities. '% n_comm)

        #new indices of nodes-------------------------------------------
        newindices = []
        for i in range(n_comm):
            newindices += comm[i]
        # print newindices
        # --------------------------------------------------------------
        M = np.zeros((N,N))  # ordered connection matrix
        for i in range(N):
            for j in range(N):
                M[i, j] = data[newindices[i], newindices[j]]


        self.G = nx.Graph()
        # whole matrix
        for i in range(N):
            for j in range(N):
                if M[i,j]!=0:
                    self.G.add_edge(i,j, weight=M[i,j], delay=delays[1])
        # in clusters
        s = 0;
        for k in range(n_comm):
            if k > 0:
                 s += lengths[k-1]
            for i in range(s,(lengths[k]+s)):
                for j in range(i+1,(lengths[k]+s)):
                    if M[i,j] != 0:
                        self.G.add_edge(i, j, weight=M[i,j], delay=delays[0])

        if plot_adj:
            # M = nx.to_numpy_matrix(self.G,weight=None)
            A = nx.attr_matrix(self.G, edge_attr='delay',  rc_order=range(N))
            self.imshow_plot(A,'0_Delay')
            B = nx.attr_matrix(self.G, edge_attr='weight',  rc_order=range(N))
            self.imshow_plot(B,'0_Weight')

        self.D   = nx.attr_matrix(self.G, edge_attr='delay',  rc_order=range(N))
        M_weight = nx.attr_matrix(self.G, edge_attr='weight', rc_order=range(N))
        self.D = tuple(np.asarray(self.D).reshape(-1))
        return   tuple(np.asarray(M_weight).reshape(-1))



    def barabasi(self, N, m, seed=123, directed=False, plot_adj=False):
        ''' returns BA graph '''

        self.N = N
        self.G = nx.barabasi_albert_graph(N,m,seed=seed)
        M = nx.to_numpy_matrix(self.G)

        if plot_adj:
            self.imshow_plot(M, "con")

        return tuple(np.asarray(M).reshape(-1))








# --------------------------------------------------------------
# --------------------------------------------------------------
def old_modular_graph(self, N, pIn, pOut, n_cluster, dIn, dOut, plot_adj=False):
        ''' returns a modular networks with :
        N : number of nodes
        pIn : conectivity of nodes insede clusters
        pOut: conectivity of nodes between clusters
        n_cluster :  number of clusters in graph
        dIn  : delay between nodes inside the clusters
        dOut : delay between nodes outside the clusters
        '''
        self.N = N
        length = int(N/n_cluster)
        M = np.zeros((N,N))
        D = np.zeros((N,N))

        for i in range(N):
            for j in range(i+1,N):
                r = np.random.rand()
                if r < pOut :
                    M[i,j] = M[j,i] = 1.0;
                    D[i,j] = D[j,i] = dOut

        s = 0;
        for k in range(n_cluster):
            if k > 0:
                 s += length
            for i in range(s,(length+s)):
                for j in range(i+1,(length+s)):
                    r = np.random.rand()
                    if r < pIn:
                        M[i,j] = M[j,i] = 1.0
                        D[i,j] = D[j,i] = dIn


        # print delay matrix
        def print_delay_matrix():
            ofi = open("delay.txt","w")
            for i in range(N):
                for j in range(N):
                    ofi.write("%2.0f"%D[i,j])
                ofi.write("\n")
            ofi.close()


        self.G = nx.from_numpy_matrix(M)
        self.D = tuple(D.reshape(-1))

        if plot_adj:
            self.imshow_plot(M)
            self.print_adjacency_matrix(self.G)

        return tuple(M.reshape(-1))
