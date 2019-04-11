import numpy as np 
import networkx as nx
import ode_solver
import pylab as pl
import jgraph as igraph
from scipy.stats.stats import pearsonr 
from mpl_toolkits.axes_grid1 import make_axes_locatable
# pl.switch_backend('agg')

adrsf="../data/f/"
adrs="../data/"

#-------------------------------------------------------------------#    
def comm_weighted(data):
    '''Community structure based on the multilevel
       algorithm of Blondel et al.'''
  
    conn_indices = np.where(data)
    weights = data[conn_indices]
    edges = zip(*conn_indices)
    G = igraph.Graph(edges=edges, directed=False)
    G.es['weight'] = weights
    comm = G.community_multilevel(weights=weights, return_levels=False)
    return comm
#-------------------------------------------------------------------#    
def comm_unweighted(data):
    '''Community structure based on the multilevel
       algorithm of Blondel et al.'''

    conn_indices = np.where(data)
    edgs = zip(*conn_indices)
    G = igraph.Graph(edges=edgs, directed=False)
    comm = G.community_multilevel(weights=None, return_levels=False)
    return comm
#-------------------------------------------------------------------#    
def calculate_NMI(comm1,comm2):
    '''Compares two community structures using normalized 
    mutual information as defined by Danon et al (2005)'''

    nmi = igraph.compare_communities(comm1, comm2, method='nmi', remove_none=False)
    return nmi

#--------------------------------------------------------------#
def imshow_plot(data):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    fig = pl.figure(100,figsize=(6,6))
    pl.clf()
    ax = pl.subplot(111)
    im = ax.imshow(data,interpolation='nearest', cmap='afmhot') # , cmap=pl.cm.ocean
    ax.invert_yaxis()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pl.colorbar(im, cax=cax)
#--------------------------------------------------------------#
def drop_array(x, xcut):
        for count, elem in enumerate(x):
            if elem > xcut:
                return count

#---------------------------------------------------------#
def find_dominant_frequency(xf,yf, fcut=0.0, 
                            xflimit=1.e-5, yflimit=0.09):

    index = drop_array(xf, fcut)
    xf = xf[index:]
    yf = yf[index:]
    # a = yf > yflimit
    # f = xf > xflimit
    # xf = xf[a & f]
    # yf = yf[a & f]
    try:
        freq = xf[np.argmax(yf)]
        pxx = np.max(yf)
    except:
        freq, pxx = [0, 0]

    return freq, pxx

#---------------------------------------------------------#
def find_frequencies(x, y, fs, method='FFT', scaling='density'):
    '''
    find power spectrum of a signal and return 
    frequencies and power spectrum
    '''
    if method == 'FFT':
        from scipy.fftpack import fft
        N = len(x)
        dt0 = x[2] - x[1]
        xf = np.linspace(0.0, 1.0 / (2.0 * dt0), N / 2)
        # xf = np.linspace(0, (N/2+1)*fs/N, N/2)
        yf = fft(y)
        yfplot = 2.0 / N * np.abs(yf[0:N // 2])
        return xf[1:], yfplot[1:]
    elif method == 'welch':
        from scipy.signal import welch
        f, pxx = welch(y, fs=fs, nperseg=len(y), scaling=scaling) # scaling='density', scaling='spectrum'
        # pl.semilogy(f, pxx)
        return f, pxx
#--------------------------------------------------------------#
def binarize(data, threshold):
        data = np.asarray(data)
        upper, lower = 1,0
        data = np.where(data>=threshold,upper,lower)
        return data
#--------------------------------------------------------------#
def calculate_correlations(sol, cor_step, dt, g, tcut=100, 
                           binarizing=True, threshold=0.95, 
                           savefig=False, savetext=False):
    ''' 
    calculate Kuramoto correlation matrices at times
    determined by cor_step

    '''

    def f_savetext():
        np.savetxt(adrsf+"Cor-"+str("%.2ff"%g)+'-'+\
                   str("%08.2f"%(step*dt))+\
                   ".txt", Cor, fmt='%15.9f')
    def f_savefig():
        # pl.savefig(adrsf+"Cor-"+str("%09.6f"%g)+'-'+str("%09.6f"%tau)+\
        #            '-'+str("%08.2f"%(step*dt))+'-'+str("%09.6f"%omega)+'-'+\
        #            str("%02d"%sim)+".png")
        pl.savefig(adrsf+"Cor-"+str("%.2f"%g)+'-'+\
                   '-'+str("%08.2f"%(step*dt))+'-'+".png")
        
        pl.close()

    rows ,colms = sol.shape
    # print rows, colms
    n = int(cor_step/dt);
    ntcut = int(tcut/dt)
    # print "n, ntcut" , n, ntcut
    if n==0 :
        n +=1
    cor_ave = 0.0
    
    counter=0
    for  step in range(ntcut, rows):
        if (step%n == 0) & (step!=0):
            x   = sol[step,:]
            Cor = ode_solver.kuramoto_correlation(tuple(x))
            # cor_ave += np.asarray(Cor)
            counter +=1
            if binarizing:
                Cor = binarize(Cor, threshold)
            if savetext:
                f_savetext()    
            if savefig:
                imshow_plot(Cor)
                f_savefig()
        # cor_ave /= float(counter)
        # np.savetxt(adrsf+'cor_ave.txt', cor_ave, fmt="%15.6f")
        # imshow_plot(cor_ave)
        # pl.savefig(adrsf+"cor_ave.png")
    inv_counter = 1.0/counter                
    # return np.asarray(cor_ave*inv_counter)
    return np.asarray(Cor)
#--------------------------------------------------------------#
def save_plot_xy(x1, y1, name, axs, SAVETXT=False, SAVEFIG=False):
    "save and plot x-y data"
    fig, ax = pl.subplots(1, figsize=(15, 5))
    # ax[0].plot(x2, y2, c="k", lw=0.5)
    # ax[0].set_xlim(np.min(x2), np.max(x2))
    ax.plot(x1,y1, c="r", lw=2)
    ax.set_xlim(axs[0], axs[1])
    ax.set_ylim(axs[2], axs[3])
    pl.savefig(name+".png")
    if SAVETXT:
        np.savetxt(name+".txt", 
                   zip(x1, y1),
                   fmt="%15.9f")
    pl.close()


#--------------------------------------------------------------#
def plot_R(R, Y, X, name="R0"):
    ''' plot time average of order parameter in 2D plane of 
    X and Y axises '''

    np.savetxt(adrsf+name+".txt", R, fmt="%15.9f")
    x_step = X[1]- X[0] 
    y_step = Y[1]  - Y[0]
    

    f, ax       = pl.subplots(1,figsize=(10,10))
    im          = ax.imshow(R, interpolation='nearest', cmap='afmhot')
    ax.invert_yaxis()

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pl.colorbar(im, cax=cax)
    step = 10
    ax.set_xticks(np.arange(0, len(X), step));
    ax.set_xticklabels(str("%.1f"%i)for i in X[0::step]);
    ax.set_yticks(np.arange(0, len(Y), step));
    ax.set_yticklabels(str("%.1f"%i)for i in Y[::step]);
    ax.set_ylabel(r"$g$" ,fontsize=16)
    ax.set_title(r"$R$")
    ax.set_xlabel(r"$\omega_0$",fontsize=16)

    pl.savefig(adrsf+name+".png")
    pl.close()
#--------------------------------------------------------------#
def display_time(time):
    ''' 
    show real time elapsed
    '''
    hour = int(time/3600);
    minute = (int(time % 3600))/60;
    second = time-(3600.*hour+60.*minute);
    print("Done in %d hours %d minutes %09.6f seconds"%(hour,minute,second))
#--------------------------------------------------------------#
def savetex(x,y,name):
    ofile = open(adrsf+name+".txt","w")
    row,col = y.shape
    for i in range(col):
        ofile.write("%.6f" % x[i])
        for j in range(row):
            ofile.write("%15.6f" % y[j,i])
        ofile.write("\n")
    ofile.close()
# ----------------------------------------- #
