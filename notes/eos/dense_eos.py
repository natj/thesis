from __future__ import division
from collections import defaultdict

import sys
import os


import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import AutoMinorLocator


print "Matplotlib version", matplotlib.__version__

from matplotlib import cm
import palettable as pal

#cmap = pal.cmocean.sequential.Matter_8.mpl_colormap #best so far
#cmap = pal.wesanderson.Zissou_5.mpl_colormap
cmap = pal.colorbrewer.qualitative.Set1_6.mpl_colormap
#cmap = plt.get_cmap('plasma_r')
#cmap = cm.get_cmap('inferno_r')

import numpy as np
from scipy.signal import savgol_filter
from scipy import interpolate

from units import *
from eos import *
from eoslibrary import EoS
from eoslibrary import eosLib



def main(argv):

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=7)
    

    fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
    #fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
    gs = plt.GridSpec(1, 1)

    N = 100

    axs = []
    for i in range(1):
        ax = plt.subplot(gs[i, 0])
        #ax.minorticks_on()
        #ax1.set_ylim(0.0, 1.1)
        ax.set_xlim(2.0e13, 2.0e16)
        ax.set_ylim(1.0e32, 1.0e39)
        ax.set_xscale('log')
        ax.set_yscale('log')

        axs.append( ax )

    axs[0].set_xlabel(r'Density $\rho$ (g cm$^{-3}$)')
    axs[0].set_ylabel(r'Pressure $P$ (dyne cm$^{-2}$)')


    # basic parameters
    rhos = np.logspace(-2, 14.2, N)
    T = 1.0 * kelvin_per_keV #about 10 million kelvins
    Z = 56.0 #hydrogen plasma


    #rho_ND
    rhond = 4.0e11

    #normal nuclear saturation density
    #rho_n
    rhon = 2.8e14



    print """ Basic parameters:
             ref pressure P_r: {:3.2e}
                 temp     T_r: {:3.2e}
                          T:   {:3.2e}
                          l_C: {:3.2e}

                          P_r^n {:3.2e}
                          T_r^n {:3.2e}

             -------------------------------------------------- 
                 """.format(Pr, Tr, T, lambdaC, Prn, Tn)


    con = container(N)
    for i, rho in enumerate( rhos ):
        p = plasma(rho, T*0.1, Z)
        con.append(p)
    axs[0].plot(con.d['rho'], con.d['P'], 'g-')

    con = container(N)
    for i, rho in enumerate( rhos ):
        p = plasma(rho, T*0.5, Z)
        con.append(p)
    axs[0].plot(con.d['rho'], con.d['P'], 'b-')

    #main plot
    con = container(N)
    for i, rho in enumerate( rhos ):
        p = plasma(rho, T, Z)
        con.append(p)
    axs[0].plot(con.d['rho'], con.d['P'], 'k-', linewidth=1.5, alpha=0.8)



    rho2 = np.logspace(3, 9, 100)
    Pnonrel = Ppoly(rho2, 5/3)
    Prel = Ppoly(rho2, 4/3)
    axs[0].plot(rho2, Pnonrel, "b", linestyle='dotted')
    axs[0].plot(rho2, Prel, "b", linestyle='dotted')

    #Core
    rho3 = np.logspace(13.0, 16, 100)

    alpha = 0.6
    for key, value in eosLib.iteritems():
        P1 = EoS(key, rho3)
        if value[4] == 'npem':
            axs[0].plot(rho3, P1, "k-", alpha = alpha)
        if value[4] == 'meson':
            axs[0].plot(rho3, P1, "b-", alpha = alpha)#, linestyle='dashdot')
        if value[4] == 'hyperon':
            axs[0].plot(rho3, P1, "g-", alpha = alpha)#, linestyle='dashed')
        if value[4] == 'quark':
            axs[0].plot(rho3, P1, "r-", alpha = alpha)#, linestyle='dotted')


    yarrow = 1.0e39
    col = 'lightgrey'

    #inner crust
    xmin = rhond
    xmax = 0.5*rhon
    axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e42, 1e42], facecolor=col, color=None, alpha=1.0, edgecolor=col)
    #plt.annotate(s='', xy=(xmax,yarrow), xytext=(xmin,yarrow), arrowprops=dict(arrowstyle='<->'))

    #outer core
    #xmin = 0.5*rhon
    #xmax = 2.0*rhon
    #axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e40, 1e40], facecolor=col, color=None, alpha=0.1, edgecolor=None)
    #plt.annotate(s='', xy=(xmax,yarrow), xytext=(xmin,yarrow), arrowprops=dict(arrowstyle='<->'))

    #inner core
    xmin = 2.0*rhon
    xmax = 100.0*rhon
    #axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e40, 1e40], facecolor=col, color=None, alpha=1.0, edgecolor=col)
    #plt.annotate(s='', xy=(xmax,yarrow), xytext=(xmin,yarrow), arrowprops=dict(arrowstyle='<->'))

    #different ns structures
    y_text = 1.0e38
    axs[0].text(5.0e13, y_text, 'Inner\ncrust', rotation=0, ha='center', va='center', size=8)
    axs[0].text(6.0e14, y_text, 'Core', rotation=0, ha='center', va='center', size=8)

    
    lstyle = 'dotted'
    axs[0].plot([rhon, rhon], [1.0e16, 1e40], "r", linestyle=lstyle)
    txt = axs[0].text(rhon, 2.0e36, r'$\rho_n$', rotation=90, ha='center', va='center', size=8)
    txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))

    axs[0].plot([rhond, rhond], [1.0e16, 1e40], "r", linestyle=lstyle)
    txt = axs[0].text(rhond*0.5, 2.0e36, r'$\rho_{\mathrm{ND}}$', rotation=90, ha='center', va='center', size=8)
    #txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))



if __name__ == "__main__":
    main(sys.argv)

    plt.subplots_adjust(left=0.15, bottom=0.16, right=0.98, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig('dense_eos.pdf')

