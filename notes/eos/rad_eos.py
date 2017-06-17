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



#--------------------------------------------------
def main(argv):

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=9)
    

    #fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
    fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
    gs = plt.GridSpec(1, 1)

    N = 1000

    axs = []
    for i in range(1):
        ax = plt.subplot(gs[i, 0])
        ax.minorticks_on()
        ax.set_xlim(0.0, 12.0)
        ax.set_ylim(1.0e14, 1.0e42)
        ax.set_yscale('log')
        axs.append( ax )


    #axs[0].set_xlabel(r' $\rho$ (g cm$^{-3}$)')
    axs[0].set_xlabel(r'Radius $R_{\mathrm{SLy}, 1.4}$ (km)')
    axs[0].set_ylabel(r'Pressure $P$ (dyne cm$^{-2}$)')


    # basic parameters
    rhos = np.logspace(-4, 15.0, N)
    T = 1.0 * kelvin_per_keV #about 10 million kelvins
    Z = 56.0 #hydrogen plasma

    #new radius axis
    rads = R_at_rho( rhos )



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
    axs[0].plot(rads, con.d['P'], 'g-')

    con = container(N)
    for i, rho in enumerate( rhos ):
        p = plasma(rho, T*0.5, Z)
        con.append(p)
    axs[0].plot(rads, con.d['P'], 'b-')

    #main plot
    con = container(N)
    for i, rho in enumerate( rhos ):
        p = plasma(rho, T, Z)
        con.append(p)
    axs[0].plot(rads, con.d['P'], 'k-', linewidth=1.5, alpha=0.8)



    rho2 = np.logspace(3, 9, 100)
    Pnonrel = Ppoly(rho2, 5/3)
    Prel = Ppoly(rho2, 4/3)
    #axs[0].plot(rho2, Pnonrel, "b", linestyle='dotted')
    #axs[0].plot(rho2, Prel, "b", linestyle='dotted')

    #Core
    rho3 = np.logspace(13.0, 16, 100)
    rads3 = R_at_rho(rho3)
    #for i, eoss in enumerate(['PAL6', 'SLy', 'FPS', 'APR1', 'WFF1', 'BBB2', 'ENG', 'MS1', 'PS', 'H1', 'ALF1']):
    #    P1 = EoS(eoss, rho3)
    #    axs[0].plot(rho3, P1, "k-")

    for key, value in eosLib.iteritems():
        P1 = EoS(key, rho3)
        if value[4] == 'npem':
            axs[0].plot(rads3, P1, "k-", alpha = 0.3)
        if value[4] == 'meson':
            axs[0].plot(rads3, P1, "b-", alpha = 0.3)
        if value[4] == 'hyperon':
            axs[0].plot(rads3, P1, "g-", alpha = 0.3)
        if value[4] == 'quark':
            axs[0].plot(rads3, P1, "r-", alpha = 0.3)



    #find transition zones
    rho_xr  = con.rho('xr', 1.0)
    rho_tr  = con.rho('tr', 1.0)
    rho_thr = con.rho('thr', 1.0)
    rho_phir= con.rho('phir', 1.0)
    rho_Ge  = con.rho('Ge', 1.0)
    rho_Gi  = con.rho('Gi', 1.0)

    
    print """Transition zones:
            xr:   {:3.2e}
            tr:   {:3.2e}
            thr:  {:3.2e}
            phir: {:3.2e}
            Ge:   {:3.2e}
            Gi:   {:3.2e}
            """.format(rho_xr, rho_tr, rho_thr, rho_phir, rho_Ge, rho_Gi)


    #highlight different regions
    xmin = R_at_rho(1.0e1)
    xmax = R_at_rho(1.0e4)
    yarrow = 1.0e39
    col = 'lightgrey'
    axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e42, 1e42], facecolor=col, color=None, alpha=1.0, edgecolor=col)
    #plt.annotate(s='', xy=(xmax,yarrow), xytext=(xmin,yarrow), arrowprops=dict(arrowstyle='<->'))

    #outer crust
    #xmin = 1.0e4
    #xmax = rhond
    #axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e40, 1e40], facecolor=col, color=None, alpha=0.1, edgecolor=col)
    #plt.annotate(s='', xy=(xmax,yarrow), xytext=(xmin,yarrow), arrowprops=dict(arrowstyle='<->'))

    #inner crust
    xmin = R_at_rho(rhond)
    xmax = R_at_rho(0.5*rhon)
    axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e42, 1e42], facecolor=col, color=None, alpha=1.0, edgecolor=col)
    #plt.annotate(s='', xy=(xmax,yarrow), xytext=(xmin,yarrow), arrowprops=dict(arrowstyle='<->'))

    #outer core
    #xmin = 0.5*rhon
    #xmax = 2.0*rhon
    #axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e40, 1e40], facecolor=col, color=None, alpha=0.1, edgecolor=None)
    #plt.annotate(s='', xy=(xmax,yarrow), xytext=(xmin,yarrow), arrowprops=dict(arrowstyle='<->'))

    #inner core
    xmin = R_at_rho(2.0*rhon)
    xmax = R_at_rho(100.0*rhon)
    #axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e40, 1e40], facecolor=col, color=None, alpha=0.1, edgecolor=None)
    #plt.annotate(s='', xy=(xmax,yarrow), xytext=(xmin,yarrow), arrowprops=dict(arrowstyle='<->'))

    #different ns structures
    #y_text = 1.0e40
    #axs[0].text(400.0,  y_text, 'Atmosphere', rotation=0, ha='center', va='center', size=10)
    #axs[0].text(1.0e8,  y_text, 'Outer crust', rotation=0, ha='center', va='center', size=10)
    #axs[0].text(6.0e12, y_text, 'Inner\ncrust', rotation=0, ha='center', va='center', size=10)
    #axs[0].text(6.0e15, y_text, 'Core', rotation=0, ha='center', va='center', size=10)

    
    #lstyle = 'dotted'
    #axs[0].plot([rhon, rhon], [1.0e16, 1e40], "r", linestyle=lstyle)
    #txt = axs[0].text(rhon, 2.0e24, r'$\rho_n$', rotation=90, ha='center', va='center', size=8)
    #txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))

    #axs[0].plot([rhond, rhond], [1.0e16, 1e40], "r", linestyle=lstyle)
    #txt = axs[0].text(rhond*0.5, 2.0e24, r'$\rho_{\mathrm{ND}}$', rotation=90, ha='center', va='center', size=8)
    ##txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))

    #axs[0].plot([rho_xr, rho_xr], [1.0e16, 1e38], "r", linestyle=lstyle)
    #txt = axs[0].text(rho_xr, 2.0e28, r'$x_r = 1$', rotation=90, ha='center', va='center', size=8)
    #txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))

    #axs[0].plot([rho_thr, rho_thr], [1.0e16, 1e28], "r", linestyle=lstyle)
    #txt = axs[0].text(rho_thr, 2.0e22, r'$\Theta_r = 1$', rotation=90, ha='center', va='center', size=8)
    #txt.set_bbox(dict(facecolor='lightgrey', edgecolor='none', pad=3))

    #axs[0].plot([rho_Ge, rho_Ge], [1.0e16, 1e30], "r", linestyle=lstyle)
    #txt = axs[0].text(rho_Ge, 2.0e26, r'$\Gamma_{\mathrm{e}} = 1$', rotation=90, ha='center', va='center', size=8)
    #txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))


    #axs[0].text(1.0e4, 2.0e18, r'$P \propto \rho^{5/3}$', rotation=27, ha='center', va='center', size=8)
    #axs[0].text(1.0e7, 2.0e25, r'$P \propto \rho^{4/3}$', rotation=22, ha='center', va='center', size=8)


    #validity zones
    #yloc = 5.0e34
    #axs[0].text(5.e2, yloc, r'ideal gas', rotation=90, ha='center', va='center', size=8)
    #axs[0].text(1.e5, yloc, 'non-relativistic\ndegenerate\nelectron gas', rotation=90, ha='center', va='center', size=8)
    #axs[0].text(1.e8, yloc, 'relativistic\ndegenerate\nelectron gas', rotation=90, ha='center', va='center', size=8)
    #axs[0].text(1.e12, yloc, r'neutron drip', rotation=90, ha='center', va='center', size=8)
    #axs[0].text(5.e13, yloc, 'nuclear\npasta', rotation=90, ha='center', va='center', size=8)





if __name__ == "__main__":
    main(sys.argv)
    plt.savefig('radius_eos.pdf')

