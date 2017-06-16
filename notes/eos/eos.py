from __future__ import division
from collections import defaultdict

import sys
import os


import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
from eoslibrary import EoS


class plasma:

    def __init__(self, rho, T, Z):
        self.rho = rho
        self.T = T
        self.kT= kB * T
        self.Z = Z


        #electron number density
        self.ne = rho/mp

        #electron Fermi momentum
        self.pFe = hbar*(3*np.pi**2*self.ne)**(1/3)

        #dimensionless momentum
        self.xr = self.pFe / me / c

        #Fermi energy
        self.eFe = c**2 * np.sqrt( (me*c)**2 + self.pFe**2 )

        #relativistic Fermi temperature
        self.gammar = np.sqrt( 1 + self.xr**2 )
        self.TF = Tr * (self.gammar - 1)

        self.thetar =  self.T / self.TF

        #relativistic temperature
        self.tr  = self.T / Tr

        ####

        #electron sphere radius
        self.re = (4*np.pi*self.ne/3)**(-1/3)

        #Electron Coulomb coupling
        self.Ge = eC**2/(self.re*self.kT)

        #ion sphere radius
        self.ri = self.re * self.Z**(1/3)

        #Ion Coulomb coupling
        self.Gi = self.Ge * self.Z**(5/3)


        #degenerate neutron gas
        self.nn = rho/mn
        self.pFn = hbar*(3*np.pi**2*self.nn)**(1/3)
        self.yr = self.pFn / mn / c
        self.eFn = c**2 * np.sqrt( (mn*c)**2 + self.pFn**2 )
        self.deltar = np.sqrt( 1 + self.yr**2 )
        self.TFn = Tn * (self.deltar - 1)
        self.trn =  self.T / Tn
        self.phir = self.T / self.TFn
        

    #ideal non-relativistic Boltzmann gas
    def Pideal(self):
        return self.ne * self.kT

    #degenerate relativistic electron gas
    def Pdegenerate(self):
        t1 = self.xr*( (2/3)*self.xr**2 -1)*self.gammar
        t2 = np.log( self.xr + self.gammar )
        t3 = (4*np.pi**2)/9 *self.tr**2 * self.xr * (self.gammar + 1.0/self.gammar ) 

        return (Pr/8/np.pi**2) * (t1 + t2 + t3)

    #Ion coulomb attraction using the ion-sphere model
    def PCoulomb(self):
        n = self.ne/self.Z #number density of ions
        return -0.3 *n*self.Gi * self.kT


    #Skyrme crust equation of state
    # done with piecewise polytropes from Read 2009
    def PSly(self):

        #baryon density as defined in Read 2009
        rho = 1.66e-24 * self.nn
        #rho = self.rho

        #constants
        K1 = 6.80110e-9
        K2 = 1.06186e-6
        K3 = 5.32697e1
        K4 = 3.99874e-8

        #polytropic indices
        G1 = 1.58425
        G2 = 1.28733
        G3 = 0.62223
        G4 = 1.35692

        #transition depths
        r2 = 2.44034e7
        r3 = 3.78358e11
        r4 = 2.62780e12


        K = K4
        G = G4
        if rho < r4:
            K = K3
            G = G3
            if rho < r3:
                K = K2
                G = G2
                if rho < r2:
                    K = K1
                    G = G1


        return (K*rho**G)*c**2



    #return pressure depending on proper states of the plasma
    def pressure(self):

        #use (degenerate) electron gas up to neutron drip
        if self.rho < 1.0e8:
            return self.Pdegenerate() + self.PCoulomb() #+ self.Pneutron()
        else:
            #Skyrm crust
            return self.PSly()





    def __str__(self):
        s = """plasma at  
                rho: {:3.2e} g/cm^3 
                n_e: {:3.2e} 1/cm^3

                xr:   {:3.2f}
                tr:   {:3.2e}
                thr:  {:3.2e}

                Ge:   {:3.2e}
                Gi:   {:3.2e}

                P:    {:3.2e}

                 """.format(
                self.rho,
                self.ne,
                self.xr,
                self.tr,
                self.thetar,

                self.Ge,
                self.Gi,

                self.pressure(),
                  )

        return s
            

class container:
    params = ['rho', 'ne', 'xr', 'tr', 'thr', 'P', 'Ge', 'Gi', 'phir', 'yr']
    iq = 0

    def __init__(self, N):

        #d = defaultdict(list)
        self.d = { }

        for p in self.params:
            self.d[p] = np.zeros(N)


    def append(self, p):
        self.d['rho'][self.iq] = p.rho
        self.d['ne'][self.iq]  = p.ne
        self.d['xr'][self.iq]  = p.xr
        self.d['tr'][self.iq]  = p.tr
        self.d['thr'][self.iq] = p.thetar

        self.d['Ge'][self.iq] = p.Ge
        self.d['Gi'][self.iq] = p.Gi

        self.d['P'][self.iq]   = p.pressure()

        self.d['phir'][self.iq] = p.phir
    

        self.iq += 1


    #find density where param attains val
    def rho(self, param, val):
        xx = self.d['rho']
        yy = self.d[param]

        #Find roots using spline interpolation
        f = interpolate.UnivariateSpline(xx, yy - val, s=0)
        roots = f.roots()

        if len(roots) < 1:
            roots = [0.0]
        return roots[0]



#polytropic pressure relation using the adiabatic index
# for degenerate electron gases
def Ppoly(rho, gad):
    ne = rho/mp
    pFe = hbar*(3*np.pi**2*ne)**(1/3)
    xr = pFe / me / c

    return (Pr/9/gad/np.pi**2) * xr**(3*gad)
    
#polytropic for neutron gas
def EoSpoly(rho, gad):
    nn = rho/mn
    pFn = hbar*(3*np.pi**2*nn)**(1/3)
    yr = pFn / mn / c

    return (Prn/9/gad/np.pi**2) * yr**(3*gad)





#--------------------------------------------------
def main(argv):

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=9)
    

    #fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
    fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
    gs = plt.GridSpec(1, 1)

    N = 100
    con = container(N)

    axs = []
    for i in range(1):
        ax = plt.subplot(gs[i, 0])
        ax.minorticks_on()
        #ax1.set_ylim(0.0, 1.1)
        ax.set_xlim(1.0e1, 1.0e17)
        ax.set_ylim(1.0e14, 1.0e42)
        ax.set_xscale('log')
        ax.set_yscale('log')

        axs.append( ax )

    axs[0].set_xlabel(r'Density $\rho$ (g cm$^{-3}$)')
    axs[0].set_ylabel(r'Pressure $P$ (dyne cm$^{-2}$)')


    # basic parameters
    rhos = np.logspace(-2, 14.2, N)
    T = 1.0 * kelvin_per_keV #about 10 million kelvins
    Z = 1.0 #hydrogen plasma


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


    for i, rho in enumerate( rhos ):

        p = plasma(rho, T, Z)
        #print(p)

        con.append(p)

        
    axs[0].plot(con.d['rho'], con.d['P'], 'b-')

    rho2 = np.logspace(3, 9, 100)
    Pnonrel = Ppoly(rho2, 5/3)
    Prel = Ppoly(rho2, 4/3)
    axs[0].plot(rho2, Pnonrel, "k", linestyle='dotted')
    axs[0].plot(rho2, Prel, "k", linestyle='dashed')

    #Core
    rho3 = np.logspace(14.2, 16, 100)
    for i, eoss in enumerate(['PAL6', 'SLy', 'FPS', 'APR1', 'WFF1', 'BBB2', 'ENG', 'MS1', 'PS', 'H1', 'ALF1']):
        P1 = EoS(eoss, rho3)
        axs[0].plot(rho3, P1, "k-")



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
            """.format(rho_xr, rho_tr, rho_thr, rho, rho_phir, rho_Ge, rho_Gi)

    #highlight different regions
    xmin = 1.0
    xmax = 1.0e4
    axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e42, 1e42], facecolor='grey', color=None, alpha=0.1, edgecolor=None)

    #outer crust
    #xmin = 1.0e4
    #xmax = rhond
    #axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e40, 1e40], facecolor='grey', color=None, alpha=0.1, edgecolor=None)

    #inner crust
    xmin = rhond
    xmax = 0.5*rhon
    axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e42, 1e42], facecolor='grey', color=None, alpha=0.1, edgecolor=None)

    #outer core
    #xmin = 0.5*rhon
    #xmax = 2.0*rhon
    #axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e40, 1e40], facecolor='grey', color=None, alpha=0.1, edgecolor=None)

    #inner core
    #xmin = 2.0*rhon
    #xmax = 100.0*rhon
    #axs[0].fill_between([xmin, xmax], [1e10, 1e10], [1e40, 1e40], facecolor='grey', color=None, alpha=0.1, edgecolor=None)

    #different ns structures
    y_text = 1.0e40
    axs[0].text(400.0,  y_text, 'Atmosphere', rotation=0, ha='center', va='center', size=10)
    axs[0].text(1.0e8,  y_text, 'Outer crust', rotation=0, ha='center', va='center', size=10)
    axs[0].text(6.0e12, y_text, 'Inner\ncrust', rotation=0, ha='center', va='center', size=10)
    axs[0].text(6.0e15, y_text, 'Core', rotation=0, ha='center', va='center', size=10)

    
    lstyle = 'dotted'
    axs[0].plot([rhon, rhon], [1.0e16, 1e40], "k", linestyle=lstyle)
    txt = axs[0].text(rhon, 2.0e24, r'$\rho_n$', rotation=90, ha='center', va='center', size=8)
    txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))

    axs[0].plot([rhond, rhond], [1.0e16, 1e40], "k", linestyle=lstyle)
    txt = axs[0].text(rhond, 2.0e24, r'$\rho_{\mathrm{ND}}$', rotation=90, ha='center', va='center', size=8)
    txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))

    axs[0].plot([rho_xr, rho_xr], [1.0e16, 1e38], "k", linestyle=lstyle)
    txt = axs[0].text(rho_xr, 2.0e28, r'$x_r = 1$', rotation=90, ha='center', va='center', size=8)
    txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))

    axs[0].plot([rho_thr, rho_thr], [1.0e16, 1e28], "k", linestyle=lstyle)
    txt = axs[0].text(rho_thr, 2.0e22, r'$\Theta_r = 1$', rotation=90, ha='center', va='center', size=8)
    txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))

    axs[0].plot([rho_Ge, rho_Ge], [1.0e16, 1e30], "k", linestyle=lstyle)
    txt = axs[0].text(rho_Ge, 2.0e26, r'$\Gamma_{\mathrm{e}} = 1$', rotation=90, ha='center', va='center', size=8)
    txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))


    axs[0].text(1.0e4, 2.0e18, r'$P \propto \rho^{5/3}$', rotation=22, ha='center', va='center', size=8)
    axs[0].text(1.0e7, 2.0e25, r'$P \propto \rho^{4/3}$', rotation=27, ha='center', va='center', size=8)



if __name__ == "__main__":
    main(sys.argv)

    plt.savefig('eos.pdf')


