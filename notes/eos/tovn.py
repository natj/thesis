from __future__ import division

import sys
import os

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

print "Matplotlib version", matplotlib.__version__

from matplotlib import cm
import palettable as pal

import numpy as np
from scipy import interpolate

from units import *
from eoslibrary import EoS
from eoslibrary import eosLib




#cmap = pal.cmocean.sequential.Matter_8.mpl_colormap #best so far
#cmap = pal.wesanderson.Zissou_5.mpl_colormap
cmap = pal.colorbrewer.qualitative.Set1_6.mpl_colormap
#cmap = plt.get_cmap('plasma_r')
#cmap = cm.get_cmap('inferno_r')

from scipy.integrate import odeint
from math import pi




#--------------------------------------------------
class monotrope:
    
    #transition continuity constant 
    a = 0.0
    cgsunits = c*c

    def __init__(self, K, G):
        self.K = K / self.cgsunits
        self.G = G
        self.n = 1.0/(G - 1)

    #pressure P(rho)
    def pressure(self, rho):
        return self.cgsunits * self.K * rho**self.G

    #energy density mu(rho)
    def edens(self, rho):
        return (1.0 + self.a)*rho + (self.K/(self.G - 1)) * rho**self.G

    #for inverse functions lets define rho(P)
    def rho(self, press):
        if press < 0.0:
            return 0.0
        return ( press/self.cgsunits/self.K )**(1 / self.G)

    #internal energy
    def intenerg(self, press):
        return press/(self.G - 1)/self.rho(press)

    




    #in terms of entalphy
    ##################################################

    #rho(eta)
    def eta_rho(self, eta):
        return ((eta/self.cgsunits - self.a)/(self.K*(self.n + 1)))**self.n

    #pressure(eta)
    def eta_pressure(self, eta):
        return self.cgsunits * self.K*( (eta/self.cgsunits - self.a)/(self.K*(self.n +1)) )**(self.n+1)

    #edens(eta)
    def eta_edens(self, eta):
        t1 = (self.a + self.n*eta/self.cgsunits )/(self.n + 1)
        return self.eta_rho( eta )*(1.0 + t1)
    
    #eta(rho)
    def rho_eta(self, rho):
        return self.a + self.K*(self.n + 1)*rho**(1.0/self.n)




class polytrope:
    
    def __init__(self, tropes, trans, prev_trope = None ):
        self.tropes      = tropes
        self.transitions = trans

        self.prs  = []
        self.eds  = []
        self.etas = []

        #if prev_trope == (None, None):
        #    prev_ed = 0.0
        #    prev_tr = 0.0
        #else:
        #    prev_ed = prev_trope[0]
        #    prev_tr = prev_trope[1]
        
        
        for (trope, transition) in zip(self.tropes, self.transitions):

            if not( prev_trope == None ):
                #trope.a = self._ai( prev_ed, prev_tr, trope.K, trope.G )
                trope.a = self._ai( prev_trope, trope, transition )
                print "continuity constant a=", trope.a
            else:
                print "no continuity constant a=", trope.a
                transition = 0.0


            ed = trope.edens(transition) 
            pr = trope.pressure(transition)
            eta= trope.rho_eta(transition)

            self.prs.append( pr )
            self.eds.append( ed )
            self.etas.append(eta)

            prev_ed = ed
            prev_tr = transition
            #first_done = True
            prev_trope = trope


    #def _ai(self, ed_im1, rho_im1, K, G):
    #    return ed_im1/rho_im1 - 1- (K/(G-1))*rho_im1**(G - 1)

    def _ai(self, pm, m, tr):
        return pm.a + (pm.K/(pm.G - 1))*tr**(pm.G-1) - (m.K/(m.G - 1))*tr**(m.G-1)



    def _find_interval_given_density(self, rho):
        if rho <= self.transitions[0]:
            return self.tropes[0]

        for q in range( len(self.transitions) - 1 ):
            if self.transitions[q] <= rho < self.transitions[q+1]:
                return self.tropes[q]

        return self.tropes[-1]


    #inverted equations as a function of pressure
    def _find_interval_given_pressure(self, press):
        if press <= self.prs[0]:
            return self.tropes[0]

        for q in range( len(self.prs) - 1):
            if self.prs[q] <= press < self.prs[q+1]:
                return self.tropes[q]

        return self.tropes[-1]

    #inverted equations as a function of enthalpy
    def _find_interval_given_enthalpy(self, eta):
        if eta <= self.etas[0]:
            return self.tropes[0]

        for q in range( len(self.etas) - 1):
        #for q in range(1, len(self.etas)):
            #if self.etas[q-1] < eta < self.etas[q]:
            if self.etas[q] <= eta < self.etas[q+1]:
                return self.tropes[q]

        return self.tropes[-1]


    ################################################## 
    def pressure(self, rho):
        trope = self._find_interval_given_density(rho)
        return trope.pressure(rho)

    #vectorized version
    def pressures(self, rhos):
        press = []
        for rho in rhos:
            pr = self.pressure(rho)
            press.append( pr )
        return press

    #def edens_given_density(self, rho):
    #    trope = self._find_interval_given_density(rho)
    #    return trope.edens(rho)

    #def edens_given_densities(self, rhos):
    #    eden = []
    #    for rho in rhos:
    #        ed = self.edens_given_density(rho)
    #        eden.append( ed ) 
    #    return eden

    def edens_inv(self, press):
        trope = self._find_interval_given_pressure(press)
        rho = trope.rho(press)
        return trope.edens(rho)

    def edenses(self, pressures):
        eden = []
        for press in pressures:
            ed = self.edens_inv(press)
            eden.append( ed ) 
        return eden

    def rho(self, press):
        trope = self._find_interval_given_pressure(press)
        return trope.rho(press)


    ################################################## 
    def eta_edens(self, eta):
        trope = self._find_interval_given_enthalpy(eta)
        return trope.eta_edens(eta)

    def eta_pressure(self, eta):
        trope = self._find_interval_given_enthalpy(eta)
        return trope.eta_pressure(eta)

    def eta_rho(self, eta):
        trope = self._find_interval_given_enthalpy(eta)
        return trope.eta_rho(eta)

    def rho_eta(self, rho):
        trope = self._find_interval_given_density(rho)
        return trope.rho_eta(rho)



##################################################
#SLy (Skyrme) crust
KSLy = [6.80110e-9, 1.06186e-6, 5.32697e1, 3.99874e-8] #Scaling constants
GSLy = [1.58425, 1.28733, 0.62223, 1.35692] #polytropic indices
RSLy = [1.e4, 2.44034e7, 3.78358e11, 2.62780e12 ] #transition depths

#KSLy = [6.80110e-9, 1.06186e-6, 5.32697e1, 3.99874e-8] #Scaling constants
#GSLy = [1.6, 1.3, 1.1, 1.4] #polytropic indices
#RSLy = [2.4e7, 3.78358e11, 2.62780e12, 1.4e14] #transition depths


tropes = []
trans = []

pm = None
for (K, G, r) in zip(KSLy, GSLy, RSLy):
    m = monotrope(K*c*c, G)
    #m = monotrope(K, G)
    tropes.append( m )

    #correct transition depths to avoid jumps
    if not(pm == None):
        rho_tr = (m.K / pm.K )**( 1.0/( pm.G - m.G ) )
        #print rho_tr, np.log10(rho_tr), r, rho_tr/r
    else:
        rho_tr = r
    pm = m

    #trans.append(r)
    trans.append(rho_tr)
SLyCrust = polytrope(tropes, trans)


#EoS class using dense eos from Read et al (2009) 
def read_eos(key):

    #read eos table for parameters
    ll = eosLib[ key ]

    #dense eos starting pressure
    p1 = 10**ll[0]

    #polytrope indices
    g1 = ll[1] 
    g2 = ll[2]
    g3 = ll[3]

    #transition densities
    r1 = 2.8e14
    r2 = 10**14.7
    r3 = 10**15.0

    #scaling constants
    #K1 = p1/(r2**g1)
    #K2 = (K1 * r2**g1)/r2**g2
    #K3 = (K2 * r3**g2)/r3**g3

    K1 = p1/(r2**g1)
    K2 = K1 * r2**(g1-g2)
    K3 = K2 * r3**(g2-g3)

    cgsunits = 1.0/c*c
    tropes = [monotrope(K1*cgsunits, g1),
              monotrope(K2*cgsunits, g2),
              monotrope(K3*cgsunits, g3) ]
    trans = [r1, r2, r3]

    dense_eos = polytrope(tropes, trans)

    return dense_eos


# Smoothly glue core to SLy crust
# for polytropic eos we can just unpack 
# and repack the piecewise presentation
def crust_and_core(crust, core):

    #unpack crust and core
    tropes_crust = crust.tropes
    trans_crust  = crust.transitions

    tropes_core = core.tropes
    trans_core  = core.transitions

    #find transition depth
    rho_tr = (tropes_core[0].K / tropes_crust[-1].K )**( 1.0/( tropes_crust[-1].G - tropes_core[0].G ) )
    print "Transition from core to crust at", rho_tr, np.log10(rho_tr), crust.edens_inv( crust.pressure( rho_tr ) )/GeVfm_per_dynecm
    trans_core[0] = rho_tr
    #trans_crust[-1] = rho_tr

    #repack
    tropes = tropes_crust + tropes_core
    trans  = trans_crust  + trans_core

    for trope in tropes:
        trope.a = 0.0

    eos = polytrope( tropes, trans )
    return eos


def simple_eos():

    print "building simple eos"
    K1 = 3.99873692e-8 
    G1 = 1.35692395 
    r1 = 0

    K2 = 2.23872092e-10
    G2 = 3 
    #a = 0.010350691 * c *c
    r2 = 1.4172900e14

    m1 = monotrope(K1*c*c, G1)
    m2 = monotrope(K2, G2)
    eos = polytrope([m1, m2], [r1, r2])

    print "a0", eos.tropes[0].a
    print "a1", eos.tropes[1].a, eos.tropes[1].a

    #eos.tropes[1].a = 0.010350691

    return eos

def simple_eos2():
    
    hbar = 1.05457266e-27
    mn   = 1.6749286e-24

    Gamma0 = 5.0/3.0 
    K0 = (3.0*pi**2)**(2.0/3.0)*hbar**2/(5.0*mn**(8.0/3.0))
    Gamma1 = 3.0 # high densities: stiffer equation of state
    #Gamma1 = 2.5 # high densities: softer equation of state
    rho1 = 5e14
    P1 = K0*rho1**Gamma0
    K1 = P1/rho1**Gamma1

    m1 = monotrope(K0, Gamma0)
    m2 = monotrope(K1, Gamma1)
    eos = polytrope([m1, m2], [0.0, rho1])

    return eos


def simpler_eos():
    K = 1.982e-6
    G = 2.75
    m = monotrope(K, G)
    eos = polytrope([m], [0.0])
    return eos




class tov:

    def __init__(self, peos):
        self.physical_eos = peos

    def tov(self, y, r):
        P, m = y
        #rho = self.physical_eos.rho( P )
        eden = self.physical_eos.edens_inv( P )

        c  = 2.99792458e10
        G  = 6.6730831e-8

        #dPdr = -G*(rho + P/c**2)*(m + 4.0*pi*r**3*P/c**2)
        #dPdr = dPdr/(r*(r - 2.0*G*m/c**2))
        #dmdr = 4.0*pi*r**2*rho

        #with energy density
        dPdr = -G*(eden + P/c**2)*(m + 4.0*pi*r**3*P/c**2)
        dPdr = dPdr/(r*(r - 2.0*G*m/c**2))
        dmdr = 4.0*pi*r**2*eden

        return [dPdr, dmdr]

    def tovsolve(self, rhoc):

        N = 500
        r = np.linspace(1e0, 18e5, N)
        P = self.physical_eos.pressure( rhoc )
        eden = self.physical_eos.edens_inv( P )
        #m = 4.0*pi*r[0]**3*rhoc
        m = 4.0*pi*r[0]**3*eden

        psol = odeint(self.tov, [P, m], r, rtol=1.0e-4, atol=1.0e-4)
        return r, psol[:,0], psol[:,1]

    def mass_radius(self):
        N = 100
        mcurve = np.zeros(N)
        rcurve = np.zeros(N)
        rhocs = np.logspace(14.0, 16.0, N)

        mass_max = 0.0
        j = 0
        for rhoc in rhocs:
            rad, press, mass = self.tovsolve(rhoc)

            rad /= 1.0e5
            mass/= Msun

            mstar = mass[-1]
            rstar = rad[-1]
            for i, p in enumerate(press):
                if p > 0.0:
                    mstar = mass[i]
                    rstar = rad[i]
            mcurve[j] = mstar
            rcurve[j] = rstar

            j += 1
            if mass_max < mstar:
                mass_max = mstar
            else:
                break

        return mcurve[:j], rcurve[:j], rhocs[:j]


    def m_given_r(self, r):
        p = self.pr(18.0 - r)
        self.rp.roots

        return self.mp(p)

    def r_given_m(self, m):
        p = self.pm(m)
        return self.rp(p)

    def mr_curve(self):
        mcurve, rcurve, rhocs = self.mass_radius()
        lrhoc = np.log10(rhocs)
        
        print mcurve
        print rcurve
        print lrhoc

        #self.pm = interpolate.UnivariateSpline(mcurve, lrhoc, s=0)
        #self.pr = interpolate.UnivariateSpline(18.0 - rcurve, lrhoc, s=0)

        self.mp = interpolate.UnivariateSpline(lrhoc, mcurve, s=0)
        self.rp = interpolate.UnivariateSpline(lrhoc, rcurve, s=0)

        #print "R_1.4: {}".format( self.r_given_m(1.4) )



    def eos_eta(self, eta):

        press = self.physical_eos.eta_pressure(eta)
        eden  = self.physical_eos.eta_edens(eta)

        #wrap physical eos to numerical units
        #pres *= GeVfm_per_dynecm
        #ed = self.physical_eos.edens( pres )
        #return ed/GeVfm_per_dynecm
        #press *= GeVfm_per_dynecm
        #eden  *= GeVfm_per_dynecm

        return press, eden


    def eos_dim(self, pres):

        if pres<=0: 
            pres=1.e-12
        #wrap physical eos to numerical units
        pres *= GeVfm_per_dynecm
        ed = self.physical_eos.edens( pres )

        return ed/GeVfm_per_dynecm


    def eos_simple(self, pres):
        # this is the fit to the non-interacting neutron matter pressure 
        anr = 3.9294
        ar = 2.86663 

        # this is proton, neutron, and electron matter in beta-equilibrium
        #anr = 4.27675893
        #ar = 2.84265221

        if pres<=0: 
            pres=1.e-12

        # results are in GeV/fm^3
        return (anr*pres**(3./5.) + ar*pres) #*GeVfm_per_dynecm


    #G=c=Msun=1 dimensionsionality
    def tovrhs_dim(self, initial,x):
        pres=initial[0]
        mass=initial[1]
        alfa = 41.325   # this is the ratio of (M_solar/R_s^3) in units of GeV/fm^3
        eden = self.eos_dim(pres) 
        one = -0.5*eden*mass*(1+(pres/eden))*(1 + (4*pi/alfa)*(pres/mass)*x**3)/(x**2-mass*x)
        two = 4*pi*x**2*eden/alfa
        f = [one,two]
        return f 

    #eos in cgs units
    def tovrhs_cgs(self, initial, rad):
        press = initial[0]
        mass  = initial[1]

        eden = self.physical_eos.edens( press )
        rho  = self.physical_eos.rho( press )
        
        G  = 6.6730831e-8
        c  = 2.99792458e10

        one = -G*(rho*(1 + eden/c/c) + press/c/c)*((mass + 4*pi*rad**3*press/c/c)/(rad*(rad - 2*G*mass/c/c)))
        two = 4*pi*rad**2*rho*(1+eden/c/c)

        return [one, two]


    def tovrhs_eta(self, initial, eta):
        rad  = initial[0]
        mass = initial[1]


        press, eden = self.eos_eta(eta)
        print "eta: {}    rad: {} mass: {} | press: {} eden: {}".format(eta, rad, mass, press, eden)

        #alfa = 41.325 #Msun/Rs^3 to GeV/fm^3

        #transform to dimensionless units
        G  = 6.6730831e-8

        #dr/d\eta
        #one = - rad * (rad - 2*mass)/( mass + 4*pi*rad**3*press )/(eta + 1)

        relcorr = (1-2*G*mass/rad/c/c)
        masst   = 1 + (4*pi*rad**3*press/mass/c/c)
        etat    = 1/(eta + c*c)
        prefac  = c*c*rad**2/mass/G
        one     = - prefac * relcorr * etat / masst
        print "relcorr: {} masst: {} etat: {} prefac: {} = {}".format(relcorr, masst, etat, prefac, one)

        #dm/d\eta
        two = 4*pi*(rad**2)*eden*one

        #print "one:", one
        #print "two:", two

        return [one, two]



    def tsolve_dim(self, pcent, xfinal):
        alfa = 41.325   # this is the ratio of (M_solar/R_s^3) in units of GeV/fm^3
        #alfa = 41.257272820415594
        eden = self.eos_dim(pcent)
        dx = 1.0e-3
        initial = pcent, 4*pi*eden*dx**3/(3.*alfa) 
        x= np.arange(dx,xfinal,dx)

        psol = odeint(self.tovrhs_dim,initial,x )#, printmessg=True, rtol=1.e-14, atol=1e-10 ,hmax=1.e-4 )
        rstar = 0. 
        mstar = 0.
        count = 0
        for i in psol[:,0]:
            if i>1.e-6:
            #if i>0:
                #print i, count, psol[count, 1]
                rstar = rstar + dx*2.95
                mstar = psol[count,1]
                count = count + 1    
        return rstar,mstar


    def mass_radius_dim(self, pmin, pmax): 
        imax = 30
        pc = np.zeros(imax)
        mass   = np.zeros(imax)
        radius = np.zeros(imax)

        #for i, p in enumerate( np.linspace(pmin, pmax, imax) ):
        for i, p in enumerate( np.logspace( np.log10(pmin), np.log10(pmax), imax) ):
            pc[i] = p
            radius[i],mass[i] = self.tsolve(pc[i], 50.0/2.95 )
        return pc,radius,mass


    def tsolve_dim_psol(self, pcent):
        xfinal = 50.0/2.95

        alfa = 41.325   # this is the ratio of (M_solar/R_s^3) in units of GeV/fm^3
        eden = self.eos_dim(pcent)
        dx = 1.0e-3
        initial = pcent, 4*pi*eden*dx**3/(3.*alfa) 
        x= np.arange(dx,xfinal,dx)

        psol = odeint(self.tovrhs_dim, initial, x) #, tcrit=[0.0] )#, printmessg=True, rtol=1.e-14, atol=1e-10 ,hmax=1.e-4 )
        return x, psol

    def tsolve_cgs(self, press_cent):
        x = np.linspace(1.0e-10, 20.0, 100)*1.0e5

        eden = self.eos_cgs(press_cent)
        initial = press_cent, 4*pi*eden*dx**3/3

        psol = odeint(self.tovrhs_cgs, initial, x, tcrit=[0.0])
        #print psol[:,0]
        #print psol[:,1]

        return x, psol


    def tsolve_eta(self, rhoc):
        etac = self.physical_eos.rho_eta( rhoc )
        print "rhoc: {} and eta: {}".format(rhoc, etac)

        x = np.logspace(np.log10(etac), -2, 20)

        #x = np.linspace(etac, 0.0, 20)
        print "x:", x

        initial = 1e-8, 1e-8
        #initial = 0.0, 0.0

        psol = odeint(self.tovrhs_eta, initial, x)

        return x, psol



#--------------------------------------------------
def main(argv):

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=7)
    

    #fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
    fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
    gs = plt.GridSpec(1, 2)

    ax = plt.subplot(gs[0, 0])
    ax.minorticks_on()
    #ax.set_xlim(8.0, 18.0)
    #ax.set_ylim(0.0, 2.5)

    ax.set_xlabel(r'Radius $R$ (km)')
    ax.set_ylabel(r'Mass $M$ (M$_{\odot}$)')

    #extra axis
    ax2 = plt.subplot(gs[0, 1])
    ax2.minorticks_on()
    ax2.set_xlabel(r'$\epsilon$ (GeV fm$^{-3}$)')
    ax2.set_ylabel(r'$P$ (GeV fm$^{-3}$)')
    ax2.set_yscale('log')
    ax2.set_xscale('log')



    key = 'SLy'
    #key = 'BGN1H1'
    #key = 'ALF2'
    #key = 'ENG'
    #key = 'PAL6'
    rho = np.logspace(4.0, 16, 100)

    dense_eos = read_eos(key)
    eos = crust_and_core( SLyCrust, dense_eos )
    #eos = simple_eos()
    #eos = simpler_eos()
    #eos = dense_eos
    #eos = SLyCrust

    #P = np.array( eos.pressures( rho ) )
    #eds = np.array( eos.edenses(P) )
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.plot( eds/GeVfm_per_dynecm, P/GeVfm_per_dynecm )

    t = tov(eos)
    #t = tov(dense_eos)
    #t = tov( simple_eos2() )





    #XXX extra plot
    if False:
        rhos = np.logspace(4, 16, 100)
        y1 = []
        y2 = []
        eos2 = simple_eos2()
        for rho in rhos:
            y1.append( eos.pressure(rho) )
            y2.append( eos2.pressure(rho) )
        ax2.plot(rhos, y1, "b-")
        ax2.plot(rhos, y2, "r--")

        press = np.logspace(19, 38, 100)
        x1 = []
        for pres in press:
            x1.append( eos.rho( pres) )
        ax2.plot(x1, press, "k--")
    if False:
        rhos = np.logspace(-4, 16, 100)
        press = []
        for rho in rhos:
            press.append( eos.pressure(rho) )

        ax2.plot(rhos, press)
    if False:
        edt = []
        edt2 = []
        prg = np.logspace(-10, 2, 1000)

        for press in prg:
            ed = t.eos_dim(press)
            ed2= t.eos_dim(press)
            edt.append( ed )
            edt2.append( ed2 )
        ax2.plot( edt, prg, "b-")
        ax2.plot( edt2, prg, "r--")
    if False: #eta - P
        xx = []
        xx2 = []
        #yy = np.logspace(-6, 2, 1000) #eta
        yy = np.logspace(2, 22, 1000) #eta in cc

        #eos2 = simpler_eos()
        eos2 = simple_eos()
        for y in yy:
            x = eos.eta_rho(y) #rho(eta)
            #x = eos.eta_pressure(y) #p(eta)
            #x = eos.eta_rho(y)
            xx.append( x )
            xx2.append( eos2.eta_rho(y) )
            #xx2.append( eos2.eta_pressure(y) )

        ax2.plot(yy, xx)
        ax2.plot(yy, xx2, "r--")

    if True: #eden - P
        xx = []
        xx2 = []
        yy = []
        yy2 = []
        rhos = np.logspace(-2, 17, 1000) #eta
        eos2 = simpler_eos()
        #eos2 = simple_eos()
        for rho in rhos:
            x = eos.edens_inv( rho )
            y = eos.pressure( rho )
            xx.append( x )
            yy.append( y )

            x = eos2.edens_inv( rho )
            y = eos2.pressure( rho )
            xx2.append( x )
            yy2.append( y )

        #ax2.plot(rhos, yy, "b-")
        #ax2.plot(rhos, yy2, "r-")
        ax2.plot(xx, yy, "b-")
        ax2.plot(xx2, yy2, "r-")

    
    if False:
        #x, psol = t.tsolve_eta(1.0e16)
        x, psol = t.tsolve_dim_psol(0.8)
        #x, psol = t.mass_radius_uno(1.0e37)
        rad = x
        press = psol[:,0]
        mass  = psol[:,1] #/Msun
        #ax.plot(x, psol[:,0], "b-")
        #ax.plot(x, psol[:,1], "g-")
        #ax.plot(rad, press)
        ax.plot(rad, mass)

    if False:
        rhoc = 1.0e16/c/c
        etac = t.physical_eos.rho_eta( rhoc )
        print "rhoc: {} and eta: {}".format(rhoc, etac)
        xx = np.logspace(np.log10(etac), -2, 20)

        yy1 = []
        yy2 = []

        initial = 1e-10, 1e-10
        for x in xx:
            print "initial: ", initial
            print " eta: ", x
            one, two = t.tovrhs_eta(initial, x)
            yy1.append(one)
            yy2.append(two)
        ax.plot(xx, yy1)

    if False:
        rhoc = 0.9e15
        rad, press, mass = t.tovsolve(rhoc)
        rad  /= 1.0e5
        mass /= 1.988435e33
        #ax.plot(rad, press)
        #ax.set_yscale('log')
        ax.plot(rad, mass)
    if False:
        mass, rad, rho = t.mass_radius()
        ax.plot(rad, mass)



    #pc, rad, mass = t.mass_radius(1.0e-4, 0.5)
    #pc, rad, mass = t.mass_radius(1.0e-6, 0.5)
    #ax.plot(rad/2.95, mass/2.95, "b-")
    #ax.plot(rad, mass, "b-") #XXX
    

    #for key, value in eosLib.iteritems():
    #    P1 = EoS(key, rho3)
    #    if value[4] == 'npem':
    #        axs[0].plot(rho3, P1, "k-", alpha = 0.3)
    #    if value[4] == 'meson':
    #        axs[0].plot(rho3, P1, "b-", alpha = 0.3)
    #    if value[4] == 'hyperon':
    #        axs[0].plot(rho3, P1, "g-", alpha = 0.3)
    #    if value[4] == 'quark':
    #        axs[0].plot(rho3, P1, "r-", alpha = 0.3)

    if True:
        for key, value in eosLib.iteritems():

            dense_eos = read_eos(key)
            eos = crust_and_core( SLyCrust, dense_eos )
            t = tov(eos)
            mass, rad, rhoc = t.mass_radius()

            t.mr_curve()

            if value[4] == 'npem':
                ax.plot(rad, mass, "k-", alpha = 0.9)
            if value[4] == 'meson':
                ax.plot(rad, mass, "b-", alpha = 0.9)
            if value[4] == 'hyperon':
                ax.plot(rad, mass, "g-", alpha = 0.9)
            if value[4] == 'quark':
                ax.plot(rad, mass, "r-", alpha = 0.9)



if __name__ == "__main__":
    main(sys.argv)
    plt.subplots_adjust(left=0.15, bottom=0.16, right=0.98, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig('mr_edens.pdf')



