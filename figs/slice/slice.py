import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#import units as cgs
from math import pi
#from polytropes import monotrope, polytrope
#from crust import SLyCrust
#from eoslib import get_eos, glue_crust_and_core, eosLib
#from scipy.integrate import odeint
#from label_line import label_line

import math

import matplotlib.patches as patches
from scipy.optimize import minimize_scalar

from pylab import Polygon
from curved_text import CurvedText

#from matplotlib.path import Path
#from matplotlib.patches import PathPatch
#from matplotlib import cm

import palettable as pal
cmap = pal.colorbrewer.qualitative.Set1_6.mpl_colormap

pi = np.pi


#--------------------------------------------------
def arc(r, a1=0.0, a2=pi/3):

    N = 100
    x = np.zeros(N) 
    y = np.zeros(N) 
    for i, a in enumerate(np.linspace(a1, a2, N)):
        x[i] = r*np.sin(a) 
        y[i] = r*np.cos(a) 

    return x, y

def slice(r1, r2, a1=0.0, a2=pi/3):

    if r1 == 0.0:
        x1 = np.array([0.0, 0.0])
        y1 = np.array([0.0, 0.0])
    else:
        x1, y1 = arc(r1, a1=a1, a2=a2)

    x2, y2 = arc(r2, a1=a1, a2=a2)

    #combine (numpy is stupid about array + scalar merging so we do it in parts)
    x = np.append( x1[0:1], x2 )
    x = np.append(x, np.flipud(x1) )
    y = np.append( y1[0:1], y2 )
    y = np.append(y, np.flipud(y1) )

    return x, y

def fill_slice(ax, r1, r2, fmt):
    xx, yy = slice(r1, r2)
    ax.add_patch(Polygon(zip(xx,yy),
                      closed=True, 
                      **fmt))

    return None


def fill_box(ax, r1, r2, fmt):
    xmax = 1.0e36
    xmin = 1.0e-2
    xx = [xmax, xmin, xmin, xmax]
    yy = [r2, r2, r1, r1]

    ax.add_patch(Polygon(zip(xx,yy),
                      closed=True, 
                      **fmt))



plt.rc('font', family='serif')
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)
plt.rc('axes', labelsize=7)


#fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
gs = plt.GridSpec(1, 2, wspace=0.0)

ax = plt.subplot(gs[0, 1])
ax.minorticks_on()
ax.set_xlim(0.0, 12.0)
ax.set_ylim(0.0, 12.0)

ax.set_ylabel(r'Radius $R$ (km)')
ax.spines['top'].set_visible(None)
ax.spines['right'].set_visible(None)
ax.spines['bottom'].set_visible(None)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.axes.get_xaxis().set_visible(False)


ax.set_yticklabels([])
ax.set_xticklabels([])



#-------------------------------------------------- 

#x, y = arc(10.0)
#ax.plot( x, y, "k-")

rad_surf  = 11.73
rad_atmos = 11.72
rad_crust = 11.4
rad_core  = 11.0
rad_coreb = 6.0

#draw full star
x, y = slice(0.0, rad_surf)
ax.plot( x, y, "k-")

#atmos boundary
#x, y = arc(rad_atmos)
#ax.plot( x, y, "k--")
x, y = arc(rad_surf)
ax.plot( x, y, "r-", linewidth=1.5, alpha=0.7)
x, y = arc(11.8)
text = CurvedText(x, y, text='           Atmosphere', va='bottom', axes=ax)

#crust boundary
x, y = arc(rad_crust)
#ax.plot( x, y, "k--")

x, y = arc(rad_core)
text = CurvedText(x, y, text='            Crust', va='bottom', axes=ax)

#core boundary
x, y = arc(rad_core)
ax.plot( x, y, "k-")

x, y = arc(rad_core-1.0)
text = CurvedText(x, y, text='        Outer core', va='top', axes=ax)
#x, y = arc(rad_core-2.0)
#text = CurvedText(x, y, text='      n p e m', va='top', axes=ax, size=10)

#inner core boundary
x, y = arc(rad_coreb)
ax.plot( x, y, "k--")

x, y = arc(rad_coreb-1.0)
text = CurvedText(x, y, text=' Inner core', va='top', axes=ax)

#annotate densities
#labs = [r'                  $(3-9) \rho_0$', r'  $\approx 2\rho_0$', r'$\approx 0.5 \rho_0$']
labs = [r'            $\approx 3\rho_0$', r'  $\approx 2\rho_0$', r'$\approx 0.5 \rho_0$']
rads = [0.0, 6.0, 10.5]
for rad, lab in zip(rads, labs):
    x, y = arc(rad, a1=0.0, a2=pi/2.75)
    ax.text(x[-1], y[-1], lab, rotation=-50.0, size=8, va='center', ha='center')



#fill_slice(ax, rad_crust, rad_surf, fmt={'facecolor':'k', 'alpha':0.1,})
#fill_slice(ax, 6.0, rad_core, fmt={'facecolor':'k', 'alpha':0.1,})

fill_slice(ax, rad_crust, rad_surf, fmt={'facecolor':'k', 'alpha':0.1,})
fill_slice(ax, rad_core, rad_crust, fmt={'facecolor':'k', 'alpha':0.15,})
fill_slice(ax, rad_coreb, rad_core, fmt={'facecolor':'k', 'alpha':0.3,})
fill_slice(ax, 0.0, rad_coreb, fmt={'facecolor':'k', 'alpha':0.4,})





#--------------------------------------------------

ax2 = plt.subplot(gs[0, 0])
ax2.minorticks_on()
#ax2.set_xlim(1.0e-2, 1.0e10)
ax2.set_xlim(1.0e36, 1.0e24)
ax2.set_xscale('log')
ax2.set_ylim(0.0, 12.0)

#ax2.set_ylabel(r'Radius $R$ (km)')
ax2.set_xlabel(r'Pressure $P$ (dyne cm$^{-2}$)', color='blue')
ax2.spines['top'].set_visible(None)
ax2.spines['left'].set_visible(None)
ax2.yaxis.set_ticks_position('right')
ax2.xaxis.set_ticks_position('bottom')

#ticklab = ax2.yaxis.get_ticklabels()[0]
#trans = ticklab.get_transform()
#ax2.yaxis.set_label_coords(1.0e25, 6, transform=trans)

ax2.text(1.0e25, 6.0, r'Radius $R$ (km)', rotation=90, size=7, ha='center', va='center')

shadefac = 0.2
fill_box(ax2, rad_crust, rad_surf, fmt={'facecolor':'k', 'alpha':0.1*shadefac, 'edgecolor':None})
fill_box(ax2, rad_core, rad_crust, fmt={'facecolor':'k', 'alpha':0.2*shadefac, 'edgecolor':None})
fill_box(ax2, rad_coreb, rad_core, fmt={'facecolor':'k', 'alpha':0.3*shadefac, 'edgecolor':None})
fill_box(ax2, 0.0,      rad_coreb, fmt={'facecolor':'k', 'alpha':0.4*shadefac, 'edgecolor':None})

#ax2.tick_params(axis='y', direction='in')
ax2.tick_params(axis='y',which="major", direction="in", 
                 left=False, right=True, reset=False, labelright=False)

#place y axis ticklabels by hand
ylabs = ['','2','4','6','8','10','12']
for i, yloc in enumerate([0.0, 2, 4, 6, 8, 10, 12]):
    ax2.text(3e24, yloc, ylabs[i], size=7, va='center', ha='center')


# map R to rho & P
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tov'))
from crust import SLyCrust
from eoslib import get_eos, glue_crust_and_core, eosLib
from tov import tov



# other axis
ax2b = ax2.twiny()
ax2b.spines["top"].set_position(("outward", 15.0))
ax2b.minorticks_on()
ax2b.set_xlim(1.5, -0.07)
ax2b.set_ylim(0.0, 12.0)
ax2b.set_xlabel(r'Mass $M$ (M$_{\odot}$)', color='red')

#ax2b.spines['top'].set_visible(None)
ax2b.spines['left'].set_visible(None)
#ax2b.yaxis.set_ticks_position('right')
#ax2b.xaxis.set_ticks_position('bottom')



dense_eos = get_eos('SLy')
eos = glue_crust_and_core( SLyCrust, dense_eos )
t = tov(eos)
#mass, rad, rho = t.mass_radius()
#for i in range(len(mass)):
#    print mass[i], rad[i], rho[i]

#for key, value in eosLib.iteritems():
for key in ['SLy']:
    dense_eos = get_eos(key)
    eos = glue_crust_and_core( SLyCrust, dense_eos )
    t = tov(eos)
    mass, rad, rho = t.mass_radius()
    rhoc = 0.0
    for i in range(len(mass)):
        #print mass[i], rad[i], rho[i]
        if mass[i] >= 1.4:
            break
        else:
            rhoc = rho[i]

    print "central density=", rhoc
    rads, press, masses =  t.tovsolve(rhoc)
    rads /= 1.0e5
    ax2.plot(press, rads, "b-")

    masses /= 1.988e33
    ax2b.plot(masses, rads, "r-")

    for i in range(len(rads)):
        print "r: {} m: {} p: {}".format(rads[i], masses[i]/masses[-1], press[i])


# Remove the last ytick label
fig.show()
fig.canvas.draw()
labels = [tick.get_text() for tick in ax2.get_xticklabels()]
labels[2] = ''
ax2.set_xticklabels(labels)







plt.subplots_adjust(left=0.1, bottom=0.12, right=0.98, top=0.85, wspace=0.1, hspace=0.1)
plt.savefig('slice.pdf')

