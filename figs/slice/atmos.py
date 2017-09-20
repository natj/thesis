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

def fill_slice(ax, r1, r2, a1, a2, fmt):
    xx, yy = slice(r1, r2, a1=a1, a2=a2)
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
fig = plt.figure(figsize=(7.48, 3.5))  #two column figure
gs = plt.GridSpec(1, 1, wspace=0.0)

ax = plt.subplot(gs[0, 0])
ax.minorticks_on()
ax.set_xlim(-11, 12.0)
ax.set_ylim(0.0, 13.0)
ax.axis('off')

#ax.set_yticklabels([])
#ax.set_xticklabels([])


#-------------------------------------------------- 

p1 = -pi/5
p2 = pi/5

def pol2xy(r, ang):
    x = r*np.sin(ang)
    y = r*np.cos(ang)
    return (x,y)


rads = [20.0, 12.0, 10.5, 9.0, 7.5, 6.0, 0.0]
labels = ['               C O R O N A',
        '        A T M O S P H E R E',
        '           O C E A N',
        '  O U T E R  C R U S T',
        'I N N E R  C R U S T',
        '   C O R E',
        ]

lstyles = ['dotted', 'dashed', 'solid', 'dashed', 'solid', None]

#for i in range(len(rads)-1):
for i in range(1,len(rads)-1):
    alpha = 0.02 + i*0.1
    r1 = rads[i]
    r2 = rads[i+1]
    fill_slice(ax, r2, r1, p1, p2, fmt={'linewidth':0.0, 'facecolor':'k', 'alpha':alpha,})

    r = 0.5*(r1+r2)
    if i == 0:
        r = 12.6
    if i == 5:
        r = 0.5*(r1 + 4.0) 

    x, y = arc(r2, p1, p2)
    lstyle = lstyles[i]
    ax.plot(x, y, linestyle=lstyle, color='black', linewidth=0.5)

    x, y = arc(r, -pi/5, pi/5)

    label = labels[i]
    text = CurvedText(x, y, text=label, va='center', axes=ax)




textang = pi/5
x, y = pol2xy(15.0, textang)
#ax.text(x,y, r'$e^{-}$ $e^{+}$', size = 10)

xloc = 8.0

#ax.text(xloc, 10.0,  'electrons, ions,\n atoms, molecules \n (gas/liquid)', size = 9, ha='center', va='center')

strings = ['electrons and positrons \n (plasma)',
        'electrons, ions, \n atoms, molecules \n (gas/liquid)',
        'electrons and nuclei \n (Coulomb liquid)',
        'electrons and nuclei \n (Coulomb crystal)',
        'electrons, neutrons, \n neutron-rich nuclei \n (Coulomb crystal)',
        'neutron, electrons, \n protons \n (liquid)',
        ]

arrowang = pi/5.0
textang = pi/4
rads = [13.5, 11.3, 9.6, 8.0, 6.0, 3.2]
ylocs = [10.0, 8.0, 6.0, 4.0, 2.0]


#for i in [0,1,2,3,4,5]:
for i in [1,2,3,4,5]:
    xa, ya = pol2xy(rads[i], arrowang)
    x, y = pol2xy(rads[i]-1.0, textang)
    string = strings[i]
    #ax.annotate(string, xy=(xa, ya), xytext=(x+2.5, y),
    #            ha='center', va='center',
    #            arrowprops=dict(facecolor='black', shrink=0.05, width=0.5, headwidth=5.0, headlength=5),
    #            size=7.0,
    #            )
    ax.text(xa+2.8, ya, string, size=9.0, ha='center', va='center')


#x, y = pol2xy(10.5, textang)
#ax.text(x,y, r'$e^{-}$ $Z$', size = 10)
#x, y = pol2xy(7.0, textang)
#ax.text(x,y, r'$e^{-}$ $Z$ $n$', size = 10)


strings2=['$-$3 +3        ~1 cm',
          '6...8  ~10-100 m',
          '$\\approx$11.6      ~0.3 km',
          '$\\approx$14           ~1 km'
          ]

rads = [10.5, 9.0, 7.5, 6.0, 4.0]
ang = -pi/5
for i in range(4):
    
    x,y = pol2xy(rads[i], ang)
    ax.text(x, y, strings2[i], size=9, rotation=38, ha='right', va='top')
    

x, y = pol2xy(13.0, -pi/5)
ax.text(x,y, 'depth', rotation=-52, ha='right', va='top', size=9)

x, y = pol2xy(13.0, -pi/4.3)
ax.text(x,y, 'log$_{10}$ $\\rho$ \n [g cm$^{-3}$]', rotation=-52, ha='right', va='top', size=9)


plt.subplots_adjust(left=0.1, bottom=0.0, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('atmos_1.pdf', bbox_inches='tight')

