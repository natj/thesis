# http://www.ster.kuleuven.be/~pieterd/python/html/plotting/mayavi_example.html

import numpy as np
from numpy import sin,cos,pi,sqrt # makes the code more readable
from scipy.optimize import newton

import sys

#mayavi 3d
#import pylab as plt
#from mayavi import mlab # or from enthought.mayavi import mlab

#matplotlib 3d
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(8.0, 6.0)) #single column fig
#fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
#fig = plt.figure(figsize=(7.48, 2.5))  #two column figure
ax = fig.add_subplot(111, projection='3d')


plt.rc('font', family='serif')
#plt.rc('xtick', labelsize=5)
#plt.rc('ytick', labelsize=5)
#plt.rc('axes', labelsize=5)


# Roche potential 
##################################################
def roche_pot_xy(x, y, q):
    r   = np.hypot(x,y)
    phi = np.arctan2(y, x)
    return roche_pot(r, phi, q)

def roche_pot(r,phi,q):
    theta = pi/2 #equatorial plane
    lamr,nu = r*cos(phi)*sin(theta),cos(theta)
    return -( 1./r  + q*( 1./sqrt(1. - 2*lamr + r**2) - lamr)  + 0.5*(q+1) * r**2 * (1-nu**2) )


#polar grid
#phi = np.linspace(0, 2*pi, 100)
#rad = np.linspace(0, 3.0,    100)
#x = np.outer( rad, np.sin(phi) )
#y = np.outer( rad, np.cos(phi) )

#Cartesian grid
x, y = np.mgrid[ -2.0:2.5:100j, -2.0:2.0:100j ]

q = 0.25
z = roche_pot_xy(x, y, q)

np.clip(z, -4.0, None, out=z)


ax.set_xlabel(r'$x$', size=18)
ax.set_ylabel(r'$y$', size=18)

ax.zaxis.set_rotate_label(False) 
ax.set_zlabel(r'$\Phi_{\mathrm{R}}$', size=18, rotation=0)
ax.minorticks_on()

ax.set_zlim(-4, -1.5)

#ax.plot_wireframe(x, y, z, )
ax.view_init(elev=60., azim=-120)


ax.plot_surface(x, y, z, 
        rstride=1,
        cstride=1,
        cmap = cm.plasma,
        linewidth=0.5,
        antialiased=False,
        )


#plt.show()
#plt.subplots_adjust(left=0.1, bottom=0.0, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('roche.pdf') #, bbox_inches='tight')

sys.exit()




def roche(r,theta,phi,pot,q):
    lamr,nu = r*cos(phi)*sin(theta),cos(theta)
    return (pot - (1./r  + q*( 1./sqrt(1. - 2*lamr + r**2) - lamr)  + 0.5*(q+1) * r**2 * (1-nu**2) ))



theta,phi = np.mgrid[0:np.pi:75j,-0.5*pi:1.5*np.pi:150j]


pot1,pot2 = 2.88,10.
q = 0.5

r_init = 1e-5

print "Newton iteration..."
r1 = [newton(roche,r_init,args=(th,ph,pot1,q)) for th,ph in zip(theta.ravel(),phi.ravel())]
r2 = [newton(roche,r_init,args=(th,ph,pot2,1./q)) for th,ph in zip(theta.ravel(),phi.ravel())]


r1 = np.array(r1).reshape(theta.shape)
r2 = np.array(r2).reshape(theta.shape)


#change to cartesian
x1 = r1*sin(theta)*cos(phi)
y1 = r1*sin(theta)*sin(phi)
z1 = r1*cos(theta)


x2 = r2*np.sin(theta)*np.cos(phi)
y2 = r2*np.sin(theta)*np.sin(phi)
z2 = r2*np.cos(theta)

#plot
print "plotting..."
#mlab.figure()
#mlab.mesh(x1,y1,z1,scalars=r1)


#move secondary to proper place
#rot_angle = pi
#Rz = np.array([[cos(rot_angle),-sin(rot_angle),0],
#               [sin(rot_angle), cos(rot_angle),0],
#               [0,             0,              1]])
#B = np.dot(Rz,np.array([x2,y2,z2]).reshape((3,-1))) # we need to have a 3x3 times 3xN array
#x2,y2,z2 = B.reshape((3,x2.shape[0],x2.shape[1])) # but we want our original shape back
#x2 += 1 # simple translation
#
#
#mlab.figure()
#mlab.mesh(x1,y1,z1,scalars=r1)
#mlab.mesh(x2,y2,z2,scalars=r2)


#ax.plot_surface(x1, y1, z2, r1)
#ax.plot_wireframe(x1, y1, z2)


ax.plot_wireframe(x1, y1, r1)

plt.show()
plt.savefig("roche.pdf")




#plt.imshow(r1)
#plt.colorbar()
#plt.figure()
#plt.show()
#plt.imshow(r2)
#plt.colorbar()


