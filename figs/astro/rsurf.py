# http://www.ster.kuleuven.be/~pieterd/python/html/plotting/mayavi_example.html

import numpy as np
from numpy import sin,cos,pi,sqrt # makes the code more readable
from scipy.optimize import newton

import sys

#mayavi 3d
import pylab as plt
from mayavi import mlab # or from enthought.mayavi import mlab



#matplotlib 3d
#import matplotlib.pyplot as plt
#from matplotlib import cm
#from mpl_toolkits.mplot3d import Axes3D

#fig = plt.figure(figsize=(8.0, 6.0)) #single column fig
#fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
#fig = plt.figure(figsize=(7.48, 2.5))  #two column figure
#ax = fig.add_subplot(111, projection='3d')

#ax.set_xlim(-1.5, 1.5)
#ax.set_ylim(-1.5, 1.5)
#ax.set_zlim(-1.5, 1.5)

#plt.rc('font', family='serif')
#plt.rc('xtick', labelsize=5)
#plt.rc('ytick', labelsize=5)
#plt.rc('axes', labelsize=5)


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
mlab.figure()
mlab.mesh(x1,y1,z1,scalars=r1)




#move secondary to proper place
rot_angle = pi
Rz = np.array([[cos(rot_angle),-sin(rot_angle),0],
               [sin(rot_angle), cos(rot_angle),0],
               [0,             0,              1]])
B = np.dot(Rz,np.array([x2,y2,z2]).reshape((3,-1))) # we need to have a 3x3 times 3xN array
x2,y2,z2 = B.reshape((3,x2.shape[0],x2.shape[1])) # but we want our original shape back
x2 += 1 # simple translation


#mlab.figure()
#mlab.mesh(x1,y1,z1,scalars=r1)
#mlab.mesh(x2,y2,z2,scalars=r2)


#ax.plot_surface(x1, y1, z2, r1)
#ax.plot_wireframe(x1, y1, z2)



#plt.show()
#plt.savefig("roche_surf.pdf")

#plt.imshow(r1)
#plt.colorbar()
#plt.figure()
#plt.show()
#plt.imshow(r2)
#plt.colorbar()


