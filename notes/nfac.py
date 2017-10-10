import numpy as np
import matplotlib.pyplot as plt


pi  = np.pi
c   = 2.99792458e10
mn  = 1.6749286e-24
mp  = 1.6726231e-24
dm  = mn - mp
rho = 1.0e14
h   = 6.6260755e-27
C   = ((3*h**3)/(8*pi))**(2/3)

print dm
print dm*c**2
print dm*c**2 * 6.242e11 / 1.0e6


#print "dm/m_n= {}".format(dm/mn)

#x = 0.0
#xn = x

#xn = (2./C)*( (dm * c**2 * mn**(2./3.))/(rho**(2./3.)) * (1. + x*(mp/mn)**(2./3.)))*( 1./mp + 2.*c )**(-1)
#print xn**(3./2.)
#--------------------------------------------------


def pf(n):
    return h*( (3.0*n)/(8*pi) )**(1./3.)

def kene(n, m):
    return (pf(n)**2)/(2.0*m)

def pene(m):
    return m*c**2

def uene(n):
    return pf(n)*c

def nnd(x, rho):
    return rho/ (mn*(1.0 + x*(mp/mn)) )

def relkene(n, m):
    return np.sqrt( uene(n)**2 + pene(m)**2 )


#non-rel massive particles
def ebalance0(x, rho):
    nn = nnd(x, rho)

    np = x*nn
    ne = np
    RH = pene(mn) + kene(nn, mn)
    LH = pene(mp) + kene(np, mp) + uene(ne)

    return RH-LH


#rel everything
def ebalance1(x, rho):
    nn = nnd(x, rho)

    np = x*nn
    ne = np

    RH = relkene(nn, mn)
    LH = relkene(np, mp) + uene(ne)

    return RH-LH


#approx x = 0 for rho
def ebalance2(x, rho):
    nn = nnd(0.0, rho)

    np = x*nn
    ne = np

    RH = relkene(nn, mn)
    LH = relkene(np, mp) + uene(ne)

    return RH-LH

def ebalanceAppr(rho):

    t1 = (2.0/C)
    t2 = (dm*c*c*mn**(2./3.))/(rho**(2./3.))
    t3 = 1.0/mn
    t4 = (1.0/mp + 2*c)
    return ( (t1*t2+t3)/t4 )**(3./2.)

#print ebalance(0.0)
#print ebalance(1.0)
#print ebalance(1./200.)

from scipy.optimize import brentq

N = 20
xs = np.zeros(N)
xs1 = np.zeros(N)
xs2 = np.zeros(N)
xs3 = np.zeros(N)

rhos = np.logspace(8, 15, N)
i = 0
for rho in rhos:
    xs[i]  = 1.0 / brentq(ebalance0, 0, 1.0, args=(rho))
    xs1[i] = 1.0 / brentq(ebalance1, 0, 1.0, args=(rho))
    xs2[i] = 1.0 / brentq(ebalance2, 0, 1.0, args=(rho))
    xs3[i] = 1.0 / ebalanceAppr(rho)

    print "rho:", np.log10(rho), xs[i], xs1[i], xs2[i], xs3[i]
    i += 1


plt.plot(np.log10(rhos), np.log10(xs), "k-")
plt.plot(np.log10(rhos), np.log10(xs1),"r-")
plt.plot(np.log10(rhos), np.log10(xs2),"b--")
plt.plot(np.log10(rhos), np.log10(xs3),"g--")
plt.show()


