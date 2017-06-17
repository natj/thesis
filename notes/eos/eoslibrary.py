import numpy as np
from units import mp
from units import c

eosLib = {
#    'PAL6'  :[ 34.380,  2.227,  2.189,  2.159, 'npem' ],
    'SLy'   :[ 34.384,  3.005,  2.988,  2.851, 'npem' ],
#    'APR1'  :[ 33.943,  2.442,  3.256,  2.908, 'npem' ],
#    'APR2'  :[ 34.126,  2.643,  3.014,  2.945, 'npem' ],
    'APR3'  :[ 34.392,  3.166,  3.573,  3.281, 'npem' ],
    'APR4'  :[ 34.269,  2.830,  3.445,  3.348, 'npem' ],
#    'FPS'   :[ 34.283,  2.985,  2.863,  2.600, 'npem' ],
    'WFF1'  :[ 34.031,  2.519,  3.791,  3.660, 'npem' ],
    'WFF2'  :[ 34.233,  2.888,  3.475,  3.517, 'npem' ],
#    'WFF3'  :[ 34.283,  3.329,  2.952,  2.589, 'npem' ],
#    'BBB2'  :[ 34.331,  3.418,  2.835,  2.832, 'npem' ],
#    'BPAL12':[ 34.358,  2.209,  2.201,  2.176, 'npem' ],
    'ENG'   :[ 34.437,  3.514,  3.130,  3.168, 'npem' ],
    'MPA1'  :[ 34.495,  3.446,  3.572,  2.887, 'npem' ],
    'MS1'   :[ 34.858,  3.224,  3.033,  1.325, 'npem' ],
#    'MS2'   :[ 34.605,  2.447,  2.184,  1.855, 'npem' ],
    'MS1b'  :[ 34.855,  3.456,  3.011,  1.425, 'npem' ],
#    'PS'    :[ 34.671,  2.216,  1.640,  2.365, 'meson' ],
#    'GS1a'  :[ 34.504,  2.350,  1.267,  2.421, 'meson' ],
#    'GS2a'  :[ 34.642,  2.519,  1.571,  2.314, 'meson' ],
#    'BGN1H1':[ 34.623,  3.258,  1.472,  2.464, 'hyperon' ],
#    'GNH3'  :[ 34.648,  2.664,  2.194,  2.304, 'hyperon' ],
#    'H1'    :[ 34.564,  2.595,  1.845,  1.897, 'hyperon' ],
#    'H2'    :[ 34.617,  2.775,  1.855,  1.858, 'hyperon' ],
#    'H3'    :[ 34.646,  2.787,  1.951,  1.901, 'hyperon' ],
    'H4'    :[ 34.669,  2.909,  2.246,  2.144, 'hyperon' ],
#    'H5'    :[ 34.609,  2.793,  1.974,  1.915, 'hyperon' ],
#    'H6a'   :[ 34.593,  2.637,  2.121,  2.064, 'hyperon' ],
#    'H7'    :[ 34.559,  2.621,  2.048,  2.006, 'hyperon' ],
#    'PCL2'  :[ 34.507,  2.554,  1.880,  1.977, 'hyperon' ],
#    'ALF1'  :[ 34.055,  2.013,  3.389,  2.033, 'quark' ],
    'ALF2'  :[ 34.616,  4.070,  2.411,  1.890, 'quark' ],
#    'ALF3'  :[ 34.283,  2.883,  2.653,  1.952, 'quark' ],
#    'ALF4'  :[ 34.314,  3.009,  3.438,  1.803, 'quark' ]
}


def EoS(name, rhos):

    rhos *= (1.66e-24/mp)

    ll = eosLib[ name ]

    p1 = 10**ll[0]
    g1 = ll[1]
    g2 = ll[2]
    g3 = ll[3]

    r1 = 2.8e14
    #r1 = 10**14.4
    r2 = 10**14.7
    r3 = 10**15.0

    #construct normalization constants starting from nuclear saturation depth
    K1 = p1/(r2**g1)
    K2 = (K1 * r2**g1)/r2**g2
    K3 = (K2 * r3**g2)/r3**g3


    Ps = np.zeros(len(rhos))
    for i, rho in enumerate(rhos):
        K = K3
        g = g3
        if rho < r3:
            K = K2
            g = g2
            if rho < r2:
                K = K1
                g = g1
        Ps[i] = K * rho**g

    #join smoothly to SLy crust
    for i, rho in enumerate(rhos):
        if rho < r1:
            Kc = 3.99874e-8
            Gc = 1.35692

            Pcrust = (Kc*rho**Gc)*c**2
            if Pcrust > Ps[i]:
                Ps[i] = Pcrust
            else:
                break


    return Ps 

#print EoS('PAL6', 1e14)
#print EoS('SLy', 1e14)
#print EoS('FPS', 1e14)




