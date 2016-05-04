#!/usr/bin/env python3

# written by Grey Christoforo <first name [at] last name [not] net>
# a place for playing around to try to get the best possible fits

import mpmath.libmp
assert mpmath.libmp.BACKEND == 'gmpy'
import numpy as np
from numpy import nan
from numpy import inf
from numpy import exp

from scipy.io import savemat

from scipy.special import lambertw

from scipy import optimize

#import matplotlib.pyplot as plt
#plt.switch_backend("Qt5Agg")

v = [-0.19994928, -0.14574878, -0.09155122, -0.03743893,  0.01663798,
       0.07081506,  0.12496336,  0.17912959,  0.23336065,  0.28753638,
       0.34162751,  0.39578977,  0.44999278,  0.50416487,  0.55837834,
       0.61247802,  0.6666159 ,  0.72084767,  0.77501357,  0.82919377,
       0.88335657,  0.93745446,  0.99165738,  1.04589999,  1.1000551 ]

i = [ 0.00631131,  0.00633053,  0.0062757 ,  0.00631208,  0.00625346,
        0.00628661,  0.00627005,  0.00628372,  0.00624534,  0.00622571,
        0.00622486,  0.00623647,  0.00620945,  0.00620305,  0.00621628,
        0.00608415,  0.00577381,  0.0052973 ,  0.00453748,  0.00347691,
        0.00208499,  0.00020583, -0.00227922, -0.00524322, -0.00856056]

#v = [np.complex128(x) for x in v]
#i = [np.complex128(x) for x in i]

cellTemp = 29 #degC all analysis is done assuming the cell is at 29 degC
T = 273.15 + cellTemp #cell temp in K
K = 1.3806488e-23 #boltzman constant
q = 1.60217657e-19 #electron charge
thermalVoltage = K*T/q #thermal voltage ~26mv

# find the sum of the square of errors for a fit to some data given the fit function, the fit parameters and the x and y data
def sse(fun,params,x,y):
    return sum([(fun(X, *params)-Y)**2 for X,Y in zip(x,y)])
    
# the lambertW function
def w(x):
    return np.real_if_close(lambertw(x, k=0, tol=1e-15))

# here's the function we want to fit to
def optimizeThis (x, I0, Iph, Rs, Rsh, n):
    #return (Rs*(I0*Rsh + Iph*Rsh - x) - thermalVoltage*n*(Rs + Rsh)*lambertw(I0*Rs*Rsh*np.exp((Rs*(I0*Rsh + Iph*Rsh - x) + x*(Rs + Rsh))/(thermalVoltage*n*(Rs + Rsh)))/(thermalVoltage*n*(Rs + Rsh)),tol=lambertWTol))/(Rs*(Rs + Rsh))
    return (Rs*(I0*Rsh + Iph*Rsh - x) - thermalVoltage*n*(Rs + Rsh)*w(I0*Rs*Rsh*exp((Rs*(I0*Rsh + Iph*Rsh - x) + x*(Rs + Rsh))/(thermalVoltage*n*(Rs + Rsh)))/(thermalVoltage*n*(Rs + Rsh))))/(Rs*(Rs + Rsh))

# [I0, I_L, R_s, R_sh, n]
guess = [7.974383037191594e-11, 0.00627619846736794, 12.743239329693433, 5694.842341863068, 2.0]
#guess = [np.complex128(x) for x in guess]

myXtol = np.float(1e-30)
myFtol = np.float(1e-30)
myGtol = np.float(1e-30)
myMax_nfev = 12000

I0_bounds = [0, inf]
I_L_bounds = [0, inf]
R_s_bounds = [0, inf]
R_sh_bounds = [0, inf]
n_bounds = [0, 4]

#I0_bounds = [-inf, inf]
#I_L_bounds = [-inf, inf]
#R_s_bounds = [-inf, inf]
#R_sh_bounds = [-inf, inf]
#n_bounds = [-inf, inf]

myBounds=([I0_bounds[0], I_L_bounds[0], R_s_bounds[0], R_sh_bounds[0], n_bounds[0]], [I0_bounds[1], I_L_bounds[1], R_s_bounds[1], R_sh_bounds[1], n_bounds[1]])
myMethod = 'trf'
#myMethod = 'dogbox'
#myMethod = 'lm'

#savemat("mat.mat",{'v':v,'i':i,'p0':guess})
#fitReturn = optimize.curve_fit(optimizeThis, v, i, p0=guess, bounds=myBounds, method=myMethod, jac ='cs', x_scale="jac", max_nfev=myMax_nfev, loss="soft_l1", tr_options={'tr_solver': "lsmr"}, verbose=2)

# the best i can do (constrained is SSE=6.70630326688e-07)
fitReturn = optimize.curve_fit(optimizeThis, v, i, p0=guess, bounds=myBounds, x_scale="jac", jac ='2-point', xtol=0, ftol=1e-17, gtol=1e-15, verbose=2,max_nfev=1200000)

# the best i can do (unconstrained is SSE=7.00172395151e-08)
#fitReturn = optimize.curve_fit(optimizeThis, v, i, p0=guess, xtol=1e-15, ftol=1e-15, gtol=1e-15,full_output=True)

fitParams = fitReturn[0]
SSE = sse(optimizeThis, fitParams, v, i)

print(fitReturn)
print("The sum of the square of the errors (SSE) is:")
print(SSE)
