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
#from mpmath import lambertw
from scipy.special import lambertw
from scipy import optimize
from scipy import interpolate

import sympy
I0, Iph, Rs, Rsh, n, I, V, Vth, V_I0, I_I0, V_n, I_n = sympy.symbols('I0 Iph Rs Rsh n I V Vth V_I0 I_I0 V_n I_n')
#this stuff is from http://dx.doi.org/10.1016/j.solmat.2003.11.018
#symbolic representation for solar cell equation:
lhs = I
rhs = Iph-((V+I*Rs)/Rsh)-I0*(sympy.exp((V+I*Rs)/(n*Vth))-1)
charEqn = sympy.Eq(lhs,rhs)

#isolate current term in solar cell equation
#current = sympy.solve(charEqn,I)

#isolate voltage term in solar cell equation
#voltage = sympy.solve(charEqn,V)

#isolate I0  in solar cell equation
eyeNot = sympy.solve(charEqn,I0)

#import matplotlib.pyplot as plt
#plt.switch_backend("Qt5Agg")

# the data to be fit
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
#v = [1e3*x for x in v]
#i = [1e3*x for x in i]

#savemat("mat.mat",{'v':v,'i':i,'p0':guess})

def makeAGuess(VV,II):
    VV = np.array(VV)
    II = np.array(II)
    #data point selection:
    #lowest voltage (might be same as Isc)
    V_start_n = VV[0]
    I_start_n = II[0]
    
    #highest voltage
    V_end_n = VV[-1]
    I_end_n = II[-1]
    
    #Isc
    iFit = interpolate.interp1d(VV,II)
    V_sc_n = 0
    try:
        I_sc_n = float(iFit(V_sc_n))
    except:
        return([[nan,nan,nan,nan,nan], [nan,nan,nan,nan,nan], nan, "hard fail", 10])
    
    #mpp
    VVcalc = VV-VV[0]
    IIcalc = II-min(II)
    Pvirtual= np.array(VVcalc*IIcalc)
    vMaxIndex = Pvirtual.argmax()
    V_vmpp_n = VV[vMaxIndex]
    I_vmpp_n = II[vMaxIndex]
    
    #Vp: half way in voltage between vMpp and the start of the dataset:
    V_vp_n = (V_vmpp_n-V_start_n)/2 +V_start_n
    try:
        I_vp_n = float(iFit(V_vp_n))
    except:
        return([[nan,nan,nan,nan,nan], [nan,nan,nan,nan,nan], nan, "hard fail", 10])
    
    #Ip: half way in current between vMpp and the end of the dataset:
    I_ip_n = (I_vmpp_n-I_end_n)/2 + I_end_n
    iFit2 = interpolate.interp1d(VV,II-I_ip_n)
    try:
        V_ip_n = optimize.brentq(iFit2, VV[0], VV[-1])
    except:
        return([[nan,nan,nan,nan,nan], [nan,nan,nan,nan,nan], nan, "hard fail", 10])
    
    diaplayAllGuesses = False
    def evaluateGuessPlot(dataX, dataY, myguess):
        myguess = [float(x) for x in myguess]
        print("myguess:")
        print(myguess)
        vv=np.linspace(min(dataX),max(dataX),1000)
        ii=vectorizedCurrent(vv,myguess[0],myguess[1],myguess[2],myguess[3],myguess[4])
        plt.title('Guess and raw data')
        plt.plot(vv,ii)
        plt.scatter(dataX,dataY)
        plt.grid(b=True)
        plt.draw()
        plt.show()
    
    # phase 1 guesses:
    I_L_initial_guess = I_sc_n
    R_sh_initial_guess = 1e6
    
    # compute intellegent guesses for Iph, Rsh by forcing the curve through several data points and numerically solving the resulting system of eqns
    newRhs = rhs - I
    aLine = Rsh*V+Iph-I
    eqnSys1 = aLine.subs([(V,V_start_n),(I,I_start_n)])
    eqnSys2 = aLine.subs([(V,V_vp_n),(I,I_vp_n)])
    
    eqnSys = (eqnSys1,eqnSys2)
    
    try:
        nGuessSln = sympy.nsolve(eqnSys,(Iph,Rsh),(I_L_initial_guess,R_sh_initial_guess),maxsteps=10000)
    except:
        return([[nan,nan,nan,nan,nan], [nan,nan,nan,nan,nan], nan, "hard fail", 10])
    
    I_L_guess = nGuessSln[0]
    R_sh_guess = -1*1/nGuessSln[1]
    R_s_guess = -1*(V_end_n-V_ip_n)/(I_end_n-I_ip_n)
    n_initial_guess = 2 #TODO: maybe a more intelegant guess for n can be found using http://pvcdrom.pveducation.org/CHARACT/IDEALITY.HTM
    I0_initial_guess = eyeNot[0].evalf(subs={Vth:thermalVoltage,Rs:R_s_guess,Rsh:R_sh_guess,Iph:I_L_guess,n:n_initial_guess,I:I_ip_n,V:V_ip_n})                         
    
    initial_guess = [I0_initial_guess, I_L_guess, R_s_guess, R_sh_guess, n_initial_guess]
    if diaplayAllGuesses:
        evaluateGuessPlot(VV, II, initial_guess)
    
    # let's try the fit now, if it works great, we're done, otherwise we can continue
    #try:
        #guess = initial_guess
        #fitParams, fitCovariance, infodict, errmsg, ier = optimize.curve_fit(optimizeThis, VV, II,p0=guess,full_output = True,xtol=1e-13,ftol=1e-15)
        #return(fitParams, fitCovariance, infodict, errmsg, ier)
    #except:
        #pass        
    
    #refine guesses for I0 and Rs by forcing the curve through several data points and numerically solving the resulting system of eqns
    eqnSys1 = newRhs.subs([(Vth,thermalVoltage),(Iph,I_L_guess),(V,V_ip_n),(I,I_ip_n),(n,n_initial_guess),(Rsh,R_sh_guess)])
    eqnSys2 = newRhs.subs([(Vth,thermalVoltage),(Iph,I_L_guess),(V,V_end_n),(I,I_end_n),(n,n_initial_guess),(Rsh,R_sh_guess)])
    eqnSys = (eqnSys1,eqnSys2)
    
    try:
        nGuessSln = sympy.nsolve(eqnSys,(I0,Rs),(I0_initial_guess,R_s_guess),maxsteps=10000)
    except:
        return([[nan,nan,nan,nan,nan], [nan,nan,nan,nan,nan], nan, "hard fail", 10])
    
    I0_guess = nGuessSln[0]
    R_s_guess = nGuessSln[1]
    
    #Rs_initial_guess = RsEqn[0].evalf(subs={I0:I0_initial_guess,Vth:thermalVoltage,Rsh:R_sh_guess,Iph:I_L_guess,n:n_initial_guess,I:I_end_n,V:V_end_n})
    #I0_guess = I0_initial_guess
    #R_s_guess = Rs_initial_guess
    
    guess = [I0_guess, I_L_guess, R_s_guess, R_sh_guess, n_initial_guess]
    
    if diaplayAllGuesses:
        evaluateGuessPlot(VV, II, guess)
    
    #nidf
    
    #give 5x weight to data around mpp
    #nP = II*VV
    #maxIndex = np.argmax(nP)
    #weights = np.ones(len(II))
    #halfRange = (V_ip_n-VV[vMaxIndex])/2
    #upperTarget = VV[vMaxIndex] + halfRange
    #lowerTarget = VV[vMaxIndex] - halfRange
    #lowerTarget = 0
    #upperTarget = V_oc_n
    #lowerI = np.argmin(abs(VV-lowerTarget))
    #upperI = np.argmin(abs(VV-upperTarget))
    #weights[range(lowerI,upperI)] = 3
    #weights[maxnpi] = 10
    #todo: play with setting up "key points"
    
    guess = [float(x) for x in guess]
    #VV = [np.float(x) for x in VV]
    #II = [np.float(x) for x in II]
    
    #odrMod = odr.Model(odrThing)
    #myData = odr.Data(VV,II)
    #myodr = odr.ODR(myData, odrMod, beta0=guess,maxit=5000,sstol=1e-20,partol=1e-20)#
    #myoutput = myodr.run()
    #myoutput.pprint()
    #see http://docs.scipy.org/doc/external/odrpack_guide.pdf    
    return guess


# some constants
cellTemp = 29 #degC all analysis is done assuming the cell is at 29 degC
T = 273.15 + cellTemp #cell temp in K
K = 1.3806488e-23 #boltzman constant
q = 1.60217657e-19 #electron charge
thermalVoltage = K*T/q #thermal voltage ~26mv

# find the sum of the square of errors for a fit to some data
# given the fit function, the fit parameters and the x and y data
def sse(fun,params,x,y):
    return sum([(fun(X, *params)-Y)**2 for X,Y in zip(x,y)])
    
# the 0th branch of the lambertW function
def w(x):
    return np.real(lambertw(x, k=0, tol=np.finfo(float).eps))
    #return lambertw(x, k=0, tol=1e-15)

# here's the function we want to fit to
def optimizeThis (x, I0, Iph, Rs, Rsh, n):
    #return np.real((Rs*(I0*Rsh + Iph*Rsh - x) - thermalVoltage*n*(Rs + Rsh)*w(I0*Rs*Rsh*exp((Rs*(I0*Rsh + Iph*Rsh - x) + x*(Rs + Rsh))/(thermalVoltage*n*(Rs + Rsh)))/(thermalVoltage*n*(Rs + Rsh))))/(Rs*(Rs + Rsh)))
    return (Rs*(I0*Rsh + Iph*Rsh - x) - thermalVoltage*n*(Rs + Rsh)*w(I0*Rs*Rsh*exp((Rs*(I0*Rsh + Iph*Rsh - x) + x*(Rs + Rsh))/(thermalVoltage*n*(Rs + Rsh)))/(thermalVoltage*n*(Rs + Rsh))))/(Rs*(Rs + Rsh))

#def make_dictionary(max_length=10, **entries):
#    return dict([(key, entries[key]) for i, key in enumerate(entries.keys()) if i < max_length])

# my guess for the fit parameters
# [I0, I_L, R_s, R_sh, n]
guess = [7.974383037191594e-11, 0.00627619846736794, 12.743239329693433, 5694.842341863068, 2.0]
#guess = {'I0': 1e-15, 'Iph': 8.6889261535882785, 'Rs': 0.5155903798902628, 'Rsh': 788.22883822837741, 'n': 40.0}
#guess = [guess['I0'],guess['Iph'],guess['Rs'],guess['Rsh'],guess['n']]

#guess = [np.complex128(x) for x in guess]

# bounds on parameters
I0_bounds = [0, inf]
I_L_bounds = [0, inf]
R_s_bounds = [0, inf]
R_sh_bounds = [0, inf]
n_bounds = [0, 40]
myBounds=([I0_bounds[0], I_L_bounds[0], R_s_bounds[0], R_sh_bounds[0], n_bounds[0]], [I0_bounds[1], I_L_bounds[1], R_s_bounds[1], R_sh_bounds[1], n_bounds[1]])

# the best I can do with unconstrained, 'lm' fit method is SSE=7.00172395151e-08
fitReturn = optimize.curve_fit(optimizeThis, v, i, p0=guess, xtol=np.finfo(float).eps, method="lm", full_output=True, check_finite=False)

# the best I can do with constrained 'trf' method is SSE=7.001723218560001e-08 :-)
#fitReturn = optimize.curve_fit(optimizeThis, v, i, p0=guess, xtol=np.finfo(float).eps, bounds=myBounds, method="trf", x_scale="jac", jac ='cs', verbose=2, max_nfev=1200000)

fitParams = fitReturn[0]
SSE = sse(optimizeThis, fitParams, v, i)

finalOnes = guess = {'I0': fitParams[0], 'Iph': fitParams[1], 'Rs': fitParams[2], 'Rsh': fitParams[3], 'n': fitParams[4]}
print(finalOnes)
print("The sum of the square of the errors (SSE) is:")
print(SSE)
