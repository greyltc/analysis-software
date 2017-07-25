from batch_iv_analysis.batch_iv_analysis_UI import Ui_batch_iv_analysis

# for performance tuning
#import cProfile, pstats, io 
#pr = cProfile.Profile()

# needed for file watching
import time

#TODO: make area editable

# to speed this up
import concurrent.futures

from collections import OrderedDict
from io import StringIO
import os, sys, inspect, csv

from PyQt5.QtCore import QSettings, Qt, QSignalMapper, QFileSystemWatcher, QDir, QFileInfo, QObject, pyqtSignal, QRunnable
from PyQt5.QtWidgets import QApplication, QMainWindow, QDialog, QFileDialog, QTableWidgetItem, QCheckBox, QPushButton

import mpmath.libmp
assert mpmath.libmp.BACKEND == 'gmpy'
import numpy as np
import sympy
import math
from numpy import nan
from numpy import inf
from numpy import exp

import functools
import scipy

import scipy.io as sio

from scipy import odr
from scipy import interpolate
from scipy import optimize
from scipy import special
from scipy.stats.distributions import t #needed for confidence interval calculation #TODO: remove this in favor of uncertainties
import matplotlib.pyplot as plt
plt.switch_backend("Qt5Agg")
#from uncertainties import ufloat #TODO: switch to using this for the error calcs

# some global variables because I'm lazy
Voc = None
Isc = None
P_prime = None
I_eqn = None
slns = None
electricalModelVarsOnly = None
I0, Iph, Rs, Rsh, n, I, V, Vth = (None,)*8

sqcmpersqm = 10000 #cm^2 per m^2
stdIrridance = 1000 #[W/m^2] standard reporting irridance
mWperW = 1000 # mW per W

def main(args=None):
  """# a tool for analysing solar cell i-v curves
  # written by Grey Christoforo <first name [at] last name [not] net>
  # please cite our work if you can!
  # DOI: 10.3390/photonics2041101
  """

  if args is None:
    args = sys.argv[1:]

  # Do argument parsing here (eg. with argparse)

  app = QApplication(sys.argv)
  analysis = MainWindow()
  analysis.show()
  sys.exit(app.exec_())

def doSymbolicManipulations(fastAndSloppy=False):
  global Voc_eqn
  global Isc_eqn
  global P_prime
  global I_eqn
  global slns
  global electricalModelVarsOnly
  global I0, Iph, Rs, Rsh, n, I, V, Vth
  
  print("Hang tight, we're doing the one-time symbolic manipulations now...")

  # let's define some variables we'll use to do some symbolic equaiton manipulation
  modelSymbols = sympy.symbols('I0 Iph Rs Rsh n I V Vth', real=True, positive=True)
  I0, Iph, Rs, Rsh, n, I, V, Vth = modelSymbols
  modelConstants = (Vth,)
  modelVariables = tuple(set(modelSymbols)-set(modelConstants))

  # calculate values for our model's constants now
  cellTemp = 29 #degC all analysis is done assuming the cell is at 29 degC
  T = 273.15 + cellTemp #cell temp in K
  K = 1.3806488e-23 #boltzman constant
  q = 1.60217657e-19 #electron charge
  thermalVoltage = K*T/q #thermal voltage ~26mv
  valuesForConstants = (thermalVoltage,)

  # define cell circuit model here
  lhs = I
  rhs = Iph-((V+I*Rs)/Rsh)-I0*(sympy.exp((V+I*Rs)/(n*Vth))-1)
  electricalModel = sympy.Eq(lhs,rhs)
  electricalModelVarsOnly = electricalModel.subs(zip(modelConstants,valuesForConstants))

  # symbolically isolate each variable in our characteristic equation
  # then make substitutions for constants
  # then make the solutions ready for use numerically
  # NOTE: this is actually pretty computationally intense;
  # some solutions might contain the Lambert W "function"
  symSolutions = {} # with constants substituted in
  symSolutionsNoSubs = {} # all the symbols preserved
  slns = {} # solutions that are ready to use numerically

  # here we define any function substitutions we'll need for lambdification later
  if fastAndSloppy:
    # for fast and inaccurate math
    functionSubstitutions = {"LambertW" : scipy.special.lambertw, "exp" : np.exp}
  else:
    # this is a massive slowdown (forces a ton of operations into mpmath)
    # but gives _much_ better accuracy and aviods overflow warnings/errors...
    functionSubstitutions = {"LambertW" : mpmath.lambertw, "exp" : mpmath.exp}

  for symbol in modelSymbols:
    symSolutionsNoSubs[str(symbol)] = sympy.solve(electricalModel,symbol)[0]
    #symSolutionsNoSubs[str(symbol)] = sympy.solveset(electricalModel,symbol,domain=sympy.S.Reals).args[0] #solveset doesn't work here (yet)
    symSolutions[str(symbol)] = symSolutionsNoSubs[str(symbol)].subs(zip(modelConstants,valuesForConstants))
    remainingVariables = list(set(modelVariables)-set([symbol]))
    slns[str(symbol)] = sympy.lambdify(remainingVariables,symSolutions[str(symbol)],functionSubstitutions,dummify=False)

  # analytical solution for Voc:
  Voc_eqn = symSolutions['V'].subs(I,0)
  Voc_eqn = sympy.lambdify((I0,Rsh,Iph,n),Voc_eqn,functionSubstitutions,dummify=False)

  # analytical solution for Isc:
  Isc_eqn = symSolutions['I'].subs(V,0)
  Isc_eqn = sympy.lambdify((I0,Rsh,Rs,Iph,n),Isc_eqn,functionSubstitutions,dummify=False)

  # analytical solution for Pmax: 
  P = symSolutions['I']*V
  P_prime = sympy.diff(P,V)
  #V_max = sympy.solve(P_prime,V,check=False,implicit=True)
  #V_max = sympy.solveset(P_prime, V, domain=sympy.S.Reals) #TODO: this is not working, but it would be cool...
  #P_max = P.subs(V,V_max)
  #P_max = sympy.lambdify((I0,Rsh,Rs,Iph,n),P_max,functionSubstitutions,dummify=False)
  # since we can't do this analytically (yet) let's try numerically

  #sympy.pprint(V_max,use_unicode=True,wrap_line=False)
  #sys.exit(0)

  # this puts the symbolic solution for I from above into a format needed for curve_fit
  I_eqn = lambda x,a,b,c,d,e: np.array([slns['I'](I0=a, Iph=b, Rs=c, Rsh=d, n=e, V=v) for v in x]).astype(complex)
  
  return fastAndSloppy

# tests if string is a number
def isNumber(s):
  try:
    float(s)
  except ValueError:
    return False
  return True

# yanked from https://github.com/randlet/to-precision
def to_precision(x,p):
  """
  returns a string representation of x formatted with a precision of p

  Based on the webkit javascript implementation taken from here:
  https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
  """

  if x is None: # catch none
    return str(x)

  if not np.isfinite(x): # catch nan and inf
    return str(x)

  if x == 0.:
    return "0." + "0"*(p-1)

  out = []

  if x < 0:
    out.append("-")
    x = -x

  e = int(math.log10(x))
  tens = math.pow(10, e - p + 1)
  n = math.floor(x/tens)

  if n < math.pow(10, p - 1):
    e = e -1
    tens = math.pow(10, e - p+1)
    n = math.floor(x / tens)

  if abs((n + 1.) * tens - x) <= abs(n * tens -x):
    n = n + 1

  if n >= math.pow(10,p):
    n = n / 10.
    e = e + 1


  m = "%.*g" % (p, n)

  if e < -2 or e >= p:
    out.append(m[0])
    if p > 1:
      out.append(".")
      out.extend(m[1:p])
    out.append('e')
    if e > 0:
      out.append("+")
    out.append(str(e))
  elif e == (p -1):
    out.append(m)
  elif e >= 0:
    out.append(m[:e+1])
    if e+1 < len(m):
      out.append(".")
      out.extend(m[e+1:])
  else:
    out.append("0.")
    out.extend(["0"]*-(e+1))
    out.append(m)

  return "".join(out)


# needed for findKnotsAndCoefs belowfrom multiprocessing import Pool
def _compute_u(p, D, dydx, dx, dx1, n):
  if p is None or p != 0:
    data = [dx[1:n - 1], 2 * (dx[:n - 2] + dx[1:n - 1]), dx[:n - 2]]
    R = scipy.sparse.spdiags(data, [-1, 0, 1], n - 2, n - 2)

  if p is None or p < 1:
    Q = scipy.sparse.spdiags(
      [dx1[:n - 2], -(dx1[:n - 2] + dx1[1:n - 1]), dx1[1:n - 1]],
          [0, -1, -2], n, n - 2)
    QDQ = (Q.T * D * Q)
    if p is None or p < 0:
      # Estimate p
      p = 1. / \
        (1. + QDQ.diagonal().sum() /
               (100. * R.diagonal().sum() ** 2))

    if p == 0:
      QQ = 6 * QDQ
    else:
      QQ = (6 * (1 - p)) * (QDQ) + p * R
  else:
    QQ = R

  # Make sure it uses symmetric matrix solver
  ddydx = np.diff(dydx, axis=0)
  #sp.linalg.use_solver(useUmfpack=True)
  u = 2 * scipy.sparse.linalg.spsolve((QQ + QQ.T), ddydx)
  return u.reshape(n - 2, -1), p        

# calculates breakpoints and coefficents for a smoothed piecewise polynomial interpolant for 1d data
def findBreaksAndCoefs(xx, yy, p=None):
  var=1
  x, y = np.atleast_1d(xx, yy)
  x = x.ravel()
  dx = np.diff(x)
  must_sort = (dx < 0).any()
  if must_sort:
    ind = x.argsort()
    x = x[ind]
    y = y[..., ind]
    dx = np.diff(x)

  n = len(x)

  #ndy = y.ndim
  szy = y.shape

  nd =  np.prod(szy[:-1]).astype(np.int)
  ny = szy[-1]

  if n < 2:
    raise ValueError('There must be >=2 data points.')
  elif (dx <= 0).any():
    raise ValueError('Two consecutive values in x can not be equal.')
  elif n != ny:
    raise ValueError('x and y must have the same length.')

  dydx = np.diff(y) / dx

  if (n == 2):  # % straight line
    coefs = np.vstack([dydx.ravel(), y[0, :]])
  else:

    dx1 = 1. / dx
    D = scipy.sparse.spdiags(var * np.ones(n), 0, n, n)  # The varianceStringIO

    u, p = _compute_u(p, D, dydx, dx, dx1, n)
    dx1.shape = (n - 1, -1)
    dx.shape = (n - 1, -1)
    zrs = np.zeros(nd)
    if p < 1:
      # faster than yi-6*(1-p)*Q*u
      ai = (y - (6 * (1 - p) * D *
                 np.diff(np.vstack([zrs,
                                          np.diff(np.vstack([zrs, u, zrs]), axis=0) * dx1,
                                          zrs]), axis=0)).T).T
    else:
      ai = y.reshape(n, -1)

    ci = np.vstack([zrs, 3 * p * u])
    di = (np.diff(np.vstack([ci, zrs]), axis=0) * dx1 / 3)
    bi = (np.diff(ai, axis=0) * dx1 - (ci + di * dx) * dx)
    ai = ai[:n - 1, ...]
    if nd > 1:
      di = di.T
      ci = ci.T
      ai = ai.T
    if not any(di):
      if not any(ci):
        coefs = np.vstack([bi.ravel(), ai.ravel()])
      else:
        coefs = np.vstack([ci.ravel(), bi.ravel(), ai.ravel()])
    else:
      coefs = np.vstack(
        [di.ravel(), ci.ravel(), bi.ravel(), ai.ravel()])

  return coefs, x

# this function makes ultra-super-intelligent guesses for all the parameters
# of the equation that w're about to attempt to fit so that
# the (relatively dumb and fragile) final optimization/curve fitting routine
# has the best chance of giving good results
def makeAReallySmartGuess(VV,II,isDarkCurve):
  #data point selection:
  #lowest voltage (might be same as Isc)
  nPoints = len(VV)
  start_i = 0
  V_start_n = VV[start_i]
  I_start_n = II[start_i]

  #highest voltage
  end_i = nPoints-1
  V_end_n = VV[end_i]
  I_end_n = II[end_i]

  # let's start out with really dumb guesses for all our variables
  guess = {'I0':1e-9, 'Iph':I_start_n, 'Rs':5, 'Rsh':1e6, 'n':1.0}    

  # try to refine guess for Iph
  try:
    # interpolate to find short circuit current estimate
    iFit = interpolate.interp1d(VV,II)
    guess['Iph'] = iFit(0).item()
  except:
    print("Warning. You really should have some negative voltages in your data...")

  # find the curve "knee"
  if isDarkCurve:
    absDeltaCurrent = abs(np.ediff1d(II))
    nums, bins = np.histogram(absDeltaCurrent,bins=30)
    iBinMax = nums.argmax() + 2
    largestBin = bins[iBinMax]
    knee_i = np.argmax(absDeltaCurrent>largestBin)
  else:
    Pvirtual = VV*II
    knee_i = Pvirtual.argmax()
  V_vmpp_n = VV[knee_i]
  I_vmpp_n = II[knee_i]        


  #key point Vp: half way in voltage between vMPP and the start of the dataset:
  try:
    V_vp_n = (V_vmpp_n+V_start_n)*3/4
    vp_i = np.searchsorted(VV, V_vp_n)+1
    dummy = VV[vp_i]# this is to check if we have a valid index
  except:
    vp_i = round(knee_i*3/4)
    print("Warning. Major issue encountered while making guesses for fit parameters.")
  V_vp_n = VV[vp_i]
  I_vp_n = II[vp_i]


  #key point Ip: half way in current between vMPP and the end of the dataset:
  try:
    I_ip_n = (I_vmpp_n+I_end_n)/4
    ip_i = np.argmax(II<I_ip_n)
  except:
    ip_i = round((nPoints+knee_i)/4)
    print("Warning. Major issue encountered while making guesses for fit parameters.")
  V_ip_n = VV[ip_i]
  I_ip_n = II[ip_i]

  # make some guess for slopes (parasitic resistances)
  # using half-way points and end points:
  if (I_ip_n == I_end_n):
    guess['Rs'] = 10e9
  else:
    guess['Rs'] = -1*(V_end_n-V_ip_n)/(I_end_n-I_ip_n)
  if (I_start_n == I_vp_n):
    guess['Rsh'] = 10e9
  else:
    guess['Rsh'] = -1*(V_start_n-V_vp_n)/(I_start_n-I_vp_n)
  # using mpp and end points:
  #guess['Rs'] = -1*(V_end_n-V_vmpp_n)/(I_end_n-I_vmpp_n)
  #guess['Rsh'] = -1*(V_start_n-V_vmpp_n)/(I_start_n-I_vmpp_n)      

  # try to further refine guesses for Iph and Rsh
  aLine = lambda x,T,Y: [-y + -1/x[0]*t + x[1] for t,y in zip(T,Y)]
  x0 = np.array([guess['Rsh'],guess['Iph']]) # initial guess vector
  optimizeResult = scipy.optimize.least_squares(aLine, x0, jac='2-point', bounds=(u'-inf', u'inf'), 
                                                method='trf', 
                                                  ftol=1e-08, 
                            xtol=1e-08, 
                            gtol=1e-08, 
                            x_scale=1.0, 
                            loss='linear', 
                            f_scale=1.0, 
                            diff_step=None, 
                            tr_solver=None, 
                            tr_options={}, 
                            jac_sparsity=None, 
                            max_nfev=None, 
                            verbose=0, args=(VV[start_i:vp_i],II[start_i:vp_i]), 
                            kwargs={})
  if optimizeResult.success:
    guess['Rsh'] = optimizeResult.x[0]
    guess['Iph'] = optimizeResult.x[1]

  # try to further refine guesses for I0 and Rs
  aLine = lambda x,T,Y: [-y + -1/x[0]*t + x[1] for t,y in zip(T,Y)]
  x0 = np.array([guess['Rs'],V_end_n/guess['Rs']]) # initial guess vector
  optimizeResult = scipy.optimize.least_squares(aLine, x0,
                                                jac='2-point',
    bounds=(u'-inf', u'inf'), 
    method='trf', 
    ftol=1e-08, 
    xtol=1e-08, 
    gtol=1e-08, 
    x_scale=1.0, 
    loss='linear', 
    f_scale=1.0, 
    diff_step=None, 
    tr_solver=None, 
    tr_options={}, 
    jac_sparsity=None, 
    max_nfev=None, 
    verbose=0, args=(VV[ip_i:end_i],II[ip_i:end_i]), 
    kwargs={})
  if optimizeResult.success:
    guess['Rs'] = optimizeResult.x[0]
    slope = optimizeResult.x[1]

  guess['I0'] = float(slns['I0'](Iph=guess['Iph'],Rs=guess['Rs'],Rsh=guess['Rsh'],n=guess['n'],I=I_ip_n,V=V_ip_n))

  # if you'd like to see how good/bad our initial guesses are
  visualizeGuess = False
  if visualizeGuess:
    print("My guesses are",guess)
    vv=np.linspace(min(VV),max(VV),1000)
    ii=np.array([slns['I'](I0=guess['I0'], Iph=guess['Iph'], Rs=guess['Rs'], Rsh=guess['Rsh'], n=guess['n'], V=v).real for v in vv])
    ii2=np.array(aLine([guess['Rs'],slope],vv,np.zeros(len(vv)))) # Rs fit line
    ii3=np.array(aLine([guess['Rsh'],guess['Iph']],vv,np.zeros(len(vv)))) # Rsh fit line
    plt.title('Guess and raw data')
    plt.plot(vv,ii) # char eqn
    plt.plot(vv,ii2) # Rs fit line
    plt.plot(vv,ii3) # Rsh fit line
    plt.plot(V_ip_n,I_ip_n,'+r',markersize=10)
    plt.plot(V_vp_n,I_vp_n,'+r',markersize=10)
    plt.plot(V_vmpp_n,I_vmpp_n,'+r',markersize=10)
    plt.scatter(VV,II)
    plt.grid(b=True)
    yRange = max(II) - min(II)
    plt.ylim(min(II)-yRange*.1,max(II)+yRange*.1)
    plt.draw()
    plt.show()
    plt.pause(1)

  return guess

# here we attempt to fit the input data to the characteristic equation
def doTheFit(VV,II,guess,bounds,windowObj):
  x0 = [guess['I0'],guess['Iph'],guess['Rs'],guess['Rsh'],guess['n']]
  #x0 = [7.974383037191593e-06, 627.619846736794, 0.00012743239329693432, 0.056948423418631065, 2.0]
  #residuals = lambda x,T,Y: [-y + float(slns['I'](I0=x[0], Iph=x[1], Rs=x[2], Rsh=x[3], n=x[4], V=t).real) for t,y in zip(T,Y)]
  residuals = lambda x,T,Y: np.abs([np.real_if_close(slns['I'](I0=x[0], Iph=x[1], Rs=x[2], Rsh=x[3], n=x[4], V=t)) - y for t,y in zip(T,Y)])
  #residuals = lambda x,T,Y: np.array([-y + slns['I'](I0=x[0], Iph=x[1], Rs=x[2], Rsh=x[3], n=x[4], V=t) for t,y in zip(T,Y)]).astype('complex')
  fitArgs = (residuals,x0)
  fitKwargs = {}
  #fitKwargs['jac'] = '3-point'
  fitKwargs['jac'] = '2-point'
  #fitKwargs['jac'] = 'cs'
  fitKwargs['ftol'] = np.finfo(float).eps
  fitKwargs['xtol'] = np.finfo(float).eps
  fitKwargs['gtol'] = np.finfo(float).eps
  #fitKwargs['x_scale'] = list(map(lambda x: x/10, x0))
  #fitKwargs['x_scale'] = list(map(lambda x: x/100, x0))
  #fitKwargs['x_scale'] = list(map(lambda x: x*1000, x0))
  #fitKwargs['x_scale'] = x0
  fitKwargs['x_scale'] = 'jac'
  fitKwargs['loss'] = 'cauchy'
  #fitKwargs['loss'] = 'arctan'
  #fitKwargs['loss'] = 'linear'
  #fitKwargs['f_scale'] = 100000.0
  #fitKwargs['diff_step'] = list(map(lambda x: x/1000000, x0))
  #fitKwargs['diff_step'] = None
  fitKwargs['tr_solver'] = 'lsmr'
  fitKwargs['tr_solver'] = None
  fitKwargs['tr_options'] = {'regularize':True}
  #fitKwargs['jac_sparsity'] = None
  fitKwargs['max_nfev'] = 100000
  fitKwargs['verbose'] = windowObj.ui.verbositySpinBox.value()
  if windowObj.ui.fitMethodComboBox.currentIndex() == 0:
    fitKwargs['method'] = 'trf'
  elif windowObj.ui.fitMethodComboBox.currentIndex() == 1:
    fitKwargs['method'] = 'dogbox'
  elif windowObj.ui.fitMethodComboBox.currentIndex() == 2:
    fitKwargs['method'] = 'lm'
    fitKwargs['loss'] = 'linear' # loss must be linear for lm method
  fitKwargs['args'] = (VV,II)
  fitKwargs['kwargs'] = {}

  # do a constrained fit when one of the bounds is not inf or -inf
  # TODO: fix me! this is not working right now!
  #if sum(sum([np.isinf(value) for key,value in bounds.items()])) != 10:
  #    residuals = lambda x,T,Y: np.array([-y + slns['I'](I0=x[0], Iph=x[1], Rs=x[2], Rsh=x[3], n=x[4], V=t) for t,y in zip(T,Y)]).astype('complex')
  #    fitKwargs['jac'] = 'cs'
  #    fitKwargs['method'] = 'trf'
  #    fitKwargs['bounds'] = [(u'-inf',u'-inf',u'-inf',u'-inf',0),(u'inf',u'inf',u'inf',u'inf',u'inf')]
  #    sqrtEPS = np.finfo(float).eps**(1/2)
  #    fitKwargs['diff_step'] = [x0[0]/10, sqrtEPS, sqrtEPS, sqrtEPS, sqrtEPS]
  #    fitKwargs['max_nfev'] = 1200

  # do the fit
  optimizeResult = scipy.optimize.least_squares(*fitArgs,**fitKwargs)

  # do the fit with curve_fit
  #tehf = lambda XX,m_I0,m_Iph,m_Rs,m_Rsh,m_n: np.real_if_close(slns['I'](I0=m_I0, Iph=m_Iph, Rs=m_Rs, Rsh=m_Rsh, n=m_n, V=XX))
  #fit_result = scipy.optimize.curve_fit(tehf,VV,II,p0=x0,method='trf',verbose=2,x_scale= list(map(lambda x: x/1000, x0)))



  #optimizeResult.success = False
  #optimize.curve_fit(I_eqn, VV, II, p0=x0, bounds=fitKwargs['bounds'], diff_step=fitKwargs['diff_step'], method="trf", x_scale="jac", jac ='cs', verbose=1, max_nfev=1200000)
  #scipy.optimize.least_squares(residuals, np.array([  1.20347834e-13,   6.28639109e+02,   1.83005279e-04, 6.49757268e-02,   1.00000000e+00]), jac='cs', bounds=[('-inf', '-inf', '-inf', '-inf', 0), ('inf', 'inf', 'inf', 'inf', 'inf')], method='trf', max_nfev=12000, x_scale='jac', verbose=1,diff_step=[1.203478342631369e-14, 1.4901161193847656e-08, 1.4901161193847656e-08, 1.4901161193847656e-08, 1.4901161193847656e-08])
  if optimizeResult.success:
    #print(optimizeResult)
    if np.any(np.isnan(optimizeResult.jac)):
      sigmas = np.empty(len(optimizeResult.x))
      sigmas[:] = nan
    else:
      # calculate covariance matrix
      # Do Moore-Penrose inverse discarding zero singular values.
      _, s, VT = scipy.linalg.svd(optimizeResult.jac, full_matrices=False)
      threshold = np.finfo(float).eps * max(optimizeResult.jac.shape) * s[0]
      s = s[s > threshold]
      VT = VT[:s.size]
      pcov = np.dot(VT.T / s**2, VT)
      # now find sigmas from covariance matrix
      error = [] 
      for i in range(len(optimizeResult.x)):
        try:
          error.append(np.absolute(pcov[i][i])**0.5)
        except:
          error.append( 0.00 )
      sigmas = np.array(error)

    return {'success':True,'optParams':optimizeResult.x,'sigmas':sigmas,'message':optimizeResult.message,'SSE':optimizeResult.cost*2}
  else:
    return {'success':False,'message':optimizeResult.message}


  #if constrainedFit:
    ## handle the case when the user sets lower bound=upper bound
    ## (take that variable out of the optimization)
    #myKwargs = {}
    #finalFitValues = {}
    #finalSigmaValues = {}
    #curve_fit_guess = np.array([])
    #curve_fit_bounds=([],[])
    #paramNames = [] # need this to keep track of where each parameter is
    #if bounds['I0'][0] == bounds['I0'][1]:
      #myKwargs['I0'] = bounds['I0'][0]
      #finalFitValues['I0'] = myKwargs['I0']
      #finalSigmaValues['I0'] = 0
    #else:
      #curve_fit_guess = np.append(curve_fit_guess,guess['I0'])
      #curve_fit_bounds[0].append(bounds['I0'][0])
      #curve_fit_bounds[1].append(bounds['I0'][1])
      #paramNames.append("I0")
    #if bounds['Iph'][0] == bounds['Iph'][1]:
      #myKwargs['Iph'] = bounds['Iph'][0]
      #finalFitValues['Iph'] = myKwargs['Iph']
      #finalSigmaValues['Iph'] = 0
    #else:
      #curve_fit_guess = np.append(curve_fit_guess,guess['Iph'])
      #curve_fit_bounds[0].append(bounds['Iph'][0])
      #curve_fit_bounds[1].append(bounds['Iph'][1])
      #paramNames.append("Iph")
    #if bounds['Rs'][0] == bounds['Rs'][1]:
      #myKwargs['Rs'] = bounds['Rs'][0]
      #finalFitValues['Rs'] = myKwargs['Rs']
      #finalSigmaValues['Rs'] = 0
    #else:
      #curve_fit_guess = np.append(curve_fit_guess,guess['Rs'])
      #curve_fit_bounds[0].append(bounds['Rs'][0])
      #curve_fit_bounds[1].append(bounds['Rs'][1])
      #paramNames.append("Rs")
    #if bounds['Rsh'][0] == bounds['Rsh'][1]:
      #myKwargs['Rsh'] = bounds['Rsh'][0]
      #finalFitValues['Rsh'] = myKwargs['Rsh']
      #finalSigmaValues['Rsh'] = 0
    #else:
      #curve_fit_guess = np.append(curve_fit_guess,guess['Rsh'])
      #curve_fit_bounds[0].append(bounds['Rsh'][0])
      #curve_fit_bounds[1].append(bounds['Rsh'][1])
      #paramNames.append("Rsh")
    #if bounds['n'][0] == bounds['n'][1]:
      #myKwargs['n'] = bounds['n'][0]
      #finalFitValues['n'] = myKwargs['n']
      #finalSigmaValues['n'] = 0
    #else:
      #curve_fit_guess = np.append(curve_fit_guess,guess['n'])
      #curve_fit_bounds[0].append(bounds['n'][0])
      #curve_fit_bounds[1].append(bounds['n'][1])
      #paramNames.append("n")

    #redirected_output = sys.stdout = StringIO()
    #redirected_error = sys.stderr = StringIO()

    #try:
      #fitParams, fitCovariance = optimize.curve_fit(I_eqn, VV, II, p0=curve_fit_guess, bounds=curve_fit_bounds, method="trf", x_scale="jac", verbose=1, max_nfev=1200000)
    #except:
      #sys.stdout = sys.__stdout__
      #sys.stderr = sys.__stderr__                    
      #return([[nan,nan,nan,nan,nan], [nan,nan,nan,nan,nan], nan, "Unexpected Error: " + str(sys.exc_info()[1]) , 10])
    #out = redirected_output.getvalue()
    #err = redirected_error.getvalue()                
    #sys.stdout = sys.__stdout__
    #sys.stderr = sys.__stderr__
    #infodict = out
    #errmsg = out.splitlines()[0]
    #ier = 0

    #sigmas = np.sqrt(np.diag(fitCovariance))

    ## need this in case we fit less than 5 parameters
    #for sigma,paramValue,paramName in zip(sigmas,fitParams,paramNames):
      #finalFitValues[paramName] = paramValue
      #finalSigmaValues[paramName] = sigma
    #fitParams = [finalFitValues['I0'],finalFitValues['Iph'],finalFitValues['Rs'],finalFitValues['Rsh'],finalFitValues['n']]
    #sigmas = [finalSigmaValues['I0'],finalSigmaValues['Iph'],finalSigmaValues['Rs'],finalSigmaValues['Rsh'],finalSigmaValues['n']]
  #else: # unconstrained "l-m" fit
    #curve_fit_guess = [guess['I0'],guess['Iph'],guess['Rs'],guess['Rsh'],guess['n']]    
    #try:
      #fitParams, fitCovariance, infodict, errmsg, ier = optimize.curve_fit(lambda *args: I_eqn(*args).astype(float), VV, II, p0=curve_fit_guess, method="lm", full_output=True)
    #except:
      #return([[nan,nan,nan,nan,nan], [nan,nan,nan,nan,nan], nan, "Unexpected Error: " + str(sys.exc_info()[1]) , 10])
    #sigmas = np.sqrt(np.diag(fitCovariance))
  #return(fitParams, sigmas, infodict, errmsg, ier)

# this routine analyzes/quantifies the goodness of our fit
def analyzeGoodness(VV,II,params,guess,msg):
  # sum of square of differences between data and fit [A^2]
  print("fit:")
  print(params)                
  print("guess:")
  print(guess)
  print("message:")
  print(msg)

  vv=np.linspace(VV[0],VV[-1],1000)
  ii=np.array([slns['I'](I0=guess['I0'], Iph=guess['Iph'], Rs=guess['Rs'], Rsh=guess['Rsh'], n=guess['n'], V=v).real for v in vv])
  ii2=np.array([slns['I'](I0=params['I0'], Iph=params['Iph'], Rs=params['Rs'], Rsh=params['Rsh'], n=params['n'], V=v).real for v in vv])
  plt.title('Fit analysis')
  p1, = plt.plot(vv,ii, label='Guess',ls='--')
  p2, = plt.plot(vv,ii2, label='Fit')
  p3, = plt.plot(VV,II,ls='None',marker='o', label='Data')
  #p4, = plt.plot(VV[range(lowerI,upperI)],II[range(lowerI,upperI)],ls="None",marker='o', label='5x Weight Data')
  ax = plt.gca()
  handles, labels = ax.get_legend_handles_labels()
  ax.legend(handles, labels, loc=3)
  plt.grid(b=True)
  plt.draw()
  plt.show()

class col:
  header = ''
  position = 0
  tooltip = ''

class WorkerSignals(QObject):
  result = pyqtSignal(dict)

# this will process a file
class Worker(QRunnable):

  def __init__(self,mainWindow,fileName):
    super(Worker, self).__init__()
    self.mainWindow = mainWindow
    self.fileName = fileName
    self.signals = WorkerSignals()

  def run(self):
    result = self.mainWindow.processFile(self.fileName)
    self.signals.result.emit(result)

##class doSymbolicMathSignals(QObject):
  ##result = pyqtSignal(bool)

### this will process a file
##class doSymbolicMath(QRunnable):

  ##def __init__(self,mainWindow):
    ##super(doSymbolicMath, self).__init__()
    ##self.mainWindow = mainWindow
    ##self.signals = doSymbolicMathSignals()

  ##def run(self):
    ##sloppy = doSymbolicManipulations(fastAndSloppy=self.mainWindow.ui.doFastAndSloppyMathCheckBox.isChecked())
    ##self.signals.result.emit(sloppy)

class MainWindow(QMainWindow):
  workingDirectory = ''
  fileNames = []
  supportedExtensions = ['*.csv','*.tsv','*.txt','*.liv1','*.liv2','*.div1','*.div2']
  bounds = {}
  bounds['I0'] = [0, inf] 
  bounds['Iph'] = [0, inf]
  bounds['Rs'] = [0, inf]
  bounds['Rsh'] = [0, inf]
  bounds['n'] = [0, inf]
  symbolCalcsNotDone = True
  upperVLim = float('inf')
  lowerVLim = float('-inf')

  # for table
  rows = 0 #this variable keepss track of how many rows there are in the results table
  cols = OrderedDict()    

  def __init__(self):
    QMainWindow.__init__(self)

    self.settings = QSettings("greyltc", "batch-iv-analysis")

    # populate column headers
    thisKey = 'plotBtn'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'Draw Plot'
    self.cols[thisKey].tooltip = 'Click this button to draw a plot for that row'        

    thisKey = 'exportBtn'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'Export'
    self.cols[thisKey].tooltip = 'Click this button to export\ninterpolated data points from fits'        

    thisKey = 'file'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'File'
    self.cols[thisKey].tooltip = 'File name\nHover to see header from data file'

    thisKey = 'SSE'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'SSE\n[mA^2]'
    self.cols[thisKey].tooltip = 'Sum of the square of the errors between the data points and the fit to the char. eqn. (a measure of fit goodness)'

    thisKey = 'pce_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'PCE\n[%]'
    self.cols[thisKey].tooltip = 'Power conversion efficiency as found from spline fit'

    thisKey = 'pmax_a_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'P_max\n[mW/cm^2]'
    self.cols[thisKey].tooltip = 'Maximum power density as found from spline fit'

    thisKey = 'jsc_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'J_sc\n[mA/cm^2]'
    self.cols[thisKey].tooltip = 'Short-circuit current density as found from spline spline fit V=0 crossing'

    thisKey = 'voc_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'V_oc\n[mV]'
    self.cols[thisKey].tooltip = 'Open-circuit voltage as found from spline fit I=0 crossing'

    thisKey = 'ff_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'FF\n[%]'
    self.cols[thisKey].tooltip = 'Fill factor as found from spline fit'

    thisKey = 'rs_a'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'R_s\n[ohm*cm^2]'
    self.cols[thisKey].tooltip = 'Specific series resistance as found from characteristic equation fit'

    thisKey = 'rsh_a'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'R_sh\n[ohm*cm^2]'
    self.cols[thisKey].tooltip = 'Specific shunt resistance as found from characteristic equation fit'

    thisKey = 'jph'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'J_ph\n[mA/cm^2]'
    self.cols[thisKey].tooltip = 'Photogenerated current density as found from characteristic equation fit'

    thisKey = 'j0'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'J_0\n[nA/cm^2]'
    self.cols[thisKey].tooltip = 'Reverse saturation current density as found from characteristic equation fit'

    thisKey = 'n'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'n'
    self.cols[thisKey].tooltip = 'Diode ideality factor as found from characteristic equation fit'

    thisKey = 'vmax_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'V_max\n[mV]'
    self.cols[thisKey].tooltip = 'Voltage at maximum power point as found from spline fit'

    thisKey = 'area'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'Area\n[cm^2]'
    self.cols[thisKey].tooltip = 'Device area'

    thisKey = 'suns'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'Suns\n'
    self.cols[thisKey].tooltip = 'Illumination intensity'        

    thisKey = 'pmax_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'P_max\n[mW]'
    self.cols[thisKey].tooltip = 'Maximum power as found from spline fit'

    thisKey = 'pce_fit'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'PCE_fit\n[%]'
    self.cols[thisKey].tooltip = 'Power conversion efficiency as found from characteristic equation fit'

    thisKey = 'pmax_fit'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'P_max_fit\n[mW]'
    self.cols[thisKey].tooltip = 'Maximum power as found from characteristic equation fit'

    thisKey = 'pmax_a_fit'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'P_max_fit\n[mW/cm^2]'
    self.cols[thisKey].tooltip = 'Maximum power density as found from characteristic equation fit'

    thisKey = 'vmax_fit'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'V_max_fit\n[mV]'
    self.cols[thisKey].tooltip = 'Voltage at maximum power point as found from characteristic equation fit'

    thisKey = 'voc_fit'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'V_oc_fit\n[mV]'
    self.cols[thisKey].tooltip = 'Open-circuit voltage as found from characteristic equation fit I=0 crossing'

    thisKey = 'ff_fit'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'FF_fit\n[%]'
    self.cols[thisKey].tooltip = 'Fill factor as found from characteristic equation fit'        

    thisKey = 'isc_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'I_sc\n[mA]'
    self.cols[thisKey].tooltip = 'Short-circuit current as found from spline V=0 crossing'

    thisKey = 'isc_fit'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'I_sc_fit\n[mA]'
    self.cols[thisKey].tooltip = 'Short-circuit current as found from characteristic equation fit V=0 crossing'

    thisKey = 'jsc_fit'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'J_sc_fit\n[mA/cm^2]'
    self.cols[thisKey].tooltip = 'Short-circuit current density as found from characteristic equation fit V=0 crossing'        

    thisKey = 'iph'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'I_ph\n[mA]'
    self.cols[thisKey].tooltip = 'Photogenerated current as found from characteristic equation fit'

    thisKey = 'jph'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'J_ph\n[mA/mc^2]'
    self.cols[thisKey].tooltip = 'Photogenerated current density as found from characteristic equation fit'

    thisKey = 'i0'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'I_0\n[nA]'
    self.cols[thisKey].tooltip = 'Reverse saturation current as found from characteristic equation fit'

    thisKey = 'j0'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'J_0\n[nA/cm^2]'
    self.cols[thisKey].tooltip = 'Reverse saturation current density as found from characteristic equation fit'

    thisKey = 'rs'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'R_s\n[ohm]'
    self.cols[thisKey].tooltip = 'Series resistance as found from characteristic equation fit'

    thisKey = 'rsh'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'R_sh\n[ohm]'
    self.cols[thisKey].tooltip = 'Shunt resistance as found from characteristic equation fit'

    #how long status messages show for
    self.messageDuration = 2500#ms

    # Set up the user interface from Designer.
    self.ui = Ui_batch_iv_analysis()
    self.ui.setupUi(self)

    # load setting for computation threads
    if not self.settings.contains('threads'):
      self.threads = 4
      self.settings.setValue('threads',self.threads)
    else:
      self.threads = int(self.settings.value('threads'))
    self.settings.setValue('threads',8)

    # load setting for lower voltage cuttoff
    if not self.settings.contains('lowerVoltageCutoff'):
      self.ui.lowerVoltageCutoffLineEdit.setText('-inf')
      self.settings.setValue('lowerVoltageCutoff','-inf')
    else:
      self.ui.lowerVoltageCutoffLineEdit.setText(self.settings.value('lowerVoltageCutoff'))
      self.lowerVLim=float(self.settings.value('lowerVoltageCutoff'))
    self.ui.lowerVoltageCutoffLineEdit.editingFinished.connect(self.handleLowerLimChange)

    # load setting for upper voltage cuttoff
    if not self.settings.contains('upperVoltageCutoff'):
      self.ui.upperVoltageCutoffLineEdit.setText('inf')
      self.settings.setValue('upperVoltageCutoff','inf')
    else:
      self.ui.upperVoltageCutoffLineEdit.setText(self.settings.value('upperVoltageCutoff'))
      self.upperVLim=float(self.settings.value('upperVoltageCutoff'))
    self.ui.upperVoltageCutoffLineEdit.editingFinished.connect(self.handleUpperLimChange)

    # load setting for fast vs accurate calculations
    if not self.settings.contains('fastAndSloppy'):
      self.ui.doFastAndSloppyMathCheckBox.setChecked(False)
      self.settings.setValue('fastAndSloppy',False)
    else:
      self.ui.doFastAndSloppyMathCheckBox.setChecked(self.settings.value('fastAndSloppy') == 'true')
    self.ui.doFastAndSloppyMathCheckBox.stateChanged.connect(self.handleMathChange)

    # load setting for fitting eqn or not
    if not self.settings.contains('fitToEqn'):
      self.ui.attemptCharEqnFitCheckBox.setChecked(False)
      self.settings.setValue('fitToEqn',False)
    else:
      self.ui.attemptCharEqnFitCheckBox.setChecked(self.settings.value('fitToEqn') == 'true')
    self.ui.attemptCharEqnFitCheckBox.stateChanged.connect(self.handleEqnFitChange)

    # set defaults
    I0_lb_string = "0" if not self.settings.contains('I0_lb') else self.settings.value('I0_lb')
    Iph_lb_string = "0" if not self.settings.contains('Iph_lb') else self.settings.value('Iph_lb')
    Rs_lb_string = "0" if not self.settings.contains('Rs_lb') else self.settings.value('Rs_lb')
    Rsh_lb_string = "0" if not self.settings.contains('Rsh_lb') else self.settings.value('Rsh_lb')
    n_lb_string = "0" if not self.settings.contains('n_lb') else self.settings.value('n_lb')

    I0_ub_string = "inf" if not self.settings.contains('I0_ub') else self.settings.value('I0_ub')
    Iph_ub_string = "inf" if not self.settings.contains('Iph_ub') else self.settings.value('Iph_ub')
    Rs_ub_string = "inf" if not self.settings.contains('Rs_ub') else self.settings.value('Rs_ub')
    Rsh_ub_string = "inf" if not self.settings.contains('Rsh_ub') else self.settings.value('Rsh_ub')
    n_ub_string = "inf" if not self.settings.contains('n_ub') else self.settings.value('n_ub')

    if self.settings.contains('fitMethod'):
      self.ui.fitMethodComboBox.setCurrentIndex(int(self.settings.value('fitMethod')))
    else:
      self.settings.setValue('fitMethod',self.ui.fitMethodComboBox.currentIndex())

    if self.settings.contains('verbosity'):
      self.ui.verbositySpinBox.setValue(int(self.settings.value('verbosity')))
    else:
      self.settings.setValue('verbosity',self.ui.verbositySpinBox.value())

    self.bounds['I0'][0] = np.float(I0_lb_string)
    self.bounds['Iph'][0] = np.float(Iph_lb_string)
    self.bounds['Rs'][0] = np.float(Rs_lb_string)
    self.bounds['Rsh'][0] = np.float(Rsh_lb_string)
    self.bounds['n'][0] = np.float(n_lb_string)

    self.bounds['I0'][1] = np.float(I0_ub_string)
    self.bounds['Iph'][1] = np.float(Iph_ub_string)
    self.bounds['Rs'][1] = np.float(Rs_ub_string)
    self.bounds['Rsh'][1] = np.float(Rsh_ub_string)
    self.bounds['n'][1] = np.float(n_ub_string)

    self.ui.I0_lb.setText(I0_lb_string)
    self.ui.Iph_lb.setText(Iph_lb_string)
    self.ui.Rs_lb.setText(Rs_lb_string)
    self.ui.Rsh_lb.setText(Rsh_lb_string)
    self.ui.n_lb.setText(n_lb_string)

    self.ui.I0_ub.setText(I0_ub_string)
    self.ui.Iph_ub.setText(Iph_ub_string)
    self.ui.Rs_ub.setText(Rs_ub_string)
    self.ui.Rsh_ub.setText(Rsh_ub_string)
    self.ui.n_ub.setText(n_ub_string)

    # connect the bounds change handler
    self.ui.I0_lb.editingFinished.connect(self.handleConstraintsChange)
    self.ui.Iph_lb.editingFinished.connect(self.handleConstraintsChange)
    self.ui.Rs_lb.editingFinished.connect(self.handleConstraintsChange)
    self.ui.Rsh_lb.editingFinished.connect(self.handleConstraintsChange)
    self.ui.n_lb.editingFinished.connect(self.handleConstraintsChange)

    self.ui.I0_ub.editingFinished.connect(self.handleConstraintsChange)
    self.ui.Iph_ub.editingFinished.connect(self.handleConstraintsChange)
    self.ui.Rs_ub.editingFinished.connect(self.handleConstraintsChange)
    self.ui.Rsh_ub.editingFinished.connect(self.handleConstraintsChange)
    self.ui.n_ub.editingFinished.connect(self.handleConstraintsChange)

    self.ui.fitMethodComboBox.currentIndexChanged.connect(self.handleFitMethodChange)

    self.ui.resetSettingsButton.clicked.connect(self.resetDefaults)

    self.ui.verbositySpinBox.valueChanged.connect(self.handleVerbosityChange)

    #insert cols
    for item in self.cols:
      blankItem = QTableWidgetItem()
      thisCol = list(self.cols.keys()).index(item)
      self.ui.tableWidget.insertColumn(thisCol)
      blankItem.setToolTip(self.cols[item].tooltip)
      blankItem.setText(self.cols[item].header)
      self.ui.tableWidget.setHorizontalHeaderItem(thisCol,blankItem)

    #file system watcher
    self.watcher = QFileSystemWatcher(self)
    self.watcher.directoryChanged.connect(self.handleWatchUpdate)

    #connect signals generated by gui elements to proper functions 
    self.ui.actionOpen.triggered.connect(self.openCall)
    self.ui.actionEnable_Watching.triggered.connect(self.watchCall)
    self.ui.actionSave.triggered.connect(self.handleSave)
    self.ui.actionWatch_2.triggered.connect(self.handleWatchAction)
    self.ui.statusbar.messageChanged.connect(self.statusChanged)

    self.ui.actionClear_Table.triggered.connect(self.clearTableCall)

    #override showMessage for the statusbar
    self.oldShowMessage = self.ui.statusbar.showMessage
    self.ui.statusbar.showMessage = self.myShowMessage

    # this pool holds the workers
    self.pool = concurrent.futures.ProcessPoolExecutor(max_workers=self.threads)
    
    # do symbolic calcs now if needed
    if self.ui.attemptCharEqnFitCheckBox.isChecked():
      doSymbolicManipulations(fastAndSloppy=self.ui.doFastAndSloppyMathCheckBox.isChecked())
       
  def handleMathFinished(self,sloppy):
    tp = QThreadPool.globalInstance()
    tp.setMaxThreadCount(self.threads)
    self.symbolCalcsNotDone = False
    print("One-time symbolic manipulations done! Fast and sloppy mode =", sloppy)

  def processFitResult(self,result):
    print('Got new fit result...')
    print(result)
    tp = QThreadPool.globalInstance()
    activeThreads = tp.activeThreadCount()
    #print('Active threads: '+str(activeThreads))
    if activeThreads == 0:
      tp.waitForDone()
    #self.ui.statusbar.showMessage('Active threads: '+str(tp.activeThreadCount()),500)

  def resetDefaults(self):
    self.ui.attemptCharEqnFitCheckBox.setChecked(True)
    self.ui.doFastAndSloppyMathCheckBox.setChecked(True)
    self.ui.lowerVoltageCutoffLineEdit.setText('-inf')
    self.ui.lowerVoltageCutoffLineEdit.editingFinished.emit()
    self.ui.upperVoltageCutoffLineEdit.setText('inf')
    self.ui.upperVoltageCutoffLineEdit.editingFinished.emit()
    self.ui.fitMethodComboBox.setCurrentIndex(2)
    self.ui.verbositySpinBox.setValue(0)


  # let's make sure to print messages for the statusbar also in the console    
  def myShowMessage(*args, **kwargs):
    print('Menubar Message:',args[1])
    return args[0].oldShowMessage(*args[1:], **kwargs)

  def exportInterp(self,row):
    thisGraphData = self.ui.tableWidget.item(row,list(self.cols.keys()).index('plotBtn')).data(Qt.UserRole)
    fitX = thisGraphData["fitX"]
    modelY = thisGraphData["modelY"]
    splineY = thisGraphData["splineY"]
    a = np.asarray([fitX, modelY, splineY])
    a = np.transpose(a).astype(float)
    destinationFolder = os.path.join(self.workingDirectory,'exports')
    QDestinationFolder = QDir(destinationFolder)
    if not QDestinationFolder.exists():
      QDir().mkdir(destinationFolder)
    saveFile = os.path.join(destinationFolder,str(self.ui.tableWidget.item(row,list(self.cols.keys()).index('file')).text())+'.csv')
    header = 'Voltage [V],CharEqn Current [mA/cm^2],Spline Current [mA/cm^2]'
    try:
      np.savetxt(saveFile, a, delimiter=",",header=header)
      self.goodMessage()
      self.ui.statusbar.showMessage("Exported " + saveFile,5000)
    except:
      self.badMessage()
      self.ui.statusbar.showMessage("Could not export " + saveFile,self.messageDuration)

  def handleUpperLimChange(self):
    lineEdit = self.sender()
    try:
      self.upperVLim = float(lineEdit.text())
      self.settings.setValue('upperVoltageCutoff',lineEdit.text())
    except:
      pass

  def handleFitMethodChange(self):
    comboBox = self.sender()
    self.settings.setValue('fitMethod',comboBox.currentIndex())

  def handleVerbosityChange(self):
    spinBox = self.sender()
    self.settings.setValue('verbosity',spinBox.value())

  def handleConstraintsChange(self):
    lineEdit = self.sender()
    name = lineEdit.objectName()
    nameSplit = name.split('_')

    try:
      text = lineEdit.text()
      value = float(text)
      if nameSplit[1] == 'lb':
        self.bounds[nameSplit[0]][0] = value
      else: # upper bound
        self.bounds[nameSplit[0]][1] = value
      self.settings.setValue(name,text)
    except:
      pass


  def handleLowerLimChange(self):
    lineEdit = self.sender()
    try:
      self.lowerVLim = float(lineEdit.text())
      self.settings.setValue('lowerVoltageCutoff',lineEdit.text())
    except:
      pass    

  def handleMathChange(self):
    checkBox = self.sender()
    self.settings.setValue('fastAndSloppy',checkBox.isChecked())
    self.symbolCalcsNotDone = True
    if self.ui.attemptCharEqnFitCheckBox.isChecked():
      tp = QThreadPool.globalInstance()
      tp.setMaxThreadCount(1)
      worker = doSymbolicMath(self)
      worker.setAutoDelete(True)
      worker.signals.result.connect(self.handleMathFinished)
      tp.start(worker)

  def handleEqnFitChange(self):
    checkBox = self.sender()
    self.settings.setValue('fitToEqn',checkBox.isChecked())
    self.symbolCalcsNotDone = True
    if checkBox.isChecked():
      tp = QThreadPool.globalInstance()
      tp.setMaxThreadCount(1)
      worker = doSymbolicMath(self)
      worker.setAutoDelete(True)
      worker.signals.result.connect(self.handleMathFinished)
      tp.start(worker)    

  def handleButton(self):
    btn = self.sender()
    #kinda hacky:
    row = self.ui.tableWidget.indexAt(btn.pos()).row()
    col = self.ui.tableWidget.indexAt(btn.pos()).column()
    if col == 0:
      self.rowGraph(row)
    if col == 1:
      self.exportInterp(row)

  def rowGraph(self,row):
    thisGraphData = self.ui.tableWidget.item(row,list(self.cols.keys()).index('plotBtn')).data(Qt.UserRole)
    filename = str(self.ui.tableWidget.item(row,list(self.cols.keys()).index('file')).text())

    v = thisGraphData["v"]
    i = thisGraphData["i"]
    if not thisGraphData["vsTime"]:
      plt.plot(v, i, c='b', marker='o', ls="None",label='I-V Data')
      plt.scatter(thisGraphData["Vmax"], thisGraphData["Imax"], c='g',marker='x',s=100)
      plt.scatter(thisGraphData["Voc"], 0, c='g',marker='x',s=100)
      plt.scatter(0, thisGraphData["Isc"], c='g',marker='x',s=100)
      fitX = thisGraphData["fitX"]
      modelY = thisGraphData["modelY"]
      splineY = thisGraphData["splineY"]
      if not mpmath.isnan(modelY[0]):
        plt.plot(fitX, modelY.astype(complex),c='k', label='CharEqn Best Fit')
      plt.plot(fitX, splineY,c='g', label='Spline Fit')
      plt.autoscale(axis='x', tight=True)
      plt.grid(b=True)
      ax = plt.gca()
      handles, labels = ax.get_legend_handles_labels()
      ax.legend(handles, labels, loc=3)

      plt.annotate(
        thisGraphData["Voc"].__format__('0.4f')+ ' V', 
              xy = (thisGraphData["Voc"], 0), xytext = (40, 20),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

      plt.annotate(
        float(thisGraphData["Isc"]).__format__('0.4f') + ' mA/cm^2', 
              xy = (0,thisGraphData["Isc"]), xytext = (40, 20),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

      plt.annotate(
        float(thisGraphData["Imax"]*thisGraphData["Vmax"]).__format__('0.4f') + '% @(' + float(thisGraphData["Vmax"]).__format__('0.4f') + ',' + float(thisGraphData["Imax"]).__format__('0.4f') + ')', 
              xy = (thisGraphData["Vmax"],thisGraphData["Imax"]), xytext = (80, 40),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))		

      plt.ylabel('Current [mA/cm^2]')
      plt.xlabel('Voltage [V]')
    else: #vs time
      tData = thisGraphData["time"]

      fig, ax1 = plt.subplots()
      ax1.plot(tData, v, 'b-',label='Voltage [V]')
      ax1.set_xlabel('Time [s]')
      # Make the y-axis label and tick labels match the line color.
      ax1.set_ylabel('Voltage [V]', color='b')
      for tl in ax1.get_yticklabels():
        tl.set_color('b')
      #fdsf
      ax2 = ax1.twinx()
      ax2.plot(tData, i, 'r-')
      ax2.set_ylabel('Current [mA/cm^2]', color='r')
      for tl in ax2.get_yticklabels():
        tl.set_color('r')            

    plt.title(filename)
    plt.draw()
    plt.show()       

  # this is how we save the table data to a .csv or .mat file
  def handleSave(self):
    if self.settings.contains('lastFolder'):
      saveDir = self.settings.value('lastFolder')
    else:
      saveDir = '.'
    path = QFileDialog.getSaveFileName(self, caption='Set Export File',filter="Comma separated values (*.csv);;MATLAB formatted data (*.mat)", directory=saveDir)
    if str(path[0]) == '':
      return
    elif '.csv' in str(path[1]): # let's write a .csv
      fullPath = str(path[0])
      if not fullPath.endswith('.csv'):
        fullPath = fullPath + '.csv'            
      with open(fullPath, 'w') as stream:
        writer = csv.writer(stream)
        rowdata = []
        for column in range(self.ui.tableWidget.columnCount()):
          item = self.ui.tableWidget.horizontalHeaderItem(column)
          if item is not None:
            rowdata.append(str(item.text()).replace('\n',' '))
          else:
            rowdata.append(b'')
        writer.writerow(rowdata[2:])                
        for row in range(self.ui.tableWidget.rowCount()):
          rowdata = []
          for column in range(self.ui.tableWidget.columnCount()):
            item = self.ui.tableWidget.item(row, column)
            if item is not None:
              rowdata.append(str(item.text()))
            else:
              rowdata.append('')
          writer.writerow(rowdata[2:])
        stream.close()
        print('Table data successfully written to', fullPath)
    elif '.mat' in str(path[1]):# let's write a .mat file
      fullPath = str(path[0])
      if not fullPath.endswith('.mat'):
        fullPath = fullPath + '.mat'
      #let's make a dict out of the table:
      tableDict = {}

      fieldsToInclude= ('pce_spline','pmax_spline','voc_spline','isc_spline','ff_spline','vmax_spline','SSE','pce_fit','pmax_fit','voc_fit','isc_fit','ff_fit','vmax_fit','rs','rsh','iph','i0','n','area','suns')

      #how many padding zeros should we use for the MATLAB variable names?
      ndigits = str(len(str(self.ui.tableWidget.rowCount()))) 

      for row in range(self.ui.tableWidget.rowCount()):
        rowDict = {}
        rowDict['file'] = self.ui.tableWidget.item(row, list(self.cols.keys()).index('file')).data(Qt.DisplayRole)
        for field in fieldsToInclude:
          rowDict[field] = self.ui.tableWidget.item(row, list(self.cols.keys()).index(field)).data(Qt.UserRole)
        rowDict['i'] = self.ui.tableWidget.item(row, list(self.cols.keys()).index('plotBtn')).data(Qt.UserRole)['i']/1000*rowDict['area']
        rowDict['v'] = self.ui.tableWidget.item(row, list(self.cols.keys()).index('plotBtn')).data(Qt.UserRole)['v']
        tableDict['thing'+format(row, '0'+ndigits)] = rowDict

      # save our dict as a .mat file
      sio.savemat(fullPath, tableDict)
      print('Table data successfully written to', fullPath)

  def formatTableRowForDisplay(self,row):      
    ignoreCols = ['plotBtn','exportBtn','file']
    cols = list(self.cols.keys())
    for coli in range(len(cols)):
      thisCol = cols[coli]
      if thisCol not in ignoreCols:
        value = self.ui.tableWidget.item(row,coli).data(Qt.UserRole)
        if value is not None:
          if thisCol == 'SSE':
            value = value*mWperW**2 # A^2 to mA^2
          elif thisCol in ['ff_spline','ff_fit']:
            value = value*100 # to percent
          elif thisCol in ['jsc_spline','isc_spline','voc_spline','voc_fit','jsc','isc','jph','iph','vmax_spline','vmax_fit','pmax_spline','pmax_fit','pmax_a_spline','pmax_a_fit']:
            value = value*1e3 # to milli-
          elif thisCol in ['area']:
            value = value*1e2 # to centi-
          elif thisCol in ['i0','j0']:
            value = value*1e9 # to nano-

        self.ui.tableWidget.item(row,coli).setData(Qt.DisplayRole,to_precision(value,4))

  def clearTableCall(self):
    for ii in range(self.rows):
      self.ui.tableWidget.removeRow(0)
    self.ui.tableWidget.clearContents()
    self.rows = 0
    self.fileNames = []

  def processFile(self,fullPath):
    result = {}
    logMessages = StringIO()
    result['fullPath'] = fullPath


    fileName, fileExtension = os.path.splitext(fullPath)
    fileName = os.path.basename(fullPath)
    self.fileNames.append(fileName)
    isSnaithFile = False
    if fileExtension == '.csv':
      delimiter = ','
    elif fileExtension == '.tsv':
      delimiter = '\t'        
    else:
      delimiter = None

    print("Processing:", fileName, file = logMessages)

    #wait here for the file to be completely written to disk and closed before trying to read it
    fi = QFileInfo(fullPath)
    while (not fi.isWritable()):
      time.sleep(0.001)
      fi.refresh()

    fp = open(fullPath, mode='r')
    fileBuffer = fp.read()
    fp.close()
    if len(fileBuffer) < 25:
      print('Could not read' + fileName +'. This file is less than 25 characters long.', file = logMessages)
      return
    first10 = fileBuffer[0:10]
    last25 = fileBuffer[-26:-1]

    isMcFile = False #true if this is a McGehee iv file format
    isSnaithFile = False # true if this is a Snaith iv file format
    #mcFile test:
    if (not first10.__contains__('#')) and (first10.__contains__('/')) and (first10.__contains__('\t')):#the first line is not a comment
      nMcHeaderLines = 25 #number of header lines in mcgehee IV file format
      #the first 8 chars do not contain comment symbol and do contain / and a tab, it's safe to assume mcgehee iv file format
      isMcFile = True
      #comment out the first 25 rows here
      fileBuffer = '#'+fileBuffer
      fileBuffer = fileBuffer.replace('\n', '\n#',nMcHeaderLines-1)
    #snaithFile test:
    elif last25.__contains__('suns:\t'):
      nSnaithFooterLines = 11 #number of footer lines in snaith IV file format
      isSnaithFile = True
      delimiter = '\t'
      if (fileExtension == '.liv1') or (fileExtension == '.div1'):
        snaithReverse = True
      if (fileExtension == '.liv2') or (fileExtension == '.div2'):
        snaithReverse = False
      fileBuffer = fileBuffer[::-1] # reverse the buffer
      fileBuffer = fileBuffer.replace('\n', '#\n',nSnaithFooterLines+1) # comment out the footer lines
      fileBuffer = fileBuffer[::-1] # un-reverse the buffer
      fileBuffer = fileBuffer[:-3] # remove the last (extra) '\r\n#'

    splitBuffer = fileBuffer.splitlines(True)

    suns = 1
    area = 1 # in cm^2
    noArea = True
    noIntensity = True
    vsTime = False #this is not an i,v vs t data file
    #extract comments lines and search for area and intensity
    comments = []
    for line in splitBuffer:
      if line.startswith('#'):
        comments.append(line)
        if line.__contains__('Area'):
          numbersHere = [float(s) for s in line.split() if isNumber(s)]
          if len(numbersHere) is 1:
            area = numbersHere[0]
            noArea = False
        elif line.__contains__('I&V vs t'):
          if float(line.split(' ')[5]) == 1:
            vsTime = True
        elif line.__contains__('Number of suns:'):
          numbersHere = [float(s) for s in line.split() if isNumber(s)]
          if len(numbersHere) is 1:
            suns = numbersHere[0]
            noIntensity = False

    jScaleFactor = 1000/area #for converstion to current density[mA/cm^2]

    c = StringIO(fileBuffer) # makes string look like a file 

    #read in data
    try:
      data = np.loadtxt(c,delimiter=delimiter)
    except:
      print('Could not read' + fileName +'. Prepend # to all non-data lines and try again', file = logMessages)
      return
    VV = data[:,0]
    II = data[:,1]
    if isMcFile or isSnaithFile: # convert from current density to amps through soucemeter
      II = II/jScaleFactor

    if vsTime:
      tData = data[:,2]
      # store off the time data in special vectors
      VVt = VV
      IIt = II
      newOrder = tData.argsort()
      VVt=VVt[newOrder]
      IIt=IIt[newOrder]
      tData=tData[newOrder]
      tData=tData-tData[0]#start time at t=0            

    # prune data points that share the same voltage
    u, indices = np.unique(VV, return_index=True)
    VV = VV[indices]
    II = II[indices]

    # sort data by ascending voltage
    newOrder = VV.argsort()
    VV=VV[newOrder]
    II=II[newOrder]

    # trim data to voltage range
    vMask = (VV>self.lowerVLim) & (VV<self.upperVLim)
    VV=VV[vMask]
    II=II[vMask]

    # the task now is to figure out how this data was collected so that we can fix it
    # this is important because the guess and fit algorithms below expect the data to be
    # in a certian way
    # essentially, we want to know if current or voltage has the wrong sign
    # the goal here is that the curve "knee" ends up in quadrant 1
    # (for a light curve, and quadrant 4 for a dark curve)
    smoothingParameter = 1-1e-3
    coefs, brks = findBreaksAndCoefs(VV, II, smoothingParameter)        
    superSmoothSpline = scipy.interpolate.PPoly(coefs,brks)
    superSmoothSplineD1 = superSmoothSpline.derivative(1) # first deravive
    superSmoothSplineD2 = superSmoothSpline.derivative(2) # second deravive

    vv=np.linspace(min(VV),max(VV),1000)
    if vv[np.abs(superSmoothSplineD2(vv)).argmax()] < 0: # fix flipped voltage sign
      VV = VV * -1
      newOrder = VV.argsort()
      II=II[newOrder]
      VV=VV[newOrder]
      vv=np.linspace(min(VV),max(VV),1000)
      print("Flipping voltage sign.", file = logMessages)
    if II[0] < II[-1]:
      II = II * -1
      print("Flipping current sign.", file = logMessages)

    # now let's do a spline fit for our data
    # this prevents measurment noise from impacting our results
    # and allows us to interpolate between data points

    smoothingParameter = 1-2e-6
    #0 -> LS-straight line
    #1 -> cubic spline interpolant

    coefs, brks = findBreaksAndCoefs(VV, II, smoothingParameter)
    smoothSpline = scipy.interpolate.PPoly(coefs,brks)        

    pCoefs, pBrks = findBreaksAndCoefs(VV, II*VV, smoothingParameter)
    powerSpline = scipy.interpolate.PPoly(pCoefs,pBrks)
    powerSplineD1 = powerSpline.derivative(1)



    ##newCoefs = np.vstack((coefs,np.zeros(coefs.shape[1])))
    ##k=4
    ##for i in range(coefs.shape[1]-1):
    ##    newCoefs[4,i+1] = sum(newCoefs[m, i] * (brks[i+1] - brks[i])**(k-m) for m in range(k+1))

    ##newCoefs = np.vstack((coefs,np.zeros(coefs.shape[1]),np.zeros(coefs.shape[1])))
    ##newCoefs = np.vstack((coefs,S))
    ##newCoefs = np.vstack((newCoefs,np.zeros(coefs.shape[1])))
    ##powerSpline = scipy.interpolate.PPoly(newCoefs,brks) # multiply smoothSpline by x        

    ##splineB = scipy.interpolate.UnivariateSpline(VV, II, s=1e-7)
    ##splineC = scipy.interpolate.Rbf(VV, II, smooth=0.05)
    ##coeffsA = scipy.signal.cspline1d(II, lamb=0.5)
    ##coeffsB = scipy.signal.cspline1d(II, smooth=0.1)
    ##splineE = scipy.interpolate.Rbf(VV, II,function="inverse",smooth=2e-5)
    ##II2 = scipy.signal.medfilt(II, kernel_size=5)
    ##splineF = scipy.interpolate.Rbf(VV, II2,function='cubic', smooth=2e-8)
    ##splineG = scipy.interpolate.UnivariateSpline(VV, II2, s=1e-7)

    ##coefs, brks = _compute_coefs(VV, II, smoothingParameter, 1)
    ##pthing = scipy.interpolate.PPoly(coefs,brks)        


    #iiA = smoothSpline(vv)
    #iiB = powerSpline(vv)
    #iiC = smoothSpline(vv)*vv
    #iiD = powerSplineD1(vv)
    ##iiB = tehPPoly2(vv)
    ##iiC = ppder1(vv)
    ##iiD = ppder2(vv)
    ###iiC = scipy.interpolate.splev(vv, tck)
    ###iiC = splineC(vv)
    ###iiD = scipy.signal.cspline1d_eval(coeffsA, vv,dx=dx)
    ###iiE = splineE(vv)
    ###iiF = splineF(vv)
    ###iiG = splineG(vv)
    ###iiH = pthing(vv)

    #plt.title('Spline analysis')
    #p1, = plt.plot(VV,II,ls='None',marker='o', label='Data')
    #p2, = plt.plot(vv,iiA, label='smoothSpline')
    #p3, = plt.plot(vv,iiB, label='powerSpline')
    #p4, = plt.plot(vv,iiC, label='cheating')
    #p5, = plt.plot(vv,iiD, label='powerSplineD1')
    ##p6, = plt.plot(vv,iiD, label='der2')
    ###p7, = plt.plot(vv,iiE, label='Rbf cubic')
    ###p8, = plt.plot(vv,iiF, label='Rbf with medfilt')
    ###p8, = plt.plot(vv,iiG, label='univariat with medfilt')
    ###p8, = plt.plot(vv,iiH, label='rippedOut')
    #ax = plt.gca()
    #handles, labels = ax.get_legend_handles_labels()
    #ax.legend(handles, labels, loc=3)
    #plt.grid(b=True)
    #plt.draw()
    #plt.show()   
    #plt.pause(500)

    #print("done")

    #if any(II>

    # catch and fix flipped current sign
    # The philosoply in use here is that energy producers have positive current defined as flowing out of the positive terminal
    #if II[0] < II[-1]:
    #    self.ui.statusbar.showMessage("Incorrect current convention detected. I'm fixing that for you.",500)
    #    II = II * -1

    #VV = VV * -1# TODO: remove this hack and properly detect inverted devices!
    #II = II * -1
    #indexInQuad1 = np.logical_and(VV>0,II>0)
    #if any(indexInQuad1): # enters statement if there is at least one datapoint in quadrant 1
      #isDarkCurve = False
    #else:
      ## pick out data points in each quadrant
      #indexInQuad2 = np.logical_and(VV<0,II>0) 
      #indexInQuad3 = np.logical_and(VV<0,II<0)
      #indexInQuad4 = np.logical_and(VV>0,II<0)
      ## find the largest powers in each quad
      #if any(indexInQuad2):
        #PP2 = np.min(VV[indexInQuad2]*II[indexInQuad2])
      #else:
        #PP2 = 0
      #if any(indexInQuad3):
        #PP3 = np.max(VV[indexInQuad3]*II[indexInQuad3])
      #else:
        #PP3 = 0
      #if any(indexInQuad4):
        #PP4 = np.max(VV[indexInQuad4]*II[indexInQuad4])
      #else:
        #PP4 = 0

      ## catch and fix flipped voltage polarity(!)
      #if (PP4<(PP2-PP3)):
        #self.ui.statusbar.showMessage("Dark curve detected",500)
        #isDarkCurve = True
      #else:
        ## TODO: dark curves of this messed up nature will likely not be caught
        #self.ui.statusbar.showMessage("Inverted I-V convention detected: I'm fixing that for you.",500)
        #II = II * -1
        #VV = VV * -1
        #newOrder = VV.argsort()
        #VV=VV[newOrder]
        #II=II[newOrder]                
        #isDarkCurve = False

    isDarkCurve = False
    Vmpp = powerSplineD1.roots(extrapolate=False,discontinuity=False)
    VmppSize = Vmpp.size
    if VmppSize is 0:
      Pmpp = nan
      Impp = nan
      Vmpp = nan
      isDarkCurve = True
    elif VmppSize is 1:
      Vmpp = float(Vmpp)
      Impp = float(smoothSpline(Vmpp))
      Pmpp = Impp*Vmpp
    else: # there are more than one local power maxima
      Impp = smoothSpline(Vmpp)
      Pmpp = Impp*Vmpp
      arg = Pmpp.argmax()
      Pmpp = float(Pmpp[arg])
      Vmpp = float(Vmpp[arg])
      Impp = float(Impp[arg])

    Voc = smoothSpline.roots(extrapolate=True,discontinuity=False)
    VocSize = Voc.size
    isDarkCurve = False
    abortTheFit = True
    if VocSize is 0:
      Voc = nan
      isDarkCurve = True
    elif VocSize is 1:
      Voc = float(Voc)
    else: # got too many answers
      valid = np.logical_and(Voc > 0, Voc < max(VV)+0.05)
      nVocs = sum(valid)
      if nVocs !=1:
        print("Warning: we found",nVocs,"values for Voc, using the last one.", file = logMessages)
        Voc = Voc[-1]
        #Voc = nan
      else:
        Voc = float(Voc[valid][0])

    if isDarkCurve:
      print("Dark curve detected.", file = logMessages)

    Isc = float(smoothSpline(0))
    FF = Pmpp/(Voc*Isc)

    # here's how we'll discretize our fits
    plotPoints = 1000
    if min(VV) > 0: # check for only positive voltages
      vvMin = -0.05 # plot at least 50 mV below zero
    else:
      vvMin = min(VV)

    if max(VV) < Voc: # check for data beyond Voc
      vvMax = Voc + 0.05 # plot at least 50 mV above Voc
    else:
      vvMax = max(VV)

    vv = np.linspace(vvMin,vvMax,plotPoints)

    splineY = smoothSpline(vv)*jScaleFactor
    modelY = [nan]
    result['graphData'] = {'vsTime':vsTime,'origRow':self.rows,'fitX':vv,'modelY':modelY,'splineY':splineY,'i':II*jScaleFactor,'v':VV,'Voc':Voc,'Isc':Isc*jScaleFactor,'Vmax':Vmpp,'Imax':Impp*jScaleFactor}

    # put items in table
    ##self.ui.tableWidget.insertRow(self.rows)
    ##for ii in range(len(self.cols)):
    ##  self.ui.tableWidget.setItem(self.rows,ii,QTableWidgetItem())

    # here's how we put data into the table
    ##insert = lambda colName,value: self.ui.tableWidget.item(self.rows,list(self.cols.keys()).index(colName)).setData(Qt.UserRole,float(np.real(value)))
    result['insert'] = {}
    result['insert']['pce_spline'] = (Pmpp/area)/(stdIrridance*suns/sqcmpersqm)*100
    result['insert']['pmax_spline'] = Pmpp/area
    result['insert']['pmax_a_spline'] = Pmpp
    result['insert']['isc_spline'] = Isc
    result['insert']['jsc_spline'] = Isc/area
    result['insert']['voc_spline'] = Voc
    result['insert']['ff_spline'] = FF
    result['insert']['vmax_spline'] = Vmpp
    result['insert']['area'] = area
    result['insert']['suns'] = suns

    #insert('pce_spline',(Pmpp/area)/(stdIrridance*suns/sqcmpersqm)*100)
    #insert('pmax_spline',Pmpp/area)
    #insert('pmax_a_spline',Pmpp)
    #insert('isc_spline',Isc)
    #insert('jsc_spline',Isc/area)
    #insert('voc_spline',Voc)
    #insert('ff_spline',FF)
    #insert('vmax_spline',Vmpp)
    #insert('area',area)
    #insert('suns',suns)

    if not vsTime:
      if not self.ui.attemptCharEqnFitCheckBox.isChecked():
        print("Not attempting fit to characteristic equation.", file = logMessages)
      else:
        # set bounds on the fit variables
        # if upper=lower bound, then that variable will be taken out of the optimization
        #bounds ={}
        #bounds['I0'] = [0, inf] 
        #bounds['Iph'] = [0, inf]
        #bounds['Rs'] = [0, inf]
        #bounds['Rsh'] = [0, inf]
        #bounds['n'] = [0, inf]
        localBounds = self.bounds

        # take a guess at what the fit parameters will be
        #pr.enable()
        #tnot = time.time()
        try:
          guess = makeAReallySmartGuess(VV,II,isDarkCurve)
        except:
          print("Warning: makeAReallySmartGuess() function failed!", file = logMessages)
          guess = {'I0':1e-9, 'Iph':II[0], 'Rs':5, 'Rsh':1e6, 'n':1.0}

        #print (time.time()-tnot)
        #print(len(VV),len(II),VV.mean(),II.mean())
        #pr.disable()
        #s = io.StringIO()
        #sortby = 'cumulative'
        #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        #ps.print_stats()
        #lines = s.getvalue()
        #print(lines.splitlines()[0])
        #print(lines)            


        #localBounds['I0'] = [x*currentScaleFactor for x in localBounds['I0']]
        #localBounds['Iph'] = [x*currentScaleFactor for x in localBounds['Iph']]
        #localBounds['Rs'] = [x/currentScaleFactor for x in localBounds['Rs']]
        #localBounds['Rsh'] = [x/currentScaleFactor for x in localBounds['Rsh']]

        # let's make sure we're not guessing outside the bounds
        if guess['I0'] < localBounds['I0'][0]:
          guess['I0'] = localBounds['I0'][0]
        elif guess['I0'] > localBounds['I0'][1]:
          guess['I0'] = localBounds['I0'][1]

        if guess['Iph'] < localBounds['Iph'][0]:
          guess['Iph'] = localBounds['Iph'][0]
        elif guess['Iph'] > localBounds['Iph'][1]:
          guess['Iph'] = localBounds['Iph'][1]

        if guess['Rs'] < localBounds['Rs'][0]:
          guess['Rs'] = localBounds['Rs'][0]
        elif guess['Rs'] > localBounds['Rs'][1]:
          guess['Rs'] = localBounds['Rs'][1]

        if guess['Rsh'] < localBounds['Rsh'][0]:
          guess['Rsh'] = localBounds['Rsh'][0]
        elif guess['Rsh'] > localBounds['Rsh'][1]:
          guess['Rsh'] = localBounds['Rsh'][1]

        if guess['n'] < localBounds['n'][0]:
          guess['n'] = localBounds['n'][0]
        elif guess['n'] > localBounds['n'][1]:
          guess['n'] = localBounds['n'][1]

        # scale the current up so that the curve fit algorithm doesn't run into machine precision issues
        currentScaleFactor = 1/abs(II.mean())
        #currentScaleFactor = 1e5
        #currentScaleFactor = 1
        guess['I0'] = guess['I0']*currentScaleFactor
        guess['Iph'] = guess['Iph']*currentScaleFactor
        guess['Rs'] = guess['Rs']/currentScaleFactor
        guess['Rsh'] = guess['Rsh']/currentScaleFactor            
        II = II*currentScaleFactor

        #pr.enable()
        try:
          result['fitResult'] = doTheFit(VV,II,guess,localBounds,self)
        except:
          result['fitResult'] = {'success': False, 'message': 'Warning: doTheFit() function crashed!'}

        #fitParams, sigmas, errmsg, status = doTheFit(VV,II,guess,localBounds)
        #{'success':True,'optParams':optimizeResult.x,'sigmas':sigmas,'message':optimizeResult.message}
        #pr.disable()
        #s = io.StringIO()
        #sortby = 'cumulative'
        #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        #ps.print_stats()
        #lines = s.getvalue()
        #print(lines.splitlines()[0])
        #print(lines)

        # unscale the things
        guess['I0'] = guess['I0']/currentScaleFactor
        guess['Iph'] = guess['Iph']/currentScaleFactor
        guess['Rs'] = guess['Rs']*currentScaleFactor
        guess['Rsh'] = guess['Rsh']*currentScaleFactor
        II = II/currentScaleFactor

        if (result['fitResult']['success']):
          print("Good fit because: " + result['fitResult']['message'], file = logMessages)

          params = {}
          params['I0'] = result['fitResult']['optParams'][0]/currentScaleFactor
          params['Iph'] = result['fitResult']['optParams'][1]/currentScaleFactor
          params['Rs'] = result['fitResult']['optParams'][2]*currentScaleFactor
          params['Rsh'] = result['fitResult']['optParams'][3]*currentScaleFactor
          params['n'] = result['fitResult']['optParams'][4]
          SSE = result['fitResult']['SSE']/currentScaleFactor**2


          # TODO: the sigmas are messed up (by scaling?) when doing a l-m fit
          # TODO: re-do uncertanties
          #sigmas = fitResult['sigmas']
          #sigmas[0] = sigmas[0]/currentScaleFactor
          #sigmas[1] = sigmas[1]/currentScaleFactor
          #sigmas[2] = sigmas[2]*currentScaleFactor
          #sigmas[3] = sigmas[3]*currentScaleFactor
          #fitParams[0] = params['I0']

          # this will produce an evaluation of how well the fit worked
          doVerboseAnalysis = False
          if doVerboseAnalysis:
            analyzeGoodness(VV,II,params,guess,result['fitResult']['message'])

          #do error estimation:
          #alpha = 0.05 # 95% confidence interval = 100*(1-alpha)

          #nn = len(VV)    # number of data points
          #p = len(sigmas) # number of parameters

          #dof = max(0, nn - p) # number of degrees of freedom

          # student-t value for the dof and confidence level
          #tval = t.ppf(1.0-alpha/2., dof) 

          #lowers = []
          #uppers = []
          #calculate 95% confidence interval
          #for a, p, sigma in zip(list(range(nn)), fitParams, sigmas):
          #    lower = p - sigma*tval
          #    upper = p + sigma*tval
          #    lowers.append(lower)
          #    uppers.append(upper)
          uppers = [nan,nan,nan,nan,nan]
          lowers = [nan,nan,nan,nan,nan]                

          # force parameter
          #params['Iph'] = 0.00192071

          # find mpp
          VmppGuess = VV[np.array(VV*II).argmax()]
          mppFound = False
          try:
            Vmpp_charEqn = np.complex(sympy.nsolve(P_prime.subs(zip([I0,Iph,Rsh,Rs,n],[params['I0'],params['Iph'],params['Rsh'],params['Rs'],params['n']])), VmppGuess))
            mppFound = True
          except:
            try: # try again with a differnt starting point
              Vmpp_guess = Vmpp_guess-0.1
              Vmpp_charEqn = np.complex(sympy.nsolve(P_prime.subs(zip([I0,Iph,Rsh,Rs,n],[params['I0'],params['Iph'],params['Rsh'],params['Rs'],params['n']])), VmppGuess))
              mppFound = True
            except: # two failures means we're done
              Vmpp_charEqn = nan
          if mppFound:
            Impp_charEqn = slns['I'](I0=params['I0'],Iph=params['Iph'],Rsh=params['Rsh'],Rs=params['Rs'],n=params['n'],V=Vmpp_charEqn)
            Pmpp_charEqn = Impp_charEqn*Vmpp_charEqn
          else:
            Impp_charEqn = nan
            Pmpp_charEqn = nan

          # find Voc
          try:
            Voc_charEqn = Voc_eqn(I0=params['I0'],Iph=params['Iph'],Rsh=params['Rsh'],n=params['n'])
          except:
            Voc_charEqn = nan

          Isc_charEqn = Isc_eqn(I0=params['I0'],Iph=params['Iph'],Rsh=params['Rsh'],Rs=params['Rs'],n=params['n'])
          FF_charEqn = Pmpp_charEqn/(Voc_charEqn*Isc_charEqn)
          result['graphData']['modelY'] = np.array([slns['I'](I0=params['I0'],Iph=params['Iph'],Rsh=params['Rsh'],Rs=params['Rs'],n=params['n'],V=x) for x in vv])*jScaleFactor

          result['insert']['SSE'] = SSE
          result['insert']['rs_a'] = params['Rs']*area
          result['insert']['rs'] = params['Rs']
          result['insert']['rsh_a'] = params['Rsh']*area
          result['insert']['rsh'] = params['Rsh']
          result['insert']['jph'] = params['Iph']/area
          result['insert']['iph'] = params['Iph']
          result['insert']['j0'] = params['I0']/area
          result['insert']['i0'] = params['I0']
          result['insert']['n'] = params['n']
          result['insert']['vmax_fit'] = Vmpp_charEqn
          result['insert']['pmax_fit'] = Pmpp_charEqn
          result['insert']['pmax_a_fit'] = Pmpp_charEqn/area
          result['insert']['pce_fit'] = (Pmpp_charEqn/area)/(stdIrridance*suns/sqcmpersqm)*100
          result['insert']['voc_fit'] = Voc_charEqn
          result['insert']['ff_fit'] = FF_charEqn
          result['insert']['isc_fit'] = Isc_charEqn
          result['insert']['jsc_fit'] = Isc_charEqn/area

          #insert('SSE',SSE)
          #insert('rs_a',params['Rs']*area)
          #insert('rs',params['Rs'])
          #insert('rsh_a',params['Rsh']*area)
          #insert('rsh',params['Rsh'])
          #insert('jph',params['Iph']/area)
          #insert('iph',params['Iph'])
          #insert('j0',params['I0']/area)
          #insert('i0',params['I0'])
          #insert('n',params['n'])
          #insert('vmax_fit',Vmpp_charEqn)
          #insert('pmax_fit',Pmpp_charEqn)
          #insert('pmax_a_fit',Pmpp_charEqn/area)
          #insert('pce_fit',(Pmpp_charEqn/area)/(stdIrridance*suns/sqcmpersqm)*100) 
          #insert('voc_fit',Voc_charEqn)
          #insert('ff_fit',FF_charEqn)
          #insert('isc_fit',Isc_charEqn)
          #insert('jsc_fit',Isc_charEqn/area)          

        else: # fit failure
          print("Bad fit because: " + result['fitResult']['message'],file = logMessages)
          #modelY = np.empty(plotPoints)*nan

    else:#vs time
      print('This file contains time data.')
      result['fitResult']['graphData'] = {'vsTime':vsTime,'origRow':self.rows,'time':tData,'i':IIt*jScaleFactor,'v':VVt}

    ###export button
    ##exportBtn = QPushButton(self.ui.tableWidget)
    ##exportBtn.setText('Export')
    ##exportBtn.clicked.connect(self.handleButton)
    ##self.ui.tableWidget.setCellWidget(self.rows,list(self.cols.keys()).index('exportBtn'), exportBtn)        

    ###file name
    ##self.ui.tableWidget.item(self.rows,list(self.cols.keys()).index('file')).setText(fileName)
    ##self.ui.tableWidget.item(self.rows,list(self.cols.keys()).index('file')).setToolTip(''.join(comments))          

    ###plot button
    ##plotBtn = QPushButton(self.ui.tableWidget)
    ##plotBtn.setText('Plot')
    ##plotBtn.clicked.connect(self.handleButton)
    ##self.ui.tableWidget.setCellWidget(self.rows,list(self.cols.keys()).index('plotBtn'), plotBtn)
    ##self.ui.tableWidget.item(self.rows,list(self.cols.keys()).index('plotBtn')).setData(Qt.UserRole,graphData)


    ##self.formatTableRowForDisplay(self.rows)
    ##self.ui.tableWidget.resizeColumnsToContents()

    ##self.rows = self.rows + 1
    logMessages.seek(0)
    result['logMessages'] = logMessages.read()
    print(result)
    return result

  def openCall(self):
    #remember the last path the user opened
    if self.settings.contains('lastFolder'):
      openDir = self.settings.value('lastFolder')
    else:
      openDir = '.'

    fileNames = QFileDialog.getOpenFileNames(self, directory = openDir, caption="Select one or more files to open", filter = '(*.csv *.tsv *.txt *.liv1 *.liv2 *.div1 *.div2);;Folders (*)')

    if len(fileNames[0])>0:#check if user clicked cancel
      self.workingDirectory = os.path.dirname(str(fileNames[0][0]))
      self.settings.setValue('lastFolder',self.workingDirectory)
      for fullPath in fileNames[0]:
        fullPath = str(fullPath)
        #worker = Worker(self, fullPath)
        #worker.setAutoDelete(True)
        #worker.signals.result.connect(self.processFitResult)
        self.pool.submit(self.processFile,fullPath)
        #self.processFile(fullPath)

      if self.ui.actionEnable_Watching.isChecked():
        watchedDirs = self.watcher.directories()
        self.watcher.removePaths(watchedDirs)
        self.watcher.addPath(self.workingDirectory)
        self.handleWatchUpdate(self.workingDirectory)

  #user chose file --> watch
  def handleWatchAction(self):
    #remember the last path th user opened
    if self.settings.contains('lastFolder'):
      openDir = self.settings.value('lastFolder')
    else:
      openDir = '.'

    myDir = QFileDialog.getExistingDirectory(self,directory = openDir, caption="Select folder to watch")

    if len(myDir)>0:#check if user clicked cancel
      self.workingDirectory = str(myDir)
      self.settings.setValue('lastFolder',self.workingDirectory)
      self.ui.actionEnable_Watching.setChecked(True)
      watchedDirs = self.watcher.directories()
      self.watcher.removePaths(watchedDirs)
      self.watcher.addPath(self.workingDirectory)
      self.handleWatchUpdate(self.workingDirectory)

  #user toggeled Tools --> Enable Watching
  def watchCall(self):
    watchedDirs = self.watcher.directories()
    self.watcher.removePaths(watchedDirs)
    if self.ui.actionEnable_Watching.isChecked():
      if (self.workingDirectory != ''):
        self.watcher.addPath(self.workingDirectory)
        self.handleWatchUpdate(self.workingDirectory)

  def handleWatchUpdate(self,path):
    myDir = QDir(path)
    myDir.setNameFilters(self.supportedExtensions)
    allFilesNow = myDir.entryList()
    allFilesNow = list(allFilesNow)
    allFilesNow = [str(item) for item in allFilesNow]

    differentFiles = list(set(allFilesNow) ^ set(self.fileNames))
    if differentFiles != []:
      for aFile in differentFiles:
        if self.fileNames.__contains__(aFile):
          #TODO: delete the file from the table
          self.ui.statusbar.showMessage('Removed' + aFile,2500)
        else:
          #process the new file
          fullPath = os.path.join(self.workingDirectory,aFile)
          worker = Worker(self, fullPath)
          tp = QThreadPool.globalInstance()
          worker.setAutoDelete(True)
          tp.start(worker)
          #self.processFile(fullPath)

  def statusChanged(self,args):
    if not args:
      # reset the statusbar background
      self.ui.statusbar.setStyleSheet("QStatusBar{padding-left:8px;background:rgba(0,0,0,0);color:black;font-weight:bold;}")

  def goodMessage(self):
    self.ui.statusbar.setStyleSheet("QStatusBar{padding-left:8px;background:rgba(0,128,0,255);color:black;font-weight:bold;}")

  def badMessage(self):
    self.ui.statusbar.setStyleSheet("QStatusBar{padding-left:8px;background:rgba(255,0,0,255);color:black;font-weight:bold;}")