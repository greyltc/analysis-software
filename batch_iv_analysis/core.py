





<<<<<<< HEAD
import mpmath.libmp
assert mpmath.libmp.BACKEND == 'gmpy'
import numpy as np
import sympy
import math
from numpy import nan
from numpy import inf
from numpy import exp
=======
>>>>>>> multiprocessing















<<<<<<< HEAD
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


# needed for findKnotsAndCoefs below
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
    D = scipy.sparse.spdiags(var * np.ones(n), 0, n, n)  # The variance

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
  fitKwargs['max_nfev'] = 20000
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

class MainWindow(QMainWindow):
  workingDirectory = ''
  fileNames = []
  supportedExtensions = ['*.csv','*.tsv','*.txt','*.liv1','*.liv2','*.div1','*.div2']
  bounds ={}
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
=======
>>>>>>> multiprocessing
  

