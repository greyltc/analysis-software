# needed for waiting for initial setup
import time

import os
import sys

import h5py

# to speed this up, we'll use a process pool from here
import concurrent.futures

import itertools
import mpmath.libmp
assert mpmath.libmp.BACKEND == 'gmpy'  # note sagemath can't be installed!
import numpy as np
import sympy

from packaging import version
from numpy import nan
from numpy import exp
from numpy import inf

import warnings

import copy

#sys.setrecursionlimit(10000000)
import dill
#import pickle
#import cloudpickle

import scipy

#to visualize guess
import matplotlib.pyplot as plt

import scipy.io as sio
from io import StringIO

from scipy import odr
from scipy import interpolate
from scipy import optimize
from scipy import special
from scipy.stats.distributions import t #needed for confidence interval calculation #TODO: remove this in favor of uncertainties
#from uncertainties import ufloat #TODO: switch to using this for the error calcs

from lmfit import Model

class Object(object):
  pass

class ivAnalyzer:
  isFastAndSloppy = None # is number crunching faster and less accurate
  symSolutions = None # symbolic solutions for solar cell parameters
  modelSymbols = None
  modelVariables = None
  dillPickle = None # dill pickled solutions for numerical usage later
  slns = None

  sqcmpersqm = 10000 #cm^2 per m^2
  stdIrridance = 1000 #[W/m^2] standard reporting irridance
  mWperW = 1000 # mW per W
  
  multiprocess = None
  readyForAnalysis = False
  
  pool = None
  poolWorkers = None
    
  def __init__(self, beFastAndSloppy=False, poolWorkers=0):
    self.__dict__['poolWorkers'] = poolWorkers

    if poolWorkers == 0:
      self.__dict__['multiprocess'] = False
    else:
      self.__dict__['multiprocess'] = True
    
    self.__dict__['isFastAndSloppy'] = beFastAndSloppy
  
  def __setattr__(self, attr, value):
    if attr == 'isFastAndSloppy':
      self.__dict__[attr] = value
      self.numericalize()
    elif attr == 'poolWorkers':
      self.__dict__[attr] = value
      if (self.multiprocess == True) and (self.pool != None):
        self.buildAPool()
    elif attr == 'multiprocess':
      changed = (self.__dict__[attr]) != value
      if changed:
        self.__dict__[attr] = value
        self.readyForAnalysis = False
        self.numericalize()
        if value == True:
          self.buildAPool()
      if not value:
        try:
          self.pool.shutdown()
        except:
          pass
        try:
          del(self.pool)
        except:
          pass
        self.pool = None
    else:
      self.__dict__[attr] = value
      
  def setup(self):
    print("Multiprocess mode = ",self.multiprocess)
    if self.multiprocess:
      self.buildAPool()
      submission = self.pool.submit(ivAnalyzer.doSymbolicManipulations,self.isFastAndSloppy)
      submission.add_done_callback(self.symbolsDone)
    else:
      results = ivAnalyzer.doSymbolicManipulations(self.isFastAndSloppy)
      self.symbolsDone(results)
      
  def buildAPool(self):
    import concurrent.futures
    if self.pool is not None:
      try:
        self.pool.shutdown()
      except:
        pass
      del(self.pool)
    self.pool = concurrent.futures.ProcessPoolExecutor(max_workers = self.poolWorkers)
      
  def getPoolStatusString(self):
    qDJobs = len(self.pool._pending_work_items)
    processes = len(self.pool._processes)
    activeJobs = len(self.pool._pending_work_items)-self.pool._work_ids.qsize()-self.pool._call_queue.qsize()
    poolStatusString = '[ Pending jobs: ' + str(qDJobs) + ' ]   [ Active jobs: ' + str(activeJobs)+'/' + str(processes) + ' ]'
    return poolStatusString
    
  def symbolsDone(self,results):
    if type(results) is concurrent.futures._base.Future:
      results = results.result()
      
    #print('results type is',type(results))
    #if self.multiprocess:
    #  results = results.result()
    self.symSolutions = results['symSolutions']
    self.modelSymbols = results['modelSymbols']
    self.modelVariables = results['modelVariables']
    self.isFastAndSloppy = results['beFastAndSloppy']
    
  def doSymbolicManipulations(beFastAndSloppy):
    print("Hang tight, we're doing the one-time symbolic manipulations now...")
    
    # let's define some variables we'll use to do some symbolic equaiton manipulation
    modelSymbols = sympy.symbols('I0 Iph Rs Rsh n Vth', rational=True, extended_nonnegative=True) + sympy.symbols('I V', rational=True)
    I0, Iph, Rs, Rsh, n, Vth, I, V = modelSymbols
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
    rhs_no_rs = rhs.subs(Rs, 0)  # with perfect series resistance 
    rhs_no_rsh = rhs.subs(Rsh, sympy.core.numbers.oo)  # with perfect sheet resistance
    rhs_no_r = rhs.subs([(Rsh,sympy.core.numbers.oo),(Rs,0)])  # with no resistive losses at all
    electricalModel = sympy.Eq(lhs,rhs)
    electricalModel_no_rs = sympy.Eq(lhs,rhs_no_rs)
    electricalModel_no_rsh = sympy.Eq(lhs,rhs_no_rsh)
    electricalModel_no_r = sympy.Eq(lhs,rhs_no_r)

    #electricalModelVarsOnly= electricalModel.subs(zip(modelConstants,valuesForConstants))
  
    # symbolically isolate each variable in our characteristic equation
    # then make substitutions for constants
    # then make the solutions ready for use numerically
    # NOTE: this is actually pretty computationally intense;
    # some solutions might contain the Lambert W "function"
    symSolutionsNoSubs = {} # all the symbols preserved
    symSolutionsNoSubs_no_rs = {} # all the symbols preserved
    symSolutionsNoSubs_no_rsh = {} # all the symbols preserved
    symSolutionsNoSubs_no_r = {} # all the symbols preserved
    
    
    solveForThese = [I, I0, V, n]
    for symbol in solveForThese:
      symSolutionsNoSubs[str(symbol)] = sympy.solve(electricalModel,symbol)[0]
      symSolutionsNoSubs_no_rs[str(symbol)] = sympy.solve(electricalModel_no_rs,symbol)[0]
      symSolutionsNoSubs_no_rsh[str(symbol)] = sympy.solve(electricalModel_no_rsh,symbol)[0]
      symSolutionsNoSubs_no_r[str(symbol)] = sympy.solve(electricalModel_no_r,symbol)[0]
      #symSolutionsNoSubs[str(symbol)] = sympy.solveset(electricalModel,symbol,domain=sympy.S.Reals).args[0] # solveset doesn't work here (yet) https://github.com/sympy/sympy/issues/12243
      #symSolutions[str(symbol)] = symSolutionsNoSubs[str(symbol)].subs(zip(modelConstants,valuesForConstants))
      ##remainingVariables = list(set(modelVariables)-set([symbol]))
      ##slns[str(symbol)] = sympy.lambdify(remainingVariables,symSolutions[str(symbol)],functionSubstitutions,dummify=False)
    
    # now we'll solve for some useful device parameters
    #self.symSolutions = symSolutions # analytical solution for all variables
    Voc_eqn = symSolutionsNoSubs['V'].subs(I,0) # analytical solution for Voc
    Voc_eqn_no_rs = symSolutionsNoSubs_no_rs['V'].subs(I,0)
    Voc_eqn_no_rsh = symSolutionsNoSubs_no_rsh['V'].subs(I,0)
    Voc_eqn_no_r = symSolutionsNoSubs_no_r['V'].subs(I,0)
    #Voc_eqn = Voc_eqn.subs(zip(modelConstants,valuesForConstants))
    ##Voc_eqn = sympy.lambdify((I0,Rsh,Iph,n),Voc_eqn,functionSubstitutions,dummify=False)
    Isc_eqn = symSolutionsNoSubs['I'].subs(V,0) # analytical solution for Isc
    Isc_eqn_no_rs = symSolutionsNoSubs_no_rs['I'].subs(V,0)
    Isc_eqn_no_rh = symSolutionsNoSubs_no_rsh['I'].subs(V,0)
    Isc_eqn_no_r = symSolutionsNoSubs_no_r['I'].subs(V,0)
    #Isc_eqn = Isc_eqn.subs(zip(modelConstants,valuesForConstants))
    ##self.Isc_eqn = sympy.lambdify((I0,Rsh,Rs,Iph,n),Isc_eqn,functionSubstitutions,dummify=False)
    PA = symSolutionsNoSubs['I']*V # analytical solution for power (voltage as independant variable)
    P_primeA = sympy.diff(PA,V) # first derivative of power (WRT V)
    PB = symSolutionsNoSubs['V']*I # analytical solution for power (current as independant variable)
    P_primeB = sympy.diff(PB,I) # first derivative of power (WRT I)
    #V_max = sympy.solve(P_prime,V,check=False,implicit=True)[0] # analytical solution for voltage at max power
    
    #V_max = sympy.solveset(P_prime, V, domain=sympy.S.Reals) #TODO: this is not working, but it would be cool...
    #P_max = P.subs(V,V_max)
    #P_max = sympy.lambdify((I0,Rsh,Rs,Iph,n),P_max,functionSubstitutions,dummify=False)
    # since we can't do this analytically (yet) let's try numerically
  
    #sympy.pprint(V_max,use_unicode=True,wrap_line=False)
    #sys.exit(0)
  
    # this puts the symbolic solution for I from above into a format needed for curve_fit
    #I_eqn = autowrap(slns['I'])
    ##I_eqn = lambda x,a,b,c,d,e: np.array([slns['I'](I0=a, Iph=b, Rs=c, Rsh=d, n=e, V=v) for v in x]).astype(complex)
    
    symSolutions = {}
    symSolutions['Isc'] = Isc_eqn.subs(zip(modelConstants,valuesForConstants))
    symSolutions['Voc'] = Voc_eqn.subs(zip(modelConstants,valuesForConstants))
    symSolutions['P_prime'] = P_primeB.subs(zip(modelConstants,valuesForConstants))
    #symSolutions['V_max'] = V_max.subs(zip(modelConstants,valuesForConstants))
    #symSolutions['P_max'] = P_max.subs(zip(modelConstants,valuesForConstants))
    symSolutions['I'] = symSolutionsNoSubs['I'].subs(zip(modelConstants,valuesForConstants))
    symSolutions['I0'] = symSolutionsNoSubs['I0'].subs(zip(modelConstants,valuesForConstants))
    symSolutions['n'] = symSolutionsNoSubs['n'].subs(zip(modelConstants,valuesForConstants))
    symSolutions['V'] = symSolutionsNoSubs['V'].subs(zip(modelConstants,valuesForConstants))
    
    results = {}
    results['symSolutions'] = symSolutions
    results['modelSymbols'] = modelSymbols
    results['modelVariables'] = modelVariables
    results['beFastAndSloppy'] = beFastAndSloppy
    return results
  
  # go from symbolic to numerical domain
  def numericalize (self):
    I0, Iph, Rs, Rsh, n, I, V, Vth = self.modelSymbols
    
    # here we define any function substitutions we'll need for lambdification later
    if self.isFastAndSloppy:
      # for fast and inaccurate math
      functionSubstitutions = {"LambertW" : scipy.special.lambertw, "exp" : np.exp, "log" : np.log}
      #functionSubstitutions = {"LambertW" : scipy.special.lambertw, "exp" : bigfloat.exp}
    else:
      # this is a massive slowdown (forces a ton of operations into mpmath)
      # but gives _much_ better accuracy and aviods overflow warnings/errors...
      functionSubstitutions = {"LambertW" : mpmath.lambertw, "exp" : mpmath.exp, "log" : mpmath.log}
    
    slns = {}
    solveForThese = [I, I0, V, n]
    for symbol in solveForThese:
        remainingVariables = list(set(self.modelVariables)-set([symbol]))
        slns[str(symbol)] = sympy.lambdify(remainingVariables,self.symSolutions[str(symbol)],functionSubstitutions,dummify=False)
        #slns[str(symbol)] = ufuncify(remainingVariables,self.symSolutions[str(symbol)],helpers=[['LambertW', sympy.LambertW(x), [x]]])  
        #slns[str(symbol)] = functools.partial(tmp) 
    
    slns['Voc'] = sympy.lambdify([I0,Rsh,Iph,n],self.symSolutions['Voc'],functionSubstitutions,dummify=False)
    slns['P_prime'] = sympy.lambdify([I0,Rsh,Iph,n,Rs,I],self.symSolutions['P_prime'],functionSubstitutions,dummify=False)
    slns['Isc'] = sympy.lambdify([I0,Rsh,Iph,n,Rs],self.symSolutions['Isc'],functionSubstitutions,dummify=False)
    
    #if not self.isFastAndSloppy:
    #  for key, value in slns.items():
    #    slns[key] = lambda **kw: np.float(value(**kw))
    
    #self.slns = {}
    #self.slns['I'] = functools.partial(tmp['I'],V=V,Iph=Iph,I0=I0,Rsh=Rsh,Rs=Rs,n=n)
    #self.slns['I0'] = functools.partial(tmp['I0'],V=V,I=I,Iph=Iph,Rsh=Rsh,Rs=Rs,n=n)
    #def I (Iph,I0,Rsh,Rs,v,n):
    #  return np.real_if_close(tmp['I'](V=v,Rsh=Rsh,n=n,Iph=Iph,Rs=Rs,I0=I0))
    #setattr(self.sols,'I',I)
    #self.I = I
    #
    #def I0 (i,Iph,Rsh,Rs,v,n):
    #  return np.real_if_close(slns['I0'](Rsh=Rsh,n=n,Iph=Iph,Rs=Rs,I=i,V=v))
    #self.I0 = I0
    
    if self.multiprocess:
      if not self.isFastAndSloppy:
        print("Error: I don't know how to pickle mpmath things.")
        self.readyForAnalysis = False
        return
      else:
        self.dillPickle = dill.dumps(slns)
    else:
      self.slns = slns
    self.readyForAnalysis = True
    print('Ready for analysis. F&S mode =',self.isFastAndSloppy)
    #return slns

  #def I (self,Iph,I0,Rsh,Rs,v,n):
  #  return np.real_if_close(slns['I'](Rsh=Rsh,n=n,Iph=Iph,Rs=Rs,I0=I0,V=v))

  #def I0 (self,i,Iph,Rsh,Rs,v,n):
  #  return np.real_if_close(slns['I0'](Rsh=Rsh,n=n,Iph=Iph,Rs=Rs,I=i,V=v))

  def printResults(results):
    if type(results) is concurrent.futures._base.Future:
      results = results.result()
    print(results)
  
  # paths = list of paths to files to be processed
  # params = list of associated analysis parameters
  # returnCall = function handle to call when the analysis is done
  # this function should be ready to be passed one argument
  # that argument could bethe analysis result dict
  # or a future object where the analysis result dict can be recovered with .result()
  def processFiles(self, paths, returnCall, prepCall):
    if type(paths) is not list:
      paths = [paths]
      
    tic = time.time()
    while not self.readyForAnalysis:
      time.sleep(0.1)
      if (time.time() - tic) > 10:
        print("Error: 10 seconds have passed and we're not ready yet. Aborting.")
        return
    
    futures = []
    for fullPath in paths:
      fileDatas = ivAnalyzer._loadFile(fullPath)
      
      # turn none into empty list
      if fileDatas == None:
        fileDatas = []
      
      #TODO: possibly load multiple curves per file here
      
      for fileData in fileDatas:
        VV = fileData.VV
        II = fileData.II
        vsTime = fileData.vsTime
        thisParams =  prepCall(fullPath, fileData)
        
        if self.multiprocess:
          futures.append(self.pool.submit(ivAnalyzer.processCurve, VV, II, vsTime, thisParams, self.dillPickle, fullPath))
          futures[-1].add_done_callback(returnCall)
        else:
          result = ivAnalyzer.processCurve(VV, II, vsTime, thisParams, self.slns, fullPath)
          returnCall(result)
    
  def _loadFile(fullPath):

    logMessages = StringIO()
    fileName, fileExtension = os.path.splitext(fullPath)
    basename = os.path.basename(fullPath)
    
    # what we'll return
    ret_list = []
    ret = Object()
    
    isMcFile = False #true if this is a McGehee iv file format
    isSnaithLegacyFile = False # true if this is a Snaith legacy iv file format
    isSnaithFile = False # true if this is a Snaith iv file format
    isMyFile = False # true if this is a custom solar sim iv file format
    isH5 = False  #true when this is an hdf5 file
    isNextTsv = False
    
    print("Processing:", fileName, file = logMessages)
    if h5py.is_hdf5(fullPath):  #  hdf5 file processing route
      isH5 = True
      try:
        h5 = h5py.File(fullPath, 'r')
      except Exception as inst:
        print("Failed to open {:} with error:".format(fullPath))
        print(inst)
        return
      
      fr = h5.attrs['Format Revision']
      if hasattr(fr, 'decode'):
        h5rev = h5.attrs['Format Revision'].decode()
      else:
        h5rev = str(h5.attrs['Format Revision'])
      print("Found HDF5 solar sim data format revision {:s} data file".format(h5rev))
      
      this_ver = version.parse(h5rev)
      
      for substrate_str in list(h5.keys()):
        substrate = h5['/'+substrate_str]
        for pixel_str in list(substrate.keys()):
          ret = Object()
          pixel = substrate[pixel_str]
          
          ret.substrate = substrate_str
          ret.pixel = pixel_str
          if 'Intensity [suns]' in h5.attrs:
            ret.suns = h5.attrs['Intensity [suns]']
          else:
            ret.sunsA = h5.attrs['Diode 1 intensity [suns]']
            ret.sunsB = h5.attrs['Diode 2 intensity [suns]']
            ret.suns =  (ret.sunsA + ret.sunsB)/2 # TODO: use the correct diode intensity for specific pixels instead of averaging the two diodes
          
          if this_ver <= version.parse('1.1.0'):
            ret.area = float(pixel.attrs['area']) / 1e4 # in m^2
          else:
            ret.area = pixel.attrs['area']
          ret.vsTime = False
          
          if 'ssPmax' in pixel.attrs:
            ret.ssPmax = pixel.attrs['ssPmax']
          
          if 'Voc' in pixel.attrs:
            ret.Voc = pixel.attrs['Voc']

          if 'Isc' in pixel.attrs:
            ret.Isc = pixel.attrs['Isc']

          if 'Vmpp' in pixel.attrs:
            ret.Vmpp = pixel.attrs['Vmpp']

          if 'Impp' in pixel.attrs:
            ret.Impp = pixel.attrs['Impp']
            
          # this is all the i-v data
          if 'all_measurements' in pixel:
            iv_data = pixel['all_measurements']
          else:
            print('WARNING: Could not find any data in {:}'.format(fileName), file = logMessages)
            return
          
          if 'V_oc dwell' in iv_data.attrs:
            ret.ssVoc = iv_data[iv_data.attrs['V_oc dwell']]

          if 'I_sc dwell' in iv_data.attrs:
            ret.ssIsc = iv_data[iv_data.attrs['I_sc dwell']]
            
          if 'MPPT' in iv_data.attrs:
            ret.mppt = iv_data[iv_data.attrs['MPPT']]       

          # now we pick out regions of interest from the big i-v data set
          if 'Snaith' in iv_data.attrs:
            snaith_region = iv_data[iv_data.attrs['Snaith']]  # I_sc --> V_oc sweep
            ret.VV = np.array([e[0] for e in snaith_region])
            ret.II = np.array([e[1] for e in snaith_region])
            ret.reverseSweep = False
            ret_list.append(copy.deepcopy(ret))
          if 'Sweep' in iv_data.attrs:
            sweep_region = iv_data[iv_data.attrs['Sweep']]  #  V_oc --> I_sc sweep
            if ret_list != []:
              ret_list.append(copy.deepcopy(ret_list[-1]))
            else:
              ret_list.append(copy.deepcopy(ret))
            ret_list[-1].VV = np.array([e[0] for e in sweep_region])
            ret_list[-1].II = np.array([e[1] for e in sweep_region])
            ret_list[-1].reverseSweep = True
      
    else:  # (legacy) non-h5py file
      ret.reverseSweep = False
      if fileExtension == '.csv':
        delimiter = ','
      elif fileExtension == '.tsv':
        delimiter = '\t'
      else:
        delimiter = None

      fp = open(fullPath, mode='r')
      fileBuffer = fp.read()
      fp.close()
      min_length = 800  # in chars
      if len(fileBuffer) < min_length:
        print('Could not read' + fileName +'. This file is less than {:} characters long.'.format(min_length), file = logMessages)
        return
      head_size = 10  # in chars
      head = fileBuffer[0:head_size]
      tail = fileBuffer[-min_length:-1]

      splitlines = fileBuffer.splitlines(True)
    
      #mcFile test:
      if (not head.__contains__('#')) and (head.__contains__('/')) and (head.__contains__('\t')):#the first line is not a comment
        nMcHeaderLines = 25 #number of header lines in mcgehee IV file format
        #the first 8 chars do not contain comment symbol and do contain / and a tab, it's safe to assume mcgehee iv file format
        isMcFile = True
        #comment out the first 25 rows here
        fileBuffer = '#'+fileBuffer
        fileBuffer = fileBuffer.replace('\n', '\n#',nMcHeaderLines-1)
      #snaithFile test:
      elif 'suns:\t' in tail:
        nSnaithFooterLines = 11 #number of footer lines in snaith IV file format
        isSnaithLegacyFile = True
        delimiter = '\t'
        if (fileExtension == '.liv1') or (fileExtension == '.div1'):
          ret.reverseSweep = True
        if (fileExtension == '.liv2') or (fileExtension == '.div2'):
          ret.reverseSweep = False
        fileBuffer = fileBuffer[::-1] # reverse the buffer
        fileBuffer = fileBuffer.replace('\n', '#\n',nSnaithFooterLines+1) # comment out the footer lines
        fileBuffer = fileBuffer[::-1] # un-reverse the buffer
        fileBuffer = fileBuffer[:-3] # remove the last (extra) '\r\n#'
      elif 'i-v file' in head:
        isMyFile = True
      elif ('(# suns)\t' in tail) or ('(# suns):\t' in tail):
        isSnaithFile = True
        footerLines = 21
        i_col = 2
        v_col = 1
        if 'concurrent\t' in tail:  # newest snaith
          footerLines = 45
          i_col = 3
          v_col = 2
        delimiter = '\t'
        if (fileExtension == '.liv1') or (fileExtension == '.div1'):
          ret.reverseSweep = True
        if (fileExtension == '.liv2') or (fileExtension == '.div2'):
          ret.reverseSweep = False
        fileBuffer = '#' + fileBuffer
        fileBuffer = fileBuffer[::-1] # reverse the buffer
        fileBuffer = fileBuffer.replace('\n', '#\n',footerLines+1) # comment out the footer lines
        fileBuffer = fileBuffer[::-1] # un-reverse the buffer
        fileBuffer = fileBuffer[:-3] # remove the last (extra) '\r\n#'
      elif 'status' in splitlines[0] and '.liv' in basename:  # new light tsv
        isNextTsv = True
        i_col = 1
        v_col = 0
      else:
        print(f"Warning: Couldn't parse file: {fileName}")
        return
    
      splitBuffer = fileBuffer.splitlines(True)
    
      ret.substrate = '?'
      ret.pixel = '?'      
      ret.suns = 1
      ret.area = 1 * 1e-4 # in m^2
      ret.vsTime = False #this is not an i,v vs t data file
      #extract comments lines and search for area and intensity
      comments = []
      for line in splitBuffer:
        if line.startswith('#'):
          comments.append(line)
          if 'Area' in line:
            numbersHere = [float(s) for s in line.split() if ivAnalyzer.isNumber(s)]
            if len(numbersHere) == 1:
              ret.area = numbersHere[0] * 1e-4
          elif 'I&V vs t' in line:
            if float(line.split(' ')[5]) == 1:
              ret.vsTime = True
          elif 'Number of suns:' in line or '(# suns)' in line:
            numbersHere = [float(s) for s in line.split() if ivAnalyzer.isNumber(s)]
            if len(numbersHere) == 1:
              ret.suns = numbersHere[0]
          elif line.startswith('#Pixel'):
            splitted = line.split('\t')
            ret.pixel = splitted[1].strip()
          elif line.startswith('#Position'):
            splitted = line.split('\t')
            ret.substrate = splitted[1].upper().strip()

    
      jScaleFactor = 1000/ (ret.area*1e4) #for converstion to current density[mA/cm^2]
    
      c = StringIO(fileBuffer) # makes string look like a file

      skiprows = 0
      if isNextTsv == True:
        skiprows = 1
    
      #read in data
      try:
        data = np.loadtxt(c,delimiter=delimiter, skiprows=skiprows)
      except:
        print('Could not read' + fileName +'. Prepend # to all non-data lines and try again', file = logMessages)
        return
      if data.size == 0:
        print('WARNING: Could not find any data in {:}'.format(fileName), file = logMessages)        
        return
      if isMyFile:
        ret.VV = data[:,2]
        ret.II = data[:,3]
      elif isSnaithFile or isNextTsv:
        ret.VV = data[:,v_col]
        ret.II = data[:,i_col]
      else:
        ret.VV = data[:,0]
        ret.II = data[:,1]
      if isMcFile or isSnaithLegacyFile: # convert from current density to amps through soucemeter
        ret.II = ret.II/jScaleFactor

      if isNextTsv:
        try:
          cur = data[1,1]
          dens = data[1,4]
          ret.area = cur/dens/10 # in m^2
          fns = basename.split('_')
          ret.pixel = fns[-2][-1]
          ret.substrate = fns[-4]
        except Exception:
          print(f"Warning: Couldn't b_parse file: {fileName}")

      ret_list.append(ret)
      
    
    for i in range(len(ret_list)):
      # prune data points that share the same voltage
      u, indices = np.unique(ret_list[i].VV, return_index=True)
      ret_list[i].VV = ret_list[i].VV[indices]
      ret_list[i].II = ret_list[i].II[indices]
    
      # sort data by ascending voltage
      newOrder = ret_list[i].VV.argsort()
      ret_list[i].VV=ret_list[i].VV[newOrder]
      ret_list[i].II=ret_list[i].II[newOrder]
    
    return ret_list
  
  def _doSplineStuff(VV,II):
    logMessages = StringIO()    
    # the task now is to figure out how this data was collected so that we can fix it
    # this is important because the guess and fit algorithms below expect the data to be
    # in a certian way
    # essentially, we want to know if current or voltage has the wrong sign
    # the goal here is that the curve "knee" ends up in quadrant 1
    # (for a light curve, and quadrant 4 for a dark curve)
    smoothingParameter = 1-1e-3
    coefs, brks = ivAnalyzer.findBreaksAndCoefs(VV, II, smoothingParameter)        
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
  
    coefs, brks = ivAnalyzer.findBreaksAndCoefs(VV, II, smoothingParameter)
    smoothSpline = scipy.interpolate.PPoly(coefs,brks)        
  
    pCoefs, pBrks = ivAnalyzer.findBreaksAndCoefs(VV, II*VV, smoothingParameter)
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
    if VmppSize == 0:
      Pmpp = nan
      Impp = nan
      Vmpp = nan
      isDarkCurve = True
    elif VmppSize == 1:
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
    abortTheFit = True
    if VocSize == 0: # never crosses zero, must be dark curve
      Voc = nan
      isDarkCurve = True
    elif VocSize == 1:
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
    #FF = Pmpp/(Voc*Isc)
  
    # here's how we'll discretize our fits
    plotPoints = 1000
    if min(VV) >= 0: # check for non-negative
      vvMin = -0.05 # plot at least 50 mV below zero
    else:
      vvMin = min(VV)
  
    if max(VV) < Voc: # check for data beyond Voc
      vvMax = Voc + 0.05 # plot at least 50 mV above Voc
    else:
      vvMax = max(VV)
  
    vv = np.linspace(vvMin,vvMax,plotPoints)
  
    #splineY = smoothSpline(vv)*jScaleFactor
    ##result['graphData'] = {'vsTime':vsTime,'fitX':vv,'modelY':modelY,'splineY':splineY,'i':II*jScaleFactor,'v':VV,'Voc':Voc,'Isc':Isc*jScaleFactor,'Vmax':Vmpp,'Imax':Impp*jScaleFactor}
  
    # put items in table
    ##self.ui.tableWidget.insertRow(self.rows)
    ##for ii in range(len(self.cols)):
    ##  self.ui.tableWidget.setItem(self.rows,ii,QTableWidgetItem())
  
    # here's how we put data into the table
    ##insert = lambda colName,value: self.ui.tableWidget.item(self.rows,list(self.cols.keys()).index(colName)).setData(Qt.UserRole,float(np.real(value)))
    #result['insert'] = {}
    #result['insert']['pce_spline'] = (Pmpp/area)/(self.stdIrridance*suns/self.sqcmpersqm)*100
    #result['insert']['pmax_spline'] = Pmpp/area
    #result['insert']['pmax_a_spline'] = Pmpp
    #result['insert']['isc_spline'] = Isc
    #result['insert']['jsc_spline'] = Isc/area
    #result['insert']['voc_spline'] = Voc
    #result['insert']['ff_spline'] = FF
    #result['insert']['vmax_spline'] = Vmpp
    #result['insert']['area'] = area
    #result['insert']['suns'] = suns
  
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
    
    ret = Object()
    ret.voltageData = VV # raw voltage data points
    ret.currentData = II # raw current data points
    ret.analyticalVoltage = vv # voltage vlaues for analytical purposes
    ret.splineCurrent = smoothSpline(vv) # current values from spline fit
    ret.Pmpp = Pmpp # power at max power point
    ret.Vmpp = Vmpp # voltage of maximum power point
    ret.Isc = Isc # short circuit current
    ret.Voc = Voc # open circuit voltage
    ret.isDarkCurve = isDarkCurve # dark curve detection flag
    logMessages.seek(0)
    ret.logMessages = logMessages.read() 
    return ret
  
  #def processFile(fullPath, params, s):    
    #fileData = ivAnalyzer._loadFile(fullPath)
    
    #VV = fileData.VV
    #II = fileData.II
    #vsTime = fileData.vsTime
    #suns = fileData.suns
    #area = fileData.area
    
    #return ivAnalyzer.processCurve(VV, II, vsTime, suns, area, params, s, fullPath)
    
  def processCurve(VV, II, vsTime, params, s, fullPath):
    result = {}
    ret = Object()
    logMessages = StringIO()
    result['fullPath'] = fullPath
    ret.params = params
    
    fileName, fileExtension = os.path.splitext(fullPath)
    fileName = os.path.basename(fullPath)
    result['fileName'] = fileName    

    splineData = ivAnalyzer._doSplineStuff(VV, II)
    VV = splineData.voltageData
    II = splineData.currentData
    vv = splineData.analyticalVoltage
    isDarkCurve = splineData.isDarkCurve
    
    # trim data to voltage range
    vMask = (VV > params['lowerVLim']) & (VV < params['upperVLim'])
    VV = VV[vMask]
    II = II[vMask]    

    if not vsTime:
      if not params['doFit']:
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
        
  
        # take a guess at what the fit parameters will be
        #pr.enable()
        #tnot = time.time()
        if type(s) is not dict:
          s = dill.loads(s) # multiprocess case
          
        try:
          guess = ivAnalyzer.makeAReallySmartGuess(VV, II, isDarkCurve, s['I'], s['I0'], s['n'])
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
        
        localBounds = params['bounds']
  
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
        #currentScaleFactor = 1/abs(II.mean())
        #currentScaleFactor = 1e5
        currentScaleFactor = 1
        guess['I0'] = guess['I0']*currentScaleFactor
        guess['Iph'] = guess['Iph']*currentScaleFactor
        guess['Rs'] = guess['Rs']/currentScaleFactor
        guess['Rsh'] = guess['Rsh']/currentScaleFactor            
        II = II*currentScaleFactor
  
        #pr.enable()
        try:
          result['fitResult'] = ivAnalyzer.doTheFit(VV, II,s['I'], guess, localBounds, method = params['method'], verbose = params['verbose'])
        except:      
          result['fitResult'] = {'success': False, 'message': 'Warning: doTheFit() function crashed: '+ str(sys.exc_info()[0])}
  
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
  
          p = {}
          p['I0'] = result['fitResult']['optParams']['I0']/currentScaleFactor
          p['Iph'] = result['fitResult']['optParams']['Iph']/currentScaleFactor
          p['Rs'] = result['fitResult']['optParams']['Rs']*currentScaleFactor
          p['Rsh'] = result['fitResult']['optParams']['Rsh']*currentScaleFactor
          p['n'] = result['fitResult']['optParams']['n']
          SSE = result['fitResult']['chi-square']
  
  
          # TODO: the sigmas are messed up (by scaling?) when doing a l-m fit
          # TODO: re-do uncertanties
          #sigmas = fitResult['sigmas']
          #sigmas[0] = sigmas[0]/currentScaleFactor
          #sigmas[1] = sigmas[1]/currentScaleFactor
          #sigmas[2] = sigmas[2]*currentScaleFactor
          #sigmas[3] = sigmas[3]*currentScaleFactor
          #fitParams[0] = p['I0']
  
          # this will produce an evaluation of how well the fit worked
          doVerboseAnalysis = False
          if doVerboseAnalysis:
            ivAnalyzer.analyzeGoodness(VV,II, s['I'],p,guess,result['fitResult']['message'])
  
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
          #p['Iph'] = 0.00192071
          
          # find Voc
          try:
            Voc_charEqn = np.real_if_close(np.float(s['Voc'](I0=p['I0'],Iph=p['Iph'],Rsh=p['Rsh'],n=p['n'])))
            if not np.isfinite(Voc_charEqn):
              Inew = lambda V: np.float(np.real(s['I'](n=p['n'],I0=p['I0'],Iph=p['Iph'],Rsh=p['Rsh'],Rs=p['Rs'],V=V)))
              sol = scipy.optimize.root(Inew,1)
              if sol.success:
                Voc_charEqn = sol.x[0]
              else:
                Voc_charEqn = nan
          except:
            Voc_charEqn = nan
            
          # find Isc
          try:
            Isc_charEqn = np.real_if_close(np.float(s['Isc'](I0=p['I0'],Iph=p['Iph'],Rsh=p['Rsh'],Rs=p['Rs'],n=p['n'])))
            if not np.isfinite(Isc_charEqn):
              Vnew = lambda I: np.float(np.real(s['V'](n=p['n'],I0=p['I0'],Iph=p['Iph'],Rsh=p['Rsh'],Rs=p['Rs'],I=I)))
              sol = scipy.optimize.root(Vnew,1e-3)
              if sol.success:
                Isc_charEqn = sol.x[0]
              else:
                Isc_charEqn = nan
          except:
            Isc_charEqn = nan
            
          def findMPP(ImppGuess,P_prime,p):
            try:
              Pnew = lambda I: np.float(np.real(P_prime(n=p['n'],I0=p['I0'],Iph=p['Iph'],Rsh=p['Rsh'],Rs=p['Rs'],I=I)))
              sol = scipy.optimize.root(Pnew,ImppGuess)            
              if sol.success:
                Impp = sol.x[0]
              else:
                Impp = nan
            except:
              Impp = nan
            
            return(Impp)
              
            
  
          # find mpp
          ImppGuess = II[np.array(VV*II).argmax()]
          Impp_charEqn = findMPP(ImppGuess, s['P_prime'], p)
          if np.isnan(Impp_charEqn):
            Impp_charEqn = findMPP(ImppGuess-1e-4, s['P_prime'], p)
          
          if not np.isnan(Impp_charEqn):
            Vmpp_charEqn = np.float(np.real_if_close(s['V'](n=p['n'],I0=p['I0'],Iph=p['Iph'],Rsh=p['Rsh'],Rs=p['Rs'],I=Impp_charEqn)))
            Pmpp_charEqn = Impp_charEqn*Vmpp_charEqn
          else:
            Vmpp_charEqn = nan
            Pmpp_charEqn = nan            

          FF_charEqn = Pmpp_charEqn/(Voc_charEqn*Isc_charEqn)
          
          ret.eqnCurrent = np.array([np.real_if_close(s['I'](n=p['n'],I0=p['I0'],Iph=p['Iph'],Rsh=p['Rsh'],Rs=p['Rs'],V=x)) for x in vv])
          ret.sse = SSE
          ret.n = p['n']
          ret.rs = p['Rs']
          ret.rsh = p['Rsh']
          ret.i0 = p['I0']
          ret.iph = p['Iph']
          ret.pmax_fit = Pmpp_charEqn
          ret.isc_fit = Isc_charEqn
          ret.voc_fit = Voc_charEqn
          ret.vmax_fit = Vmpp_charEqn
          
          #result['graphData'] = {}
          #result['graphData']['modelY'] = *jScaleFactor
          
          #result['graphData']['modelY'] = np.array([slns['I'](I0=p['I0'],Iph=p['Iph'],Rsh=p['Rsh'],Rs=p['Rs'],n=p['n'],V=x) for x in vv])*jScaleFactor
  
          result['insert'] = {}
          result['insert']['SSE'] = SSE
          
          result['insert']['rs'] = p['Rs']
          result['insert']['rsh'] = p['Rsh']
          result['insert']['iph'] = p['Iph']
          result['insert']['i0'] = p['I0']
          result['insert']['n'] = p['n']
          result['insert']['vmax_fit'] = Vmpp_charEqn
          result['insert']['pmax_fit'] = Pmpp_charEqn
          result['insert']['voc_fit'] = Voc_charEqn
          result['insert']['ff_fit'] = FF_charEqn
          result['insert']['isc_fit'] = Isc_charEqn
  
          #insert('SSE',SSE)
          #insert('rs_a',p['Rs']*area)
          #insert('rs',p['Rs'])
          #insert('rsh_a',p['Rsh']*area)
          #insert('rsh',p['Rsh'])
          #insert('jph',p['Iph']/area)
          #insert('iph',p['Iph'])
          #insert('j0',p['I0']/area)
          #insert('i0',p['I0'])
          #insert('n',p['n'])
          #insert('vmax_fit',Vmpp_charEqn)
          #insert('pmax_fit',Pmpp_charEqn)
          #insert('pmax_a_fit',Pmpp_charEqn/area)
          #insert('pce_fit',(Pmpp_charEqn/area)/(stdIrridance*suns/sqcmpersqm)*100) 
          #insert('voc_fit',Voc_charEqn)
          #insert('ff_fit',FF_charEqn)
          #insert('isc_fit',Isc_charEqn)
          #insert('jsc_fit',IscPmpp_charEqn/area)          
  
        else: # fit failure
          print("Bad fit because: " + result['fitResult']['message'],file = logMessages)
          #modelY = np.empty(plotPoints)*nan
  
    else:#vs time
      print('This file contains time data.')
      result['fitResult']['graphData'] = {'vsTime':vsTime,'time':tData,'i':IIt,'v':VVt}
  
    ###export button
    ##exportBtn = QPushButton(self.ui.tableWidget)
    ##exportBtn.setText('Export')
    ##exportBtn.clicked.connect(self.handleButton)
    ##self.ui.tableWidget.setCellWidget(self.rows,list(self.cols.keys()).index('exportBtn'), exportBtn)        
  
    ###file name
    ##self.ui.tableWidget.item(self.rows,list(self.cols.keys()).index('file')).setText(fileName)
    ##self.ui.tableWidget.item(self.rows,list(self.cols.keys()).index('file')).setToolTip(''.join(comments))          
  

  
    ##self.formatTableRowForDisplay(self.rows)
    ##self.ui.tableWidget.resizeColumnsToContents()
  
    ##self.rows = self.rows + 1
    
    #print(result)pmax_a_spline
    
    #ret.pce = (splineData.Pmpp/fileData.area)/(ivAnalyzer.stdIrridance*fileData.suns/ivAnalyzer.sqcmpersqm)
    ret.pmpp = splineData.Pmpp # power at max power point
    ret.vmpp = splineData.Vmpp
    ret.isc = splineData.Isc
    ret.voc = splineData.Voc
    

    ret.v = splineData.voltageData
    ret.i = splineData.currentData
    ret.x = splineData.analyticalVoltage
    ret.splineCurrent = splineData.splineCurrent
    logMessages.seek(0)
    ret.logMessages = logMessages.read()    
    return ret
    #return result
  
  # tests if string is a number
  def isNumber(s):
    try:
      float(s)
    except ValueError:
      return False
    return True


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
      D = scipy.sparse.spdiags(var * np.ones(n), 0, n, n)  # The varianceStringIO
  
      u, p = ivAnalyzer._compute_u(p, D, dydx, dx, dx1, n)
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
  
  def lineFit (xData, yData, mGuess, bGuess):
    lineResiduals = lambda p,dataX,dataY: np.square([p[0]*x + p[1] - y for x,y in zip(dataX,dataY)])
    p0 = [mGuess, bGuess]
    fitArgs = (lineResiduals,p0)
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
    #fitKwargs['x_scale'] = 'jac'
    fitKwargs['loss'] = 'soft_l1'
    #fitKwargs['loss'] = 'arctan'
    #fitKwargs['loss'] = 'linear'
    #fitKwargs['f_scale'] = 100000.0
    #fitKwargs['diff_step'] = list(map(lambda x: x/1000000, x0))
    #fitKwargs['diff_step'] = None
    #fitKwargs['tr_solver'] = 'lsmr'
    fitKwargs['tr_solver'] = None
    #fitKwargs['tr_options'] = {'regularize':True}
    #fitKwargs['jac_sparsity'] = None
    fitKwargs['max_nfev'] = 20000
    fitKwargs['verbose'] = 0
    #fitKwargs['method'] = method
    fitKwargs['method'] = 'trf'
    if fitKwargs['method'] == 'lm':
      fitKwargs['loss'] = 'linear' # loss must be linear for lm method
    fitKwargs['args'] = (xData,yData)
    fitKwargs['kwargs'] = {}
  
    # do the fit
    optimizeResult = scipy.optimize.least_squares(*fitArgs,**fitKwargs)
    
    result={}
    if optimizeResult.success:
      return optimizeResult.x
    else:
      print('Warning: Line fit fail!')
      return p0
  
  # this function makes ultra-super-intelligent guesses for all the parameters
  # of the equation that w're about to attempt to fit so that
  # the (relatively dumb and fragile) final optimization/curve fitting routine
  # has the best chance of giving good results
  def makeAReallySmartGuess(VV, II, isDarkCurve, fI, fI0, fn):
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
      vZeroCurrent = iFit(0).item()
    except:
      print("Warning. You really should have some negative voltages in your data...")
      vZeroCurrent = I_start_n
    guess['Iph'] = vZeroCurrent      
  
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
    xData = VV[start_i:vp_i]
    yData = II[start_i:vp_i]
    result = ivAnalyzer.lineFit(xData, yData, -1/guess['Rsh'], guess['Iph'])
    guess['Rsh'] = -1/result[0]
    guess['Iph'] = result[1]
  
    # try to further refine guess for Rs
    xData = VV[ip_i:end_i]
    yData = II[ip_i:end_i]
    result = ivAnalyzer.lineFit(xData, yData, -1/guess['Rs'], V_end_n/guess['Rs'])
    guess['Rs'] = -1/result[0]
    RsYInter = result[1]
      
    # try to refine guess for n
    #starti = ip_i
    starti = int(round((vp_i+knee_i)/2))
    endi = int(round((ip_i+knee_i)/2))
    yData = np.log(np.abs(II[starti:endi]-vZeroCurrent))
    xData = VV[starti:endi]
    goodis = np.isfinite(yData)
    yData = yData[goodis]
    xData = xData[goodis]
    result = ivAnalyzer.lineFit(xData, yData, 1/guess['n'], 1)
    guess['n'] = 1/result[0]
  
    # take a stab at a guess for I0
    guess['I0'] = float(np.real_if_close(fI0(n=guess['n'],V=V_ip_n,Iph=guess['Iph'],I=I_ip_n,Rs=guess['Rs'],Rsh=guess['Rsh'])))
    
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
    
    moreSyms = sympy.symbols('i1 i2 v1 v2', real=True)
    i1, i2, v1, v2 = moreSyms
    
    
    ceqn = sympy.Eq(lhs,rhs)
    
    
    I0_alone = sympy.solve(ceqn,I0)[0]
    n_alone = sympy.solve(ceqn,n)[0]
    
    n_alone_1 = n_alone.subs([(I,i1),(V,v1)])
    I0_alone_2 = I0_alone.subs([(I,i2),(V,v2)])

    
    #refine guesses for I0 and Rs by forcing the curve through several data points and numerically solving the resulting system of eqns
    #eqnSys1 = zero.subs([(Vth,thermalVoltage),(Iph,guess['Iph']),(V,V_vmpp_n),(I,I_vmpp_n),(Rs,guess['Rs']),(Rsh,guess['Rsh'])])
    #eqnSys2 = zero.subs([(Vth,thermalVoltage),(Iph,guess['Iph']),(V,V_ip_n),(I,I_ip_n),(Rs,guess['Rs']),(Rsh,guess['Rsh'])])
    I0_zero = I0_alone_2.subs(n,n_alone_1) - I0
    
    I0_zero_subs = I0_zero.subs([(Vth,thermalVoltage),(Iph,guess['Iph']),(v1,V_vmpp_n),(i1,I_vmpp_n),(v2,V_ip_n),(i2,I_ip_n),(Rs,guess['Rs']),(Rsh,guess['Rsh'])])
    
    guess['I0'] = float(sympy.nsolve(I0_zero_subs,1e-9))
    guess['n'] = float(n_alone.subs([(Vth,thermalVoltage),(V,V_vmpp_n),(I,I_vmpp_n),(Rs,guess['Rs']),(Rsh,guess['Rsh']),(Iph,guess['Iph']),(I0,guess['I0'])]))
    
    #eqnSys2 = n_alone - n
    
    #eqnSys1n = eqnSys1.subs([(Vth,thermalVoltage),(Iph,guess['Iph']),(V,V_vmpp_n),(I,I_vmpp_n),(Rs,guess['Rs']),(Rsh,guess['Rsh'])])
    #eqnSys2n = eqnSys2.subs([(Vth,thermalVoltage),(Iph,guess['Iph']),(V,V_ip_n),(I,I_ip_n),(Rs,guess['Rs']),(Rsh,guess['Rsh'])])
    
    
    #ceqn1 = ceqn.subs([(Vth,thermalVoltage),(Iph,guess['Iph']),(V,V_vmpp_n),(I,I_vmpp_n),(Rs,guess['Rs']),(Rsh,guess['Rsh'])])
    #ceqn2 = ceqn.subs([(Vth,thermalVoltage),(Iph,guess['Iph']),(V,V_ip_n),(I,I_ip_n),(Rs,guess['Rs']),(Rsh,guess['Rsh'])])
    
    #ceqn1 = ceqn1.rhs - ceqn1.lhs
    #ceqn2 = ceqn2.rhs - ceqn2.lhs
    
    #eqnSysn = (eqnSys1n,eqnSys2n)
    #eqnSys = (eqnSys1,eqnSys2)
    
    #n_collect = n_alone.subs(I0,I0_alone) - n
    #I0_collect = I0_alone.subs(n,n_alone) - I0
    
    #sympy.nsolve(eqnSys,(I0,n),(guess['I0'],guess['n']),maxsteps=10000)
    
    
  
    #try:
    #  sln = sympy.nsolve(eqnSys,(I0,n),(guess['I0'],guess['n']),maxsteps=10000)
    #except:
    #  return([[nan,nan,nan,nan,nan], [nan,nan,nan,nan,nan], nan, "hard fail", 10])
  
    #guess['I0'] = 0.0000000000010296906156466855596021832205111
    #guess['n'] = 48.353157516330718582509891126366
    # now reevaluate n at mpp:
    #guess['n'] = float(np.real_if_close(fn(I0=guess['I0'],V=V_vmpp_n,Iph=guess['Iph'],I=I_vmpp_n,Rs=guess['Rs'],Rsh=guess['Rsh'])))
    
    #ivAnalyzer.visualizeGuess(VV,II,guess,fI,RsYInter,V_ip_n,I_ip_n,V_vp_n,I_vp_n,V_vmpp_n,I_vmpp_n)
    return guess
  
  # so you'd like to see how smart/dumb our really smart guess was...
  def visualizeGuess(VV,II,guess,fI,RsYInter,V_ip_n,I_ip_n,V_vp_n,I_vp_n,V_vmpp_n,I_vmpp_n):
    tehRange = max(VV)-min(VV)
    vv = np.linspace(min(VV)-0.1*tehRange,max(VV)+0.1*tehRange,1000)
    aLine = lambda x,T,Y: [-y + -1/x[0]*t + x[1] for t,y in zip(T,Y)]
    print("My guesses are",guess)
    ii=np.real_if_close(np.array([fI(n=guess['n'],I0=guess['I0'],Iph=guess['Iph'],Rsh=guess['Rsh'],Rs=guess['Rs'],V=v) for v in vv]))
    ii2=np.array(aLine([guess['Rs'],RsYInter],vv,np.zeros(len(vv)))) # Rs fit line
    ii3=np.array(aLine([guess['Rsh'],guess['Iph']],vv,np.zeros(len(vv)))) # Rsh fit line
    plt.title('Guess and raw data')
    plt.plot(vv,ii,'k',label='Char. Eqn.') # char eqn
    plt.plot(vv,ii2,'g',label='R_s') # Rs fit line
    plt.plot(vv,ii3,'r',label='R_sh') # Rsh fit line
    plt.plot(V_ip_n,I_ip_n,'Xc',markersize=10,label='ip')
    plt.plot(V_vp_n,I_vp_n,'Dc',markersize=10,label='vp')
    plt.plot(V_vmpp_n,I_vmpp_n,'+c',markersize=10,label='mpp')
    plt.scatter(VV,II,label='Data')
    plt.grid(b=True)
    plt.legend()
    yRange = max(II) - min(II)
    plt.ylim(min(II)-yRange*.1,max(II)+yRange*.1)
    plt.draw()
    plt.show()
    plt.pause(1)
    
  # here we attempt to fit the input data to the characteristic equation
  def doTheFit(VV, II, fI, guess, bounds, method = 'trf',verbose = 0):
    x0 = [guess['I0'],guess['Iph'],guess['Rs'],guess['Rsh'],guess['n']]
    #x0 = [7.974383037191593e-06, 627.619846736794, 0.00012743239329693432, 0.056948423418631065, 2.0]
    residuals = lambda x,T,Y: np.real(np.array([-y + fI(n=x[4],I0=x[0],Iph=x[1],Rsh=x[3],Rs=x[2],V=t) for t,y in zip(T,Y)]))
    #residuals = lambda x,T,Y: np.abs([np.real_if_close(np.float(fI(n=x[4],I0=x[0],Iph=x[1],Rsh=x[3],Rs=x[2],V=t))) - y for t,y in zip(T,Y)])
    
    #residuals = lambda x,T,Y: np.abs([np.float(fI(n=x[4],I0=x[0],Iph=x[1],Rsh=x[3],Rs=x[2],V=t).real) - y for t,y in zip(T,Y)])
    
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
    fitKwargs['x_scale'] = x0
    #fitKwargs['x_scale'] = 'jac'
    #fitKwargs['loss'] = 'soft_l1'
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
    fitKwargs['verbose'] = verbose
    #fitKwargs['verbose'] = 2
    fitKwargs['method'] = method
    #fitKwargs['method'] = 'lm'
    if fitKwargs['method'] == 'lm':
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
    
    #cellModel = Model(fI, nan_policy='raise', independent_vars=['V'])
    def fcn(VVV,n,I0,Iph,Rsh,Rs):
      try:
        send_up = np.array([np.float(np.real(fI(V=x,n=n,I0=I0,Iph=Iph,Rsh=Rsh,Rs=Rs))) for x in VVV], dtype=float)
        if any(np.isinf(send_up)):
          print('woopsie2')
          send_up = np.full_like(VVV, np.nan)
      except:
        print('woopsie')
        send_up = np.full_like(VVV, np.nan)
      return send_up

    #fI2 = lambda V,n,I0,Iph,Rsh,Rs: np.real(np.abs(fI(V=V,n=n,I0=I0,Iph=Iph,Rsh=Rsh,Rs=Rs)))
    cellModel = Model(fcn, nan_policy='omit')
    cellModel.set_param_hint('Rsh',value=guess['Rsh'])
    cellModel.set_param_hint('n',value=guess['n'])
    cellModel.set_param_hint('Rs',value=guess['Rs'])
    cellModel.set_param_hint('Iph',value=guess['Iph'])
    cellModel.set_param_hint('I0',value=guess['I0'])
    
    cellModel.set_param_hint('Rsh',min=1)
    cellModel.set_param_hint('Rs',min=1)
    
    #fitResult = cellModel.fit(II,cellParams, V=VV)
    #print(fitResult.fit_report())
    
    #cellParams['Rsh'].set(min=0)
    #cellParams['n'].set(min=0)
    #cellParams['Rs'].set(min=0)
    #cellParams['Iph'].set(min=0)
    #cellParams['I0'].set(min=0)
    
    #cellParams['Rsh'].set(expr='Rsh>=0')
    #cellParams['Rs'].set(expr='Rs>=0')
    #cellParams['Rsh'].set(max=1e12)
    #cellParams['n'].set(max=1e5)
    #cellParams['Rs'].set(max=1e5)
    #cellParams['Iph'].set(max=2e4)
    #cellParams['I0'].set(max=1)    
    #fitResult = cellModel.fit(II,cellParams, V=VV)
    #print("FIT2")
    #print(fitResult.fit_report())

    cmodel = Model(fcn, nan_policy='omit')
    params = cmodel.make_params(n=guess['n'], I0=guess['I0'], Iph=guess['Iph'], Rs=guess['Rs'], Rsh=guess['Rsh'])
    params['n'].min=0.1
    params['Rsh'].min=0.1
    params['Rs'].min=0.1
    params['Iph'].min=0
    params['I0'].min=0
    
    #fitResult = cellModel.fit(II, V=VV, method='powell',fit_kws={'options':{'xtol':1e-6,'ftol':1e-6}})
    fitResult = cmodel.fit(II, params, VVV=VV, fit_kws={'maxfev':24000})
    #cellModel.set_param_hint('Rsh',min=-np.inf)
    #cellModel.set_param_hint('Rs',min=-np.inf)
    params['Rsh'].min = float('-inf')
    params['Rs'].min = float('-inf')
    fitResult = cmodel.fit(II, params, VVV=VV, method='nelder',fit_kws={'reduce_fcn':'neglogcauchy','options':{'xatol':1e-14,'maxfev':44000,'fatol':1e-14,'disp':True}})
    print(fitResult.fit_report())
    
  
    # do the fit
    #optimizeResult = scipy.optimize.least_squares(*fitArgs,**fitKwargs)
  
    # do the fit with curve_fit
    #tehf =  lambda XX,m_I0,m_Iph,m_Rs,m_Rsh,m_n: np.real_if_close(np.array([fI(I0=m_I0, Iph=m_Iph, Rs=m_Rs, Rsh=m_Rsh, n=m_n, V=t) for t in XX]))
    #popt, pcov = scipy.optimize.curve_fit(tehf,VV,II,p0=x0,method='trf',verbose=2,xtol=np.finfo(float).eps)
    #ivAnalyzer.analyzeGoodness(VV, II, fI, {'I0':popt[0],'Iph':popt[1],'Rs':popt[2],'Rsh':popt[3],'n':popt[4]}, guess, 'tool')
    
  
  
    #optimizeResult.success = False
    #optimize.curve_fit(I_eqn, VV, II, p0=x0, bounds=fitKwargs['bounds'], diff_step=fitKwargs['diff_step'], method="trf", x_scale="jac", jac ='cs', verbose=1, max_nfev=1200000)
    #scipy.optimize.least_squares(residuals, np.array([  1.20347834e-13,   6.28639109e+02,   1.83005279e-04, 6.49757268e-02,   1.00000000e+00]), jac='cs', bounds=[('-inf', '-inf', '-inf', '-inf', 0), ('inf', 'inf', 'inf', 'inf', 'inf')], method='trf', max_nfev=12000, x_scale='jac', verbose=1,diff_step=[1.203478342631369e-14, 1.4901161193847656e-08, 1.4901161193847656e-08, 1.4901161193847656e-08, 1.4901161193847656e-08])
    if fitResult.success:
     
      return {'success':True,'optParams':fitResult.best_values,'message':fitResult.message,'chi-square':fitResult.chisqr}
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
  def analyzeGoodness(VV, II, fI, params, guess, msg):
    # sum of square of differences between data and fit [A^2]
    print("fit:")
    print(params)                
    print("guess:")
    print(guess)
    print("message:")
    print(msg)
  
    vv=np.linspace(VV[0],VV[-1],1000)
    ii=np.real_if_close(np.array([fI(n=guess['n'],I0=guess['I0'],Iph=guess['Iph'],Rsh=guess['Rsh'],Rs=guess['Rs'],V=v) for v in vv]))
    ii2=np.real_if_close(np.array([fI(n=params['n'],I0=params['I0'],Iph=params['Iph'],Rsh=params['Rsh'],Rs=params['Rs'],V=v) for v in vv]))
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
    plt.pause(1)


if __name__ == '__main__':
  # execute only if run as a script
  main()
