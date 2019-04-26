#!/usr/bin/env python

from lmfit import Model
import argparse
import os
from scipy import special
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

import mpmath.libmp
assert mpmath.libmp.BACKEND == 'gmpy'

parser = argparse.ArgumentParser(description='Fit solar cell data')
parser.add_argument('--draw-plots', dest='drawPlots', action='store_true', default=False, help="Draw analysis plots or each file processed")
parser.add_argument('--database', default=':memory:', help="Save/append analysis data to this sqlite database file")
parser.add_argument('--freeze-file', dest='freezeObj', type=argparse.FileType('x'), help="Freeze/draw data to a .csv file")
parser.add_argument('-v', '--verbose', action='store_true', default=False, help="Dump everything to the terminal during analysis")
parser.add_argument('input', type=argparse.FileType('r'), nargs='*', help="File(s) to process")
pArgs = parser.parse_args()

cellTemp = 29 #degC all analysis is done assuming the cell is at 29 degC
T = 273.15 + cellTemp #cell temp in K
K = 1.3806488e-23 #boltzman constant
q = 1.60217657e-19 #electron charge
thermalVoltage = K*T/q #thermal voltage ~26mv
Vth = thermalVoltage
Vth = mpmath.convert(Vth)

tol = 1e-15

# setup our functions
#exp = math.exp
#exp = np.exp
exp = mpmath.exp

# returns I from V and params
def cellEqn(V,n,Rs,Rsh,I0,Iph):
  return (Rs*(I0*Rsh + Iph*Rsh - V) - Vth*n*(Rs + Rsh)*special.lambertw(I0*Rs*Rsh*exp((Rs*(I0*Rsh + Iph*Rsh - V)/(Rs + Rsh) + V)/(Vth*n))/(Vth*n*(Rs + Rsh)),tol=tol))/(Rs*(Rs + Rsh))

# returns V from I and params
def cellEqnV(I,n,Rs,Rsh,I0,Iph):
  I = mpmath.matrix(I)
  n = mpmath.convert(n)
  Rs = mpmath.convert(Rs)
  Rsh = mpmath.convert(Rsh)
  I0 = mpmath.convert(I0)
  Iph = mpmath.convert(Iph)
  #if n < 0:
    #n = 500
  #if Rs < 0:
    #Rs = 0.1
  #if Rsh < 0:
    #Rsh = 0.1
  #if I0 < 0:
    #I0 = 100
  #if Iph < 0:
    #Iph = 100
  inExp = Rsh*(-I + I0 + Iph)/(Vth*n)
  tehExp = pool.map(exp,inExp)
  tehExpList = list(tehExp) # expensive
  tehExpMP = mpmath.matrix(tehExpList)
  inW = I0*Rsh*tehExpMP/(Vth*n)
  tehW = pool.map(mpmath.lambertw,inW)
  tehWList = list(tehW) # expensive
  tehWMP = mpmath.matrix(tehWList)
  evaluation = -I*Rs - I*Rsh + I0*Rsh + Iph*Rsh - Vth*n*tehWMP
  #evaluation = -I*Rs - I*Rsh + I0*Rsh + Iph*Rsh - Vth*n*special.lambertw(I0*Rsh*exp(Rsh*(-I + I0 + Iph)/(Vth*n))/(Vth*n),tol=tol)
  evaluation = np.array(evaluation.tolist(), dtype=np.float64)
  retVal = evaluation
  #retVal = np.real_if_close(evaluation)
  print(n,Rs,Rsh,I0,Iph)
  #print(len(retVal))
  return retVal

def SSE(data, fit_dv):
  r = data - fit_dv
  return (r*r).sum()

# loop through each file in the input
for f in pArgs.input:
  df = pd.read_table(f,sep=',',skiprows=2)
  vv = df.voltage.as_matrix()
  ii = df.current.as_matrix()
  ii = ii*-1
  
  cellModel = Model(cellEqn,nan_policy='omit')
  cellModel.set_param_hint('n',value=1)
  cellModel.set_param_hint('Rs',value=6)
  cellModel.set_param_hint('Rsh',value=1e5)
  cellModel.set_param_hint('Iph',value=20e-3)
  cellModel.set_param_hint('I0',value=1e-9)
  
  #cellModel.set_param_hint('n',min=0)
  #cellModel.set_param_hint('Rs',min=0)
  #cellModel.set_param_hint('Rsh',min=0)
  #cellModel.set_param_hint('Iph',min=0)
  #cellModel.set_param_hint('I0',min=0)  
  
  cellModelV = Model(cellEqnV,nan_policy='omit')
  cellModelV.set_param_hint('n',value=1)
  cellModelV.set_param_hint('Rs',value=6)
  cellModelV.set_param_hint('Rsh',value=1e5)
  cellModelV.set_param_hint('Iph',value=20e-3)
  cellModelV.set_param_hint('I0',value=1e-9)
  
  cellModelV.set_param_hint('n',min=1)
  cellModelV.set_param_hint('Rs',min=1)
  cellModelV.set_param_hint('Rsh',min=1)
  #cellModelV.set_param_hint('Iph',min=0)
  cellModelV.set_param_hint('I0',min=0)  
  
  #fitResult = cellModel.fit(ii, V=vv)
  #resultParams = fitResult.params
  #print(fitResult.fit_report())
  #print(fitResult.message)
  
  fitResult = cellModelV.fit(vv, I=ii,method='nelder')
  resultParams = fitResult.params
  print(fitResult.fit_report())
  print(fitResult.message)
  exit(code=0)
  
  #resultParams['Rsh'].value = resultParams['Rsh'].value + 1000
  #resultParams['Iph'].value = resultParams['Iph'].value - 0.09e-3
  #resultParams['n'].value = resultParams['n'].value + 0.2
  
  fig, ax = plt.subplots()
  fitVals = cellModel.eval(params=resultParams,V=vv)
  ax.plot(vv,fitVals,label='Fit')
  ax.plot(vv,ii,'.',label='Data')
  ax.legend()
  fig.tight_layout()
  
  
  print("SSE={:}".format(SSE(fitVals,ii)))
  
  plt.show()
  
  #cellModelInv = Model(cellEqnV)
  #print(f.name)
  