from batch_iv_analysis_UI import Ui_batch_iv_analysis

from collections import OrderedDict

import os, sys, inspect

from PyQt4.QtCore import QString, QThread, pyqtSignal, QTimer, QSettings, Qt
from PyQt4.QtGui import QApplication, QDialog, QMainWindow, QFileDialog, QTableWidgetItem, QCheckBox

import numpy as np
from numpy import exp, log
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize, root
from scipy.special import lambertw
from scipy import interpolate
from scipy import optimize

#all analysis is done assuming the cell is at 29 degC
cellTemp = 29; #degC 
T = 273.15 + cellTemp; #temp in K
K = 1.3806488e-23; #boltzman constant
q = 1.60217657e-19; #electron charge
Vth = K*T/q; #thermal voltage

#this is ghetto, should not have to use global here
bounds = []

#maximum voltage at maximum power point
def Vmax (guess,I0, Iph, Rs, Rsh, n):
	#myArgs = [I0, Iph, Rs, Rsh, n]
	res =  minimize(PowWrap, guess, args=(I0, Iph, Rs, Rsh, n))
	return res.x[0]

#power as a function of voltage and cell parameters
def PowWrap(x, I0, Iph, Rs, Rsh, n):
	return charEqnI(x[0], I0, Iph, Rs, Rsh, n) * x[0] *-1

#open circuit voltage
def Voc (guess,I0, Iph, Rs, Rsh, n):
	res = root(charEqnI, guess, args=(I0, Iph, Rs, Rsh, n))
	return res.x[0]

#short circuit current
def Isc(I0, Iph, Rs, Rsh, n):
	return charEqnI(0, I0, Iph, Rs, Rsh, n)

#current for nonideal solar cell equation as a function of voltage and cell parameters
def charEqnI(x, I0, Iph, Rs, Rsh, n):
	return np.real((Rsh*(I0+Iph)/(Rs+Rsh))-(x/(Rs+Rsh))-lambertw((Rs*I0*Rsh)/(n*Vth*(Rs+Rsh))*exp((Rsh*(Rs*Iph+Rs*I0+x))/(n*Vth*(Rs+Rsh))))*n*Vth/Rs);

#a wrapper around the nonideal solar cell equation to limit the fit variables to upper and lower bounds
def charWrap(x, I0, Iph, Rs, Rsh, n):
	lowerBounds, upperBounds = bounds #this is ghetto, should not have to use global here
	variables = [I0, Iph, Rs, Rsh, n]
	I0, Iph, Rs, Rsh, n = sigmoidAll(variables, lowerBounds, upperBounds)
	return charEqnI(x, I0, Iph, Rs, Rsh, n)

#restricts input, x to the interval [a, b]
def sigmoid(x,a,b):
	return a+(b-a)/(1+exp(-x))

#inverse sigmoid
def invSigmoid(x,a,b):
	return -log((a - b)/(a - x) - 1)

def sigmoidAll (variables,lowerBounds, upperBounds):
	wat = range(len(variables))
	i=0
	for v in variables:
		wat[i] = sigmoid(variables[i], lowerBounds[i], upperBounds[i])
		i = i+1
	return wat

def invSigmoidAll (variables,lowerBounds,upperBounds):
	wat = range(len(variables))
	i=0
	for v in variables:
		wat[i] = invSigmoid(variables[i], lowerBounds[i], upperBounds[i])
		i = i+1
	return wat

class col:
	header = ''
	position = 0
	tooltip = ''

class MainWindow(QMainWindow):
	def __init__(self):
		QMainWindow.__init__(self)
		self.rows = 0 #keep track of how many rows there are in the table
		
		self.cols = OrderedDict()
		thisKey = 'file'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'File'
		self.cols[thisKey].tooltip = 'File name'
		
		#thisKey = 'pceSpline'
		#self.cols[thisKey] = col()
		#self.cols[thisKey].header = 'PCE (spline)\n[%]'
		#self.cols[thisKey].tooltip = 'Power conversion efficiency found from spline interpolation'		
		
		thisKey = 'pce'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'PCE\n[%]'
		self.cols[thisKey].tooltip = 'Power conversion efficiency'
		
		thisKey = 'pmax'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'P_max\n[mW/cm^2]'
		self.cols[thisKey].tooltip = 'Maximum power density'
		
		thisKey = 'jsc'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'J_sc\n[mA/cm^2]'
		self.cols[thisKey].tooltip = 'Short-circuit current density'
		
		thisKey = 'voc'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'V_oc\n[mV]'
		self.cols[thisKey].tooltip = 'Open-circuit voltage'
		
		thisKey = 'ff'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'FF'
		self.cols[thisKey].tooltip = 'Fill factor'
		
		thisKey = 'rs'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'R_s\n[ohm*cm^2]'
		self.cols[thisKey].tooltip = 'Specific series resistance'
		
		thisKey = 'rsh'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'R_sh\n[ohm*cm^2]'
		self.cols[thisKey].tooltip = 'Specific shunt resistance'
		
		thisKey = 'jph'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'J_ph\n[mA/cm^2]'
		self.cols[thisKey].tooltip = 'Photogenerated current density'
		
		thisKey = 'j0'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'J_0\n[nA/cm^2]'
		self.cols[thisKey].tooltip = 'Reverse saturation current density'
		
		thisKey = 'n'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'n'
		self.cols[thisKey].tooltip = 'Diode ideality factor'
		
		thisKey = 'Vmax'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'V_max\n[mV]'
		self.cols[thisKey].tooltip = 'Voltage at maximum power point'
		
		thisKey = 'area'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'Area\n[cm^2]'
		self.cols[thisKey].tooltip = 'Device area'
		
		thisKey = 'pmax2'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'P_max\n[mW]'
		self.cols[thisKey].tooltip = 'Maximum power'
		
		thisKey = 'isc'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'I_sc\n[mA]'
		self.cols[thisKey].tooltip = 'Short-circuit current'
		
		thisKey = 'iph'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'I_ph\n[mA]'
		self.cols[thisKey].tooltip = 'Photogenerated current'
		
		thisKey = 'i0'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'I_0\n[nA]'
		self.cols[thisKey].tooltip = 'Reverse saturation current'
		
		thisKey = 'rs2'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'R_s\n[ohm]'
		self.cols[thisKey].tooltip = 'Series resistance'
		
		thisKey = 'rsh2'
		self.cols[thisKey] = col()
		self.cols[thisKey].header = 'R_sh\n[ohm]'
		self.cols[thisKey].tooltip = 'Shunt resistance'		

		self.graphData = []
		
		#how long status messages show for
		self.messageDuration = 1000#ms

		# Set up the user interface from Designer.
		self.ui = Ui_batch_iv_analysis()
		self.ui.setupUi(self)
		
		#insert cols
		for item in self.cols:
			blankItem = QTableWidgetItem()
			thisCol = self.cols.keys().index(item)
			self.ui.tableWidget.insertColumn(thisCol)
			blankItem.setToolTip(self.cols[item].tooltip)
			blankItem.setText(self.cols[item].header)
			self.ui.tableWidget.setHorizontalHeaderItem(thisCol,blankItem)
		
		#connect signals generated by gui elements to proper functions 
		self.ui.actionOpen.triggered.connect(self.openCall)
		self.ui.tableWidget.cellDoubleClicked.connect(self.rowGraph)
		
		self.ui.actionClear_Table.triggered.connect(self.clearTableCall)
		
	def rowGraph(self,row):
		v = self.graphData[row]['v']
		i = self.graphData[row]['i']
		plt.scatter(v, i,c='r',marker='o')
		plt.scatter(self.graphData[row]['Vmax'], self.graphData[row]['Imax'], c='g',marker='x',s=100)
		plt.scatter(self.graphData[row]['Voc'], 0, c='g',marker='x',s=100)
		plt.scatter(0, self.graphData[row]['Isc'], c='g',marker='x',s=100)
		fitX = self.graphData[row]['fitX']
		fitY = self.graphData[row]['fitY']
		plt.plot(fitX, fitY)
		plt.autoscale(axis='x', tight=True)
		plt.grid(b=True)
		
		plt.annotate(
			self.graphData[row]['Voc'].__format__('0.4f')+ 'A', 
			xy = (self.graphData[row]['Voc'], 0), xytext = (40, 20),
			textcoords = 'offset points', ha = 'right', va = 'bottom',
			bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
			arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
		
		plt.annotate(
			float(self.graphData[row]['Isc']).__format__('0.4f') + 'V', 
			xy = (0,self.graphData[row]['Isc']), xytext = (40, 20),
			textcoords = 'offset points', ha = 'right', va = 'bottom',
			bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
			arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
		
		plt.annotate(
		        '(' + float(self.graphData[row]['Vmax']).__format__('0.4f') + ',' + float(self.graphData[row]['Imax']).__format__('0.4f') + ')', 
		        xy = (self.graphData[row]['Vmax'],self.graphData[row]['Imax']), xytext = (40, 20),
		        textcoords = 'offset points', ha = 'right', va = 'bottom',
		        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
		        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))		
		
		plt.ylabel('Current [A]')
		plt.xlabel('Voltage [V]')
		
		plt.draw()
		plt.show()
		             
	def clearTableCall(self):
		self.ui.tableWidget.clearContents()
		self.rows = 0

	def openCall(self):
		global bounds #this is ghetto, should not have to use global here
		fileNames = QFileDialog.getOpenFileNamesAndFilter(caption="Select one or more files to open", filter = '*.txt')
		for thisFile in fileNames[0]:				
			thisPath = str(thisFile)
			fileName = os.path.split(thisPath)[-1]
			self.ui.tableWidget.insertRow(self.rows)
			print "computing: "+ fileName
			for ii in range(len(self.cols)):
				self.ui.tableWidget.setItem(self.rows,ii,QTableWidgetItem())			
			
			#get header
			header = []
			fp = open(thisPath, mode='r', buffering=1)
			for ii, line in enumerate(fp):
				header.append(line)
				if ii == 21:
					break
			fp.close()
			
			area = float(header[14].split(' ')[3])	
			
			v, i = np.loadtxt(thisPath,skiprows=25,unpack=True)
			i = i * -1 /1000*area 
			
			maxVoltage = np.max(v);
			minVoltage = np.min(v);
			vGuess = (maxVoltage+minVoltage)/2; #voltage guess for max power and Voc
			
			#initial guess for solar cell parameters
			guess = [7e-9,2e-3,16,1e8,1e0]
			lowerBound = [0,-1e-1,0,   0,   -2]#lower bounds for variables
			upperBound = [1e-3, 2e-1, 400, 1e9, 3]#lower bounds for variables
			bounds = (lowerBound,upperBound) #this is ghetto, should not have to use global here
			
			
			#print invSigmoidAll(guess, lowerBound, upperBound)
			guess = invSigmoidAll(guess, lowerBound, upperBound)
			#do the fit here
			#guess = [0,0,0,0,0]
			try:
				errmsg = "curve_fit hard crash"
				fitParams, fitCovariance, infodict, errmsg, ier = curve_fit(charWrap, v, i,p0=guess,full_output = True)
				fitParams = sigmoidAll(fitParams,lowerBound,upperBound) #unwrap here
				I0 = fitParams[0]
				Iph = fitParams[1]
				Rs = fitParams[2]
				Rsh = fitParams[3]
				n = fitParams[4]				
			except:
				print "Fit failed: " + errmsg
				I0 = np.nan
				Iph = np.nan
				Rs = np.nan
				Rsh = np.nan
				n = np.nan				

			#Iscf = Isc(I0, Iph, Rs, Rsh, n)
			#Vmaxf = Vmax(vGuess,I0, Iph, Rs, Rsh, n)
			#Vocf = Voc(vGuess,I0, Iph, Rs, Rsh, n)
			#Imax = charEqnI(Vmaxf,I0, Iph, Rs, Rsh, n)
			#Pmax = Vmaxf*Imax;
			#FF = Pmax/(Iscf*Vocf)
			v=v[::-1]#need to re-order v for spline calculations
			i=i[::-1]#need to re-order v for spline calculations
			xr = np.linspace(min(v),max(v),100000)

			#interpolation types:
			smoothingFactor = 0 #zero sends spline through all datapoints, 
			currentInterp = interpolate.UnivariateSpline(v,i,s=smoothingFactor)
			#currentInterp = interpolate.interp1d(v,i)
			
			def invPowerFunction (pv):
				return -1*currentInterp(pv)*pv
			
			powerInterp_min_res = minimize(invPowerFunction,(max(v)-min(v))/2)
			Vmax2 = powerInterp_min_res.x[0]
			print Vmax2
			
			pMaxIndex=np.argmax(currentInterp(xr)*xr)
			Vmax = xr[pMaxIndex]
			Pmax = currentInterp(Vmax)*Vmax
			#fit in power space test
			#powerInterp = interpolate.InterpolatedUnivariateSpline(v,i*v)
			#print max(powerInterp(xr))				
			#print Pmax
			
			Isc=currentInterp(0)
			Voc=optimize.brentq(currentInterp, v[0], v[-1])
			FF = Pmax/(Isc*Voc)
			
			fitX = xr
			fitY = currentInterp(xr)
			self.graphData.append({'fitX':fitX,'fitY':fitY,'i':i,'v':v,'Voc':Voc,'Isc':Isc,'Vmax':Vmax,'Imax':currentInterp(Vmax)})			
			
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('file')).setText(fileName)
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('file')).setToolTip(''.join(header))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('pce')).setData(Qt.DisplayRole,float(Pmax/area*1e3).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('pmax')).setData(Qt.DisplayRole,float(Pmax/area*1e3).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('jsc')).setData(Qt.DisplayRole,float(Isc/area*1e3).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('voc')).setData(Qt.DisplayRole,float(Voc*1e3).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('ff')).setData(Qt.DisplayRole,float(FF).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('rs')).setData(Qt.DisplayRole,float(Rs*area).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('rsh')).setData(Qt.DisplayRole,float(Rsh*area).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('jph')).setData(Qt.DisplayRole,float(Iph/area*1e3).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('j0')).setData(Qt.DisplayRole,float(I0/area*1e9).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('n')).setData(Qt.DisplayRole,float(n).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('Vmax')).setData(Qt.DisplayRole,float(Vmax*1e3).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('area')).setData(Qt.DisplayRole,float(area).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('pmax2')).setData(Qt.DisplayRole,float(Pmax*1e3).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('isc')).setData(Qt.DisplayRole,float(Isc*1e3).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('iph')).setData(Qt.DisplayRole,float(Iph*1e3).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('i0')).setData(Qt.DisplayRole,float(I0*1e9).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('rs2')).setData(Qt.DisplayRole,float(Rs).__format__('.3f'))
			self.ui.tableWidget.item(self.rows,self.cols.keys().index('rsh2')).setData(Qt.DisplayRole,float(Rsh).__format__('.3f'))
			
			self.rows = self.rows + 1
		self.ui.tableWidget.setVisible(False)
		self.ui.tableWidget.resizeColumnsToContents()
		self.ui.tableWidget.setVisible(True)
	

if __name__ == "__main__":
	app = QApplication(sys.argv)
	analysis = MainWindow()
	analysis.show()
	sys.exit(app.exec_())
