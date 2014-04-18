from batch_iv_analysis_UI import Ui_batch_iv_analysis

#TODO: intigrate QFileSystemWatcher
#TODO: enable table sorting
#TODO: make area editable
#TODO: handle area from custom input file

from interpolate import SmoothSpline
#cite:
#References
#----------
#.. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
    #Data by Simplified Least Squares Procedures. Analytical
    #Chemistry, 1964, 36 (8), pp 1627-1639.
#.. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
    #W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
    #Cambridge University Press ISBN-13: 9780521880688

from collections import OrderedDict

import os, sys, inspect, csv

from PyQt4.QtCore import QString, QThread, pyqtSignal, QTimer, QSettings, Qt, QSignalMapper, QTemporaryFile, QFileSystemWatcher, QDir, QStringList
from PyQt4.QtGui import QApplication, QDialog, QMainWindow, QFileDialog, QTableWidgetItem, QCheckBox, QPushButton, QWidget

import platform
if not platform.system()=='Windows':
    #these next two lines ensure numerical operation under linux is fast (i can't get gmpy to work on windows, so it's much slower):
    import mpmath.libmp
    assert mpmath.libmp.BACKEND == 'gmpy'
import numpy as np
import sympy
from numpy import nan

from scipy.special import lambertw

from scipy import odr
from scipy import interpolate
from scipy import optimize
from scipy.stats.distributions import  t #needed for confidence interval calculation
import matplotlib.pyplot as plt
#from uncertainties import ufloat

I0, Iph, Rs, Rsh, n, I, V, Vth, V_I0, I_I0, V_n, I_n = sympy.symbols('I0 Iph Rs Rsh n I V Vth V_I0 I_I0 V_n I_n')

cellTemp = 29; #degC all analysis is done assuming the cell is at 29 degC
T = 273.15 + cellTemp; #cell temp in K
K = 1.3806488e-23; #boltzman constant
q = 1.60217657e-19; #electron charge
thermalVoltage = K*T/q; #thermal voltage ~26mv		

#this stuff is from http://dx.doi.org/10.1016/j.solmat.2003.11.018
#symbolic representation for solar cell equation:
lhs = I
rhs = Iph-((V+I*Rs)/Rsh)-I0*(sympy.exp((V+I*Rs)/(n*Vth))-1)
charEqn = sympy.Eq(lhs,rhs)

#isolate current term in solar cell equation
current = sympy.solve(charEqn,I)

#isolate voltage term in solar cell equation
voltage = sympy.solve(charEqn,V)

#isolate I0  in solar cell equation
eyeNot = sympy.solve(charEqn,I0)
#isolate n  in solar cell equation
nidf = sympy.solve(charEqn,n)
RshEqn = sympy.solve(charEqn,Rsh)

#solve for I0 in terms of n:
I0_in_n = eyeNot[0].subs(V,V_I0).subs(I,I_I0)

#use I0 to find a guess value for n:
nEqn = sympy.Eq(n,nidf[0].subs(I0,I0_in_n).subs(V,V_n).subs(I,I_n))
#forNguess = sympy.solve(nEqn,n)


Isc = current[0].subs(V,0)
#Isc_n = sympy.lambdify((I0,Iph,Rs,Rsh,n,Vth),Isc)
Voc = voltage[0].subs(I,0)
#Voc_n = sympy.lambdify((I0,Iph,Rs,Rsh,n,Vth),Voc)

#numeric substitution for thermalVoltage
current = current[0].subs(Vth,thermalVoltage)

#get symbolic solution for current ready for fast number crunching
current_n = sympy.lambdify((V,I0,Iph,Rs,Rsh,n),current)

def odrThing(B,x):
    I0, Iph, Rs, Rsh, n = B
    return np.real((Rs*(I0*Rsh + Iph*Rsh - x) - thermalVoltage*n*(Rs + Rsh)*lambertw(I0*Rs*Rsh*np.exp((Rs*(I0*Rsh + Iph*Rsh - x) + x*(Rs + Rsh))/(thermalVoltage*n*(Rs + Rsh)))/(thermalVoltage*n*(Rs + Rsh))))/(Rs*(Rs + Rsh)))    

def optimizeThis (x,I0, Iph, Rs, Rsh, n):
    return np.real((Rs*(I0*Rsh + Iph*Rsh - x) - thermalVoltage*n*(Rs + Rsh)*lambertw(I0*Rs*Rsh*np.exp((Rs*(I0*Rsh + Iph*Rsh - x) + x*(Rs + Rsh))/(thermalVoltage*n*(Rs + Rsh)))/(thermalVoltage*n*(Rs + Rsh))))/(Rs*(Rs + Rsh)))

#allow current solutoin to operate on vectors of voltages (needed for curve fitting)
def vectorizedCurrent(vVector, I0_n, Iph_n, Rsn_n, Rsh_n, n_n):
    if hasattr(vVector, '__iter__'):
        return [float(sympy.re(current_n(x, I0_n,Iph_n,Rsn_n,Rsh_n,n_n))) for x in vVector]
    else:
        return float(sympy.re(current_n(vVector, I0_n,Iph_n,Rsn_n,Rsh_n,n_n)))
    #TODO: this is VERY bad practice to allow global use of current_n in this function 

#residual function for leastsq fit
def residual(p, vVector, y, weights):
    I0_n, Iph_n, Rsn_n, Rsh_n, n_n = p
    weights = np.asarray(weights)
    return (optimizeThis(vVector, I0_n,Iph_n,Rsn_n,Rsh_n,n_n) -y)*weights
    #return ([float(sympy.re(current_n(x, I0_n,Iph_n,Rsn_n,Rsh_n,n_n))) for x in vVector] - y)*weights
    #TODO: this is VERY bad practice to allow global use of current_n in this function

#residual function for leastsq fit with negative penalty
def residual2(p, vVector, y, weights):
    negativePenalty = 2
    #negativePenalty = 1e-10
    if np.any(p<0):
        return negativePenalty * np.linalg.norm(p[p<0]) * residual(p, vVector, y, weights)
    else:
        return residual(p, vVector, y, weights)

class col:
    header = ''
    position = 0
    tooltip = ''

class MainWindow(QMainWindow):
    workingDirectory = ''
    fileNames = []
    supportedExtensions = QStringList(('*.txt','*.csv'))
    def __init__(self):
        QMainWindow.__init__(self)

        self.settings = QSettings("greyltc", "batch-iv-analysis")

        self.rows = 0 #keep track of how many rows there are in the table

        self.cols = OrderedDict()

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

        thisKey = 'pce'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'PCE\n[%]'
        self.cols[thisKey].tooltip = 'Power conversion efficiency as found from spline fit\nHover for value from characteristic equation fit'

        thisKey = 'pmax'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'P_max\n[mW/cm^2]'
        self.cols[thisKey].tooltip = 'Maximum power density as found from spline fit\nHover for value from characteristic equation fit'

        thisKey = 'jsc'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'J_sc\n[mA/cm^2]'
        self.cols[thisKey].tooltip = 'Short-circuit current density as found from spline fit\nHover for value from characteristic equation fit'

        thisKey = 'voc'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'V_oc\n[mV]'
        self.cols[thisKey].tooltip = 'Open-circuit voltage as found from spline fit\nHover for value from characteristic equation fit'

        thisKey = 'ff'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'FF'
        self.cols[thisKey].tooltip = 'Fill factor as found from spline fit\nHover for value from characteristic equation fit'

        thisKey = 'rs'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'R_s\n[ohm*cm^2]'
        self.cols[thisKey].tooltip = 'Specific series resistance as found from characteristic equation fit\nHover for 95% confidence interval'

        thisKey = 'rsh'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'R_sh\n[ohm*cm^2]'
        self.cols[thisKey].tooltip = 'Specific shunt resistance as found from characteristic equation fit\nHover for 95% confidence interval'

        thisKey = 'jph'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'J_ph\n[mA/cm^2]'
        self.cols[thisKey].tooltip = 'Photogenerated current density as found from characteristic equation fit\nHover for 95% confidence interval'

        thisKey = 'j0'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'J_0\n[nA/cm^2]'
        self.cols[thisKey].tooltip = 'Reverse saturation current density as found from characteristic equation fit\nHover for 95% confidence interval'

        thisKey = 'n'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'n'
        self.cols[thisKey].tooltip = 'Diode ideality factor as found from characteristic equation fit\nHover for 95% confidence interval'

        thisKey = 'Vmax'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'V_max\n[mV]'
        self.cols[thisKey].tooltip = 'Voltage at maximum power point as found from spline fit\nHover for value from characteristic equation fit'

        thisKey = 'area'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'Area\n[cm^2]'
        self.cols[thisKey].tooltip = 'Device area'

        thisKey = 'pmax2'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'P_max\n[mW]'
        self.cols[thisKey].tooltip = 'Maximum power as found from spline fit\nHover for value from characteristic equation fit'

        thisKey = 'isc'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'I_sc\n[mA]'
        self.cols[thisKey].tooltip = 'Short-circuit current as found from characteristic equation fit\nHover for 95% confidence interval'

        thisKey = 'iph'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'I_ph\n[mA]'
        self.cols[thisKey].tooltip = 'Photogenerated current as found from characteristic equation fit\nHover for 95% confidence interval'

        thisKey = 'i0'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'I_0\n[nA]'
        self.cols[thisKey].tooltip = 'Reverse saturation current as found from characteristic equation fit\nHover for 95% confidence interval'

        thisKey = 'rs2'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'R_s\n[ohm]'
        self.cols[thisKey].tooltip = 'Series resistance as found from characteristic equation fit\nHover for 95% confidence interval'

        thisKey = 'rsh2'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'R_sh\n[ohm]'
        self.cols[thisKey].tooltip = 'Shunt resistance as found from characteristic equation fit\nHover for 95% confidence interval'		

        self.graphData = []

        #how long status messages show for
        self.messageDuration = 2500#ms

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
        
        #file system watcher
        self.watcher = QFileSystemWatcher(self)
        self.watcher.directoryChanged.connect(self.handleWatchUpdate)

        #connect signals generated by gui elements to proper functions 
        self.ui.actionOpen.triggered.connect(self.openCall)
        self.ui.actionEnable_Watching.triggered.connect(self.watchCall)
        self.ui.actionSave.triggered.connect(self.handleSave)
        self.ui.actionWatch_2.triggered.connect(self.handleWatchAction)
        

        self.ui.actionClear_Table.triggered.connect(self.clearTableCall)

    def exportInterp(self,row):
        fitX = self.graphData[row]['fitX']
        modelY = self.graphData[row]['modelY']
        splineY = self.graphData[row]['splineY']
        a = np.asarray([fitX, modelY, splineY])
        a = np.transpose(a)
        destinationFolder = os.path.join(self.workingDirectory,'exports')
        QDestinationFolder = QDir(destinationFolder)
        if not QDestinationFolder.exists():
            QDir().mkdir(destinationFolder)
        saveFile = os.path.join(destinationFolder,str(self.ui.tableWidget.item(row,self.cols.keys().index('file')).text())+'.csv')
        header = 'Voltage [V],CharEqn Current [mA/cm^2],Spline Current [mA/cm^2]'
        try:
            np.savetxt(saveFile, a, delimiter=",",header=header)
            self.ui.statusbar.showMessage("Exported " + saveFile,self.messageDuration)
        except:
            self.ui.statusbar.showMessage("Could not export " + saveFile,self.messageDuration)
        


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
        filename = str(self.ui.tableWidget.item(row,self.cols.keys().index('file')).text())
        plt.title(filename)
        v = self.graphData[row]['v']
        i = self.graphData[row]['i']
        plt.plot(v, i, c='b', marker='o', ls="None",label='I-V Data')
        plt.scatter(self.graphData[row]['Vmax'], self.graphData[row]['Imax'], c='g',marker='x',s=100)
        plt.scatter(self.graphData[row]['Voc'], 0, c='g',marker='x',s=100)
        plt.scatter(0, self.graphData[row]['Isc'], c='g',marker='x',s=100)
        fitX = self.graphData[row]['fitX']
        modelY = self.graphData[row]['modelY']
        splineY = self.graphData[row]['splineY']
        if not np.isnan(modelY[0]):
            plt.plot(fitX, modelY,c='k', label='CharEqn Best Fit')
        plt.plot(fitX, splineY,c='g', label='Spline Fit')
        plt.autoscale(axis='x', tight=True)
        plt.grid(b=True)
        ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=3)

        plt.annotate(
            self.graphData[row]['Voc'].__format__('0.4f')+ ' V', 
            xy = (self.graphData[row]['Voc'], 0), xytext = (40, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        plt.annotate(
            float(self.graphData[row]['Isc']).__format__('0.4f') + ' mA/cm^2', 
            xy = (0,self.graphData[row]['Isc']), xytext = (40, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        plt.annotate(
            float(self.graphData[row]['Imax']*self.graphData[row]['Vmax']).__format__('0.4f') + '% @(' + float(self.graphData[row]['Vmax']).__format__('0.4f') + ',' + float(self.graphData[row]['Imax']).__format__('0.4f') + ')', 
            xy = (self.graphData[row]['Vmax'],self.graphData[row]['Imax']), xytext = (80, 40),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))		

        plt.ylabel('Current [mA/cm^2]')
        plt.xlabel('Voltage [V]')

        plt.draw()
        plt.show()

    def handleSave(self):
        if self.settings.contains('lastFolder'):
            self.settings.value('lastFolder')
            saveDir = self.settings.value('lastFolder').toString()
        else:
            saveDir = '.'        
        path = QFileDialog.getSaveFileName(self, caption='Save File', directory=saveDir)
        if not str(path[0]) == '':
            with open(unicode(path), 'wb') as stream:
                writer = csv.writer(stream)
                rowdata = []
                for column in range(self.ui.tableWidget.columnCount()):
                    item = self.ui.tableWidget.horizontalHeaderItem(column)
                    if item is not None:
                        rowdata.append(unicode(item.text()).encode('utf8').replace('\n',' '))
                    else:
                        rowdata.append('')
                writer.writerow(rowdata)                
                for row in range(self.ui.tableWidget.rowCount()):
                    rowdata = []
                    for column in range(self.ui.tableWidget.columnCount()):
                        item = self.ui.tableWidget.item(row, column)
                        if item is not None:
                            rowdata.append(unicode(item.text()).encode('utf8'))
                        else:
                            rowdata.append('')
                    writer.writerow(rowdata)
                stream.close()

    def clearTableCall(self):
        for ii in range(self.rows):
            self.ui.tableWidget.removeRow(0)
        self.ui.tableWidget.clearContents()
        self.rows = 0
        self.graphData = []

    def processFile(self,fullPath):
        fileName, fileExtension = os.path.splitext(fullPath)
        fileName = os.path.basename(fullPath)
        self.fileNames.append(fileName)
        if fileExtension == '.csv':
            delimiter = ','
        else:
            delimiter = None

        self.ui.statusbar.showMessage("processing: "+ fileName,2500)
        
        fp = open(fullPath, mode='r')
        fileBuffer = fp.read()
        fp.close()
        isMcFileFormat = False
        if not fileBuffer[0] == '#':#the first line is not a comment
            if fileBuffer[0:3].find('/') != -1: #the first line is probably a date (assume mcgehee iv file format)
                isMcFileFormat = True
                #comment out the first 25 rows here
                fileBuffer = '#'+fileBuffer
                fileBuffer = fileBuffer.replace('\n', '\n#',24)

        splitBuffer = fileBuffer.splitlines(True)

        #extract header lines
        header = []
        for line in splitBuffer:
            if line.startswith('#'):
                header.append(line)
            else:
                break

        if isMcFileFormat:
            area = float(header[14].split(' ')[3])
            noArea = False
        else:
            area = 1
            noArea = True

        tempFile = QTemporaryFile()
        tempFile.open()
        tempFile.writeData(fileBuffer)
        tempFile.flush()

        #read in data
        try:
            data = np.loadtxt(str(tempFile.fileName()),delimiter=delimiter)
            VV = data[:,0]
            II = data[:,1]
        except:
            self.ui.statusbar.showMessage('Could not read' + fileName +'. Prepend # to all non-data lines and try again',2500)
            VV = nan
            II = nan
            return
        tempFile.close()
        tempFile.remove()

        #assumes current data has the wrong sign and it is in mA
        II = II * -1 /1000*area #flip current and scale it to absolute amps

        #sort data by ascending voltage 
        newOrder = VV.argsort()
        VV=VV[newOrder]
        II=II[newOrder]

        #catch and fix flipped current sign:
        if II[0] < II[-1]:
            II = II * -1       

        indexInQuad1 = np.logical_and(VV>0,II>0)
        if any(indexInQuad1): #enters statement if there is at least one datapoint in quadrant 1
            isDarkCurve = False
        else:
            self.ui.statusbar.showMessage("Dark curve detected",500)
            isDarkCurve = True        

        fitParams, fitCovariance, infodict, errmsg, ier = self.bestEffortFit(VV,II)

        alwaysShowRecap = False
        if  alwaysShowRecap:
            vv=np.linspace(VV[0],VV[-1],1000)
            print "fit:"
            print fitParams                
            print "guess:"
            print guess
            print ier
            print errmsg
            ii=vectorizedCurrent(vv,guess[0],guess[1],guess[2],guess[3],guess[4])
            ii2=vectorizedCurrent(vv,fitParams[0],fitParams[1],fitParams[2],fitParams[3],fitParams[4])
            plt.title('Fit analysis for' + fileName)
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


        I0_fit = fitParams[0]
        Iph_fit = fitParams[1]
        Rs_fit = fitParams[2]
        Rsh_fit = fitParams[3]
        n_fit = fitParams[4]

        #smoothingDegree = 3 #must be <=5 int, default 3
        #smoothingFactor = 1.5e-7 #zero sends spline through all datapoints, 
        #iFitSpline = interpolate.UnivariateSpline(VV,II,s=smoothingFactor,k=smoothingDegree)
        iFitSpline = SmoothSpline(VV, II, p=1-1e-5)

        def cellModel(voltageIn):
            #voltageIn = np.array(voltageIn)
            return vectorizedCurrent(voltageIn, I0_fit, Iph_fit, Rs_fit, Rsh_fit, n_fit)

        def invCellPowerSpline(voltageIn):
            if voltageIn < 0:
                return 0
            else:
                return -1*voltageIn*iFitSpline(voltageIn)

        def invCellPowerModel(voltageIn):
            if voltageIn < 0:
                return 0
            else:
                return -1*voltageIn*cellModel(voltageIn)

        if not isDarkCurve:
            VVq1 = VV[indexInQuad1]
            IIq1 = II[indexInQuad1]
            vMaxGuess = VVq1[np.array(VVq1*IIq1).argmax()]
            powerSearchResults = optimize.minimize(invCellPowerSpline,vMaxGuess)
            #catch a failed max power search:
            if not powerSearchResults.status == 0:
                print "power search exit code = " + str(powerSearchResults.status)
                print powerSearchResults.message
                vMax = nan
                iMax = nan
                pMax = nan
            else:
                vMax = powerSearchResults.x[0]
                iMax = iFitSpline([vMax])[0]
                pMax = vMax*iMax                

            #only do this stuff if the char eqn fit was good
            if ier < 5:
                powerSearchResults_charEqn = optimize.minimize(invCellPowerModel,vMaxGuess)
                #catch a failed max power search:
                if not powerSearchResults_charEqn.status == 0:
                    print "power search exit code = " + str(powerSearchResults_charEqn.status)
                    print powerSearchResults_charEqn.message
                    vMax_charEqn = nan
                else:
                    vMax_charEqn = powerSearchResults_charEqn.x[0]
                #dude
                try:
                    Voc_nn_charEqn=optimize.brentq(cellModel, VV[0], VV[-1])
                except:
                    Voc_nn_charEqn = nan
            else:
                Voc_nn_charEqn = nan
                vMax_charEqn = nan


            try:
                Voc_nn = optimize.brentq(iFitSpline, VV[0], VV[-1])
            except:
                Voc_nn = nan

        else:
            Voc_nn = nan
            vMax = nan
            iMax = nan
            pMax = nan
            Voc_nn_charEqn = nan
            vMax_charEqn = nan
            iMax_charEqn = nan
            pMax_charEqn = nan



        if ier < 5:
            dontFindBounds = False
            iMax_charEqn = cellModel([vMax_charEqn])[0]
            pMax_charEqn = vMax_charEqn*iMax_charEqn
            Isc_nn_charEqn = cellModel(0)
            FF_charEqn = pMax_charEqn/(Voc_nn_charEqn*Isc_nn_charEqn)
        else:
            dontFindBounds = True
            iMax_charEqn = nan
            pMax_charEqn = nan
            Isc_nn_charEqn = nan
            FF_charEqn = nan

        #there is a maddening bug in SmoothingSpline: it can't evaluate 0 alone, so I have to do this:
        try:
            Isc_nn = iFitSpline([0,1e-55])[0]
        except:
            Isc_nn = nan
        #Isc_nn_charEqn = vectorizedCurrent(0,I0_fit,Iph_fit,Rs_fit,Rsh_fit,n_fit)

        #Voc_nn_charEqn = Voc_n(I0_fit, Iph_fit, Rs_fit, Rsh_fit, n_fit, thermalVoltage)
        #Isc_nn_charEqn = Isc_n(I0_fit, Iph_fit, Rs_fit, Rsh_fit, n_fit, thermalVoltage)
        FF = pMax/(Voc_nn*Isc_nn)

        if (ier != 7) and (ier != 6) and (not dontFindBounds):
            #error estimation:
            alpha = 0.05 # 95% confidence interval = 100*(1-alpha)

            nn = len(VV)    # number of data points
            p = len(fitParams) # number of parameters

            dof = max(0, nn - p) # number of degrees of freedom

            # student-t value for the dof and confidence level
            tval = t.ppf(1.0-alpha/2., dof) 

            lowers = []
            uppers = []
            #calculate 95% confidence interval
            for a, p,var in zip(range(nn), fitParams, np.diag(fitCovariance)):
                sigma = var**0.5
                lower = p - sigma*tval
                upper = p + sigma*tval
                lowers.append(lower)
                uppers.append(upper)

        else:
            uppers = [nan,nan,nan,nan,nan]
            lowers = [nan,nan,nan,nan,nan]

        plotPoints = 1000
        fitX = np.linspace(VV[0],VV[-1],plotPoints)
        outputScaleFactor = np.array(1000/area)
        if ier < 5:
            modelY = cellModel(fitX)*outputScaleFactor
        else:
            modelY = np.empty(plotPoints)*nan
        splineY = iFitSpline(fitX)*outputScaleFactor
        self.graphData.append({'origRow':self.rows,'fitX':fitX,'modelY':modelY,'splineY':splineY,'i':II*outputScaleFactor,'v':VV,'Voc':Voc_nn,'Isc':Isc_nn*outputScaleFactor,'Vmax':vMax,'Imax':iMax*outputScaleFactor})			

        self.ui.tableWidget.insertRow(self.rows)
        for ii in range(len(self.cols)):
            self.ui.tableWidget.setItem(self.rows,ii,QTableWidgetItem())

        #plot button
        plotBtn = QPushButton(self.ui.tableWidget)
        plotBtn.setText('Plot')
        plotBtn.clicked.connect(self.handleButton)
        self.ui.tableWidget.setCellWidget(self.rows,self.cols.keys().index('plotBtn'), plotBtn)

        #export button
        exportBtn = QPushButton(self.ui.tableWidget)
        exportBtn.setText('Export')
        exportBtn.clicked.connect(self.handleButton)
        self.ui.tableWidget.setCellWidget(self.rows,self.cols.keys().index('exportBtn'), exportBtn)

        self.ui.tableWidget.item(self.rows,self.cols.keys().index('file')).setText(fileName)
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('file')).setToolTip(''.join(header))            
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('pce')).setData(Qt.DisplayRole,float(pMax/area*1e3).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('pce')).setToolTip(str(float(pMax_charEqn/area*1e3).__format__('.3f')))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('pmax')).setData(Qt.DisplayRole,float(pMax/area*1e3).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('pmax')).setToolTip(str(float(pMax_charEqn/area*1e3).__format__('.3f')))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('jsc')).setData(Qt.DisplayRole,float(Isc_nn/area*1e3).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('jsc')).setToolTip(str(float(Isc_nn_charEqn/area*1e3).__format__('.3f')))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('voc')).setData(Qt.DisplayRole,float(Voc_nn*1e3).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('voc')).setToolTip(str(float(Voc_nn_charEqn*1e3).__format__('.3f')))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('ff')).setData(Qt.DisplayRole,float(FF).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('ff')).setToolTip(str(float(FF_charEqn).__format__('.3f')))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('rs')).setData(Qt.DisplayRole,float(Rs_fit*area).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('rs')).setToolTip('[{0}  {1}]'.format(lowers[2]*area, uppers[2]*area))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('rsh')).setData(Qt.DisplayRole,float(Rsh_fit*area).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('rsh')).setToolTip('[{0}  {1}]'.format(lowers[3]*area, uppers[3]*area))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('jph')).setData(Qt.DisplayRole,float(Iph_fit/area*1e3).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('jph')).setToolTip('[{0}  {1}]'.format(lowers[1]/area*1e3, uppers[1]/area*1e3))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('j0')).setData(Qt.DisplayRole,float(I0_fit/area*1e9).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('j0')).setToolTip('[{0}  {1}]'.format(lowers[0]/area*1e9, uppers[0]/area*1e9))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('n')).setData(Qt.DisplayRole,float(n_fit).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('n')).setToolTip('[{0}  {1}]'.format(lowers[4], uppers[4]))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('Vmax')).setData(Qt.DisplayRole,float(vMax*1e3).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('Vmax')).setToolTip(str(float(vMax_charEqn*1e3).__format__('.3f')))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('area')).setData(Qt.DisplayRole,float(area).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('pmax2')).setData(Qt.DisplayRole,float(pMax*1e3).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('pmax2')).setToolTip(str(float(pMax_charEqn*1e3).__format__('.3f')))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('isc')).setData(Qt.DisplayRole,float(Isc_nn*1e3).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('isc')).setToolTip(str(float(Isc_nn_charEqn*1e3).__format__('.3f')))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('iph')).setData(Qt.DisplayRole,float(Iph_fit*1e3).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('iph')).setToolTip('[{0}  {1}]'.format(lowers[1]*1e3, uppers[1]*1e3))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('i0')).setData(Qt.DisplayRole,float(I0_fit*1e9).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('i0')).setToolTip('[{0}  {1}]'.format(lowers[0]*1e9, uppers[0]*1e9))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('rs2')).setData(Qt.DisplayRole,float(Rs_fit).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('rs2')).setToolTip('[{0}  {1}]'.format(lowers[2], uppers[2]))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('rsh2')).setData(Qt.DisplayRole,float(Rsh_fit).__format__('.3f'))
        self.ui.tableWidget.item(self.rows,self.cols.keys().index('rsh2')).setToolTip('[{0}  {1}]'.format(lowers[3], uppers[3]))

        self.rows = self.rows + 1
        self.ui.tableWidget.setVisible(False)
        self.ui.tableWidget.resizeColumnsToContents()
        self.ui.tableWidget.setVisible(True)        

    def bestEffortFit(self,VV,II):

        #splineTestVV=np.linspace(VV[0],VV[-1],1000)
        #splineTestII=iFitSpline(splineTestVV)
        #p1, = plt.plot(splineTestVV,splineTestII)
        #p3, = plt.plot(VV,II,ls='None',marker='o', label='Data')
        #plt.draw()
        #plt.show()            

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
            print "myguess:"
            print myguess
            vv=np.linspace(min(dataX),max(dataX),1000)
            ii=vectorizedCurrent(vv,myguess[0],myguess[1],myguess[2],myguess[3],myguess[4])
            plt.title('Guess and raw data')
            plt.plot(vv,ii)
            plt.scatter(dataX,dataY)
            plt.grid(b=True)
            plt.draw()
            plt.show()

        #phase 1 guesses:
        I_L_initial_guess = I_sc_n
        R_sh_initial_guess = 1e6

        #compute intellegent guesses for Iph, Rsh by forcing the curve through several data points and numerically solving the resulting system of eqns
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
        n_initial_guess = 2
        I0_initial_guess = eyeNot[0].evalf(subs={Vth:thermalVoltage,Rs:R_s_guess,Rsh:R_sh_guess,Iph:I_L_guess,n:n_initial_guess,I:I_ip_n,V:V_ip_n})                         

        if diaplayAllGuesses:
            initial_guess = [I0_initial_guess, I_L_guess, R_s_guess, R_sh_guess, n_initial_guess]
            evaluateGuessPlot(VV, II, initial_guess)                

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
        guess = [I0_guess, I_L_guess, R_s_guess, R_sh_guess, n_initial_guess]
        if diaplayAllGuesses:
            evaluateGuessPlot(VV, II, guess)

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

        #odrMod = odr.Model(odrThing)
        #myData = odr.Data(VV,II)
        #myodr = odr.ODR(myData, odrMod, beta0=guess,maxit=5000,sstol=1e-20,partol=1e-20)#
        #myoutput = myodr.run()
        #myoutput.pprint()
        #see http://docs.scipy.org/doc/external/odrpack_guide.pdf


        try:
            #myoutput = myodr.run()
            #fitParams = myoutput.beta
            #print myoutput.stopreason
            #print myoutput.info
            #ier = 1
            return(optimize.curve_fit(optimizeThis, VV, II,p0=guess,full_output = True,xtol=1e-12,ftol=1e-14))
            #fitParams, fitCovariance, infodict, errmsg, ier = optimize.leastsq(func=residual, args=(VV, II, np.ones(len(II))),x0=guess,full_output=1,xtol=1e-12,ftol=1e-14)#,xtol=1e-12,ftol=1e-14,maxfev=12000
            #fitParams, fitCovariance, infodict, errmsg, ier = optimize.leastsq(func=residual, args=(VV, II, weights),x0=fitParams,full_output=1,ftol=1e-15,xtol=0)#,xtol=1e-12,ftol=1e-14
        except:
            return([[nan,nan,nan,nan,nan], [nan,nan,nan,nan,nan], nan, "hard fail", 10])



    def openCall(self):
        #remember the last path th user opened
        if self.settings.contains('lastFolder'):
            self.settings.value('lastFolder')
            openDir = self.settings.value('lastFolder').toString()
        else:
            openDir = '.'

        fileNames = QFileDialog.getOpenFileNamesAndFilter(directory = openDir, caption="Select one or more files to open", filter = '(*.txt *.csv);;Folders (*)')       
        #fileNames = QFileDialog.getExistingDirectory(directory = openDir, caption="Select one or more files to open")       
        
        if len(fileNames[0])>0:#check if user clicked cancel
            self.workingDirectory = os.path.dirname(str(fileNames[0][0]))
            self.settings.setValue('lastFolder',self.workingDirectory)
            for fullPath in fileNames[0]:
                fullPath = str(fullPath)
                self.processFile(fullPath)
        
            if self.ui.actionEnable_Watching.isChecked():
                self.handleWatchUpdate(self.workingDirectory)
                watchedDirs = self.watcher.directories()
                self.watcher.removePaths(watchedDirs)
                self.watcher.addPath(self.workingDirectory)
    
    def handleWatchAction(self):
        #remember the last path th user opened
        if self.settings.contains('lastFolder'):
            self.settings.value('lastFolder')
            openDir = self.settings.value('lastFolder').toString()
        else:
            openDir = '.'
        
        myDir = QFileDialog.getExistingDirectory(directory = openDir, caption="Select folder to watch")
        
        if len(myDir)>0:#check if user clicked cancel
            self.workingDirectory = str(myDir)
            self.settings.setValue('lastFolder',self.workingDirectory)
            self.ui.actionEnable_Watching.setChecked(True)
            watchedDirs = self.watcher.directories()
            self.watcher.removePaths(watchedDirs)
            self.watcher.addPath(self.workingDirectory)
            self.handleWatchUpdate(self.workingDirectory)
            
    
    def watchCall(self):
        watchedDirs = self.watcher.directories()
        if self.ui.actionEnable_Watching.isChecked():
            if (self.workingDirectory != ''):
                self.watcher.addPath(self.workingDirectory)
                self.handleWatchUpdate(self.workingDirectory)
        else:
            self.watcher.removePaths(watchedDirs)
            
    def handleWatchUpdate(self,path):
        myDir = QDir(path)
        myDir.setNameFilters(self.supportedExtensions)
        allFilesNow = myDir.entryList()
        allFilesNow = list(allFilesNow)
        allFilesNow = [str(item) for item in allFilesNow]
        #print allFilesNow
        #print self.fileNames
        #print set(allFilesNow) & set(self.fileNames)
        differentFiles = list(set(allFilesNow) ^ set(self.fileNames))
        if differentFiles != []:
            for aFile in differentFiles:
                if self.fileNames.__contains__(aFile):
                    #TODO: delete the file from the table
                    self.ui.statusbar.showMessage('Removed' + aFile,2500)
                else:
                    #process the new file
                    self.processFile(os.path.join(self.workingDirectory,aFile))
        
        


if __name__ == "__main__":
    app = QApplication(sys.argv)
    analysis = MainWindow()
    analysis.show()
    sys.exit(app.exec_())
