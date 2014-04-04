from batch_iv_analysis_UI import Ui_batch_iv_analysis

from collections import OrderedDict

import os, sys, inspect

from PyQt4.QtCore import QString, QThread, pyqtSignal, QTimer, QSettings, Qt
from PyQt4.QtGui import QApplication, QDialog, QMainWindow, QFileDialog, QTableWidgetItem, QCheckBox

import platform
if not platform.system()=='Windows':
    #these next two lines ensure numerical operation under linux is fast (i can't get gmpy to work on windows, so it's much slower):
    import mpmath.libmp
    assert mpmath.libmp.BACKEND == 'gmpy'
import numpy as np
import sympy
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

#allow current solutoin to operate on vectors of voltages (needed for curve fitting)
def vectorizedCurrent(vVector, I0_n, Iph_n, Rsn_n, Rsh_n, n_n):
    return [float(sympy.re(current_n(x, I0_n,Iph_n,Rsn_n,Rsh_n,n_n))) for x in vVector]
    #TODO: this is VERY bad practice to allow global use of current_n in this function 

class col:
    header = ''
    position = 0
    tooltip = ''

class MainWindow(QMainWindow):
    

    def __init__(self):
        QMainWindow.__init__(self)
        
        self.settings = QSettings("greyltc", "batch-iv-analysis")
        
        self.rows = 0 #keep track of how many rows there are in the table

        self.cols = OrderedDict()
        thisKey = 'file'
        self.cols[thisKey] = col()
        self.cols[thisKey].header = 'File'
        self.cols[thisKey].tooltip = 'File name'

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
            self.graphData[row]['Voc'].__format__('0.4f')+ 'V', 
            xy = (self.graphData[row]['Voc'], 0), xytext = (40, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        plt.annotate(
            float(self.graphData[row]['Isc']).__format__('0.4f') + 'A', 
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
        #remember the last path th user opened
        if self.settings.contains('lastFolder'):
            self.settings.value('lastFolder')
            openDir = self.settings.value('lastFolder').toString()
        else:
            openDir = '.'
        
        fileNames = QFileDialog.getOpenFileNamesAndFilter(directory = openDir, caption="Select one or more files to open", filter = '*.txt')       
        
        self.settings.setValue('lastFolder',os.path.dirname(str(fileNames[0][0])))
        for thisFile in fileNames[0]:
            
            thisPath = str(thisFile)
            fileName = os.path.split(thisPath)[-1]
            self.ui.tableWidget.insertRow(self.rows)
            print "computing: "+ fileName
            for ii in range(len(self.cols)):
                self.ui.tableWidget.setItem(self.rows,ii,QTableWidgetItem())
            
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('file')).setText(fileName)
                      

            #do a thing here:
            #grab the area out of the header section of the data file
            header = []
            fp = open(thisFile, mode='r', buffering=1)
            for ii, line in enumerate(fp):
                header.append(line)
                if ii == 21:
                    break
            fp.close()
            area = float(header[14].split(' ')[3])
            
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('file')).setToolTip(''.join(header))
            

            #read in data
            VV, II = np.loadtxt(str(thisFile),skiprows=25,unpack=True)
            II = II * -1 /1000*area #flip current and scale it to absolute amps


            #sort data by ascending voltage 
            newOrder = VV.argsort()
            VV=VV[newOrder]
            II=II[newOrder]

            maxVoltage = VV[-1];
            minVoltage = VV[0];

            #smoothingDegree = 3 #must be <=5 int, default 3
            #smoothingFactor = 0 #zero sends spline through all datapoints, 
            #iFit = interpolate.UnivariateSpline(VV,II,s=smoothingFactor,k=smoothingDegree)

            #use linear interpolation to find approximation for photocurrent
            iFit = interpolate.interp1d(VV,II)
            I_L_guess = float(iFit(0))

            #we'll use a system of two equations and two unknowns to make intelligent guesses for I0 and n
            #for this we'll pick two data points and force the curve through both of them
            #one of the two points is the data point with the highest voltage:
            V_I0_n = VV[-1] #voltage value to use to calculate I0 guess
            I_I0_n = II[-1] #current value to use to calculate I0 guess

            #the other is the "virtual max power point"
            VVcalc = VV+abs(minVoltage)
            IIcalc = II+abs(min(II))
            Pvirtual= np.array(VVcalc*IIcalc)
            maxIndex = Pvirtual.argmax()
            V_n_n = VV[maxIndex]
            I_n_n = II[maxIndex]
            
            V_p2=optimize.brentq(iFit, minVoltage, maxVoltage)
            I_p2=iFit(V_p2)            

            #these guesses will likely break fitting for curves that don't show rectification
            R_s_guess = 10
            R_sh_guess = 1e6
            
            #compute intellegent guesses for n and I0
            n_initial_guess = 1
            n_guess = n_initial_guess            
            #n_guess = sympy.nsolve(nEqn.subs([(Vth,thermalVoltage),(V_n,V_n_n),(I_n,I_n_n),(Rs,R_s_guess),(Rsh,R_sh_guess),(Iph,I_L_guess),(I_I0,I_I0_n),(V_I0,V_I0_n)]),n_initial_guess)
            #I0_guess = I0_in_n.evalf(subs={Vth:thermalVoltage,Rs:R_s_guess,Rsh:R_sh_guess,Iph:I_L_guess,n:n_guess,I_I0:I_I0_n,V_I0:V_I0_n})
            
            #dumb guess:
            I0_guess= eyeNot[0].evalf(subs={Vth:thermalVoltage,Rs:R_s_guess,Rsh:R_sh_guess,Iph:I_L_guess,n:n_guess,I:I_p2,V:V_p2})
            
            guess = [float(I0_guess), I_L_guess, R_s_guess, R_sh_guess, float(n_guess)]
            fitParams, fitCovariance, infodict, errmsg, ier = optimize.curve_fit(vectorizedCurrent, VV, II,p0=guess,full_output = True)
            
            #catch a failed fit attempt:
            if not ier == 1:
                vv=np.linspace(minVoltage,maxVoltage,1000)
                print "fit:"
                print fitParams                
                print "guess:"
                print guess
                print ier
                print errmsg
                ii=vectorizedCurrent(vv,guess[0],guess[1],guess[2],guess[3],guess[4])
                ii2=vectorizedCurrent(vv,fitParams[0],fitParams[1],fitParams[2],fitParams[3],fitParams[4])
                plt.title('The fit did not terminate properly. This was the guess.')
                plt.plot(vv,ii)
                plt.scatter(VV,II)
                plt.draw()
                plt.show()


            I0_fit = fitParams[0]
            Iph_fit = fitParams[1]
            Rs_fit = fitParams[2]
            Rsh_fit = fitParams[3]
            n_fit = fitParams[4]

            def cellModel(voltageIn):
                voltageIn = np.array(voltageIn)
                return vectorizedCurrent(voltageIn, I0_fit, Iph_fit, Rs_fit, Rsh_fit, n_fit)

            def invCellPower(voltageIn):
                voltageIn = np.array(voltageIn)
                return -1*voltageIn*cellModel(voltageIn)

            vMaxGuess = VV[np.array(VV*II).argmax()]
            powerSearchResults = optimize.minimize(invCellPower,vMaxGuess)
            print "Model fit exit code = " + str(ier) + " power search exit code = " + str(powerSearchResults.status)
            
            #catch a failed max power search:
            if not powerSearchStatus == 0:
                print "power search exit code = " + str(powerSearchResults.status)
                print powerSearchResults.message
            
            vMax = powerSearchResults.x[0]
            iMax = cellModel([vMax])[0]
            pMax = vMax*iMax
            
            Voc_nn=optimize.brentq(cellModel, minVoltage, maxVoltage)#more robust than symbolic?
            #Voc_nn = Voc.evalf(subs={I0:I0_fit, Iph:Iph_fit, Rs:Rs_fit, Rsh:Rsh_fit, n:n_fit, Vth:thermalVoltage})
            Isc_nn = Isc.evalf(subs={I0:I0_fit, Iph:Iph_fit, Rs:Rs_fit, Rsh:Rsh_fit, n:n_fit, Vth:thermalVoltage})
            #Voc_nn = Voc_n(I0_fit, Iph_fit, Rs_fit, Rsh_fit, n_fit, thermalVoltage)
            #Isc_nn = Isc_n(I0_fit, Iph_fit, Rs_fit, Rsh_fit, n_fit, thermalVoltage)
            FF = pMax/(Voc_nn*Isc_nn)

            #error estimation:
            alpha = 0.05 # 95% confidence interval = 100*(1-alpha)

            nn = len(VV)    # number of data points
            p = len(fitParams) # number of parameters

            dof = max(0, nn - p) # number of degrees of freedom

            # student-t value for the dof and confidence level
            tval = t.ppf(1.0-alpha/2., dof) 

            bounds = []
            #calculate 95% confidence interval
            for a, p,var in zip(range(nn), fitParams, np.diag(fitCovariance)):
                sigma = var**0.5
                bs = '[{0}  {1}]'.format(p - sigma*tval, p + sigma*tval)
                bounds.append(bs)

            fitX = np.linspace(minVoltage,maxVoltage,1000)
            fitY = cellModel(fitX)
            self.graphData.append({'fitX':fitX,'fitY':fitY,'i':II,'v':VV,'Voc':Voc_nn,'Isc':Isc_nn,'Vmax':vMax,'Imax':iMax})			


            self.ui.tableWidget.item(self.rows,self.cols.keys().index('pce')).setData(Qt.DisplayRole,float(pMax/area*1e3).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('pmax')).setData(Qt.DisplayRole,float(pMax/area*1e3).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('jsc')).setData(Qt.DisplayRole,float(Isc_nn/area*1e3).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('voc')).setData(Qt.DisplayRole,float(Voc_nn*1e3).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('ff')).setData(Qt.DisplayRole,float(FF).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('rs')).setData(Qt.DisplayRole,float(Rs_fit*area).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('rsh')).setData(Qt.DisplayRole,float(Rsh_fit*area).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('jph')).setData(Qt.DisplayRole,float(Iph_fit/area*1e3).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('j0')).setData(Qt.DisplayRole,float(I0_fit/area*1e9).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('n')).setData(Qt.DisplayRole,float(n_fit).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('n')).setToolTip(bounds[4])
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('Vmax')).setData(Qt.DisplayRole,float(vMax*1e3).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('area')).setData(Qt.DisplayRole,float(area).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('pmax2')).setData(Qt.DisplayRole,float(pMax*1e3).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('isc')).setData(Qt.DisplayRole,float(Isc_nn*1e3).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('iph')).setData(Qt.DisplayRole,float(Iph_fit*1e3).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('iph')).setToolTip(bounds[1])
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('i0')).setData(Qt.DisplayRole,float(I0_fit*1e9).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('i0')).setToolTip(bounds[0])
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('rs2')).setData(Qt.DisplayRole,float(Rs_fit).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('rs2')).setToolTip(bounds[2])
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('rsh2')).setData(Qt.DisplayRole,float(Rsh_fit).__format__('.3f'))
            self.ui.tableWidget.item(self.rows,self.cols.keys().index('rsh2')).setToolTip(bounds[3])

            self.rows = self.rows + 1
        self.ui.tableWidget.setVisible(False)
        self.ui.tableWidget.resizeColumnsToContents()
        self.ui.tableWidget.setVisible(True)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    analysis = MainWindow()
    analysis.show()
    sys.exit(app.exec_())