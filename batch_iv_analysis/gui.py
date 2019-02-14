from batch_iv_analysis_UI import Ui_batch_iv_analysis
from ivAnalyzer import ivAnalyzer

# needed for file watching
import time

# for performance tuning
#import cProfile, pstats, io 
#pr = cProfile.Profile()

import math

#TODO: make area editable

from collections import OrderedDict

import os, sys, inspect, csv

from numpy import inf
import numpy as np

from PyQt5.QtCore import QSettings, Qt, QSignalMapper, QFileSystemWatcher, QDir, QFileInfo, QObject, pyqtSignal, QRunnable
from PyQt5.QtWidgets import QApplication, QMainWindow, QDialog, QFileDialog, QTableWidgetItem, QCheckBox, QPushButton, QItemDelegate

import matplotlib.pyplot as plt
plt.switch_backend("Qt5Agg")

class Object(object):
  pass

def runGUI(analyzer,args):
  app = QApplication(sys.argv)
  analysis = MainWindow(analyzer)
  analysis.show()
  ret = app.exec_()
  sys.exit(ret)

class customSignals(QObject):
  newFitResult = pyqtSignal(object)
  #populateRow = pyqtSignal(object)
  #analysisResult = pyqtSignal(dict)
  #sloppy = pyqtSignal(bool)

#mySignals = customSignals()
  
class col:
  header = ''
  position = 0
  tooltip = ''
  
class FloatDelegate(QItemDelegate):
  def __init__(self, sigFigs, parent=None):
    QItemDelegate.__init__(self, parent=parent)
    self.sigFigs = sigFigs
  def paint(self, painter, option, index):
    value = index.model().data(index, Qt.DisplayRole)
    try:
      number = float(value)
      painter.drawText(option.rect, Qt.AlignLeft|Qt.AlignVCenter, MainWindow.to_precision(value,self.sigFigs))
    except :
      QItemDelegate.paint(self, painter, option, index)  

class MainWindow(QMainWindow):
  workingDirectory = ''
  fileNames = []
  supportedExtensions = ['*.csv','*.tsv','*.txt','*.liv1','*.liv2','*.div1','*.div2', '*.h5']
  bounds = {}
  bounds['I0'] = [0, inf] 
  bounds['Iph'] = [0, inf]
  bounds['Rs'] = [0, inf]
  bounds['Rsh'] = [0, inf]
  bounds['n'] = [0, inf]
  symbolCalcsNotDone = True
  upperVLim = float('inf')
  lowerVLim = float('-inf')
  analyzer = None
  uid = 0 # unique identifier associated with each file

  # for table
  #rows = 0 #this variable keepss track of how many rows there are in the results table
  cols = OrderedDict()
  #nextRow = 0
  
  def closeEvent(self, event):
    pass
    #self.pool.shutdown(wait=False)

  def __init__(self, analyzer):
    QMainWindow.__init__(self)
    
    self.settings = QSettings("greyltc", "batch-iv-analysis")
    self.analyzer = analyzer


    #how long status messages show for
    self.messageDuration = 2500#ms

    # Set up the user interface from Designer.
    self.ui = Ui_batch_iv_analysis()
    self.ui.setupUi(self)
    self.ui.tableWidget.setItemDelegate(FloatDelegate(4))
    
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
    
    thisKey = 'substrate'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'Subs'
    self.cols[thisKey].tooltip = 'Substrate position'    
    
    thisKey = 'pixel'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'Pix'
    self.cols[thisKey].tooltip = 'Pixel number'
    
    thisKey = 'ssPCE'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'ssPCE\n[%]'
    self.cols[thisKey].tooltip = 'Final value taken during max power point tracking stage'
    
    thisKey = 'ssVoc'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'ssV_oc\n[mV]'
    self.cols[thisKey].tooltip = 'Final value taken during V_oc dwell stage'
    
    thisKey = 'ssJsc'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'ssJ_sc\n[mV]'
    self.cols[thisKey].tooltip = 'Final value taken during J_sc dwell stage'
    
    thisKey = 'ssff'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'ssFF\n[%]'
    self.cols[thisKey].tooltip = 'Fill factor as found from the "steady state" Mpp, V_oc and I_sc'

    thisKey = 'direction'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'Dir'
    self.cols[thisKey].tooltip = 'Scan direction'
    
    thisKey = 'pce_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'PCE\n[%]'
    self.cols[thisKey].tooltip = 'Power conversion efficiency as found from spline fit'
    
    thisKey = 'pmax_a_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'P_max\n[mW/cm^2]'
    self.cols[thisKey].tooltip = 'Maximum power density as found from spline fit'

    thisKey = 'voc_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'V_oc\n[mV]'
    self.cols[thisKey].tooltip = 'Open-circuit voltage as found from spline fit I=0 crossing'
    
    thisKey = 'jsc_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'J_sc\n[mA/cm^2]'
    self.cols[thisKey].tooltip = 'Short-circuit current density as found from spline spline fit V=0 crossing'
  
    thisKey = 'ff_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'FF\n[%]'
    self.cols[thisKey].tooltip = 'Fill factor as found from spline fit'
  
    thisKey = 'area'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'Area\n[cm^2]'
    self.cols[thisKey].tooltip = 'Device area'
  
    thisKey = 'suns'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'Suns\n'
    self.cols[thisKey].tooltip = 'Illumination intensity'        
  
    thisKey = 'vmax_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'V_max\n[mV]'
    self.cols[thisKey].tooltip = 'Voltage at maximum power point as found from spline fit'
  
    thisKey = 'isc_spline'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'I_sc\n[mA]'
    self.cols[thisKey].tooltip = 'Short-circuit current as found from spline V=0 crossing'    
  
    thisKey = 'SSE'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'SSE\n[mA^2]'
    self.cols[thisKey].tooltip = 'Sum of the square of the errors between the data points and the fit to the char. eqn. (a measure of fit goodness)'
    
    thisKey = 'n'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'n'
    self.cols[thisKey].tooltip = 'Diode ideality factor as found from characteristic equation fit'
  
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
  
    thisKey = 'pce_fit'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'PCE_fit\n[%]'
    self.cols[thisKey].tooltip = 'Power conversion efficiency as found from characteristic equation fit'
  
    # thisKey = 'pmax_fit'
    # self.cols[thisKey] = col()
    # self.cols[thisKey].header = 'P_max_fit\n[mW]'
    # self.cols[thisKey].tooltip = 'Maximum power as found from characteristic equation fit'
  
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
    
    thisKey = 'jsc_fit'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'J_sc_fit\n[mA/cm^2]'
    self.cols[thisKey].tooltip = 'Short-circuit current density as found from characteristic equation fit V=0 crossing'        

    thisKey = 'isc_fit'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'I_sc_fit\n[mA]'
    self.cols[thisKey].tooltip = 'Short-circuit current as found from characteristic equation fit V=0 crossing'

    thisKey = 'iph'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'I_ph\n[mA]'
    self.cols[thisKey].tooltip = 'Photogenerated current as found from characteristic equation fit'
  
    thisKey = 'jph'
    self.cols[thisKey] = col()
    self.cols[thisKey].header = 'J_ph\n[mA/cm^2]'
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
      self.ui.doFastAndSloppyMathCheckBox.setChecked(True)
      self.settings.setValue('fastAndSloppy',True)
    else:
      self.ui.doFastAndSloppyMathCheckBox.setChecked(self.settings.value('fastAndSloppy') == 'true')
    self.ui.doFastAndSloppyMathCheckBox.stateChanged.connect(self.handleMathChange)
    
    # load setting for multiprocessing
    if not self.settings.contains('multiprocessing'):
      self.ui.useMultithreadingModeCheckBox.setChecked(True)
      self.settings.setValue('multiprocessing',True)
    else:
      value = self.settings.value('multiprocessing') == 'true'
      self.ui.useMultithreadingModeCheckBox.setChecked(value)
      self.ui.analysisThreadsSpinBox.setEnabled(value)
    self.ui.useMultithreadingModeCheckBox.stateChanged.connect(self.handleMultiprocessingChange)
    
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

    if self.settings.contains('threads'):
      self.ui.analysisThreadsSpinBox.setValue(int(self.settings.value('threads')))
    else:
      self.settings.setValue('threads',self.ui.analysisThreadsSpinBox.value())
    self.ui.analysisThreadsSpinBox.valueChanged.connect(self.handleNThreadChange)

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
    
    self.mySignals = customSignals()
    self.mySignals.newFitResult.connect(self._processFitResult)
    #self.mySignals.populateRow.connect(self.populateRow)
    #mySignals.sloppy.connect(self.handleMathFinished)
    #mySignals.analysisResult.connect(self.processFitResult)
    
    if self.analyzer.isFastAndSloppy is None:
      self.analyzer.__dict__['isFastAndSloppy'] = self.ui.doFastAndSloppyMathCheckBox.isChecked()
    if self.analyzer.poolWorkers is None:
      self.analyzer.__dict__['poolWorkers'] = self.ui.analysisThreadsSpinBox.value()
    if self.analyzer.multiprocess is None:
      self.analyzer.__dict__['multiprocess'] = self.ui.useMultithreadingModeCheckBox.isChecked()
      
    self.analyzer.setup()
    #self.analyzer = ivAnalyzer(beFastAndSloppy=beFastAndSloppy, multiprocess=multiprocess, poolWorkers=poolWorkers)
     
    # do symbolic calcs now if needed
    #if self.ui.attemptCharEqnFitCheckBox.isChecked():
    #  if self.multiprocess:
    #    submission = self.pool.submit(ivAnalyzer)
    #    #submission = self.pool.submit(self.analyzer.doSymbolicManipulations,fastAndSloppy=self.ui.doFastAndSloppyMathCheckBox.isChecked())
    #    submission.add_done_callback(self.handleMathFinished)
    #    #self.analyzer.doSymbolicManipulations(fastAndSloppy=self.ui.doFastAndSloppyMathCheckBox.isChecked())
    #    #doSymbolicManipulations(fastAndSloppy=self.ui.doFastAndSloppyMathCheckBox.isChecked())
    #  else:
    #    self.handleMathFinished(ivAnalyzer())
       
  #def handleMathFinished(self,submission):
  #def handleMathFinished(self,thing):
    #self.symbolCalcsNotDone = False
    #if self.multiprocess:
    #  self.analyzer = thing.result()
    #else:
    #  self.analyzer = thing
    #print("One-time symbolic manipulations done!")
    #self.analyzer.numericalize(beFastAndSloppy=self.ui.doFastAndSloppyMathCheckBox.isChecked())
    #print("Fast and sloppy mode =", self.analyzer.isFastAndSloppy)
    #print(self.analyzer)
    #self.analyzer.I_eqn = submission.result()['I_eqn']
    #self.analyzer.P_prime = submission.result()['P_prime']
    #self.analyzer.slns = submission.result()['slns']
    #self.analyzer.electricalModelVarsOnly = submission.result()['electricalModelVarsOnly']
    #print(self.analyzer.I_eqn)
        
  def distillAnalysisParams(self):
    analysisParams = {}
    analysisParams['lowerVLim'] = self.lowerVLim
    analysisParams['upperVLim'] = self.upperVLim
    analysisParams['doFit'] = self.ui.attemptCharEqnFitCheckBox.isChecked()
    analysisParams['bounds'] = self.bounds
    analysisParams['uid'] = self.uid # unique identifier
    self.uid = self.uid + 1
    
    if self.ui.fitMethodComboBox.currentIndex() == 0:
      analysisParams['method'] = 'trf'
    elif self.ui.fitMethodComboBox.currentIndex() == 1:
      analysisParams['method'] = 'dogbox'
    elif self.ui.fitMethodComboBox.currentIndex() == 2:
      analysisParams['method'] = 'lm'
    
    analysisParams['verbose'] = self.ui.verbositySpinBox.value()
    
    return analysisParams
  
  def updatePoolStatus(self):
    self.myShowMessage(self.analyzer.getPoolStatusString())
    
  def resetDefaults(self):
    self.ui.attemptCharEqnFitCheckBox.setChecked(True)
    self.ui.doFastAndSloppyMathCheckBox.setChecked(True)
    self.ui.lowerVoltageCutoffLineEdit.setText('-inf')
    self.ui.lowerVoltageCutoffLineEdit.editingFinished.emit()
    self.ui.upperVoltageCutoffLineEdit.setText('inf')
    self.ui.upperVoltageCutoffLineEdit.editingFinished.emit()
    self.ui.fitMethodComboBox.setCurrentIndex(2)
    self.ui.verbositySpinBox.setValue(0)
    self.ui.analysisThreadsSpinBox.setValue(8)
    self.ui.analysisThreadsSpinBox.setEnabled(True)
    self.ui.useMultithreadingModeCheckBox.setChecked(True)

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

  def handleNThreadChange(self):
    spinBox = self.sender()
    value = spinBox.value()
    self.settings.setValue('threads',value)
    self.analyzer.poolWorkers = value

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
    self.analyzer.isFastAndSloppy = checkBox.isChecked()
    #self.analyzer.numericalize(beFastAndSloppy=checkBox.isChecked())
    #print("Fast and sloppy mode =", self.analyzer.isFastAndSloppy)
    
  def handleEqnFitChange(self):
    checkBox = self.sender()
    self.settings.setValue('fitToEqn',checkBox.isChecked())
    
  def handleMultiprocessingChange(self):
    checkBox = self.sender()
    value = checkBox.isChecked()
    self.settings.setValue('multiprocessing',value)
    self.ui.analysisThreadsSpinBox.setEnabled(value)
    self.analyzer.multiprocess = value

  def handleButton(self):
    btn = self.sender()
    #kinda hacky:
    row = self.ui.tableWidget.indexAt(btn.pos()).row()
    col = self.ui.tableWidget.indexAt(btn.pos()).column()
    if col == 0:
      self.rowGraph(row)
    elif col == 1:
      self.exportInterp(row)
    elif col == self.getCol('ssVoc'):
      self.ssVocGraph(row)
    elif col == self.getCol('ssJsc'):
      self.ssJscGraph(row)
    elif col == self.getCol('ssPCE'):
      self.mpptGraph(row)
      
  def ssVocGraph(self, row):
    thisGraphData = self.ui.tableWidget.item(row, self.getCol('plotBtn')).data(Qt.UserRole)
    filename = str(self.ui.tableWidget.item(row, self.getCol('file')).text())
    substrate = str(self.ui.tableWidget.item(row, self.getCol('substrate')).text())
    pixel = str(self.ui.tableWidget.item(row, self.getCol('pixel')).text())
    
    measurements =  thisGraphData['ssVoc']
    v = np.array([e[0] for e in measurements])
    i = np.array([e[1] for e in measurements])
    t = np.array([e[2] for e in measurements])
    s = np.array([int(e[3]) for e in measurements])
    
    x = t - t[0]
    y = abs(v * 1000)
    
    plt.plot(x, y, c='b', marker='o', ls="None",label='Voc')
    
    plt.title("{:}, Pixel {:}{:}".format(filename, substrate, pixel))
    plt.ylabel('Open-circuit voltage [mV]')
    plt.xlabel('Time [s]')
    plt.grid()
    plt.show()
    
  def ssJscGraph(self, row):
    thisGraphData = self.ui.tableWidget.item(row, self.getCol('plotBtn')).data(Qt.UserRole)
    filename = str(self.ui.tableWidget.item(row, self.getCol('file')).text())
    substrate = str(self.ui.tableWidget.item(row, self.getCol('substrate')).text())
    pixel = str(self.ui.tableWidget.item(row, self.getCol('pixel')).text())
    area = self.ui.tableWidget.item(row, self.getCol('area')).data(Qt.UserRole)
    areacm = area * 1e4
    
    measurements =  thisGraphData['ssIsc']
    v = np.array([e[0] for e in measurements])
    i = np.array([e[1] for e in measurements])
    t = np.array([e[2] for e in measurements])
    s = np.array([int(e[3]) for e in measurements])
    
    x = t - t[0]
    y = abs(i * 1000) / areacm
    
    plt.plot(x, y, c='b', marker='o', ls="None",label='Jsc')
    
    plt.title("{:}, Pixel {:}{:}".format(filename, substrate, pixel))
    plt.ylabel('Short-circuit current density [mA/cm^2]')
    plt.xlabel('Time [s]')
    plt.grid()
    plt.show()
    
  def mpptGraph(self, row):
    thisGraphData = self.ui.tableWidget.item(row, self.getCol('plotBtn')).data(Qt.UserRole)
    filename = str(self.ui.tableWidget.item(row, self.getCol('file')).text())
    substrate = str(self.ui.tableWidget.item(row, self.getCol('substrate')).text())
    pixel = str(self.ui.tableWidget.item(row, self.getCol('pixel')).text())
    area = self.ui.tableWidget.item(row, self.getCol('area')).data(Qt.UserRole)
    areacm = area * 1e4
    
    measurements =  thisGraphData['mppt']
    v = np.array([e[0] for e in measurements])
    i = np.array([e[1] for e in measurements])
    t = np.array([e[2] for e in measurements])
    s = np.array([int(e[3]) for e in measurements])
    
    x = t - t[0]
    y = abs(i*v * 1000) / areacm
    
    plt.plot(x, y, c='b', marker='o', ls="None",label='mppt')
    
    plt.title("{:}, Pixel {:}{:}".format(filename, substrate, pixel))
    plt.ylabel('Power Density [mW/cm^2]')
    plt.xlabel('Time [s]')
    plt.grid()
    plt.show()
    

  def rowGraph(self,row):
    thisGraphData = self.ui.tableWidget.item(row, self.getCol('plotBtn')).data(Qt.UserRole)
    filename = str(self.ui.tableWidget.item(row, self.getCol('file')).text())
    substrate = str(self.ui.tableWidget.item(row, self.getCol('substrate')).text())
    pixel = str(self.ui.tableWidget.item(row, self.getCol('pixel')).text())
    direction = str(self.ui.tableWidget.item(row, self.getCol('direction')).text())

    v = thisGraphData["v"]
    i = thisGraphData["i"]

    if direction == 'Fwd.':
      plt.plot(v, i, c='b', marker='o', ls="None",label='I-V Data (Fwd.)')
    else:
      plt.plot(v, i, c='r', marker='o', ls="None",label='I-V Data (Rev.)')
    plt.scatter(thisGraphData["Vmax"], thisGraphData["Imax"], c='g',marker='x',s=100)
    plt.scatter(thisGraphData["Voc"], 0, c='g',marker='x',s=100)
    plt.scatter(0, thisGraphData["Isc"], c='g',marker='x',s=100)
    fitX = thisGraphData["fitX"]
    modelY = thisGraphData["modelY"]
    modelY = np.array(thisGraphData["modelY"]).astype(complex)
    splineY = thisGraphData["splineY"]
    if not np.isnan(modelY[0]):
      plt.plot(fitX, modelY,c='k', label='CharEqn Best Fit')
    plt.plot(fitX, splineY,c='g', label='Spline Fit')
    plt.autoscale(axis='x', tight=True)
    plt.grid()

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
      float(thisGraphData["Imax"]*thisGraphData["Vmax"]).__format__('0.4f') + 'mW/cm^2 @(' + float(thisGraphData["Vmax"]).__format__('0.4f') + ',' + float(thisGraphData["Imax"]).__format__('0.4f') + ')', 
            xy = (thisGraphData["Vmax"],thisGraphData["Imax"]), xytext = (80, 40),
              textcoords = 'offset points', ha = 'right', va = 'bottom',
              bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
              arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))		

    plt.ylabel('Current [mA/cm^2]')
    plt.xlabel('Voltage [V]')

    plt.title("{:}, Pixel {:}{:}".format(filename, substrate, pixel))
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc=3)    
    # ax.grid()
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
      with open(fullPath, 'w',newline='') as stream:
        writer = csv.writer(stream, dialect="excel")
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

      fieldsToInclude= ('pce_spline','pmax_a_spline','voc_spline','isc_spline','ff_spline','vmax_spline','SSE','pce_fit','pmax_a_fit','voc_fit','isc_fit','ff_fit','vmax_fit','rs','rsh','iph','i0','n','area','suns')

      #how many padding zeros should we use for the MATLAB variable names?
      ndigits = str(len(str(self.ui.tableWidget.rowCount()))) 

      for row in range(self.ui.tableWidget.rowCount()):
        rowDict = {}
        rowDict['file'] = self.ui.tableWidget.item(row, list(self.cols.keys()).index('file')).data(Qt.DisplayRole)
        for field in fieldsToInclude:
          rowDict[field] = self.ui.tableWidget.item(row, list(self.cols.keys()).index(field)).data(Qt.UserRole)
        rowDict['i'] = self.ui.tableWidget.item(row, list(self.cols.keys()).index('plotBtn')).data(Qt.UserRole)['i']/rowDict['area']
        rowDict['v'] = self.ui.tableWidget.item(row, list(self.cols.keys()).index('plotBtn')).data(Qt.UserRole)['v']
        tableDict['thing'+format(row, '0'+ndigits)] = rowDict

      # save our dict as a .mat file
      sio.savemat(fullPath, tableDict)
      print('Table data successfully written to', fullPath)

  # takes cell data and modifies it for display
  def sanitizeRow(self,row):      
    ignoreCols = ['plotBtn','exportBtn','file']
    cols = list(self.cols.keys())
    for coli in range(len(cols)):
      thisCol = cols[coli]      
      if thisCol not in ignoreCols:
        thisTableItem = self.ui.tableWidget.item(row,coli)
        if thisTableItem is not None:
          value = thisTableItem.data(Qt.UserRole)
          if value is not None and not np.isnan(value):
            saneValue = float(np.real(value))
            if thisCol == 'SSE':
              displayValue = saneValue*1000**2 # A^2 to mA^2
            elif thisCol in ['ff_spline','ff_fit','pce_spline','ssPCE','ssff']:
              displayValue = saneValue*100 # to percent
            elif thisCol in ['ssVoc','voc_spline','voc_fit','vmax_spline','vmax_fit','isc_spline','isc','iph']:
              displayValue = saneValue*1e3 # to milli-
            elif thisCol in ['jsc_spline','jsc','ssJsc','jph','pmax_a_spline','pmax_a_fit']:
              displayValue = saneValue*1e3*1e-4 # to milli- per cm^2
            elif thisCol in ['area']:
              displayValue = saneValue*1e4 # to cm^2
            elif thisCol in ['i0','j0']:
              displayValue = saneValue*1e9 # to nano-
            else:
              displayValue = saneValue
            displayValue = MainWindow.to_precision(displayValue,4)
            self.ui.tableWidget.item(row,coli).setData(Qt.DisplayRole,float(displayValue))
            self.ui.tableWidget.resizeColumnToContents(coli)
            self.ui.tableWidget.viewport().update()

  # returns table column number given name
  def getCol(self,colName):
    return list(self.cols.keys()).index(colName)
  
  # returns row number associated with a unique identifier
  def getRowByUID(self,uid):
    nRows = self.ui.tableWidget.rowCount()
    fileCol = self.getCol('file')
    row = None
    for i in range(nRows):
      thisCellItem = self.ui.tableWidget.item(i,fileCol)
      if thisCellItem.data(Qt.UserRole) == uid:
        row = i
        break
    return row
  
  def clearTableCall(self):
    for ii in range(self.ui.tableWidget.rowCount()):
      self.ui.tableWidget.removeRow(0)

    
    #self.ui.tableWidget.clear()
    #self.ui.tableWidget.clearContents()
    self.fileNames = []
    
  def newFiles(self, fullPaths):
    self.analyzer.processFiles(fullPaths, self.processFitResult, self.primeRow)

  def primeRow(self, fullPath, fileData):
    """primes a new row in the table"""
    #analysisParams = []
    
    #for i in range(len(fullPaths)):
      # grab settings from gui
      #analysisParams.append(self.distillAnalysisParams())
    params = self.distillAnalysisParams()
      
    #wait here for the file to be completely written to disk and closed before trying to read it
    fi = QFileInfo(fullPath)
    while (not fi.isWritable()):
      time.sleep(0.01)
      fi.refresh()
      
    # insert filename into table immediately
    thisRow = self.ui.tableWidget.rowCount()
    self.ui.tableWidget.setSortingEnabled(False) # fix strange sort behavior
    self.ui.tableWidget.insertRow(thisRow)
    for ii in range(self.ui.tableWidget.columnCount()):
      self.ui.tableWidget.setItem(thisRow,ii,QTableWidgetItem())
    fileName = os.path.basename(fullPath)
      
    self.tableInsert(thisRow,'file', fileName, role=Qt.DisplayRole)
    self.tableInsert(thisRow,'file', params['uid'])
    
    self.tableInsert(thisRow,'substrate', fileData.substrate, role=Qt.DisplayRole)
    self.tableInsert(thisRow,'pixel', fileData.pixel, role=Qt.DisplayRole)
    self.tableInsert(thisRow,'direction', 'Rev.' if fileData.reverseSweep else 'Fwd.', role=Qt.DisplayRole)
    
    graphData = {}
    if hasattr(fileData, 'mppt'):
      graphData['mppt'] = fileData.mppt
    if hasattr(fileData, 'ssVoc'):
      graphData['ssVoc'] = fileData.ssVoc
    if hasattr(fileData, 'ssIsc'):
      graphData['ssIsc'] = fileData.ssIsc
    

    self.tableInsert(thisRow,'plotBtn', graphData)
    
    self.tableInsert(thisRow,'suns', fileData.suns)
    self.tableInsert(thisRow,'area', fileData.area)  # in m^2
    
    #if hasattr(fileData, 'Impp'):
      #self.tableInsert(thisRow,'???', fileData.Impp)
    #if hasattr(fileData, 'Vmpp'):
      #self.tableInsert(thisRow,'???', fileData.Vmpp)
    if hasattr(fileData, 'Voc'):
      self.tableInsert(thisRow,'ssVoc', fileData.Voc)
    if hasattr(fileData, 'ssPmax'):
      self.tableInsert(thisRow,'ssPCE', fileData.ssPmax / fileData.area / ivAnalyzer.stdIrridance / fileData.suns)
    if hasattr(fileData, 'Isc'):
      self.tableInsert(thisRow,'ssJsc', fileData.Isc / fileData.area)
    
    if hasattr(fileData, 'Isc') and hasattr(fileData, 'Voc') and hasattr(fileData, 'ssPmax'):
      self.tableInsert(thisRow,'ssff', abs(fileData.ssPmax/(fileData.Isc*fileData.Voc)))

    self.ui.tableWidget.setSortingEnabled(True) # fix strange sort behavior
    self.fileNames.append(fileName)
      
    return params

  def tableInsert(self,thisRow,colName,value,role=Qt.UserRole):
    thisCol = self.getCol(colName)
    thisItem = self.ui.tableWidget.item(thisRow,thisCol)
    thisItem.setData(role,value)
    self.ui.tableWidget.resizeColumnToContents(thisCol)
    
  def processFitResult(self,result):
    try:# this handles the multiprocessing case
      if result.done():
        exception = result.exception(timeout=0)
        if exception is None:
          result = result.result()
        else:
          print('Error during file processing:', result.exception(timeout=0))
          return
      else:
        print("Somehow the future isn't 'done'")
        return
    except:
      pass
    
    self.mySignals.newFitResult.emit(result)
    #self._processFitResult(result)
    
  def _processFitResult(self,result):
    if self.ui.useMultithreadingModeCheckBox.isChecked():
      self.updatePoolStatus()
    uid = result.params['uid']
    thisRow = self.getRowByUID(uid)
    #print('Got new fit result, UID:',result['params']['uid'])
    #print(result['fitResult'])
    
    #
    
    #thisItem = QTableWidgetItem()
    #thisItem.setData
    
    fitData = Object()
    fitData.pmax_spline = result.pmpp
    fitData.vmax_spline = result.vmpp
    fitData.isc_spline = result.isc
    fitData.voc_spline = result.voc
    fitData.row = thisRow
    fitData.v = result.v
    fitData.i = result.i
    fitData.x = result.x
    fitData.splineCurrent = result.splineCurrent
    
    fitData.SSE = result.sse if hasattr(result,'sse') else np.nan
    fitData.eqnCurrent = result.eqnCurrent if hasattr(result,'eqnCurrent') else np.array([np.nan])
    fitData.n = result.n if hasattr(result,'n') else np.nan
    fitData.rs = result.rs if hasattr(result,'rs') else np.nan
    fitData.rsh = result.rsh if hasattr(result,'rsh') else np.nan
    fitData.i0 = result.i0 if hasattr(result,'i0') else np.nan
    fitData.iph = result.iph if hasattr(result,'iph') else np.nan
    fitData.pmax_fit = result.pmax_fit if hasattr(result,'pmax_fit') else np.nan
    fitData.isc_fit = result.isc_fit if hasattr(result,'isc_fit') else np.nan
    fitData.voc_fit = result.voc_fit if hasattr(result,'voc_fit') else np.nan
    fitData.vmax_fit = result.vmax_fit if hasattr(result,'vmax_fit') else np.nan
    
    #print('got new fit result:',uid)
    self.populateRow(fitData)
    #self.mySignals.populateRow.emit(rowData)
    
    
    #self.tableInsert(thisRow, 'pce_spline', getattr(result, 'pce'))
    #item = self.ui.tableWidget.item(thisRow,self.getCol(thisThing))
    #thisItem = QTableWidgetItem()
    #value = result['insert'][thisThing]
    #role = Qt.UserRole
    #item.setData(role,value)
    
    #insert = lambda colName,value: self.ui.tableWidget.item(thisRow,self.getCol(colName)).setData(Qt.UserRole,value)
    #thisThing = 'pce_spline'
    #insert(thisThing,result['insert'][thisThing])

  def populateRow(self,fitData):
    self.ui.tableWidget.setSortingEnabled(False) # fix strange sort behavior
    
    # add in the export button
    exportBtn = QPushButton(self.ui.tableWidget)
    exportBtn.setText('Export')
    exportBtn.clicked.connect(self.handleButton)
    exportCol = self.getCol('exportBtn')
    self.ui.tableWidget.setCellWidget(fitData.row, exportCol, exportBtn)
    
    # add in the plot button
    plotBtn = QPushButton(self.ui.tableWidget)
    plotBtn.setText('Plot')
    plotBtn.clicked.connect(self.handleButton)
    plotCol = self.getCol('plotBtn')
    self.ui.tableWidget.setCellWidget(fitData.row, plotCol, plotBtn)
    
    if self.ui.tableWidget.item(fitData.row, plotCol).data(Qt.UserRole) == None:
      graphData = {}
    else:
      graphData = self.ui.tableWidget.item(fitData.row, plotCol).data(Qt.UserRole)
      
    # copy fit data over to row data
    rowData = fitData
    
    #  retrieve area and intensity from the table
    area = self.ui.tableWidget.item(fitData.row, self.getCol('area')).data(Qt.UserRole)  # in m^2
    suns = self.ui.tableWidget.item(fitData.row, self.getCol('suns')).data(Qt.UserRole)
    areacm = area * 1e4  #  area in cm^2
    
    # derived row data values:
    rowData.pce_spline = rowData.pmax_spline / area / ivAnalyzer.stdIrridance / suns
    rowData.pmax_a_spline = rowData.pmax_spline / area
    rowData.ff_spline = rowData.pmax_spline / (rowData.isc_spline*rowData.voc_spline)
    rowData.jsc_spline = rowData.isc_spline / area
    rowData.rs_a = rowData.rs*area
    rowData.rsh_a = rowData.rsh/area
    rowData.jph = rowData.iph/area
    rowData.j0 = rowData.i0/area
    rowData.pce_fit = rowData.pmax_fit / area / ivAnalyzer.stdIrridance / suns
    rowData.ff_fit = rowData.pmax_fit/(rowData.isc_fit*rowData.voc_fit)
    rowData.jsc_fit = rowData.isc_fit/area
    rowData.pmax_a_fit = rowData.pmax_fit / area
    
    graphData["v"] = rowData.v
    graphData["i"] = rowData.i/areacm * 1000  # in mA/cm^2
    graphData["vsTime"] = False
    graphData["Vmax"] = rowData.vmax_spline
    graphData["Imax"] = rowData.pmax_spline/rowData.vmax_spline * 1000  # in mA
    graphData["Voc"] = rowData.voc_spline
    graphData["Isc"] = rowData.isc_spline * 1000  # in mA
    graphData["fitX"] = rowData.x
    graphData["modelY"] = rowData.eqnCurrent/areacm * 1000  # in mA/cm^2
    graphData["splineY"] = rowData.splineCurrent/areacm * 1000  # in mA/cm^2
    self.ui.tableWidget.item(rowData.row, plotCol).setData(Qt.UserRole, graphData)
    
    for key,value in rowData.__dict__.items():
      colName = key
      if key not in ['row','i','v','vsTime','x','splineCurrent','eqnCurrent', 'area', 'suns', 'pmax_spline', 'pmax_fit']:
        self.tableInsert(rowData.row, key, value)
        
    # add in the Voc button
    thisGraphData = self.ui.tableWidget.item(rowData.row, self.getCol('plotBtn')).data(Qt.UserRole)
    if 'ssVoc' in thisGraphData:
      vocBtn = QPushButton(self.ui.tableWidget)
      vocCol = self.getCol('ssVoc')
      voc = abs(self.ui.tableWidget.item(fitData.row, vocCol).data(Qt.UserRole))
      vocBtn.setText("{:}".format(MainWindow.to_precision(voc*1000,4)))
      vocBtn.clicked.connect(self.handleButton)
      self.ui.tableWidget.setCellWidget(rowData.row, vocCol, vocBtn)
    
    # add in the Jsc button
    if 'ssIsc' in thisGraphData:
      jscBtn = QPushButton(self.ui.tableWidget)
      jscCol = self.getCol('ssJsc')
      jsc = abs(self.ui.tableWidget.item(fitData.row, jscCol).data(Qt.UserRole))
      jscBtn.setText("{:}".format(MainWindow.to_precision(jsc*1000*1e-4,4)))
      jscBtn.clicked.connect(self.handleButton)
      
      self.ui.tableWidget.setCellWidget(rowData.row, jscCol, jscBtn)
    
    # add in the PCE button
    if 'mppt' in thisGraphData:
      pceBtn = QPushButton(self.ui.tableWidget)
      pceCol = self.getCol('ssPCE')
      pce = self.ui.tableWidget.item(fitData.row, pceCol).data(Qt.UserRole)
      pceBtn.setText("{:}".format(MainWindow.to_precision(pce*100,4)))
      pceBtn.clicked.connect(self.handleButton)
      self.ui.tableWidget.setCellWidget(rowData.row, pceCol, pceBtn)      
    
    self.sanitizeRow(rowData.row)
    self.ui.tableWidget.setSortingEnabled(True)
    
  def openCall(self):
    #remember the last path the user opened
    if self.settings.contains('lastFolder'):
      openDir = self.settings.value('lastFolder')
    else:
      openDir = '.'

    fileNames = QFileDialog.getOpenFileNames(self, directory = openDir, caption="Select one or more files to open", filter = '(*.csv *.tsv *.txt *.liv1 *.liv2 *.div1 *.div2 *.h5);;Folders (*)')

    if len(fileNames[0])>0:#check if user clicked cancel
      self.workingDirectory = os.path.dirname(str(fileNames[0][0]))
      self.settings.setValue('lastFolder',self.workingDirectory)
      
      fullPaths = fileNames[0]
      self.newFiles(fullPaths)
      
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

    newFiles = list(set(allFilesNow) - set(self.fileNames))
    if newFiles != []:
      # prepend full path
      for i in range(len(newFiles)):
        newFiles[i] = os.path.join(self.workingDirectory,newFiles[i])
      # process all the new files
      self.newFiles(newFiles)
      
  def statusChanged(self,args):
    if not args:
      # reset the statusbar background
      self.ui.statusbar.setStyleSheet("QStatusBar{padding-left:8px;background:rgba(0,0,0,0);color:black;font-weight:bold;}")

  def goodMessage(self):
    self.ui.statusbar.setStyleSheet("QStatusBar{padding-left:8px;background:rgba(0,128,0,255);color:black;font-weight:bold;}")

  def badMessage(self):
    self.ui.statusbar.setStyleSheet("QStatusBar{padding-left:8px;background:rgba(255,0,0,255);color:black;font-weight:bold;}")

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
