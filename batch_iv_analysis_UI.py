# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'batch-iv-analysis.ui'
#
# Created: Mon Mar 31 20:19:50 2014
#      by: PyQt4 UI code generator 4.10
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_batch_iv_analysis(object):
    def setupUi(self, batch_iv_analysis):
        batch_iv_analysis.setObjectName(_fromUtf8("batch_iv_analysis"))
        batch_iv_analysis.resize(1238, 588)
        self.centralwidget = QtGui.QWidget(batch_iv_analysis)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.tableWidget = QtGui.QTableWidget(self.centralwidget)
        self.tableWidget.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.tableWidget.setAlternatingRowColors(True)
        self.tableWidget.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.tableWidget.setSelectionBehavior(QtGui.QAbstractItemView.SelectItems)
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setObjectName(_fromUtf8("tableWidget"))
        self.tableWidget.setRowCount(0)
        self.tableWidget.horizontalHeader().setCascadingSectionResizes(False)
        self.tableWidget.horizontalHeader().setDefaultSectionSize(92)
        self.tableWidget.horizontalHeader().setSortIndicatorShown(False)
        self.tableWidget.verticalHeader().setVisible(False)
        self.tableWidget.verticalHeader().setDefaultSectionSize(30)
        self.tableWidget.verticalHeader().setSortIndicatorShown(False)
        self.gridLayout.addWidget(self.tableWidget, 0, 0, 1, 1)
        batch_iv_analysis.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(batch_iv_analysis)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1238, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuTools = QtGui.QMenu(self.menubar)
        self.menuTools.setObjectName(_fromUtf8("menuTools"))
        batch_iv_analysis.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(batch_iv_analysis)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        batch_iv_analysis.setStatusBar(self.statusbar)
        self.actionQuit = QtGui.QAction(batch_iv_analysis)
        self.actionQuit.setObjectName(_fromUtf8("actionQuit"))
        self.actionOpen = QtGui.QAction(batch_iv_analysis)
        self.actionOpen.setObjectName(_fromUtf8("actionOpen"))
        self.actionClear_Graph = QtGui.QAction(batch_iv_analysis)
        self.actionClear_Graph.setObjectName(_fromUtf8("actionClear_Graph"))
        self.actionClear_Table = QtGui.QAction(batch_iv_analysis)
        self.actionClear_Table.setObjectName(_fromUtf8("actionClear_Table"))
        self.actionFsadf = QtGui.QAction(batch_iv_analysis)
        self.actionFsadf.setObjectName(_fromUtf8("actionFsadf"))
        self.actionSet_Bounds = QtGui.QAction(batch_iv_analysis)
        self.actionSet_Bounds.setObjectName(_fromUtf8("actionSet_Bounds"))
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionClear_Graph)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionQuit)
        self.menuTools.addAction(self.actionSet_Bounds)
        self.menuTools.addAction(self.actionClear_Table)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuTools.menuAction())

        self.retranslateUi(batch_iv_analysis)
        QtCore.QObject.connect(self.actionQuit, QtCore.SIGNAL(_fromUtf8("triggered()")), batch_iv_analysis.close)
        QtCore.QMetaObject.connectSlotsByName(batch_iv_analysis)

    def retranslateUi(self, batch_iv_analysis):
        batch_iv_analysis.setWindowTitle(_translate("batch_iv_analysis", "batch-iv-analysis", None))
        self.tableWidget.setSortingEnabled(False)
        self.menuFile.setTitle(_translate("batch_iv_analysis", "File", None))
        self.menuTools.setTitle(_translate("batch_iv_analysis", "Tools", None))
        self.actionQuit.setText(_translate("batch_iv_analysis", "Quit", None))
        self.actionOpen.setText(_translate("batch_iv_analysis", "Open", None))
        self.actionClear_Graph.setText(_translate("batch_iv_analysis", "Save", None))
        self.actionClear_Table.setText(_translate("batch_iv_analysis", "Clear Table", None))
        self.actionFsadf.setText(_translate("batch_iv_analysis", "fsadf", None))
        self.actionSet_Bounds.setText(_translate("batch_iv_analysis", "Set Bounds", None))

