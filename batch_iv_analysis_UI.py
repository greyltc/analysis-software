# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'batch-iv-analysis.ui'
#
# Created: Wed Mar 26 13:49:40 2014
#      by: PyQt4 UI code generator 4.10.2
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
        batch_iv_analysis.resize(1238, 599)
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
        self.tableWidget.horizontalHeader().setCascadingSectionResizes(True)
        self.tableWidget.horizontalHeader().setDefaultSectionSize(92)
        self.tableWidget.horizontalHeader().setSortIndicatorShown(False)
        self.tableWidget.verticalHeader().setDefaultSectionSize(30)
        self.tableWidget.verticalHeader().setSortIndicatorShown(True)
        self.gridLayout.addWidget(self.tableWidget, 0, 0, 1, 1)
        batch_iv_analysis.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(batch_iv_analysis)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1238, 27))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        batch_iv_analysis.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(batch_iv_analysis)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        batch_iv_analysis.setStatusBar(self.statusbar)
        self.actionQuit = QtGui.QAction(batch_iv_analysis)
        self.actionQuit.setObjectName(_fromUtf8("actionQuit"))
        self.actionOpen = QtGui.QAction(batch_iv_analysis)
        self.actionOpen.setObjectName(_fromUtf8("actionOpen"))
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionQuit)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(batch_iv_analysis)
        QtCore.QObject.connect(self.actionQuit, QtCore.SIGNAL(_fromUtf8("triggered()")), batch_iv_analysis.close)
        QtCore.QMetaObject.connectSlotsByName(batch_iv_analysis)

    def retranslateUi(self, batch_iv_analysis):
        batch_iv_analysis.setWindowTitle(_translate("batch_iv_analysis", "batch-iv-analysis", None))
        self.tableWidget.setSortingEnabled(False)
        self.menuFile.setTitle(_translate("batch_iv_analysis", "File", None))
        self.actionQuit.setText(_translate("batch_iv_analysis", "Quit", None))
        self.actionOpen.setText(_translate("batch_iv_analysis", "Open", None))

