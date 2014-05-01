# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'batch-iv-analysis.ui'
#
# Created: Wed Apr 30 20:20:49 2014
#      by: PyQt4 UI code generator 4.9.6
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
        self.tableWidget.horizontalHeader().setSortIndicatorShown(True)
        self.tableWidget.verticalHeader().setVisible(False)
        self.tableWidget.verticalHeader().setDefaultSectionSize(30)
        self.tableWidget.verticalHeader().setSortIndicatorShown(True)
        self.gridLayout.addWidget(self.tableWidget, 0, 0, 1, 1)
        batch_iv_analysis.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(batch_iv_analysis)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1238, 21))
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
        self.actionSave = QtGui.QAction(batch_iv_analysis)
        self.actionSave.setObjectName(_fromUtf8("actionSave"))
        self.actionClear_Table = QtGui.QAction(batch_iv_analysis)
        self.actionClear_Table.setEnabled(True)
        self.actionClear_Table.setObjectName(_fromUtf8("actionClear_Table"))
        self.actionFsadf = QtGui.QAction(batch_iv_analysis)
        self.actionFsadf.setObjectName(_fromUtf8("actionFsadf"))
        self.actionSet_Bounds = QtGui.QAction(batch_iv_analysis)
        self.actionSet_Bounds.setObjectName(_fromUtf8("actionSet_Bounds"))
        self.actionWatch = QtGui.QAction(batch_iv_analysis)
        self.actionWatch.setObjectName(_fromUtf8("actionWatch"))
        self.actionEnable_Watching = QtGui.QAction(batch_iv_analysis)
        self.actionEnable_Watching.setCheckable(True)
        self.actionEnable_Watching.setChecked(False)
        self.actionEnable_Watching.setObjectName(_fromUtf8("actionEnable_Watching"))
        self.actionWatch_2 = QtGui.QAction(batch_iv_analysis)
        self.actionWatch_2.setObjectName(_fromUtf8("actionWatch_2"))
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionWatch_2)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionSave)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionQuit)
        self.menuTools.addAction(self.actionClear_Table)
        self.menuTools.addAction(self.actionEnable_Watching)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuTools.menuAction())

        self.retranslateUi(batch_iv_analysis)
        QtCore.QObject.connect(self.actionQuit, QtCore.SIGNAL(_fromUtf8("triggered()")), batch_iv_analysis.close)
        QtCore.QMetaObject.connectSlotsByName(batch_iv_analysis)

    def retranslateUi(self, batch_iv_analysis):
        batch_iv_analysis.setWindowTitle(_translate("batch_iv_analysis", "batch-iv-analysis", None))
        self.tableWidget.setSortingEnabled(True)
        self.menuFile.setTitle(_translate("batch_iv_analysis", "File", None))
        self.menuTools.setTitle(_translate("batch_iv_analysis", "Tools", None))
        self.actionQuit.setText(_translate("batch_iv_analysis", "Quit", None))
        self.actionQuit.setShortcut(_translate("batch_iv_analysis", "Ctrl+Q", None))
        self.actionOpen.setText(_translate("batch_iv_analysis", "Open", None))
        self.actionOpen.setShortcut(_translate("batch_iv_analysis", "Ctrl+O", None))
        self.actionSave.setText(_translate("batch_iv_analysis", "Export", None))
        self.actionSave.setShortcut(_translate("batch_iv_analysis", "Ctrl+S", None))
        self.actionClear_Table.setText(_translate("batch_iv_analysis", "Clear Table", None))
        self.actionClear_Table.setShortcut(_translate("batch_iv_analysis", "Ctrl+Backspace", None))
        self.actionFsadf.setText(_translate("batch_iv_analysis", "fsadf", None))
        self.actionSet_Bounds.setText(_translate("batch_iv_analysis", "Set Bounds", None))
        self.actionWatch.setText(_translate("batch_iv_analysis", "Watch", None))
        self.actionWatch.setShortcut(_translate("batch_iv_analysis", "Ctrl+W", None))
        self.actionEnable_Watching.setText(_translate("batch_iv_analysis", "Enable Watching", None))
        self.actionEnable_Watching.setShortcut(_translate("batch_iv_analysis", "Ctrl+E", None))
        self.actionWatch_2.setText(_translate("batch_iv_analysis", "Watch", None))
        self.actionWatch_2.setShortcut(_translate("batch_iv_analysis", "Ctrl+W", None))

