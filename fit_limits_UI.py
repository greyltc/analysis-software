# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'fit-limits.ui'
#
# Created: Wed Mar 26 19:51:09 2014
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

class Ui_fit_limits(object):
    def setupUi(self, fit_limits):
        fit_limits.setObjectName(_fromUtf8("fit_limits"))
        fit_limits.resize(400, 300)
        self.buttonBox = QtGui.QDialogButtonBox(fit_limits)
        self.buttonBox.setGeometry(QtCore.QRect(30, 250, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.label = QtGui.QLabel(fit_limits)
        self.label.setGeometry(QtCore.QRect(10, 50, 81, 21))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_2 = QtGui.QLabel(fit_limits)
        self.label_2.setGeometry(QtCore.QRect(10, 90, 81, 21))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_3 = QtGui.QLabel(fit_limits)
        self.label_3.setGeometry(QtCore.QRect(10, 130, 81, 21))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_4 = QtGui.QLabel(fit_limits)
        self.label_4.setGeometry(QtCore.QRect(10, 170, 81, 21))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.label_5 = QtGui.QLabel(fit_limits)
        self.label_5.setGeometry(QtCore.QRect(10, 210, 41, 21))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.label_6 = QtGui.QLabel(fit_limits)
        self.label_6.setGeometry(QtCore.QRect(100, 10, 81, 21))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.label_7 = QtGui.QLabel(fit_limits)
        self.label_7.setGeometry(QtCore.QRect(250, 10, 111, 21))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.lineEdit = QtGui.QLineEdit(fit_limits)
        self.lineEdit.setGeometry(QtCore.QRect(100, 40, 113, 33))
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        self.lineEdit_2 = QtGui.QLineEdit(fit_limits)
        self.lineEdit_2.setGeometry(QtCore.QRect(250, 40, 113, 33))
        self.lineEdit_2.setObjectName(_fromUtf8("lineEdit_2"))
        self.lineEdit_3 = QtGui.QLineEdit(fit_limits)
        self.lineEdit_3.setGeometry(QtCore.QRect(100, 80, 113, 33))
        self.lineEdit_3.setObjectName(_fromUtf8("lineEdit_3"))
        self.lineEdit_4 = QtGui.QLineEdit(fit_limits)
        self.lineEdit_4.setGeometry(QtCore.QRect(100, 120, 113, 33))
        self.lineEdit_4.setObjectName(_fromUtf8("lineEdit_4"))
        self.lineEdit_5 = QtGui.QLineEdit(fit_limits)
        self.lineEdit_5.setGeometry(QtCore.QRect(100, 160, 113, 33))
        self.lineEdit_5.setObjectName(_fromUtf8("lineEdit_5"))
        self.lineEdit_6 = QtGui.QLineEdit(fit_limits)
        self.lineEdit_6.setGeometry(QtCore.QRect(250, 80, 113, 33))
        self.lineEdit_6.setObjectName(_fromUtf8("lineEdit_6"))
        self.lineEdit_7 = QtGui.QLineEdit(fit_limits)
        self.lineEdit_7.setGeometry(QtCore.QRect(250, 120, 113, 33))
        self.lineEdit_7.setObjectName(_fromUtf8("lineEdit_7"))
        self.lineEdit_8 = QtGui.QLineEdit(fit_limits)
        self.lineEdit_8.setGeometry(QtCore.QRect(250, 160, 113, 33))
        self.lineEdit_8.setObjectName(_fromUtf8("lineEdit_8"))
        self.lineEdit_9 = QtGui.QLineEdit(fit_limits)
        self.lineEdit_9.setGeometry(QtCore.QRect(100, 200, 113, 33))
        self.lineEdit_9.setObjectName(_fromUtf8("lineEdit_9"))
        self.lineEdit_10 = QtGui.QLineEdit(fit_limits)
        self.lineEdit_10.setGeometry(QtCore.QRect(250, 200, 113, 33))
        self.lineEdit_10.setObjectName(_fromUtf8("lineEdit_10"))

        self.retranslateUi(fit_limits)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), fit_limits.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), fit_limits.reject)
        QtCore.QMetaObject.connectSlotsByName(fit_limits)

    def retranslateUi(self, fit_limits):
        fit_limits.setWindowTitle(_translate("fit_limits", "Set Parameter Bounds", None))
        self.label.setText(_translate("fit_limits", "I_0 [A]", None))
        self.label_2.setText(_translate("fit_limits", "I_ph [A]", None))
        self.label_3.setText(_translate("fit_limits", "R_s [Ohm]", None))
        self.label_4.setText(_translate("fit_limits", "R_sh [Ohm]", None))
        self.label_5.setText(_translate("fit_limits", "n", None))
        self.label_6.setText(_translate("fit_limits", "Minimum", None))
        self.label_7.setText(_translate("fit_limits", "Maximum", None))
        self.lineEdit.setText(_translate("fit_limits", "0", None))
        self.lineEdit_2.setText(_translate("fit_limits", "1e-3", None))
        self.lineEdit_3.setText(_translate("fit_limits", "-1e-1", None))
        self.lineEdit_4.setText(_translate("fit_limits", "0", None))
        self.lineEdit_5.setText(_translate("fit_limits", "0", None))
        self.lineEdit_6.setText(_translate("fit_limits", "2e-1", None))
        self.lineEdit_7.setText(_translate("fit_limits", "400", None))
        self.lineEdit_8.setText(_translate("fit_limits", "1e7", None))
        self.lineEdit_9.setText(_translate("fit_limits", "-2", None))
        self.lineEdit_10.setText(_translate("fit_limits", "7", None))

