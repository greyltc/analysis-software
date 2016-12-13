#!/usr/bin/env python

import sys
from batch_iv_analysis import MainWindow
from PyQt5.QtWidgets import QApplication

def main(args=None):
    """# a tool for analysing solar cell i-v curves
	# written by Grey Christoforo <first name [at] last name [not] net>
	# please cite our work if you can!
	# DOI: 10.3390/photonics2041101
	"""
	
    if args is None:
        args = sys.argv[1:]

    app = QApplication(sys.argv)
    analysis = MainWindow()
    analysis.show()
    sys.exit(app.exec_())

    # Do argument parsing here (eg. with argparse) and anything else
    # you want your project to do.

if __name__ == "__main__":
    main()