import os, sys
sys.path.append(os.path.dirname(__file__))

from gui import runGUI
from cli import runCLI
import ivAnalyzer
import argparse

def main(args=None):
  """# a tool for analysing solar cell i-v curves
  # written by Grey Christoforo <first name [at] last name [not] net>
  # please cite our work if you can!
  # DOI: 10.3390/photonics2041101
  """
  a = ivAnalyzer.ivAnalyzer()

  if args is None:
    runGUI(a)
  else:
    runCLI(a,args)
    
def parseTehArgs():
  parser = argparse.ArgumentParser(description='Process some iv data.')
  
  parser.add_argument('-f', '--files',
                      help='File(s) to analyze.')
  args = parser.parse_args()
  
  if args.files is None:
    return None
  else:
    return args

if __name__ == '__main__':
  main(args=parseTehArgs())