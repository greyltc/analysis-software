import os, sys
sys.path.append(os.path.dirname(__file__))

import ivAnalyzer
import argparse

def main():
  """# a tool for analyzing solar cell i-v curves
  # written by Grey Christoforo <first name [at] last name [not] net>
  # please cite our work if you can!
  # DOI: 10.3390/photonics2041101
  """
  
  parser = argparse.ArgumentParser(description='Process some iv data.')
  
  parser.add_argument('-f', '--files', type=argparse.FileType('r'), help='File(s) to analyze.')
  parser.add_argument('-c', '--cli', dest='cli', action='store_true', help='Run without GUI.')
  parser.add_argument('-s', '--sloppy', dest='sloppyMath', action='store_true', help='Do sloppy math (faster).')
  parser.add_argument('-w', '--workers', dest='workers', type=int, help='Multiprocessing control. w > 1 enables it with w workers.')
  
  args = parser.parse_args()    
  
  if args.cli is True:
    import cli
    if args.sloppyMath is None:
      args.sloppyMath = False
    if args.workers is None:
      args.workers = 1
      multiprocess = False
    elif args.workers is 1:
      multiprocess = False
    else:
      multiprocess = True
    a = ivAnalyzer.ivAnalyzer(beFastAndSloppy=args.sloppyMath, multiprocess=multiprocess, poolWorkers=args.workers)
    cli.runCLI(a,args)
  else:
    from PyQt5.QtCore import QSettings
    settings = QSettings("greyltc", "batch-iv-analysis") # gui mode reads QSettings
    
    if args.sloppyMath is None:
      if settings.contains('fastAndSloppy'):
        args.sloppyMath = self.settings.value('fastAndSloppy') == 'true'
      else:
        args.sloppyMath = False
    if args.workers is None:
      args.workers = 1
      
    if args.workers is None:
      if settings.contains('multiprocessing'):
        multiprocess = self.settings.value('multiprocessing') == 'true'
      else:
        multiprocess = False
      if settings.contains('threads'):
        args.workers = int(self.settings.value('threads'))
    elif args.workers is 1:
      multiprocess = False
    else:
      multiprocess = True

    a = ivAnalyzer.ivAnalyzer(beFastAndSloppy=args.sloppyMath, multiprocess=multiprocess, poolWorkers=args.workers)
    import gui
    gui.runGUI(a,args)
    


if __name__ == '__main__':
  main()