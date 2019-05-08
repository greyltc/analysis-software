import os, sys

"""
a tool for analyzing solar cell i-v curves
written by Grey Christoforo <grey@mutovis.com>
"""
def main():
  if ('MUTOVIS_CLI_ANALYSIS' in os.environ) or ('-cli' in sys.argv[0]):
    import batch_iv_analysis.cli as cli
    cli.handle_cli()
  else:
    import batch_iv_analysis.gui as gui
    from batch_iv_analysis.ivAnalyzer import ivAnalyzer
    a = ivAnalyzer()
    gui.runGUI(a,None)

if __name__ == "__main__":
    main()
