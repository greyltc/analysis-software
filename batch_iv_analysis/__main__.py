import os, sys
import batch_iv_analysis.cli as cli

"""
a tool for analyzing solar cell i-v curves
written by Grey Christoforo <grey@mutovis.com>
"""

if ('MUTOVIS_CLI_ANALYSIS' in os.environ) or ('-cli' in sys.argv[0]):
  cli.handle_cli(sys.argv)
else:
  cli.handle_cli(sys.argv+["--gui"])
