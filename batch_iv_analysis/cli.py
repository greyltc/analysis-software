from batch_iv_analysis.ivAnalyzer import ivAnalyzer
import argparse
import sys


def runCLI(analyzer:ivAnalyzer, args):
    analyzer.setup()
    print("Got args:", args)


def handle_cli(sysargs):
    parser = argparse.ArgumentParser(description="Process some iv data.")

    parser.add_argument("-f", "--files", default=None, type=argparse.FileType("r"), help="File(s) to analyze.")
    parser.add_argument("-g", "--gui", default=False, action="store_true", help="Run with GUI.")
    parser.add_argument("-s", "--no-sloppy", dest="sloppyMath", default=True, action="store_false", help="Don't do sloppy math (slower).")
    parser.add_argument("-w", "--workers", default=4, type=int, help="Multiprocessing control. w=0 disables it. w>0 enables it with w workers.")

    args = parser.parse_args(sysargs[1:])

    if args.gui == False and args.files == None:
        raise (ValueError("Command line interface (cli) mode needs files to process"))

    a = ivAnalyzer(beFastAndSloppy=args.sloppyMath, poolWorkers=args.workers)
    if args.gui == False:
        runCLI(a, args)
    else:
        import batch_iv_analysis.gui as gui

        gui.runGUI(a, args)


if __name__ == "__main__":
    handle_cli(sys.argv)
