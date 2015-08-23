#!/usr/bin/env python

# ROOT command line tools: rootmv
# Author: Julien Ripoche
# Mail: julien.ripoche@u-psud.fr
# Date: 20/08/15

"""Command line to move objects from ROOT files to another"""

import cmdLineUtils
import sys

# Help strings
COMMAND_HELP = "Move objects from ROOT files to another"

EPILOG = """Examples:
- roomv source.root:hist* dest.root
  Move all histograms whose named starts with 'hist' from 'source.root' to 'dest.root'.

- roomv source1.root:hist1 source2.root:hist2 dest.root
  Move histograms 'hist1' from 'source1.root' and 'hist2' from 'source2.root' to 'dest.root'.

- roomv --recreate source.root:hist dest.root
  Recreate the destination file 'dest.root' and move the histogram named 'hist' from 'source.root' into it.

- roomv -c 1 source.root:hist dest.root
  Change the compression level of the destination file 'dest.root' and move the histogram named 'hist' from 'source.root' into it. For more information about compression settings of ROOT file, please look at the reference guide available on the ROOT site.
"""

def execute():
    # Collect arguments with the module argparse
    parser = cmdLineUtils.getParserSourceDest(COMMAND_HELP, EPILOG)
    parser.add_argument("-c","--compress", type=int, help=cmdLineUtils.COMPRESS_HELP)
    parser.add_argument("-i","--interactive", help=cmdLineUtils.INTERACTIVE_HELP, action="store_true")
    parser.add_argument("--recreate", help=cmdLineUtils.RECREATE_HELP, action="store_true")

    # Put arguments in shape
    sourceList, destFileName, destPathSplit, optDict = cmdLineUtils.getSourceDestListOptDict(parser)

    # Process rootMv
    return cmdLineUtils.rootMv(sourceList, destFileName, destPathSplit, \
                               compress=optDict["compress"], interactive=optDict["interactive"], \
                               recreate=optDict["recreate"])

sys.exit(execute())
