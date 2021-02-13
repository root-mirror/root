#-----------------------------------------------------------------------------
#  Author: Danilo Piparo <Danilo.Piparo@cern.ch> CERN
#  Author: Enric Tejedor <enric.tejedor.saavedra@cern.ch> CERN
#-----------------------------------------------------------------------------

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from IPython.core.magic import (Magics, magics_class, cell_magic)
from IPython.core.magic_arguments import (argument, magic_arguments, parse_argstring)
from JupyROOT.helpers import utils


@magics_class
class CppMagics(Magics):
    @cell_magic
    @magic_arguments()
    @argument('-a', '--aclic', action="store_true", help='Compile code with ACLiC.')
    @argument('-d', '--declare', action="store_true", help='Declare functions and/or classes.')
    def cpp(self, line, cell):
        '''Executes the content of the cell as C++ code.'''
        args = parse_argstring(self.cpp, line)
        if args.aclic:
            utils.invokeAclic(cell)
        elif args.declare:
            utils.declareCppCode(cell)
        else:
            utils.processMagicCppCode(cell)

def load_ipython_extension(ipython):
    ipython.register_magics(CppMagics)
