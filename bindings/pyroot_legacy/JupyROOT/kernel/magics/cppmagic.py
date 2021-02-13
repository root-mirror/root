# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------------
#  Authors: Omar Zapata <Omar.Zapata@cern.ch> http://oproject.org
#           Danilo Piparo <Danilo.Piparo@cern.ch> CERN
#           Enric Tejedor enric.tejedor.saavedra@cern.ch> CERN
#-----------------------------------------------------------------------------

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from metakernel import Magic, option

import sys

#NOTE:actually JupyROOT is not capturing the error on %%cpp -d if the function is wrong
class CppMagics(Magic):
    def __init__(self, kernel):
        super(CppMagics, self).__init__(kernel)
    @option(
        '-a', '--aclic', action='store', default="default", help='Compile code with ACLiC.'
    )
    @option(
        '-d', '--declare', action='store', default=None, help='Declare functions and/or classes.'
    )
    def cell_cpp(self, args):
        '''Executes the content of the cell as C++ code.'''
        if self.code.strip():
             execFunc = None
             if args == '-a':
                 execFunc = self.kernel.ACLiC
             elif args == '-d':
                 execFunc = self.kernel.Declarer.Run
             else: # normal flow
                self.kernel.do_execute_direct(str(self.code))
                self.evaluate = False
                return

             self.kernel.ioHandler.Clear()
             self.kernel.ioHandler.InitCapture()
             execFunc(self.code.encode('utf8'))
             self.kernel.ioHandler.EndCapture()
             self.kernel.print_output(self.kernel.ioHandler)

        self.evaluate = False

def register_magics(kernel):
    kernel.register_magics(CppMagics)

