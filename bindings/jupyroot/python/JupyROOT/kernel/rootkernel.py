#!/usr/bin/env python
# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------------
#  Authors: Omar Zapata <Omar.Zapata@cern.ch> http://oproject.org
#           Danilo Piparo <Danilo.Piparo@cern.ch> CERN
#           Enric Tejedor <enric.tejedor.saavedra@cern.ch> CERN
#-----------------------------------------------------------------------------

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from __future__ import print_function

import sys

try:
    from metakernel import MetaKernel
    from metakernel.display import HTML
except ImportError:
    raise Exception("Error: package metakernel not found.(install it running 'pip install metakernel')")

import ROOT

from JupyROOT.helpers.utils import setStyle, invokeAclic, GetDrawers
from JupyROOT.helpers.handlers import RunAsyncAndPrint, Display
from JupyROOT.helpers.cppcompleter import CppCompleter
from JupyROOT.kernel.utils import GetIOHandler, GetPoller, GetExecutor, GetDeclarer, GetDisplayer, MagicLoader

import IPython

# We want iPython to take over the graphics
ROOT.gROOT.SetBatch()

def Debug(msg):
     print('Kernel main: %r' % msg, file=sys.__stderr__)

class ROOTKernel(MetaKernel):
    # These two regexes are considered by the parser of the metakernel
    # there is no need to create one explicitly
    identifier_regex = r'(?:\w(?:\w|\.|->|::|\d)*)'
    func_call_regex = r'(?:\w(?:(?:\w|\.|->|::|\d))*)\([^\)\()]*\Z'

    implementation = 'ROOT'
    implementation_version = '1.0'
    language = 'c++'
    language_version = '0.1'
    language_info = {'name': 'c++',
                     'codemirror_mode': 'text/x-c++src',
                     'mimetype': ' text/x-c++src',
                     'file_extension': '.C'}
    banner = "ROOT Kernel"

    def __init__(self,**kwargs):

        MetaKernel.__init__(self,**kwargs)
        setStyle()
        self.ioHandler = GetIOHandler()
        self.Poller    = GetPoller()
        self.Executor  = GetExecutor(self.Poller)
        self.Declarer  = GetDeclarer(self.Poller) #required for %%cpp -d magic
        self.Displayer = GetDisplayer(self.Poller)
        self.ACLiC     = invokeAclic
        self.magicloader = MagicLoader(self)
        self.completer = CppCompleter()
        self.completer.activate()

    def __del__(self):
        self.Poller.Stop()

    def get_completions(self, info):
        return self.completer._completeImpl(info['code'])

    def print_output(self, handler):
        streamDicts = handler.GetStreamsDicts()
        for streamDict in filter(lambda d: None != d, streamDicts):
            self.send_response(self.iopub_socket, 'stream', streamDict)

    def do_execute_direct(self, code, silent=False):

        if not code.strip():
            return

        status = 'ok'
        try:
            RunAsyncAndPrint(self.Executor,
                             code,
                             self.ioHandler,
                             self.print_output,
                             silent,
                             .1)

            Display(self.Displayer, self.Display)

        except KeyboardInterrupt:
            ROOT.gROOT.SetInterrupt()
            status = 'interrupted'
            self.ioHandler.EndCapture()
        if not silent:
            self.print_output(self.ioHandler)

        traceback = None
        reply = {'status': status,
                'execution_count': self.execution_count,
                'payload': [],
                'user_expressions': {},
                }

        if status == 'interrupted':
            pass
        elif status == 'error':
            err = {
                'ename': 'ename',
                'evalue': 'evalue',
                'traceback': traceback,
            }
            self.send_response(self.iopub_socket, 'error', err)
            reply.update(err)
        elif status == 'ok':
            pass
        else:
            raise ValueError("Invalid status: %r" % status)


def main():
    """launch a root kernel"""
    try:
        from ipykernel.kernelapp import IPKernelApp
    except ImportError:
        from IPython.kernel.zmq.kernelapp import IPKernelApp
    IPKernelApp.launch_instance(kernel_class=ROOTKernel)

if __name__ == '__main__':
    main()
