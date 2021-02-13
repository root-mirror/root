# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------------
#  Authors: Danilo Piparo
#           Omar Zapata <Omar.Zapata@cern.ch> http://oproject.org
#-----------------------------------------------------------------------------

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from threading import Thread
from time import sleep as timeSleep
from sys import platform
from os import path
import sys
if sys.hexversion >= 0x3000000:
    import queue
else:
    import Queue as queue

from JupyROOT import helpers

# import libJupyROOT with Python version number
import importlib
major, minor = sys.version_info[0:2]
libjupyroot_mod_name = 'libJupyROOT{}_{}'.format(major, minor)
_lib = importlib.import_module(libjupyroot_mod_name)


class IOHandler(object):
    r'''Class used to capture output from C/C++ libraries.
    >>> import sys
    >>> h = IOHandler()
    >>> h.GetStdout()
    ''
    >>> h.GetStderr()
    ''
    >>> h.GetStreamsDicts()
    (None, None)
    >>> del h
    '''
    def __init__(self):
        _lib.JupyROOTExecutorHandler_Ctor()

    def __del__(self):
        _lib.JupyROOTExecutorHandler_Dtor()

    def Clear(self):
        _lib.JupyROOTExecutorHandler_Clear()

    def Poll(self):
        _lib.JupyROOTExecutorHandler_Poll()

    def InitCapture(self):
        _lib.JupyROOTExecutorHandler_InitCapture()

    def EndCapture(self):
        _lib.JupyROOTExecutorHandler_EndCapture()

    def GetStdout(self):
       return _lib.JupyROOTExecutorHandler_GetStdout()

    def GetStderr(self):
       return _lib.JupyROOTExecutorHandler_GetStderr()

    def GetStreamsDicts(self):
       out = self.GetStdout()
       err = self.GetStderr()
       outDict = {'name': 'stdout', 'text': out} if out != "" else None
       errDict = {'name': 'stderr', 'text': err} if err != "" else None
       return outDict,errDict

class Poller(Thread):
    def __init__(self):
        Thread.__init__(self, group=None, target=None, name="JupyROOT Poller Thread")
        self.poll = True
        self.is_running = False
        self.queue = queue.Queue()

    def run(self):
        while self.poll:
            work_item = self.queue.get()
            if work_item is not None:
                function, argument = work_item
                self.is_running = True
                function(argument)
                self.is_running = False
            else:
                self.poll = False

    def Stop(self):
        if self.is_alive():
            self.queue.put(None)
            self.join() 

class Runner(object):
    ''' Asynchrously run functions
    >>> import time
    >>> def f(code):
    ...    print(code)
    >>> p = Poller(); p.start()
    >>> r= Runner(f, p)
    >>> r.Run("ss")
    ss
    >>> r.AsyncRun("ss");time.sleep(1)
    ss
    >>> def g(msg):
    ...    time.sleep(.5)
    ...    print(msg)
    >>> r= Runner(g, p)
    >>> r.AsyncRun("Asynchronous");print("Synchronous");time.sleep(1)
    Synchronous
    Asynchronous
    >>> r.AsyncRun("Asynchronous"); print(r.HasFinished())
    False
    >>> time.sleep(1)
    Asynchronous
    >>> print(r.HasFinished())
    True
    >>> p.Stop()
    '''
    def __init__(self, function, poller):
        self.function = function
        self.poller = poller

    def Run(self, argument):
        return self.function(argument)

    def AsyncRun(self, argument):
        self.poller.is_running = True
        self.poller.queue.put((self.Run, argument))

    def Wait(self):
        while self.poller.is_running:
            timeSleep(.1)

    def HasFinished(self):
        return not self.poller.is_running

class JupyROOTDeclarer(Runner):
    ''' Asynchrously execute declarations
    >>> import ROOT
    >>> p = Poller(); p.start()
    >>> d = JupyROOTDeclarer(p)
    >>> d.Run("int f(){return 3;}")
    1
    >>> ROOT.f()
    3
    >>> p.Stop()
    '''
    def __init__(self, poller):
       super(JupyROOTDeclarer, self).__init__(_lib.JupyROOTDeclarer, poller)

class JupyROOTExecutor(Runner):
    r''' Asynchrously execute process lines
    >>> import ROOT
    >>> p = Poller(); p.start()
    >>> d = JupyROOTExecutor(p)
    >>> d.Run('cout << "Here am I" << endl;')
    1
    >>> p.Stop()
    '''
    def __init__(self, poller):
       super(JupyROOTExecutor, self).__init__(_lib.JupyROOTExecutor, poller)

def display_drawables(displayFunction):
    drawers = helpers.utils.GetDrawers()
    for drawer in drawers:
        for dobj in drawer.GetDrawableObjects():
            displayFunction(dobj)

class JupyROOTDisplayer(Runner):
    ''' Display all canvases'''
    def __init__(self, poller):
       super(JupyROOTDisplayer, self).__init__(display_drawables, poller)

def RunAsyncAndPrint(executor, code, ioHandler, printFunction, displayFunction, silent = False, timeout = 0.1):
    ioHandler.Clear()
    ioHandler.InitCapture()
    executor.AsyncRun(code)
    while not executor.HasFinished():
        ioHandler.Poll()
        if not silent:
            printFunction(ioHandler)
            ioHandler.Clear()
        if executor.HasFinished(): break
        timeSleep(.1)
    executor.Wait()
    ioHandler.EndCapture()

def Display(displayer, displayFunction):
    displayer.AsyncRun(displayFunction)
    displayer.Wait()

