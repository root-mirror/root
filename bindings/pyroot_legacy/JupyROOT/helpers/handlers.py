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

from ctypes import CDLL, c_char_p
from threading import Thread
from time import sleep as timeSleep
from sys import platform
from os import path
import weakref
import sys
if sys.hexversion >= 0x3000000:
    import queue
else:
    import Queue as queue

_lib = CDLL(path.join(path.dirname(path.dirname(path.dirname(__file__))), 'libJupyROOT.so'))

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
        for cfunc in [_lib.JupyROOTExecutorHandler_GetStdout,
                      _lib.JupyROOTExecutorHandler_GetStderr]:
           cfunc.restype = c_char_p
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

    def Decode(self, obj):
        import sys
        if sys.version_info >= (3, 0):
            return obj.decode('utf-8')
        else:
            return obj

    def GetStdout(self):
       return self.Decode(_lib.JupyROOTExecutorHandler_GetStdout())

    def GetStderr(self):
       return self.Decode(_lib.JupyROOTExecutorHandler_GetStderr())

    def GetStreamsDicts(self):
       out = self.GetStdout()
       err = self.GetStderr()
       outDict = {'name': 'stdout', 'text': out} if out != "" else None
       errDict = {'name': 'stderr', 'text': err} if err != "" else None
       return outDict,errDict

class Poller(Thread):
    def __init__(self, runner_obj, name):
        Thread.__init__(self, group=None, target=None, name=name)
        self.poll = True
        self.ro_ref = weakref.ref(runner_obj)
    def run(self):
        while self.poll:
            work_item_argument = self.ro_ref().argument_queue.get()
            if work_item_argument is not None:
                self.ro_ref().is_running = True
                self.ro_ref().Run(work_item_argument)
                self.ro_ref().is_running = False
            else:
                self.poll = False
        return

class Runner(object):
    ''' Asynchrously run functions
    >>> import time
    >>> def f(code):
    ...    print(code)
    >>> r= Runner(f)
    >>> r.Run("ss")
    ss
    >>> r.AsyncRun("ss");time.sleep(1)
    ss
    >>> def g(msg):
    ...    time.sleep(.5)
    ...    print(msg)
    >>> r= Runner(g)
    >>> r.AsyncRun("Asynchronous");print("Synchronous");time.sleep(1)
    Synchronous
    Asynchronous
    >>> r.AsyncRun("Asynchronous"); print(r.HasFinished())
    False
    >>> time.sleep(1)
    Asynchronous
    >>> print(r.HasFinished())
    True
    >>> r.Stop()
    '''
    def __init__(self, function):
        self.function = function
        self.is_running = False
        self.argument_queue = queue.Queue()
        self.poller = Poller(runner_obj=self, name = "JupyROOT Runner Thread")
        self.poller.start()

    def __del__(self):
        if self.poller.is_alive():
            self.Stop()
        self.poller.join()

    def Run(self, argument):
        return self.function(argument)

    def AsyncRun(self, argument):
        self.is_running = True
        self.argument_queue.put(argument)

    def Wait(self):
        while self.is_running: pass

    def HasFinished(self):
        if self.is_running: return False
        return True

    def Stop(self):
        self.Wait()
        self.argument_queue.put(None)
        self.Wait()


class JupyROOTDeclarer(Runner):
    ''' Asynchrously execute declarations
    >>> import ROOT
    >>> d = JupyROOTDeclarer()
    >>> d.Run("int f(){return 3;}".encode("utf-8"))
    1
    >>> ROOT.f()
    3
    >>> d.Stop()
    '''
    def __init__(self):
       super(JupyROOTDeclarer, self).__init__(_lib.JupyROOTDeclarer)

class JupyROOTExecutor(Runner):
    r''' Asynchrously execute process lines
    >>> import ROOT
    >>> d = JupyROOTExecutor()
    >>> d.Run('cout << "Here am I" << endl;'.encode("utf-8"))
    1
    >>> d.Stop()
    '''
    def __init__(self):
       super(JupyROOTExecutor, self).__init__(_lib.JupyROOTExecutor)

def RunAsyncAndPrint(executor, code, ioHandler, printFunction, silent = False, timeout = 0.1):
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
