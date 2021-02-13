# Author: Enric Tejedor CERN  04/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

import sys
import time

from cppyy.gbl import gSystem, gInterpreter

from libROOTPythonizations import InitApplication, InstallGUIEventInputHook


class PyROOTApplication(object):
    """
    Application class for PyROOT.
    Configures the interactive usage of ROOT from Python.
    """

    def __init__(self, config, is_ipython):
        # Construct a TApplication for PyROOT
        InitApplication(config.IgnoreCommandLineOptions)

        self._is_ipython = is_ipython

    @staticmethod
    def _ipython_config():
        # Integrate IPython >= 5 with ROOT's event loop
        # Check for new GUI events until there is some user input to process

        from IPython import get_ipython
        from IPython.terminal import pt_inputhooks

        def inputhook(context):
            while not context.input_is_ready():
                gSystem.ProcessEvents()
                time.sleep(0.01)

        pt_inputhooks.register('ROOT', inputhook)

        ipy = get_ipython()
        if ipy:
            get_ipython().run_line_magic('gui', 'ROOT')

    @staticmethod
    def _inputhook_config():
        # PyOS_InputHook-based mechanism
        # Point to a function which will be called when Python's interpreter prompt
        # is about to become idle and wait for user input from the terminal
        InstallGUIEventInputHook()

    @staticmethod
    def _set_display_hook():
        # Set the display hook

        orig_dhook = sys.displayhook

        def displayhook(v):
            # sys.displayhook is called on the result of evaluating an expression entered
            # in an interactive Python session.
            # Therefore, this function will call EndOfLineAction after each interactive
            # command (to update display etc.)
            gInterpreter.EndOfLineAction()
            return orig_dhook(v)

        sys.displayhook = displayhook

    def init_graphics(self):
        """Configure ROOT graphics to be used interactively"""

        # Note that we only end up in this function if gROOT.IsBatch() is false
        import __main__
        if self._is_ipython and 'IPython' in sys.modules and sys.modules['IPython'].version_info[0] >= 5:
            # ipython and notebooks, register our event processing with their hooks
            self._ipython_config()
        elif sys.flags.interactive == 1 or not hasattr(__main__, '__file__') or gSystem.InheritsFrom('TMacOSXSystem'):
            # Python in interactive mode, use the PyOS_InputHook to call our event processing
            # - sys.flags.interactive checks for the -i flags passed to python
            # - __main__ does not have the attribute __file__ if the Python prompt is started directly
            # - MacOS does not allow to run a second thread to process events, fall back to the input hook
            self._inputhook_config()
        else:
            # Python in script mode, start a separate thread for the event processing
            def _process_root_events(self):
                while self.keep_polling:
                    gSystem.ProcessEvents()
                    time.sleep(0.01)
            import threading
            self.keep_polling = True # Used to shut down the thread safely at teardown time
            update_thread = threading.Thread(None, _process_root_events, None, (self,))
            self.process_root_events = update_thread # The thread is joined at teardown time
            update_thread.daemon = True
            update_thread.start()

        self._set_display_hook()
