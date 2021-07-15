# Authors:
# * Jonas Rembser 06/2021
# * Harshal Shende 06/2021

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################


r"""
/**
\class RooDataHist
\brief \parblock \endparblock
\htmlonly
<div class="pyrootbox">
\endhtmlonly

## PyROOT

Constructor of RooDataHist takes a RooCmdArg as argument also supports keyword arguments.
For example, the following code is equivalent in PyROOT:
\code{.py}
# Directly passing a RooCmdArg:
dh = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(x), ROOT.RooFit.Import("SampleA", histo))

# With keyword arguments:
dh = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(x), Import=("SampleA", histo))
\endcode

\htmlonly
</div>
\endhtmlonly
*/
"""

from ._utils import _kwargs_to_roocmdargs, _dict_to_std_map


class RooDataHist(object):
    def __init__(self, *args, **kwargs):
        # Redefinition of `RooDataHist` constructor for keyword arguments and converting python dict to std::map.
        # The keywords must correspond to the CmdArg of the constructor function.
        # The instances in dict must correspond to the template argument in std::map of the constructor.
        if len(args) > 4 and isinstance(args[4], dict):
            args = list(args)
            args[4] = _dict_to_std_map(args[4])

        args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
        self._init(*args, **kwargs)

    def plotOn(self, *args, **kwargs):
        # Redefinition of `RooDataHist.plotOn` for keyword arguments.
        # The keywords must correspond to the CmdArg of the `plotOn` function.
        args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
        return self._plotOn(*args, **kwargs)
