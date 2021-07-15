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
\class RooCategory
\brief \parblock \endparblock
\htmlonly
<div class="pyrootbox">
\endhtmlonly

## PyROOT

Constructor of RooCategory takes a map as an argument also supports python dictionaries.
For example, the following code is equivalent in PyROOT:
\code{.py}
# Default bindings :
mixState = ROOT.RooCategory("mixState", "B0/B0bar mixing state")
mixState.defineType("mixed", -1)
mixState.defineType("unmixed", 1)

# With pythonization :
mixState = ROOT.RooCategory("mixState", "B0/B0bar mixing state", {"mixed" : -1, "unmixed" : 1})

\endcode

\htmlonly
</div>
\endhtmlonly
*/
"""

from ._utils import _dict_to_std_map


class RooCategory(object):
    def __init__(self, *args, **kwargs):
        # Redefinition of `RooCategory` constructor for converting python dict to std::map.
        # The instances in the dict must correspond to the template argument in std::map of the constructor.
        if len(args) > 2 and isinstance(args[2], dict):
            args = list(args)
            args[2] = _dict_to_std_map(args[2])

        self._init(*args, **kwargs)
