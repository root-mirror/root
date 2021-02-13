# Author: Danilo Piparo, Stefan Wunsch CERN  08/2018

################################################################################
# Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

r'''
/**
\class TDirectory
\brief \parblock \endparblock
\htmlonly
<div class="pyrootbox">
\endhtmlonly
## PyROOT

From Python, it is possible to inspect the content of a TDirectory object
as if the subdirectories and objects it contains were its attributes.
Moreover, once a subdirectory or object is accessed for the first time,
it is cached for later use.
For example, assuming `d` is a TDirectory instance:
\code{.py}
# Access a subdirectory
d.subdir

# We can go further down in the hierarchy of directories
d.subdir.subsubdir

# Access an object (e.g. a histogram) in the directory
d.obj

# ... or in a subdirectory
d.subdir.obj

# Wrong attribute: raises AttributeError
d.wrongAttr
\endcode

Furthermore, TDirectory implements a `WriteObject` Python method which relies
on TDirectory::WriteObjectAny. This method is a no-op for TDirectory objects,
but it is useful for objects of TDirectory subclasses such as TDirectoryFile
and TFile, which inherit it. Please refer to the documentation of those classes
for more information.
\htmlonly
</div>
\endhtmlonly
*/
'''

from libROOTPythonizations import AddDirectoryGetAttrPyz, AddDirectoryWritePyz
from ROOT import pythonization
import cppyy

# This pythonization must be set as not lazy, otherwise the mechanism cppyy uses
# to pythonize classes will not be able to be triggered on this very core class.
# The pythonization does not have arguments since it is not fired by cppyy but
# manually upon import of the ROOT module.
@pythonization(lazy = False)
def pythonize_tdirectory():
    klass = cppyy.gbl.TDirectory
    AddDirectoryGetAttrPyz(klass)
    AddDirectoryWritePyz(klass)
    return True
