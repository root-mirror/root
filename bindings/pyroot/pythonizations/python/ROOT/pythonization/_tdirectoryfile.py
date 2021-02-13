# Author: Danilo Piparo, Stefan Wunsch, Massimiliano Galli CERN  08/2018

################################################################################
# Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

r'''
/**
\class TDirectoryFile
\brief \parblock \endparblock
\htmlonly
<div class="pyrootbox">
\endhtmlonly
## PyROOT

In the same way as for TDirectory, it is possible to inspect the content of a
TDirectoryFile object from Python as if the subdirectories and objects it
contains were its attributes. For more information, please refer to the
TDirectory documentation.

In addition to the attribute syntax, one can inspect a TDirectoryFile in Python
via the `Get` method. In this case, the subdirectory/object name is specified
as a string:
\code{.py}
# Access a subdirectory
d.Get('subdir')

# We can go further down in the hierarchy of directories
d.Get('subdir/subsubdir')

# Access an object (e.g. a histogram) in the directory
d.Get('obj')

# ... or in a subdirectory
d.Get('subdir/obj')

# Wrong attribute: returns null
x = d.Get('wrongAttr')  # x points to null
\endcode

Furthermore, TDirectoryFile inherits a `WriteObject` Python method from
TDirectory. Such method allows to write an object into a TDirectoryFile
with the following syntax:
\code{.py}
# Write object obj with identifier 'keyName'
d.WriteObject(obj, 'keyName')
\endcode
\htmlonly
</div>
\endhtmlonly
*/
'''

from libROOTPythonizations import AddTDirectoryFileGetPyz
from ROOT import pythonization

# Pythonizor function
@pythonization()
def pythonize_tdirectoryfile(klass, name):
    """
    TDirectoryFile inherits from TDirectory the pythonized attr syntax (__getattr__)
    and WriteObject method.
    On the other side, the Get() method is pythonised only in TDirectoryFile.
    Thus, the situation is now the following:

    1) __getattr__ : TDirectory --> TDirectoryFile --> TFile
        1.1) caches the returned object for future attempts
        1.2) raises AttributeError if object not found

    2) Get() : TDirectoryFile --> TFile
        2.1) does not cache the returned object
        2.2 returns nullptr if object not found
    """

    if name == 'TDirectoryFile':
        AddTDirectoryFileGetPyz(klass)

    return True
