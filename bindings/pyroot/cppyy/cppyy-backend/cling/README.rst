cppyy-cling
===========

A repackaging of Cling, the LLVM-based interactive C++ interpreter, as a
library for use as the backend to cppyy. This version of Cling is patched for
improved performance and better use with Python.

Wheels are available for the major platforms, but if you have to build from
source, building of LLVM will take a long time. By default, all cores will be
used, but it is also recommended to add the verbose flag to see progress:

  $ python -m pip install --verbose cppyy-cling

For further details, see cppyy's installation instructions:
  https://cppyy.readthedocs.io/en/latest/installation.html


Cling documentation is here:
  https://root.cern.ch/cling

----

Full cppyy documentation is here:
  http://cppyy.readthedocs.io/

Change log:
  https://cppyy.readthedocs.io/en/latest/changelog.html

Bug reports/feedback:
  https://bitbucket.org/wlav/cppyy/issues?status=new&status=open
