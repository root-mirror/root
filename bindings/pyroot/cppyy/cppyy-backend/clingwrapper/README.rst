cppyy-backend
=============

C/C++ wrapper around Cling, the LLVM-based interactive C++ interpreter, for
use by cppyy, providing stable C and C++ Reflection APIs.

The compilation of cppyy-backend is very fast, but it will pull in
cppyy-cling, which takes a long time to install if there is no matching wheel
for your platform, forcing a build from source. By default, all cores will be
used, but it is also recommended to add the verbose flag to see progress:

  $ python -m pip install --verbose cppyy-backend

For further details, see cppyy's installation instructions:
  https://cppyy.readthedocs.io/en/latest/installation.html


Cling documentation is here:
  https://root.cern.ch/cling

----

Find the cppyy documentation here:
  http://cppyy.readthedocs.io

Change log:
  https://cppyy.readthedocs.io/en/latest/changelog.html

Bug reports/feedback:
  https://bitbucket.org/wlav/cppyy/issues?status=new&status=open
