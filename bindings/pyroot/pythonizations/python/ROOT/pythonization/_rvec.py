# Author: Stefan Wunsch CERN  08/2018

################################################################################
# Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

r'''
/**
\class ROOT::VecOps::RVec
\brief \parblock \endparblock
\htmlonly
<div class="pyrootbox">
\endhtmlonly
## PyROOT

The ROOT::RVec class has additional features in Python, which allow to adopt memory
from Numpy arrays and vice versa. The purpose of these features is the copyless
interfacing of Python and C++ using their most common data containers, Numpy arrays
and RVec with a std::vector interface.

### Conversion of RVecs to Numpy arrays

RVecs of fundamental types (int, float, ...) have in Python the `__array_interface__`
attribute attached. This information allows Numpy to adopt the memory of RVecs without
copying the content. You can find further documentation regarding the Numpy array interface
[here](https://numpy.org/doc/stable/reference/arrays.interface.html). The following code example
demonstrates the memory adoption mechanism using `numpy.asarray`.

\code{.py}
rvec = ROOT.RVec('double')((1, 2, 3))
print(rvec) # { 1.0000000, 2.0000000, 3.0000000 }

npy = numpy.asarray(rvec)
print(npy) # [1. 2. 3.]

rvec[0] = 42
print(npy) # [42. 2. 3.]
\endcode

### Conversion of Numpy arrays to RVecs

Data owned by Numpy arrays with fundamental types (int, float, ...) can be adopted by RVecs. To
create an RVec from a Numpy array, ROOT offers the facility ROOT.VecOps.AsRVec, which performs
a similar operation to `numpy.asarray`, but vice versa. A code example demonstrating the feature and
the adoption of the data owned by the Numpy array is shown below.

\code{.py}
npy = numpy.array([1.0, 2.0, 3.0])
print(npy) # [1. 2. 3.]

rvec = ROOT.VecOps.AsRVec(npy)
print(rvec) # { 1.0000000, 2.0000000, 3.0000000 }

npy[0] = 42
print(rvec) # { 42.000000, 2.0000000, 3.0000000 }
\endcode

\htmlonly
</div>
\endhtmlonly
*/
'''

from ROOT import pythonization
from libROOTPythonizations import GetEndianess, GetDataPointer, GetSizeOfType


_array_interface_dtype_map = {
    "float": "f",
    "double": "f",
    "int": "i",
    "long": "i",
    "Long64_t": "i",
    "unsigned int": "u",
    "unsigned long": "u",
    "ULong64_t": "u",
}


def get_array_interface(self):
    cppname = type(self).__cpp_name__
    for dtype in _array_interface_dtype_map:
        if cppname.endswith("<{}>".format(dtype)):
            dtype_numpy = _array_interface_dtype_map[dtype]
            dtype_size = GetSizeOfType(dtype)
            endianess = GetEndianess()
            size = self.size()
            # Numpy breaks for data pointer of 0 even though the array is empty.
            # We set the pointer to 1 but the value itself is arbitrary and never accessed.
            if self.empty():
                pointer = 1
            else:
                pointer = GetDataPointer(self, cppname, "data")
            return {
                "shape": (size, ),
                "typestr": "{}{}{}".format(endianess, dtype_numpy, dtype_size),
                "version": 3,
                "data": (pointer, False)
            }


def add_array_interface_property(klass, name):
    if True in [
            name.endswith("<{}>".format(dtype)) for dtype in _array_interface_dtype_map
    ]:
        klass.__array_interface__ = property(get_array_interface)


@pythonization()
def pythonize_rvec(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class

    if name.startswith("ROOT::VecOps::RVec<"):
        # Add numpy array interface
        add_array_interface_property(klass, name)

    return True
