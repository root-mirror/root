# Author: Stefan Wunsch CERN  07/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from ROOT import pythonization
from libROOTPythonizations import GetEndianess, GetDataPointer, GetSizeOfType
from ROOT.pythonization._rvec import _array_interface_dtype_map
import cppyy


def get_array_interface(self):
    """
    Return the array interface dictionary

    Parameters:
        self: RTensor object
    Returns:
        Dictionary following the Numpy array interface specifications
    """
    cppname = type(self).__cpp_name__
    idx1 = cppname.find("RTensor<")
    idx2 = cppname.find(",", idx1)
    dtype = cppname[idx1 + 8:idx2]
    dtype_numpy = _array_interface_dtype_map[dtype]
    dtype_size = GetSizeOfType(dtype)
    endianess = GetEndianess()
    shape = self.GetShape()
    strides = self.GetStrides()
    # Numpy breaks for data pointer of 0 even though the array is empty.
    # We set the pointer to 1 but the value itself is arbitrary and never accessed.
    pointer = GetDataPointer(self, cppname, "GetData")
    if pointer == 0:
        pointer == 1
    return {
        "shape": tuple(s for s in shape),
        "strides": tuple(s * dtype_size for s in strides),
        "typestr": "{}{}{}".format(endianess, dtype_numpy, dtype_size),
        "version": 3,
        "data": (pointer, False)
    }


def add_array_interface_property(klass, name):
    """
    Attach the array interface as property if the data-type of the RTensor
    elements is one of the supported basic types

    Parameters:
        klass: class to be pythonized
        name: string containing the name of the class
    """
    if True in [
            name.find("RTensor<{},".format(dtype)) != -1 for dtype in _array_interface_dtype_map
    ]:
        klass.__array_interface__ = property(get_array_interface)


def RTensorGetitem(self, idx):
    """
    Implementation of the __getitem__ special function for RTensor

    Parameters:
        self: RTensor object
        idx: Indices passed to RTensor[indices] operator
    Returns:
        New RTensor object if indices represent a slice or the requested element
    """
    # Make single index iterable and convert to list
    if not hasattr(idx, "__len__"):
        idx = [idx]
    idx = list(idx)

    # Check shape
    shape = self.GetShape()
    if shape.size() != len(idx):
        raise Exception("RTensor with rank {} got {} indices.".format(shape.size(), len(idx)))

    # Convert negative indices and Nones
    isSlice = False
    for i, x in enumerate(idx):
        if type(x) == slice:
            isSlice = True
            start = 0 if x.start is None else x.start
            stop = shape[i] if x.stop is None else x.stop
            if stop < 0:
                stop += shape[i]
            if x.step != None:
                raise Exception("RTensor does not support slices with step size unequal 1.")
            idx[i] = slice(start, stop, None)
        else:
            if x < 0:
                idx[i] += shape[i]

    # If a slice is requested, return a new RTensor
    if isSlice:
        idxVec = cppyy.gbl.std.vector("vector<size_t>")(len(idx))
        for i, x in enumerate(idx):
            idxVec[i].resize(2)
            if type(x) == slice:
                idxVec[i][0] = x.start
                idxVec[i][1] = x.stop
            else:
                idxVec[i][0] = x
                idxVec[i][1] = x + 1
        return self.Slice(idxVec)

    # Otherwise, access element by array of indices
    idxVec = cppyy.gbl.std.vector("size_t")(len(idx))
    for i, x in enumerate(idx):
        idxVec[i] = x
    return self(idxVec)


@pythonization()
def pythonize_rtensor(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class

    if name.startswith("TMVA::Experimental::RTensor<"):
        # Add numpy array interface
        add_array_interface_property(klass, name)
        # Get elements, including slices
        klass.__getitem__ = RTensorGetitem

    return True
