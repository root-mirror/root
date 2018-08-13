from ROOT import pythonization
from libROOTPython import GetEndianess, GetVectorDataPointer, GetSizeOfType

_array_interface_dtypes = [
    "float", "double", "int", "long", "unsigned int", "unsigned long"
]

_array_interface_dtype_map = {
    "float": "f",
    "double": "f",
    "int": "i",
    "long": "i",
    "unsigned int": "u",
    "unsigned long": "u"
}


def get_array_interface(self):
    cppname = type(self).__cppname__
    for dtype in _array_interface_dtypes:
        if cppname.endswith("<{}>".format(dtype)):
            dtype_numpy = _array_interface_dtype_map[dtype]
            dtype_size = GetSizeOfType(dtype)
            endianess = GetEndianess()
            size = self.size()
            pointer = GetVectorDataPointer(self, cppname)
            return {
                "shape": (size, ),
                "typestr": "{}{}{}".format(endianess, dtype_numpy, dtype_size),
                "version": 3,
                "data": (pointer, False)
            }


def add_array_interface_property(klass, name):
    if True in [
            "<{}>".format(dtype) in name for dtype in _array_interface_dtypes
    ]:
        klass.__array_interface__ = property(get_array_interface)


@pythonization
def pythonize_rvec(klass, name):
    # Parameters:
    # klass: class to be pythonized
    # name: string containing the name of the class

    # Add numpy array interface
    add_array_interface_property(klass, name)

    return True
