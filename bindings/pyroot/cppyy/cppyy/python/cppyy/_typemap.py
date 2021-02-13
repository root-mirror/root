""" Externally provided types: get looked up if all else fails, e.g.
    for typedef-ed C++ builtin types.
"""

import sys

def _create_mapper(cls, extra_dct=None):
    def mapper(name, scope):
        if scope:
            cppname = scope+'::'+name
            modname = 'cppyy.gbl.'+scope
        else:
            cppname = name
            modname = 'cppyy.gbl'
        dct = {'__cpp_name__' : cppname, '__module__' : modname}
        if extra_dct: dct.update(extra_dct)
        return type(name, (cls,), dct)
    return mapper

# from six.py ---
# Copyright (c) 2010-2017 Benjamin Peterson
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

def with_metaclass(meta, *bases):
    """Create a base class with a metaclass."""
    # This requires a bit of explanation: the basic idea is to make a dummy
    # metaclass for one level of class instantiation that replaces itself with
    # the actual metaclass.
    class metaclass(type):

        def __new__(cls, name, this_bases, d):
            return meta(name, bases, d)

        @classmethod
        def __prepare__(cls, name, this_bases):
            return meta.__prepare__(name, bases)
    return type.__new__(metaclass, 'temporary_class', (), {})
# --- end from six.py

class _BoolMeta(type):
    def __call__(self, val = bool()):
        if val: return True
        else: return False

class _Bool(with_metaclass(_BoolMeta, object)):
    pass


def initialize(backend):
    if not hasattr(backend, 'type_map'):
        return

    tm = backend.type_map

    # boolean type (builtin type bool can nog be subclassed)
    tm['bool'] = _create_mapper(_Bool)

    # char types
    str_tm = _create_mapper(str)
    for tp in ['char', 'unsigned char', 'signed char']:
        tm[tp] = str_tm
    if sys.hexversion < 0x3000000:
        tm['wchar_t'] = _create_mapper(unicode)
    else:
        tm['wchar_t'] = str_tm

    # integer types
    int_tm = _create_mapper(int)
    for tp in ['short', 'unsigned short', 'int']:
        tm[tp] = int_tm

    if sys.hexversion < 0x3000000:
        long_tm = _create_mapper(long)
    else:
        long_tm = tm['int']
    for tp in ['unsigned int', 'long', 'unsigned long', 'long long', 'unsigned long long']:
        tm[tp] = long_tm

    # floating point types
    float_tm = _create_mapper(float)
    for tp in ['float', 'double', 'long double']:
        tm[tp] = float_tm

    # void*
    import ctypes
    def voidp_init(self, arg=0):
        import cppyy, ctypes
        if arg == cppyy.nullptr: arg = 0
        ctypes.c_void_p.__init__(self, arg)
    tm['void*'] = _create_mapper(ctypes.c_void_p, {'__init__' : voidp_init})
