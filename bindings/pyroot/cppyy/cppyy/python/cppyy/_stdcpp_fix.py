import sys

# It may be that the interpreter (wether python or pypy-c) was not linked
# with C++; force its loading before doing anything else (note that not
# linking with C++ spells trouble anyway for any C++ libraries ...)
if 'linux' in sys.platform and 'GCC' in sys.version:
    # TODO: check executable to see whether linking indeed didn't happen
    import ctypes
    try:
        stdcpp = ctypes.CDLL('libstdc++.so', ctypes.RTLD_GLOBAL)
    except Exception:
        pass
# TODO: what if Linux/clang and what if Mac?
