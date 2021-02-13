from cppyy.interactive import *

# namespace at the global level
assert std

# cppyy functions
assert cppdef
assert include

try:
    import __pypy__
  # 'cppyy.gbl' bound to 'g'
    assert g
    assert g.std
except ImportError:
 # full lazy lookup available
    assert gInterpreter
