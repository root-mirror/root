import py, pytest, os, sys, math, warnings
from support import setup_make


setup_make("runvectorDict.so")

currpath = py.path.local(__file__).dirpath()
test_dct = str(currpath.join("runvectorDict.so"))

import cppyy
cppyy.load_reflection_info(test_dct)

all_configs = [('cppyy', 'cppyy.gbl')]

preamble = "@pytest.mark.benchmark(group=group, warmup=True)"


try:
    import __pypy__
    import py_runvector
    all_configs.append(('py', 'py_runvector'))    # too slow to run on CPython
except ImportError:
    try:
        import py11_runvector
        all_configs.append(('py11', 'py11_runvector'))
        py11 = True
    except ImportError:
        warnings.warn('pybind11 tests disabled')
        py11 = False

    try:
        import swig_runvector
        all_configs.append(('swig', 'swig_runvector.cvar'))
        swig = True
    except ImportError:
        warnings.warn('swig tests disabled')
        swig = False

all_benches = []


#- group: stl-vector ---------------------------------------------------------
all_benches.append(('stl-vector', (
"""
def test_{0}_stl_vector(benchmark):
    benchmark(sum, {1}.global_vector)
""",
)))


#- actual creation of all benches --------------------------------------------
for group, benches in all_benches:
    for bench in benches:
        for label, modname in all_configs:
            exec(preamble+bench.format(label, modname))
