import py, pytest, os, sys, math, warnings
from support import setup_make


setup_make("functioncallsDict.so")

currpath = py.path.local(__file__).dirpath()
test_dct = str(currpath.join("functioncallsDict.so"))

import cppyy
cppyy.load_reflection_info(test_dct)

import py_functioncalls

all_configs = [('py', 'py_functioncalls'), ('cppyy', 'cppyy.gbl')]

N = 10000
preamble = "@pytest.mark.benchmark(group=group, warmup=True)"
looprange = range
if sys.hexversion < 0x3000000:
    looprange = xrange

try:
    import __pypy__
except ImportError:
    try:
        import py11_functioncalls
        all_configs.append(('py11', 'py11_functioncalls'))
        py11 = True
    except ImportError:
        warnings.warn('pybind11 tests disabled')
        py11 = False

    try:
        import swig_functioncalls
        all_configs.append(('swig', 'swig_functioncalls'))
        swig = True
    except ImportError:
        warnings.warn('swig tests disabled')
        swig = False

all_benches = []


#- group: empty-free ---------------------------------------------------------
all_benches.append(('empty-free', (
"""
def test_{0}_free_empty_call(benchmark):
    benchmark({1}.empty_call)
""",
)))

#- group: empty-inst ---------------------------------------------------------
def call_instance_empty(inst):
    for i in looprange(N):
        inst.empty_call()

all_benches.append(('empty-inst', (
"""
def test_{0}_inst_empty_call(benchmark):
    inst = {1}.EmptyCall()
    benchmark(call_instance_empty, inst)
""",
)))


#- group: builtin-args-free --------------------------------------------------
all_benches.append(('builtin-args-free', (
"""
def test_{0}_free_take_an_int(benchmark):
    benchmark({1}.take_an_int, 1)
""",
"""
def test_{0}_free_take_a_double(benchmark):
    benchmark({1}.take_a_double, 1.)
""",
"""
def test_{0}_free_take_a_struct(benchmark):
    benchmark({1}.take_a_struct, {1}.Value())
""",
)))

#- group: builtin-args-inst --------------------------------------------------
def call_instance_take_an_int(inst, val):
    for i in looprange(N):
        inst.take_an_int(1)

def call_instance_take_a_double(inst, val):
    for i in looprange(N):
        inst.take_a_double(val)

def call_instance_take_a_struct(inst, val):
    for i in looprange(N):
        inst.take_a_struct(val)


all_benches.append(('builtin-args-inst', (
"""
def test_{0}_inst_take_an_int(benchmark):
    inst = {1}.TakeAValue()
    benchmark(call_instance_take_an_int, inst, 1)
""",
"""
def test_{0}_inst_take_a_double(benchmark):
    inst = {1}.TakeAValue()
    benchmark(call_instance_take_a_double, inst, 1.)
""",
"""
def test_{0}_inst_take_a_struct(benchmark):
    inst = {1}.TakeAValue()
    benchmark(call_instance_take_a_struct, inst, {1}.Value())
""",
)))

#- group: builtin-args-pass --------------------------------------------------
def call_instance_pass_int(inst, val):
    for i in looprange(N):
        inst.pass_int(val)

all_benches.append(('builtin-args-pass', (
"""
def test_{0}_inst_pass_int(benchmark):
    inst = {1}.TakeAValue()
    benchmark(call_instance_pass_int, inst, 1)
""",
)))


#- group: do_work-free -------------------------------------------------------
all_benches.append(('do_work-free', (
"""
def test_{0}_free_do_work(benchmark):
    benchmark({1}.do_work, 1.)
""",
)))

#- group: do_work-inst -------------------------------------------------------
def call_instance_do_work(inst):
    for i in looprange(N):
        inst.do_work(1.)

all_benches.append(('do_work-inst', (
"""
def test_{0}_inst_do_work(benchmark):
    inst = {1}.DoWork()
    benchmark(call_instance_do_work, inst)
""",
)))


#- group: overload-inst ------------------------------------------------------
def call_instance_overload(inst):
    for i in looprange(N):
        inst.add_it(1.)

all_benches.append(('overload-inst', (
"""
def test_{0}_inst_overload(benchmark):
    inst = {1}.OverloadedCall()
    benchmark(call_instance_overload, inst)
""",
)))


#- actual creation of all benches --------------------------------------------
for group, benches in all_benches:
    for bench in benches:
        for label, modname in all_configs:
            exec(preamble+bench.format(label, modname))
