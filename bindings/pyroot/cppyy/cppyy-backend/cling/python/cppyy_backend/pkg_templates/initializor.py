"""
Support utilities for bindings.
"""
import glob
import json
from distutils.command.clean import clean
from distutils.util import get_platform
from setuptools.command.build_py import build_py
from wheel.bdist_wheel import bdist_wheel
import gettext
import inspect
import os
import re
import setuptools
import subprocess
import sys
try:
    #
    # Python2.
    #
    from imp import load_source
except ImportError:
    #
    # Python3.
    #
    import importlib.util

    def load_source(module_name, file_path, add_to_sys=False):
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        # Optional; only necessary if you want to be able to import the module
        # by name later.
        if add_to_sys:
            sys.modules[module_name] = module
        return module

import cppyy


gettext.install(__name__)

# Keep PyCharm happy.
_ = _


PRIMITIVE_TYPES = re.compile(r"\b(bool|char|short|int|unsigned|long|float|double)\b")


def add_pythonizations(py_files, noisy=False):
    for py_file in py_files:
        if noisy:
            print('check', py_file, 'for pythonizors')
        if not os.path.basename(py_file).startswith('pythonize'):
            continue
        module_name = inspect.getmodulename(py_file)
        module      = load_source(module_name, py_file)
        funcs       = inspect.getmembers(module, predicate=inspect.isroutine)

        pythonizors = {}
        for name, func in funcs:
            if not name.startswith('pythonize'):
                continue
            tokens = name.split('_')
            if len(tokens) > 1:
                namespace = tokens[1]
                if noisy:
                    print('added pythonization', func, namespace)
                if namespace == 'gbl':
                    cppyy.py.add_pythonization(func)
                else:
                    cppyy.py.add_pythonization(func, namespace)



def initialise(pkg, lib_file, map_file, noisy=False):
    """
    Initialise the bindings module.

    :param pkg:             The bindings package.
    :param __init__py:      Base __init__.py file of the bindings.
    :param cmake_shared_library_prefix:
                            ${cmake_shared_library_prefix}
    :param cmake_shared_library_suffix:
                            ${cmake_shared_library_suffix}
    """
    def add_to_pkg(file, keyword, simplenames, children):
        def map_operator_name(name):
            """
            Map the given C++ operator name on the python equivalent.
            """
            CPPYY__idiv__ = "__idiv__"
            CPPYY__div__  = "__div__"
            gC2POperatorMapping = {
                "[]": "__getitem__",
                "()": "__call__",
                "/": CPPYY__div__,
                "%": "__mod__",
                "**": "__pow__",
                "<<": "__lshift__",
                ">>": "__rshift__",
                "&": "__and__",
                "|": "__or__",
                "^": "__xor__",
                "~": "__inv__",
                "+=": "__iadd__",
                "-=": "__isub__",
                "*=": "__imul__",
                "/=": CPPYY__idiv__,
                "%=": "__imod__",
                "**=": "__ipow__",
                "<<=": "__ilshift__",
                ">>=": "__irshift__",
                "&=": "__iand__",
                "|=": "__ior__",
                "^=": "__ixor__",
                "==": "__eq__",
                "!=": "__ne__",
                ">": "__gt__",
                "<": "__lt__",
                ">=": "__ge__",
                "<=": "__le__",
            }

            op = name[8:]
            result = gC2POperatorMapping.get(op, None)
            if result:
                return result

            bTakesParams = 1
            if op == "*":
                # dereference v.s. multiplication of two instances
                return "__mul__" if bTakesParams else "__deref__"
            elif op == "+":
                # unary positive v.s. addition of two instances
                return "__add__" if bTakesParams else "__pos__"
            elif op == "-":
                # unary negative v.s. subtraction of two instances
                return "__sub__" if bTakesParams else "__neg__"
            elif op == "++":
                # prefix v.s. postfix increment
                return "__postinc__" if bTakesParams else "__preinc__"
            elif op == "--":
                # prefix v.s. postfix decrement
                return "__postdec__" if bTakesParams else "__predec__"
            # might get here, as not all operator methods are handled (new, delete, etc.)
            return name

        #
        # Add level 1 objects to the pkg namespace.
        #
        if len(simplenames) > 1:
            return
        #
        # Ignore some names based on heuristics.
        #
        simplename = simplenames[0]
        if simplename in ('void', 'sizeof', 'const'):
            return
        if simplename[0] in '0123456789':
            #
            # Don't attempt to look up numbers (i.e. non-type template parameters).
            #
            return
        if PRIMITIVE_TYPES.search(simplename):
            return
        if simplename.startswith("operator"):
            simplename = map_operator_name(simplename)
        #
        # Classes, variables etc.
        #
        try:
            entity = getattr(cppyy.gbl, simplename)
        except AttributeError as e:
            if noisy:
                print(_("Unable to lookup {}:{} cppyy.gbl.{} ({})").format(file, keyword, simplename, children))
            #raise
        else:
            if getattr(entity, "__module__", None) == "cppyy.gbl":
                setattr(entity, "__module__", pkg)
            setattr(pkg_module, simplename, entity)

    pkg_dir = os.path.dirname(__file__)
    if "." in pkg:
        pkg_namespace, pkg_simplename = pkg.rsplit(".", 1)
    else:
        pkg_namespace, pkg_simplename = "", pkg
    pkg_module = sys.modules[pkg]
    #
    # Load the library.
    #
    cppyy.load_reflection_info(os.path.join(pkg_dir, lib_file))

    #
    # Load pythonizations
    #
    try:
        pythonization_files = glob.glob(os.path.join(pkg_dir, '**/pythonize*.py'), recursive=True)
    except TypeError:
        # versions older than 3.5 do not support 'recursive'
        # TODO: below is good enough for most cases, but not recursive
        pythonization_files = glob.glob(os.path.join(pkg_dir, 'pythonize*.py'))
    add_pythonizations(pythonization_files, noisy=noisy)

    #
    # Parse the map file.
    #
    with open(os.path.join(pkg_dir, map_file), 'r') as map_file:
        files = json.load(map_file)

    #
    # Iterate over all the items at the top level of each file, and add them
    # to the pkg.
    #
    for file in files:
        for child in file["children"]:
            if not child["kind"] in ('class', 'var', 'namespace', 'typedef'):
                continue
            simplenames = child["name"].split('::')
            add_to_pkg(file["name"], child["kind"], simplenames, child)
