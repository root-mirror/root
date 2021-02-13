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

    def load_source(module_name, file_path):
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        # Optional; only necessary if you want to be able to import the module
        # by name later.
        sys.modules[module_name] = module
        return module

import cppyy


gettext.install(__name__)

# Keep PyCharm happy.
_ = _


PRIMITIVE_TYPES = re.compile(r"\b(bool|char|short|int|unsigned|long|float|double)\b")


def initialise(pkg, __init__py, cmake_shared_library_prefix, cmake_shared_library_suffix):
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
            print(children)

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
            print(_("Unable to lookup {}:{} cppyy.gbl.{} ({})").format(file, keyword, simplename, children))
            #raise
        else:
            if getattr(entity, "__module__", None) == "cppyy.gbl":
                setattr(entity, "__module__", pkg)
            setattr(pkg_module, simplename, entity)

    pkg_dir = os.path.dirname(__init__py)
    if "." in pkg:
        pkg_namespace, pkg_simplename = pkg.rsplit(".", 1)
    else:
        pkg_namespace, pkg_simplename = "", pkg
    pkg_module = sys.modules[pkg]
    lib_name = pkg_namespace + pkg_simplename + "Cppyy"
    lib_file = cmake_shared_library_prefix + lib_name + cmake_shared_library_suffix
    map_file = os.path.join(pkg_dir, pkg_simplename + ".map")
    #
    # Load the library.
    #
    cppyy.load_reflection_info(os.path.join(pkg_dir, lib_file))
    #
    # Parse the map file.
    #
    with open(map_file, 'rU') as map_file:
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
    #
    # Load any customisations.
    #
    extra_pythons = glob.glob(os.path.join(pkg_dir, "extra_*.py"))
    for extra_python in extra_pythons:
        extra_module = os.path.basename(extra_python)
        extra_module = pkg + "." + os.path.splitext(extra_module)[0]
        #
        # Deleting the modules after use runs the risk of GC running on
        # stuff we are using, such as ctypes.c_int.
        #
        extra = load_source(extra_module, extra_python)
        #
        # Valid customisations are routines named "c13n_<something>".
        #
        fns = inspect.getmembers(extra, predicate=inspect.isroutine)
        fns = {fn[0]: fn[1] for fn in fns if fn[0].startswith("c13n_")}
        for fn in sorted(fns):
            fns[fn](pkg_module)


def setup(pkg, setup_py, cmake_shared_library_prefix, cmake_shared_library_suffix, extra_pythons,
          pkg_version, author, author_email, url, license):
    """
    Wrap setuptools.setup for some bindings.

    :param pkg:             Name of the bindings.
    :param setup_py:        Base __init__.py file of the bindings.
    :param cmake_shared_library_prefix:
                            ${cmake_shared_library_prefix}
    :param cmake_shared_library_suffix:
                            ${cmake_shared_library_suffix}
    :param extra_pythons:   Semicolon-separated list of customisation code.
    :param pkg_version:     The version of the bindings.
    :param author:          The name of the library author.
    :param author_email:    The email address of the library author.
    :param url:             The home page for the library.
    :param license:         The license.
    """
    pkg_dir = os.path.dirname(setup_py)
    if "." in pkg:
        pkg_namespace, pkg_simplename = pkg.rsplit(".", 1)
    else:
        pkg_namespace, pkg_simplename = "", pkg
    lib_name = pkg_namespace + pkg_simplename + "Cppyy"
    lib_file = cmake_shared_library_prefix + lib_name + cmake_shared_library_suffix
    long_description = """Bindings for {}.
These bindings are based on https://cppyy.readthedocs.io/en/latest/, and can be
used as per the documentation provided via the cppyy.cgl namespace. The environment
variable LD_LIBRARY_PATH must contain the path of the {}.rootmap file. Use
"import cppyy; from cppyy.gbl import <some-C++-entity>".

Alternatively, use "import {}". This convenience wrapper supports "discovery" of the
available C++ entities using, for example Python 3's command line completion support.
""".replace("{}", pkg)

    class my_build_py(build_py):
        def run(self):
            #
            # Base build.
            #
            build_py.run(self)
            #
            # Custom build.
            #
            #
            # Move CMake output to self.build_lib.
            #
            pkg_subdir = pkg.replace(".", os.path.sep)
            if pkg_namespace:
                #
                # Implement a pkgutil-style namespace package as per the guidance on
                # https://packaging.python.org/guides/packaging-namespace-packages.
                #
                namespace_init = os.path.join(pkg_namespace, "__init__.py")
                with open(namespace_init, "w") as f:
                    f.write("__path__ = __import__('pkgutil').extend_path(__path__, __name__)\n")
                self.copy_file(namespace_init, os.path.join(self.build_lib, namespace_init))
            for f in self.package_data[pkg]:
                self.copy_file(os.path.join(pkg_dir, pkg_subdir, f), os.path.join(self.build_lib, pkg_subdir, f))

    class my_clean(clean):
        def run(self):
            #
            # Custom clean.
            # TODO: There is no way to reliably clean the "dist" directory.
            #
            #
            #  Base clean.
            #
            clean.run(self)

    class my_bdist_wheel(bdist_wheel):
        def finalize_options(self):
            #
            # This is a universal (Python2/Python3), but platform-specific (has
            # compiled parts) package; a combination that wheel does not recognize,
            # thus simply fool it.
            #
            self.plat_name = get_platform()
            bdist_wheel.finalize_options(self)
            self.root_is_pure = True

    package_data = [lib_file, pkg_simplename + '.rootmap', pkg_simplename + '_rdict.pcm', pkg_simplename + ".map"]
    package_data += [ep for ep in extra_pythons.split(";") if ep]
    setuptools.setup(
        name=pkg,
        version=pkg_version,
        author=author,
        author_email=author_email,
        url=url,
        license=license,
        description='Bindings for ' + pkg,
        long_description=long_description,
        platforms=['any'],
        package_data={pkg: package_data},
        packages=[pkg],
        zip_safe=False,
        cmdclass={
            'build_py': my_build_py,
            'clean': my_clean,
            'bdist_wheel': my_bdist_wheel,
        },
    )


def find_pips():
    """
    What pip versions do we have?

    :return: [pip_program]
    """
    possible_pips = ['pip', 'pip2', 'pip3']
    pips = {}
    for pip in possible_pips:
        try:
            #
            # The command 'pip -V' returns a string of the form:
            #
            #   pip 9.0.1 from /usr/lib/python2.7/dist-packages (python 2.7)
            #
            version = subprocess.check_output([pip, '-V'])
        except subprocess.CalledProcessError:
            pass
        else:
            version = version.rsplit('(', 1)[-1]
            version = version.split()[-1]
            #
            # All pip variants that map onto a given Python version are de-duped.
            #
            pips[version] = pip
    #
    # We want the pip names.
    #
    assert len(pips), 'No viable pip versions found'
    return pips.values()
