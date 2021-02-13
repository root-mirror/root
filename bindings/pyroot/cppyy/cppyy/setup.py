import codecs, glob, os, sys, re
from setuptools import setup, find_packages, Extension
from distutils import log

from setuptools.dist import Distribution
from setuptools.command.install import install as _install

force_bdist = False
if '--force-bdist' in sys.argv:
    force_bdist = True
    sys.argv.remove('--force-bdist')

add_pkg = ['cppyy']
try:
    import __pypy__, sys
    version = sys.pypy_version_info
    requirements = ['cppyy-backend']
    if version[0] == 5:
        if version[1] <= 9:
            requirements = ['cppyy-cling<6.12', 'cppyy-backend<0.3']
            add_pkg += ['cppyy_compat']
        elif version[1] <= 10:
            requirements = ['cppyy-cling<=6.15', 'cppyy-backend<0.4']
    elif version[0] == 6:
        if version[1] <= 0:
            requirements = ['cppyy-cling<=6.15', 'cppyy-backend<1.1']
    elif version[0] == 7:
        if version[1] <= 2:
            requirements = ['cppyy-cling<=6.18.2.3', 'cppyy-backend<=1.10']
        else:
            requirements = ['cppyy-cling<=6.18.2.7', 'cppyy-backend<=1.10']
except ImportError:
    # CPython
    requirements = ['cppyy-cling==6.18.2.7', 'cppyy-backend==1.10.8', 'CPyCppyy==1.10.2']

setup_requirements = ['wheel']
if 'build' in sys.argv or 'install' in sys.argv:
    setup_requirements += requirements

here = os.path.abspath(os.path.dirname(__file__))
with codecs.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# https://packaging.python.org/guides/single-sourcing-package-version/
def read(*parts):
    with codecs.open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


#
# platform-dependent helpers
#
def is_manylinux():
    try:
       for line in open('/etc/redhat-release').readlines():
           if 'CentOS release 6.10 (Final)' in line:
               return True
    except (OSError, IOError):
        pass
    return False


#
# customized commands
#
class my_install(_install):
    def run(self):
        # base install
        _install.run(self)

        # force build of the .pch underneath the cppyy package if not available yet
        install_path = os.path.join(os.getcwd(), self.install_libbase, 'cppyy')

        try:
            import cppyy_backend as cpb
            if not os.path.exists(os.path.join(cpb.__file__, 'etc', 'allDict.cxx.pch')):
                log.info("installing pre-compiled header in %s", install_path)
                cpb.loader.set_cling_compile_options(True)
                cpb.loader.ensure_precompiled_header(install_path, 'allDict.cxx.pch')
        except (ImportError, AttributeError):
            # ImportError may occur with wrong pip requirements resolution (unlikely)
            # AttributeError will occur with (older) PyPy as it relies on older backends
            pass

    def get_outputs(self):
        outputs = _install.get_outputs(self)
        # pre-emptively add allDict.cxx.pch, which may or may not be created; need full
        # path to make sure the final relative path is correct
        outputs.append(os.path.join(os.getcwd(), self.install_libbase, 'cppyy', 'allDict.cxx.pch'))
        return outputs


cmdclass = {
        'install': my_install }


#
# customized distribition to disable binaries
#
class MyDistribution(Distribution):
    def run_commands(self):
        # pip does not resolve dependencies before building binaries, so unless
        # packages are installed one-by-one, on old install is used or the build
        # will simply fail hard. The following is not completely quiet, but at
        # least a lot less conspicuous.
        if not is_manylinux() and not force_bdist:
            disabled = set((
                'bdist_wheel', 'bdist_egg', 'bdist_wininst', 'bdist_rpm'))
            for cmd in self.commands:
                if not cmd in disabled:
                    self.run_command(cmd)
                else:
                    log.info('Command "%s" is disabled', cmd)
                    cmd_obj = self.get_command_obj(cmd)
                    cmd_obj.get_outputs = lambda: None
        else:
            return Distribution.run_commands(self)


setup(
    name='cppyy',
    version=find_version('python', 'cppyy', '_version.py'),
    description='Cling-based Python-C++ bindings',
    long_description=long_description,

    url='http://cppyy.readthedocs.org',

    # Author details
    author='Wim Lavrijsen',
    author_email='WLavrijsen@lbl.gov',

    license='LBNL BSD',

    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Developers',

        'Topic :: Software Development',
        'Topic :: Software Development :: Interpreters',

        'License :: OSI Approved :: BSD License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: C',
        'Programming Language :: C++',

        'Natural Language :: English'
    ],

    setup_requires=setup_requirements,
    install_requires=requirements,

    keywords='C++ bindings data science calling language integration',

    package_dir={'': 'python'},
    packages=find_packages('python', include=add_pkg),

    cmdclass=cmdclass,
    distclass=MyDistribution,

    zip_safe=False,
)
