import codecs, glob, os, sys, subprocess
from setuptools import setup, find_packages, Extension
from distutils import log

from setuptools.dist import Distribution
from setuptools.command.install import install as _install
from distutils.command.build_ext import build_ext as _build_ext
from distutils.command.clean import clean as _clean
from distutils.dir_util import remove_tree
try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
    has_wheel = True
except ImportError:
    has_wheel = False
from distutils.errors import DistutilsSetupError

force_bdist = False
if '--force-bdist' in sys.argv:
    force_bdist = True
    sys.argv.remove('--force-bdist')

requirements = ['cppyy-cling<6.18.3', 'cppyy-cling>=6.18.2.4']
setup_requirements = ['wheel']
if 'build' in sys.argv or 'install' in sys.argv:
    setup_requirements += requirements

here = os.path.abspath(os.path.dirname(__file__))
with codecs.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


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

def _get_linker_options():
    if 'win32' in sys.platform:
        link_libraries = ['libCore', 'libThread', 'libRIO', 'libCling']
        import cppyy_backend
        link_dirs = [os.path.join(os.path.dirname(cppyy_backend.__file__), 'lib')]
    else:
        link_libraries = None
        link_dirs = None
    return link_libraries, link_dirs

def _get_config_exec():
    return [sys.executable, '-m', 'cppyy_backend._cling_config']

def get_include_path():
    config_exec_args = _get_config_exec()
    config_exec_args.append('--incdir')
    cli_arg = subprocess.check_output(config_exec_args)
    return cli_arg.decode("utf-8").strip()

def get_cflags():
    config_exec_args = _get_config_exec()
    config_exec_args.append('--auxcflags')
    cli_arg = subprocess.check_output(config_exec_args)
    return cli_arg.decode("utf-8").strip()


#
# customized commands
#
class my_build_cpplib(_build_ext):
    def build_extension(self, ext):
        include_dirs = ext.include_dirs + [get_include_path()]
        log.info('checking for %s', self.build_temp)
        if not os.path.exists(self.build_temp):
            log.info('creating %s', self.build_temp)
            os.makedirs(self.build_temp)
        extra_postargs = ['-O2']+get_cflags().split()
        if 'win32' in sys.platform:
        # /EHsc and sometimes /MT are hardwired in distutils, but the compiler/linker will
        # let the last argument take precedence
            extra_postargs += ['/GR', '/EHsc-', '/MD']
        objects = self.compiler.compile(
            ext.sources,
            output_dir=self.build_temp,
            include_dirs=include_dirs,
            debug=self.debug,
            extra_postargs=extra_postargs)

        ext_path = self.get_ext_fullpath(ext.name)
        output_dir = os.path.dirname(ext_path)
        libname_base = 'libcppyy_backend'
        libname = libname_base+self.compiler.shared_lib_extension
        extra_postargs = list()
        if 'linux' in sys.platform:
            extra_postargs.append('-Wl,-Bsymbolic-functions')
        elif 'win32' in sys.platform:
            # force the export results in the proper directory.
            extra_postargs.append('/IMPLIB:'+os.path.join(output_dir, libname_base+'.lib'))
            import platform
            if '64' in platform.architecture()[0]:
                extra_postargs.append('/EXPORT:?__type_info_root_node@@3U__type_info_node@@A')

        log.info("now building %s", libname)
        link_libraries, link_dirs = _get_linker_options()
        self.compiler.link_shared_object(
            objects, libname,
            libraries=link_libraries, library_dirs=link_dirs,
            # export_symbols=[], # ie. all (hum, that puts the output in the wrong directory)
            build_temp=self.build_temp,
            output_dir=output_dir,
            debug=self.debug,
            target_lang='c++',
            extra_postargs=extra_postargs)

class my_clean(_clean):
    def run(self):
        # Custom clean. Clean everything except that which the base clean
        # (see below) or create_src_directory.py is responsible for.
        topdir = os.getcwd()
        if self.all:
            # remove build directories
            for directory in (os.path.join(topdir, "dist"),
                              os.path.join(topdir, "python", "cppyy_backend.egg-info")):
                if os.path.exists(directory):
                    remove_tree(directory, dry_run=self.dry_run)
                else:
                    log.warn("'%s' does not exist -- can't clean it",
                             directory)
        # Base clean.
        _clean.run(self)

class my_install(_install):
    def _get_install_path(self):
        # depending on goal, copy over pre-installed tree
        if hasattr(self, 'bdist_dir') and self.bdist_dir:
            install_path = self.bdist_dir
        else:
            install_path = self.install_lib
        return install_path

    def run(self):
        # base install
        _install.run(self)

        # custom install of backend
        log.info('Now installing cppyy_backend')
        builddir = self.build_lib
        if not os.path.exists(builddir):
            raise DistutilsSetupError('Failed to find build dir!')

        install_path = self._get_install_path()
        log.info('Copying installation to: %s ...', install_path)
        self.copy_tree(builddir, install_path)

        log.info('Install finished')

    def get_outputs(self):
        outputs = _install.get_outputs(self)
        #outputs.append(os.path.join(self._get_install_path(), 'cppyy_backend'))
        return outputs


cmdclass = {
        'build_ext': my_build_cpplib,
        'clean': my_clean,
        'install': my_install }
if has_wheel:
    class my_bdist_wheel(_bdist_wheel):
        def finalize_options(self):
         # this is a universal, but platform-specific package; a combination
         # that wheel does not recognize, thus simply fool it
            from distutils.util import get_platform
            self.plat_name = get_platform()
            self.universal = True
            _bdist_wheel.finalize_options(self)
            self.root_is_pure = True
    cmdclass['bdist_wheel'] = my_bdist_wheel


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
    name='cppyy-backend',
    description='C/C++ wrapper for Cling',
    long_description=long_description,
    url='http://pypy.org',

    # Author details
    author='PyPy Developers',
    author_email='pypy-dev@python.org',

    version='1.10.8',

    license='LBNL BSD',

    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',

        'Topic :: Software Development',
        'Topic :: Software Development :: Interpreters',

        'License :: OSI Approved :: BSD License',

        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',

        'Programming Language :: C',
        'Programming Language :: C++',

        'Natural Language :: English'
    ],

    keywords='C++ bindings data science',

    setup_requires=setup_requirements,
    install_requires=requirements,

    ext_modules=[Extension(os.path.join('cppyy_backend', 'lib', 'libcppyy_backend'),
        sources=glob.glob(os.path.join('src', 'clingwrapper.cxx')))],

    cmdclass=cmdclass,
    distclass=MyDistribution,

    zip_safe=False,
)
