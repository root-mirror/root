import codecs, glob, os, sys, subprocess
from setuptools import setup, find_packages, Extension
from distutils import log

from setuptools.dist import Distribution
from distutils.command.build_ext import build_ext as _build_ext
try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
    has_wheel = True
except ImportError:
    has_wheel = False

force_bdist = False
if '--force-bdist' in sys.argv:
    force_bdist = True
    sys.argv.remove('--force-bdist')

requirements = ['cppyy-cling==6.18.2.*', 'cppyy-backend==1.10.*']
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

def _get_link_libraries():
    if 'win32' in sys.platform:
        return ['libcppyy_backend', 'libCore']
    return []

def _get_link_dirs():
    if 'win32' in sys.platform:
        try:
            import cppyy_backend
            return [os.path.join(os.path.dirname(cppyy_backend.__file__), 'lib')]
        except ImportError:       # happens during egg_info and other non-build/install commands
            pass
    return []

def _get_config_exec():
    return [sys.executable, '-m', 'cppyy_backend._cling_config']

def get_cflags():
    config_exec_args = _get_config_exec()
    config_exec_args.append('--auxcflags')
    cli_arg = subprocess.check_output(config_exec_args)
    return cli_arg.decode("utf-8").strip()


#
# customized commands
#
class my_build_extension(_build_ext):
    def build_extension(self, ext):
        ext.extra_compile_args += ['-O2']+get_cflags().split()
        if ('linux' in sys.platform) or ('darwin' in sys.platform):
            if 'clang' in self.compiler.compiler_cxx[0]:
                ext.extra_compile_args += \
                   ['-Wno-bad-function-cast']    # clang for same
            elif 'g++' in self.compiler.compiler_cxx[0]:
                ext.extra_compile_args += \
                   ['-Wno-cast-function-type',   # g++ >8.2, complaint of CPyFunction cast
                    '-Wno-unknown-warning']      # since clang/g++ don't have the same options
            ext.extra_compile_args += \
                ['-Wno-register',                # C++17, Python headers
                 '-Wno-strict-aliasing']         # not all Pythons provide this
        if 'linux' in sys.platform:
            ext.extra_link_args += ['-Wl,-Bsymbolic-functions']
        elif 'win32' in sys.platform:
        # /EHsc and sometimes /MT are hardwired in distutils, but the compiler/linker will
        # let the last argument take precedence
            ext.extra_compile_args += ['/GR', '/EHsc-', '/MD']
            ext.extra_link_args += ['/EXPORT:_Init_thread_abort', '/EXPORT:_Init_thread_epoch',
                '/EXPORT:_Init_thread_footer', '/EXPORT:_Init_thread_header', '/EXPORT:_tls_index']
        return _build_ext.build_extension(self, ext)


cmdclass = {
    'build_ext': my_build_extension }


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
    name='CPyCppyy',
    version='1.10.2',
    description='Cling-based Python-C++ bindings for CPython',
    long_description=long_description,

    url='http://cppyy.readthedocs.io/',

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
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: C',
        'Programming Language :: C++',

        'Natural Language :: English'
    ],

    setup_requires=setup_requirements,
    install_requires=requirements,

    keywords='C++ bindings data science',

    ext_modules=[Extension('libcppyy',
        sources=glob.glob(os.path.join('src', '*.cxx')),
        include_dirs=['include'],
        libraries=_get_link_libraries(),
        library_dirs=_get_link_dirs())],

    headers=glob.glob(os.path.join('include', 'CPyCppyy', '*.h')),

    cmdclass=cmdclass,
    distclass=MyDistribution,

    zip_safe=False,
)
