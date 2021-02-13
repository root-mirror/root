from __future__ import print_function
""" cppyy_backend loader
"""

__all__ = [
    'load_cpp_backend',           # load libcppyy_backend
    'set_cling_compile_options',  # set EXTRA_CLING_ARGS envar
    'ensure_precompiled_header'   # build precompiled header as necessary
]

import os, sys, ctypes, subprocess, sysconfig, warnings

if 'win32' in sys.platform:
    soext = '.dll'
else:
    soext = '.so'

soabi = sysconfig.get_config_var("SOABI")


def _load_helper(bkname):
 # normal load, allowing for user overrides of LD_LIBRARY_PATH
    try:
        return ctypes.CDLL(bkname, ctypes.RTLD_GLOBAL)
    except OSError:
         pass

 # failed ... load dependencies explicitly
    try:
        pkgpath = os.path.dirname(bkname)
        if not pkgpath:
            pkgpath = os.path.dirname(__file__)
        elif os.path.basename(pkgpath) in ['lib', 'bin']:
            pkgpath = os.path.dirname(pkgpath)
        for dep in ['liblzma', 'libCore', 'libThread', 'libRIO', 'libCling']:
            for loc in ['lib', 'bin']:
                fpath = os.path.join(pkgpath, loc, dep+soext)
                if os.path.exists(fpath):
                    ldtype = ctypes.RTLD_GLOBAL
                    if dep == 'libCling': ldtype = ctypes.RTLD_LOCAL
                    ctypes.CDLL(fpath, ldtype)
                    break
        return ctypes.CDLL(os.path.join(pkgpath, 'lib', bkname), ctypes.RTLD_GLOBAL)
    except OSError:
        pass

    return None


_precompiled_header_ensured = False
def load_cpp_backend():
    set_cling_compile_options()

    if not _precompiled_header_ensured:
     # the precompiled header of standard and system headers is not part of the
     # distribution as there are too many varieties; create it now if needed
        ensure_precompiled_header()

    altbkname = None
    try:
        bkname = os.environ['CPPYY_BACKEND_LIBRARY']
        if bkname.rfind(soext) < 0:
            bkname += soext
    except KeyError:
        bkname = 'libcppyy_backend'+soext
        if soabi:
            altbkname = 'libcppyy_backend.'+soabi+soext

    c = _load_helper(bkname)
    if not c and altbkname is not None:
        c = _load_helper(altbkname)

    if not c:
        raise RuntimeError("could not load cppyy_backend library")

    return c


def set_cling_compile_options(add_defaults = False):
 # extra optimization flags for Cling
    if not 'EXTRA_CLING_ARGS' in os.environ:
        CURRENT_ARGS = ''
        add_defaults = True
    else:
        CURRENT_ARGS = os.environ['EXTRA_CLING_ARGS']

    if add_defaults:
        has_avx = False
        try:
            with open('/proc/cpuinfo', 'r') as ci:
                for line in ci:
                    if 'avx' in line:
                        has_avx = True
                        break
        except Exception:
            try:
                cli_arg = subprocess.check_output(['sysctl', 'machdep.cpu.features'])
                has_avx = 'avx' in cli_arg.decode("utf-8").strip().lower()
            except Exception:
                pass
        CURRENT_ARGS += ' -O2'
        if has_avx: CURRENT_ARGS += ' -mavx'
        os.putenv('EXTRA_CLING_ARGS', CURRENT_ARGS)
        os.environ['EXTRA_CLING_ARGS'] = CURRENT_ARGS


def _disable_pch():
    os.putenv('CLING_STANDARD_PCH', 'none')

def _warn_no_pch(msg, pchname=None):
    if pchname is None or not os.path.exists(pchname):
        _disable_pch()
        warnings.warn('No precompiled header available (%s); this may impact performance.' % msg)
    else:
        warnings.warn('Precompiled header may be out of date (%s).' % msg)

def _is_uptodate(pchname, incpath):
  # test whether the pch is older than the include directory
    try:
        return os.stat(pchname).st_mtime >= os.stat(incpath).st_mtime
    except Exception:
        if not os.path.exists(incpath):
            return True     # no point in updating as it will fail
    return False

def ensure_precompiled_header(pchdir = '', pchname = ''):
  # the precompiled header of standard and system headers is not part of the
  # distribution as there are too many varieties; create it now if needed
     global _precompiled_header_ensured
     _precompiled_header_ensured = True     # only ever call once, even if failed

     olddir = os.getcwd()
     try:
         pkgpath = os.path.abspath(os.path.dirname(__file__))
         if 'CLING_STANDARD_PCH' in os.environ:
             stdpch = os.environ['CLING_STANDARD_PCH']
             if stdpch.lower() == 'none':   # magic keyword to disable pch
                 _disable_pch()
                 os.chdir(olddir)
                 return                     # quiet
             pchdir = os.path.dirname(stdpch)
             pchname = os.path.basename(stdpch)
         else:
             if not pchdir:
                 pchdir = os.path.join(pkgpath, 'etc')
             if not pchname:
                 pchname = 'allDict.cxx.pch'

         os.chdir(pkgpath)
         full_pchname = os.path.join(pchdir, pchname)
         incpath = os.path.join(pkgpath, 'include')
         is_uptodate = _is_uptodate(full_pchname, incpath)
         if not os.path.exists(full_pchname) or not is_uptodate:
             if os.access(pchdir, os.R_OK|os.W_OK):
                 print('(Re-)building pre-compiled headers (options:%s); this may take a minute ...' % os.environ.get('EXTRA_CLING_ARGS', ' none'))
                 makepch = os.path.join(pkgpath, 'etc', 'dictpch', 'makepch.py')
                 if subprocess.call([sys.executable, makepch, full_pchname, '-I'+incpath]) != 0:
                     _warn_no_pch('failed to build', full_pchname)
             else:
                 _warn_no_pch('%s not writable' % pchdir, full_pchname)

     except Exception as e:
         _warn_no_pch(str(e))
     finally:
         os.chdir(olddir)
