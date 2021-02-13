#!/usr/bin/env python
from __future__ import print_function

import os, sys, subprocess
import shutil, tarfile
try:
    import urllib2
except ModuleNotFoundError:
    import urllib.request as urllib2  # p3

if 'win32' in sys.platform:
    def rename(f1, f2):
        os.remove(f2)
        os.rename(f1, f2)
else:
    rename = os.rename

def is_manylinux():
    _is_manylinux = False
    try:
        for line in open('/etc/redhat-release').readlines():
            if 'CentOS release 6.10 (Final)' in line:
                _is_manylinux = True
                break
    except (OSError, IOError):
        pass
    return _is_manylinux


DEBUG_TESTBUILD = False

TARBALL_CACHE_DIR = 'releases'

ROOT_KEEP = ['build', 'cmake', 'config', 'core', 'etc', 'interpreter',
             'io', 'LICENSE', 'LGPL2_1.txt', 'Makefile', 'CMakeLists.txt', 'math',
             'main', # main only needed in more recent root b/c of rootcling
             'builtins']
ROOT_CORE_KEEP = ['CMakeLists.txt', 'base', 'clib', 'clingutils', 'cont',
                  'dictgen', 'foundation', 'macosx', 'meta',
                  'metacling', 'metautils', 'rootcling_stage1', 'textinput',
                  'thread', 'unix', 'utils', 'winnt', 'zip']
ROOT_BUILTINS_KEEP = ['openssl', 'pcre', 'xxhash', 'zlib']
ROOT_IO_KEEP = ['CMakeLists.txt', 'io', 'rootpcm']
ROOT_MATH_KEEP = ['CMakeLists.txt', 'mathcore']
ROOT_ETC_KEEP = ['Makefile.arch', 'class.rules', 'cmake', 'dictpch',
                 'gdb-backtrace.sh', 'gitinfo.txt', 'helgrind-root.supp',
                 'hostcert.conf', 'plugins', 'system.plugins-ios',
                 'valgrind-root-python.supp', 'valgrind-root.supp', 'vmc']
ROOT_PLUGINS_KEEP = ['TVirtualStreamerInfo']

ROOT_EXPLICIT_REMOVE = [os.path.join('core', 'base', 'v7'),
                        os.path.join('math', 'mathcore', 'v7'),
                        os.path.join('io', 'io', 'v7')]


ERR_RELEASE_NOT_FOUND = 2


def get_root_version(try_recover=True):
    import pkg_resources
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'python')
    dists = pkg_resources.find_distributions(path)
    try:
        cppyy_cling = [d for d in dists if d.key == 'cppyy-cling'][0]
        version = cppyy_cling.version
    except IndexError:
        if try_recover and os.path.exists('setup.py'):
            print('No egg_info ... running "python setup.py egg_info"')
            if subprocess.call(['python', 'setup.py', 'egg_info']) != 0:
                print('ERROR: creation of egg_info failed ... giving up')
                sys.exit(2)
            return get_root_version(False)
        else:
            print('ERROR: cannot determine version. Please run "python setup.py egg_info" first.')
            sys.exit(1)
    #
    parts = version.split('.', 3)
    major, minor, patch = map(int, parts[:3])
    root_version = '%d.%02d.%02d' % (major, minor, patch)
    return root_version


ROOT_VERSION = get_root_version()

#
## ROOT source pull and cleansing
#
def clean_directory(directory, keeplist, trim_cmake=True):
    removed_entries = []
    for entry in os.listdir(directory):
        if entry[0] == '.' or entry in keeplist:
            continue
        removed_entries.append(entry)
        entry = os.path.join(directory, entry)
        print('now removing', entry)
        if os.path.isdir(entry):
            shutil.rmtree(entry)
        else:
            os.remove(entry)

    if not trim_cmake:
        return

    # now take the removed entries out of the CMakeLists.txt
    inp = os.path.join(directory, 'CMakeLists.txt')
    if removed_entries and os.path.exists(inp):
        print('trimming', inp)
        outp = inp+'.new'
        new_cml = open(outp, 'w')
        for line in open(inp).readlines():
            if ('add_subdirectory' in line) or\
               ('COMMAND' in line and 'copy' in line) or\
               ('ROOT_ADD_TEST_SUBDIRECTORY' in line) or\
               ('install(DIRECTORY' in line):
                for sub in removed_entries:
                    if sub in line:
                        line = '#'+line
                        break
            new_cml.write(line)
        new_cml.close()
        rename(outp, inp)
    else:
        print('reusing existing %s/CMakeLists.txt' % (directory,))
 

if not os.path.exists(TARBALL_CACHE_DIR):
    os.mkdir(TARBALL_CACHE_DIR)

fn = 'root_v%s.source.tar.gz' % ROOT_VERSION
addr = 'https://root.cern.ch/download/'+fn
if not os.path.exists(os.path.join(TARBALL_CACHE_DIR, fn)):
    try:
        print('retrieving', fn)
        if sys.hexversion < 0x3000000:
            output_fn = fn
        else:
            output_fn = bytes(fn, 'utf-8')
        resp = urllib2.urlopen(addr, output_fn)
        out = open(os.path.join(TARBALL_CACHE_DIR, fn), 'wb')
        out.write(resp.read())
        out.close()
    except urllib2.HTTPError:
        print('release %s not found' % ROOT_VERSION)
        sys.exit(ERR_RELEASE_NOT_FOUND)
else:
    print('reusing', fn, 'from local directory')


fn = os.path.join(TARBALL_CACHE_DIR, fn)
pkgdir = os.path.join('root-' + ROOT_VERSION)
if not os.path.exists(pkgdir):
    print('now extracting', ROOT_VERSION)
    tf = tarfile.TarFile.gzopen(fn)
    tf.extractall()
    tf.close()
else:
    print('reusing existing directory', pkgdir)

# remove everything except for the listed set of libraries
try:
    shutil.rmtree('src')
except OSError:
    pass
os.chdir(pkgdir)
clean_directory(os.path.curdir,                  ROOT_KEEP)
clean_directory('core',                          ROOT_CORE_KEEP)
clean_directory('builtins',                      ROOT_BUILTINS_KEEP)
clean_directory('etc',                           ROOT_ETC_KEEP, trim_cmake=False)
clean_directory(os.path.join('etc', 'plugins'),  ROOT_PLUGINS_KEEP, trim_cmake=False)
clean_directory('io',                            ROOT_IO_KEEP)
clean_directory('math',                          ROOT_MATH_KEEP)


# trim main (only need rootcling)
print('trimming main')
for entry in os.listdir(os.path.join('main', 'src')):
    if entry != 'rootcling.cxx':
        os.remove(os.path.join('main', 'src', entry))
inp = os.path.join('main', 'CMakeLists.txt')
outp = inp+'.new'
new_cml = open(outp, 'w')
now_stripping = False
for line in open(inp).readlines():
    if ('ROOT_EXECUTABLE' in line or 'root.exe' in line or \
        'SET_TARGET_PROPERTIES' in line) and\
       not 'rootcling' in line:
        line = '#'+line
    elif '#---CreateHadd' == line[0:14]:
        now_stripping = True
    elif now_stripping and line == ')\n':
        line = '#'+line
        now_stripping = False
    if now_stripping:
        line = '#'+line
    new_cml.write(line)
new_cml.close()
rename(outp, inp)

# trim core (remove lz4 and lzma, and gGLManager which otherwise crashes on call)
print('trimming core/CMakeLists.txt')
inp = os.path.join('core', 'CMakeLists.txt')
outp = inp+'.new'
new_cml = open(outp, 'w')
for line in open(inp).readlines():
    if ('Lzma' in line or 'Lz4' in line):
        line = '#'+line
    else:
        line = line.replace(' ${LZMA_LIBRARIES}', '')
        line = line.replace(' LZ4::LZ4', '')
        line = line.replace(' LZMA', '')
    new_cml.write(line)
new_cml.close()
rename(outp, inp)

print('trimming core/base')
os.remove(os.path.join('core', 'base', 'src', 'TVirtualGL.cxx'))
os.remove(os.path.join('core', 'base', 'inc', 'TVirtualGL.h'))
inp = os.path.join('core', 'base', 'CMakeLists.txt')
outp = inp+'.new'
now_stripping = False
new_cml = open(outp, 'w')
for line in open(inp).readlines():
    if 'if(cxx14 OR cxx17 OR root7)' == line[0:27]:   # get rid of v7 stuff
        now_stripping = True
    elif 'if(root7)' == line[0:9]:
        now_stripping = False
    elif 'RLogger' in line or \
         ('BASE_HEADER_DIRS' in line and 'v7' in line) or \
         'TVirtualGL' in line or \
         ('argparse' in line and 'Man' in line):
        line = '#'+line
    if now_stripping:
        line = '#'+line
    new_cml.write(line)
new_cml.close()
rename(outp, inp)

# do not copy wchar.h & friends b/c the pch should be generated at install time,
# so preventing conflict
print('trimming clingutils')
inp = os.path.join('core', 'clingutils', 'CMakeLists.txt')
outp = inp+'.new'
now_stripping = False
new_cml = open(outp, 'w')
for line in open(inp).readlines():
    if '# Capture their build-time' == line[0:26]:
        now_stripping = True
    elif 'set(stamp_file' == line[0:14]:
        now_stripping = False
    if now_stripping:
        line = '#'+line
    new_cml.write(line)
new_cml.close()
rename(outp, inp)

# remove extraneous external software explicitly
print('trimming externals')
for cmf in ['AfterImage', 'FTGL']:
    fname = os.path.join('cmake', 'modules', 'Find%s.cmake' % (cmf,))
    if os.path.exists(fname):
        os.remove(fname)
inp = os.path.join('cmake', 'modules', 'SearchInstalledSoftware.cmake')
outp = inp+'.new'
now_stripping = False
new_cml = open(outp, 'w')
for line in open(inp).readlines():
    if '#---Check for ftgl if needed' == line[0:28] or\
       '#---Check for AfterImage' == line[0:24] or\
       '#---Check for Freetype' == line[0:22] or\
       '#---Check for LZMA' == line[0:18] or\
       '#---Check for LZ4' == line[0:17] or\
       '#-------' == line[0:8]:   # openui5 (doesn't follow convention)
        now_stripping = True
    elif '#---Check' == line[0:9] or\
         '#---Report' == line[0:10]:
        now_stripping = False
    if now_stripping or 'builtin_freetype' in line:
        line = '#'+line
    new_cml.write(line)
new_cml.close()
rename(outp, inp)

inp = os.path.join('cmake', 'modules', 'RootBuildOptions.cmake')
outp = inp+'.new'
new_cml = open(outp, 'w')
for line in open(inp).readlines():
    if 'ROOT_BUILD_OPTION(builtin_ftgl' in line or\
       'ROOT_BUILD_OPTION(builtin_afterimage' in line:
        line = '#'+line
    new_cml.write(line)
new_cml.close()
rename(outp, inp)

# remove testing, examples, and notebook
print('trimming testing')
inp = 'CMakeLists.txt'
outp = inp+'.new'
now_stripping = False
new_cml = open(outp, 'w')
for line in open(inp).readlines():
    if '#---Configure Testing using CTest' == line[0:33] or\
       '#---hsimple.root' == line[0:16]:
        now_stripping = True
    elif '#---Packaging' == line[0:13] or\
         '#---version' == line[0:11]:
        now_stripping = False
    if now_stripping:
        line = '#'+line
    elif 'root_kernel' in line:
        line = '#'+line
    new_cml.write(line)
new_cml.close()
rename(outp, inp)

print('trimming RootCPack')
inp = os.path.join('cmake', 'modules', 'RootCPack.cmake')
outp = inp+'.new'
new_cml = open(outp, 'w')
for line in open(inp):
    if 'README.txt' in line:
        line = '#'+line
    new_cml.write(line)
new_cml.close()
rename(outp, inp)


# some more explicit removes:
for dir_to_remove in ROOT_EXPLICIT_REMOVE:
    try:
        shutil.rmtree(dir_to_remove)
    except OSError:
        pass


# special fixes
for inp in [os.path.join('core', 'unix', 'src', 'TUnixSystem.cxx'),
            os.path.join('core', 'winnt', 'src', 'TWinNTSystem.cxx')]:
    outp = inp+'.new'
    new_cml = open(outp, 'w')
    for line in open(inp):
        if '#include "TSocket.h"' == line[0:20]:
            line = """//#include "TSocket.h"
enum ESockOptions {
   kSendBuffer,        // size of send buffer
   kRecvBuffer,        // size of receive buffer
   kOobInline,         // OOB message inline
   kKeepAlive,         // keep socket alive
   kReuseAddr,         // allow reuse of local portion of address 5-tuple
   kNoDelay,           // send without delay
   kNoBlock,           // non-blocking I/O
   kProcessGroup,      // socket process group (used for SIGURG and SIGIO)
   kAtMark,            // are we at out-of-band mark (read only)
   kBytesToRead        // get number of bytes to read, FIONREAD (read only)
};

enum ESendRecvOptions {
   kDefault,           // default option (= 0)
   kOob,               // send or receive out-of-band data
   kPeek,              // peek at incoming message (receive only)
   kDontBlock          // send/recv as much data as possible without blocking
};
"""
        new_cml.write(line)
    new_cml.close()
    rename(outp, inp)

print('trimming mathcore')
os.remove(os.path.join('math', 'mathcore', 'src', 'triangle.h'))
os.remove(os.path.join('math', 'mathcore', 'src', 'triangle.c'))
os.remove(os.path.join('math', 'mathcore', 'inc', 'Math', 'Delaunay2D.h'))
os.remove(os.path.join('math', 'mathcore', 'src', 'Delaunay2D.cxx'))

inp = os.path.join('math', 'mathcore', 'CMakeLists.txt')
outp = inp+'.new'
now_stripping = False
new_cml = open(outp, 'w')
for line in open(inp).readlines():
    if 'Delaunay2D' in line or 'triangle' in line:
        new_cml.write('#')
    new_cml.write(line)
new_cml.close()
rename(outp, inp)

# done
os.chdir(os.path.pardir)

# debugging: run a test build
if DEBUG_TESTBUILD:
    print('running a debug test build')
    tb = "test_builddir"
    if os.path.exists(tb):
        shutil.rmtree(tb)
    os.mkdir(tb)
    os.chdir(tb)
    os.system('cmake %s -DCMAKE_INSTALL_PREFIX=%s -Dminimal=ON -Dasimage=OFF' % \
                 (os.path.join(os.pardir, pkgdir), os.path.join(os.pardir, 'install')))
    os.system('make -j 32')


#
## package creation
#
countdown = 0

print('adding src ... ')
if not os.path.exists('src'):
    os.mkdir('src')
for entry in os.listdir(pkgdir):
    fullp = os.path.join(pkgdir, entry)
    if entry[0] == '.':
        continue
    dest = os.path.join('src', entry)
    if os.path.isdir(fullp):
        if not os.path.exists(dest):
            shutil.copytree(fullp, dest)
    else:
        if not os.path.exists(dest):
            shutil.copy2(fullp, dest)

#
## apply patches (in order)
#
try:
    import patch
except ImportError:
    class patch(object):
        @staticmethod
        def fromfile(fdiff):
            return patch(fdiff)

        def __init__(self, fdiff):
            self.fdiff = fdiff

        def apply(self):
            res = os.system('patch -p1 < ' + self.fdiff)
            return res == 0

for fdiff in ('cleanup_tstring', 'scanner', 'scanner_2', 'faux_typedef', 'classrules', 'template_fwd',
              'dep_template', 'no_long64_t', 'using_decls', 'sfinae', 'typedef_of_private',
              'optlevel2_forced', 'silence', 'explicit_template', 'alias_template', 'lambda', 'templ_ops',
              'private_type_args', 'incomplete_types', 'clang_printing', 'resolution',
              'stdfunc_printhack', 'anon_union', 'no_inet', 'signaltrycatch', 'nofastmath', 'pch',
              'stackoverflow', 'stdvalue_type', 'strip_lz4_lzma',
              'msvc', 'win64rtti', 'win64', 'win64s2'):
    fpatch = os.path.join('patches', fdiff+'.diff')
    print(' ==> applying patch:', fpatch)
    pset = patch.fromfile(fpatch)
    if not pset or not pset.apply():
        print("Failed to apply patch:", fdiff)
        sys.exit(2)

#
## manylinux1 specific patch, as there a different, older, compiler is used
#
if is_manylinux():
    print(' ==> applying patch:', 'manylinux1')
    patch.fromfile(os.path.join('patches', 'manylinux1.diff')).apply()

#
## finally, remove the ROOT source directory, as it can not be reused
#
print("removing", pkgdir)
shutil.rmtree(pkgdir)

# done!
