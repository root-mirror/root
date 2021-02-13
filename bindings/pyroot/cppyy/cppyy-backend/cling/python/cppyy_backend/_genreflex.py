from __future__ import print_function
import os, sys, subprocess

MYHOME = os.path.dirname(__file__)

def main():
    if len(sys.argv) == 2 and \
           (sys.argv[1] == '--cflags' or sys.argv[1] == '--cppflags'):
        print('-I%s/include' % (MYHOME,))
        return 0
    os.environ['LD_LIBRARY_PATH'] = os.path.join(MYHOME, 'lib')
    genreflex = os.path.join(MYHOME, 'bin', 'genreflex')
    if 'win32' in sys.platform:
        genreflex += '.exe'
    if not os.path.exists(genreflex):
        raise RuntimeError("genreflex not installed in standard location")

    from ._get_cppflags import get_cppflags
    extra_flags = get_cppflags()
    if extra_flags is not None and 1 < len(sys.argv):
      # genreflex is picky about order ...
       args = [sys.argv[1], '--cxxflags', extra_flags] + sys.argv[2:]
    else:
       args = sys.argv[1:]

    return subprocess.call([genreflex] + args)


if __name__ == "__main__":
    sys.exit(main())
