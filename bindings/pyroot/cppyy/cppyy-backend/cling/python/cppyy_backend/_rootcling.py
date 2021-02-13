import os, sys, subprocess

MYHOME = os.path.dirname(__file__)

def main():
    os.environ['LD_LIBRARY_PATH'] = os.path.join(MYHOME, 'lib')
    rootcling = os.path.join(MYHOME, 'bin', 'rootcling')
    if 'win32' in sys.platform:
        rootcling += '.exe'
    if not os.path.exists(rootcling):
        raise RuntimeError("rootcling not installed in standard location")

    from ._get_cppflags import get_cppflags
    extra_flags = get_cppflags()

    rc_idx = 0
    for arg in sys.argv:
        if 'rootcling' in arg:
            break
        rc_idx += 1

    linkdef = None
    if sys.argv and 'linkdef.h' in sys.argv[-1].lower() or \
            sys.argv[-1][-4:] == '.xml':
        linkdef = sys.argv[-1]
        sysargs = sys.argv[:-1]
    else:
        sysargs = sys.argv[:]

    if extra_flags is not None:
      # rootcling is picky about order ...
        args = list()
        try:
            check_idx = rc_idx+1
            if sysargs[check_idx].find('-v', 0, 2) == 0:
                args.append(sysargs[check_idx])
                check_idx += 1
            if sysargs[check_idx] == '-f':
                args += sysargs[check_idx:check_idx+2]
                check_idx += 2
          # skip past either an output file or input header
            while sysargs[check_idx][0] != '-':
                args.append(sysargs[check_idx])
                check_idx += 1
        except IndexError:
            pass
        if args:
            args = args + ['-cxxflags', extra_flags] + sysargs[check_idx:]
        else:
            args = sysargs[rc_idx+1:] + ['-cxxflags', extra_flags]
        if linkdef:
            args.append(linkdef)
    else:
       args = sys.argv[1:]

    return subprocess.call([rootcling] + args)


if __name__ == "__main__":
    sys.exit(main())
