from __future__ import print_function
import py, sys, subprocess

currpath = py.path.local(__file__).dirpath()


def setup_make(targetname):
    if sys.platform == 'win32':
        raise OSError("win32 not supported yet")
    popen = subprocess.Popen(["make", targetname], cwd=str(currpath),
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, _ = popen.communicate()
    if popen.returncode:
        raise OSError("'make' failed:\n%s" % (stdout,))
