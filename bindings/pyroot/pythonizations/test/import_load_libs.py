import unittest
import re
import os


class ImportLoadLibs(unittest.TestCase):
    """
    Test which libraries are loaded during importing ROOT
    """

    # The whitelist is a list of regex expressions that mark wanted libraries
    # Note that the regex has to result in an exact match with the library name.
    known_libs = [
            # libCore and dependencies
            'libCore',
            'libm',
            'liblz4',
            'libxxhash',
            'liblzma',
            'libzstd',
            'libz',
            'libpthread',
            'libc',
            'libdl',
            'libpcre',
            # libCling and dependencies
            'libCling.*',
            'librt',
            'libncurses.*',
            'libtinfo', # by libncurses (on some older platforms)
            # libTree and dependencies
            'libTree',
            'libThread',
            'libRIO',
            'libNet',
            'libImt',
            'libMathCore',
            'libMultiProc',
            'libssl',
            'libcrypt.*', # by libssl
            'libtbb',
            'liburing', # by libRIO if uring option is enabled
            # On centos7 libssl links against kerberos pulling in all dependencies below, removed with libssl1.1.0
            'libgssapi_krb5',
            'libkrb5',
            'libk5crypto',
            'libkrb5support',
            'libselinux',
            'libkeyutils',
            'libcom_err',
            'libresolv',
            # cppyy and Python libraries
            'libcppyy.*',
            'libROOTPythonizations.*',
            'libpython.*',
            'libutil.*',
            '.*cpython.*',
            '_.*',
            '.*module',
            'operator',
            'cStringIO',
            'binascii',
            'libbz2',
            'libexpat',
            'ISO8859-1',
            # System libraries and others
            'libnss_.*',
            'ld.*',
            'libffi',
            ]

    # Verbose mode of the test
    verbose = False

    def test_import(self):
        """
        Test libraries loaded after importing ROOT
        """
        import ROOT
        libs = str(ROOT.gSystem.GetLibraries())

        if self.verbose:
            print("Initial output from ROOT.gSystem.GetLibraries():\n" + libs)

        # Split paths
        libs = libs.split(' ')

        # Get library name without full path and .so* suffix
        libs = [os.path.basename(l).split('.so')[0] for l in libs \
                if not l.startswith('-l') and not l.startswith('-L')]

        # Check that the loaded libraries are white listed
        bad_libs = []
        good_libs = []
        matched_re = []
        for l in libs:
            matched = False
            for r in self.known_libs:
                m = re.match(r, l)
                if m:
                    if m.group(0) == l:
                        matched = True
                        good_libs.append(l)
                        matched_re.append(r)
                        break
            if not matched:
                bad_libs.append(l)

        if self.verbose:
            print('Found whitelisted libraries after importing ROOT with the shown regex match:')
            for l, r in zip(good_libs, matched_re):
                print(' - {} ({})'.format(l, r))
            import sys
            sys.stdout.flush()

        if bad_libs:
            raise Exception('Found not whitelisted libraries after importing ROOT:' \
                    + '\n - ' + '\n - '.join(bad_libs) \
                    + '\nIf the test fails with a library that is loaded on purpose, please add it to the whitelist.')


if __name__ == '__main__':
    unittest.main()
