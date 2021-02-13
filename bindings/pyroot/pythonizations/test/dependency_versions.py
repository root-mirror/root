import unittest
import pkg_resources
import subprocess
import os
import sys


# Compile list of packages to be ignored in the test
ignore = []

if sys.version_info[0] == 2 and "ROOTTEST_IGNORE_NUMBA_PY2" in os.environ or \
   sys.version_info[0] == 3 and "ROOTTEST_IGNORE_NUMBA_PY3" in os.environ:
       ignore += ['numba', 'cffi']

if sys.version_info[0] == 2 and "ROOTTEST_IGNORE_JUPYTER_PY2" in os.environ or \
   sys.version_info[0] == 3 and "ROOTTEST_IGNORE_JUPYTER_PY3" in os.environ:
       ignore += ['notebook', 'metakernel']


class DependencyVersions(unittest.TestCase):
    def test_versions(self):
        '''
        Test the versions of the installed packages versus the
        requirements file in ROOT
        '''
        # For implementation details see
        # https://stackoverflow.com/questions/16294819/check-if-my-python-has-all-required-packages/45474387#45474387

        # Get source directory with requirements.txt
        p = subprocess.Popen(['root-config', '--srcdir'], stdout=subprocess.PIPE)
        r, _ = p.communicate()
        rootsrc = r.decode('UTF-8').strip()

        # Check each requirement separately
        path = os.path.join(rootsrc, 'requirements.txt')
        f = open(path)
        requirements = pkg_resources.parse_requirements(f)
        errors = []
        for requirement in requirements:
            requirement_str = str(requirement)
            name = requirement.project_name
            if name in ignore:
                print('Ignore dependency {}'.format(requirement_str))
                continue
            try:
                pkg_resources.require(requirement_str)
            except Exception as e:
                errors.append(e)
        f.close()
        if errors:
            print()
            print('Full path to requirements.txt: {}'.format(path))
            print('Details about not matched dependencies:')
            print('\n'.join([' - ' + e.report() for e in errors]))
            raise Exception('Found not matched dependencies declared in the requirements.txt, see test output for details')


if __name__ == '__main__':
    unittest.main()
