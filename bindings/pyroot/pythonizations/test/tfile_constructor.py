import unittest

import ROOT


class TFileConstructor(unittest.TestCase):
    """
    Test for the TFile constructor
    """

    def test_oserror(self):
        # check that an OSError is raised when the string passed as argument
        # refers to an inexistent file
        self.assertRaises(OSError, ROOT.TFile, 'inexistent_file.root')


if __name__ == '__main__':
    unittest.main()
