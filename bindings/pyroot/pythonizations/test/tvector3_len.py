import unittest

import ROOT


class TVector3Len(unittest.TestCase):
    """
    Test for the pythonization that allows to get the size of a
    TVector3 (always 3) by calling `len` on it.
    """

    # Tests
    def test_len(self):
        v = ROOT.TVector3(1., 2., 3.)
        self.assertEqual(len(v), 3)


if __name__ == '__main__':
    unittest.main()
