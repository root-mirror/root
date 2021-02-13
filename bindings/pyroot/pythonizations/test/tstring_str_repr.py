import unittest

import ROOT


class TStringStrRepr(unittest.TestCase):
    """
    Test for the pythonizations that provide a string representation
    for instances of TString (__str__, __repr__).
    """

    # Tests
    def test_str(self):
        s = 'test'
        ts = ROOT.TString(s)
        self.assertEqual(str(ts), s)

    def test_repr(self):
        s = 'test'
        ts = ROOT.TString(s)
        self.assertEqual(repr(ts), repr(s))


if __name__ == '__main__':
    unittest.main()
