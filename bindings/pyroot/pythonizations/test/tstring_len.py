import unittest

import ROOT


class TStringLen(unittest.TestCase):
    """
    Test for the pythonization that provides the length of a
    TString instance `s` via `len(s)`.
    """

    # Tests
    def test_len(self):
        s = 'test'
        ts = ROOT.TString(s)
        self.assertEqual(len(ts), len(s))


if __name__ == '__main__':
    unittest.main()
