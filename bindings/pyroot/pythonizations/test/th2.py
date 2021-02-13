import unittest

import ROOT


class TH2UsingDecls(unittest.TestCase):
    """
    Test that the method pulled in via using decls from TH1 are accessible
    """

    # Tests
    def test_GetBinError(self):
        h = ROOT.TH2F("h", "h", 1, 1, 0, 1, 0, 1)
        for _ in range(4):
            h.Fill(1, 1, 1)
        self.assertEqual(h.GetBinErrorUp(1, 1), 2)
        self.assertEqual(h.GetBinErrorLow(1, 1), 2)


if __name__ == '__main__':
    unittest.main()
