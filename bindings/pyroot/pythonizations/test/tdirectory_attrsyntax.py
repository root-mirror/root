import unittest

import ROOT
from libcppyy import SetOwnership


class TDirectoryReadWrite(unittest.TestCase):
    """
    Test for the attr syntax of TDirectory.
    """

    nbins = 8
    xmin = 0
    xmax = 4

    # Setup
    @classmethod
    def setUpClass(cls):
        cls.dir0 = ROOT.TDirectory("dir0", "dir0")
        h = ROOT.TH1F("h", "h", cls.nbins, cls.xmin, cls.xmax)
        SetOwnership(h, False)
        # this must be there otherwise the histogram is not attached to dir0
        h.SetDirectory(cls.dir0)

        dir1 = cls.dir0.mkdir("dir1")
        dir1.cd()
        h1 = ROOT.TH1F("h1", "h1", cls.nbins, cls.xmin, cls.xmax)
        SetOwnership(h1, False)

        dir2 = dir1.mkdir("dir2")
        dir2.cd()
        h2 = ROOT.TH1F("h2", "h2", cls.nbins, cls.xmin, cls.xmax)
        SetOwnership(h2, False)

    def checkHisto(self, h):
        xaxis = h.GetXaxis()
        self.assertEqual(self.nbins, h.GetNbinsX())
        self.assertEqual(self.xmin, xaxis.GetXmin())
        self.assertEqual(self.xmax, xaxis.GetXmax())

    # Tests
    def test_readHisto_attrsyntax(self):
        self.checkHisto(self.dir0.h)
        self.checkHisto(self.dir0.dir1.h1)
        self.checkHisto(self.dir0.dir1.dir2.h2)

    def test_caching_getattr(self):
        # check that __dict__ of self.dir_caching is initially empty
        self.assertFalse(self.dir0.__dict__)
        self.dir0.h
        # check that after call __dict__ is not empty anymore
        self.assertTrue(self.dir0.__dict__)
        # check that __dict__ has only one entry
        self.assertEqual(len(self.dir0.__dict__), 1)
        # check that the value in __dict__ is actually the object
        # inside the directory
        self.assertEqual(self.dir0.__dict__['h'], self.dir0.h)


if __name__ == '__main__':
    unittest.main()
