import unittest

import ROOT
from libcppyy import SetOwnership


class TFileOpenReadWrite(unittest.TestCase):
    """
    Test for the TFile.Open factory like creation of TFile
    """

    filename  = 'tfileopenreadwrite.root'
    nbins = 8
    xmin = 0
    xmax = 4

    # Setup
    @classmethod
    def setUpClass(cls):
        f = ROOT.TFile.Open(cls.filename, "RECREATE")
        h = ROOT.TH1F("h", "h", cls.nbins, cls.xmin, cls.xmax)
        SetOwnership(h, False)
        f.WriteObject(h, "h")

        dir1 = f.mkdir("dir1")
        dir1.cd()
        h1 = ROOT.TH1F("h1", "h1", cls.nbins, cls.xmin, cls.xmax)
        SetOwnership(h1, False)
        h1.Write()

        dir2 = dir1.mkdir("dir2")
        dir2.cd()
        h2 = ROOT.TH1F("h2", "h2", cls.nbins, cls.xmin, cls.xmax)
        SetOwnership(h2, False)
        h2.Write()

        f.Close()

    def checkHisto(self, h):
        xaxis = h.GetXaxis()
        self.assertEqual(self.nbins, h.GetNbinsX())
        self.assertEqual(self.xmin, xaxis.GetXmin())
        self.assertEqual(self.xmax, xaxis.GetXmax())

    # Tests
    def test_readHisto_attrsyntax(self):
        f = ROOT.TFile.Open(self.filename)
        self.checkHisto(f.h)
        self.checkHisto(f.dir1.h1)
        self.checkHisto(f.dir1.dir2.h2)

    def test_readHisto(self):
        f = ROOT.TFile.Open(self.filename)
        self.checkHisto(f.Get("h"))
        self.checkHisto(f.Get("dir1/h1"))
        self.checkHisto(f.Get("dir1/dir2/h2"))

    def test_caching_getattr(self):
        f = ROOT.TFile.Open(self.filename)
        # check that __dict__ of self.dir_caching is initially empty
        self.assertFalse(f.__dict__)
        f.h
        # check that after call is not empty anymore
        self.assertTrue(f.__dict__)
        # check that __dict__ has only one entry
        self.assertEqual(len(f.__dict__), 1)
        # check that the value in __dict__ is actually the object
        # inside the directory
        self.assertEqual(f.__dict__['h'], f.h)

    def test_oserror(self):
        # check that an OSError is raised when an inexistent file is opened
        # both with a string and an instance of TFileOpenHandle as arguments
        self.assertRaises(OSError, ROOT.TFile.Open, 'inexistent_file.root')
        handle = ROOT.TFile.AsyncOpen("inexistent_file.root")
        self.assertRaises(OSError, ROOT.TFile.Open, handle)


if __name__ == '__main__':
    unittest.main()
