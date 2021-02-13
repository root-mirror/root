import unittest

import ROOT


class RooAbsCollectionLen(unittest.TestCase):
    """
    Test for the pythonization that allows to access the number of elements of a
    RooAbsCollection (or subclass) by calling `len` on it.
    """

    num_elems = 3
    rooabsarg_list = [ ROOT.RooStringVar(str(i),'','') for i in range(num_elems) ]
    tlist = ROOT.TList()

    # Setup
    @classmethod
    def setUpClass(cls):
        for elem in cls.rooabsarg_list:
            cls.tlist.Add(elem)

    @classmethod
    def tearDownClass(cls):
        # Clear TList before Python list deletes the objects
        cls.tlist.Clear()

    # Helpers
    def check_len(self, c):
        self.assertEqual(len(c), self.num_elems)
        self.assertEqual(len(c), c.getSize())

    # Tests
    def test_rooarglist(self):
        self.check_len(ROOT.RooArgList(self.tlist))

    def test_rooargset(self):
        self.check_len(ROOT.RooArgSet(self.tlist))


if __name__ == '__main__':
    unittest.main()
