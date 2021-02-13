import unittest

import ROOT
from libcppyy import SetOwnership


class TObjectContains(unittest.TestCase):
    """
    Test for the __contains__ pythonisation of TObject and subclasses.
    Such pythonisation relies on TObject::FindObject, which is redefined
    in some of its subclasses, such as TCollection.
    Thanks to this pythonisation, we can use the syntax `obj in col`
    to know if col contains obj.
    """

    num_elems = 3

    # Helpers
    def create_tlist(self):
        l = ROOT.TList()
        for _ in range(self.num_elems):
            o = ROOT.TObject()
            # Prevent immediate deletion of C++ TObjects
            SetOwnership(o, False)
            l.Add(o)

        return l

    # Tests
    def test_contains(self):
        l = self.create_tlist()

        for elem in l:
            self.assertTrue(elem in l)
            # Make sure it does not work just because of __iter__
            self.assertTrue(l.__contains__(elem))

        o = ROOT.TObject()
        self.assertFalse(o in l)
        self.assertFalse(l.__contains__(o))


if __name__ == '__main__':
    unittest.main()
