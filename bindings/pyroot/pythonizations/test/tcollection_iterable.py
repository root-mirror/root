import unittest

import ROOT
from libcppyy import SetOwnership


class TCollectionIterable(unittest.TestCase):
    """
    Test for the pythonization that makes instances of TCollection subclasses
    iterable in Python.
    For example, this allows to do:
    `for elem in collection:`
         `...`
    """

    num_elems = 3

    # Helpers
    def create_tcollection(self):
        c = ROOT.TList()
        for _ in range(self.num_elems):
            o = ROOT.TObject()
            # Prevent immediate deletion of C++ TObjects
            SetOwnership(o, False)
            c.Add(o)

        return c

    # Tests
    def test_iterable(self):
        c = self.create_tcollection()

        itc = ROOT.TIter(c)
        for elem in c:
            self.assertEqual(elem, itc.Next())


if __name__ == '__main__':
    unittest.main()

