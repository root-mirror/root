import unittest

import ROOT
from libcppyy import SetOwnership


class TIterIterator(unittest.TestCase):
    """
    Test for the pythonization that allows instances of TIter to
    behave as Python iterators.
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
        # Check that TIter instances are iterable
        c = self.create_tcollection()

        itc = ROOT.TIter(c)
        # An iterator of an iterator is itself
        self.assertEqual(itc, iter(itc))

    def test_iterator(self):
        # Check that TIter instances are iterators
        c = self.create_tcollection()

        itc1 = ROOT.TIter(c)
        itc2 = ROOT.TIter(c)
        for _ in range(c.GetEntries()):
            self.assertEqual(next(itc1), itc2.Next())

    def test_for_loop_syntax(self):
        # Somehow redundant, but good to test with real syntax
        c = self.create_tcollection()

        itc1 = ROOT.TIter(c)
        itc2 = ROOT.TIter(c)
        for elem1, elem2 in zip(itc1, itc2):
            self.assertEqual(elem1, elem2)


if __name__ == '__main__':
    unittest.main()
