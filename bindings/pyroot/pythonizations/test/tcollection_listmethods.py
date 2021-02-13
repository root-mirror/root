import unittest

import ROOT
from libcppyy import SetOwnership


class TCollectionListMethods(unittest.TestCase):
    """
    Test for the Python-list-like methods added to TCollection (and subclasses):
    append, remove, extend, count
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
    def test_append(self):
        c = self.create_tcollection()

        o = ROOT.TObject()
        self.assertFalse(c.Contains(o))
        len1 = c.GetEntries()

        c.append(o)

        len2 = c.GetEntries()
        self.assertEqual(len1 + 1, len2)
        self.assertTrue(c.Contains(o))

        # Skip elements that were already there
        itc = ROOT.TIter(c)
        for _ in range(len1):
            itc.Next()

        # Check that `o` is indeed the last element
        self.assertEqual(o, itc.Next())

    def test_remove(self):
        c = ROOT.TList()

        o1 = ROOT.TObject()
        o2 = ROOT.TObject()

        c.Add(o1)
        c.Add(o2)
        c.Add(o1)

        self.assertTrue(c.Contains(o1))
        self.assertEqual(c.GetEntries(), 3)

        c.remove(o1)

        self.assertEqual(c.GetEntries(), 2)

        c.remove(o1)

        self.assertEqual(c.GetEntries(), 1)
        self.assertFalse(c.Contains(o1))

        with self.assertRaises(ValueError):
            c.remove(o1)

    def test_extend(self):
        c1 = self.create_tcollection()
        c2 = self.create_tcollection()

        len1 = c1.GetEntries()
        len2 = c2.GetEntries()

        c1.extend(c2)

        len1_final = c1.GetEntries()

        self.assertEqual(len1_final, len1 + len2)

        # Skip elements that were already there
        itc1 = ROOT.TIter(c1)
        for _ in range(len1):
            itc1.Next()

        # Compare with elements of second collection
        itc2 = ROOT.TIter(c2)
        for _ in range(len2):
            self.assertEqual(itc1.Next(), itc2.Next())

    def test_count(self):
        c = ROOT.TList()

        o1 = ROOT.TObject()
        o2 = ROOT.TObject()

        c.Add(o1)
        c.Add(o2)
        c.Add(o1)

        self.assertEqual(c.count(o1), 2)
        self.assertEqual(c.count(o2), 1)


if __name__ == '__main__':
    unittest.main()

