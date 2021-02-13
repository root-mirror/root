import unittest

import ROOT
from libcppyy import SetOwnership


class TSeqCollectionItemAccess(unittest.TestCase):
    """
    Test for the item access methods added to TSeqCollection (and subclasses):
    __getitem__, __setitem__, __delitem__.
    Both the index (l[i]) and slice (l[i:j:k]) syntaxes are tested.
    """

    num_elems = 3

    # Helpers
    def create_tseqcollection(self):
        sc = ROOT.TList()
        for _ in range(self.num_elems):
            o = ROOT.TObject()
            # Prevent immediate deletion of C++ TObjects
            SetOwnership(o, False)
            sc.Add(o)

        return sc

    # Tests
    def test_getitem(self):
        sc = self.create_tseqcollection()

        # Get items
        it = ROOT.TIter(sc)
        for i in range(self.num_elems):
            self.assertEqual(it.Next(), sc[i])

        # Get items, negative indices
        it2 = ROOT.TIter(sc)
        neg_idcs = [ -i-1 for i in reversed(range(self.num_elems)) ]
        for i in neg_idcs:
            self.assertEqual(it2.Next(), sc[i])

        # Check invalid index cases
        with self.assertRaises(IndexError):
            sc[self.num_elems]

        with self.assertRaises(IndexError):
            sc[-(self.num_elems + 1)]

        with self.assertRaises(TypeError):
            sc[1.0]

    def test_getitem_slice(self):
        sc = self.create_tseqcollection()

        # All items
        slice1 = sc[:]
        for i in range(slice1.GetEntries()):
            self.assertEqual(sc[i], slice1[i])

        # First two items
        slice2 = sc[0:2]
        self.assertEqual(sc[0], slice2[0])
        self.assertEqual(sc[1], slice2[1])

        # Last two items
        slice3 = sc[-2:]
        self.assertEqual(sc[1], slice3[0])
        self.assertEqual(sc[2], slice3[1])

        # First and third items
        slice4 = sc[0::2]
        self.assertEqual(sc[0], slice4[0])
        self.assertEqual(sc[2], slice4[1])

        # All items, reverse order
        slice5 = sc[::-1]
        for i in range(slice5.GetEntries()):
            self.assertEqual(sc[i], slice5[self.num_elems - 1 - i])

        # First and third items, reverse order
        slice6 = sc[::-2]
        self.assertEqual(sc[0], slice6[1])
        self.assertEqual(sc[2], slice6[0])

        # Step cannot be zero
        with self.assertRaises(ValueError):
            sc[::0]

    def test_setitem(self):
        sc = self.create_tseqcollection()
        l = []

        # Set items
        for i in range(self.num_elems):
            o = ROOT.TObject()
            sc[i] = o
            l.append(o)

        # Check previously set items
        it = ROOT.TIter(sc)
        for i in range(self.num_elems):
            self.assertEqual(it.Next(), l[i])

        # Set items, negative indices
        l2 = []
        neg_idcs = [ -i-1 for i in reversed(range(self.num_elems)) ]
        for i in neg_idcs:
            o = ROOT.TObject()
            sc[i] = o
            l2.append(o)

        # Check previously set items
        it2 = ROOT.TIter(sc)
        for i in range(self.num_elems):
            self.assertEqual(it2.Next(), l2[i])

        # Check invalid index cases
        with self.assertRaises(IndexError):
            sc[self.num_elems] = ROOT.TObject()

        with self.assertRaises(IndexError):
            sc[-(self.num_elems + 1)] = ROOT.TObject()

        with self.assertRaises(TypeError):
            sc[1.0] = ROOT.TObject()

    def test_setitem_slice(self):
        sc1 = self.create_tseqcollection()
        sc2 = self.create_tseqcollection()

        # Replace all items
        sc1[:] = sc2
        self.assertEqual(sc1.GetEntries(), self.num_elems)
        for i in range(self.num_elems):
            self.assertEqual(sc1[i], sc2[i])

        # Append items
        sc1 = self.create_tseqcollection()
        l = [ elem for elem in sc1 ]

        sc1[self.num_elems:] = sc2

        self.assertEqual(sc1.GetEntries(), 2 * self.num_elems)
        i = 0
        for elem in l: # first half
            self.assertEqual(sc1[i], elem)
            i += 1
        for elem in sc2: # second half
            self.assertEqual(sc1[i], elem)
            i += 1

        # Assign second item.
        # This time use a Python list as assigned value
        sc3 = self.create_tseqcollection()
        l2 = [ ROOT.TObject() ]
        l3 = [ elem for elem in sc3 ]

        sc3[1:2] = l2

        self.assertEqual(sc3.GetEntries(), self.num_elems)
        self.assertEqual(sc3[0], l3[0])
        self.assertEqual(sc3[1], l2[0])
        self.assertEqual(sc3[2], l3[2])

        # Assign second and third items to just one item.
        # This tests that the third item is removed
        sc4 = self.create_tseqcollection()
        l4 = [ ROOT.TObject() ]
        l5 = [ elem for elem in sc4 ]

        sc4[1:3] = l4

        self.assertEqual(sc4.GetEntries(), self.num_elems - 1)
        self.assertEqual(sc4[0], l5[0])
        self.assertEqual(sc4[1], l4[0])

        # Assign with step
        sc5 = self.create_tseqcollection()
        o = sc5[1]
        len6 = 2
        l6 = [ ROOT.TObject() for _ in range(len6) ]

        sc5[::2] = l6

        self.assertEqual(sc5.GetEntries(), self.num_elems)
        self.assertEqual(sc5[0], l6[0])
        self.assertEqual(sc5[1], o)
        self.assertEqual(sc5[2], l6[1])

        # Assign with step (start from end)
        sc5[::-2] = l6

        self.assertEqual(sc5.GetEntries(), self.num_elems)
        self.assertEqual(sc5[0], l6[1])
        self.assertEqual(sc5[1], o)
        self.assertEqual(sc5[2], l6[0])

        # Step cannot be zero
        sc6 = self.create_tseqcollection()
        with self.assertRaises(ValueError):
            sc6[::0] = [ ROOT.TObject() ]

    def test_delitem(self):
        sc = self.create_tseqcollection()

        self.assertEqual(sc.GetEntries(), self.num_elems)

        # Delete all elements
        for _ in range(self.num_elems):
            del sc[0]

        self.assertEqual(sc.GetEntries(), 0)

        sc = ROOT.TList()
        o1 = ROOT.TObject()
        o2 = ROOT.TObject()
        sc.Add(o1)
        sc.Add(o2)
        sc.Add(o1)

        # Delete o2
        del sc[1]

        # Only o1s should be there
        self.assertEqual(sc.GetEntries(), 2)

        it = ROOT.TIter(sc)
        for _ in range(2):
            self.assertEqual(it.Next(), o1)

        # Check invalid index cases
        with self.assertRaises(IndexError):
            del sc[2]

        with self.assertRaises(IndexError):
            del sc[-3]

        with self.assertRaises(TypeError):
            del sc[1.0]

    def test_delitem_slice(self):
        # Delete all items
        sc1 = self.create_tseqcollection()
        del sc1[:]
        self.assertEqual(sc1.GetEntries(), 0)

        # Do not delete anything (slice out of range)
        sc2 = self.create_tseqcollection()
        l2 = [ elem for elem in sc2 ]
        del sc2[self.num_elems:]
        self.assertEqual(sc2.GetEntries(), self.num_elems)
        for el1, el2 in zip(sc2, l2):
            self.assertEqual(el1, el2)

        # Delete first two items
        sc3 = self.create_tseqcollection()
        o = sc3[2]
        del sc3[0:2]
        self.assertEqual(sc3.GetEntries(), 1)
        self.assertEqual(sc3[0], o)

        # Delete first and third items
        sc4 = self.create_tseqcollection()
        o = sc4[1]
        del sc4[::2]
        self.assertEqual(sc4.GetEntries(), 1)
        self.assertEqual(sc4[0], o)

        # Delete first and third items (start from end)
        sc5 = self.create_tseqcollection()
        o = sc5[1]
        del sc5[::-2]
        self.assertEqual(sc5.GetEntries(), 1)
        self.assertEqual(sc5[0], o)

        # Step cannot be zero
        sc6 = self.create_tseqcollection()
        with self.assertRaises(ValueError):
            sc6[::0]


if __name__ == '__main__':
    unittest.main()
