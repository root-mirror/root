import unittest

import ROOT
from libcppyy import SetOwnership


class TSeqCollectionListMethods(unittest.TestCase):
    """
    Test for the Python-list-like methods added to TSeqCollection
    (and subclasses): insert, pop, reverse, sort, index
    """

    num_elems = 3

    # Helpers
    def create_tseqcollection(self):
        sc = ROOT.TList()
        for i in reversed(range(self.num_elems)):
            sc.Add(ROOT.TObjString(str(i)))

        return sc

    # Tests
    def test_insert(self):
        sc = self.create_tseqcollection()

        # Insert with positive index
        o = ROOT.TObject()
        sc.insert(1, o)
        self.assertEqual(sc.GetEntries(), self.num_elems + 1)
        self.assertEqual(sc.At(1), o)

        # Insert with negative index (starts from end)
        o2 = ROOT.TObject()
        sc.insert(-1, o2)
        self.assertEqual(sc.GetEntries(), self.num_elems + 2)
        self.assertEqual(sc.At(self.num_elems), o2)

        # Insert with index beyond lower boundary.
        # Inserts at the beginning
        o3 = ROOT.TObject()
        sc.insert(-(self.num_elems + 3), o3)
        self.assertEqual(sc.GetEntries(), self.num_elems + 3)
        self.assertEqual(sc.At(0), o3)

        # Insert with index beyond upper boundary.
        # Inserts at the end
        o4 = ROOT.TObject()
        sc.insert(self.num_elems + 4, o4)
        self.assertEqual(sc.GetEntries(), self.num_elems + 4)
        self.assertEqual(sc.At(self.num_elems + 3), o4)

    def test_pop(self):
        sc = self.create_tseqcollection()
        l = [ elem for elem in sc ]

        # No arguments, pop last item
        self.assertEqual(sc.pop(), l[-1])
        self.assertEqual(sc.GetEntries(), self.num_elems - 1)

        # Pop first item, positive index
        self.assertEqual(sc.pop(0), l[0])
        self.assertEqual(sc.GetEntries(), self.num_elems - 2)

        # Pop last item, negative index
        self.assertEqual(sc.pop(-1), l[1])
        self.assertEqual(sc.GetEntries(), self.num_elems - 3)

        # Pop from empty collection
        with self.assertRaises(IndexError):
            sc.pop()

        # Index out of range, positive
        sc2 = self.create_tseqcollection()
        with self.assertRaises(IndexError):
            sc2.pop(self.num_elems)

        # Index out of range, negative
        with self.assertRaises(IndexError):
            sc2.pop(-(self.num_elems + 1))

        # Pop with non-integer argument
        with self.assertRaises(TypeError):
            sc2.pop(1.0)

        # Pop a repeated element
        sc2.append(ROOT.TObjString('2'))
        elem = sc2.pop()
        self.assertEqual(sc2.At(0), elem)

    def test_reverse(self):
        sc = self.create_tseqcollection()
        l = [ elem for elem in sc ]

        sc.reverse()

        self.assertEqual(sc.GetEntries(), self.num_elems)
        for i,elem in zip(range(self.num_elems), sc):
            self.assertEqual(elem, l[-(i+1)])

        # Empty collection
        sc2 = ROOT.TList()
        sc2.reverse()
        self.assertEqual(sc2.GetEntries(), 0)

    def test_sort(self):
        sc = self.create_tseqcollection()
        l = [ elem for elem in sc ]

        # Regular sort, rely on TList::Sort
        sc.sort()
        # We need to set `key` until the pythonization to
        # make TObjString comparable is there
        l.sort(key = lambda s: s.GetName())

        self.assertEqual(sc.GetEntries(), self.num_elems)
        self.assertEqual(l[0], sc[0])
        for el1, el2 in zip(sc, l):
            self.assertEqual(el1, el2)

        # Python sort, key and reverse arguments.
        # Sort by hash in reverse order
        sc2 = self.create_tseqcollection()
        l2 = [ elem for elem in sc2 ]

        fsort = lambda elem: elem.Hash()
        rev = True
        sc2.sort(key = fsort, reverse = rev)
        l2.sort(key = fsort, reverse = rev)

        self.assertEqual(sc2.GetEntries(), self.num_elems)
        for el1, el2 in zip(sc2, l2):
            self.assertEqual(el1, el2)

        # Empty collection
        sc4 = ROOT.TList()
        sc4.sort()
        self.assertEqual(sc4.GetEntries(), 0)

    def test_index(self):
        sc = self.create_tseqcollection()

        # Check all elements of collection
        for i, elem in zip(range(self.num_elems), sc):
            self.assertEqual(sc.index(elem), i)

        # Check element not in collection
        o = ROOT.TObjString(str(self.num_elems))
        with self.assertRaises(ValueError):
            sc.index(o)


if __name__ == '__main__':
    unittest.main()
