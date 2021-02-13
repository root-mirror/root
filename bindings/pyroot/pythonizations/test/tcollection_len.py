import unittest

import ROOT


class TCollectionLen(unittest.TestCase):
    """
    Test for the pythonization that allows to access the number of elements of a
    TCollection (or subclass) by calling `len` on it.
    """

    num_elems = 3
    tobject_list = [ ROOT.TObject() for _ in range(num_elems) ]

    # Helpers
    def add_elems_check_len(self, c):
        for elem in self.tobject_list:
            c.Add(elem)

        self.assertEqual(len(c), self.num_elems)
        self.assertEqual(len(c), c.GetEntries())

    # Tests
    def test_tlist(self):
        self.add_elems_check_len(ROOT.TList())

    def test_tobjarray(self):
        self.add_elems_check_len(ROOT.TObjArray())

    def test_thashtable(self):
        self.add_elems_check_len(ROOT.THashTable())


if __name__ == '__main__':
    unittest.main()

