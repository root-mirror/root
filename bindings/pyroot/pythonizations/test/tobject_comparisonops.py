import unittest

import ROOT
from ROOT import TUrl


class TObjectComparisonOps(unittest.TestCase):
    """
    Test for the comparison operators of TObject and subclasses:
    __eq__, __ne__, __lt__, __le__, __gt__, __ge__.
    Such pythonisations rely on TObject::IsEqual and TObject::Compare,
    which can be reimplemented in subclasses.
    """

    num_elems = 3

    # Tests
    def test_eq(self):
        o = ROOT.TObject()

        # TObject::IsEqual compares internal address
        self.assertTrue(o == o)

        # Test comparison with no TObject
        self.assertFalse(o == 1)

        # Test comparison with None
        self.assertFalse(o == None)

    def test_ne(self):
        o = ROOT.TObject()

        # TObject::IsEqual compares internal address
        self.assertFalse(o != o)

        # Test comparison with no TObject
        self.assertTrue(o != 1)

        # Test comparison with None
        self.assertTrue(o != None)

    def test_lt(self):
        a = TUrl("a")
        b = TUrl("b")

        # TUrl::Compare compares URL strings
        self.assertTrue(a < b)
        self.assertFalse(b < a)

        # Test comparison with no TObject
        self.assertEqual(a.__lt__(1), NotImplemented)

    def test_le(self):
        a1 = TUrl("a")
        a2 = TUrl("a")
        b  = TUrl("b")

        # TUrl::Compare compares URL strings
        self.assertTrue(a1 <= a2)
        self.assertTrue(a2 <= a1)
        self.assertTrue(a1 <= b)
        self.assertFalse(b <= a1)

        # Test comparison with no TObject
        self.assertEqual(a1.__le__(1), NotImplemented)

    def test_gt(self):
        a = TUrl("a")
        b = TUrl("b")

        # TUrl::Compare compares URL strings
        self.assertFalse(a > b)
        self.assertTrue(b > a)

        # Test comparison with no TObject
        self.assertEqual(a.__gt__(1), NotImplemented)

    def test_ge(self):
        a1 = TUrl("a")
        a2 = TUrl("a")
        b  = TUrl("b")

        # TUrl::Compare compares URL strings
        self.assertTrue(a1 >= a2)
        self.assertTrue(a2 >= a1)
        self.assertTrue(b >= a1)
        self.assertFalse(a1 >= b)

        # Test comparison with no TObject
        self.assertEqual(a1.__ge__(1), NotImplemented)

    def test_list_sort(self):
        l1 = [ ROOT.TUrl(str(i)) for i in range(self.num_elems) ]
        l2 = list(reversed(l1))

        self.assertNotEqual(l1, l2)

        # Test that comparison operators enable list sorting
        l2.sort()

        self.assertEqual(l1, l2)


if __name__ == '__main__':
    unittest.main()
