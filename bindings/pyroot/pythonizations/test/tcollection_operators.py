import unittest

import ROOT
from libcppyy import SetOwnership


class TCollectionOperators(unittest.TestCase):
    """
    Test for the Python operators defined in TCollection (and subclasses):
    __add__, __mul__, __rmul__, __imul__
    """

    num_elems = 3
    factor = 2

    # Helpers
    def create_tcollection(self):
        c = ROOT.TList()
        for _ in range(self.num_elems):
            o = ROOT.TObject()
            # Prevent immediate deletion of C++ TObjects
            SetOwnership(o, False)
            c.Add(o)

        return c

    def check_mul_result(self, c, cmul):
        lenc = c.GetEntries()

        self.assertEqual(cmul.GetEntries(), lenc * self.factor)

        itmul = ROOT.TIter(cmul)
        for _ in range(self.factor):
            itc = ROOT.TIter(c)
            for _ in range(lenc):
                oc = itc.Next()
                omul = itmul.Next()
                self.assertEqual(oc, omul)

    # Tests
    def test_add(self):
        c1 = self.create_tcollection()
        c2 = self.create_tcollection()

        len1 = c1.GetEntries()
        len2 = c2.GetEntries()

        cadd = c1 + c2

        len_add = cadd.GetEntries()

        self.assertEqual(len_add, len1 + len2)

        # Compare with elements of first collection
        itc1 = ROOT.TIter(c1)
        itadd = ROOT.TIter(cadd)
        for _ in range(len1):
            oc1 = itc1.Next()
            oadd = itadd.Next()
            self.assertEqual(oc1, oadd)

        # Compare with elements of second collection
        itc2 = ROOT.TIter(c2)
        for _ in range(len2):
            oc2 = itc2.Next()
            oadd = itadd.Next()
            self.assertEqual(oc2, oadd)

    def test_mul(self):
        c = self.create_tcollection()

        cmul = c * self.factor

        self.check_mul_result(c, cmul)

    def test_rmul(self):
        c = self.create_tcollection()

        cmul = self.factor * c

        self.check_mul_result(c, cmul)

    def test_imul(self):
        c = self.create_tcollection()
        lenc = c.GetEntries()

        c *= self.factor

        self.assertEqual(c.GetEntries(), lenc * self.factor)

        it = ROOT.TIter(c)
        subc = []
        for _ in range(lenc):
            subc.append(it.Next())

        for _ in range(self.factor - 1):
            for o in subc:
                self.assertEqual(o, it.Next())


if __name__ == '__main__':
    unittest.main()

