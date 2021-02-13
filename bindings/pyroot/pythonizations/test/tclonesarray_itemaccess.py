import unittest

import ROOT


class TClonesArrayItemAccess(unittest.TestCase):
    """
    Test for the item access method added to TClonesArray:
    __setitem__.
    """

    num_elems = 3

    # Helpers
    def create_tclonesarray_and_list(self):
        ca = ROOT.TClonesArray("TObject", self.num_elems)
        l = [ ROOT.TObject() for _ in range(self.num_elems) ]
        return ca, l

    # Tests
    def test_setitem(self):
        ca, l = self.create_tclonesarray_and_list()
        for item, i in zip(l, range(self.num_elems)):
            ca[i] = item
            self.assertEqual(ca[i], item)


if __name__ == '__main__':
    unittest.main()
