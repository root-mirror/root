import unittest

import ROOT


class TArrayLen(unittest.TestCase):
    """
    Test for the pythonization that allows to access the number of elements of a
    TArray (or subclass) by calling `len` on it.
    """

    num_elems = 3

    # Helpers
    def check_len(self, a):
        self.assertEqual(len(a), self.num_elems)
        self.assertEqual(len(a), a.GetSize())

    # Tests
    def test_tarrayi(self):
        self.check_len(ROOT.TArrayI(self.num_elems))

    def test_tarrayd(self):
        self.check_len(ROOT.TArrayD(self.num_elems))


if __name__ == '__main__':
    unittest.main()

