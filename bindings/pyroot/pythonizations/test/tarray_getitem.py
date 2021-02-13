import unittest

import ROOT


class TArrayGetItem(unittest.TestCase):
    """
    Test for the pythonization that allows to: (i) get an item of a
    TArray with boundary check for the index and (ii) iterate over
    a TArray.
    """

    num_elems = 3

    # Tests
    def test_boundary_check(self):
        a = ROOT.TArrayI(self.num_elems)

        # In range
        self.assertEqual(a[0], a[0])

        # Out of range
        with self.assertRaises(IndexError):
            a[-1]

        # Out of range
        with self.assertRaises(IndexError):
            a[self.num_elems]

    def test_iterable(self):
        a = ROOT.TArrayI(self.num_elems)
        val = 1

        for i in range(self.num_elems):
            a[i] = val

        self.assertEquals(list(a), [ val for _ in range(self.num_elems) ])


if __name__ == '__main__':
    unittest.main()
