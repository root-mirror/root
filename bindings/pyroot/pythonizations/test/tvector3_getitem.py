import unittest

import ROOT


class TVector3GetItem(unittest.TestCase):
    """
    Test for the pythonization that allows to: (i) get an item of a
    TVector3 with boundary check for the index and (ii) iterate over
    a TVector3.
    """

    # Tests
    def test_boundary_check(self):
        v = ROOT.TVector3(1., 2., 3.)

        # In range
        self.assertEqual(v[0], v[0])

        # Out of range
        with self.assertRaises(IndexError):
            v[-1]

        # Out of range
        with self.assertRaises(IndexError):
            v[3]

    def test_iterable(self):
        v = ROOT.TVector3(1., 2., 3.)

        self.assertEquals(list(v), [1., 2., 3.])


if __name__ == '__main__':
    unittest.main()
