import unittest

import ROOT


class TVectorTLen(unittest.TestCase):
    """
    Test for the pythonization that allows to get the size of a
    TVectorT by calling `len` on it.
    """

    num_elems = 3

    # Tests
    def test_len(self):
        v = ROOT.TVectorT[float](self.num_elems)
        self.assertEqual(len(v), self.num_elems)


if __name__ == '__main__':
    unittest.main()
