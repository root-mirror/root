import unittest

import ROOT


class TObjStringStrRepr(unittest.TestCase):
    """
    Test for the pythonizations that provide a string representation
    for instances of TObjString (__str__, __repr__).
    """

    # Tests
    def test_str(self):
        s = 'test'
        tos = ROOT.TObjString(s)
        self.assertEqual(str(tos), s)

    def test_repr(self):
        s = 'test'
        tos = ROOT.TObjString(s)
        self.assertEqual(repr(tos), repr(s))


if __name__ == '__main__':
    unittest.main()
