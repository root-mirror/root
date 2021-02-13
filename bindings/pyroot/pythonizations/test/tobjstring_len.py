import unittest

import ROOT


class TObjStringLen(unittest.TestCase):
    """
    Test for the pythonization that provides the length of a
    TObjString instance `s` via `len(s)`.
    """

    # Tests
    def test_len(self):
        s = 'test'
        tos = ROOT.TObjString(s)
        self.assertEqual(len(tos), len(s))


if __name__ == '__main__':
    unittest.main()
