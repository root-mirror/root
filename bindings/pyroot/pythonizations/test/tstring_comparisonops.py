import unittest

from ROOT import TString

class TStringComparisonOps(unittest.TestCase):
    """
    Test for the comparison operators of TString:
    __eq__, __ne__, __lt__, __le__, __gt__, __ge__.
    """

    num_elems = 3
    test_str1 = 'test1'
    test_str2 = 'test2'

    # Tests
    def test_eq(self):
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str1)
        ts3 = TString(self.test_str2)

        # Comparison between TStrings
        self.assertTrue(ts1 == ts2)
        self.assertFalse(ts1 == ts3)

        # Comparison with Python string
        self.assertTrue(ts1 == self.test_str1)
        self.assertFalse(ts1 == self.test_str2)

        # Comparison with non-string
        self.assertFalse(ts1 == 1)

    def test_ne(self):
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str1)
        ts3 = TString(self.test_str2)

        # Comparison between TStrings
        self.assertFalse(ts1 != ts2)
        self.assertTrue(ts1 != ts3)

        # Comparison with Python string
        self.assertFalse(ts1 != self.test_str1)
        self.assertTrue(ts1 != self.test_str2)

        # Comparison with non-string
        self.assertTrue(ts1 != 1)

    def test_lt(self):
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str2)

        # Comparison between TStrings
        self.assertTrue(ts1 < ts2)
        self.assertFalse(ts2 < ts1)

        # Comparison with Python string
        self.assertTrue(ts1 < self.test_str2)
        self.assertFalse(ts2 < self.test_str1)

    def test_le(self):
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str1)
        ts3 = TString(self.test_str2)

        # Comparison between TStrings
        self.assertTrue(ts1 <= ts2)
        self.assertTrue(ts1 <= ts3)
        self.assertFalse(ts3 <= ts1)

        # Comparison with Python string
        self.assertTrue(ts1 <= self.test_str1)
        self.assertTrue(ts1 <= self.test_str2)
        self.assertFalse(ts3 <= self.test_str1)

    def test_gt(self):
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str2)

        # Comparison between TStrings
        self.assertFalse(ts1 > ts2)
        self.assertTrue(ts2 > ts1)

        # Comparison with Python string
        self.assertFalse(ts1 > self.test_str2)
        self.assertTrue(ts2 > self.test_str1)

    def test_ge(self):
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str1)
        ts3 = TString(self.test_str2)

        # Comparison between TStrings
        self.assertTrue(ts1 >= ts2)
        self.assertFalse(ts1 >= ts3)
        self.assertTrue(ts3 >= ts1)

        # Comparison with Python string
        self.assertTrue(ts1 >= self.test_str1)
        self.assertFalse(ts1 >= self.test_str2)
        self.assertTrue(ts3 >= self.test_str1)

    def test_list_sort(self):
        l1 = [ TString(str(i)) for i in range(self.num_elems) ]
        l2 = list(reversed(l1))

        self.assertNotEqual(l1, l2)

        # Test that comparison operators enable list sorting
        l2.sort()

        self.assertEqual(l1, l2)


if __name__ == '__main__':
    unittest.main()
