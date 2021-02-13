import unittest

from ROOT import TObjString, TString

class TObjStringComparisonOps(unittest.TestCase):
    """
    Test for the comparison operators of TObjString:
    __eq__, __ne__, __lt__, __le__, __gt__, __ge__.
    """

    num_elems = 3
    test_str1 = 'test1'
    test_str2 = 'test2'

    # Tests
    def test_eq(self):
        tos1 = TObjString(self.test_str1)
        tos2 = TObjString(self.test_str1)
        tos3 = TObjString(self.test_str2)

        # Comparison between TObjStrings
        self.assertTrue(tos1 == tos2)
        self.assertFalse(tos1 == tos3)

        # Comparison with TString
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str2)
        self.assertTrue(tos1 == ts1)
        self.assertFalse(tos1 == ts2)

        # Comparison with Python string
        self.assertTrue(tos1 == self.test_str1)
        self.assertFalse(tos1 == self.test_str2)

        # Comparison with non-string
        self.assertFalse(tos1 == 1)

    def test_ne(self):
        tos1 = TObjString(self.test_str1)
        tos2 = TObjString(self.test_str1)
        tos3 = TObjString(self.test_str2)

        # Comparison between TObjStrings
        self.assertFalse(tos1 != tos2)
        self.assertTrue(tos1 != tos3)

        # Comparison with TString
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str2)
        self.assertFalse(tos1 != ts1)
        self.assertTrue(tos1 != ts2)

        # Comparison with Python string
        self.assertFalse(tos1 != self.test_str1)
        self.assertTrue(tos1 != self.test_str2)

        # Comparison with non-string
        self.assertTrue(tos1 != 1)

    def test_lt(self):
        tos1 = TObjString(self.test_str1)
        tos2 = TObjString(self.test_str2)

        # Comparison between TObjStrings
        self.assertTrue(tos1 < tos2)
        self.assertFalse(tos2 < tos1)

        # Comparison with TString
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str2)
        self.assertTrue(tos1 < ts2)
        self.assertFalse(tos2 < ts1)

        # Comparison with Python string
        self.assertTrue(tos1 < self.test_str2)
        self.assertFalse(tos2 < self.test_str1)

    def test_le(self):
        tos1 = TObjString(self.test_str1)
        tos2 = TObjString(self.test_str1)
        tos3 = TObjString(self.test_str2)

        # Comparison between TObjStrings
        self.assertTrue(tos1 <= tos2)
        self.assertTrue(tos1 <= tos3)
        self.assertFalse(tos3 <= tos1)

        # Comparison with TString
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str2)
        self.assertTrue(tos1 <= ts1)
        self.assertTrue(tos1 <= ts2)
        self.assertFalse(tos3 <= ts1)

        # Comparison with Python string
        self.assertTrue(tos1 <= self.test_str1)
        self.assertTrue(tos1 <= self.test_str2)
        self.assertFalse(tos3 <= self.test_str1)

    def test_gt(self):
        tos1 = TObjString(self.test_str1)
        tos2 = TObjString(self.test_str2)

        # Comparison between TObjStrings
        self.assertFalse(tos1 > tos2)
        self.assertTrue(tos2 > tos1)

        # Comparison with TString
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str2)
        self.assertFalse(tos1 > ts2)
        self.assertTrue(tos2 > ts1)

        # Comparison with Python string
        self.assertFalse(tos1 > self.test_str2)
        self.assertTrue(tos2 > self.test_str1)

    def test_ge(self):
        tos1 = TObjString(self.test_str1)
        tos2 = TObjString(self.test_str1)
        tos3 = TObjString(self.test_str2)

        # Comparison between TObjStrings
        self.assertTrue(tos1 >= tos2)
        self.assertFalse(tos1 >= tos3)
        self.assertTrue(tos3 >= tos1)

        # Comparison with TString
        ts1 = TString(self.test_str1)
        ts2 = TString(self.test_str2)
        self.assertTrue(tos1 >= ts1)
        self.assertFalse(tos1 >= ts2)
        self.assertTrue(tos3 >= ts1)

        # Comparison with Python string
        self.assertTrue(tos1 >= self.test_str1)
        self.assertFalse(tos1 >= self.test_str2)
        self.assertTrue(tos3 >= self.test_str1)

    def test_list_sort(self):
        l1 = [ TObjString(str(i)) for i in range(self.num_elems) ]
        l2 = list(reversed(l1))

        self.assertNotEqual(l1, l2)

        # Test that comparison operators enable list sorting
        l2.sort()

        self.assertEqual(l1, l2)


if __name__ == '__main__':
    unittest.main()
