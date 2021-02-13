import unittest
import ROOT
from ROOT import TComplex

class TestTComplexOperators(unittest.TestCase):
    """
    Test for the operators of TComplex:
    __radd__, __rsub__, __rmul__, __rtruediv__/__rdiv__.
    """

    c = TComplex(4.,0)

    d = 2.

    s = 'string'

    # check the expected result for d + c and that Re(c + d) == Re(d + c)
    def test_radd(self):
        self.assertEqual((self.d + self.c).Re(), 6.0)
        self.assertEqual((self.c + self.d).Re(), (self.d + self.c).Re())

    # check the expected result for d - c and that Re(c - d) == -Re(d - c)
    def test_rsub(self):
        self.assertEqual((self.d - self.c).Re(), -2.0)
        self.assertEqual((self.c - self.d).Re(), -((self.d - self.c).Re()))

    # check the expected result for d * c and that Re(c * d) == Re(d * c)
    def test_rmul(self):
        self.assertEqual((self.d * self.c).Re(), 8.0)
        self.assertEqual((self.c * self.d).Re(), (self.d * self.c).Re())

    # check the expected result for d / c
    def test_rdiv(self):
        self.assertEqual((self.d / self.c).Re(), 0.5)
        with self.assertRaises(TypeError):
            self.s / self.c


if __name__ == '__main__':
    unittest.main()
