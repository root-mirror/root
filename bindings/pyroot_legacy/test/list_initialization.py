import unittest
import ROOT

cppcode = """
struct A {
  A(int _i, double _f) : i(_i), f(_f) {}
  A(int _i) : i(_i), f(99.99) {}
  int i;
  double f;
};
A func_by_ref(const A& a) { return a; }
A func_by_val(A a) { return a; }

"""

class ListInitialization(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ROOT.gInterpreter.Declare(cppcode)

    def test_function_call_by_ref(self):
        self.assertEqual(ROOT.func_by_ref((10,20.)).i, 10)
        self.assertEqual(ROOT.func_by_ref((10,20.)).f, 20.)
        self.assertAlmostEqual(ROOT.func_by_ref((10,)).f, 99.99)
    
    def test_invalid_constructor(self):
        with self.assertRaises(TypeError): ROOT.func_by_ref((1,2,3))
        with self.assertRaises(TypeError): ROOT.func_by_ref((1,'abc'))

    def test_function_call_by_val(self):
        self.assertEqual(ROOT.func_by_val((10,20.)).i, 10)
        self.assertEqual(ROOT.func_by_val((10,20.)).f, 20.)
        self.assertAlmostEqual(ROOT.func_by_val((10,)).f, 99.99)

if __name__ == '__main__':
    unittest.main()
