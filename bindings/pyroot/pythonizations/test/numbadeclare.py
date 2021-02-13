import unittest
import ROOT
import sys
import os
import numpy as np
import gc


# Check whether these tests should be skipped
skip = False
skip_reason = ""
if sys.version_info[:3] <= (2, 7, 5):
    skip = True
    skip_reason = "Python version <= 2.7.5"
elif sys.version_info[0] == 2 and "ROOTTEST_IGNORE_NUMBA_PY2" in os.environ:
    skip = True
    skip_reason = "Running python2 and ROOTTEST_IGNORE_NUMBA_PY2 was set"
elif sys.version_info[0] == 3 and "ROOTTEST_IGNORE_NUMBA_PY3" in os.environ:
    skip = True
    skip_reason = "Running python3 and ROOTTEST_IGNORE_NUMBA_PY3 was set"

if not skip:
    import numba as nb


# long does not exist anymore on Python 3, map it to int
if sys.version_info[0] > 2:
    long = int


default_test_inputs = [-1.0, 0.0, 100.0]


class NumbaDeclareSimple(unittest.TestCase):
    """
    Test decorator to create C++ wrapper for Python callables using numba with fundamental types
    """

    test_inputs = default_test_inputs


    # Test refcounts
    @unittest.skipIf(skip, skip_reason)
    def test_refcount_decorator(self):
        """
        Test refcount of decorator
        """
        x = ROOT.Numba.Declare(["float"], "float")
        gc.collect()
        self.assertEqual(sys.getrefcount(x), 2)

    @unittest.skipIf(skip, skip_reason)
    def test_refcount_pycallable(self):
        """
        Test refcount of decorated callable
        """
        def f1(x):
            return x
        def f2(x):
            return x
        fn0 = ROOT.Numba.Declare(["float"], "float")(f1)
        ref = nb.cfunc("float32(float32)", nopython=True)(f2)
        gc.collect()
        if sys.version_info.major == 2:
            self.assertEqual(sys.getrefcount(f1), sys.getrefcount(f2) + 1)
        else:
            self.assertEqual(sys.getrefcount(f1), sys.getrefcount(f2) + 2)

    # Test optional name
    @unittest.skipIf(skip, skip_reason)
    def test_optional_name(self):
        """
        Test optional name of wrapper function
        """
        optname = "optname2"
        @ROOT.Numba.Declare(["float"], "float", name=optname)
        def f(x):
            return x
        self.assertTrue(hasattr(ROOT.Numba, optname))

    # Test attributes
    @unittest.skipIf(skip, skip_reason)
    def test_additional_attributes(self):
        """
        Test additional attributes
        """
        @ROOT.Numba.Declare(["float"], "float")
        def fn1(x):
            return x

        gc.collect()

        self.assertTrue(hasattr(fn1, "__cpp_wrapper__"))
        self.assertTrue(type(fn1.__cpp_wrapper__) == str)
        self.assertEqual(sys.getrefcount(fn1.__cpp_wrapper__), 2)

        self.assertTrue(hasattr(fn1, "__py_wrapper__"))
        self.assertTrue(type(fn1.__py_wrapper__) == str)
        self.assertEqual(sys.getrefcount(fn1.__py_wrapper__), 2)

        self.assertTrue(hasattr(fn1, "numba_func"))
        self.assertEqual(sys.getrefcount(fn1.numba_func), 3)

    # Test cling integration
    @unittest.skipIf(skip, skip_reason)
    def test_cling(self):
        """
        Test function call in cling
        """
        @ROOT.Numba.Declare(["float"], "float")
        def fn12(x):
            return 2.0 * x
        ROOT.gInterpreter.ProcessLine("y12 = Numba::fn12(42.0);")
        self.assertEqual(fn12(42.0), ROOT.y12)

    # Test RDataFrame integration
    @unittest.skipIf(skip, skip_reason)
    def test_rdataframe(self):
        """
        Test function call as part of RDataFrame
        """
        @ROOT.Numba.Declare(["unsigned int"], "float")
        def fn13(x):
            return 2.0 * x
        df = ROOT.RDataFrame(4).Define("x", "rdfentry_").Define("y", "Numba::fn13(x)")
        mean_x = df.Mean("x")
        mean_y = df.Mean("y")
        self.assertEqual(mean_x.GetValue(), 1.5)
        self.assertEqual(mean_y.GetValue(), 3.0)

    # Test wrappings
    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_void(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare([], "float")
        def fn2n():
            return float(42)
        x1 = fn2n()
        x2 = ROOT.Numba.fn2n()
        self.assertEqual(x1, x2)
        self.assertEqual(type(x1), type(x2))

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_f(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["float"], "float")
        def fn2(x):
            return float(x)
        for v in self.test_inputs:
            x1 = fn2(v)
            x2 = ROOT.Numba.fn2(v)
            self.assertEqual(x1, x2)
            self.assertEqual(type(x1), type(x2))

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_d(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["float"], "double")
        def fn2d(x):
            return float(x)
        for v in self.test_inputs:
            x1 = fn2d(v)
            x2 = ROOT.Numba.fn2d(v)
            self.assertEqual(x1, x2)
            # NOTE: There is no double in Python because everything is a double.
            self.assertEqual(type(x1), type(x2))

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_i(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["float"], "int")
        def fn3(x):
            return int(x)
        for v in self.test_inputs:
            x1 = fn3(v)
            x2 = ROOT.Numba.fn3(v)
            self.assertEqual(x1, x2)
            self.assertEqual(type(x1), type(x2))

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_l(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["float"], "long")
        def fn4(x):
            return x
        for v in self.test_inputs:
            x1 = fn4(v)
            x2 = ROOT.Numba.fn4(v)
            self.assertEqual(x1, x2)
            self.assertEqual(long, type(x2))

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_u(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["float"], "unsigned int")
        def fn5(x):
            return abs(x)
        for v in self.test_inputs:
            x1 = fn5(v)
            x2 = ROOT.Numba.fn5(v)
            self.assertEqual(x1, x2)
            # NOTE: cppyy does not return an Python int for unsigned int but a long.
            # This could be fixed but as well should not have any impact.
            self.assertEqual(type(x2), long)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_k(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["float"], "unsigned long")
        def fn6(x):
            return abs(x)
        for v in self.test_inputs:
            x1 = fn6(v)
            x2 = ROOT.Numba.fn6(v)
            self.assertEqual(x1, x2)
            self.assertEqual(long, type(x2))

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_b(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["float"], "bool")
        def fn6b(x):
            return x>0
        for v in self.test_inputs:
            x1 = fn6b(v)
            x2 = ROOT.Numba.fn6b(v)
            self.assertEqual(x1, x2)
            self.assertEqual(type(x1), type(x2))

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_b(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["bool"], "bool")
        def fn6b2(x):
            return not x
        for v in [True, False]:
            x1 = fn6b2(v)
            x2 = ROOT.Numba.fn6b2(v)
            self.assertEqual(x1, x2)
            self.assertEqual(type(x1), type(x2))

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_i(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["int"], "float")
        def fn7i(x):
            return 2.0 * x
        for v in [-1, 0, 1, 999]:
            x1 = fn7i(v)
            x2 = ROOT.Numba.fn7i(v)
            self.assertEqual(x1, x2)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_l(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["long"], "float")
        def fn7l(x):
            return 2.0 * x
        for v in [-1, 0, 1, 999]:
            x1 = fn7l(v)
            x2 = ROOT.Numba.fn7l(v)
            self.assertEqual(x1, x2)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_ui(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["unsigned int"], "float")
        def fn7ui(x):
            return 2.0 * x
        for v in [0, 1, 999]:
            x1 = fn7ui(v)
            x2 = ROOT.Numba.fn7ui(v)
            self.assertEqual(x1, x2)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_ul(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["unsigned long"], "float")
        def fn7ul(x):
            return 2.0 * x
        for v in [0, 1, 999]:
            x1 = fn7ul(v)
            x2 = ROOT.Numba.fn7ul(v)
            self.assertEqual(x1, x2)


class NumbaDeclareArray(unittest.TestCase):
    """
    Test decorator to create C++ wrapper for Python callables using numba with RVecs
    """

    test_inputs = [default_test_inputs]

    # The global module index does not have RVec entities preloaded and
    # gInterpreter.Declare is not allowed to load libROOTVecOps for RVec.
    # Preload the library now.
    ROOT.gSystem.Load("libROOTVecOps")

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_vecf(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<float>"], "float")
        def g1(x):
            return x.sum()

        for v in self.test_inputs:
            x1 = g1(np.array(v, dtype=np.float32))
            x2 = ROOT.Numba.g1(ROOT.VecOps.RVec('float')(v))
            self.assertEqual(x1, x2)
            self.assertEqual(type(x2), float)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_vecf_vecd(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<float>", "RVec<double>"], "float")
        def g1_2vec(x, y):
            return x.sum() + y.sum()

        for v in self.test_inputs:
            x1 = g1_2vec(np.array(v, dtype=np.float32), np.array(v, dtype=np.float64))
            x2 = ROOT.Numba.g1_2vec(ROOT.VecOps.RVec('float')(v), ROOT.VecOps.RVec('double')(v))
            self.assertEqual(x1, x2)
            self.assertEqual(type(x2), float)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_vecd(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<double>"], "float")
        def g1d(x):
            return x.sum()

        for v in self.test_inputs:
            x1 = g1d(np.array(v, dtype=np.float64))
            x2 = ROOT.Numba.g1d(ROOT.VecOps.RVec('double')(v))
            self.assertEqual(x1, x2)
            self.assertEqual(type(x2), float)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_veci(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<int>"], "int")
        def g1i(x):
            return x.sum()

        for v in self.test_inputs:
            x1 = g1i(np.array(v, dtype=np.int32))
            x2 = ROOT.Numba.g1i(ROOT.VecOps.RVec('int')(int(x) for x in v))
            self.assertEqual(x1, x2)
            self.assertEqual(type(x2), int)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_vecl(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<long>"], "int")
        def g1l(x):
            return x.sum()

        for v in self.test_inputs:
            x1 = g1l(np.array(v, dtype=np.int64))
            x2 = ROOT.Numba.g1l(ROOT.VecOps.RVec('long')(int(x) for x in v))
            self.assertEqual(x1, x2)
            self.assertEqual(type(x2), int)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_vecui(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<unsigned int>"], "int")
        def g1ui(x):
            return x.sum()

        for v in [[0, 1, 999]]:
            x1 = g1ui(np.array(v, dtype=np.uint32))
            x2 = ROOT.Numba.g1ui(ROOT.VecOps.RVec('unsigned int')(v))
            self.assertEqual(x1, x2)
            self.assertEqual(type(x2), int)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_vecul(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<unsigned long>"], "int")
        def g1ul(x):
            return x.sum()

        for v in [[0, 1, 999]]:
            x1 = g1ul(np.array(v, dtype=np.uint64))
            x2 = ROOT.Numba.g1ul(ROOT.VecOps.RVec('unsigned long')(v))
            self.assertEqual(x1, x2)
            self.assertEqual(type(x2), int)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_vecb(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<bool>"], "int")
        def g1b(x):
            return x.sum()

        for v in [[True, False, True]]:
            x1 = g1b(np.array(v, dtype=np.float32))
            x2 = ROOT.Numba.g1b(ROOT.VecOps.RVec('bool')(v))
            self.assertEqual(x1, x2)
            self.assertEqual(type(x2), int)

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_vecf(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<float>"], "RVec<float>")
        def g2f(x):
            return x[::-1]

        for v in [[0, 1, 999]]:
            x1 = g2f(np.array(v, dtype=np.float32))
            x2 = ROOT.Numba.g2f(ROOT.VecOps.RVec('float')(v))
            self.assertTrue((x1 == x2).all())

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_vecd(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<double>"], "RVec<double>")
        def g2d(x):
            return x[::-1]

        for v in [[0, 1, 999]]:
            x1 = g2d(np.array(v, dtype=np.float64))
            x2 = ROOT.Numba.g2d(ROOT.VecOps.RVec('double')(v))
            self.assertTrue((x1 == x2).all())

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_veci(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<int>"], "RVec<int>")
        def g2i(x):
            return x[::-1]

        for v in [[0, 1, 999]]:
            x1 = g2i(np.array(v, dtype=np.int32))
            x2 = ROOT.Numba.g2i(ROOT.VecOps.RVec('int')(v))
            self.assertTrue((x1 == x2).all())

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_vecl(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<long>"], "RVec<long>")
        def g2l(x):
            return x[::-1]

        for v in [[0, 1, 999]]:
            x1 = g2l(np.array(v, dtype=np.int64))
            x2 = ROOT.Numba.g2l(ROOT.VecOps.RVec('long')(v))
            self.assertTrue((x1 == x2).all())

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_vecul(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<unsigned long>"], "RVec<unsigned long>")
        def g2ul(x):
            return x[::-1]

        for v in [[0, 1, 999]]:
            x1 = g2ul(np.array(v, dtype=np.uint64))
            x2 = ROOT.Numba.g2ul(ROOT.VecOps.RVec('unsigned long')(v))
            self.assertTrue((x1 == x2).all())

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_vecui(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<unsigned int>"], "RVec<unsigned int>")
        def g2ui(x):
            return x[::-1]

        for v in [[0, 1, 999]]:
            x1 = g2ui(np.array(v, dtype=np.uint32))
            x2 = ROOT.Numba.g2ui(ROOT.VecOps.RVec('unsigned int')(v))
            self.assertTrue((x1 == x2).all())

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_out_vecb(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<bool>"], "RVec<bool>")
        def g2b(x):
            return x[::-1]

        for v in [[True, False]]:
            x1 = g2b(np.array(v, dtype=np.bool))
            x2 = ROOT.Numba.g2b(ROOT.VecOps.RVec('bool')(v))
            self.assertEqual(x1[0], bool(x2[0]))
            self.assertEqual(x1[1], bool(x2[1]))

    @unittest.skipIf(skip, skip_reason)
    def test_wrapper_in_vecfb_out_vecf(self):
        """
        Test wrapper with different input/output configurations
        """
        @ROOT.Numba.Declare(["RVec<float>", "RVec<bool>"], "RVec<bool>")
        def g2fb(x, y):
            return (x > 1) | y

        for vf, vb in [[[1.0, 2.0], [True, False]]]:
            x1 = g2fb(np.array(vf, dtype=np.float32), np.array(vb, dtype=np.bool))
            x2 = ROOT.Numba.g2fb(ROOT.VecOps.RVec('float')(vf), ROOT.VecOps.RVec('bool')(vb))
            self.assertEqual(x1[0], bool(x2[0]))
            self.assertEqual(x1[1], bool(x2[1]))


if __name__ == '__main__':
    unittest.main()
