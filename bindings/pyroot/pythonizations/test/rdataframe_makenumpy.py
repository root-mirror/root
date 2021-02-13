import unittest
import ROOT
import numpy as np
import sys
import gc


class MakeNumpyDataFrame(unittest.TestCase):
    """
    Tests for the MakeNumpyDataFrame feature enabling to read numpy arrays
    with RDataFrame.
    """

    dtypes = [
        "int32", "int64", "uint32", "uint64", "float32", "float64"
    ]

    def test_dtypes(self):
        """
        Test reading different datatypes
        """
        for dtype in self.dtypes:
            data = {"x": np.array([1, 2, 3], dtype=dtype)}
            df = ROOT.RDF.MakeNumpyDataFrame(data)
            self.assertEqual(df.Mean("x").GetValue(), 2)

    def test_multiple_columns(self):
        """
        Test reading multiple columns
        """
        data = {}
        for dtype in self.dtypes:
            data[dtype] = np.array([1, 2, 3], dtype=dtype)
        df = ROOT.RDF.MakeNumpyDataFrame(data)
        colnames = df.GetColumnNames()
        # Test column names
        for dtype in colnames:
            self.assertIn(dtype, self.dtypes)
        # Test mean
        for dtype in self.dtypes:
            self.assertEqual(df.Mean(dtype).GetValue(), 2)

    def test_refcount(self):
        """
        Check refcounts of associated PyObjects
        """
        data = {"x": np.array([1, 2, 3], dtype="float32")}
        gc.collect()
        self.assertEqual(sys.getrefcount(data), 2)
        self.assertEqual(sys.getrefcount(data["x"]), 2)

        df = ROOT.RDF.MakeNumpyDataFrame(data)
        gc.collect()
        self.assertTrue(hasattr(df, "__data__"))
        self.assertEqual(sys.getrefcount(df), 2)

        self.assertEqual(sys.getrefcount(data["x"]), 3)

    def test_transformations(self):
        """
        Test the use of transformations
        """
        data = {"x": np.array([1, 2, 3], dtype="float32")}
        df = ROOT.RDF.MakeNumpyDataFrame(data)
        df2 = df.Filter("x>1").Define("y", "2*x")
        self.assertEqual(df2.Mean("x").GetValue(), 2.5)
        self.assertEqual(df2.Mean("y").GetValue(), 5)

    def test_delete_dict(self):
        """
        Test behaviour with data dictionary going out of scope
        """
        data = {"x": np.array([1, 2, 3], dtype="float32")}
        df = ROOT.RDF.MakeNumpyDataFrame(data)
        del data
        self.assertEqual(df.Mean("x").GetValue(), 2)

    def test_delete_numpy_array(self):
        """
        Test behaviour with numpy array going out of scope
        """
        x = np.array([1, 2, 3], dtype="float32")
        data = {"x": x}
        df = ROOT.RDF.MakeNumpyDataFrame(data)
        del x
        self.assertEqual(df.Mean("x").GetValue(), 2)

    def test_inplace_dict(self):
        """
        Test behaviour with inplace dictionary
        """
        df = ROOT.RDF.MakeNumpyDataFrame({"x": np.array([1, 2, 3], dtype="float32")})
        self.assertEqual(df.Mean("x").GetValue(), 2)

    def test_lifetime_numpy_array(self):
        """
        Test lifetime of numpy array
        """
        x = np.array([1, 2, 3], dtype="float32")
        gc.collect()
        ref1 = sys.getrefcount(x)

        df = ROOT.RDF.MakeNumpyDataFrame({"x": x})
        gc.collect()
        ref2 = sys.getrefcount(x)
        self.assertEqual(ref2, ref1 + 1)

        del df
        gc.collect()
        ref3 = sys.getrefcount(x)
        self.assertEqual(ref1, ref3)

    def test_lifetime_datasource(self):
        """
        Test lifetime of datasource

        Datasource survives until last node of the graph goes out of scope
        """
        x = np.array([1, 2, 3], dtype="float32")
        gc.collect()
        ref1 = sys.getrefcount(x)

        # Data source has dictionary with RVecs attached, which take a reference
        # to the numpy array
        df = ROOT.RDF.MakeNumpyDataFrame({"x": x})
        m = df.Mean("x")
        gc.collect()
        ref2 = sys.getrefcount(x)
        self.assertEqual(ref1 + 1, ref2)

        # Deleting the root node does not change anything since the datasource
        # owns the RVecs
        del df
        self.assertEqual(m.GetValue(), 2)
        gc.collect()
        ref3 = sys.getrefcount(x)
        self.assertEqual(ref1 + 1, ref3)

        # Deleting the last node releases the RVecs and releases the reference
        # to the numpy array
        del m
        gc.collect()
        ref4 = sys.getrefcount(x)
        self.assertEqual(ref1, ref4)


if __name__ == '__main__':
    unittest.main()
