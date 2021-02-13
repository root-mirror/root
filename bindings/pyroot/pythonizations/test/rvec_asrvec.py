import unittest
import ROOT
import numpy as np
import sys
import gc


def get_maximum_for_dtype(dtype):
    if np.issubdtype(dtype, np.integer):
        return np.iinfo(dtype).max
    if np.issubdtype(dtype, np.floating):
        return np.finfo(dtype).max

def get_minimum_for_dtype(dtype):
    if np.issubdtype(dtype, np.integer):
        return np.iinfo(dtype).min
    if np.issubdtype(dtype, np.floating):
        return np.finfo(dtype).min


class AsRVec(unittest.TestCase):
    """
    Tests for the AsRVec feature enabling to adopt memory of Python objects
    with an array interface member using RVec as C++ container.
    """

    # Helpers
    dtypes = [
        "int32", "int64", "uint32", "uint64", "float32", "float64"
    ]

    def check_memory_adoption(self, root_obj, np_obj):
        np_obj[0] = get_maximum_for_dtype(np_obj.dtype)
        np_obj[1] = get_minimum_for_dtype(np_obj.dtype)
        self.assertEqual(root_obj[0], np_obj[0])
        self.assertEqual(root_obj[1], np_obj[1])

    def check_size(self, expected_size, obj):
        self.assertEqual(expected_size, obj.size())

    # Tests
    def test_dtypes(self):
        """
        Test adoption of numpy arrays with different data types
        """
        for dtype in self.dtypes:
            np_obj = np.empty(2, dtype=dtype)
            root_obj = ROOT.VecOps.AsRVec(np_obj)
            self.check_memory_adoption(root_obj, np_obj)
            self.check_size(2, root_obj)

    def test_multidim(self):
        """
        Test adoption of multi-dimensional numpy arrays
        """
        np_obj = np.array([[1, 2], [3, 4]], dtype="float32")
        rvec = ROOT.VecOps.AsRVec(np_obj)
        self.assertEqual(rvec.size(), 4)

    def test_size_zero(self):
        """
        Test adoption of numpy array with size 0
        """
        np_obj = np.array([], dtype="float32")
        rvec = ROOT.VecOps.AsRVec(np_obj)
        self.assertEqual(rvec.size(), 0)

    def test_adopt_rvec(self):
        """
        Test adoption of RVecs
        """
        rvec = ROOT.VecOps.RVec("float")(1)
        rvec[0] = 42
        rvec2 = ROOT.VecOps.AsRVec(rvec)
        self.assertEqual(rvec.size(), rvec2.size())
        self.assertEqual(rvec[0], rvec2[0])
        rvec2[0] = 43
        self.assertEqual(rvec[0], rvec2[0])

    def test_ownership(self):
        """
        Test ownership of returned RVec (to be owned by Python)
        """
        np_obj = np.array([1, 2])
        rvec = ROOT.VecOps.AsRVec(np_obj)
        self.assertEqual(rvec.__python_owns__, True)

    def test_attribute_adopted(self):
        """
        Test __adopted__ attribute of returned RVecs
        """
        np_obj = np.array([1, 2])
        rvec = ROOT.VecOps.AsRVec(np_obj)
        self.assertTrue(hasattr(rvec, "__adopted__"))
        self.assertEqual(id(rvec.__adopted__), id(np_obj))

    def test_refcount(self):
        """
        Test reference count of returned RVec

        We expect a refcount of 2 for the RVec because the call to sys.getrefcount
        creates a second reference by itself.
        We attach the adopted pyobject to the RVec and increase the refcount of the
        numpy array. After deletion of the rvec, the refcount of the numpy array
        is decreased.
        """
        np_obj = np.array([1, 2])
        rvec = ROOT.VecOps.AsRVec(np_obj)
        gc.collect()
        self.assertEqual(sys.getrefcount(rvec), 2)
        self.assertEqual(sys.getrefcount(np_obj), 3)
        del rvec
        gc.collect()
        self.assertEqual(sys.getrefcount(np_obj), 2)


if __name__ == '__main__':
    unittest.main()
