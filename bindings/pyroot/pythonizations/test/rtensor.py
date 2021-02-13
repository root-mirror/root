import unittest
import ROOT
RTensor = ROOT.TMVA.Experimental.RTensor
import numpy as np


def check_shape(root_obj, np_obj):
    root_shape = tuple(root_obj.GetShape())
    np_shape = tuple(np_obj.shape)
    return root_shape == np_shape


class AsRTensor(unittest.TestCase):
    """
    Test AsRTensor adoption mechanism
    """

    # Helpers
    dtypes = [
        "int32", "int64", "uint32", "uint64", "float32", "float64"
    ]

    # Tests
    def test_dtypes(self):
        """
        Test adoption of numpy arrays with different data types
        """
        for dtype in self.dtypes:
            np_obj = np.array([[1, 2], [3, 4], [5, 6]], dtype=dtype)
            root_obj = ROOT.TMVA.Experimental.AsRTensor(np_obj)
            self.assertTrue(check_shape(root_obj, np_obj))
            np_obj[0,0] = 42
            self.assertTrue(root_obj[0,0] == 42)

    def test_memoryLayout(self):
        """
        Test adoption of the memory layout
        """
        np_obj = np.array([[1, 2], [3, 4], [5, 6]])
        root_obj = ROOT.TMVA.Experimental.AsRTensor(np_obj)
        self.assertTrue(np_obj.flags.c_contiguous)
        self.assertEqual(root_obj.GetMemoryLayout(), 1)

        np_obj2 = np_obj.T
        root_obj2 = ROOT.TMVA.Experimental.AsRTensor(np_obj2)
        self.assertTrue(np_obj2.flags.f_contiguous)
        self.assertEqual(root_obj2.GetMemoryLayout(), 2)

    def test_strides(self):
        """
        Test adoption of the strides

        Note that numpy multiplies the strides with the size of the element
        in bytes.
        """
        np_obj = np.array([[1, 2], [3, 4], [5, 6]], dtype="float32")
        root_obj = ROOT.TMVA.Experimental.AsRTensor(np_obj)

        np_strides = np_obj.strides
        root_strides = root_obj.GetStrides()
        self.assertEqual(len(np_strides), len(root_strides))
        for x, y in zip(np_strides, root_strides):
            self.assertEqual(x, y * 4)

        np_obj = np_obj.T
        root_obj = ROOT.TMVA.Experimental.AsRTensor(np_obj)
        np_strides = np_obj.strides
        root_strides = root_obj.GetStrides()
        self.assertEqual(len(np_strides), len(root_strides))
        for x, y in zip(np_strides, root_strides):
            self.assertEqual(x, y * 4)


class ArrayInterface(unittest.TestCase):
    """
    Test memory adoption of RTensor array interface.
    """

    # Helpers
    dtypes = [
        "int", "unsigned int", "long", "long long", "Long64_t", "unsigned long",
        "unsigned long long", "ULong64_t", "float", "double"
    ]

    def get_maximum_for_dtype(self, dtype):
        if np.issubdtype(dtype, np.integer):
            return np.iinfo(dtype).max
        if np.issubdtype(dtype, np.floating):
            return np.finfo(dtype).max

    def get_minimum_for_dtype(self, dtype):
        if np.issubdtype(dtype, np.integer):
            return np.iinfo(dtype).min
        if np.issubdtype(dtype, np.floating):
            return np.finfo(dtype).min

    # Tests
    def test_memoryAdoption(self):
        """
        Test correct adoption of different datatypes
        """
        shape = ROOT.std.vector("size_t")((2, 2))
        for dtype in self.dtypes:
            root_obj = RTensor(dtype)(shape)
            np_obj = np.asarray(root_obj)
            self.assertTrue(check_shape(root_obj, np_obj))
            np_obj[0,0] = 42
            self.assertTrue(root_obj[0,0] == 42)

    def test_memoryLayout(self):
        """
        Test adoption of the memory layout
        """
        shape = ROOT.std.vector("size_t")((2, 2))
        x = RTensor("float")(shape)
        y = np.asarray(x)
        self.assertTrue(y.flags.c_contiguous)

        x = x.Transpose()
        y = np.asarray(x)
        self.assertTrue(y.flags.f_contiguous)

    def test_ownData(self):
        """
        Test ownership of adopted numpy array
        """
        shape = ROOT.std.vector("size_t")((2, 2))
        x = RTensor("float")(shape)
        y = np.asarray(x)
        self.assertFalse(y.flags.owndata)

        y = np.transpose(y)
        self.assertFalse(y.flags.owndata)

        y = np.copy(y)
        self.assertTrue(y.flags.owndata)


class NumpyCompliance(unittest.TestCase):
    """
    Test compliance of the RTensor methods with the numpy interface
    """

    def test_transpose(self):
        """
        Test np.transpose vs RTensor::Transpose
        """
        shape = ROOT.std.vector("size_t")((2, 3))
        x = RTensor("float")(shape)
        y = np.asarray(x)
        self.assertEqual(x.GetMemoryLayout(), 1)
        self.assertEqual(y.flags.c_contiguous, True)

        for i, j in zip(x.GetShape(), y.shape):
            self.assertEqual(i, j)

        x = x.Transpose()
        y = np.transpose(y)
        self.assertEqual(x.GetMemoryLayout(), 2)
        self.assertEqual(y.flags.f_contiguous, True)

        for i, j in zip(x.GetShape(), y.shape):
            self.assertEqual(i, j)

    def test_expandDims(self):
        """
        Test np.expand_dims vs RTensor::ExpandDims
        """
        shape = ROOT.std.vector("size_t")((2, 2))
        x = RTensor("float")(shape)
        y = np.asarray(x)

        x1 = x.ExpandDims(0)
        y1 = np.expand_dims(y, 0)
        for i, j in zip(x1.GetShape(), y1.shape):
            self.assertEqual(i, j)

        x2 = x.ExpandDims(-1)
        y2 = np.expand_dims(y, -1)
        for i, j in zip(x2.GetShape(), y2.shape):
            self.assertEqual(i, j)

    def test_squeeze(self):
        """
        Test np.squeeze vs RTensor::Squeeze
        """
        shape = ROOT.std.vector("size_t")((1, 2))
        x = RTensor("float")(shape)
        y = np.asarray(x)

        x = x.Squeeze()
        y = np.squeeze(y)

        for i, j in zip(x.GetShape(), y.shape):
            self.assertEqual(i, j)

    def test_reshape(self):
        """
        Test np.reshape vs RTensor::Reshape
        """
        shape = ROOT.std.vector("size_t")((1, 2))
        x = RTensor("float")(shape)
        y = np.asarray(x)
        y[0,0] = 1
        y[0,1] = 2

        shape[0] = 2
        shape[1] = 1
        x = x.Reshape(shape)
        y = np.reshape(y, shape)

        for i, j in zip(x.GetShape(), y.shape):
            self.assertEqual(i, j)

        self.assertEqual(x[0,0], y[0,0])
        self.assertEqual(x[1,0], y[1,0])

        self.assertEqual(x.GetMemoryLayout(), 1)
        self.assertEqual(y.flags.c_contiguous, True)

    def test_slice(self):
        """
        Test slicing operations
        """
        shape = ROOT.std.vector("size_t")((2, 2))
        x = RTensor("float")(shape)
        y = np.asarray(x)
        y[0,0] = 1
        y[0,1] = 2
        y[1,0] = 3
        y[1,1] = 4

        x1 = x[:,0]
        y1 = y[:,0]
        for i, j in zip(x1.GetShape(), y1.shape):
            self.assertEqual(i, j)
        self.assertEqual(x1[0], y1[0])
        self.assertEqual(x1[1], y1[1])

        x2 = x[:,:]
        y2 = y[:,:]
        for i, j in zip(x2.GetShape(), y2.shape):
            self.assertEqual(i, j)
        self.assertEqual(x2[0,0], y2[0,0])
        self.assertEqual(x2[0,1], y2[0,1])
        self.assertEqual(x2[1,0], y2[1,0])
        self.assertEqual(x2[1,1], y2[1,1])

        # We know that we differ in this case since numpy
        # does only squeeze the dimensions of a slice if
        # the dimension is requested by a single index.
        x3 = x[1:2,:]
        y3 = np.squeeze(y[1:2,:])
        for i, j in zip(x3.GetShape(), y3.shape):
            self.assertEqual(i, j)
        self.assertEqual(x3[0], y3[0])
        self.assertEqual(x3[1], y3[1])

        x4 = x[:,-2]
        y4 = y[:,-2]
        for i, j in zip(x4.GetShape(), y4.shape):
            self.assertEqual(i, j)
        self.assertEqual(x4[0], y4[0])
        self.assertEqual(x4[1], y4[1])
