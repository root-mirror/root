import unittest
import ROOT
import numpy as np
import pickle


def make_tree(*dtypes):
    """
    Make a tree with branches of different data-types
    """
    tree = ROOT.TTree("test", "description")
    col_names = ["col_{}".format(d) for d in dtypes]

    col_vars = []
    for dtype in dtypes:
        if "F" in dtype:
            var = np.empty(1, dtype=np.float32)
        elif "D" in dtype:
            var = np.empty(1, dtype=np.float64)
        elif "I" in dtype:
            var = np.empty(1, dtype=np.int32)
        elif "i" in dtype:
            var = np.empty(1, dtype=np.uint32)
        elif "L" in dtype:
            var = np.empty(1, dtype=np.int64)
        elif "l" in dtype:
            var = np.empty(1, dtype=np.uint64)
        elif "S" in dtype:
            var = np.empty(1, dtype=np.int16)
        elif "s" in dtype:
            var = np.empty(1, dtype=np.uint16)
        elif "B" in dtype:
            var = np.empty(1, dtype=np.int8)
        elif "b" in dtype:
            var = np.empty(1, dtype=np.uint8)
        elif "O" in dtype:
            var = np.empty(1, dtype=np.uint8)
        else:
            raise Exception(
                "Type {} not known to create branch.".format(dtype))
        col_vars.append(var)

    for dtype, name, var in zip(dtypes, col_names, col_vars):
        tree.Branch(name, var, name + "/" + dtype)

    reference = {col: [] for col in col_names}
    for i in range(5):
        for i_var, var in enumerate(col_vars):
            var[0] = i
            reference[col_names[i_var]].append(var[0])
        tree.Fill()
    reference = {col: np.array(reference[col]) for col in reference}

    return tree, reference, dtypes, col_names, col_vars


def create_slice_in_scope():
    """
    Read-out as a numpy array and return a slice which is a view on the data
    """
    df = ROOT.ROOT.RDataFrame(4).Define("x", "(double)rdfentry_")
    npy = df.AsNumpy()
    x = npy["x"]
    x2 = x[:2]
    return x2


class RDataFrameAsNumpy(unittest.TestCase):
    """
    Testing of RDataFrame.AsNumpy pythonization
    """
    def test_branch_dtypes(self):
        """
        Test supported data-types for read-out
        """
        root_dtypes = ["S", "s", "I", "i", "L", "l", "F", "D"]
        tree, ref, _, col_names, _ = make_tree(*root_dtypes)
        df = ROOT.ROOT.RDataFrame(tree)
        npy = df.AsNumpy()
        for col in col_names:
            self.assertTrue(all(npy[col] == ref[col]))

    def test_branch_bool(self):
        """
        Test bool data-type as a special case since we cannot adopt
        the std::vector<bool> with numpy arrays
        """
        df = ROOT.RDataFrame(2).Define("x", "bool(rdfentry_)")
        npy = df.AsNumpy()
        self.assertTrue(bool(npy["x"][0]) == False)
        self.assertTrue(bool(npy["x"][1]) == True)

    def test_read_array(self):
        """
        Testing reading a std::array
        """
        ROOT.gInterpreter.Declare("""
        std::array<unsigned int, 3> create_array(unsigned int n) {
            return std::array<unsigned int, 3>({n, n, n});
        }
        """)
        df = ROOT.ROOT.RDataFrame(5).Define("x", "create_array(rdfentry_)")
        npy = df.AsNumpy()
        self.assertEqual(npy["x"].size, 5)
        self.assertEqual(list(npy["x"][0]), [0, 0, 0])
        self.assertIn("array<unsigned int,3>", str(type(npy["x"][0])))

    def test_read_th1f(self):
        """
        Testing reading a TH1F
        """
        ROOT.gInterpreter.Declare("""
        TH1F create_histo(unsigned int n) {
            const auto str = TString::Format("h%i", n);
            return TH1F(str, str, 4, 0, 1);
        }
        """)
        df = ROOT.ROOT.RDataFrame(5).Define("x", "create_histo(rdfentry_)")
        npy = df.AsNumpy()
        self.assertEqual(npy["x"].size, 5)
        self.assertIn("TH1F", str(type(npy["x"][0])))

    def test_read_vector_constantsize(self):
        """
        Testing reading a std::vector with constant size
        """
        ROOT.gInterpreter.Declare("""
        std::vector<unsigned int> create_vector_constantsize(unsigned int n) {
            return std::vector<unsigned int>({n, n, n});
        }
        """)
        df = ROOT.ROOT.RDataFrame(5).Define("x",
                                       "create_vector_constantsize(rdfentry_)")
        npy = df.AsNumpy()
        self.assertEqual(npy["x"].size, 5)
        self.assertEqual(list(npy["x"][0]), [0, 0, 0])
        self.assertIn("vector<unsigned int>", str(type(npy["x"][0])))

    def test_read_vector_variablesize(self):
        """
        Testing reading a std::vector with variable size
        """
        ROOT.gInterpreter.Declare("""
        std::vector<unsigned int> create_vector_variablesize(unsigned int n) {
            return std::vector<unsigned int>(n);
        }
        """)
        df = ROOT.ROOT.RDataFrame(5).Define("x",
                                       "create_vector_variablesize(rdfentry_)")
        npy = df.AsNumpy()
        self.assertEqual(npy["x"].size, 5)
        self.assertEqual(list(npy["x"][3]), [0, 0, 0])
        self.assertIn("vector<unsigned int>", str(type(npy["x"][0])))

    def test_read_tlorentzvector(self):
        """
        Testing reading a TLorentzVector
        """

        # The global module index does not have it preloaded and
        # gInterpreter.Declare is not allowed to load libPhysics for
        # TLorentzVector. Preload the library now.
        ROOT.gSystem.Load("libPhysics")
        ROOT.gInterpreter.Declare("""
        TLorentzVector create_tlorentzvector() {
            auto v = TLorentzVector();
            v.SetPtEtaPhiM(1, 2, 3, 4);
            return v;
        }
        """)
        df = ROOT.ROOT.RDataFrame(5).Define("x", "create_tlorentzvector()")
        npy = df.AsNumpy()
        self.assertEqual(npy["x"].size, 5)
        self.assertIn("TLorentzVector", str(type(npy["x"][0])))

    def test_read_custom_class(self):
        """
        Testing reading a custom class injected in the interpreter
        """
        ROOT.gInterpreter.Declare("""
        struct CustomClass {
            unsigned int fMember = 42;
        };
        CustomClass create_custom_class() {
            return CustomClass();
        }
        """)
        df = ROOT.ROOT.RDataFrame(5).Define("x", "create_custom_class()")
        npy = df.AsNumpy()
        self.assertEqual(npy["x"].size, 5)
        self.assertIn("CustomClass", str(type(npy["x"][0])))
        self.assertEqual(npy["x"][0].fMember, 42)

    def test_define_columns(self):
        """
        Testing reading defined columns
        """
        df = ROOT.ROOT.RDataFrame(4).Define("x", "1").Define("y", "2").Define(
            "z", "3")
        npy = df.AsNumpy(columns=["x", "y"])
        ref = {"x": np.array([1] * 4), "y": np.array([2] * 4)}
        self.assertTrue(sorted(["x", "y"]) == sorted(npy.keys()))
        self.assertTrue(all(ref["x"] == npy["x"]))
        self.assertTrue(all(ref["y"] == npy["y"]))

    def test_exclude_columns(self):
        """
        Testing excluding columns from read-out
        """
        df = ROOT.ROOT.RDataFrame(4).Define("x", "1").Define("y", "2").Define(
            "z", "3")
        npy = df.AsNumpy(exclude=["z"])
        ref = {"x": np.array([1] * 4), "y": np.array([2] * 4)}
        self.assertTrue(sorted(["x", "y"]) == sorted(npy.keys()))
        self.assertTrue(all(ref["x"] == npy["x"]))
        self.assertTrue(all(ref["y"] == npy["y"]))

        df2 = ROOT.ROOT.RDataFrame(4).Define("x", "1").Define("y", "2").Define(
            "z", "3")
        npy = df2.AsNumpy(columns=["x", "y"], exclude=["y"])
        ref = {"x": np.array([1] * 4)}
        self.assertTrue(["x"] == list(npy.keys()))
        self.assertTrue(all(ref["x"] == npy["x"]))

    def test_numpy_result_ptr(self):
        """
        Testing result pointer being attribute of returned numpy array
        """
        df = ROOT.ROOT.RDataFrame(4).Define("x", "(double)rdfentry_")
        npy = df.AsNumpy()
        x = npy["x"]
        self.assertTrue(hasattr(x, "result_ptr"))

    def test_numpy_slice(self):
        """
        Testing ownership of numpy array as owner of the data
        """
        df = ROOT.ROOT.RDataFrame(4).Define("x", "(double)rdfentry_")
        npy = df.AsNumpy()
        x = npy["x"]
        ref = np.array([0, 1, 2, 3])
        self.assertTrue(all(x == ref))

        x2 = x[:2]
        ref2 = ref[:2]
        self.assertTrue(all(x2 == ref2))

    def test_numpy_slice_in_scope(self):
        """
        Testing ownership of numpy array as view on data
        """
        x = create_slice_in_scope()
        self.assertEqual(x.flags["OWNDATA"], False)
        ref = np.array([0, 1])
        self.assertTrue(all(x == ref))
        self.assertTrue(hasattr(x, "result_ptr"))

    def test_empty_array(self):
        """
        Testing readout of empty std::vectors
        """
        df = ROOT.ROOT.RDataFrame(1).Define("x", "std::vector<float>()")
        npy = df.AsNumpy(["x"])
        self.assertEqual(npy["x"].size, 1)
        self.assertTrue(npy["x"][0].empty())

    def test_empty_selection(self):
        """
        Testing readout of empty selection
        """
        df = ROOT.ROOT.RDataFrame(10).Define("x", "1.0").Filter("x<0")
        npy = df.AsNumpy(["x"])
        self.assertEqual(npy["x"].size, 0)

    def test_pickle(self):
        """
        Testing pickling of returned numpy array
        """
        df = ROOT.ROOT.RDataFrame(10).Define("x", "1.0")
        npy = df.AsNumpy(["x"])
        arr = npy["x"]

        with open("rdataframe_asnumpy.pickle", "wb") as f:
            pickle.dump(arr, f)
        with open("rdataframe_asnumpy.pickle", "rb") as f:
            arr2 = pickle.load(f)
        self.assertTrue(all(arr == arr2))

    def test_memory_adoption_fundamental_types(self):
        """
        Testing the adoption of the memory from the C++ side for fundamental types
        """
        df = ROOT.ROOT.RDataFrame(1).Define("x", "1.0")
        npy = df.AsNumpy(["x"])
        pyarr = npy["x"]
        cpparr = pyarr.result_ptr.GetValue()
        pyarr[0] = 42
        self.assertTrue(cpparr[0] == pyarr[0])

    def test_memory_adoption_complex_types(self):
        """
        Testing the adoption of the memory from the C++ side for complex types
        """
        df = ROOT.ROOT.RDataFrame(1).Define("x", "std::vector<float>({1, 2, 3})")
        npy = df.AsNumpy(["x"])
        pyarr = npy["x"]
        cpparr = pyarr.result_ptr.GetValue()
        pyarr[0][0] = 42
        self.assertTrue(cpparr[0][0] == pyarr[0][0])


if __name__ == '__main__':
    unittest.main()
