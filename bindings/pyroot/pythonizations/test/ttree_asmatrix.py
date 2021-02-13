import unittest
import ROOT
import numpy as np
from sys import maxsize


class TTreeAsMatrix(unittest.TestCase):
    is_64bit = True if maxsize > 2**32 else False

    # Helpers
    def make_tree(self, *dtypes):
        tree = ROOT.TTree("test", "description")
        col_names = ["col{}".format(i) for i in range(len(dtypes))]

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
                if self.is_64bit:
                    var = np.empty(1, dtype=np.int64)
                else:
                    var = np.empty(1, dtype=np.int32)
            elif "l" in dtype:
                if self.is_64bit:
                    var = np.empty(1, dtype=np.uint64)
                else:
                    var = np.empty(1, dtype=np.uint32)
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

        reference = []
        for i in range(4):
            row = []
            for i_var, var in enumerate(col_vars):
                var[0] = (i + i_var + 0.5) * np.power(-1, i)
                row.append(var[0])
            reference.append(row)
            tree.Fill()

        return tree, reference, dtypes, col_names, col_vars

    def make_example(self, *dtypes):
        tree, reference, dtypes, col_names, col_vars = self.make_tree(*dtypes)
        matrix_ttree = tree.AsMatrix()
        matrix_ref = np.asarray(reference)
        return matrix_ttree, matrix_ref

    # Tests
    def test_instance(self):
        """
        Test instantiation
        """
        tree, _, _, col_names, _ = self.make_tree("F", "F")
        tree.AsMatrix(col_names)

    def test_return_labels(self):
        """
        Test returned labels
        """
        tree, _, _, col_names, _ = self.make_tree("F", "F")
        matrix, labels = tree.AsMatrix(col_names, return_labels=True)
        self.assertEqual(labels, col_names)

    def test_exclude_columns(self):
        """
        Test excluding columns
        """
        tree, reference, _, _, _ = self.make_tree("F", "F")
        matrix_ttree = tree.AsMatrix(exclude=["col0"])
        matrix_ref = np.asarray([x[1] for x in reference])
        for value_ttree, value_ref in zip(matrix_ttree, matrix_ref):
            self.assertEqual(value_ttree, value_ref)

    def test_not_supported_dtype(self):
        """
        Test error-handling of unsupported data-types
        """
        tree, _, _, col_names, _ = self.make_tree("F", "F")
        try:
            tree.AsMatrix(col_names, dtype="foo")
            self.assertFail()
        except Exception as exception:
            self.assertIn("Data-type foo is not supported, select from",
                          exception.args[0])

    def test_not_existent_column(self):
        """
        Test if a column does not exist
        """
        tree, _, _, col_names, _ = self.make_tree("F", "F")
        try:
            tree.AsMatrix(["foo"])
            self.assertFail()
        except Exception as exception:
            self.assertIn("branch not existent", exception.args[0])

    def test_no_branches_selected(self):
        """
        Test if no branch is selected by the requested columns
        """
        tree, _, _, col_names, _ = self.make_tree("F", "F")
        try:
            tree.AsMatrix(columns=["col0"], exclude=["col0"])
            self.assertFail()
        except Exception as exception:
            self.assertIn("Arguments resulted in no selected branches.",
                          exception.args[0])

    def test_shape(self):
        """
        Test shape of returned matrix
        """
        matrix_ttree, matrix_ref = self.make_example("F", "F")
        for i in range(2):
            self.assertEqual(matrix_ttree.shape[i], matrix_ref.shape[i])

    def test_take_all_columns(self):
        """
        Test taking all columns
        """
        tree, reference, _, _, _ = self.make_tree("F", "F")
        matrix_ttree = tree.AsMatrix()
        matrix_ref = np.asarray(reference)
        for i in range(matrix_ref.shape[0]):
            for j in range(matrix_ref.shape[1]):
                self.assertEqual(matrix_ttree[i, j], matrix_ref[i, j])

    def test_values(self):
        """
        Test correctness of returned values
        """
        matrix_ttree, matrix_ref = self.make_example("F", "F")
        for i in range(matrix_ref.shape[0]):
            for j in range(matrix_ref.shape[1]):
                self.assertEqual(matrix_ttree[i, j], matrix_ref[i, j])

    def test_zero_entries(self):
        """
        Test behaviour for a tree with zero entries
        """
        tree = ROOT.TTree("test", "description")
        var = np.empty(1, np.float32)
        tree.Branch("col", var, "col/F")
        try:
            tree.AsMatrix(["col"])
            self.assertFail()
        except Exception as exception:
            self.assertEqual("Tree test has no entries.", exception.args[0])

    def test_multiple_leaves(self):
        """
        Test a column with multiple leaves
        """
        tree = ROOT.TTree("test", "description")
        var = np.ones(2, np.float32)
        tree.Branch("col", var, "sub1/F:sub2/F")
        tree.Fill()
        try:
            tree.AsMatrix(["col"])
            self.assertFail()
        except Exception as exception:
            self.assertIn("branch has multiple leaves", exception.args[0])

    def test_dtype_conversion(self):
        """
        Test conversion of data-types of the columns to uniform output data-type
        """
        numpy_dtype = {
            "unsigned int": np.dtype(np.uint32),
            "int": np.dtype(np.int32),
            "unsigned long": np.dtype(np.uint64),
            "long": np.dtype(np.int64),
            "float": np.dtype(np.float32),
            "double": np.dtype(np.float64)
        }
        if not self.is_64bit:
            numpy_dtype["long"] = np.dtype(np.int32)
            numpy_dtype["unsigned long"] = np.dtype(np.uint32)
        tree, reference, _, _, _ = self.make_tree("F", "F")
        matrix_ref = np.asarray(reference, dtype=np.float64)
        for dtype in numpy_dtype:
            matrix_ttree = tree.AsMatrix(dtype=dtype)
            self.assertEqual(matrix_ttree.dtype.name, numpy_dtype[dtype].name)

    def test_branch_dtypes(self):
        """
        Test branches of different data-types
        """
        root_dtypes = ["B", "b", "S", "s", "I", "i", "L", "l", "F", "D"]
        tree, _, _, col_names, _ = self.make_tree(*root_dtypes)
        tree.AsMatrix()

        try:
            tree, _, _, col_names, _ = self.make_tree("O")
            tree.AsMatrix()
        except Exception as exception:
            self.assertIn("branch has unsupported data-type ['Bool_t']",
                          exception.args[0])

    def test_tchain(self):
        """
        Test whether the pythonization works as well for a TChain
        """
        tree1 = ROOT.TTree("test1", "description")
        tree2 = ROOT.TTree("test2", "description")
        var = np.ones(1, np.float32)
        tfile = ROOT.TFile("test_ttree_asmatrix.root", "RECREATE")
        tree1.Branch("col", var, "col/F")
        tree1.Fill()
        tree2.Branch("col", var, "col/F")
        tree2.Fill()
        tree1.Write()
        tree2.Write()
        tfile.Close()

        chain = ROOT.TChain("chain")
        chain.Add("test_ttree_asmatrix.root/test1")
        chain.Add("test_ttree_asmatrix.root/test2")
        m = chain.AsMatrix(["col"])
        self.assertTrue((m == np.ones((2, 1), dtype=m.dtype)).all())


if __name__ == '__main__':
    unittest.main()
