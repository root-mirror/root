import unittest

import ROOT


class TStringConverter(unittest.TestCase):
    """
    Tests for passing a Python string to a C++ function that expects a TString.

    This feature is not implemented by a PyROOT pythonization, but by a converter
    that was added to Cppyy to create a TString out of a Python string.
    """

    test_str = "test"

    # Helpers
    def check_type_conversion(self):
        s = ROOT.TString(self.test_str)

        # Works with TString...
        self.assertEqual(ROOT.myfun(s), self.test_str)
        # ... and Python string
        self.assertEqual(ROOT.myfun(self.test_str), self.test_str)

    # Tests
    def test_by_value(self):
        ROOT.gInterpreter.Declare("""
        const char* myfun(TString s) { return s.Data(); }
        """)

        self.check_type_conversion()

    def test_by_reference(self):
        ROOT.gInterpreter.Declare("""
        const char* myfun(TString &s) { return s.Data(); }
        """)

        self.check_type_conversion()

    def test_by_const_reference(self):
        ROOT.gInterpreter.Declare("""
        const char* myfun(const TString &s) { return s.Data(); }
        """)

        self.check_type_conversion()


if __name__ == '__main__':
    unittest.main()
