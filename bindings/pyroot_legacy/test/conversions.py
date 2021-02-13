import unittest
import ROOT

cppcode = """
void stringViewConv(std::string_view) {};
"""

class ListInitialization(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ROOT.gInterpreter.Declare(cppcode)

    def test_string_view_conv(self):
        ROOT.stringViewConv("pyString")

if __name__ == '__main__':
    unittest.main()
