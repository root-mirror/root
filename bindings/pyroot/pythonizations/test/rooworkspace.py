import unittest
import ROOT

class RooWorkspace_test(unittest.TestCase):
    """
    Test for the pythonizations of RooWorkspace.
    """

    # Setup
    def setUp(self):
       self.x = ROOT.RooRealVar("x", "x", 1.337, 0, 10)
       self.ws = ROOT.RooWorkspace("ws", "A workspace")

    # Tests
    def test_import(self):
       self.ws.Import(self.x)
       x = self.ws.var('x')
       self.assertEqual(x.GetName(), 'x')
       self.assertEqual(x.getVal(), self.x.getVal())
    
    def test_import_with_arg(self):
       # Prepare workspace with variables and a PDF
       self.ws.Import(self.x)
       rename = ROOT.RooFit.RenameAllVariables('exp')
       exp = ROOT.RooExponential("exp", "exp", self.x, self.x)
       self.ws.Import(exp, rename)
       
       # Test that rename argument has worked
       x = self.ws.var('x_exp')
       self.assertEqual(x.getVal(), self.x.getVal())
       pdf = self.ws.pdf('exp')
       self.assertGreater(pdf.getVal(), 0)

    def test_import_argset(self):
       argSet = ROOT.RooArgSet(self.x)
       self.ws.Import(argSet)
       x = self.ws.arg('x')
       self.assertEqual(x.GetName(), 'x')
       self.assertEqual(x.getVal(), self.x.getVal())

if __name__ == '__main__':
    unittest.main()

