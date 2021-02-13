import unittest
import ROOT


class PrettyPrinting(unittest.TestCase):
    # Helpers
    def _print(self, obj):
        print("print({}) -> {}".format(obj.__cppname__, obj))

    # Tests
    def test_RVec(self):
        x = ROOT.VecOps.RVec("float")(4)
        for i in range(x.size()):
            x[i] = i
        self._print(x)
        self.assertIn("{ 0", x.__str__())

    def test_STLVector(self):
        x = ROOT.std.vector("float")(4)
        for i in range(x.size()):
            x[i] = i
        self._print(x)
        self.assertIn("{ 0", x.__str__())

    def test_STLMap(self):
        x = ROOT.std.map("string", "int")()
        for i, s in enumerate(["foo", "bar"]):
            x[s] = i
        self._print(x)
        self.assertIn("foo", x.__str__())
        self.assertIn("bar", x.__str__())

    def test_STLPair(self):
        x = ROOT.std.pair("string", "int")("foo", 42)
        self._print(x)
        self.assertIn("foo", x.__str__())

    def test_TNamed(self):
        x = ROOT.TNamed("name", "title")
        self._print(x)
        self.assertEqual("Name: name Title: title", x.__str__())

    def test_TObject(self):
        x = ROOT.TObject()
        self._print(x)
        self.assertEqual("Name: TObject Title: Basic ROOT object", x.__str__())

    def test_TH1F(self):
        x = ROOT.TH1F("name", "title", 10, 0, 1)
        self._print(x)
        self.assertEqual("Name: name Title: title NbinsX: 10", x.__str__())

    def test_user_class(self):
        # Test fall-back to __repr__
        ROOT.gInterpreter.Declare('class MyClass {};')
        x = ROOT.MyClass()
        self._print(x)
        s = x.__str__()
        r = x.__repr__()
        self.assertIn("MyClass object at", s)
        self.assertEqual(s, r)


if __name__ == '__main__':
    unittest.main()
