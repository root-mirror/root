import unittest
import ROOT


class PrettyPrinting(unittest.TestCase):
    # Helpers
    def _print(self, obj):
        print("print({}) -> {}".format(repr(obj), obj))

    # Tests
    def test_RVec(self):
        x = ROOT.ROOT.VecOps.RVec("float")(4)
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

    def test_STLString(self):
         # std::string is not pythonized with the pretty printing, because:
         # 1. gInterpreter->ToString("s") returns ""s""
         # 2. cppyy already does the right thing
         s = ROOT.std.string("x")
         self.assertEqual(str(s), "x")

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

    def test_null_object(self):
        # ROOT-9935: test null proxied cpp object
        x = ROOT.MakeNullPointer("TTree")
        s = x.__str__()
        r = x.__repr__()
        self.assertIn("TTree object at", s)
        self.assertEqual(s, r)

    def test_user_class_with_str(self):
        # ROOT-10967: Respect existing __str__ method defined in C++
        ROOT.gInterpreter.Declare('struct MyClassWithStr { std::string __str__() { return "foo"; } };')
        x = ROOT.MyClassWithStr()
        self._print(x)
        s = x.__str__()
        r = x.__repr__()
        self.assertIn("MyClassWithStr object at", r)
        self.assertEqual(s, "foo")

        # Test inherited class
        ROOT.gInterpreter.Declare('struct MyClassWithStr2 : public MyClassWithStr { };')
        x2 = ROOT.MyClassWithStr2()
        self._print(x2)
        s2 = x2.__str__()
        r2 = x2.__repr__()
        self.assertIn("MyClassWithStr2 object at", r2)
        self.assertEqual(s2, "foo")


    # TNamed and TObject are not pythonized because these object are touched
    # by PyROOT before any pythonizations are added. Following, the classes
    # are not piped through the pythonizor functions again.
    """
    def test_TNamed(self):
        x = ROOT.TNamed("name", "title")
        self._print(x)
        self.assertEqual("Name: name Title: title", x.__str__())

    def test_TObject(self):
        x = ROOT.TObject()
        self._print(x)
        self.assertEqual("Name: TObject Title: Basic ROOT object", x.__str__())
    """


if __name__ == '__main__':
    unittest.main()
