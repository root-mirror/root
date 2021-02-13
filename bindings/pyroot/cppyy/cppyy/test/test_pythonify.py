import py, os, sys
from pytest import raises
from .support import setup_make, pylong

currpath = py.path.local(__file__).dirpath()
test_dct = str(currpath.join("example01Dict"))

def setup_module(mod):
    setup_make("example01")


class TestPYTHONIFY:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.example01 = cppyy.load_reflection_info(cls.test_dct)

    def test01_load_dictionary_cache(self):
        """Test whether loading a dictionary twice results in the same object"""

        import cppyy
        lib2 = cppyy.load_reflection_info(self.test_dct)
        assert self.example01 is lib2

    def test02_finding_classes(self):
        """Test the lookup of a class, and its caching"""

        import cppyy
        example01_class = cppyy.gbl.example01
        cl2 = cppyy.gbl.example01
        assert example01_class is cl2

        with raises(AttributeError):
            cppyy.gbl.nonexistingclass

    def test03_calling_static_functions(self):
        """Test calling of static methods"""

        import cppyy, sys, math
        example01_class = cppyy.gbl.example01
        res = example01_class.staticAddOneToInt(1)
        assert res == 2

        res = example01_class.staticAddOneToInt(pylong(1))
        assert res == 2
        res = example01_class.staticAddOneToInt(1, 2)
        assert res == 4
        res = example01_class.staticAddOneToInt(-1)
        assert res == 0
        maxint32 = int(2 ** 31 - 1)
        res = example01_class.staticAddOneToInt(maxint32-1)
        assert res == maxint32
        res = example01_class.staticAddOneToInt(maxint32)
        assert res == -maxint32-1

        raises(TypeError, example01_class.staticAddOneToInt, 1, [])
        raises(TypeError, example01_class.staticAddOneToInt, 1.)
        raises(TypeError, example01_class.staticAddOneToInt, maxint32+1)
        res = example01_class.staticAddToDouble(0.09)
        assert res == 0.09 + 0.01

        res = example01_class.staticAtoi("1")
        assert res == 1

        res = example01_class.staticStrcpy("aap")     # TODO: this leaks
        assert res == "aap"
        res = example01_class.staticStrcpy(u"aap")    # TODO: id.
        assert res == "aap"
        raises(TypeError, example01_class.staticStrcpy, 1.)    # TODO: id.

    def test04_constructing_and_calling(self):
        """Test object and method calls"""

        import cppyy
        example01_class = cppyy.gbl.example01
        assert example01_class.getCount() == 0
        instance = example01_class(7)
        assert example01_class.getCount() == 1
        res = instance.addDataToInt(4)
        assert res == 11
        res = instance.addDataToInt(-4)
        assert res == 3
        instance.__destruct__()
        assert example01_class.getCount() == 0
        raises(ReferenceError, instance.addDataToInt, 4)

        instance = example01_class(7)
        instance2 = example01_class(8)
        assert example01_class.getCount() == 2
        instance.__destruct__()
        assert example01_class.getCount() == 1
        instance2.__destruct__()
        assert example01_class.getCount() == 0

        instance = example01_class(13)
        res = instance.addDataToDouble(16)
        assert round(res-29, 8) == 0.
        instance.__destruct__()
        instance = example01_class(-13)
        res = instance.addDataToDouble(16)
        assert round(res-3, 8) == 0.
        instance.__destruct__() 

        instance = example01_class(42)
        assert example01_class.getCount() == 1

        res = instance.addDataToAtoi("13")
        assert res == 55

        res = instance.addToStringValue("12")    # TODO: this leaks
        assert res == "54"
        res = instance.addToStringValue("-12")   # TODO: this leaks
        assert res == "30"

        res = instance.staticAddOneToInt(pylong(1))
        assert res == 2

        instance.__destruct__()
        assert example01_class.getCount() == 0

    def test05_passing_object_by_pointer(self):
        """Pass object by pointer"""

        import cppyy
        example01_class = cppyy.gbl.example01
        payload_class = cppyy.gbl.payload

        e = example01_class(14)
        pl = payload_class(3.14)
        assert round(pl.getData()-3.14, 8) == 0

        example01_class.staticSetPayload(pl, 41.)
        assert pl.getData() == 41.
        example01_class.staticSetPayload(pl, 43.)
        assert pl.getData() == 43.
        e.staticSetPayload(pl, 45.)
        assert pl.getData() == 45.

        e.setPayload(pl)
        assert round(pl.getData()-14., 8) == 0

        pl.__destruct__()
        e.__destruct__()
        assert example01_class.getCount() == 0

    def test06_returning_object_by_pointer(self):
        """Return an object py pointer"""

        import cppyy
        example01_class = cppyy.gbl.example01
        payload_class = cppyy.gbl.payload

        pl = payload_class(3.14)
        assert round(pl.getData()-3.14, 8) == 0

        pl2 = example01_class.staticCyclePayload(pl, 38.)
        assert pl2.getData() == 38.

        e = example01_class(14)

        pl2 = e.cyclePayload(pl)
        assert round(pl2.getData()-14., 8) == 0

        pl.__destruct__()
        e.__destruct__()
        assert example01_class.getCount() == 0

    def test07_returning_object_by_value(self):
        """Return an object by value"""

        import cppyy
        example01_class = cppyy.gbl.example01
        payload_class = cppyy.gbl.payload

        pl = payload_class(3.14)
        assert round(pl.getData()-3.14, 8) == 0

        pl2 = example01_class.staticCopyCyclePayload(pl, 38.)
        assert pl2.getData() == 38.
        pl2.__destruct__()

        e = example01_class(14)

        pl2 = e.copyCyclePayload(pl)
        assert round(pl2.getData()-14., 8) == 0
        pl2.__destruct__()

        pl.__destruct__()
        e.__destruct__()
        assert example01_class.getCount() == 0

    def test08_global_functions(self):
        """Call a global function"""

        import cppyy

        assert cppyy.gbl.globalAddOneToInt(3) == 4     # creation lookup
        assert cppyy.gbl.globalAddOneToInt(3) == 4     # cached lookup

        assert cppyy.gbl.ns_example01.globalAddOneToInt(4) == 5
        assert cppyy.gbl.ns_example01.globalAddOneToInt(4) == 5

    def test09_memory(self):
        """Test proper C++ destruction by the garbage collector"""

        import cppyy, gc
        example01_class = cppyy.gbl.example01
        payload_class = cppyy.gbl.payload

        pl = payload_class(3.14)
        assert payload_class.count == 1
        assert round(pl.getData()-3.14, 8) == 0

        pl2 = example01_class.staticCopyCyclePayload(pl, 38.)
        assert payload_class.count == 2
        assert pl2.getData() == 38.
        pl2 = None
        gc.collect()
        assert payload_class.count == 1

        e = example01_class(14)

        pl2 = e.copyCyclePayload(pl)
        assert payload_class.count == 2
        assert round(pl2.getData()-14., 8) == 0
        pl2 = None
        gc.collect()
        assert payload_class.count == 1

        pl = None
        e = None
        gc.collect()
        assert payload_class.count == 0
        assert example01_class.getCount() == 0

        pl = payload_class(3.14)
        pl_a = example01_class.staticCyclePayload(pl, 66.)
        pl_a.getData() == 66.
        assert payload_class.count == 1
        pl_a = None
        pl = None
        gc.collect()
        assert payload_class.count == 0

        # TODO: need ReferenceError on touching pl_a

    def test10_default_arguments(self):
        """Test propagation of default function arguments"""

        import cppyy
        a = cppyy.gbl.ArgPasser()

        # NOTE: when called through the stub, default args are fine
        f = a.stringRef
        s = cppyy.gbl.std.string
        assert f(s("aap"), 0, s("noot")) == "aap"
        assert f(s("noot"), 1) == "default"
        assert f(s("mies")) == "mies"

        for itype in ['short', 'ushort', 'int', 'uint', 'long', 'ulong']:
            g = getattr(a, '%sValue' % itype)
            raises(TypeError, g, 1, 2, 3, 4, 6)
            assert g(11, 0, 12, 13) == 11
            assert g(11, 1, 12, 13) == 12
            assert g(11, 1, 12)     == 12
            assert g(11, 2, 12)     ==  2
            assert g(11, 1)         ==  1
            assert g(11, 2)         ==  2
            assert g(11)            == 11

        for ftype in ['float', 'double']:
            g = getattr(a, '%sValue' % ftype)
            raises(TypeError, g, 1., 2, 3., 4., 6.)
            assert g(11., 0, 12., 13.) == 11.
            assert g(11., 1, 12., 13.) == 12.
            assert g(11., 1, 12.)      == 12.
            assert g(11., 2, 12.)      ==  2.
            assert g(11., 1)           ==  1.
            assert g(11., 2)           ==  2.
            assert g(11.)              == 11.

    def test11_overload_on_arguments(self):
        """Test functions overloaded on arguments"""

        import cppyy
        e = cppyy.gbl.example01(1)

        assert e.addDataToInt(2)                 ==  3
        assert e.overloadedAddDataToInt(3)       ==  4
        assert e.overloadedAddDataToInt(4, 5)    == 10
        assert e.overloadedAddDataToInt(6, 7, 8) == 22

    def test12_typedefs(self):
        """Test access and use of typedefs"""

        import cppyy

        assert cppyy.gbl.example01 == cppyy.gbl.example01_t

    def test13_underscore_in_class_name(self):
        """Test recognition of '_' as part of a valid class name"""

        import cppyy

        assert cppyy.gbl.z_ == cppyy.gbl.z_

        z = cppyy.gbl.z_()

        assert hasattr(z, 'myint')
        assert z.gime_z_(z)

    def test14_bound_unbound_calls(self):
        """Test (un)bound method calls"""

        import cppyy

        raises(TypeError, cppyy.gbl.example01.addDataToInt, 1)

        meth = cppyy.gbl.example01.addDataToInt
        raises(TypeError, meth)
        raises(TypeError, meth, 1)

        e = cppyy.gbl.example01(2)
        assert 5 == meth(e, 3)

    def test15_installable_function(self):
       """Test installing and calling global C++ function as python method"""

       import cppyy

       cppyy.gbl.example01.fresh = cppyy.gbl.installableAddOneToInt

       e = cppyy.gbl.example01(0)
       assert 2 == e.fresh(1)
       assert 3 == e.fresh(2)

    def test16_subclassing(self):
        """A sub-class on the python side should have that class as type"""

        import cppyy, gc
        gc.collect()

        example01 = cppyy.gbl.example01

        assert example01.getCount() == 0

        o = example01()
        assert type(o) == example01
        assert example01.getCount() == 1
        o.__destruct__()
        assert example01.getCount() == 0

        class MyClass1(example01):
            def myfunc(self):
                return 1

        o = MyClass1()
        assert type(o) == MyClass1
        assert isinstance(o, example01)
        assert example01.getCount() == 1
        assert o.myfunc() == 1
        o.__destruct__()
        assert example01.getCount() == 0

        class MyClass2(example01):
            def __init__(self, what):
                example01.__init__(self)
                self.what = what

        o = MyClass2('hi')
        assert type(o) == MyClass2
        assert example01.getCount() == 1
        assert o.what == 'hi'
        o.__destruct__()

        assert example01.getCount() == 0

    def test17_chaining(self):
        """Respective return values of temporaries should not go away"""

        import cppyy

        cppyy.cppdef("""namespace Lifeline {
        struct A1 { A1(int x) : x(x) {} int x; };
        struct A2 { A2(int x) { v.emplace_back(x); } std::vector<A1> v; std::vector<A1>& get() { return v; } };
        struct A3 { A3(int x) { v.emplace_back(x); } std::vector<A2> v; std::vector<A2>& get() { return v; } };
        struct A4 { A4(int x) { v.emplace_back(x); } std::vector<A3> v; std::vector<A3>& get() { return v; } };
        struct A5 { A5(int x) { v.emplace_back(x); } std::vector<A4> v; std::vector<A4>& get() { return v; } };

        A5 gime(int i) { return A5(i); }
        }""")

        assert cppyy.gbl.Lifeline.gime(42).get()[0].get()[0].get()[0].get()[0].x == 42

    def test18_keywords(self):
        """Use of keyword arguments"""

        import cppyy

        cppyy.cppdef("""namespace KeyWords {
        struct A {
 	    A(std::initializer_list<int> vals) : fVals(vals) {}
            std::vector<int> fVals;
        };

        struct B {
	    B() = default;
	    B(const A& in_A, const A& out_A) : fVal(42), fIn(in_A), fOut(out_A) {}
            B(int val, const A& in_A, const A& out_A) : fVal(val), fIn(in_A), fOut(out_A) {}
            int fVal;
            A fIn, fOut;
        };

        int callme(int choice, int a, int b, int c) {
            if (choice == 0) return a;
            if (choice == 1) return b;
            return c;
        }

        struct C {
            int fChoice;
        };

        int callme_c(const C& o, int a, int b, int c) {
            return callme(o.fChoice, a, b, c);
        } }""")

      # constructor and implicit conversion with keywords
        A = cppyy.gbl.KeyWords.A
        B = cppyy.gbl.KeyWords.B

        def verify_b(b, val, ti, to):
            assert b.fVal              == val
            assert b.fIn.fVals.size()  == len(ti)
            assert tuple(b.fIn.fVals)  == ti
            assert b.fOut.fVals.size() == len(to)
            assert tuple(b.fOut.fVals) == to

        b = B(in_A=(256,), out_A=(32,))
        verify_b(b, 42, (256,), (32,))

        b = B(out_A=(32,), in_A=(256,))
        verify_b(b, 42, (256,), (32,))

        with raises(TypeError):
            b = B(in_B=(256,), out_A=(32,))

        b = B(17, in_A=(23,), out_A=(78,))
        verify_b(b, 17, (23,), (78,))

        with raises(TypeError):
            b = B(17, val=23, out_A=(78,))

        with raises(TypeError):
            b = B(17, out_A=(78,)) 

      # global function with keywords
        callme = cppyy.gbl.KeyWords.callme
        for i in range(3):
            assert callme(i, a=1, b=2, c=3) == i+1
            assert callme(i, b=2, c=3, a=1) == i+1
            assert callme(i, c=3, a=1, b=2) == i+1

        with raises(TypeError):
            callme(0, a=1, b=2, d=3)

        with raises(TypeError):
            callme(0, 1, a=2, c=3)

        with raises(TypeError):
            callme(0, a=1, b=2)

      # global function as method with keywords
        c = cppyy.gbl.KeyWords.C()
        cppyy.gbl.KeyWords.C.callme = cppyy.gbl.KeyWords.callme_c

        for i in range(3):
            c.fChoice = i
            assert c.callme(a=1, b=2, c=3) == i+1
            assert c.callme(b=2, c=3, a=1) == i+1
            assert c.callme(c=3, a=1, b=2) == i+1

        c.fChoice = 0
        with raises(TypeError):
            c.callme(a=1, b=2, d=3)

        with raises(TypeError):
            c.callme(1, a=2, c=3)

        with raises(TypeError):
            c.callme(a=1, b=2)


class TestPYTHONIFY_UI:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.example01 = cppyy.load_reflection_info(cls.test_dct)

    def test01_pythonizations(self):
        """Test addition of user-defined pythonizations"""

        import cppyy

        def example01a_pythonize(pyclass, pyname):
            if pyname == 'example01a':
                def getitem(self, idx):
                    return self.addDataToInt(idx)
                pyclass.__getitem__ = getitem

        cppyy.py.add_pythonization(example01a_pythonize)

        e = cppyy.gbl.example01a(1)

        assert e[0] == 1
        assert e[1] == 2
        assert e[5] == 6

    def test02_fragile_pythonizations(self):
        """Test pythonizations error reporting"""

        import cppyy

        example01_pythonize = 1
        raises(TypeError, cppyy.py.add_pythonization, example01_pythonize)

    def test03_write_access_to_globals(self):
        """Test overwritability of globals"""

        import cppyy

        oldval = cppyy.gbl.ns_example01.gMyGlobalInt
        assert oldval == 99

        proxy = cppyy.gbl.ns_example01.__class__.__dict__['gMyGlobalInt']
        cppyy.gbl.ns_example01.gMyGlobalInt = 3
        assert proxy.__get__(proxy, None) == 3

        cppyy.gbl.ns_example01.gMyGlobalInt = oldval
