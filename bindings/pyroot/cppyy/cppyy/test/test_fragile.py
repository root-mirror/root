import py, os, sys
from pytest import raises
from .support import setup_make, IS_WINDOWS

currpath = py.path.local(__file__).dirpath()
test_dct = str(currpath.join("fragileDict"))

def setup_module(mod):
    setup_make("fragile")


class TestFRAGILE:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.fragile = cppyy.load_reflection_info(cls.test_dct)

    def test01_load_failure(self):
        """Test failure to load dictionary"""

        import cppyy
        raises(RuntimeError, cppyy.load_reflection_info, "does_not_exist")

        try:
            cppyy.load_reflection_info("does_not_exist")
        except RuntimeError as e:
            assert "does_not_exist" in str(e)

    def test02_missing_classes(self):
        """Test (non-)access to missing classes"""

        import cppyy

        raises(AttributeError, getattr, cppyy.gbl, "no_such_class")

        assert cppyy.gbl.fragile is cppyy.gbl.fragile
        assert cppyy.gbl.fragile == cppyy.gbl.fragile
        fragile = cppyy.gbl.fragile

        raises(AttributeError, getattr, fragile, "no_such_class")

        assert fragile.C is fragile.C
        assert fragile.C == fragile.C
        assert fragile.C().check() == ord('C')

        assert fragile.B is fragile.B
        assert fragile.B == fragile.B
        assert fragile.B().check() == ord('B')
        assert not fragile.B().gime_no_such()

        assert fragile.C is fragile.C
        assert fragile.C == fragile.C
        assert fragile.C().check() == ord('C')
        raises(TypeError, fragile.C().use_no_such, None)

    def test03_arguments(self):
        """Test reporting when providing wrong arguments"""

        import cppyy

        assert cppyy.gbl.fragile == cppyy.gbl.fragile
        fragile = cppyy.gbl.fragile

        assert fragile.D == fragile.D
        assert fragile.D().check() == ord('D')

        d = fragile.D()
        raises(TypeError, d.overload, None)
        raises(TypeError, d.overload, None, None, None)

        d.overload('a')
        d.overload(1)

    def test04_unsupported_arguments(self):
        """Test arguments that are yet unsupported"""

        import cppyy

        assert cppyy.gbl.fragile == cppyy.gbl.fragile
        fragile = cppyy.gbl.fragile

        assert fragile.E == fragile.E
        assert fragile.E().check() == ord('E')

        e = fragile.E()
        raises(TypeError, e.overload, None)
        # allowing access to e.m_pp_no_such is debatable, but it provides a raw pointer
        # which may be useful ...
        assert e.m_pp_no_such[0] == 0xdead

    def test05_wrong_arg_addressof(self):
        """Test addressof() error reporting"""

        import cppyy

        assert cppyy.gbl.fragile == cppyy.gbl.fragile
        fragile = cppyy.gbl.fragile

        assert fragile.F == fragile.F
        assert fragile.F().check() == ord('F')

        f = fragile.F()
        o = object()

        cppyy.addressof(f)
        raises(TypeError, cppyy.addressof, o)
        raises(TypeError, cppyy.addressof, 1)

        # regression (m_int is 0 by default, but its address is not)
        assert cppyy.addressof(f, 'm_int')

        # see also test08_void_pointer_passing in test_advancedcpp.py

    def test06_wrong_this(self):
        """Test that using an incorrect self argument raises"""

        import cppyy

        assert cppyy.gbl.fragile == cppyy.gbl.fragile
        fragile = cppyy.gbl.fragile

        a = fragile.A()
        assert fragile.A.check(a) == ord('A')

        b = fragile.B()
        assert fragile.B.check(b) == ord('B')
        raises(TypeError, fragile.A.check, b)
        raises(TypeError, fragile.B.check, a)

        assert not a.gime_null()

        assert isinstance(a.gime_null(), fragile.A)
        raises(ReferenceError, fragile.A.check, a.gime_null())

    def test07_unnamed_enum(self):
        """Test that an unnamed enum does not cause infinite recursion"""

        import cppyy

        assert cppyy.gbl.fragile is cppyy.gbl.fragile
        fragile = cppyy.gbl.fragile
        assert cppyy.gbl.fragile is fragile

        g = fragile.G()

    def test08_unhandled_scoped_datamember(self):
        """Test that an unhandled scoped data member does not cause infinite recursion"""

        import cppyy

        assert cppyy.gbl.fragile is cppyy.gbl.fragile
        fragile = cppyy.gbl.fragile
        assert cppyy.gbl.fragile is fragile

        h = fragile.H()

    def test09_operator_bool(self):
        """Access to global vars with an operator bool() returning False"""

        import cppyy

        i = cppyy.gbl.fragile.I()
        assert not i

        g = cppyy.gbl.fragile.gI
        assert not g

    def test10_documentation(self):
        """Check contents of documentation"""

        import cppyy

        assert cppyy.gbl.fragile == cppyy.gbl.fragile
        fragile = cppyy.gbl.fragile

        d = fragile.D()
        try:
            d.check(None)         # raises TypeError
            assert 0
        except TypeError as e:
            assert "fragile::D::check()" in str(e)
            assert "TypeError: takes at most 0 arguments (1 given)" in str(e)
            assert "TypeError: takes at least 2 arguments (1 given)" in str(e)

        try:
            d.overload(None)      # raises TypeError
            assert 0
        except TypeError as e:
            # TODO: pypy-c does not indicate which argument failed to convert, CPython does
            # likewise there are still minor differences in descriptiveness of messages
            assert "fragile::D::overload()" in str(e)
            assert "TypeError: takes at most 0 arguments (1 given)" in str(e)
            assert "fragile::D::overload(fragile::no_such_class*)" in str(e)
            #assert "no converter available for 'fragile::no_such_class*'" in str(e)
            assert "void fragile::D::overload(char, int i = 0)" in str(e)
            #assert "char or small int type expected" in str(e)
            assert "void fragile::D::overload(int, fragile::no_such_class* p = 0)" in str(e)
            #assert "int/long conversion expects an integer object" in str(e)

        j = fragile.J()
        assert fragile.J.method1.__doc__ == j.method1.__doc__
        assert j.method1.__doc__ == "int fragile::J::method1(int, double)"

        f = fragile.fglobal
        assert f.__doc__ == "void fragile::fglobal(int, double, char)"

        try:
            o = fragile.O()       # raises TypeError
            assert 0
        except TypeError as e:
            assert "cannot instantiate abstract class 'fragile::O'" in str(e)

    def test11_dir(self):
        """Test __dir__ method"""

        import cppyy

        members = dir(cppyy.gbl.fragile)
        assert 'A' in members
        assert 'B' in members
        assert 'C' in members
        assert 'D' in members                # classes

        assert 'nested1' in members          # namespace

        # TODO: think this through ... probably want this, but interferes with
        # the (new) policy of lazy lookups
        #assert 'fglobal' in members          # function
        #assert 'gI'in members                # variable

    def test12_imports(self):
        """Test ability to import from namespace (or fail with ImportError)"""

        import cppyy

        # TODO: namespaces aren't loaded (and thus not added to sys.modules)
        # with just the from ... import statement; actual use is needed
        from cppyy.gbl import fragile

        def fail_import():
            from cppyy.gbl import does_not_exist
        raises(ImportError, fail_import)

        from cppyy.gbl.fragile import A, B, C, D
        assert cppyy.gbl.fragile.A is A
        assert cppyy.gbl.fragile.B is B
        assert cppyy.gbl.fragile.C is C
        assert cppyy.gbl.fragile.D is D

        # according to warnings, can't test "import *" ...

        from cppyy.gbl.fragile import nested1
        assert cppyy.gbl.fragile.nested1 is nested1
        assert nested1.__name__ == 'nested1'
        assert nested1.__module__ == 'cppyy.gbl.fragile'
        assert nested1.__cpp_name__ == 'fragile::nested1'

        from cppyy.gbl.fragile.nested1 import A, nested2
        assert cppyy.gbl.fragile.nested1.A is A
        assert A.__name__ == 'A'
        assert A.__module__ == 'cppyy.gbl.fragile.nested1'
        assert A.__cpp_name__ == 'fragile::nested1::A'
        assert cppyy.gbl.fragile.nested1.nested2 is nested2
        assert nested2.__name__ == 'nested2'
        assert nested2.__module__ == 'cppyy.gbl.fragile.nested1'
        assert nested2.__cpp_name__ == 'fragile::nested1::nested2'

        from cppyy.gbl.fragile.nested1.nested2 import A, nested3
        assert cppyy.gbl.fragile.nested1.nested2.A is A
        assert A.__name__ == 'A'
        assert A.__module__ == 'cppyy.gbl.fragile.nested1.nested2'
        assert A.__cpp_name__ == 'fragile::nested1::nested2::A'
        assert cppyy.gbl.fragile.nested1.nested2.nested3 is nested3
        assert nested3.__name__ == 'nested3'
        assert nested3.__module__ == 'cppyy.gbl.fragile.nested1.nested2'
        assert nested3.__cpp_name__ == 'fragile::nested1::nested2::nested3'

        from cppyy.gbl.fragile.nested1.nested2.nested3 import A
        assert cppyy.gbl.fragile.nested1.nested2.nested3.A is nested3.A
        assert A.__name__ == 'A'
        assert A.__module__ == 'cppyy.gbl.fragile.nested1.nested2.nested3'
        assert A.__cpp_name__ == 'fragile::nested1::nested2::nested3::A'

        # test writability of __module__
        nested3.__module__ = "peanut butter"
        assert nested3.__module__ == "peanut butter"

        # classes in namespace should inherit
        assert A.__module__ == 'peanut butter.nested3'
        assert 'peanut butter' in repr(A)
        assert 'class' in repr(A)
        assert 'peanut butter' in repr(nested3)
        assert 'namespace' in repr(nested3)

        # as should objects
        a = A()
        assert 'peanut butter' in repr(a)
        assert 'object' in repr(a)

    def test13_missing_casts(self):
        """Test proper handling when a hierarchy is not fully available"""

        import cppyy

        k = cppyy.gbl.fragile.K()

        assert k is k.GimeK(False)
        assert k is not k.GimeK(True)

        kd = k.GimeK(True)
        assert kd is k.GimeK(True)
        assert kd is not k.GimeK(False)

        l = k.GimeL()
        assert l is k.GimeL()

    def test14_double_enum_trouble(self):
        """Test a redefinition of enum in a derived class"""

        import cppyy

        M = cppyy.gbl.fragile.M
        N = cppyy.gbl.fragile.N

        assert M.kOnce == N.kOnce
        assert M.kTwice == N.kTwice

    def test15_const_in_name(self):
        """Make sure 'const' is not erased when part of a name"""

        import cppyy

        cppyy.cppdef("""
            struct Some0Class {}        myvar0;
            struct constSome1Class {}   myvar1;
            struct Some2Classconst {}   myvar2;
            struct Some_const_Class3 {} myvar3;
            struct SomeconstClass4 {}   myvar4;
        """)

        assert cppyy.gbl.myvar0
        assert cppyy.gbl.myvar1
        assert cppyy.gbl.myvar2
        assert cppyy.gbl.myvar3
        assert cppyy.gbl.myvar4

    def test16_opaque_handle(self):
        """Support use of opaque handles"""

        import cppyy

        assert cppyy.gbl.fragile.OpaqueType
        assert cppyy.gbl.fragile.OpaqueHandle_t

        handle = cppyy.gbl.fragile.OpaqueHandle_t(0x42)
        assert handle
        assert cppyy.addressof(handle) == 0x42

        raises(TypeError, cppyy.gbl.fragile.OpaqueType)
        assert not 'OpaqueType' in cppyy.gbl.fragile.__dict__

        handle = cppyy.gbl.fragile.OpaqueHandle_t()
        assert not handle

        addr = cppyy.gbl.fragile.create_handle(handle);
        assert addr
        assert not not handle

        assert cppyy.gbl.fragile.destroy_handle(handle, addr);

        # now define OpaqueType
        cppyy.cppdef("namespace fragile { class OpaqueType { public: int m_int; }; }")

        # get fresh (should not have been cached while incomplete)
        o = cppyy.gbl.fragile.OpaqueType()
        assert hasattr(o, 'm_int')

        assert 'OpaqueType' in cppyy.gbl.fragile.__dict__

    def test17_interactive(self):
        """Test the usage of 'from cppyy.interactive import *'"""

        import assert_interactive

    def test18_overload(self):
        """Test usage of __overload__"""

        import cppyy

        cppyy.cppdef("""struct Variable {
            Variable(double lb, double ub, double value, bool binary, bool integer, const string& name) {}
            Variable(int) {}
        };""")

        for sig in ['double, double, double, bool, bool, const string&',
                    'double,double,double,bool,bool,const string&',
                    'double lb, double ub, double value, bool binary, bool integer, const string& name']:
            assert cppyy.gbl.Variable.__init__.__overload__(sig)

    def test19_gbl_contents(self):
        """Assure cppyy.gbl is mostly devoid of ROOT thingies"""


        import cppyy

        dd = dir(cppyy.gbl)

        assert not 'TCanvasImp' in dd
        assert not 'ESysConstants' in dd
        assert not 'kDoRed' in dd


class TestSIGNALS:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.fragile = cppyy.load_reflection_info(cls.test_dct)

    def test01_abortive_signals(self):
        """Conversion from abortive signals to Python exceptions"""

        import cppyy
        import cppyy.ll

        f = cppyy.gbl.fragile

        assert issubclass(cppyy.ll.BusError,               cppyy.ll.FatalError)
        assert issubclass(cppyy.ll.SegmentationViolation,  cppyy.ll.FatalError)
        assert issubclass(cppyy.ll.IllegalInstruction,     cppyy.ll.FatalError)
        assert issubclass(cppyy.ll.AbortSignal,            cppyy.ll.FatalError)

        import os
        os.putenv('CPPYY_CRASH_QUIET', '1')

        with raises((cppyy.ll.SegmentationViolation, cppyy.ll.IllegalInstruction)):
            with cppyy.ll.signals_as_exception():
                f.segfault()

        with raises(cppyy.ll.AbortSignal):
            with cppyy.ll.signals_as_exception():
                f.sigabort()

      # can only recover once from each error on Windows, which is functionally
      # enough, but precludes further testing here
        if not IS_WINDOWS:
            cppyy.ll.set_signals_as_exception(True)
            with raises((cppyy.ll.SegmentationViolation, cppyy.ll.IllegalInstruction)):
                f.segfault()
            with raises(cppyy.ll.AbortSignal):
                f.sigabort()
            cppyy.ll.set_signals_as_exception(False)

            f.segfault.__sig2exc__ = True
            with raises((cppyy.ll.SegmentationViolation, cppyy.ll.IllegalInstruction)):
                f.segfault()

            f.sigabort.__sig2exc__ = True
            with raises(cppyy.ll.AbortSignal):
                f.sigabort()
