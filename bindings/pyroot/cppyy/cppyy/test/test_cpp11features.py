import py, os, sys
from pytest import raises
from .support import setup_make

try:
    import __pypy__
    is_pypy = True
except ImportError:
    is_pypy = False


currpath = py.path.local(__file__).dirpath()
test_dct = str(currpath.join("cpp11featuresDict"))

def setup_module(mod):
    setup_make("cpp11features")

class TestCPP11FEATURES:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.cpp11features = cppyy.load_reflection_info(cls.test_dct)

    def test01_shared_ptr(self):
        """Usage and access of std::shared_ptr<>"""

        from cppyy.gbl import TestSharedPtr, create_shared_ptr_instance

      # proper memory accounting
        assert TestSharedPtr.s_counter == 0

        ptr1 = create_shared_ptr_instance()
        assert ptr1
        assert not not ptr1
        assert TestSharedPtr.s_counter == 1

        ptr2 = create_shared_ptr_instance()
        assert ptr2
        assert not not ptr2
        assert TestSharedPtr.s_counter == 2

        del ptr2
        import gc; gc.collect()
        assert TestSharedPtr.s_counter == 1

        del ptr1
        gc.collect()
        assert TestSharedPtr.s_counter == 0

    def test02_shared_ptr_construction(self):
        """Shared pointer ctor is templated, taking special care"""

        from cppyy.gbl import std, TestSharedPtr

      # proper memory accounting
        assert TestSharedPtr.s_counter == 0

        class C(TestSharedPtr):
            pass

        c = C()
        assert TestSharedPtr.s_counter == 1
        c.__python_owns__ = False
        cc = std.shared_ptr[TestSharedPtr](c)
        assert cc.__python_owns__

        del c

        assert cc
        assert TestSharedPtr.s_counter == 1

        del cc

        import gc
        gc.collect()
        assert TestSharedPtr.s_counter == 0

    def test03_shared_ptr_memory_handling(self):
        """Test shared pointer memory ownership"""

        from cppyy.gbl import std, TestSharedPtr

      # proper memory accounting
        assert TestSharedPtr.s_counter == 0

        t = TestSharedPtr()
        assert TestSharedPtr.s_counter == 1
        assert t.__python_owns__

        tt = std.shared_ptr[TestSharedPtr](t)
        assert not t.__python_owns__

        class C(TestSharedPtr):
            pass

        c = C()
        assert TestSharedPtr.s_counter == 2
        assert c.__python_owns__

        cc = std.shared_ptr[TestSharedPtr](c)
        assert not c.__python_owns__

        del cc, tt

        import gc
        gc.collect()
        assert TestSharedPtr.s_counter == 0

    def test04_shared_ptr_passing(self):
        """Ability to pass shared_ptr<Derived> through shared_ptr<Base>"""

        from cppyy.gbl import std, TestSharedPtr, DerivedTestSharedPtr
        from cppyy.gbl import pass_shared_ptr

    # proper memory accounting
        assert TestSharedPtr.s_counter == 0

        dd = std.make_shared[DerivedTestSharedPtr](DerivedTestSharedPtr(24))
        assert TestSharedPtr.s_counter == 1

        assert pass_shared_ptr(dd) == 100

        del dd

        import gc
        gc.collect()
        assert TestSharedPtr.s_counter == 0

    def test05_nullptr(self):
        """Allow the programmer to pass NULL in certain cases"""
      
        import cppyy

      # test existence
        nullptr = cppyy.nullptr
      # assert not hasattr(cppyy.gbl, 'nullptr')

      # usage is tested in datatypes.py:test15_nullptr_passing
 
    def test06_move(self):
        """Move construction, assignment, and methods"""

        import cppyy, gc

        def moveit(T):
            assert T.s_instance_counter == 0

            from cppyy.gbl import std

          # move constructor
            i1 = T()
            assert T.s_move_counter == 0

            i2 = T(i1)  # cctor
            assert T.s_move_counter == 0

            if is_pypy or 0x3000000 <= sys.hexversion:
                i3 = T(std.move(T()))            # can't check ref-count
            else:
                i3 = T(T()) # should call move, not memoized cctor
            assert T.s_move_counter == 1

            i3 = T(std.move(T()))                # both move and ref-count
            assert T.s_move_counter == 2

            i4 = T(std.move(i1))
            assert T.s_move_counter == 3

          # move assignment
            i4.__assign__(i2)
            assert T.s_move_counter == 3

            if is_pypy or 0x3000000 <= sys.hexversion:
                i4.__assign__(std.move(T()))     # can't check ref-count
            else:
                i4.__assign__(T())
            assert T.s_move_counter == 4

            i4.__assign__(std.move(i2))
            assert T.s_move_counter == 5

            del i1, i2, i3, i4
            gc.collect()
            assert T.s_instance_counter == 0

      # order of moving and normal functions are reversed in 1, 2, for
      # overload resolution testing
        moveit(cppyy.gbl.TestMoving1)
        moveit(cppyy.gbl.TestMoving2)

      # implicit conversion and move
        assert cppyy.gbl.TestMoving1.s_instance_counter == 0
        assert cppyy.gbl.TestMoving2.s_instance_counter == 0
        cppyy.gbl.implicit_converion_move(cppyy.gbl.TestMoving1())
        cppyy.gbl.implicit_converion_move(cppyy.gbl.TestMoving2())
        gc.collect()
        assert cppyy.gbl.TestMoving1.s_instance_counter == 0
        assert cppyy.gbl.TestMoving2.s_instance_counter == 0

    def test07_initializer_list(self):
        """Initializer list construction"""

        from cppyy.gbl import std, TestData, TestData2, WithInitList

        v = std.vector[int]((1, 2, 3, 4))
        assert list(v) == [1, 2, 3, 4]

        v = std.vector['double']((1, 2, 3, 4))
        assert list(v) == [1., 2., 3., 4.]

        raises(TypeError, std.vector[int], [1., 2., 3., 4.])

        for cls in [std.vector, WithInitList]:
            for cls_arg in [TestData, TestData2]:
                l = list()
                for i in range(10):
                    l.append(cls_arg(i))

                v = cls[cls_arg](l)
                assert len(v) == len(l)
                for i in range(len(l)):
                    assert v[i].m_int == l[i].m_int

    def test08_lambda_calls(self):
        """Call (global) lambdas"""

        import cppyy

        cppyy.cppdef("auto gMyLambda = [](int a) { return 40 + a; };")

        assert cppyy.gbl.gMyLambda
        assert cppyy.gbl.gMyLambda(2)  == 42
        assert cppyy.gbl.gMyLambda(40) == 80

        if cppyy.gbl.gInterpreter.ProcessLine("__cplusplus;") >= 201402:
            cppyy.cppdef("auto gime_a_lambda1() { return []() { return 42; }; }")
            l1 = cppyy.gbl.gime_a_lambda1()
            assert l1
            assert l1() == 42

            cppyy.cppdef("auto gime_a_lambda2() { int a = 4; return [a](int b) { return 42+a+b; }; }")
            l2 = cppyy.gbl.gime_a_lambda2()
            assert l2
            assert l2(2) == 48

            cppyy.cppdef("auto gime_a_lambda3(int a ) { return [a](int b) { return 42+a+b; }; }")
            l3 = cppyy.gbl.gime_a_lambda3(4)
            assert l3
            assert l3(2) == 48

    def test09_optional(self):
        """Use of optional and nullopt"""

        import cppyy

        if 201703 <= cppyy.gbl.gInterpreter.ProcessLine("__cplusplus;"):
            assert cppyy.gbl.std.optional
            assert cppyy.gbl.std.nullopt

            cppyy.cppdef("""
                enum Enum { A = -1 };
                bool callopt(std::optional<Enum>) { return true; }
            """)

            a = cppyy.gbl.std.optional[cppyy.gbl.Enum]()
            assert cppyy.gbl.callopt(a)

            c = cppyy.gbl.nullopt
            assert cppyy.gbl.callopt(c)

    def test10_chrono(self):
        """Use of chrono and overloaded operator+"""

        import cppyy
        from cppyy.gbl import std

        t = std.chrono.system_clock.now() - std.chrono.seconds(1)
        # following used to fail with compilation error
        t = std.chrono.system_clock.now() + std.chrono.seconds(1)


    def test11_stdfunction(self):
        """Use of std::function with arguments in a namespace"""

        import cppyy
        from cppyy.gbl import FunctionNS, FNTestStruct, FNCreateTestStructFunc

        t = FNTestStruct(42)
        f = FNCreateTestStructFunc()
        assert f(t) == 42

        t = FunctionNS.FNTestStruct(13)
        f = FunctionNS.FNCreateTestStructFunc()
        assert f(t) == 13

      # and for good measure, inline
        cppyy.cppdef("""namespace FunctionNS2 {
            struct FNTestStruct { FNTestStruct(int i) : t(i) {} int t; };
            std::function<int(const FNTestStruct& t)> FNCreateTestStructFunc() { return [](const FNTestStruct& t) { return t.t; }; }
        }""")

        from cppyy.gbl import FunctionNS2

        t = FunctionNS.FNTestStruct(27)
        f = FunctionNS.FNCreateTestStructFunc()
        assert f(t) == 27

    def test12_stdhash(self):
        """Use of std::hash"""

        import cppyy
        from cppyy.gbl import StructWithHash, StructWithoutHash

        for i in range(3):   # to test effect of caching
            swo = StructWithoutHash()
            assert hash(swo) == object.__hash__(swo)
            assert hash(swo) == object.__hash__(swo)

            sw = StructWithHash()
            assert hash(sw)  == 17
            assert hash(sw)  == 17
