import py, os, sys
from pytest import raises
from .support import setup_make, pylong, maxvalue

currpath = py.path.local(__file__).dirpath()
test_dct = str(currpath.join("operatorsDict"))

def setup_module(mod):
    setup_make("operators")


class TestOPERATORS:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.operators = cppyy.load_reflection_info(cls.test_dct)
        cls.N = cppyy.gbl.N

    def teardown_method(self, meth):
        import gc
        gc.collect()

    def test01_math_operators(self):
        """Test overloading of math operators"""

        import cppyy
        number = cppyy.gbl.number

        assert (number(20) + number(10)) == number(30)
        assert (number(20) + 10        ) == number(30)
        assert (number(20) - number(10)) == number(10)
        assert (number(20) - 10        ) == number(10)
        assert (number(20) / number(10)) == number(2)
        assert (number(20) / 10        ) == number(2)
        assert (number(20) * number(10)) == number(200)
        assert (number(20) * 10        ) == number(200)
        assert (number(20) % 10        ) == number(0)
        assert (number(20) % number(10)) == number(0)
        assert (number(5)  & number(14)) == number(4)
        assert (number(5)  | number(14)) == number(15)
        assert (number(5)  ^ number(14)) == number(11)
        assert (number(5)  << 2) == number(20)
        assert (number(20) >> 2) == number(5)

    def test02_unary_math_operators(self):
        """Test overloading of unary math operators"""

        import cppyy
        number = cppyy.gbl.number

        n  = number(20)
        n += number(10)
        n -= number(10)
        n *= number(10)
        n /= number(2)
        assert n == number(100)

        nn = -n;
        assert nn == number(-100)

    def test03_comparison_operators(self):
        """Test overloading of comparison operators"""

        import cppyy
        number = cppyy.gbl.number

        assert (number(20) >  number(10)) == True
        assert (number(20) <  number(10)) == False
        assert (number(20) >= number(20)) == True
        assert (number(20) <= number(10)) == False
        assert (number(20) != number(10)) == True
        assert (number(20) == number(10)) == False

    def test04_boolean_operator(self):
        """Test implementation of operator bool"""

        import cppyy
        number = cppyy.gbl.number

        n = number(20)
        assert n

        n = number(0)
        assert not n

    def test05_exact_types(self):
        """Test converter operators of exact types"""

        import cppyy
        gbl = cppyy.gbl

        o = gbl.operator_char_star()
        assert o.m_str == 'operator_char_star'
        assert str(o)  == 'operator_char_star'

        o = gbl.operator_const_char_star()
        assert o.m_str == 'operator_const_char_star'
        assert str(o)  == 'operator_const_char_star'

        o = gbl.operator_int(); o.m_int = -13
        assert o.m_int == -13
        assert int(o)  == -13

        o = gbl.operator_long(); o.m_long = 42
        assert o.m_long  == 42
        assert pylong(o) == 42

        o = gbl.operator_double(); o.m_double = 3.1415
        assert o.m_double == 3.1415
        assert float(o)   == 3.1415

    def test06_approximate_types(self):
        """Test converter operators of approximate types"""

        import cppyy, sys
        gbl = cppyy.gbl

        o = gbl.operator_short(); o.m_short = 256
        assert o.m_short == 256
        assert int(o)    == 256

        o = gbl.operator_unsigned_int(); o.m_uint = 2147483647 + 32
        assert o.m_uint  == 2147483647 + 32
        assert pylong(o) == 2147483647 + 32

        o = gbl.operator_unsigned_long();
        o.m_ulong = maxvalue + 128
        assert o.m_ulong == maxvalue + 128
        assert pylong(o) == maxvalue + 128

        o = gbl.operator_float(); o.m_float = 3.14
        assert round(o.m_float - 3.14, 5) == 0.
        assert round(float(o) - 3.14, 5)  == 0.

    def test07_virtual_operator_eq(self):
        """Test use of virtual bool operator=="""

        import cppyy

        b1  = cppyy.gbl.v_opeq_base(1)
        b1a = cppyy.gbl.v_opeq_base(1)
        b2  = cppyy.gbl.v_opeq_base(2)
        b2a = cppyy.gbl.v_opeq_base(2)

        assert b1 == b1
        assert b1 == b1a
        assert not b1 == b2
        assert not b1 == b2a
        assert b2 == b2
        assert b2 == b2a

        d1  = cppyy.gbl.v_opeq_derived(1)
        d1a = cppyy.gbl.v_opeq_derived(1)
        d2  = cppyy.gbl.v_opeq_derived(2)
        d2a = cppyy.gbl.v_opeq_derived(2)

        # derived operator== returns opposite
        assert not d1 == d1
        assert not d1 == d1a
        assert d1 == d2
        assert d1 == d2a
        assert not d2 == d2
        assert not d2 == d2a

        # the following is a wee bit interesting due to python resolution
        # rules on the one hand, and C++ inheritance on the other: python
        # will never select the derived comparison b/c the call will fail
        # to pass a base through a const derived&
        assert b1 == d1
        assert d1 == b1
        assert not b1 == d2
        assert not d2 == b1
        
    def test08_call_to_getsetitem_mapping(self):
        """Map () to []"""

        import cppyy

        m = cppyy.gbl.YAMatrix1()
        assert m.m_val == 42
        assert m[1,2]  == 42
        assert m(1,2)  == 42
        m[1,2] = 27
        assert m.m_val == 27
        assert m[1,2]  == 27
        assert m(1,2)  == 27

        m = cppyy.gbl.YAMatrix2()
        assert m.m_val == 42
        assert m[1]    == 42
        m[1] = 27
        assert m.m_val == 27
        assert m[1]    == 27

        for cls in [cppyy.gbl.YAMatrix3, cppyy.gbl.YAMatrix4,
                    cppyy.gbl.YAMatrix5, cppyy.gbl.YAMatrix6,
                    cppyy.gbl.YAMatrix7]:
            m = cls()
            assert m.m_val == 42
            assert m[1,2]  == 42
            assert m[1]    == 42
            assert m(1,2)  == 42

            m[1,2]  = 27
            assert m.m_val == 27
            assert m[1,2]  == 27
            assert m[1]    == 27
            assert m(1,2)  == 27

            m[1]    = 83
            assert m.m_val == 83
            assert m[1,2]  == 83
            assert m[1]    == 83
            assert m(1,2)  == 83

            m.m_val = 74
            assert m.m_val == 74
            assert m[1,2]  == 74
            assert m[1]    == 74
            assert m(1,2)  == 74

    def test09_templated_operator(self):
        """Templated operator<()"""

        from cppyy.gbl import TOIClass

        assert (TOIClass() < 1)

    def test10_r_non_associative(self):
        """Use of radd/rmul with non-associative types"""

        import cppyy

        # Note: calls are repeated to test caching, if any

        a = cppyy.gbl.AssocADD(5.)
        assert 5+a == 10.
        assert a+5 == 10.
        assert 5+a == 10.
        assert a+5 == 10.

        a = cppyy.gbl.NonAssocRADD(5.)
        assert 5+a == 10.
        assert 5+a == 10.
        with raises(NotImplementedError):
            v = a+5

        a = cppyy.gbl.AssocMUL(5.)
        assert 2*a == 10.
        assert a*2 == 10.
        assert 2*a == 10.
        assert a*2 == 10.

        m = cppyy.gbl.NonAssocRMUL(5.)
        assert 2*m == 10.
        assert 2*m == 10.
        with raises(NotImplementedError):
            v = m*2

    def test11_overloaded_operators(self):
        """Overloaded operator*/+-"""

        import cppyy

        v = cppyy.gbl.MultiLookup.Vector2(1, 2)
        w = cppyy.gbl.MultiLookup.Vector2(3, 4)

        u = v*2
        assert u.x == 2.
        assert u.y == 4.

        assert v*w == 1*3 + 2*4

        u = v/2
        assert u.x == 0.5
        assert u.y == 1.0

        assert round(v/w - (1./3. + 2./4.), 8) == 0.

        u = v+2
        assert u.x == 3.
        assert u.y == 4.

        assert v+w == 1+3 + 2+4

        u = v-2
        assert u.x == -1.
        assert u.y ==  0.

        assert v-w == 1-3 + 2-4

    def test12_unary_operators(self):
        """Unary operator-+~"""

        import cppyy

        for cls in [cppyy.gbl.SomeGlobalNumber, cppyy.gbl.Unary.SomeNumber]:
            n = cls(42)

            assert (-n).i == -42
            assert (+n).i ==  42
            #assert (~n).i == ~42

    def test13_comma_operator(self):
        """Comma operator"""

        import cppyy

        c = cppyy.gbl.CommaOperator(1)
        assert c.__comma__(2).__comma__(3).fInt == 6
