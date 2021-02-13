# -*- coding: UTF-8 -*-
import py, os, sys
from pytest import raises
from .support import setup_make, pylong, pyunicode, maxvalue

try:
    import __pypy__
    is_pypy = True
except ImportError:
    is_pypy = False


currpath = py.path.local(__file__).dirpath()
test_dct = str(currpath.join("stltypesDict"))

def setup_module(mod):
    setup_make("stltypes")


# after CPython's Lib/test/seq_tests.py
def iterfunc(seqn):
    """Regular generator"""
    for i in seqn:
        yield i

class Sequence:
    """Sequence using __getitem__"""
    def __init__(self, seqn):
        self.seqn = seqn
    def __getitem__(self, i):
        return self.seqn[i]

class IterFunc:
    """Sequence using iterator protocol"""
    def __init__(self, seqn):
        self.seqn = seqn
        self.i = 0
    def __iter__(self):
        return self
    def __next__(self):
        if self.i >= len(self.seqn): raise StopIteration
        v = self.seqn[self.i]
        self.i += 1
        return v
    next = __next__ # p2.7

class IterGen:
    """Sequence using iterator protocol defined with a generator"""
    def __init__(self, seqn):
        self.seqn = seqn
        self.i = 0
    def __iter__(self):
        for val in self.seqn:
            yield val

class IterNextOnly:
    """Missing __getitem__ and __iter__"""
    def __init__(self, seqn):
        self.seqn = seqn
        self.i = 0
    def __next__(self):
        if self.i >= len(self.seqn): raise StopIteration
        v = self.seqn[self.i]
        self.i += 1
        return v
    next = __next__ # p2.7

class IterNoNext:
    """Iterator missing __next__()"""
    def __init__(self, seqn):
        self.seqn = seqn
        self.i = 0
    def __iter__(self):
        return self

class IterGenExc:
    """Test propagation of exceptions"""
    def __init__(self, seqn):
        self.seqn = seqn
        self.i = 0
    def __iter__(self):
        return self
    def __next__(self):
        3 // 0
    next = __next__ # p2.7

class IterFuncStop:
    """Test immediate stop"""
    def __init__(self, seqn):
        pass
    def __iter__(self):
        return self
    def __next__(self):
        raise StopIteration
    next = __next__ # p2.7

from itertools import chain
def itermulti(seqn):
    """Test multiple tiers of iterators"""
    return chain(map(lambda x:x, iterfunc(IterGen(Sequence(seqn)))))

class LyingTuple(tuple):
    def __iter__(self):
        yield 1

class LyingList(list):
    def __iter__(self):
        yield 1

def constructors_cpython_test(type2test):
    l0 = []
    l1 = [0]
    l2 = [0, 1]

    u = type2test()
    u0 = type2test(l0)
    u1 = type2test(l1)
    u2 = type2test(l2)

    uu = type2test(u)
    uu0 = type2test(u0)
    uu1 = type2test(u1)
    uu2 = type2test(u2)

    v = type2test(tuple(u))
    class OtherSeq:
        def __init__(self, initseq):
            self.__data = initseq
        def __len__(self):
            return len(self.__data)
        def __getitem__(self, i):
            return self.__data[i]
    s = OtherSeq(u0)
    v0 = type2test(s)
    assert len(v0) == len(s)

    # the following does not work for type-checked containers
    #s = "this is also a sequence"
    #vv = type2test(s)
    #assert len(vv) == len(s)

  # Create from various iteratables
    # as above, can not put strings in type-checked containers
    #for s in ("123", "", range(1000), ('do', 1.2), range(2000,2200,5)):
    for s in (range(1000), range(2000,2200,5)):
        for g in (Sequence, IterFunc, IterGen,
                  itermulti, iterfunc):
            assert type2test(g(s)) == type2test(s)
        assert type2test(IterFuncStop(s))  ==  type2test()
        # as above, no strings
        #assert type2test(c for c in "123") == type2test("123")
        raises(TypeError, type2test, IterNextOnly(s))
        raises(TypeError, type2test, IterNoNext(s))
        raises(ZeroDivisionError, type2test, IterGenExc(s))

  # Issue #23757 (in CPython)
    #assert type2test(LyingTuple((2,))) == type2test((1,))
    #assert type2test(LyingList([2]))   == type2test([1])

def getslice_cpython_test(type2test):
    """Detailed slicing tests from CPython"""

    l = [0, 1, 2, 3, 4]
    u = type2test(l)

    assert u[0:0]        == type2test()
    assert u[1:2]        == type2test([1])
    assert u[-2:-1]      == type2test([3])
    assert u[-1000:1000] == u
    assert u[1000:-1000] == type2test([])
    assert u[:]          == u
    assert u[1:None]     == type2test([1, 2, 3, 4])
    assert u[None:3]     == type2test([0, 1, 2])

  # Extended slices
    assert u[::]          == u
    assert u[::2]         == type2test([0, 2, 4])
    assert u[1::2]        == type2test([1, 3])
    assert u[::-1]        == type2test([4, 3, 2, 1, 0])
    assert u[::-2]        == type2test([4, 2, 0])
    assert u[3::-2]       == type2test([3, 1])
    assert u[3:3:-2]      == type2test([])
    assert u[3:2:-2]      == type2test([3])
    assert u[3:1:-2]      == type2test([3])
    assert u[3:0:-2]      == type2test([3, 1])
    assert u[::-100]      == type2test([4])
    assert u[100:-100:]   == type2test([])
    assert u[-100:100:]   == u
    assert u[100:-100:-1] == u[::-1]
    assert u[-100:100:-1] == type2test([])
    assert u[-pylong(100):pylong(100):pylong(2)] == type2test([0, 2, 4])

  # Test extreme cases with long ints
    a = type2test([0,1,2,3,4])
    # the following two fail b/c PySlice_GetIndices succeeds w/o error, while
    # returning an overflown value (list object uses different internal APIs)
    #assert a[ -pow(2,128): 3 ] == type2test([0,1,2])
    #assert a[ 3: pow(2,145) ]  == type2test([3,4])
    assert a[3::maxvalue]      == type2test([3])


class TestSTLVECTOR:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)
        cls.N = cppyy.gbl.N

    def test01_builtin_type_vector_types(self):
        """Test access to std::vector<int>/std::vector<double>"""

        import cppyy

        assert cppyy.gbl.std        is cppyy.gbl.std
        assert cppyy.gbl.std.vector is cppyy.gbl.std.vector

        assert callable(cppyy.gbl.std.vector)

        type_info = (
            ("int",     int),
            ("float",   "float"),
            ("double",  "double"),
        )

        for c_type, p_type in type_info:
            tv1 = getattr(cppyy.gbl.std, 'vector<%s>' % c_type)
            tv2 = cppyy.gbl.std.vector(p_type)
            assert tv1 is tv2
            assert tv1.iterator is cppyy.gbl.std.vector(p_type).iterator

            #----- 
            v = tv1()
            assert not v
            v += range(self.N)
            assert v
            if p_type == int:
                assert v.begin().__eq__(v.begin())
                assert v.begin() == v.begin()
                assert v.end() == v.end()
                assert v.begin() != v.end()
                assert v.end() != v.begin()

            #-----
            for i in range(self.N):
                v[i] = i
                assert v[i] == i
                assert v.at(i) == i

            assert v.size() == self.N
            assert len(v) == self.N
            assert len(v.data()) == self.N

            #-----
            v = tv1()
            for i in range(self.N):
                v.push_back(i)
                assert v.size() == i+1
                assert v.at(i) == i
                assert v[i] == i

            assert v.size() == self.N
            assert len(v) == self.N

    def test02_user_type_vector_type(self):
        """Test access to an std::vector<just_a_class>"""

        import cppyy

        assert cppyy.gbl.std        is cppyy.gbl.std
        assert cppyy.gbl.std.vector is cppyy.gbl.std.vector

        assert callable(cppyy.gbl.std.vector)

        tv1 = getattr(cppyy.gbl.std, 'vector<just_a_class>')
        tv2 = cppyy.gbl.std.vector('just_a_class')
        tv3 = cppyy.gbl.std.vector(cppyy.gbl.just_a_class)

        assert tv1 is tv2
        assert tv2 is tv3

        v = tv3()
        assert hasattr(v, 'size')
        assert hasattr(v, 'push_back')
        assert hasattr(v, '__getitem__')
        assert hasattr(v, 'begin')
        assert hasattr(v, 'end')

        for i in range(self.N):
            v.push_back(cppyy.gbl.just_a_class())
            v[i].m_i = i
            assert v[i].m_i == i

        assert len(v) == self.N
        v.__destruct__()

    def test03_empty_vector_type(self):
        """Test behavior of empty std::vector<int>"""

        import cppyy

        v = cppyy.gbl.std.vector(int)()
        assert not v
        for arg in v:
            pass
        v.__destruct__()

    def test04_vector_iteration(self):
        """Test iteration over an std::vector<int>"""

        import cppyy

        v = cppyy.gbl.std.vector(int)()

        for i in range(self.N):
            v.push_back(i)
            assert v.size() == i+1
            assert v.at(i) == i
            assert v[i] == i

        assert v.size() == self.N
        assert len(v) == self.N

        i = 0
        for arg in v:
            assert arg == i
            i += 1

        assert list(v) == [i for i in range(self.N)]

        v.__destruct__()

    def test05_push_back_iterables_with_iadd(self):
        """Test usage of += of iterable on push_back-able container"""

        import cppyy

        v = cppyy.gbl.std.vector(int)()

        v += [1, 2, 3]
        assert len(v) == 3
        assert v[0] == 1
        assert v[1] == 2
        assert v[2] == 3

        v += (4, 5, 6)
        assert len(v) == 6
        assert v[3] == 4
        assert v[4] == 5
        assert v[5] == 6

        raises(TypeError, v.__iadd__, (7, '8'))  # string shouldn't pass
        assert len(v) == 7   # TODO: decide whether this should roll-back

        v2 = cppyy.gbl.std.vector(int)()
        v2 += [8, 9]
        assert len(v2) == 2
        assert v2[0] == 8
        assert v2[1] == 9

        v += v2
        assert len(v) == 9
        assert v[6] == 7
        assert v[7] == 8
        assert v[8] == 9

    def test06_vector_indexing(self):
        """Test python-style indexing to an std::vector<int>"""

        import cppyy

        v = cppyy.gbl.std.vector(int)()

        for i in range(self.N):
            v.push_back(i)

        with raises(IndexError):
            v[self.N]
        with raises(IndexError):
            v[self.N+1]

        assert v[-1] == self.N-1
        assert v[-2] == self.N-2

        assert len(v[0:0]) == 0
        assert v[1:2][0] == v[1]

        v2 = v[2:-1]
        assert len(v2) == self.N-3     # 2 off from start, 1 from end
        assert v2[0] == v[2]
        assert v2[-1] == v[-2]
        assert v2[self.N-4] == v[-2]

    def test07_vector_bool(self):
        """Usability of std::vector<bool> which can be a specialization"""

        import cppyy

        vb = cppyy.gbl.std.vector(bool)(8)
        assert [x for x in vb] == [False]*8

        vb[0] = True
        assert vb[0]
        vb[-1] = True
        assert vb[7]

        assert [x for x in vb] == [True]+[False]*6+[True]

        assert len(vb[4:8]) == 4
        assert list(vb[4:8]) == [False]*3+[True]

    def test08_vector_enum(self):
        """Usability of std::vector<> of some enums"""

        import cppyy

        # TODO: would like to use cppyy.gbl.VecTestEnum but that's an int
        assert cppyy.gbl.VecTestEnum
        ve = cppyy.gbl.std.vector['VecTestEnum']()
        ve.push_back(cppyy.gbl.EVal1);
        assert ve[0] == 1
        ve[0] = cppyy.gbl.EVal2
        assert ve[0] == 3

        assert cppyy.gbl.VecTestEnumNS.VecTestEnum
        ve = cppyy.gbl.std.vector['VecTestEnumNS::VecTestEnum']()
        ve.push_back(cppyy.gbl.VecTestEnumNS.EVal1);
        assert ve[0] == 5
        ve[0] = cppyy.gbl.VecTestEnumNS.EVal2
        assert ve[0] == 42

    def test09_vector_of_string(self):
        """Adverse effect of implicit conversion on vector<string>"""

        import cppyy

        assert cppyy.gbl.vectest_ol1("")  == 2
        assert cppyy.gbl.vectest_ol1("a") == 2
        assert cppyy.gbl.vectest_ol2("")  == 2
        assert cppyy.gbl.vectest_ol2("a") == 2

        raises(TypeError, cppyy.gbl.std.vector["std::string"], "abc")

    def test10_vector_std_distance(self):
        """Use of std::distance with vector"""

        import cppyy
        from cppyy.gbl import std

        v = std.vector[int]([1, 2, 3])
        assert v.size() == 3
        assert std.distance(v.begin(), v.end()) == v.size()
        assert std.distance[type(v).iterator](v.begin(), v.end()) == v.size()

    def test11_vector_of_pair(self):
        """Use of std::vector<std::pair>"""

        import cppyy

      # after the original bug report
        cppyy.cppdef("""
        class PairVector {
        public:
            std::vector<std::pair<double, double>> vector_pair(const std::vector<std::pair<double, double>>& a) {
                return a;
            }
        };
        """)

        from cppyy.gbl import PairVector
        a = PairVector()
        ll = [[1., 2.], [2., 3.], [3., 4.], [4., 5.]]
        v = a.vector_pair(ll)

        assert len(v) == 4
        i = 0
        for p in v:
            p.first  == ll[i][0]
            p.second == ll[i][1]
            i += 1
        assert i == 4

      # TODO: nicer error handling for the following (current: template compilation failure trying
      # to assign a pair with <double, string> to <double, double>)
        # ll2 = ll[:]
        # ll2[2] = ll[2][:]
        # ll2[2][1] = 'a'
        # v = a.vector_pair(ll2)

        ll3 = ll[:]
        ll3[0] = 'a'
        raises(TypeError, a.vector_pair, ll3)

        ll4 = ll[:]
        ll4[1] = 'a'
        raises(TypeError, a.vector_pair, ll4)

    def test12_vector_lifeline(self):
        """Check lifeline setting on vectors of objects"""

        import cppyy

        cppyy.cppdef("""namespace Lifeline {
        static int count = 0;
        template <typename T>
        struct C {
            C(int x) : x(x) { ++count; }
            C(const C& c) : x(c.x) { ++count; }
            ~C() { --count; }
            int x;
        };
        auto foo() { return std::vector<C<int>>({C<int>(1337)}); }
        }""")

        assert cppyy.gbl.Lifeline.count == 0
        assert not cppyy.gbl.Lifeline.foo()._getitem__unchecked.__set_lifeline__
        assert cppyy.gbl.Lifeline.foo()[0].x == 1337
        raises(IndexError, cppyy.gbl.Lifeline.foo().__getitem__, 1)
        assert cppyy.gbl.Lifeline.foo()._getitem__unchecked.__set_lifeline__

        import gc
        gc.collect()
        assert cppyy.gbl.Lifeline.count == 0

    def test13_vector_smartptr_iteration(self):
        """Iteration over smart pointers"""

        import cppyy

        cppyy.cppdef("""namespace VectorOfShared {
        struct X {
            int y;
            X() : y(0) {}
            X(int y) : y(y) { }
            std::vector<std::shared_ptr<X>> gimeVec() {
                std::vector<std::shared_ptr<X>> result;
                for (int i = 0; i < 10; ++i) {
                    result.push_back(std::make_shared<X>(i));
                }
                return result;
            }
        }; }""")

        test = cppyy.gbl.VectorOfShared.X()
        result = test.gimeVec()
        assert 'shared' in type(result).__cpp_name__
        assert len(result) == 10

        for i in range(len(result)):
            assert result[i].y == i

        i = 0
        for res in result:
            assert res.y == i
            i += 1
        assert i == len(result)

    def test14_vector_of_vector_of_(self):
        """Nested vectors"""

        from cppyy.gbl.std import vector

        vv = vector[vector[int]](((1, 2), [3, 4]))

        assert len(vv) == 2
        assert list(vv[0]) == [1, 2]
        assert vv[0][0] == 1
        assert vv[0][1] == 2
        assert list(vv[1]) == [3, 4]
        assert vv[1][0] == 3
        assert vv[1][1] == 4

    def test15_vector_slicing(self):
        """Advanced test of vector slicing"""

        from cppyy.gbl.std import vector

        l = list(range(10))
        v = vector[int](range(10))

        assert list(v[2:2])    == l[2:2]
        assert list(v[2:2:-1]) == l[2:2:-1]
        assert list(v[2:5])    == l[2:5]
        assert list(v[5:2])    == l[5:2]
        assert list(v[2:5:-1]) == l[2:5:-1]
        assert list(v[5:2:-1]) == l[5:2:-1]
        assert list(v[2:5: 2]) == l[2:5: 2]
        assert list(v[5:2: 2]) == l[5:2: 2]
        assert list(v[2:5:-2]) == l[2:5:-2]
        assert list(v[5:2:-2]) == l[5:2:-2]
        assert list(v[2:5: 7]) == l[2:5: 7]
        assert list(v[5:2: 7]) == l[5:2: 7]
        assert list(v[2:5:-7]) == l[2:5:-7]
        assert list(v[5:2:-7]) == l[5:2:-7]

      # additional test from CPython's test suite
        getslice_cpython_test(vector[int])

    def test16_vector_construction(self):
        """Vector construction following CPython's sequence"""

        import cppyy

        constructors_cpython_test(cppyy.gbl.std.vector[int])


class TestSTLSTRING:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)

    def test01_string_argument_passing(self):
        """Test mapping of python strings and std::[w]string"""

        import cppyy
        std = cppyy.gbl.std

        for stp, pystp in [(std.string, str), (std.wstring, pyunicode)]:
            stringy_class = cppyy.gbl.stringy_class[stp]

            c, s = stringy_class(pystp("")), stp(pystp("test1"))

            # pass through const std::[w]string&
            c.set_string1(s)
            assert c.get_string1() == s

            c.set_string1(pystp("test2"))
            assert c.get_string1() == pystp("test2")

            # pass through std::string (by value)
            s = stp(pystp("test3"))
            c.set_string2(s)
            assert c.get_string1() == s

            c.set_string2(pystp("test4"))
            assert c.get_string1() == pystp("test4")

            # getting through std::[w]string&
            s2 = stp()
            c.get_string2(s2)
            assert s2 == "test4"

            raises(TypeError, c.get_string2, "temp string")

    def test02_string_data_access(self):
        """Test access to std::string object data members"""

        import cppyy
        std = cppyy.gbl.std

        for stp, pystp in [(std.string, str), (std.wstring, pyunicode)]:
            stringy_class = cppyy.gbl.stringy_class[stp]

            c, s = stringy_class(pystp("dummy")), stp(pystp("test string"))

            c.m_string = pystp("another test")
            assert c.m_string == pystp("another test")
            assert pystp(c.m_string) == c.m_string
            assert c.get_string1() == pystp("another test")

            c.m_string = s
            assert pystp(c.m_string) == s
            assert c.m_string == s
            assert c.get_string1() == s

    def test03_string_with_null_character(self):
        """Test that strings with NULL do not get truncated"""

        import cppyy
        std = cppyy.gbl.std
        stringy_class = cppyy.gbl.stringy_class["std::string"]

        t0 = "aap\0noot"
        assert t0 == "aap\0noot"

        c, s = stringy_class(""), std.string(t0)

        c.set_string1(s)
        assert t0 == c.get_string1()
        assert s == c.get_string1()

        assert std.string('ab\0c')       == 'ab\0c'
        assert repr(std.string('ab\0c')) == repr('ab\0c')
        assert str(std.string('ab\0c'))  == 'ab\0c'

    def test04_array_of_strings(self):
        """Access to global arrays of strings"""

        import cppyy

        assert tuple(cppyy.gbl.str_array_1) == ('a', 'b', 'c')
        str_array_2 = cppyy.gbl.str_array_2
        # fix up the size
        str_array_2.size = 4
        assert tuple(str_array_2) == ('d', 'e', 'f', 'g')
        assert tuple(str_array_2) == ('d', 'e', 'f', 'g')

        # multi-dimensional
        vals = ['a', 'b', 'c', 'd', 'e', 'f']
        str_array_3 = cppyy.gbl.str_array_3
        for i in range(3):
            for j in range(2):
                assert str_array_3[i][j] == vals[i*2+j]

        vals = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
                'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p']
        str_array_4 = cppyy.gbl.str_array_4
        for i in range(4):
            for j in range(2):
                for k in range(2):
                    assert str_array_4[i][j][k] == vals[i*4+j*2+k]

    def test05_stlstring_and_unicode(self):
        """Mixing unicode and std::string"""

        import cppyy

        uas = cppyy.gbl.UnicodeAndSTL

        actlen = len(u'ℕ'.encode(encoding='UTF-8'))
        assert uas.get_size('ℕ')    == actlen
        assert uas.get_size_cr('ℕ') == actlen
        assert uas.get_size_cc('ℕ') == actlen

        assert uas.get_size_w('ℕ')   == 1
        assert uas.get_size_wcr('ℕ') == 1

        assert str(uas.get_string('ℕ'))     == 'ℕ'
        assert str(uas.get_string_cr('ℕ'))  == 'ℕ'
        assert str(uas.get_string_cc('ℕ'))  == 'ℕ'

        assert uas.get_string_w('ℕ').encode(encoding='UTF-8')   == 'ℕ'
        assert uas.get_string_wcr('ℕ').encode(encoding='UTF-8') == 'ℕ'

        bval = u'ℕ'.encode(encoding='UTF-8')
        actlen = len(bval)
        assert uas.get_size(bval)    == actlen
        assert uas.get_size_cr(bval) == actlen
        assert uas.get_size_cc(bval) == actlen

        assert str(uas.get_string(bval))    == 'ℕ'
        assert str(uas.get_string_cr(bval)) == 'ℕ'
        assert str(uas.get_string_cc(bval)) == 'ℕ'


class TestSTLLIST:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)
        cls.N = 13

    def test01_builtin_list_type(self):
        """Test access to a list<int>"""

        import cppyy
        from cppyy.gbl import std

        type_info = (
            ("int",     int),
            ("float",   "float"),
            ("double",  "double"),
        )

        for c_type, p_type in type_info:
            tl1 = getattr(std, 'list<%s>' % c_type)
            tl2 = cppyy.gbl.std.list(p_type)
            assert tl1 is tl2
            assert tl1.iterator is cppyy.gbl.std.list(p_type).iterator

            #-----
            a = tl1()
            assert not a
            for i in range(self.N):
                a.push_back(i)
            assert a

            assert len(a) == self.N
            assert 11 < self.N
            assert 11 in a

            #-----
            ll = list(a)
            for i in range(self.N):
                assert ll[i] == i

            for val in a:
                assert ll[ll.index(val)] == val

    def test02_empty_list_type(self):
        """Test behavior of empty list<int>"""

        import cppyy
        from cppyy.gbl import std

        a = std.list(int)()
        assert not a
        for arg in a:
            pass

    def test03_replacement_of_eq(self):
        """A global templated function should work as a method"""

        import cppyy

        cppyy.cppdef("""template<class C>
            bool cont_eq(const typename C::iterator& a, const typename C::iterator& b) {
                return a != b;
            }""")

        a = cppyy.gbl.std.list[int]()
        a.begin().__class__.__eq__ = cppyy.gbl.cont_eq[cppyy.gbl.std.list[int]]
        assert not (a.begin() == a.end())

        a = cppyy.gbl.std.list[float]()
        assert a.begin() == a.end()
        assert not cppyy.gbl.cont_eq[cppyy.gbl.std.list[float]](a.begin(), a.begin())
        a.push_back(1)
        assert     cppyy.gbl.cont_eq[cppyy.gbl.std.list[float]](a.begin(), a.end())


class TestSTLMAP:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)
        cls.N = 13

    def test01_builtin_map_type(self):
        """Test access to a map<int,int>"""

        import cppyy
        std = cppyy.gbl.std

        a = std.map(int, int)()
        for i in range(self.N):
            a[i] = i
            assert a[i] == i

        assert len(a) == self.N

        for key, value in a:
            assert key == value
        assert key   == self.N-1
        assert value == self.N-1

        # add a variation, just in case
        m = std.map(int, int)()
        for i in range(self.N):
            m[i] = i*i
            assert m[i] == i*i

        for key, value in m:
            assert key*key == value
        assert key   == self.N-1
        assert value == (self.N-1)*(self.N-1)

    def test02_keyed_maptype(self):
        """Test access to a map<std::string,int>"""

        import cppyy
        std = cppyy.gbl.std

        a = std.map(std.string, int)()
        assert not a
        for i in range(self.N):
            a[str(i)] = i
            assert a[str(i)] == i
        assert a

        assert len(a) == self.N

    def test03_empty_maptype(self):
        """Test behavior of empty map<int,int>"""

        import cppyy
        std = cppyy.gbl.std

        m = std.map(int, int)()
        assert not m
        for key, value in m:
            pass

    def test04_unsignedvalue_typemap_types(self):
        """Test assignability of maps with unsigned value types"""

        import cppyy, math, sys
        std = cppyy.gbl.std

        mui = std.map(str, 'unsigned int')()
        mui['one'] = 1
        assert mui['one'] == 1
        raises(ValueError, mui.__setitem__, 'minus one', -1)

        # UInt_t is always 32b, maxvalue is sys.maxint/maxsize and follows system int
        maxint32 = int(math.pow(2,31)-1)
        mui['maxint'] = maxint32 + 3
        assert mui['maxint'] == maxint32 + 3

        mul = std.map(str, 'unsigned long')()
        mul['two'] = 2
        assert mul['two'] == 2
        mul['maxint'] = maxvalue + 3
        assert mul['maxint'] == maxvalue + 3

        raises(ValueError, mul.__setitem__, 'minus two', -2)

    def test05_bool_typemap(self):
        """Test mapping of bool type typedefs"""

        import cppyy

        cppyy.cppdef("""
        struct BoolTypeMapTest {
            typedef bool BoolType;
        };
        """)

        bt = cppyy.gbl.BoolTypeMapTest.BoolType

        assert bt.__name__ == 'BoolType'
        assert bt.__cpp_name__ == 'BoolTypeMapTest::BoolType'
        assert bt(1)
        assert bt(1) == True
        assert bt(1) != False
        assert bt(1) is True
        assert bt() == bt(0)
        assert not bt()
        assert bt() == False
        assert bt() != True
        assert bt() is False
        assert str(bt(1)) == 'True'
        assert str(bt(0)) == 'False'

    def test06_STL_like_class_indexing_overloads(self):
        """Test overloading of operator[] in STL like class"""

        import cppyy
        stl_like_class = cppyy.gbl.stl_like_class

        a = stl_like_class(int)()
        assert a["some string" ] == 'string'
        assert a[3.1415] == 'double'

    def test07_STL_like_class_iterators(self):
        """Test the iterator protocol mapping for an STL like class"""

        import cppyy

        a = cppyy.gbl.stl_like_class(int)()
        assert len(a) == 4
        for i, j in enumerate(a):
            assert i == j

        assert i == len(a)-1

        for cls in [cppyy.gbl.stl_like_class2, cppyy.gbl.stl_like_class3]:
            b = cls[float, 2]()
            b[0] = 27; b[1] = 42
            limit = len(b)+1
            for x in b:
                limit -= 1
                assert limit and "iterated too far!"
                assert x in [27, 42]
            assert x == 42
            del x, b

        for num in [4, 5, 6, 7]:
            cls = getattr(cppyy.gbl, 'stl_like_class%d' % num)
            count = 0
            for i in cls():
                count += 1
            assert count == 10


class TestSTLITERATOR:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)

    def test01_builtin_vector_iterators(self):
        """Test iterator comparison with operator== reflected"""

        import cppyy
        from cppyy.gbl import std

        v = std.vector(int)()
        v.resize(1)

        b1, e1 = v.begin(), v.end()
        b2, e2 = v.begin(), v.end()

        assert b1 == b2
        assert not b1 != b2

        assert e1 == e2
        assert not e1 != e2

        assert not b1 == e1
        assert b1 != e1

        b1.__preinc__()
        assert not b1 == b2
        assert b1 == e2
        assert b1 != b2
        assert b1 == e2


class TestSTLARRAY:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)

    def test01_array_of_basic_types(self):
        """Usage of std::array of basic types"""

        import cppyy
        from cppyy.gbl import std

        a = std.array[int, 4]()
        assert len(a) == 4
        for i in range(len(a)):
            a[i] = i
            assert a[i] == i

    def test02_array_of_pods(self):
        """Usage of std::array of PODs"""

        import cppyy
        from cppyy import gbl
        from cppyy.gbl import std

        a = std.array[gbl.ArrayTest.Point, 4]()
        assert len(a) == 4
        for i in range(len(a)):
            a[i].px = i
            assert a[i].px == i
            a[i].py = i**2
            assert a[i].py == i**2

        if is_pypy:
            raise RuntimeError("test fails with crash")
        # test assignment
        assert a[2]
        a[2] = gbl.ArrayTest.Point(6, 7)
        assert a[2].px == 6
        assert a[2].py == 7

    def test03_array_of_pointer_to_pods(self):
        """Usage of std::array of pointer to PODs"""

        import cppyy
        from cppyy import gbl
        from cppyy.gbl import std

        ll = [gbl.ArrayTest.Point() for i in range(4)]
        for i in range(len(ll)):
            ll[i].px = 13*i
            ll[i].py = 42*i

        a = std.array['ArrayTest::Point*', 4]()
        assert len(a) == 4
        if is_pypy:
            raise RuntimeError("test fails with crash")
        for i in range(len(a)):
            a[i] = ll[i]
            assert a[i].px == 13*i
            assert a[i].py == 42*i

        raises(TypeError, a.__setitem__, 1, 42)

        for i in range(len(a)):
            assert gbl.ArrayTest.get_pp_px(a.data(), i) == 13*i
            assert gbl.ArrayTest.get_pp_py(a.data(), i) == 42*i

            assert gbl.ArrayTest.get_pa_px(a.data(), i) == 13*i
            assert gbl.ArrayTest.get_pa_py(a.data(), i) == 42*i


class TestSTLSTRING_VIEW:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)

    def test01_string_through_string_view(self):
        """Usage of std::string_view as formal argument"""

        import cppyy
        if cppyy.gbl.gInterpreter.ProcessLine("__cplusplus;") <= 201402:
            # string_view exists as of C++17
            return
        countit = cppyy.gbl.StringViewTest.count
        countit_cr = cppyy.gbl.StringViewTest.count_cr

        assert countit("aap")    == 3
        assert countit_cr("aap") == 3
        s = cppyy.gbl.std.string("noot")
        assert countit(s)    == 4
        assert countit_cr(s) == 4
        v = cppyy.gbl.std.string_view(s.data(), s.size())
        assert v[0] == 'n'
        assert countit(v)    == 4
        assert countit_cr(v) == 4


class TestSTLDEQUE:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)
        cls.N = cppyy.gbl.N

    def test01_deque_byvalue_regression(self):
        """Return by value of a deque used to crash"""

        import cppyy
        assert cppyy.cppdef("""std::deque<long double> f() {
            std::deque<long double> d; d.push_back(0); return d ; }""")
        x = cppyy.gbl.f()
        assert x
        del x


class TestSTLSET:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)
        cls.N = cppyy.gbl.N

    def test01_set_iteration(self):
        """Iterate over a set"""

        import cppyy

        s = cppyy.gbl.std.set[int]()
        r = range(self.N)
        for i in r:
            s.insert(i)

        assert len(s) == len(r)
        assert sum(s) == sum(r)

        for i in s:
            assert i in s
            assert i in r

    def test02_set_iterators(self):
        """Access to set iterators and their comparisons"""

        import cppyy

        cppyy.include("iterator")
        s = cppyy.gbl.std.set[int]()

        assert s.begin()  == s.end()
        assert s.rbegin() == s.rend()

        val = 42
        s.insert(val)

        assert len(s) == 1
        assert s.begin().__deref__()  == val
        assert s.rbegin().__deref__() == val

        assert s.begin()  != s.end()
        assert s.begin().__preinc__()  == s.end()
        assert s.rbegin() != s.rend()
        assert s.rbegin().__preinc__() == s.rend()


class TestSTLTUPLE:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)
        cls.N = cppyy.gbl.N

    def test01_tuple_creation_and_access(self):
        """Create tuples and access their elements"""

        import cppyy
        std = cppyy.gbl.std

        t1 = std.make_tuple(1, 'a')
        assert t1
        assert std.get[0](t1) == 1
        assert std.get[1](t1) == 'a'

        t2 = std.make_tuple(1, 'a')
        assert t1 == t2

        t3 = std.make_tuple[int, 'char'](1, 'a')
        assert t3
        assert std.get[0](t3) == 1
        assert std.get[1](t3) == 'a'

        # assert t1 != t3     # fails to link (?!)

        t4 = std.make_tuple(7., 1, 'b')
        assert t4
        assert std.get[0](t4) == 7.
        assert std.get[1](t4) == 1
        assert std.get[2](t4) == 'b'

        v = std.vector[int](range(self.N))
        t5 = std.make_tuple(v, False)
        assert std.get[0](t5).size() == self.N
        assert not std.get[1](t5)

        # TODO: should be easy enough to add iterators over std::tuple?

    def test02_tuple_size(self):
        """Usage of tuple_size helper class"""

        import cppyy
        std = cppyy.gbl.std

        t = std.make_tuple("aap", 42, 5.)
        assert std.tuple_size(type(t)).value == 3

    def test03_tuple_iter(self):
        """Pack/unpack tuples"""

        import cppyy, ctypes
        std = cppyy.gbl.std

        t = std.make_tuple(1, '2', 5.)
        assert len(t) == 3

        a, b, c = t
        assert a == 1
        assert b == '2'
        assert c == 5.


class TestSTLPAIR:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)
        cls.N = cppyy.gbl.N

    def test01_pair_pack_unpack(self):
        """Pack/unpack pairs"""

        import cppyy
        std = cppyy.gbl.std

        p = std.make_pair(1, 2)
        a, b = p

        assert a == 1
        assert b == 2


class TestSTLEXCEPTION:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.stltypes = cppyy.load_reflection_info(cls.test_dct)

    def test01_basics(self):
        """Test behavior of std::exception derived classes"""

        import cppyy

        assert issubclass(cppyy.gbl.std.exception, BaseException)
        assert cppyy.gbl.std.exception is cppyy.gbl.std.exception

        MyError = cppyy.gbl.MyError
        assert MyError is cppyy.gbl.MyError
        assert issubclass(MyError, BaseException)
        assert issubclass(MyError, cppyy.gbl.std.exception)
        assert MyError.__name__     == 'MyError'
        assert MyError.__cpp_name__ == 'MyError'
        assert MyError.__module__   == 'cppyy.gbl'

        YourError = cppyy.gbl.YourError
        assert YourError is cppyy.gbl.YourError
        assert issubclass(YourError, MyError)
        assert YourError.__name__     == 'YourError'
        assert YourError.__cpp_name__ == 'YourError'
        assert YourError.__module__   == 'cppyy.gbl'

        MyError = cppyy.gbl.ErrorNamespace.MyError
        assert MyError is not cppyy.gbl.MyError
        assert MyError is cppyy.gbl.ErrorNamespace.MyError
        assert issubclass(MyError, BaseException)
        assert issubclass(MyError, cppyy.gbl.std.exception)
        assert MyError.__name__     == 'MyError'
        assert MyError.__cpp_name__ == 'ErrorNamespace::MyError'
        assert MyError.__module__   == 'cppyy.gbl.ErrorNamespace'

        YourError = cppyy.gbl.ErrorNamespace.YourError
        assert YourError is not cppyy.gbl.YourError
        assert YourError is cppyy.gbl.ErrorNamespace.YourError
        assert issubclass(YourError, MyError)
        assert YourError.__name__     == 'YourError'
        assert YourError.__cpp_name__ == 'ErrorNamespace::YourError'
        assert YourError.__module__   == 'cppyy.gbl.ErrorNamespace'

    def test02_raising(self):
        """Raise a C++ std::exception derived class as a Python excption"""

        import cppyy

        assert issubclass(cppyy.gbl.MyError, BaseException)

        def raiseit(cls):
            raise cls('Oops')

        with raises(Exception):
            raiseit(cppyy.gbl.MyError)

        with raises(cppyy.gbl.MyError):
            raiseit(cppyy.gbl.MyError)

        try:
            raiseit(cppyy.gbl.MyError)
        except Exception as e:
            assert e.what() == 'Oops'

        try:
            raiseit(cppyy.gbl.MyError)
        except cppyy.gbl.MyError as e:
            assert e.what() == 'Oops'

        try:
            raiseit(cppyy.gbl.YourError)
        except cppyy.gbl.MyError as e:
            assert e.what() == 'Oops'

        try:
            raiseit(cppyy.gbl.YourError)
        except cppyy.gbl.YourError as e:
            assert e.what() == 'Oops'

    def test03_memory(self):
        """Memory handling of C++ c// helper for exception base class testing"""

        import cppyy, gc

        MyError   = cppyy.gbl.MyError
        YourError = cppyy.gbl.YourError

        gc.collect()
        assert cppyy.gbl.GetMyErrorCount() == 0

        m = MyError('Oops')
        assert cppyy.gbl.GetMyErrorCount() == 1
        del m
        gc.collect()
        assert cppyy.gbl.GetMyErrorCount() == 0

        def raiseit(cls):
            raise cls('Oops')

        def run_raiseit(t1, t2):
            try:
                raiseit(t1)
            except t2 as e:
                assert e.what() == 'Oops'
                return
            assert not "should not reach this point"

        for t1, t2 in [(MyError,   Exception),
                       (MyError,   MyError),
                       (YourError, MyError),
                       (YourError, YourError)]:
            with raises(t2):
                raiseit(t1)
            gc.collect()
            assert cppyy.gbl.GetMyErrorCount() == 0

            run_raiseit(t1, t2)
            gc.collect()
            assert cppyy.gbl.GetMyErrorCount() == 0

        gc.collect()
        assert cppyy.gbl.GetMyErrorCount() == 0

    def test04_from_cpp(self):
        """Catch C++ exceptiosn from C++"""

        import cppyy, gc

        gc.collect()
        assert cppyy.gbl.GetMyErrorCount() == 0

        with raises(cppyy.gbl.MyError):
            cppyy.gbl.ErrorNamespace.throw_error(0)

        with raises(cppyy.gbl.MyError):
            cppyy.gbl.ErrorNamespace.throw_error(1)

        with raises(cppyy.gbl.YourError):
            cppyy.gbl.ErrorNamespace.throw_error(1)

        with raises(cppyy.gbl.ErrorNamespace.MyError):
            cppyy.gbl.ErrorNamespace.throw_error(2)

        with raises(cppyy.gbl.ErrorNamespace.MyError):
            cppyy.gbl.ErrorNamespace.throw_error(3)

        with raises(cppyy.gbl.ErrorNamespace.YourError):
            cppyy.gbl.ErrorNamespace.throw_error(3)

        gc.collect()
        assert cppyy.gbl.GetMyErrorCount() == 0
