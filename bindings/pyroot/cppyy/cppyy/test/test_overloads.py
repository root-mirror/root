import py, os, sys
from pytest import raises
from .support import setup_make

currpath = py.path.local(__file__).dirpath()
test_dct = str(currpath.join("overloadsDict"))

def setup_module(mod):
    setup_make("overloads")


class TestOVERLOADS:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.overloads = cppyy.load_reflection_info(cls.test_dct)

    def test01_class_based_overloads(self):
        """Test functions overloaded on different C++ clases"""

        import cppyy
        a_overload = cppyy.gbl.a_overload
        b_overload = cppyy.gbl.b_overload
        c_overload = cppyy.gbl.c_overload
        d_overload = cppyy.gbl.d_overload

        ns_a_overload = cppyy.gbl.ns_a_overload
        ns_b_overload = cppyy.gbl.ns_b_overload

        assert c_overload().get_int(a_overload()) == 42
        assert c_overload().get_int(b_overload()) == 13
        assert d_overload().get_int(a_overload()) == 42
        assert d_overload().get_int(b_overload()) == 13

        assert c_overload().get_int(ns_a_overload.a_overload()) ==  88
        assert c_overload().get_int(ns_b_overload.a_overload()) == -33

        assert d_overload().get_int(ns_a_overload.a_overload()) ==  88
        assert d_overload().get_int(ns_b_overload.a_overload()) == -33

    def test02_class_based_overloads_explicit_resolution(self):
        """Test explicitly resolved function overloads"""

        import cppyy
        a_overload = cppyy.gbl.a_overload
        b_overload = cppyy.gbl.b_overload
        c_overload = cppyy.gbl.c_overload
        d_overload = cppyy.gbl.d_overload

        ns_a_overload = cppyy.gbl.ns_a_overload

        c = c_overload()
        raises(TypeError, c.__dispatch__, 'get_int', 12)
        raises(LookupError, c.__dispatch__, 'get_int', 'does_not_exist')
        assert c.__dispatch__('get_int', 'a_overload*')(a_overload())             == 42
        assert c_overload.get_int.__overload__('a_overload*')(c, a_overload())    == 42
        assert c.__dispatch__('get_int', 'b_overload*')(b_overload())             == 13
        assert c_overload.get_int.__overload__('b_overload*')(c, b_overload())    == 13

        assert c_overload().__dispatch__('get_int', 'a_overload*')(a_overload())  == 42
        # TODO: #assert c_overload.__dispatch__('get_int', 'b_overload*')(c, b_overload()) == 13

        d = d_overload()
        assert d.__dispatch__('get_int', 'a_overload*')(a_overload())             == 42
        assert d_overload.get_int.__overload__('a_overload*')(d, a_overload())    == 42
        assert d.__dispatch__('get_int', 'b_overload*')(b_overload())             == 13
        assert d_overload.get_int.__overload__('b_overload*')(d, b_overload())    == 13

        nb = ns_a_overload.b_overload()
        raises(TypeError, nb.f, c_overload())

    def test03_fragile_class_based_overloads(self):
        """Test functions overloaded on void* and non-existing classes"""

        import cppyy
        more_overloads = cppyy.gbl.more_overloads
        aa_ol = cppyy.gbl.aa_ol
        bb_ol = cppyy.gbl.bb_ol
        cc_ol = cppyy.gbl.cc_ol
        dd_ol = cppyy.gbl.dd_ol

        assert more_overloads().call(aa_ol()) == "aa_ol"
        bb = cppyy.bind_object(cppyy.nullptr, bb_ol)
        assert more_overloads().call(bb     ) == "bb_ol"
        assert more_overloads().call(cc_ol()) == "cc_ol"
        dd = cppyy.bind_object(cppyy.nullptr, dd_ol)
        assert more_overloads().call(dd     ) == "dd_ol"

    def test04_fully_fragile_overloads(self):
        """Test that unknown* is preferred over unknown&"""

        import cppyy
        more_overloads2 = cppyy.gbl.more_overloads2

        bb = cppyy.bind_object(cppyy.nullptr, cppyy.gbl.bb_ol)
        assert more_overloads2().call(bb)    == "bb_olptr"

        dd = cppyy.bind_object(cppyy.nullptr, cppyy.gbl.dd_ol)
        assert more_overloads2().call(dd, 1) == "dd_olptr"

    def test05_array_overloads(self):
        """Test functions overloaded on different arrays"""

        import cppyy
        c_overload = cppyy.gbl.c_overload
        d_overload = cppyy.gbl.d_overload

        from array import array

        ai = array('i', [525252])
        assert c_overload().get_int(ai) == 525252
        assert d_overload().get_int(ai) == 525252

        ah = array('h', [25])
        assert c_overload().get_int(ah) == 25
        assert d_overload().get_int(ah) == 25

    def test06_double_int_overloads(self):
        """Test overloads on int/doubles"""

        import cppyy
        more_overloads = cppyy.gbl.more_overloads

        assert more_overloads().call(1)   == "int"
        assert more_overloads().call(1.)  == "double"
        assert more_overloads().call1(1)  == "int"
        assert more_overloads().call1(1.) == "double"

    def test07_mean_overloads(self):
        """Adapted test for array overloading"""

        import cppyy, array
        cmean = cppyy.gbl.calc_mean

        numbers = [8, 2, 4, 2, 4, 2, 4, 4, 1, 5, 6, 3, 7]
        mean, median = 4.0, 4.0

        for l in ['f', 'd', 'i', 'h', 'l']:
            a = array.array(l, numbers)
            assert round(cmean(len(a), a) - mean, 8) == 0
