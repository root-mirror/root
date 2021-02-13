import py, os, sys
from pytest import raises
from .support import setup_make, pylong, pyunicode, IS_WINDOWS

currpath = py.path.local(__file__).dirpath()
test_dct = str(currpath.join("datatypesDict"))

def setup_module(mod):
    setup_make("datatypes")


class TestLOWLEVEL:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.datatypes = cppyy.load_reflection_info(cls.test_dct)
        cls.N = cppyy.gbl.N

    def test01_builtin_cpp_casts(self):
        """C++ casting of builtin types"""

        from cppyy import ll

        for cast in (ll.cast, ll.static_cast):
            assert type(cast[float](1)) == float
            assert cast[float](1) == 1.

            assert type(cast[int](1.1)) == int
            assert cast[int](1.1) == 1

        assert len(ll.reinterpret_cast['int*'](0)) == 0
        raises(ReferenceError, ll.reinterpret_cast['int*'](0).__getitem__, 0)

    def test02_memory(self):
        """Memory allocation and free-ing"""

        import cppyy
        from cppyy import ll

      # regular C malloc/free
        mem = cppyy.gbl.malloc(16)
        cppyy.gbl.free(mem)

      # typed styles
        mem = cppyy.ll.malloc[int](self.N)
        assert len(mem) == self.N
        for i in range(self.N):
            mem[i] = i+1
            assert type(mem[i]) == int
            assert mem[i] == i+1
        cppyy.ll.free(mem)

      # C++ arrays
        mem = cppyy.ll.array_new[int](self.N)
        assert len(mem) == self.N
        for i in range(self.N):
            mem[i] = i+1
            assert type(mem[i]) == int
            assert mem[i] == i+1
        cppyy.ll.array_delete(mem)

    def test03_python_casts(self):
        """Casts to common Python pointer encapsulations"""

        import cppyy, cppyy.ll

        cppyy.cppdef("""namespace pycasts {
        struct SomeObject{};
        uintptr_t get_address(SomeObject* ptr) { return (intptr_t)ptr; }
        uintptr_t get_deref(void* ptr) { return (uintptr_t)(*(void**)ptr); }
        }""")

        from cppyy.gbl import pycasts

        s = pycasts.SomeObject()
        actual = pycasts.get_address(s)

        assert cppyy.ll.addressof(s) == actual
        assert cppyy.ll.as_ctypes(s).value == actual

        ptrptr = cppyy.ll.as_ctypes(s, byref=True)
        assert pycasts.get_deref(ptrptr) == actual

    def test04_array_as_ref(self):
        """Use arrays for pass-by-ref"""

        import cppyy, sys
        from array import array

        ctd = cppyy.gbl.CppyyTestData()

      # boolean type
        b = array('b', [0]); ctd.set_bool_r(b); assert b[0] == True

      # char types (as data)
        c = array('B', [0]); ctd.set_uchar_r(c); assert c[0] == ord('d')

      # integer types
        i = array('h', [0]);     ctd.set_short_r(i);  assert i[0] == -1
        i = array('H', [0]);     ctd.set_ushort_r(i); assert i[0] ==  2
        i = array('i', [0]);     ctd.set_int_r(i);    assert i[0] == -3
        i = array('I', [0]);     ctd.set_uint_r(i);   assert i[0] ==  4
        i = array('l', [0]);     ctd.set_long_r(i);   assert i[0] == -5
        i = array('L', [0]);     ctd.set_ulong_r(i);  assert i[0] ==  6
        if sys.hexversion >= 0x3000000:
            i = array('q', [0]); ctd.set_llong_r(i);  assert i[0] == -7
            i = array('Q', [0]); ctd.set_ullong_r(i); assert i[0] ==  8

      # floating point types
        f = array('f', [0]);     ctd.set_float_r(f);  assert f[0] ==  5.
        f = array('d', [0]);     ctd.set_double_r(f); assert f[0] == -5.

    def test05_ctypes_as_ref_and_ptr(self):
        """Use ctypes for pass-by-ref/ptr"""

        # See:
        #  https://docs.python.org/2/library/ctypes.html#fundamental-data-types
        #
        # ctypes type       C type                                      Python type
        # ------------------------------------------------------------------------------
        # c_bool            _Bool                                       bool (1)
        #
        # c_char            char 1-character                            string
        # c_wchar           wchar_t 1-character                         unicode string
        # c_byte            char                                        int/long
        # c_ubyte           unsigned char                               int/long
        #
        # c_short           short                                       int/long
        # c_ushort          unsigned short                              int/long
        # c_int             int                                         int/long
        # c_uint            unsigned int                                int/long
        # c_long            long                                        int/long
        # c_ulong           unsigned long                               int/long
        # c_longlong        __int64 or long long                        int/long
        # c_ulonglong       unsigned __int64 or unsigned long long      int/long
        #
        # c_float           float                                       float
        # c_double          double                                      float
        # c_longdouble      long double                                 float

        import cppyy, ctypes

        ctd = cppyy.gbl.CppyyTestData()

      ### pass by reference/pointer and set value back

        for e in ['_r', '_p']:
          # boolean type
            b = ctypes.c_bool(False);      getattr(ctd, 'set_bool'+e)(b);     assert b.value == True

          # char types
            if e == '_r':
                c = ctypes.c_char(b'\0');  getattr(ctd, 'set_char'+e)(c);     assert c.value == b'a'
                c = ctypes.c_wchar(u'\0'); getattr(ctd, 'set_wchar'+e)(c);    assert c.value == u'b'
                c = ctypes.c_byte(0);      getattr(ctd, 'set_schar'+e)(c);    assert c.value == ord('c')
            c = ctypes.c_ubyte(0);         getattr(ctd, 'set_uchar'+e)(c);    assert c.value == ord('d')

          # integer types
            i = ctypes.c_short(0);         getattr(ctd, 'set_short'+e)(i);    assert i.value == -1
            i = ctypes.c_ushort(0);        getattr(ctd, 'set_ushort'+e)(i);   assert i.value ==  2
            i = ctypes.c_int(0);           getattr(ctd, 'set_int'+e)(i);      assert i.value == -3
            i = ctypes.c_uint(0);          getattr(ctd, 'set_uint'+e)(i);     assert i.value ==  4
            i = ctypes.c_long(0);          getattr(ctd, 'set_long'+e)(i);     assert i.value == -5
            i = ctypes.c_ulong(0);         getattr(ctd, 'set_ulong'+e)(i);    assert i.value ==  6
            i = ctypes.c_longlong(0);      getattr(ctd, 'set_llong'+e)(i);    assert i.value == -7
            i = ctypes.c_ulonglong(0);     getattr(ctd, 'set_ullong'+e)(i);   assert i.value ==  8

          # floating point types
            f = ctypes.c_float(0);         getattr(ctd, 'set_float'+e)(f);    assert f.value ==  5.
            f = ctypes.c_double(0);        getattr(ctd, 'set_double'+e)(f);   assert f.value == -5.
            f = ctypes.c_longdouble(0);    getattr(ctd, 'set_ldouble'+e)(f);  assert f.value == 10.

      ### pass by pointer and set value back, now using byref (not recommended)

        cb = ctypes.byref

      # boolean type
        b = ctypes.c_bool(False);     ctd.set_bool_p(cb(b));     assert b.value == True

      # char types
        c = ctypes.c_ubyte(0);        ctd.set_uchar_p(cb(c));    assert c.value == ord('d')

      # integer types
        i = ctypes.c_short(0);        ctd.set_short_p(cb(i));    assert i.value == -1
        i = ctypes.c_ushort(0);       ctd.set_ushort_p(cb(i));   assert i.value ==  2
        i = ctypes.c_int(0);          ctd.set_int_p(cb(i));      assert i.value == -3
        i = ctypes.c_uint(0);         ctd.set_uint_p(cb(i));     assert i.value ==  4
        i = ctypes.c_long(0);         ctd.set_long_p(cb(i));     assert i.value == -5
        i = ctypes.c_ulong(0);        ctd.set_ulong_p(cb(i));    assert i.value ==  6
        i = ctypes.c_longlong(0);     ctd.set_llong_p(cb(i));    assert i.value == -7
        i = ctypes.c_ulonglong(0);    ctd.set_ullong_p(cb(i));   assert i.value ==  8

      # floating point types
        f = ctypes.c_float(0);        ctd.set_float_p(cb(f));    assert f.value ==  5.
        f = ctypes.c_double(0);       ctd.set_double_p(cb(f));   assert f.value == -5.

      ### pass by ptr/ptr with allocation (ptr/ptr is ambiguous in it's task, so many
        # types are allowed to pass; this tests allocation into the pointer)

        from ctypes import POINTER

      # boolean type
        b = POINTER(ctypes.c_bool)();     ctd.set_bool_ppa(b);
        assert b[0] == True; assert b[1] == False; assert b[2] == True
        cppyy.ll.array_delete(b)

      # char types
        c = POINTER(ctypes.c_ubyte)();    ctd.set_uchar_ppa(c)
        assert c[0] == ord('k'); assert c[1] == ord('l'); assert c[2] == ord('m')
        cppyy.ll.array_delete(c)

      # integer types
        i = POINTER(ctypes.c_short)();        ctd.set_short_ppa(i)
        assert i[0] ==  -1; assert i[1] ==  -2; assert i[2] ==  -3
        cppyy.ll.array_delete(i)
        i = POINTER(ctypes.c_ushort)();       ctd.set_ushort_ppa(i)
        assert i[0] ==   4; assert i[1] ==   5; assert i[2] ==   6
        cppyy.ll.array_delete(i)
        i = POINTER(ctypes.c_int)();          ctd.set_int_ppa(i)
        assert i[0] ==  -7; assert i[1] ==  -8; assert i[2] ==  -9
        cppyy.ll.array_delete(i)
        i = POINTER(ctypes.c_uint)();         ctd.set_uint_ppa(i)
        assert i[0] ==  10; assert i[1] ==  11; assert i[2] ==  12
        cppyy.ll.array_delete(i)
        i = POINTER(ctypes.c_long)();         ctd.set_long_ppa(i)
        assert i[0] == -13; assert i[1] == -14; assert i[2] == -15
        cppyy.ll.array_delete(i)
        i = POINTER(ctypes.c_ulong)();        ctd.set_ulong_ppa(i)
        assert i[0] ==  16; assert i[1] ==  17; assert i[2] ==  18
        cppyy.ll.array_delete(i)
        i = POINTER(ctypes.c_longlong)();     ctd.set_llong_ppa(i)
        assert i[0] == -19; assert i[1] == -20; assert i[2] == -21
        cppyy.ll.array_delete(i)
        i = POINTER(ctypes.c_ulonglong)();    ctd.set_ullong_ppa(i)
        assert i[0] ==  22; assert i[1] ==  23; assert i[2] ==  24
        cppyy.ll.array_delete(i)

      # floating point types
        f = POINTER(ctypes.c_float)();        ctd.set_float_ppa(f)
        assert f[0] ==   5; assert f[1] ==  10; assert f[2] ==  20
        cppyy.ll.array_delete(f)
        f = POINTER(ctypes.c_double)();       ctd.set_double_ppa(f)
        assert f[0] ==  -5; assert f[1] == -10; assert f[2] == -20
        cppyy.ll.array_delete(f)
        f = POINTER(ctypes.c_longdouble)();   ctd.set_ldouble_ppa(f)
        assert f[0] ==   5; assert f[1] ==  10; assert f[2] ==  20
        cppyy.ll.array_delete(f)

    def test06_ctypes_pointer_types(self):
        """Use ctypes for pass-by-ptr/ptr-ptr"""

        # See:
        #  https://docs.python.org/2/library/ctypes.html#fundamental-data-types
        #
        # ctypes type       C type                                      Python type
        # ------------------------------------------------------------------------------
        # c_char_p          char* (NULL terminated)                     string or None
        # c_wchar_p         wchar_t* (NULL terminated)                  unicode or None
        # c_void_p          void*                                       int/long or None

        import cppyy, ctypes

        ctd = cppyy.gbl.CppyyTestData()

        ptr = ctypes.c_char_p()
        for meth in ['char', 'cchar']:
            val = getattr(ctd, 'set_'+meth+'_ppm')(ptr)
            assert ctd.freeit(ptr) == val

        ptr = ctypes.c_wchar_p()
        for meth in ['wchar', 'cwchar']:
            val = getattr(ctd, 'set_'+meth+'_ppm')(ptr)
            assert ctd.freeit(ptr) == val

        ptr = ctypes.c_void_p()
        val = ctd.set_void_ppm(ptr)
        assert ctd.freeit(ptr) == val

    def test07_ctypes_type_correctness(self):
        """If types don't match with ctypes, expect exceptions"""

        import cppyy, ctypes

        ctd = cppyy.gbl.CppyyTestData()

        meth_types = ['bool', 'double']
        if not IS_WINDOWS:
            meth_types.append('long')

        i = ctypes.c_int(0);
        for ext in ['_r', '_p']:
            for meth in meth_types:
                with raises(TypeError): getattr(ctd, 'set_'+meth+ext)(i)
