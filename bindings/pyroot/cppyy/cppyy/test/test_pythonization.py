import py, os, sys
from pytest import raises
from .support import setup_make, pylong

currpath = py.path.local(__file__).dirpath()
test_dct = str(currpath.join("pythonizablesDict"))

def setup_module(mod):
    setup_make("pythonizables")


class TestClassPYTHONIZATION:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.pyzables = cppyy.load_reflection_info(cls.test_dct)

    def test00_api(self):
        """Test basic semantics of the pythonization API"""

        import cppyy

        raises(TypeError, cppyy.py.add_pythonization, 1)

        def pythonizor1(klass, name):
            pass

        def pythonizor2(klass, name):
            pass

        pythonizor3 = pythonizor1

        cppyy.py.add_pythonization(pythonizor1)
        assert cppyy.py.remove_pythonization(pythonizor2) == False
        assert cppyy.py.remove_pythonization(pythonizor3) == True

        def pythonizor(klass, name):
            if name == 'pyzables::SomeDummy1':
                klass.test = 1

        cppyy.py.add_pythonization(pythonizor)
        assert cppyy.gbl.pyzables.SomeDummy1.test == 1

        def pythonizor(klass, name):
            if name == 'SomeDummy2':
                klass.test = 2
        cppyy.py.add_pythonization(pythonizor, 'pyzables')

        def pythonizor(klass, name):
            if name == 'pyzables::SomeDummy2':
                klass.test = 3
        cppyy.py.add_pythonization(pythonizor)

        assert cppyy.gbl.pyzables.SomeDummy2.test == 2

        def root_pythonizor(klass, name):
            if name == 'TString':
                klass.__len__ = klass.Length

        cppyy.py.add_pythonization(root_pythonizor)

        assert len(cppyy.gbl.TString("aap")) == 3

    def test01_size_mapping(self):
        """Use composites to map GetSize() onto buffer returns"""

        import cppyy

        def set_size(self, buf):
            buf.reshape((self.GetN(),))
            return buf

        cppyy.py.add_pythonization(
            cppyy.py.compose_method('NakedBuffers$', 'Get[XY]$', set_size), 'pyzables')

        bsize, xval, yval = 3, 2, 5
        m = cppyy.gbl.pyzables.NakedBuffers(bsize, xval, yval)

        x = m.GetX()
        assert len(x) == bsize
        assert list(x) == list(map(lambda x: x*xval, range(bsize)))

        y = m.GetY()
        assert len(y) == bsize
        assert list(y) == list(map(lambda x: x*yval, range(bsize)))

    def test02_size_mapping_of_templated_method(self):
        """Use composites to map GetSize() onto buffer returns"""

        import cppyy

        def set_size(self, buf):
            buf.reshape((self.GetN(),))
            return buf

        cppyy.py.add_pythonization(
            cppyy.py.compose_method('NakedBuffers2.*Vector.*', 'Get[XY]$', set_size), 'pyzables')

        bsize, xval, yval = 3, 2, 5
        m = cppyy.gbl.pyzables.NakedBuffers2[cppyy.gbl.pyzables.Vector](bsize, xval, yval)

        x = m.GetX()
        assert len(x) == bsize
        assert list(x) == list(map(lambda x: x*xval, range(bsize)))

        y = m.GetY()
        assert len(y) == bsize
        assert list(y) == list(map(lambda x: x*yval, range(bsize)))

    def test03_type_pinning(self):
        """Verify pinnability of returns"""

        import cppyy

        cppyy.gbl.pyzables.GimeDerived.__creates__ = True

        result = cppyy.gbl.pyzables.GimeDerived()
        assert type(result) == cppyy.gbl.pyzables.MyDerived

        cppyy.py.pin_type(cppyy.gbl.pyzables.MyBase)
        assert type(result) == cppyy.gbl.pyzables.MyDerived


    def test04_transparency(self):
        """Transparent use of smart pointers"""

        import cppyy

        Countable = cppyy.gbl.pyzables.Countable
        mine = cppyy.gbl.pyzables.mine

        assert type(mine) == Countable
        assert mine.m_check == 0xcdcdcdcd
        assert type(mine.__smartptr__()) == cppyy.gbl.std.shared_ptr(Countable)
        assert mine.__smartptr__().get().m_check == 0xcdcdcdcd
        assert mine.say_hi() == "Hi!"

    def test05_converters(self):
        """Smart pointer argument passing"""

        import cppyy

        pz = cppyy.gbl.pyzables
        mine = pz.mine

        assert 0xcdcdcdcd == pz.pass_mine_rp_ptr(mine)
        assert 0xcdcdcdcd == pz.pass_mine_rp_ref(mine)
        assert 0xcdcdcdcd == pz.pass_mine_rp(mine)

        assert 0xcdcdcdcd == pz.pass_mine_sp_ptr(mine)
        assert 0xcdcdcdcd == pz.pass_mine_sp_ref(mine)

        assert 0xcdcdcdcd == pz.pass_mine_sp_ptr(mine.__smartptr__())
        assert 0xcdcdcdcd == pz.pass_mine_sp_ref(mine.__smartptr__())

        assert 0xcdcdcdcd == pz.pass_mine_sp(mine)
        assert 0xcdcdcdcd == pz.pass_mine_sp(mine.__smartptr__())

        # TODO:
        # cppyy.gbl.mine = mine
        pz.renew_mine()

    def test06_executors(self):
        """Smart pointer return types"""

        import cppyy

        pz = cppyy.gbl.pyzables
        Countable = pz.Countable

        mine = pz.gime_mine_ptr()
        assert type(mine) == Countable
        assert mine.m_check == 0xcdcdcdcd
        assert type(mine.__smartptr__()) == cppyy.gbl.std.shared_ptr(Countable)
        assert mine.__smartptr__().get().m_check == 0xcdcdcdcd
        assert mine.say_hi() == "Hi!"

        mine = pz.gime_mine_ref()
        assert type(mine) == Countable
        assert mine.m_check == 0xcdcdcdcd
        assert type(mine.__smartptr__()) == cppyy.gbl.std.shared_ptr(Countable)
        assert mine.__smartptr__().get().m_check == 0xcdcdcdcd
        assert mine.say_hi() == "Hi!"

        mine = pz.gime_mine()
        assert type(mine) == Countable
        assert mine.m_check == 0xcdcdcdcd
        assert type(mine.__smartptr__()) == cppyy.gbl.std.shared_ptr(Countable)
        assert mine.__smartptr__().get().m_check == 0xcdcdcdcd
        assert mine.say_hi() == "Hi!"

    def test07_creates_flag(self):
        """Effect of creates flag on return type"""

        import cppyy, gc

        pz = cppyy.gbl.pyzables
        Countable = pz.Countable

        gc.collect()
        oldcount = Countable.sInstances     # there's eg. one global variable

        pz.gime_naked_countable.__creates__ = True
        for i in range(10):
            cnt = pz.gime_naked_countable()
            gc.collect()
            assert Countable.sInstances == oldcount + 1
        del cnt
        gc.collect()

        assert Countable.sInstances == oldcount

    def test08_base_class_pythonization(self):
        """Derived class should not re-pythonize base class pythonization"""

        import cppyy

        d = cppyy.gbl.pyzables.IndexableDerived()

        assert d[0]  == 42
        assert d[-1] == 42
        # skip the IndexErorr test: pythonization for __getitem__[index] < size()
        # can not be applied strict enough (instead of an index, this could be an
        # associative  container, with 'index' a key, not a counter
        #raises(IndexError, d.__getitem__, 1)

    def test09_cpp_side_pythonization(self):
        """Use of C++ side pythonizations"""

        import cppyy

      # explicit pythonization
        for kls in [cppyy.gbl.pyzables.WithCallback1, cppyy.gbl.pyzables.WithCallback2]:
            w = kls(42)
            assert hasattr(w, 'GetInt')
            assert not hasattr(w, 'get_int')
            assert w.GetInt() == 42

            assert hasattr(w, 'SetInt')
            assert not hasattr(w, 'set_int')
            w.SetInt(17)
            assert w.GetInt() == 17

            assert kls.klass_name == kls.__cpp_name__

      # up-the-hierarchy pythonization
        w = cppyy.gbl.pyzables.WithCallback3(42)
        assert hasattr(w, 'GetInt')
        assert not hasattr(w, 'get_int')
        assert w.GetInt() == 2*42

        assert hasattr(w, 'SetInt')
        assert not hasattr(w, 'set_int')
        w.SetInt(17)
        assert w.GetInt() == 4*17

        assert cppyy.gbl.pyzables.WithCallback2.klass_name == 'pyzables::WithCallback3'


## actual test run
if __name__ == '__main__':
    result = run_pytest(__file__)
    sys.exit(result)
