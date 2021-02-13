import py, os, sys
from pytest import raises
from .support import setup_make

currpath = py.path.local(__file__).dirpath()
test_dct = str(currpath.join("std_streamsDict"))

def setup_module(mod):
    setup_make("std_streams")


class TestSTDStreams:
    def setup_class(cls):
        cls.test_dct = test_dct
        import cppyy
        cls.streams = cppyy.load_reflection_info(cls.test_dct)

    def test01_std_ostream(self):
        """Test availability of std::ostream"""

        import cppyy

        assert cppyy.gbl.std is cppyy.gbl.std
        assert cppyy.gbl.std.ostream is cppyy.gbl.std.ostream

        assert callable(cppyy.gbl.std.ostream)

    def test02_std_cout(self):
        """Test access to std::cout"""

        import cppyy

        assert not (cppyy.gbl.std.cout is None)
