import py, os, sys
from pytest import raises
from .support import setup_make


class TestBOOSTANY:
    def setup_class(cls):
        cls.disable = False
        import cppyy
        # TODO: better error handling
        cppyy.include('boost/any.hpp')
        if not hasattr(cppyy.gbl, 'boost'):
            cls.disable = True

    def test01_any_class(self):
        """Availability of boost::any"""

        if self.disable:
            import warnings
            warnings.warn("no boost/any.hpp found, skipping test01_any_class")
            return

        import cppyy
        assert cppyy.gbl.boost.any

        from cppyy.gbl import std
        from cppyy.gbl.boost import any

        assert std.list(any)

    def test02_any_usage(self):
        """boost::any assignment and casting"""

        if self.disable:
            import warnings
            warnings.warn("no boost/any.hpp found, skipping test02_any_usage")
            return

        import cppyy
        assert cppyy.gbl.boost

        from cppyy.gbl import std, boost

        val = boost.any()
        # test both by-ref and by rvalue
        v = std.vector[int]()
        val.__assign__(v)
        val.__assign__(std.move(std.vector[int](range(100))))
        assert val.type() == cppyy.typeid(std.vector[int])

        extract = boost.any_cast[std.vector[int]](val)
        assert type(extract) is std.vector[int]
        assert len(extract) == 100
        extract += range(100)
        assert len(extract) == 200

        val.__assign__(std.move(extract))   # move forced
        #assert len(extract) == 0      # not guaranteed by the standard

        # TODO: we hit boost::any_cast<int>(boost::any* operand) instead
        # of the reference version which raises
        boost.any_cast.__useffi__ = False
        try:
          # raises(Exception, boost.any_cast[int], val)
            assert not boost.any_cast[int](val)
        except Exception:
          # getting here is good, too ...
            pass

        extract = boost.any_cast[std.vector[int]](val)
        assert len(extract) == 200
