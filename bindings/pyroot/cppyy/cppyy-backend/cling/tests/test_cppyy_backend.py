"""
Pytest/nosetest tests.
"""
import json
import logging
import os
import tempfile

from cppyy_backend import _cppyy_generator


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


class TestCppyyGenerator(object):
    """
    Test cppyy_generator.
    """
    @classmethod
    def setup_class(klass):
        pass

    @classmethod
    def teardown_class(klass):
        pass

    def setUp(self):
        '''This method is run once before _each_ test method is executed'''

    def teardown(self):
        '''This method is run once after _each_ test method is executed'''

    def test_generator(self):
        cpp_h = os.path.join(SCRIPT_DIR, "test_cppyy_backend.h")
        mapfile = tempfile.NamedTemporaryFile().name
        #
        # Create mapping.
        #
        logging.info("mapfile is {}".format(mapfile))
        result = _cppyy_generator.main([
            "",
            "-v",
            "--flags=\\-fvisibility=hidden;\\-D__PIC__;\\-Wno-macro-redefined;\\-std=c++14",
            mapfile, cpp_h
        ])
        assert result == 0
        #
        # Read mapping.
        #
        with open(mapfile, "rU") as mapfile:
            mapping = json.load(mapfile)
            assert mapping[0]["children"][0]["name"] == "Outer"
