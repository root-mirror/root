import unittest


class ROOTModule(unittest.TestCase):
    """
    Testing features of the ROOT module implemented in the ROOT module facade
    """

    def test_import(self):
        """
        Test import
        """
        import ROOT
        self.assertEqual(ROOT.__name__, "ROOT")


    def test_relative_import(self):
        """
        Test relative import
        """
        from ROOT import TH1F
        self.assertEqual(TH1F.__name__, "TH1F")
        from ROOT import TH1F as Foo
        self.assertEqual(Foo.__name__, "TH1F")
        self.assertEqual(TH1F, Foo)


    def test_version(self):
        """
        Test __version__ property
        """
        import ROOT
        v = ROOT.__version__
        self.assertTrue(type(v) == str)
        self.assertIn("/", v)
        self.assertIn(".", v)
        self.assertEqual(v, ROOT.gROOT.GetVersion())


    def test_ignore_cmdline_options(self):
        """
        Test module flag to ignore command line options
        """
        import ROOT
        self.assertEqual(ROOT.PyConfig.IgnoreCommandLineOptions, True)
        ROOT.PyConfig.IgnoreCommandLineOptions = False
