import unittest

import ROOT


class RooDataHistPlotOn(unittest.TestCase):
    """
    Initially, this was a test for the pythonization that allowed
    RooDataHist to use the overloads of plotOn defined in RooAbsData.
    Currently, such functionality is automatically provided by Cppyy
    and ROOT meta: the overloads obtained with 'using' declarations
    are taken into account when calling a method.
    We keep this test to check that the aforementioned functionality
    works properly in a case that is important for RooFit.
    """

    # Helpers
    def create_hist_and_frame(self):
        # Inspired by the code of rf402_datahandling.py

        x = ROOT.RooRealVar('x', 'x', -10, 10)
        y = ROOT.RooRealVar('y', 'y', 0, 40)
        x.setBins(10)
        y.setBins(10)

        d = ROOT.RooDataSet('d', 'd', ROOT.RooArgSet(x, y))
        for i in range(10):
            x.setVal(i / 2)
            y.setVal(i)
            d.add(ROOT.RooArgSet(x, y))

        dh = ROOT.RooDataHist('dh', 'binned version of d', ROOT.RooArgSet(x, y), d)

        yframe = ROOT.RooPlot('yplot', 'Operations on binned datasets', y, 0, 40, 10)

        return dh, yframe

    # Tests
    def test_overload1(self):
        dh, yframe = self.create_hist_and_frame()

        # Overload in RooDataHist
        # RooPlot* RooDataHist::plotOn(RooPlot* frame, RooAbsData::PlotOpt o)
        res = dh.plotOn(yframe, ROOT.RooAbsData.PlotOpt())
        self.assertEqual(type(res), ROOT.RooPlot)

    def test_overload2(self):
        dh, yframe = self.create_hist_and_frame()

        # Overload taken from RooAbsData
        # RooPlot* RooAbsData::plotOn(RooPlot* frame, const RooCmdArg& arg1 = RooCmdArg::none(),
        # const RooCmdArg& arg2 = RooCmdArg::none(), const RooCmdArg& arg3 = RooCmdArg::none(),
        # const RooCmdArg& arg4 = RooCmdArg::none(), const RooCmdArg& arg5 = RooCmdArg::none(),
        # const RooCmdArg& arg6 = RooCmdArg::none(), const RooCmdArg& arg7 = RooCmdArg::none(),
        # const RooCmdArg& arg8 = RooCmdArg::none())
        res = dh.plotOn(yframe)
        self.assertEqual(type(res), ROOT.RooPlot)

    def test_overload3(self):
        dh, yframe = self.create_hist_and_frame()

        # Overload taken from RooAbsData
        # RooPlot* RooAbsData::plotOn(RooPlot* frame, const RooLinkedList& cmdList)
        res = dh.plotOn(yframe, ROOT.RooLinkedList())
        self.assertEqual(type(res), ROOT.RooPlot)


if __name__ == '__main__':
    unittest.main()
