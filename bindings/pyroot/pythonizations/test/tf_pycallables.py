"""
Tests for passing Python callables when constructing TFX classes.

This feature is not implemented by a PyROOT pythonization, but by a converter of
Cppyy that creates a C++ wrapper to invoke the Python callable.
"""

import unittest
import math

import ROOT


def pyf_tf1_identity(x, p):
    return x[0]


def pyf_tf1_params(x, p):
    return p[0] * x[0] + p[1]


class pyf_tf1_callable:
    def __call__(self, x, p):
        return p[0] * x[0] + p[1]


def pyf_tf1_gauss(x, p):
    return p[0] * 1.0 / math.sqrt(2.0 * math.pi * p[2]**2) * math.exp(-(x[0] - p[1])**2 / 2.0 / p[2]**2)


class TF1(unittest.TestCase):
    """
    Test passing Python callables to ROOT::TF1
    """

    def test_identity(self):
        """
        Test simple function without parameters
        """
        f = ROOT.TF1("tf1_identity", pyf_tf1_identity, 0.0, 1.0)
        for x in [0.0, -1.0, 42.0]:
            self.assertEqual(f.Eval(x), x)

    def test_params(self):
        """
        Test function with parameters
        """
        npars = 2
        f = ROOT.TF1("tf1_params", pyf_tf1_params, 0.0, 1.0, npars)
        par1 = 2.0
        par2 = -1.0
        f.SetParameter(0, par1)
        f.SetParameter(1, par2)
        for x in [0.0, -1.0, 42.0]:
            self.assertEqual(f.Eval(x), pyf_tf1_params([x], [par1, par2]))

    def test_callable(self):
        """
        Test function provided as callable
        """
        npars = 2
        pycallable = pyf_tf1_callable()
        f = ROOT.TF1("tf1_callable", pycallable, 0.0, 1.0, npars)
        par1 = 2.0
        par2 = -1.0
        f.SetParameter(0, par1)
        f.SetParameter(1, par2)
        for x in [0.0, -1.0, 42.0]:
            self.assertEqual(f.Eval(x), pycallable([x], [par1, par2]))


    def test_fitgauss(self):
        """
        Test fitting a histogram to a Python function
        """
        # Gaus function
        f = ROOT.TF1("tf1_fitgauss", pyf_tf1_gauss, -4, 4, 3)
        f.SetParameter(0, 10.0) # scale
        f.SetParameter(1, -1.0) # mean
        f.SetParameter(2, 2.0) # standard deviation

        # Sample gauss in histogram
        h = ROOT.TH1F("h", "test", 100, -4, 4)
        h.FillRandom("gaus", 100000)
        h.Scale(1.0 / 100000.0 * 100.0 / 8.0) # Normalize as density

        # Fit to histogram and get parameters
        h.Fit( f, "0Q" )
        scale = f.GetParameter(0)
        mean = f.GetParameter(1)
        std = f.GetParameter(2)

        self.assertAlmostEqual(scale, 1.0, 2)
        self.assertAlmostEqual(mean, 0.0, 2)
        self.assertAlmostEqual(abs(std), 1.0, 2)


def pyf_tf2_params(x, p):
    return p[0] * x[0] + p[1] * x[1] + p[2]


class TF2(unittest.TestCase):
    """
    Test passing Python callables to ROOT::TF2
    """

    def test_params(self):
        """
        Test function with parameters
        """
        npars = 3
        f = ROOT.TF1("tf2_params", pyf_tf2_params, 0.0, 1.0, npars)
        par1 = 2.0
        par2 = -1.0
        par3 = 1.0
        f.SetParameter(0, par1)
        f.SetParameter(1, par2)
        f.SetParameter(2, par3)
        for x in [(0.0, 0.0), (-1.0, 1.0), (42.0, 0.0)]:
            self.assertEqual(f.Eval(*x), pyf_tf2_params(x, [par1, par2, par3]))


def pyf_tf3_params(x, p):
    return p[0] * x[0] + p[1] * x[1] + p[2] * x[2] + p[3]


class TF3(unittest.TestCase):
    """
    Test passing Python callables to ROOT::TF3
    """

    def test_params(self):
        """
        Test function with parameters
        """
        npars = 4
        f = ROOT.TF1("tf2_params", pyf_tf2_params, 0.0, 1.0, npars)
        par1 = 2.0
        par2 = -1.0
        par3 = 1.0
        par4 = 3.0
        f.SetParameter(0, par1)
        f.SetParameter(1, par2)
        f.SetParameter(2, par3)
        f.SetParameter(3, par4)
        for x in [(0.0, 0.0, 0.0), (-1.0, 1.0, 2.0), (42.0, 0.0, -10.0)]:
            self.assertEqual(f.Eval(*x), pyf_tf2_params(x, [par1, par2, par3, par4]))


if __name__ == '__main__':
    unittest.main()

