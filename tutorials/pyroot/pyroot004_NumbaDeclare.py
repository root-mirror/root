## \file
## \ingroup tutorial_pyroot
## \notebook -nodraw
## This tutorial illustrates how PyROOT supports declaring C++ callables from
## Python callables making them, for example, usable with RDataFrame. The feature
## uses the numba Python package for just-in-time compilation of the Python callable
## and supports fundamental types and ROOT::RVec thereof.
##
## \macro_code
## \macro_output
##
## \date March 2020
## \author Stefan Wunsch

import ROOT

# To mark a Python callable to be used from C++, you have to use the decorator
# provided by PyROOT passing the C++ types of the input arguments and the return
# value.
@ROOT.Numba.Declare(['float', 'int'], 'float')
def pypow(x, y):
    return x**y

# The Python callable is now available from C++ in the Numba namespace.
# For example, we can use it from the interpreter.
ROOT.gInterpreter.ProcessLine('cout << "2^3 = " << Numba::pypow(2, 3) << endl;')

# Or we can use the callable as well within a RDataFrame workflow.
data = ROOT.RDataFrame(4).Define('x', '(float)rdfentry_')\
                         .Define('x_pow3', 'Numba::pypow(x, 3)')\
                         .AsNumpy()

print('pypow({}) = {}'.format(data['x'], data['x_pow3']))

# ROOT uses the numba Python package to create C++ functions from python ones.
# We support as input and return types of the callable fundamental types and
# ROOT::RVec thereof. See the following callable computing the power of the
# elements in an array.
@ROOT.Numba.Declare(['RVec<float>', 'int'], 'RVec<float>')
def pypowarray(x, y):
    return x**y

ROOT.gInterpreter.ProcessLine('''
ROOT::RVec<float> x = {0, 1, 2, 3};
cout << "pypowarray(" << x << ") =  " << Numba::pypowarray(x, 3) << endl;
''')
