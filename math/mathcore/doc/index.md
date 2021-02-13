\defgroup MathCore  MathCore
\ingroup Math
\brief The Core Mathematical Library of %ROOT.

**MathCore** provides a collection of functions and C++ classes for HEP numerical computing.
This library provides only the basic mathematical functions and algorithms and not all the
functionality required by the HEP community. More advanced mathematical functionalities is
provided by the \ref MathMore. The current set includes classes and functions for:

*   \ref SpecFunc "Basic special functions" used in HEP like the gamma, beta and error function;
*   \ref StatFunc : mathematical functions used in statistics, such as the probability density
functions and the cumulative distributions functions (lower and upper integral of the pdf's)
for continuous and discrete distributions.;
*   \ref CppFunctions :
    *   \ref GenFunc, including helper class to wrap free (static) and non-static member functions
    *   \ref ParamFunc
*   Numerical algorithms: user classes with (in some cases) basic implementations for:
    *   \ref Integration
    *   \ref Deriv
    *   \ref RootFinders
    *   \ref Min1D and \ref MultiMin
*   \ref Fit :classes for fitting and parameter estimation from a given data set.

Note that in this latest release the \ref GenVector "GenVector" (physics and geometry vectors
for 2,3 and 4 dimensions with their transformations) is not anymore part of MathCore, but is
built as a separate library. MathCore contains instead now classes which were originally part
of _libCore_. These include:

*   **TMath** : namespace with mathematical functions and basic function algorithms.
*   **TComplex**: class for complex numbers.
*   Random classes: base class **TRandom** and the derived classes TRandom1, TRandom2 and
    TRandom3, implementing the pseudo-random number generators.
*   Other classes, such as
   *   TKDTree for partitioning the data using a kd-Tree and TKDTreeBinning for binning data using a kdTree
   *   ROOT::Math::GoFTest for goodness of fit tests
