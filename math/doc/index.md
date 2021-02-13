\defgroup Math Math
\brief The %ROOT Mathematical Libraries.

They consist of the following components:

- \ref MathCore "MathCore":   a self-consistent minimal set of tools required for the basic numerical computing.
  It provides the major mathematical functions in the namespaces ROOT::Math and TMath,
  classes for random number generators, TRandom, class for complex numbers, TComplex,
  common interfaces for function evaluation and numerical algorithms.
  Basic implementations of some of the numerical algorithms such as integration or derivation, are also provided by MathCore.
  together with the core classes needed to fit any generic data set.

- \ref MathMore "MathMore": a package incorporating advanced numerical functionality and dependent on external libraries like the GNU Scientific Library ([GSL](http://www.gnu.org/software/gsl/)). It complements the MathCore library by providing a more complete sets of special mathematical functions and implementations of the numerical algorithms interfaces defined in MathCore using GSL.

- **Minimization and Fitting Libraries**
 Libraries required for numerical minimization and fitting. The minimization libraries include the numerical methods for solving the fitting problem by finding minimum of multi-dimensional
  function. The current common interface for minimization is the class ROOT::Math::Minimizer and implemented by derived classes in the minimization and fitting libraries. The fitting in %ROOT is
  organized in fitting classes present in MathCore in the (ROOT::Fit namespace) for providing the fitting functionality and the use the minimization libraries via the common interface (ROOT::Math::Minimizer). In detail the minimization libraries, implementing all the new and old minimization interface, include:

   -  \ref MinuitOld "Minuit": library providing via a class TMinuit an implementation of the popular MINUIT minimization package. In addition the library contains also an implementation of the linear fitter (class TLinearFitter), for solving linear least square fits.
   - \subpage Minuit2Page "Minuit2": new object-oriented implementation of MINUIT, with the same minimization algorithms (such as Migrad or Simplex). In addition it provides a new implementation of the Fumili algorithm, a specialized method for finding the minimum of a standard least square or likelihood functions.
   - **Fumili**: library providing the implementation of the original Fumili fitting algorithm (class TFumili).

- **Linear algebra**. Two libraries are contained in %ROOT for describing linear algebra matrices and vector classes:

   - Matrix: general matrix package providing matrix classes (TMatrixD and TMatrixF)  and vector classes (TVectorD and TVectorF) and the complete environment to perform linear algebra calculations, like equation solving and eigenvalue decompositions.
   - \subpage SMatrixPage "SMatrix": package optimized for high performances matrix and vector computations of small and fixed size. It is based on expression templates to achieve an high level optimization.


- **Physics Vectors**: Classes for describing vectors in 2, 3 and 4 dimensions (relativistic vectors) and their rotation and transformation algorithms. Two package exist in %ROOT:

   - Physics: library with the TVector3 and TLorentzVector classes.
   - GenVector: new library providing generic class templates for modeling the vectors. See the \ref GenVector "GenVector" page.

- \ref Unuran "UNURAN": Package with universal algorithms for generating non-uniform pseudo-random numbers, from a large classes of continuous or discrete distributions in one or multi-dimensions.

- **Foam**  Multi-dimensional general purpose Monte Carlo event generator (and integrator). It generates randomly points (vectors) according to an arbitrary probability distribution  in n dimensions.

- **FFTW** Library with implementation of the fast Fourier transform (FFT) using the FFTW package. It requires a previous installation of [FFTW](http://www.fftw.org).

- **MLP** Library with the neural network class, TMultiLayerPerceptron based on the NN algorithm from the mlpfit package.

- **Quadp** Optimization library with linear and quadratic programming methods. It is based on the Matrix package.


Further information is available at the following links:

- [The Math page in the manual](https://root.cern/manual/math)
- [The Linear Algebra section in the manual](https://root.cern/manual/math/#linear-algebra-packages)
- [The Fitting histograms page in the manual](https://root.cern/manual/fitting/)
- [Inventory of Math functions and algorithms] (http://project-mathlibs.web.cern.ch/project-mathlibs/mathTable.html)

