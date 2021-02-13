\defgroup SMatrixGroup SMatrix Package
\ingroup  Math


**SMatrix** is a C++ package for high performance vector and matrix computations. It can be
used only in problems when the size of the matrices is known at compile time, like in the
tracking reconstruction of HEP experiments. It is based on a C++ technique, called expression
templates, to achieve an high level optimization. The C++ templates can be used to implement
vector and matrix expressions such that these expressions can be transformed at compile time
to code which is equivalent to hand optimized code in a low-level language like FORTRAN or
C (see for example ref. 1)

The SMatrix has been developed initially by T. Glebe of the Max-Planck-Institut, Heidelberg,
as part of the HeraB analysis framework. A subset of the original package has been now incorporated
in the %ROOT distribution, with the aim to provide to the LHC experiments a stand-alone and
high performant matrix package for reconstruction. The API of the current package differs
from the original one, in order to be compliant to the %ROOT coding conventions.

SMatrix contains generic \ref SMatrixSVector to describe matrix and vector of arbitrary
dimensions and of arbitrary type. The classes are templated on the scalar type and on the
size of the matrix (number of rows and columns) or the vector. Therefore, the size has to
be known at compile time. Since the release 5.10, SMatrix supports symmetric matrices using
a storage class (ROOT::Math::MatRepSym) which contains only the N*(N+1)/2 independent element
of a NxN symmetric matrix.
It is not in the mandate of this package to provide a complete linear algebra functionality
for these classes. What is provided are basic \ref MatrixFunctions and \ref VectFunction,
such as the matrix-matrix, matrix-vector, vector-vector operations, plus some extra
functionality for square matrices, like inversion, which is based on the optimized Cramer
method for squared matrices of size up to 6x6, and determinant calculation.
For a more detailed descriptions and usage examples see:

*   \ref SVectorDoc
*   \ref SMatrixDoc
*   \ref MatVecFunctions

The SMatrix package contains only header files. Normally one does not need to build any library.
In the %ROOT distribution a library, _libSmatrix_ is produced with the C++ dictionary information
for vectors, symmetric and squared matrices for double, float types up to dimension 7.
The current version of SMatrix can be downloaded from [here](../SMatrix.tar.gz). If you want
to install the header files or run the test _configure_ script and then _make install_ or
_make check_ to build the tests. No dictionary library is built in this case.

## References

1.  T. Veldhuizen, [_Expression Templates_](http://osl.iu.edu/~tveldhui/papers/Expression-Templates/exprtmpl.html),
    C++ Report, 1995.
2.  T. Glebe, _SMatrix - A high performance library for Vector/Matrix calculation and Vertexing_,
    HERA-B Software Note 01-134, December 2, 2003 ([pdf](http://seal.web.cern.ch/seal/documents/mathlib/smatrix_herab.pdf))
3.  L. Moneta, %ROOT Math proposal for Linear Algebra, [presentation](http://seal.cern.ch/documents/mathlib/aa_matrix_nov05.pdf)
    at the LCG Application Area meeting, November 23, 2005

* * *

@authors the %ROOT Math Library Team, T. Glebe (original SMatrix author) and J. Palacios (LHCb)
