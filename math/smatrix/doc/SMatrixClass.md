\page SMatrixDoc SMatrix Class Properties

The template ROOT::Math::SMatrix class has 4 template parameters which define, at compile
time, its properties. These are:

*   type of the contained elements, T, for example _float_ or _double_;
*   number of rows;
*   number of columns;
*   representation type (\ref MatRep). This is a class describing the underlined storage
    model of the Matrix. Presently exists only two types of this class:
    1.  ROOT::Math::MatRepStd for a general nrows x ncols matrix. This class is itself a
        template on the contained type T, the number of rows and the number of columns. Its
        data member is an array T[nrows*ncols] containing the matrix data. The data are
        stored in the row-major C convention. For example, for a matrix, M, of size 3x3,
        the data \f$ \left[a_0,a_1,a_2,.......,a_7,a_8 \right] \f$ are stored in the following
        order: \f[ M = \left( \begin{array}{ccc} a_0 & a_1 & a_2 \\ a_3 & a_4 & a_5 \\ a_6 & a_7 & a_8 \end{array} \right) \f]
    2.  ROOT::Math::MatRepSym for a symmetric matrix of size NxN. This class is a template
        on the contained type and on the symmetric matrix size, N. It has as data member an
        array of type T of size N*(N+1)/2, containing the lower diagonal block of the matrix.
        The order follows the lower diagonal block, still in a row-major convention. For
        example for a symmetric 3x3 matrix the order of the 6 elements
        \f$ \left[a_0,a_1.....a_5 \right]\f$ is: \f[ M = \left( \begin{array}{ccc} a_0 & a_1 & a_3 \\ a_1 & a_2 & a_4 \\ a_3 & a_4 & a_5 \end{array} \right) \f]

### Creating a matrix

The following constructors are available to create a matrix:

*   Default constructor for a zero matrix (all elements equal to zero).
*   Constructor of an identity matrix.
*   Copy constructor (and assignment) for a matrix with the same representation, or from a
    different one when possible, for example from a symmetric to a general matrix.
*   Constructor (and assignment) from a matrix expression, like D = A*B + C. Due to the
    expression template technique, no temporary objects are created in this operation. In
    the case of an operation like A = A*B + C, a temporary object is needed and it is created
    automatically to store the intermediary result in order to preserve the validity of
    this operation.
*   Constructor from a generic STL-like iterator copying the data referred by the iterator,
    following its order. It is both possible to specify the _begin_ and _end_ of the iterator
    or the _begin_ and the size. In case of a symmetric matrix, it is required only the
    triangular block and the user can specify whether giving a block representing the lower
    (default case) or the upper diagonal part.
*   Constructor of a symmetric matrix NxN passing a ROOT::Math::SVector with dimension
    N*(N+1)/2 containing the lower (or upper) block data elements.

Here are some examples on how to create a matrix. We use _typedef's_ in the following examples
to avoid the full C++ names for the matrix classes. Notice that for a general matrix the
representation has the default value, ROOT::Math::MatRepStd, and it is not needed to be
specified. Furthermore, for a general square matrix, the number of column may be as well omitted.

~~~ {.cpp}
// typedef definitions used in the following declarations
typedef ROOT::Math::SMatrix<double,3>                                       SMatrix33;
typedef ROOT::Math::SMatrix<double,2>                                       SMatrix22;
typedef ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> >    SMatrixSym3;
typedef ROOT::Math::SVector<double,2>                                       SVector2;
typedef ROOT::Math::SVector<double,3>                                       SVector3;
typedef ROOT::Math::SVector<double,6>                                       SVector6;

SMatrix33   m0;                         // create a zero 3x3 matrix
// create an 3x3 identity matrix
SMatrix33   i = ROOT::Math::SMatrixIdentity();
double   a[9] = {1,2,3,4,5,6,7,8,9};    // input matrix data
SMatrix33   m(a,9);                     // create a matrix using the a[] data
// this will produce the 3x3 matrix
//    (  1    2    3
//       4    5    6
//       7    8    9  )
~~~


Example to create a symmetric matrix from an _std::vector_:

~~~ {.cpp}
std::vector<double> v(6);
for (int i = 0; i<6; ++i) v[i] = double(i+1);
SMatrixSym3  s(v.begin(),v.end())
// this will produce the symmetric  matrix
//    (  1    2    4
//       2    3    5
//       4    5    6  )

// create a a general matrix from a symmetric matrix. The opposite will not compile
SMatrix33    m2 = s;
~~~


Example to create a symmetric matrix from a ROOT::Math::SVector contining the lower/upper data block:

~~~ {.cpp}
ROOT::Math::SVector<double, 6> v(1,2,3,4,5,6);
SMatrixSym3 s1(v);  // lower block (default)
// this will produce the symmetric  matrix
//    (  1    2    4
//       2    3    5
//       4    5    6  )

SMatrixSym3 s2(v,false);  // upper block
// this will produce the symmetric  matrix
//    (  1    2    3
//       2    4    5
//       3    5    6  )
~~~


### Accessing and Setting Methods

The matrix elements can be set using the _operator()(irow,icol)_, where irow and icol are
the row and column indexes or by using the iterator interface. Notice that the indexes start
from zero and not from one as in FORTRAN. All the matrix elements can be set also by using
the ROOT::Math::SetElements function passing a generic iterator.
The elements can be accessed by these same methods and also by using the
ROOT::Math::SMatrix::apply function. The _apply(i)_ function has exactly the same behavior
for general and symmetric matrices, in contrast to the iterator access methods which behave
differently (it follows the data order).

~~~ {.cpp}
SMatrix33   m;
m(0,0)  = 1;                               // set the element in first row and first column
*(m.begin()+1) = 2;                    // set the second element (0,1)
double d[9]={1,2,3,4,5,6,7,8,9};
m.SetElements(d,d+9);                      // set the d[] values in m

double x = m(2,1);                         // return the element in third row and first column
x = m.apply(7);                        // return the 8-th element (row=2,col=1)
x = *(m.begin()+7);                    // return the 8-th element (row=2,col=1)
// symmetric matrices (note the difference in behavior between apply and the iterators)
x = *(m.begin()+4)                     // return the element (row=2,col=1).
x = m.apply(7);                        // returns again the (row=2,col=1) element
~~~


There are methods to place and/or retrieve ROOT::Math::SVector objects as rows or columns
in (from) a matrix. In addition one can put (get) a sub-matrix as another
ROOT::Math::SMatrix object in a matrix. If the size of the the sub-vector or sub-matrix are
larger than the matrix size a static assert ( a compilation error) is produced. The non-const

~~~ {.cpp}


SMatrix33            m;
SVector2       v2(1,2);
// place a vector of size 2 in the first row starting from element (0,1) : m(0,1) = v2[0]
m.Place_in_row(v2,0,1);
// place the vector in the second column from (0,1) : m(0,1) = v2[0]
m.Place in_col(v2,0,1);
SMatrix22           m2;
// place the sub-matrix m2 in m starting from the element (1,1) : m(1,1) = m2(0,0)
m.Place_at(m2,1,1);
SVector3     v3(1,2,3);
// set v3 as the diagonal elements of m  : m(i,i) = v3[i] for i=0,1,2
m.SetDiagonal(v3)
~~~


The const methods retrieving contents (getting slices of a matrix) are:

~~~ {.cpp}
a = {1,2,3,4,5,6,7,8,9};
SMatrix33       m(a,a+9);
SVector3 irow = m.Row(0);            // return as vector the first matrix row
SVector3 jcol = m.Col(1);            // return as vector the second matrix column
// return a slice of the first row from element (0,1) : r2[0] = m(0,1); r2[1] = m(0,2)
SVector2 r2   =  m.SubRow<SVector2> (0,1);
// return a slice of the second column from element (0,1) : c2[0] = m(0,1); c2[1] = m(1,1);
SVector2 c2   =  m.SubCol<SVector2> (1,0);
// return a sub-matrix 2x2 with the upper left corner at the values (1,1)
SMatrix22 subM = m.Sub<SMatrix22>   (1,1);
// return the diagonal element in a SVector
SVector3  diag = m.Diagonal();
// return the upper(lower) block of the matrix m
SVector6 vub = m.UpperBlock();        //  vub = [ 1, 2, 3, 5, 6, 9 ]
SVector6 vlb = m.LowerBlock();        //  vlb = [ 1, 4, 5, 7, 8, 9 ]
~~~


### Linear Algebra Functions

Only limited linear algebra functionality is available for SMatrix. It is possible
for squared matrices NxN, to find the inverse or to calculate the determinant.
Different inversion algorithms are used if the matrix is smaller than 6x6 or if it
is symmetric. In the case of a small matrix, a faster direct inversion is used.
For a large (N > 6) symmetric matrix the Bunch-Kaufman diagonal pivoting method
is used while for a large (N > 6) general matrix an LU factorization is performed
using the same algorithm as in the CERNLIB routine
[dinv](https://cern-tex.web.cern.ch/cern-tex/shortwrupsdir/f010/top.html).

~~~ {.cpp}
//  Invert a NxN matrix. The inverted matrix replace the existing one and returns if the result is successful
bool ret = m.Invert()
// return the inverse matrix of m. If the inversion fails ifail is different than zero
int ifail = 0;
mInv = m.Inverse(ifail);
~~~


The determinant of a square matrix can be obtained as follows:

~~~ {.cpp}
double det;
// calculate the determinant modifying the matrix content. Returns if the calculation was successful
bool ret = m.Det(det);
// calculate the determinant using a temporary matrix but preserving the matrix content
bool ret = n.Det2(det);
~~~


For additional Matrix functionality see the \ref MatVecFunctions page

