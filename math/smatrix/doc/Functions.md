\page MatVecFunctions Matrix and Vector Operators and Functions


## Matrix and Vector Operators

The ROOT::Math::SVector and ROOT::Math::SMatrix classes defines the following operators
described below. The _m1,m2,m3_ are vectors or matrices of the same type (and size)
and _a_ is a scalar value:

~~~ {.cpp}
m1 == m2           // returns whether m1 is equal to m2 (element by element comparison)
m1 != m2           // returns whether m1 is NOT equal to m2 (element by element comparison)
m1 < m2            // returns whether m1 is less than m2 (element wise comparison)
m1 > m2            // returns whether m1 is greater than m2 (element wise comparison)
// in the following m1 and m3 can be general and m2 symmetric, but not vice-versa
m1 += m2           // add m2 to m1
m1 -= m2           // subtract m2 to m1
m3 = m1 + m2       // addition
m1 - m2            // subtraction

// Multiplication and division via a scalar value a
m3 = a*m1; m3 = m1*a; m3 = m1/a;
~~~


### Vector-Vector multiplication

The _operator *_ defines an element by element multiplication between vectors. For the
standard vector-vector multiplication, \f$ a = v^T v \f$, (dot product) one must use the
ROOT::Math::Dot function. In addition, the Cross (only for vector sizes of 3),
ROOT::Math::Cross, and the Tensor product, ROOT::Math::TensorProd, are defined.

### Matrix - Vector multiplication

The _operator *_ defines the matrix-vector multiplication, \f$ y_i = \sum_{j} M_{ij} x_j\f$:

~~~ {.cpp}
// M is a  N1xN2 matrix, x is a N2 size vector, y is a N1 size vector
y = M * x
~~~


It compiles only if the matrix and the vectors have the right sizes.

### Matrix - Matrix multiplication

The _operator *_ defines the matrix-matrix multiplication, \f$ C_{ij} = \sum_{k} A_{ik} B_{kj}\f$:

~~~ {.cpp}
// A is a N1xN2 matrix, B is a N2xN3 matrix and C is a N1xN3 matrix
C = A * B
~~~


The operation compiles only if the matrices have the right size. In the case that A and B
are symmetric matrices, C is a general one, since their product is not guaranteed to be symmetric.

#### Special note on using the C++ auto keyword

Special care must be taken when using the C++ ``auto`` keyword with expression templates.
Some expression can lead to temporary objects that the compiler might remove. One example
is when dealing with an expression like:

~~~ {.cpp}
auto D = (A * B) * C;
~~~

while instead declaring directly the matrix as

~~~ {.cpp}
SMatrix<double, N, N> D = (A * B) * C;
~~~

will be fine, because it will force the evaluation of the expression template when
constructing the matrix D.
This is a limitation of the package, see [ROOT-6731](https://sft.its.cern.ch/jira/browse/ROOT-6371)
and present in other similar libraries such as Eigen.

### Matrix and Vector Functions

The most used matrix functions are:

*   **ROOT::Math::Transpose**(M) : return the transpose matrix, \f$ M^T \f$
*   **ROOT::Math::Similarity**( v, M) : returns the scalar value resulting from the matrix- vector product \f$ v^T M v \f$
*   **ROOT::Math::Similarity**( U, M) : returns the matrix resulting from the product \f$ U M U^T \f$. If M is symmetric, the returned resulting matrix is also symmetric
*   **ROOT::Math::SimilarityT**( U, M) : returns the matrix resulting from the product \f$ U^T M U \f$. If M is symmetric, the returned resulting matrix is also symmetric

See \ref MatrixFunctions for the documentation of all existing matrix functions in the package.
The major Vector functions are:

*   **ROOT::Math::Dot**( v1, v2) : returns the scalar value resulting from the vector dot product
*   **ROOT::Math::Cross**( v1, v2) : returns the vector cross product for two vectors of size 3\. Note that the Cross product is not defined for other vector sizes
*   **ROOT::Math::Unit**( v) : returns unit vector. One can use also the _v.Unit()_ method.
*   **ROOT::Math::TensorProd**(v1,v2) : returns a general matrix M of size N1xN2 resulting from the [Tensor Product](http://en.wikipedia.org/wiki/Tensor_product) between the vector v1 of size N1) and v2 of size N2

See \ref VectFunction for the list and documentation of all of them.

### Matrix and Vector I/O

One can print (or write in an output stream) Vectors (and also Matrices) using the Print method or the << operator, like:

~~~ {.cpp}
v.Print(std::cout);
std::cout << v << std::endl;
~~~


In the ROOT distribution, the CINT dictionary is generated for SMatrix and SVector for double types and sizes up to 5\. This allows the storage of them in a ROOT file.

