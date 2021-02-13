\page SVectorDoc SVector Class Properties

The template ROOT::Math::SVector class has 2 template parameters which define, at compile
time, its properties. These are:

*   type of the contained elements, for example _float_ or double.
*   size of the vector.

### Creating a Vector

The following constructors are available to create a vector:

*   Default constructor for a zero vector (all elements equal to zero)
*   Constructor (and assignment) from a vector expression, like v = p*q + w. Due to the
    expression template technique, no temporary objects are created in this operation.
*   Construct a vector passing directly the elements. This is possible only for vector up to size 10.
*   Constructor from an iterator copying the data referred by the iterator. It is possible
    to specify the _begin_ and _end_ of the iterator or the _begin_ and the size. Note that
    for the Vector the iterator is not generic and must be of type _T*,_ where T is the type
    of the contained elements.

Here are some examples on how to create a vector. In the following we assume that we are
using the namespace ROOT::Math.

~~~ {.cpp}
SVector<double,N>  v;                         // create a vector of size N, v[i]=0
SVector<double,3>  v(1,2,3);                  // create a vector of size 3, v[0]=1,v[1]=2,v[2]=3
double   a[9] = {1,2,3,4,5,6,7,8,9};          // input data
SVector<double,9>  v(a,9);                    // create a vector using the a[] data
~~~


### Accessing and Setting Methods

The single vector elements can be set or retrieved using the _operator[i]_ , _operator(i)_
or the iterator interface. Notice that the index starts from zero and not from one as in
FORTRAN. Also no check is performed on the passed index. Furthermore, all the matrix
elements can be set also by using the ROOT::SVector::SetElements function passing a generic
iterator. The elements can be accessed also by using the ROOT::Math::SVector::apply(i) function.

~~~ {.cpp}
v[0]  = 1;                     // set the first element
v(1)  = 2;                     // set the second element
*(v.begin()+3) = 3;        // set the third element
// set vector elements from a std::vector<double>::iterator</double>
std::vector <double> w(3);
v.SetElements(w.begin(),w.end());

double x = m(i);                     // return the i-th element
x = m.apply(i);                  // return the i-th element
x = *(m.begin()+i);              // return the i-th element
~~~


In addition there are methods to place a sub-vector in a vector. If the size of the the
sub-vector is larger than the vector size a static assert ( a compilation error) is produced.

~~~ {.cpp}
SVector<double,N>  v;
SVector<double,M>  w;          // M <= N otherwise a compilation error is obtained later
// place a vector of size M starting from element ioff,  v[ioff + i] = w[i]
v.Place_at(w,ioff);
// return a sub-vector of size M starting from v[ioff]:  w[i] = v[ioff + i]
w = v.Sub < SVector>double,M> > (ioff);
~~~


For additional Vector functionality see the \ref MatVecFunctions page
