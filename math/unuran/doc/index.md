\defgroup Unuran Unuran
\ingroup Math
\brief Universal Non Uniform Random number generator for generating non uniform pseudo-random numbers.


UNU.RAN, (Universal Non Uniform Random number generator for generating non uniform pseudo-random numbers)
is an ANSI C library licensed under GPL.<br>
It contains universal (also called automatic or black-box) algorithms that can generate random numbers from
large classes of continuous or discrete distributions, and also from practically all standard distributions.
An extensive online documentation are available at the (UNU.RAN Web Site)[http://statistik.wu-wien.ac.at/unuran/]

New classes have been introduced to use the UNU.RAN C library from ROOT and C++ from ROOT and using C++ objects.
To use UNU.RAN one needs always an instance of the class **TUnuran**.
It can then be used in two distinct ways:
- using the UNU.RAN native string API for pre-defined distributions (see<a href="http://statistik.wu-wien.ac.at/unuran/doc/unuran.html#StringAPI"> UNU.RAN documentation</a> for the string API):

~~~{.cpp}
          TUnuran unr;
          //initialize unuran to generate normal random numbers using a "arou" method
          unr.Init("normal()","method=arou");
          //......
          // sample distributions N times (generate N random numbers)
          for (int i = 0; i &lt; N; ++i)
                double x = unr.Sample();

~~~


- Using a distribution object. We have then the following cases depending on the dimension and the distribution object.

- For 1D distribution the class **TUnuranContDist** must be used.
    - A **TUnuranContDist** object can be created from a function
    providing the pdf (probability density function) and optionally one providing the derivative of the pdf.
   - If the derivative is not provided and the generation method requires it, then it is estimated numerically.
    - The user can optionally provide the
     - cdf (cumulative distribution function) via the **TUnuranContDist::SetCdf** function,
     - the mode via **TUnuranContDist::SetMode**,
     - the domain via **TUnuranContDist::SetDomain** for generating numbers in a restricted region,
     - the area below the pdf via **TUnuranContDist::SetPdfArea**.

Some of this information is required depending on the chosen UNURAN generation method.

~~~~{.cpp}
       //1D case: create a distribution from two TF1 object pointers pdfFunc
       TUnuranContDist  dist( pdfFunc);
       //initialize unuran passing the distribution and a string defining the method
       unr.Init(dist, "method=hinv");
       // sample distribution  N times (generate N random numbers)
       for (int i = 0; i &lt; N; ++i)
          double x = unr.Sample();
~~~~

- For multi-dimensional distribution the class **TUnuranMultiContDist** must be used.
In this case only the multi-dimensional pdf is required

~~~~{.cpp}
      //Multi-Dim case from a TF1 (or TF2 or TF3) object describing a multi-dimensional function
      TUnuranMultiContDist  dist( pdfFuncMulti);
      // the recommended method for multi-dimensional function is "hitro"
      unr.Init(dist, "method=hitro");
      // sample distribution  N times (generate N random numbers)
      double x[NDIM];
      for (int i = 0; i &lt; N; ++i)
      unr.SampleMulti(x);
~~~~

- For discrete distribution the class **TUnuranDiscrDist** must be used.
   The distribution can be initialized from a TF1 or from a vector of probabilities.

~~~~{.cpp}
      // create distribution from a vector of probabilities
         double pv[NSize] = {0.1,0.2,.......};
      TUnuranDiscrDist  dist(pv, pv+NSize);
      // the recommended method for discrete distribution is
      unr.Init(dist, "method=dgt");
      // sample N times (generate N random numbers)
      for (int i = 0; i &lt; N; ++i)
         int k = unr.SampleDiscr();
~~~~

- For empirical distribution the class **TUnuranEmpDist** must be used.
  In this case one can generate random numbers from a set of data (un-binned)  in one or multi-dimension or
  from a set of binned data in one dimension (similar to TH1::GetRandom() ).
  - For unbin data the parent distribution is estimated by UNU.RAN using a gaussian kernel smoothing algorithm.
   One can create the distribution class directly from a vector of data or from the buffer of TH1.

~~~~{.cpp}
      // create distribution from a set of data 1D
      // vdata is an std::vector containing the data
      TUnuranEmpDist  dist( vdata.begin(),vdata.end());
      unr.Init(dist);
         // sample N times (generate N random numbers)
      for (int i = 0; i &lt; N; ++i)
         double x = unr.Sample();
~~~~

- In the case of multi-dimension empirical distributions one needs to pass in addition to the iterators, the data dimension. It is assumed that the data are stored in the vector in this order : `(x0,y0,...),(x1,y1,....)`.

- For binned data (only one dimensional data are supported) one uses directly the histogram

~~~{.cpp}
      // create an empirical distribution from an histogram
      // if the histogram has a buffer one must use TUnuranEmpDist(h1,false)
      TH1 * h1 = ... // histogram pointer
      TUnuranEmpDist  binDist( h1);
      unr.Init(binDist);
         // sample N times (generate N random numbers)
      for (int i = 0; i &lt; N; ++i)
         double x = unr.Sample();
~~~

- This is equivalent to TH1::GetRandom(), but sampling is faster, therefore, since it requires some initialization time,
  it becomes convenient when generating a large sample of random numbers.

Functionality is also provided via the C++ classes for using a different random number generator by passing a
TRandom pointer when constructing the TUnuran class (by default the ROOT gRandom is passed to UNURAN).

The (UNU.RAN documentation)[http://statistik.wu-wien.ac.at/unuran/doc/unuran.html#Top] provides a detailed
description of all the available methods and the possible options which one can pass to UNU.RAN for the various distributions.
