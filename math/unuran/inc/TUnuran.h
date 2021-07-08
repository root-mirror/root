// @(#)root/unuran:$Id$
// Author: L. Moneta Tue Sep 26 16:25:09 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class TUnuran

#ifndef ROOT_TUnuran
#define ROOT_TUnuran

#include <string>

#include "TUnuranBaseDist.h"


class TUnuranContDist;
class TUnuranDiscrDist;
class TUnuranMultiContDist;
class TUnuranEmpDist;

#include <memory>


/**

   \class TUnuran
   \ingroup Unuran

   TUnuran class.
   Interface to the UNU.RAN package for generating non uniform random
   numbers. This class wraps the UNU.RAN calls in C++ methods.
   It provides methods for initializing Unuran and then to sample the
   desired distribution.
   It provides support for initializing UNU.RAN in these following way (various signatures
   for TUnuran::Init)
   - with string API via TUnuran::Init passing the distribution type and the method
   - using a one-dimensional distribution object defined by TUnuranContDist
   - using a multi-dimensional distribution object defined by TUnuranMultiContDist
   - using a discrete one-dimensional distribution object defined by TUnuranDiscrDist
   - using an empirical distribution defined by TUnuranEmpDist
   - using pre-defined distributions. Presently only support for Poisson (TUnuran::InitPoisson)
     and Binomial (TUnuran::InitBinomial) are provided. Other distributions can however be generated
     using the previous methods (in particular via the string API)

   The sampling is provided via these methods:
    - TUnuran::Sample()   returns a double for all one-dimensional distribution
    - TUnuran::SampleDiscr()  returns an integer for one-dimensional discrete distribution
    - TUnuran::Sample(double *) sample a multi-dimensional distribution. A pointer to a vector with
      size at least equal to the distribution dimension must be passed

   In addition is possible to set the random number generator in the constructor of the class, its seed
   via the TUnuran::SetSeed() method.
*/



//class TUnuranGenerator;
struct unur_gen;
typedef struct unur_gen UNUR_GEN;

// struct unur_urng_generic;
// typedef struct unur_urng_generic UNUR_URNG;

struct unur_distr;
typedef struct unur_distr UNUR_DISTR;

struct unur_urng;
typedef struct unur_urng  UNUR_URNG;


class TRandom;
class TH1;

class TUnuran {

public:

   /**
      Constructor with a generator instance and given level of log output
   */
   TUnuran (TRandom * r = 0, unsigned int log = 0);


   /**
      Destructor
   */
   ~TUnuran ();

private:
   // usually copying is non trivial, so we make this unaccessible

   /**
      Copy constructor
   */
   TUnuran(const TUnuran &);

   /**
      Assignment operator
   */
   TUnuran & operator = (const TUnuran & rhs);

public:
   /**
      Initialize with Unuran string API interface.
      See  https://statmath.wu.ac.at/unuran/doc/unuran.html#StringAPI

      @param distr   : UNU.RAN distribution string
      @param  method : UNU.RAN  method string

      Here is an example using the string API:
      ```
      Tunuran unr;
      unr.Init("normal(3.,0.75); domain = (0,inf)", "method = tdr; c = 0");
      ```
   */
   bool Init(const std::string & distr, const std::string  & method);


   /**
      Initialize method for continuous one-dimensional distribution.
      User must provide a distribution object (which is copied inside) and a string for a method.
      For the list of available method for 1D cont. distribution see the
      <A href="http://statmath.wu-wien.ac.at/unuran/doc/unuran.html#Methods_005ffor_005fCONT">UnuRan doc</A>.
      A re-initialization is needed whenever distribution parameters have been changed.
      Note that the method string can contain in addition to the method name all the specific method 
      parameters specified using the UNURAN methiod string API. 
      For example a valid string can be `"method=arou; max_segments=1000; max_sqhratio = 0.9"`
   */
   bool Init(const TUnuranContDist & distr, const std::string & method = "auto");

   /**
      Initialize method for continuous multi-dimensional distribution.
      User must provide a distribution object (which is copied inside) and a string for a method.
      For the list of available method for multivariate cont. distribution see the
      <A href="http://statmath.wu-wien.ac.at/unuran/doc/unuran.html#Methods_005ffor_005fCVEC">UnuRan doc</A>
      A re-initialization is needed whenever distribution parameters have been changed.

   */
   bool Init(const TUnuranMultiContDist & distr, const std::string & method = "hitro");


   /**
      Initialize method for continuous one-dimensional discrete distribution.
      User must provide a distribution object (which is copied inside) and a string for a method.
      For the list of available method for 1D discrete distribution see the
      <A href="http://statmath.wu-wien.ac.at/unuran/doc/unuran.html#Methods_005ffor_005fDISCR">UnuRan doc</A>
      A re-initialization is needed whenever distribution parameters have been changed.

   */
   bool Init(const TUnuranDiscrDist & distr, const std::string & method = "auto");


   /**
      Initialize method for continuous empirical distribution.
      User must provide a distribution object (which is copied inside) and a string for a method.
      The distribution object can represent binned (only 1D) or unbinned (1D or multi-dim) data
      The method for the unbinned empirical distribution are based on the kernel smoothing, see
      <A href="http://statmath.wu-wien.ac.at/software/unuran/doc/unuran.html#EMPK">UnuRan doc</A>
      A re-initialization is needed whenever distribution parameters have been changed.

   */
   bool Init(const TUnuranEmpDist & distr, const std::string & method = "empk");


   /**
      Initialize method for the Poisson distribution.
      Used to generate poisson numbers for a constant parameter mu of the Poisson distribution.
      Use after the method TUnuran::SampleDiscr to generate the numbers.
      The flag reinit perform a fast re-initialization when only the distribution parameters
      are changed in the subsequent calls.
      If the same TUnuran object is used to generate with other distributions it cannot be used.
   */
   bool InitPoisson(double mu, const std::string & method = "dstd");

   /**
      Initialize method for the Binomial distribution.
      Used to generate poisson numbers for a constant parameters (n,p) of the Binomial distribution.
      Use after the method TUnuran::SampleDiscr to generate the numbers.
      The flag reinit perform a fast re-initialization when only the distribution parameters
      are changed in the subsequent calls.
      If the same TUnuran object is used to generate with other distributions it cannot be used.
   */
   bool InitBinomial(unsigned int ntot, double prob, const std::string & method = "dstd");

   /**
      Reinitialize UNURAN by changing the distribution parameters but mantaining same distribution and method.
      It is implemented now only for predefined discrete distributions like the poisson or the binomial
   */
   bool ReInitDiscrDist(unsigned int npar, double * params);

   /**
      Sample 1D distribution.
      User is responsible for having previously correctly initialized with TUnuran::Init
   */
   double Sample();

   /**
      Sample multidimensional distributions.
      User is responsible for having previously correctly initialized with TUnuran::Init
   */
   bool SampleMulti(double * x);

   /**
      Sample discrete distributions.
      User is responsible for having previously correctly initialized with TUnuran::Init
   */
   int SampleDiscr();

   /**
      Set the random engine.
      Must be called before init to have effect
    */
   void SetRandom(TRandom * r) {
      fRng = r;
   }

   /**
      Return instance of the random engine used.
    */
   TRandom * GetRandom() {
      return fRng;
   }

   /**
      Return an information string about used unuran generator method.
      @param extended : if true return some help about the possible options for the the method.
   */
   std::string GetInfo(bool extended = false);

   /**
      Return an ID string about the unuran generator method.
   */
   std::string GetGenId() const;

   /**
      Return the dimension of unuran generator method.
      For 1D method returns 1 and for the multi-dimensional case 
      must be equal to the distribution dimension.
   */
   int GetDimension() const;

   /**
      Return the type of the distribution.
      See documentation of `unuran_distr_get_type` for the possible
      types of distributions.
   */
   int GetDistType() const; 

   /// Return true for a univariate continous distribution.
   bool IsDistCont() const; 
   /// Return true for a multivariate continous distribution.
   bool IsDistMultiCont() const;  
    /// Return true for a discrete distribution.
   bool IsDistDiscrete() const; 
    /// Return true for an empirical distribution.
   bool IsDistEmpirical() const; 
   


   /**
      set the seed for the random number generator
    */
   void SetSeed(unsigned int seed);

   /**
      set log level
   */
   bool SetLogLevel(unsigned int iflag = 1);

   /**
      set stream for log and error (not yet implemented)
   */
   bool SetLogStream() { return false;}

   /**
      used Unuran method
    */
   const std::string & MethodName() const { return fMethod; }

protected:


   bool SetRandomGenerator();

   bool SetContDistribution(const TUnuranContDist & dist );

   bool SetMultiDistribution(const TUnuranMultiContDist & dist );

   bool SetDiscreteDistribution(const TUnuranDiscrDist & dist );

   bool SetEmpiricalDistribution(const TUnuranEmpDist & dist );

   /**
      change the method and initialize Unuran with the previously given distribution
    */
   bool SetMethodAndInit();



// private:

   UNUR_GEN * fGen;                      //pointer to the UnuRan C generator struct
   UNUR_DISTR * fUdistr;                 //pointer to the UnuRan C distribution struct
   UNUR_URNG  * fUrng;                   // pointer to Unuran C random generator struct
   std::unique_ptr<TUnuranBaseDist> fDist; // pointer for distribution wrapper
   TRandom * fRng;                       //pointer to ROOT random number generator
   std::string fMethod;                  //string representing the method

};


#endif /* ROOT_Math_TUnuran */
