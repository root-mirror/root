// @(#)root/mathcore:$Id$
// Authors: L. Moneta    8/2015

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015 , ROOT MathLib Team                             *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for random class
//
//
// Created by: Lorenzo Moneta  : Tue 4 Aug 2015
//
//
#ifndef ROOT_Math_RandomFunctions
#define ROOT_Math_RandomFunctions


#include <type_traits>
#include <cmath>
#include "RtypesCore.h"
#include "TMath.h"
#include <cassert>
#include <vector>

#include "TRandomEngine.h"


namespace ROOT {
namespace Math {


//___________________________________________________________________________________


   // class DefaultEngineType {};

   
   /**
       Documentation for the RandomFunction class 

       @ingroup  Random
   */


   typedef TRandomEngine DefaultEngineType;
   //class DefaultEngineType {};  // for generic types



      /**
      Definition of the generic impelmentation class for the RandomFunctions.
      Needs to have specialized implementations on the different type of engines      
    */
   template <class EngineBaseType> 
   class  RandomFunctionsImpl {
   public:
      void SetEngine(void *) {}
   };

   /**
      Implementation class for the RandomFunction for all the engined that derives from 
      TRandomEngine class, which defines an interface which has TRandomEngine::Rndm()
      In this way we can have a common implementation for the RandomFunctions
    */

   template<>
   class RandomFunctionsImpl<TRandomEngine> { 

   public: 

      /// class constructor 
      RandomFunctionsImpl() : fBaseEngine(0) {}

      void SetEngine(void *r) {
         fBaseEngine = static_cast<TRandomEngine*>(r);
         assert(fBaseEngine);  // to be sure the static cast works 
      }
      

      ///Generate binomial numbers
      int Binomial(int ntot, double prob);

      /// Return a number distributed following a BreitWigner function with mean and gamma.
      double BreitWigner(double mean, double gamma);

      /// Generates random vectors, uniformly distributed over a circle of given radius.
      ///   Input : r = circle radius
      ///   Output: x,y a random 2-d vector of length r
      void Circle(double &x, double &y, double r);

      /// Returns an exponential deviate.
      ///    exp( -t/tau )
      double  Exp(double tau);
      
      /// generate Gaussian number using Box-Muller method
      double GausBM( double mean, double sigma);

      /// generate random numbers according to the Accemptance-Complemet-Ratio method
      double GausACR( double mean, double sigma);

      /// Generate a random number following a Landau distribution
      /// with location parameter mu and scale parameter sigma:
      ///      Landau( (x-mu)/sigma )
      double Landau(double mu, double sigma);

      /// Generates a random integer N according to a Poisson law.
      /// Prob(N) = exp(-mean)*mean^N/Factorial(N)
      int Poisson(double mean);
      double PoissonD(double mean);

      /// Generate numbers distributed following a gaussian with mean=0 and sigma=1.
      /// Using the Box-Muller method 
      void Rannor(double &a, double  &b);

      /// Generates random vectors, uniformly distributed over the surface
      /// of a sphere of given radius.
      void Sphere(double &x, double &y, double &z, double r);

      /// generate random numbers following a Uniform distribution in the [a,b] interval
      double Uniform(double a, double b);
      double Uniform(double a);

   protected:
      TRandomEngine* fBaseEngine;

   private:
      // Internal method used by the functions 
      double Rndm() { return fBaseEngine->Rndm(); }
      // for internal usage
      double Gaus(double mean, double sigma) { return GausACR(mean,sigma); }


   };


   template < class Engine, class EngineBaseType>
   class RandomFunctions { //: public RandomFunctionsImpl<EngineBaseType> {


   public:

      //RandomFunctions() {} 

      RandomFunctions(Engine & rng) : fEngine(&rng) {
         fImpl.SetEngine(&rng);
      }

      /// destructor (no op) we do not mantain the engine)
      ~RandomFunctions() {}


      /// non-virtual method
      inline double operator() () { return (*fEngine)(); }


      ///Generate binomial numbers
      int Binomial(int ntot, double prob) {
         return fImpl.Binomial(ntot,prob); 
      }

      /// Return a number distributed following a BreitWigner function with mean and gamma.
      double BreitWigner(double mean, double gamma) {
         return fImpl.BreitWigner(mean,gamma);
      }

      /// Generates random vectors, uniformly distributed over a circle of given radius.
      ///   Input : r = circle radius
      ///   Output: x,y a random 2-d vector of length r
      void Circle(double &x, double &y, double r) {
         return fImpl.Circle(x,y,r);
      }

      /// Returns an exponential deviate.
      ///    exp( -t/tau )
      double  Exp(double tau) {
         return fImpl.Exp(tau); 
      }
      
      /// generate Gaussian number using Box-Muller method
      double GausBM( double mean, double sigma) {
         return fImpl.GausBM(mean,sigma);
      }

      /// generate random numbers according to the Accemptance-Complemet-Ratio method
      double GausACR( double mean, double sigma) {
         return fImpl.GausACR(mean, sigma); 
      }

      /// Generate a random number following a Landau distribution
      /// with location parameter mu and scale parameter sigma:
      ///      Landau( (x-mu)/sigma )
      double Landau(double mu, double sigma) {
         return fImpl.Landau(mu,sigma); 
      }

      /// Generates a random integer N according to a Poisson law.
      /// Prob(N) = exp(-mean)*mean^N/Factorial(N)
      int Poisson(double mean) { return fImpl.Poisson(mean); }
      double PoissonD(double mean) { return fImpl.PoissonD(mean); }

      /// Generate numbers distributed following a gaussian with mean=0 and sigma=1.
      /// Using the Box-Muller method 
      void Rannor(double &a, double  &b) {
         return fImpl.Rannor(a,b);
      }

      /// Generates random vectors, uniformly distributed over the surface
      /// of a sphere of given radius.
      void Sphere(double &x, double &y, double &z, double r) {
         return fImpl.Sphere(x,y,z,r);
      }

      /// generate random numbers following a Uniform distribution in the [a,b] interval
      double Uniform(double a, double b) {
         return (b-a) * Rndm_impl() + a; 
      }
     
      /// generate random numbers following a Uniform distribution in the [0,a] interval
      double Uniform(double a) {
         return a * Rndm_impl() ; 
      }


      /// generate Gaussian number using defqault method
      inline double Gaus( double mean, double sigma) {
         return fImpl.GausACR(mean,sigma);
      }


      // /// re-implement Gaussian 
      // double GausBM2(double mean, double sigma) {
      //    double y =  Rndm_impl();
      //    double z =  Rndm_impl();
      //    double x = z * 6.28318530717958623;
      //    double radius = std::sqrt(-2*std::log(y));
      //    double g = radius * std::sin(x);
      //    return mean + g * sigma; 
      // }


      /// methods which are only for GSL random generators 
      

      /// Gamma functions (not implemented here, requires a GSL random engine)
      double Gamma( double , double ) {
         //r.Error("Error: Gamma() requires a GSL Engine type"); 
         static_assert(std::is_fundamental<Engine>::value,"Error: Gamma() requires a GSL Engine type");
         return 0;
      }
      double Beta( double , double ) {
         static_assert(std::is_fundamental<Engine>::value,"Error: Beta() requires a GSL Engine type");
         return 0;
      }
      double LogNormal(double, double) {
         static_assert(std::is_fundamental<Engine>::value,"Error: LogNormal() requires a GSL Engine type");
         return 0;
      }
      double ChiSquare(double) {
         static_assert(std::is_fundamental<Engine>::value,"Error: ChiSquare() requires a GSL Engine type");
         return 0;
      }
      double Rayleigh( double ) {
         static_assert(std::is_fundamental<Engine>::value,"Error: Rayleigh() requires a GSL Engine type");
         return 0;
      }
      double Logistic( double ) {
         static_assert(std::is_fundamental<Engine>::value,"Error: Logistic() requires a GSL Engine type");
         return 0;
      }
      double Pareto( double , double ) {
         static_assert(std::is_fundamental<Engine>::value,"Error: Pareto() requires a GSL Engine type");
         return 0;
      }
      double FDist(double, double) {
         static_assert(std::is_fundamental<Engine>::value,"Error: FDist() requires a GSL Engine type");
         return 0;
      }
      double tDist(double) {
         static_assert(std::is_fundamental<Engine>::value,"Error: tDist() requires a GSL Engine type");
         return 0;
      }
      unsigned int NegativeBinomial(double , double ) {
         static_assert(std::is_fundamental<Engine>::value,"Error: NegativeBinomial() requires a GSL Engine type");
         return 0;
      }
      std::vector<unsigned int> MultiNomial(unsigned int, const std::vector<double> &){
         static_assert(std::is_fundamental<Engine>::value,"Error: MultiNomial() requires a GSL Engine type");
         return std::vector<unsigned int>();
      }


   protected:

      Engine & Rng() { assert(fEngine); return *fEngine; }

      /// Internal impelmentation to return random number
      /// Since this one is not a virtual function is faster than Rndm
      inline double Rndm_impl() { return (*fEngine)(); }


   private:    

      Engine * fEngine;   //! random number generator engine
      RandomFunctionsImpl<EngineBaseType> fImpl;   //! instance of the class implementing the functions
      

  };




} // namespace Math
} // namespace ROOT

#endif /* ROOT_Math_RandomFunctions */
