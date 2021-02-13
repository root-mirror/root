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
#ifndef ROOT_Math_Random
#define ROOT_Math_Random

/**
@defgroup Random Interface class for Random number generation
*/

#include "Math/RandomFunctions.h"

#include <string>
#include <vector>


namespace ROOT {
namespace Math {


//___________________________________________________________________________________
   /**
       Documentation for the Random class

       @ingroup  Random
   */

   template < class Engine>
   class Random {

   public:

      typedef typename Engine::BaseType EngineBaseType;
      typedef RandomFunctions<Engine, EngineBaseType> RndmFunctions;

      Random() :
         fEngine(),
         fFunctions(fEngine)
      {}

      explicit Random(unsigned int seed) :
         fEngine(),
         fFunctions(fEngine)
      {
         fEngine.SetSeed(seed);
      }

      double Rndm() {
         return fEngine();
      }

      /**
         Generate an array of random numbers between ]0,1]
         0 is excluded and 1 is included
         Function to preserve ROOT Trandom compatibility
      */
      void RndmArray(int n, double * array) {
         fEngine.RandomArray(array, array+n);
      }

      /**
         Return the type (name) of the used generator
      */
      std::string Type() const {
         return fEngine.Name();
      }

      /**
       Return the size of the generator state
     */
      unsigned int EngineSize() const {
         return fEngine.Size();
      }


      double operator() (){
         return fEngine();
      }

      uint64_t Integer() {
         return fEngine.IntRndm();
      }

      static uint64_t MaxInt()  {
         return Engine::Max();
      }

      Engine & Rng() {
         return fEngine;
      }

      /// Exponential distribution
      double Exp(double tau) {
         return fFunctions.Exp(tau);
      }

      double Gaus(double mean = 0, double sigma = 1) {
         return fFunctions.Gaus(mean,sigma);
      }

      /// Gamma distribution
      double Gamma(double a, double b) {
         return fFunctions.Gamma(a,b);
      }

      /// Beta distribution
      double Beta(double a, double b) {
         return fFunctions.Beta(a,b);
      }

      ///Log-normal distribution
      double LogNormal(double zeta, double sigma) {
         return fFunctions.LogNormal(zeta,sigma);
      }

      /// chi-square
      double  ChiSquare(double nu) {
         return fFunctions.ChiSquare(nu);
      }

      /// Rayleigh distribution
      double  Rayleigh(double sigma) {
         return fFunctions.Rayleigh(sigma);
      }

      /// Logistic distribution
      double  Logistic(double a) {
         return fFunctions.Logistic(a);
      }

      /// Pareto distribution
      double  Pareto(double a, double b) {
         return fFunctions.Pareto(a, b);
      }

      ///F-distribution
      double FDist(double nu1, double nu2) {
         return fFunctions.FDist(nu1,nu2);
      }

      ///  t student distribution
      double tDist(double nu) {
         return fFunctions.tDist(nu);
      }

      /// Landau distribution
      double Landau(double m = 0, double s = 1) {
         return fFunctions.Landau(m,s);
      }
     ///  Breit Wigner distribution
      double BreitWigner(double mean = 0., double gamma = 1) {
         return fFunctions.BreitWigner(mean,gamma);
      }

     ///  generate random numbers in a 2D circle of radious 1
      void Circle(double &x, double &y, double r = 1) {
         fFunctions.Circle(x,y,r);
      }

      ///  generate random numbers in a 3D sphere of radious 1
      void Sphere(double &x, double &y, double &z,double r = 1) {
         fFunctions.Sphere(x,y,z,r);
      }


      ///discrete distributions

      /// Binomial distribution
      unsigned int Binomial(unsigned int ntot, double prob) {
         return fFunctions.Binomial(prob,ntot);
      }


      ///   Poisson distribution
      unsigned int Poisson(double mu)  {
         return fFunctions.Poisson(mu);
      }

      /// Negative Binomial distribution
      ///  First parameter is n, second is probability
      ///  To be consistent with Random::Binomial
     unsigned int NegativeBinomial(double n, double prob) {
      return fFunctions.NegativeBinomial(prob,n);
     }

      ///     Multinomial distribution
      std::vector<unsigned int> Multinomial( unsigned int ntot, const std::vector<double> & p ) {
         return fFunctions.Multinomial(ntot,p);
      }



      double Uniform(double a, double b) {
         return fFunctions.Uniform(a,b);
      }
      double Uniform(double a = 1.0) {
         return fFunctions.Uniform(a);
      }
      double Uniform2(double a, double b) {
         return fFunctions.UniformBase(a,b);
      }


      RandomFunctions<Engine,EngineBaseType> & Functions() {
         return fFunctions;
      }

      void SetSeed(int seed) { fEngine.SetSeed(seed);}

   private:

      Engine fEngine;             //  random generator engine
      RndmFunctions fFunctions;   //! random functions object


  };




} // namespace Math
} // namespace ROOT

#include "Math/MixMaxEngine.h"
#include "Math/MersenneTwisterEngine.h"
#include "Math/StdEngine.h"

namespace ROOT {
namespace Math {

   /// Useful typedef definitions

   typedef   Random<ROOT::Math::MixMaxEngine<240,0>>            RandomMixMax;
   typedef   Random<ROOT::Math::MersenneTwisterEngine>   RandomMT19937;
   typedef   Random<ROOT::Math::StdEngine<std::mt19937_64>> RandomMT64;
   typedef   Random<ROOT::Math::StdEngine<std::ranlux48>> RandomRanlux48;

} // namespace Math
} // namespace ROOT


#endif /* ROOT_Math_Random */
