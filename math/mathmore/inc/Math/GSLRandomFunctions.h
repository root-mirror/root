// @(#)root/mathmore:$Id$
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
#ifndef ROOT_Math_GSLRandomFunctions
#define ROOT_Math_GSLRandomFunctions


//#include <type_traits>

#include "Math/RandomFunctions.h"

#include "Math/GSLRndmEngines.h"

#include <vector>

namespace ROOT {
namespace Math {


//___________________________________________________________________________________
   /**
       Specialized implementation of the Random functions based on the GSL library. 
       These will work onlmy with a GSLRandomEngine type 

       @ingroup  Random
   */


   template <class EngineType >
   class RandomFunctions<EngineType, ROOT::Math::GSLRandomEngine> : public RandomFunctions<EngineType, DefaultEngineType> {
      //class RandomFunctions<Engine, ROOT::Math::GSLRandomEngine>  {

      //typdef TRandomEngine DefaulEngineType;
      
   public:

      RandomFunctions() {} 

      RandomFunctions(EngineType & rng) : RandomFunctions<EngineType, DefaultEngineType>(rng) {}


      inline EngineType & Engine() { return  RandomFunctions<EngineType,DefaultEngineType>::Rng(); }
      
      double GausZig(double mean, double sigma) {
         return Engine().GaussianZig(sigma) + mean;
      }
      // double GausRatio(double mean, double sigma) {
      //    auto & r =  RandomFunctions<Engine,DefaultEngineType>::Rng();
      //    return r.GaussianRatio(sigma) + mean; 
      // }
      
      /**
         Gaussian distribution. Default method (use Ziggurat)
      */
      double Gaus(double mean = 0, double sigma = 1) {
         return mean + Engine().GaussianZig(sigma);
      }

      /**
         Gaussian distribution (Box-Muller method)
      */
      double GausBM(double mean = 0, double sigma = 1) {
         return mean + Engine().Gaussian(sigma);
      }
      
      /**
         Gaussian distribution (Ratio Method)
      */
      double GausR(double mean = 0, double sigma = 1) {
         return mean + Engine().GaussianRatio(sigma);
      }
      
      /**
         Gaussian Tail distribution
      */
      double GaussianTail(double a, double sigma = 1) {
         return Engine().GaussianTail(a,sigma);
      }
      
      /**
         Bivariate Gaussian distribution with correlation
      */
      void Gaussian2D(double sigmaX, double sigmaY, double rho, double &x, double &y) {
         Engine().Gaussian2D(sigmaX, sigmaY, rho, x, y);
      }

      /**
         Exponential distribution
      */
      double Exp(double tau) {
         return Engine().Exponential(tau);
      }
      /**
         Breit Wigner distribution
      */
      double BreitWigner(double mean = 0., double gamma = 1) {
         return mean + Engine().Cauchy( gamma/2.0 );
      }
      
      /**
         Landau distribution
      */
      double Landau(double mean = 0, double sigma = 1) {
         return mean + sigma*Engine().Landau();
      }

      /**
         Gamma distribution
      */
      double Gamma(double a, double b) {
         return Engine().Gamma(a,b);
      }

      /**
         Beta distribution
      */
      double Beta(double a, double b) {
         return Engine().Beta(a,b);
      }

      /**
         Log Normal distribution
      */
      double LogNormal(double zeta, double sigma) {
         return Engine().LogNormal(zeta,sigma);
      }

      /**
         Chi square distribution
      */
      double ChiSquare(double nu) {
         return Engine().ChiSquare(nu);
      }

      /**
         F distrbution
      */
      double FDist(double nu1, double nu2) {
         return Engine().FDist(nu1,nu2);
      }

      /**
         t student distribution
      */
      double tDist(double nu) {
         return Engine().tDist(nu);
      }
      /**
         Rayleigh distribution
      */
      double Rayleigh(double sigma)  {
         return Engine().Rayleigh(sigma);
      }

      /**
         Logistic distribution
      */
      double Logistic(double a) {
         return Engine().Logistic(a);
      }

      /**
         Pareto distribution
      */
      double Pareto(double a, double b)  {
         return Engine().Pareto(a,b);
      }

      /**
         generate random numbers in a 2D circle of radious 1
      */
      void Circle(double &x, double &y, double r = 1) {
         Engine().Dir2D(x,y);
         x *= r;
         y *= r;
      }

      /**
         generate random numbers in a 3D sphere of radious 1
      */
      void Sphere(double &x, double &y, double &z,double r = 1) {
         Engine().Dir3D(x,y,z);
         x *= r;
         y *= r;
         z *= r;
      }

      /**
         Poisson distribution
      */
      unsigned int Poisson(double mu) {
         return Engine().Poisson(mu);
      }

      /**
         Binomial distribution
      */
      unsigned int Binomial(unsigned int ntot, double prob) {
         return Engine().Binomial(prob,ntot);
      }

      /**
         Negative Binomial distribution
         First parameter is n, second is probability
         To be consistent with Random::Binomial
      */
      unsigned int NegativeBinomial(double n, double prob) {
         return Engine().NegativeBinomial(prob,n);
      }

      /**
         Multinomial distribution
      */
      std::vector<unsigned int> Multinomial( unsigned int ntot, const std::vector<double> & p ) {
         return Engine().Multinomial(ntot,p);
      }



   };




} // namespace Math
} // namespace ROOT

#endif /* ROOT_Math_GSLRandomFunctions */
