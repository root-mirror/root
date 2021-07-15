// @(#)root/mathcore:$Id$
// Author: L. Moneta Fri Sep 22 15:06:47 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class DistSampler

#ifndef ROOT_Math_DistSampler
#define ROOT_Math_DistSampler

#include "Math/IFunctionfwd.h"

#include "Math/WrappedFunction.h"

#include <vector>
#include <cassert>

class TRandom;

namespace ROOT {

   namespace Fit {

      class DataRange;
      class BinData;
      class UnBinData;
   }

   namespace Math {

      class DistSamplerOptions;

/**
   @defgroup Random Random Classes

   Pseudo-random numbers generator classes and for generation of random number distributions.
   These classes implement several pseudo-random number generators and method for generation of random numbers
   according to arbitrary distributions

   @ingroup MathCore

*/

/**
   Interface class for generic sampling of a distribution,
   i.e. generating random numbers according to arbitrary distributions

   @ingroup Random
*/


class DistSampler {

public:

   /// default constructor
   DistSampler() : fOwnFunc(false), fRange(0), fFunc(0) {}


   /// virtual destructor
   virtual ~DistSampler();



   /// set the parent function distribution to use for sampling (generic case)
   template<class Function>
   void SetFunction(Function & func, unsigned int dim) {
      WrappedMultiFunction<Function &> wf(func, dim);
      fData.resize(dim);
      // need to clone to avoid temporary
      DoSetFunction(wf,true);
   }

   /// set the parent function distribution to use for random sampling (one dim case)
   virtual void SetFunction(const ROOT::Math::IGenFunction & func)  {
      SetFunction<const ROOT::Math::IGenFunction>(func, 1);
   }


   /// set the parent function distribution to use for random sampling (multi-dim case)
   virtual void SetFunction(const ROOT::Math::IMultiGenFunction & func)  {
      DoSetFunction(func,false);
   }

   /// return the dimension of the parent distribution (and the data)
   unsigned int NDim() const { return fData.size(); }


   /**
      Initialize the sampling generator with the given algorithm.
      Implemented by the derived classes who needs it
      (like UnuranSampler).
      If nothing is specified use default algorithm
      from DistSamplerOptions::SetDefaultAlgorithm
   */
   virtual bool Init(const char * =""/* algorithm */) { return true;}

   /**
      Initialize the generators with the given DistSamplerOption object.
      The string will include the algorithm and in case additional options
      which can be interpreted by a re-implemented method in the derived class.
      The default implementation just calls the above method
      passing just the algorithm name
   */
   virtual bool Init(const DistSamplerOptions & opt );


   /**
       Set the random engine to be used.
       To be implemented by the derived classes who provides
       random sampling
   */
   virtual void SetRandom(TRandom *  ) {}

   /**
       Set the random seed for the TRandom instances used by the sampler
       classes.
       To be implemented by the derived classes who provides random sampling
   */
   virtual void SetSeed(unsigned int /*seed*/ ) {}

   /**
      Get the random engine used by the sampler.
      To be implemented by the derived classes who needs it
      Returns zero by default
    */
   virtual TRandom * GetRandom() { return 0; }

   /// Set the range in a given dimension.
   void SetRange(double xmin, double xmax, int icoord = 0);

   /// Set the range for all dimensions.
   void SetRange(const double * xmin, const double * xmax);
   /// Set the range for all dimensions (use std::vector)
   void SetRange(const std::vector<double> & xmin, const std::vector<double> & xmax){ 
      assert(xmin.size() >= NDim() && xmax.size() >= NDim());
      SetRange(xmin.data(),xmax.data());
   }

   /// Set the range using the ROOT::Fit::DataRange class.
   void SetRange(const ROOT::Fit::DataRange & range);

   /// Set the mode of the distribution (1D case).
   /// It could be useful or needed by some sampling methods.
   /// It is implemented by derived classes if needed (e.g. TUnuranSampler)
   virtual void SetMode(double  ) {}

   /// Set the mode of the distribution (Multi-dim case).
   virtual void SetMode(const std::vector<double> &) {}

   /// Set the normalization area of distribution.
   /// Implemented by derived classes if needed
   virtual void SetArea(double) {}

   /// Use the log of the provided pdf.
   /// Implemented by the derived classes 
   virtual void SetUseLogPdf(bool = true) {}

   /// Set usage of Derivative of PDF. 
   /// Can be implemented by derived class 
   virtual void SetDPdf(const ROOT::Math::IGenFunction & ) {}

   /// Set usage of Cumulative of PDF.
   /// Can be implemented by derived class
   virtual void SetCdf(const ROOT::Math::IGenFunction &) {}

   /// Get the parent distribution function (must be called after setting the function).
   const ROOT::Math::IMultiGenFunction & ParentPdf() const {
      return *fFunc;
   }

   /// Check if there is a parent distribution defined.
   bool HasParentPdf() const { return fFunc != nullptr; }

   /**
      Sample one event in one dimension.
      Specialized implementation could be provided by the derived classes
   */
   virtual double Sample1D() {
      Sample(&fData[0]);
      return fData[0];
   }

   /**
      Sample one event and return an array x with 
      sample coordinates values.
    */
   const double *  Sample() {
      Sample(&fData[0]);
      return &fData.front();
   }

   /**
      Sample one event in multi-dimension by filling the given array.
      Return false if the sampling failed.
      Abstract method to be re-implmented by the derived classes
   */
   virtual bool Sample(double * x) = 0;

   /**
      Sample one bin given an estimated of the pdf in the bin.
      (this can be function value at the center or its integral in the bin
      divided by the bin width)
      By default do not do random sample, just return the function values
      Typically Poisson statistics will be used
    */
   virtual bool SampleBin(double prob, double & value, double * error = 0) {
      value = prob;
      if (error) *error = 0;
      return true;
   }
   /**
      Sample a set of bins given a vector of probabilities
      Typically multinomial statistics will be used and the sum of the probabilities
      will be equal to the total number of events to be generated
      For sampling the bins indipendently, SampleBin should be used
    */
   virtual bool SampleBins(unsigned int n, const double * prob, double * values, double * errors  = 0)  {
      std::copy(prob,prob+n, values); // default impl returns prob values (Asimov data)
      if (errors) std::fill(errors,errors+n,0);
      return true;
   }


   /**
      Generate a un-binned data sets by fill the given data set.
      If dataset is not empty, append the new data.
   */
   virtual bool Generate(unsigned int nevt, ROOT::Fit::UnBinData & data);
   /**
      Generate a vector of events by fillling the passed data vector.
      The flag eventRow indicates how the events are arrenged in the multi-dim case. 
      The can be arranged in rows or in columns.  
      With eventRow=false events are the columns in data: {x1,x2,.....,xn},{y1,....yn} 
      With eventRow=true  events are rows in data: {x1,y1},{x2,y2},.....{xn,yn} 
   */
   virtual bool Generate(unsigned int nevt, double * data, bool eventRow = false);

   /**
      Generate a binned data set.
      A range must have been set before (otherwise inf is returned)
      and the bins are equidinstant in the previously defined range
      bin center values must be present in given data set
      If the sampler is implemented by a random one, the entries
      will be binned according to the Poisson distribution
      It is assumed the distribution is normalized, otherwise the nevt must be scaled
      accordingly. The expected value/bin nexp  = f(x_i) * binArea/ nevt
      Extend control if use a fixed (i.e. multinomial statistics) or floating total number of events
   */
   virtual bool Generate(unsigned int nevt, const  int * nbins, ROOT::Fit::BinData & data, bool extend = true, bool expErr = true);
   /**
      Same as before but passing the range in case of 1 dim data.
    */
   bool Generate(unsigned int nevt, int nbins, double xmin, double xmax, ROOT::Fit::BinData & data, bool extend = true, bool expErr = true ) {
      SetRange(xmin,xmax);
      int nbs[1]; nbs[0] = nbins;
      return Generate(nevt, nbs, data, extend, expErr);
   }


protected:

   // internal method to set the function
   virtual void DoSetFunction(const ROOT::Math::IMultiGenFunction & func, bool copy);
   // internal method to set the dimension
   virtual void DoSetDimension(unsigned int ndim);
   // check if generator have been initialized correctly and one can start generating
   bool IsInitialized() ;
   /// return the data range of the Pdf . Must be called after setting the function
   const ROOT::Fit::DataRange & PdfRange() const {
      assert(fRange);
      return *fRange;
   }

private:

   // private methods

   bool fOwnFunc;                         /// flag to indicate if the function is owned
   mutable std::vector<double> fData;     ///! internal array used to cached the sample data
   ROOT::Fit::DataRange    *   fRange;    /// data range
   const ROOT::Math::IMultiGenFunction * fFunc; /// internal function (ND)


};

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_DistSampler */
