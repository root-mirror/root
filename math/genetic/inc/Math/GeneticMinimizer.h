// @(#)root/mathcore:$Id

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class GeneticMinimizer

#ifndef ROOT_Math_GeneticMinimizer
#define ROOT_Math_GeneticMinimizer

#include "Math/Minimizer.h"

#include "RtypesCore.h"

#include <vector>
#include <string>

namespace TMVA {
   class IFitterTarget;
   class Interval;
}

namespace ROOT {
   namespace Math {


//_______________________________________________________________________________
/*
  structure containing the parameters of the genetic minimizer
 */
struct GeneticMinimizerParameters {

   Int_t fPopSize;
   Int_t fNsteps;
   Int_t fCycles;
   Int_t fSC_steps;
   Int_t fSC_rate;
   Double_t fSC_factor;
   Double_t fConvCrit;
   Int_t fSeed;


   // constructor with default value
   GeneticMinimizerParameters();
};



//_______________________________________________________________________________
/**
   GeneticMinimizer

   Minimizer class based on the Gentic algorithm implemented in TMVA

   @ingroup MultiMin
*/
class GeneticMinimizer: public ROOT::Math::Minimizer {

public:

   //GeneticMinimizer (int = 0);
   GeneticMinimizer (int i = 0);
   virtual ~GeneticMinimizer ();

   virtual void Clear();
   using ROOT::Math::Minimizer::SetFunction;
   virtual void SetFunction(const ROOT::Math::IMultiGenFunction & func);

   virtual bool SetLimitedVariable(unsigned int , const std::string& , double , double , double, double);
   virtual bool SetVariable(unsigned int ivar, const std::string & name, double val, double step);
   virtual bool SetFixedVariable(unsigned int ivar  , const std::string & name , double val);

   virtual  bool Minimize();
   virtual double MinValue() const;
   virtual double Edm() const;
   virtual const double *  X() const;
   virtual const double *  MinGradient() const;
   virtual unsigned int NCalls() const;

   virtual unsigned int NDim() const;
   virtual unsigned int NFree() const;

   virtual bool ProvidesError() const;
   virtual const double * Errors() const;

   virtual double CovMatrix(unsigned int i, unsigned int j) const;

   void SetParameters(const GeneticMinimizerParameters & params );

   void SetRandomSeed(int seed) { fParameters.fSeed = seed; }

   const GeneticMinimizerParameters & MinimizerParameters() const { return fParameters; }

   virtual ROOT::Math::MinimizerOptions Options() const;

   virtual void SetOptions(const ROOT::Math::MinimizerOptions & opt);

protected:

   void GetGeneticOptions(ROOT::Math::MinimizerOptions & opt) const;

   std::vector<TMVA::Interval*> fRanges;
   TMVA::IFitterTarget* fFitness;
   double fMinValue;
   std::vector<double> fResult;

   GeneticMinimizerParameters fParameters;

};


   } // end namespace Math
} // end namespace ROOT

#endif /* ROOT_Math_GeneticMinimizer */
