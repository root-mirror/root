// @(#)root/minuit:$Id$
// Author: L. Moneta Wed Oct 25 16:28:55 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class TLinearMinimizer

#ifndef ROOT_TLinearMinimizer
#define ROOT_TLinearMinimizer

#include "Math/Minimizer.h"

#include "Rtypes.h"

#include <vector>
#include <string>

class TLinearFitter;




/**
   TLinearMinimizer class: minimizer implementation based on TMinuit.
*/
class TLinearMinimizer  : public ROOT::Math::Minimizer {

public:

   /**
      Default constructor
   */
   TLinearMinimizer (int type = 0);

   /**
      Constructor from a char * (used by PM)
   */
   TLinearMinimizer ( const char * type );

   /**
      Destructor (no operations)
   */
   virtual ~TLinearMinimizer ();

private:
   // usually copying is non trivial, so we make this unaccessible

   /**
      Copy constructor
   */
   TLinearMinimizer(const TLinearMinimizer &);

   /**
      Assignment operator
   */
   TLinearMinimizer & operator = (const TLinearMinimizer & rhs);

public:

   /// set the fit model function
   virtual void SetFunction(const ROOT::Math::IMultiGenFunction & func);

   /// set the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGradFunction & func);

   /// set free variable (dummy impl. )
   virtual bool SetVariable(unsigned int , const std::string & , double , double ) { return false; }

   /// set fixed variable (override if minimizer supports them )
   virtual bool SetFixedVariable(unsigned int /* ivar */, const std::string & /* name */, double /* val */);

   /// method to perform the minimization
   virtual  bool Minimize();

   /// return minimum function value
   virtual double MinValue() const { return fMinVal; }

   /// return expected distance reached from the minimum
   virtual double Edm() const { return 0; }

   /// return  pointer to X values at the minimum
   virtual const double *  X() const { return &fParams.front(); }

   /// return pointer to gradient values at the minimum
   virtual const double *  MinGradient() const { return 0; } // not available in Minuit2

   /// number of function calls to reach the minimum
   virtual unsigned int NCalls() const { return 0; }

   /// this is <= Function().NDim() which is the total
   /// number of variables (free+ constrained ones)
   virtual unsigned int NDim() const { return fDim; }

   /// number of free variables (real dimension of the problem)
   /// this is <= Function().NDim() which is the total
   virtual unsigned int NFree() const { return fNFree; }

   /// minimizer provides error and error matrix
   virtual bool ProvidesError() const { return true; }

   /// return errors at the minimum
   virtual const double * Errors() const { return  (fErrors.empty()) ? 0 : &fErrors.front(); }

   /** return covariance matrices elements
       if the variable is fixed the matrix is zero
       The ordering of the variables is the same as in errors
   */
   virtual double CovMatrix(unsigned int i, unsigned int j) const {
      return (fCovar.empty()) ? 0 : fCovar[i + fDim* j];
   }

   /// return covariance matrix status
   virtual int CovMatrixStatus() const {
      if (fCovar.size() == 0) return 0;
      return (fStatus ==0) ? 3 : 1;
   }

   /// return reference to the objective function
   ///virtual const ROOT::Math::IGenFunction & Function() const;




protected:

private:

   bool fRobust;
   unsigned int fDim;
   unsigned int fNFree;
   double fMinVal;
   std::vector<double> fParams;
   std::vector<double> fErrors;
   std::vector<double> fCovar;

   const ROOT::Math::IMultiGradFunction * fObjFunc;
   TLinearFitter * fFitter;

   ClassDef(TLinearMinimizer,1)  //Implementation of the Minimizer interface using TLinearFitter

};



#endif /* ROOT_TLinearMinimizer */
