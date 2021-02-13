// @(#)root/fumili:$Id$
// Author: L. Moneta Wed Oct 25 16:28:55 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class TFumiliMinimizer

#ifndef ROOT_TFumiliMinimizer
#define ROOT_TFumiliMinimizer

#include "Math/Minimizer.h"

#include "Math/FitMethodFunction.h"

#include "Rtypes.h"
#include <vector>
#include <string>

class TFumili;



// namespace ROOT {

//    namespace Math {

//       class BasicFitMethodFunction<ROOT::Math::IMultiGenFunction>;
//       class BasicFitMethodFunction<ROOT::Math::IMultiGradFunction>;

//    }
// }



/**
   TFumiliMinimizer class: minimizer implementation based on TFumili.
*/
class TFumiliMinimizer  : public ROOT::Math::Minimizer {

public:

   /**
      Default constructor (an argument is needed by plug-in manager)
   */
   TFumiliMinimizer (int dummy=0 );


   /**
      Destructor (no operations)
   */
   ~TFumiliMinimizer ();

private:
   // usually copying is non trivial, so we make this unaccessible

   /**
      Copy constructor
   */
   TFumiliMinimizer(const TFumiliMinimizer &);

   /**
      Assignment operator
   */
   TFumiliMinimizer & operator = (const TFumiliMinimizer & rhs);

public:

   /// set the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGenFunction & func);

   /// set the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGradFunction & func);

   /// set free variable
   virtual bool SetVariable(unsigned int ivar, const std::string & name, double val, double step);

   /// set upper/lower limited variable (override if minimizer supports them )
   virtual bool SetLimitedVariable(unsigned int ivar , const std::string & name , double val , double step , double /* lower */, double /* upper */);

#ifdef LATER
   /// set lower limit variable  (override if minimizer supports them )
   virtual bool SetLowerLimitedVariable(unsigned int  ivar , const std::string & name , double val , double step , double lower );
   /// set upper limit variable (override if minimizer supports them )
   virtual bool SetUpperLimitedVariable(unsigned int ivar , const std::string & name , double val , double step , double upper );
#endif

   /// set fixed variable (override if minimizer supports them )
   virtual bool SetFixedVariable(unsigned int /* ivar */, const std::string & /* name */, double /* val */);

   /// set the value of an existing variable
   virtual bool SetVariableValue(unsigned int ivar, double val );

   /// method to perform the minimization
   virtual  bool Minimize();

   /// return minimum function value
   virtual double MinValue() const { return fMinVal; }

   /// return expected distance reached from the minimum
   virtual double Edm() const { return fEdm; }

   /// return  pointer to X values at the minimum
   virtual const double *  X() const { return &fParams.front(); }

   /// return pointer to gradient values at the minimum
   virtual const double *  MinGradient() const { return 0; } // not available

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
   virtual const double * Errors() const { return  &fErrors.front(); }

   /** return covariance matrices elements
       if the variable is fixed the matrix is zero
       The ordering of the variables is the same as in errors
   */
   virtual double CovMatrix(unsigned int i, unsigned int j) const {
      return fCovar[i + fDim* j];
   }

   /*
     return covariance matrix status
   */
   virtual int CovMatrixStatus() const {
      if (fCovar.size() == 0) return 0;
      return (fStatus ==0) ? 3 : 1;
   }





protected:

   /// implementation of FCN for Fumili
   static void Fcn( int &, double * , double & f, double * , int);
   /// implementation of FCN for Fumili when user provided gradient is used
   //static void FcnGrad( int &, double * g, double & f, double * , int);

   /// static function implementing the evaluation of the FCN (it uses static instance fgFumili)
   static double EvaluateFCN(const double * x, double * g);

private:


   unsigned int fDim;
   unsigned int fNFree;
   double fMinVal;
   double fEdm;
   std::vector<double> fParams;
   std::vector<double> fErrors;
   std::vector<double> fCovar;

   TFumili * fFumili;

   // need to have a static copy of the function
   //NOTE: This is NOT thread safe.
   static ROOT::Math::FitMethodFunction * fgFunc;
   static ROOT::Math::FitMethodGradFunction * fgGradFunc;

   static TFumili * fgFumili; // static instance (used by fcn function)

   ClassDef(TFumiliMinimizer,1)  //Implementation of Minimizer interface using TFumili

};



#endif /* ROOT_TFumiliMinimizer */
