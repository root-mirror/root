// @(#)root/mathcore:$Id$
// Author: L. Moneta Fri Sep 22 15:06:47 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class Minimizer

#ifndef ROOT_Math_Minimizer
#define ROOT_Math_Minimizer

#include "Math/IFunction.h"

#include "Math/MinimizerOptions.h"

#include "Math/Util.h"

#include "Math/Error.h"


#include <string>
#include <limits>
#include <cmath>


namespace ROOT {

   namespace Fit {
      class ParameterSettings;
   }


   namespace Math {

/**
   @defgroup MultiMin Multi-dimensional Minimization
   @ingroup NumAlgo

   Classes implementing algorithms for multi-dimensional minimization
 */



//_______________________________________________________________________________
/**
   Abstract Minimizer class, defining  the interface for the various minimizer
   (like Minuit2, Minuit, GSL, etc..)
   Plug-in's exist in ROOT to be able to instantiate the derived classes like
   ROOT::Math::GSLMinimizer or ROOT::Math::Minuit2Minimizer via the
   plug-in manager.

   Provides interface for setting the function to be minimized.
   The function must  implemente the multi-dimensional generic interface
   ROOT::Math::IBaseFunctionMultiDim.
   If the function provides gradient calculation
   (implements the ROOT::Math::IGradientFunctionMultiDim interface) this will be
   used by the Minimizer.

   It Defines also interface for setting the initial values for the function variables (which are the parameters in
   of the model function in case of solving for fitting) and especifying their limits.

   It defines the interface to set and retrieve basic minimization parameters
   (for specific Minimizer parameters one must use the derived classes).

   Then it defines the interface to retrieve the result of minimization ( minimum X values, function value,
   gradient, error on the mimnimum, etc...)

   @ingroup MultiMin
*/

class Minimizer {

public:

   /**
      Default constructor
   */
   Minimizer () :
      fValidError(false),
      fStatus(-1)
   {}

   /**
      Destructor (no operations)
   */
   virtual ~Minimizer ()  {}




private:
   // usually copying is non trivial, so we make this unaccessible

   /**
      Copy constructor
   */
   Minimizer(const Minimizer &) {}

   /**
      Assignment operator
   */
   Minimizer & operator = (const Minimizer & rhs)  {
      if (this == &rhs) return *this;  // time saving self-test
      return *this;
   }

public:

   /// reset for consecutive minimizations - implement if needed
   virtual void Clear() {}

   /// set the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGenFunction & func) = 0;

   /// set a function to minimize using gradient
   virtual void SetFunction(const ROOT::Math::IMultiGradFunction & func)
   {
      SetFunction(static_cast<const ::ROOT::Math::IMultiGenFunction &> (func));
   }


   /// add variables  . Return number of variables successfully added
   template<class VariableIterator>
   int SetVariables(const VariableIterator & begin, const VariableIterator & end) {
      unsigned int ivar = 0;
      for ( VariableIterator vitr = begin; vitr != end; ++vitr) {
         bool iret = false;
         if (vitr->IsFixed() )
            iret = SetFixedVariable(ivar,  vitr->Name(), vitr->Value() );
         else if (vitr->IsDoubleBound() )
            iret = SetLimitedVariable(ivar,  vitr->Name(), vitr->Value(), vitr->StepSize(), vitr->LowerLimit(), vitr->UpperLimit() );
         else if (vitr->HasLowerLimit() )
            iret = SetLowerLimitedVariable(ivar,  vitr->Name(), vitr->Value(), vitr->StepSize(), vitr->LowerLimit() );
         else if (vitr->HasUpperLimit() )
            iret = SetUpperLimitedVariable(ivar,  vitr->Name(), vitr->Value(), vitr->StepSize(), vitr->UpperLimit() );
         else
            iret = SetVariable( ivar, vitr->Name(), vitr->Value(), vitr->StepSize() );

         if (iret) ivar++;

         // an error message should be eventually be reported in the virtual single SetVariable methods
      }
      return ivar;
   }
   /// set a new free variable
   virtual bool SetVariable(unsigned int ivar, const std::string & name, double val, double step) = 0;
   /// set a new lower limit variable  (override if minimizer supports them )
   virtual bool SetLowerLimitedVariable(unsigned int  ivar , const std::string & name , double val , double step , double lower ) {
      return SetLimitedVariable(ivar, name, val, step, lower, std::numeric_limits<double>::infinity() );
   }
   /// set a new upper limit variable (override if minimizer supports them )
   virtual bool SetUpperLimitedVariable(unsigned int ivar , const std::string & name , double val , double step , double upper ) {
      return SetLimitedVariable(ivar, name, val, step, - std::numeric_limits<double>::infinity(), upper );
   }
   /// set a new upper/lower limited variable (override if minimizer supports them ) otherwise as default set an unlimited variable
   virtual bool SetLimitedVariable(unsigned int ivar  , const std::string & name  , double val  , double  step ,
                                   double lower , double  upper ) {
      MATH_WARN_MSG("Minimizer::SetLimitedVariable","Setting of limited variable not implemented - set as unlimited");
      MATH_UNUSED(lower); MATH_UNUSED(upper);
      return SetVariable(ivar, name, val, step);
   }
   /// set a new fixed variable (override if minimizer supports them )
   virtual bool SetFixedVariable(unsigned int  ivar  , const std::string &  name , double val  ) {
      MATH_ERROR_MSG("Minimizer::SetFixedVariable","Setting of fixed variable not implemented");
      MATH_UNUSED(ivar); MATH_UNUSED(name); MATH_UNUSED(val);
      return false;
   }
   /// set the value of an already existing variable
   virtual bool SetVariableValue(unsigned int ivar , double value) {
      MATH_ERROR_MSG("Minimizer::SetVariableValue","Set of a variable value not implemented");
      MATH_UNUSED(ivar); MATH_UNUSED(value);
      return false;
   }
   /// set the values of all existing variables (array must be dimensioned to the size of the existing parameters)
   virtual bool SetVariableValues(const double * x) {
      bool ret = true;
      unsigned int i = 0;
      while ( i <= NDim() && ret) {
         ret &= SetVariableValue(i,x[i] ); i++;
      }
      return ret;
   }
   /// set the step size of an already existing variable
   virtual bool SetVariableStepSize(unsigned int ivar, double value ) {
      MATH_ERROR_MSG("Minimizer::SetVariableStepSize","Setting an existing variable step size not implemented");
      MATH_UNUSED(ivar); MATH_UNUSED(value);
      return false;
   }
   /// set the lower-limit of an already existing variable
   virtual bool SetVariableLowerLimit(unsigned int ivar, double lower) {
      MATH_ERROR_MSG("Minimizer::SetVariableLowerLimit","Setting an existing variable limit not implemented");
      MATH_UNUSED(ivar); MATH_UNUSED(lower);
      return false;
   }
   /// set the upper-limit of an already existing variable
   virtual bool SetVariableUpperLimit(unsigned int ivar, double upper) {
      MATH_ERROR_MSG("Minimizer::SetVariableUpperLimit","Setting an existing variable limit not implemented");
      MATH_UNUSED(ivar); MATH_UNUSED(upper);
      return false;
   }
   /// set the limits of an already existing variable
   virtual bool SetVariableLimits(unsigned int ivar, double lower, double upper) {
      return SetVariableLowerLimit(ivar,lower) && SetVariableUpperLimit(ivar,upper);
   }
   /// fix an existing variable
   virtual bool FixVariable(unsigned int ivar) {
      MATH_ERROR_MSG("Minimizer::FixVariable","Fixing an existing variable not implemented");
      MATH_UNUSED(ivar);
      return false;
   }
   /// release an existing variable
   virtual bool ReleaseVariable(unsigned int ivar) {
      MATH_ERROR_MSG("Minimizer::ReleaseVariable","Releasing an existing variable not implemented");
      MATH_UNUSED(ivar);
      return false;
   }
   /// query if an existing variable is fixed (i.e. considered constant in the minimization)
   /// note that by default all variables are not fixed
   virtual bool IsFixedVariable(unsigned int ivar) const {
      MATH_ERROR_MSG("Minimizer::IsFixedVariable","Quering an existing variable not implemented");
      MATH_UNUSED(ivar);
      return false;
   }
   /// get variable settings in a variable object (like ROOT::Fit::ParamsSettings)
   virtual bool GetVariableSettings(unsigned int ivar, ROOT::Fit::ParameterSettings & pars) const {
      MATH_ERROR_MSG("Minimizer::GetVariableSettings","Quering an existing variable not implemented");
      MATH_UNUSED(ivar); MATH_UNUSED(pars);
      return false;
   }


   /// set the initial range of an existing variable
   virtual bool SetVariableInitialRange(unsigned int /* ivar */, double /* mininitial */, double /* maxinitial */) {
     return false;
   }

   /// method to perform the minimization
   virtual  bool Minimize() = 0;

   /// return minimum function value
   virtual double MinValue() const = 0;

   /// return  pointer to X values at the minimum
   virtual const double *  X() const = 0;

   /// return expected distance reached from the minimum (re-implement if minimizer provides it
   virtual double Edm() const { return -1; }

   /// return pointer to gradient values at the minimum
   virtual const double *  MinGradient() const { return NULL; }

   /// number of function calls to reach the minimum
   virtual unsigned int NCalls() const { return 0; }

   /// number of iterations to reach the minimum
   virtual unsigned int NIterations() const { return NCalls(); }

   /// this is <= Function().NDim() which is the total
   /// number of variables (free+ constrained ones)
   virtual unsigned int NDim() const = 0;

   /// number of free variables (real dimension of the problem)
   /// this is <= Function().NDim() which is the total
   /// (re-implement if minimizer supports bounded parameters)
   virtual unsigned int NFree() const { return NDim(); }

   /// minimizer provides error and error matrix
   virtual bool ProvidesError() const { return false; }

   /// return errors at the minimum
   virtual const double * Errors() const { return NULL; }

   /** return covariance matrices element for variables ivar,jvar
       if the variable is fixed the return value is zero
       The ordering of the variables is the same as in the parameter and errors vectors
   */
   virtual double CovMatrix(unsigned int  ivar , unsigned int jvar ) const {
      MATH_UNUSED(ivar); MATH_UNUSED(jvar);
      return 0;
   }

   /**
       Fill the passed array with the  covariance matrix elements
       if the variable is fixed or const the value is zero.
       The array will be filled as cov[i *ndim + j]
       The ordering of the variables is the same as in errors and parameter value.
       This is different from the direct interface of Minuit2 or TMinuit where the
       values were obtained only to variable parameters
   */
   virtual bool GetCovMatrix(double * covMat) const {
      MATH_UNUSED(covMat);
      return false;
   }

   /**
       Fill the passed array with the Hessian matrix elements
       The Hessian matrix is the matrix of the second derivatives
       and is the inverse of the covariance matrix
       If the variable is fixed or const the values for that variables are zero.
       The array will be filled as h[i *ndim + j]
   */
   virtual bool GetHessianMatrix(double * hMat) const {
      MATH_UNUSED(hMat);
      return false;
   }


   ///return status of covariance matrix
   /// using Minuit convention {0 not calculated 1 approximated 2 made pos def , 3 accurate}
   /// Minimizer who implements covariance matrix calculation will re-implement the method
   virtual int CovMatrixStatus() const {
      return 0;
   }

   /**
      return correlation coefficient between variable i and j.
      If the variable is fixed or const the return value is zero
    */
   virtual double Correlation(unsigned int i, unsigned int j ) const {
      double tmp = CovMatrix(i,i) * CovMatrix(j,j);
      return ( tmp < 0) ? 0 : CovMatrix(i,j) / std::sqrt( tmp );
   }

   /**
      return global correlation coefficient for variable i
      This is a number between zero and one which gives
      the correlation between the i-th parameter  and that linear combination of all
      other parameters which is most strongly correlated with i.
      Minimizer must overload method if implemented
    */
   virtual double GlobalCC(unsigned int ivar) const {
      MATH_UNUSED(ivar);
      return -1;
   }

   /**
      minos error for variable i, return false if Minos failed or not supported
      and the lower and upper errors are returned in errLow and errUp
      An extra flag  specifies if only the lower (option=-1) or the upper (option=+1) error calculation is run
   */
   virtual bool GetMinosError(unsigned int ivar , double & errLow, double & errUp, int option = 0) {
      MATH_ERROR_MSG("Minimizer::GetMinosError","Minos Error not implemented");
      MATH_UNUSED(ivar); MATH_UNUSED(errLow); MATH_UNUSED(errUp); MATH_UNUSED(option);
      return false;
   }

   /**
      perform a full calculation of the Hessian matrix for error calculation
    */
   virtual bool Hesse() {
      MATH_ERROR_MSG("Minimizer::Hesse","Hesse not implemented");
      return false;
   }

   /**
      scan function minimum for variable i. Variable and function must be set before using Scan
      Return false if an error or if minimizer does not support this functionality
    */
   virtual bool Scan(unsigned int ivar , unsigned int & nstep , double * x , double * y ,
                     double xmin = 0, double xmax = 0) {
      MATH_ERROR_MSG("Minimizer::Scan","Scan not implemented");
      MATH_UNUSED(ivar); MATH_UNUSED(nstep); MATH_UNUSED(x); MATH_UNUSED(y);
      MATH_UNUSED(xmin); MATH_UNUSED(xmax);
      return false;
   }

   /**
      find the contour points (xi, xj) of the function for parameter ivar and jvar around the minimum
      The contour will be find for value of the function = Min + ErrorUp();
    */
   virtual bool Contour(unsigned int ivar , unsigned int jvar, unsigned int & npoints,
                        double *  xi , double * xj ) {
      MATH_ERROR_MSG("Minimizer::Contour","Contour not implemented");
      MATH_UNUSED(ivar); MATH_UNUSED(jvar); MATH_UNUSED(npoints);
      MATH_UNUSED(xi); MATH_UNUSED(xj);
      return false;
   }

   /// return reference to the objective function
   ///virtual const ROOT::Math::IGenFunction & Function() const = 0;

   /// print the result according to set level (implemented for TMinuit for mantaining Minuit-style printing)
   virtual void PrintResults() {}

   /// get name of variables (override if minimizer support storing of variable names)
   /// return an empty string if variable is not found
   virtual std::string VariableName(unsigned int ivar) const {
      MATH_UNUSED(ivar);
      return std::string(); // return empty string
   }

   /// get index of variable given a variable given a name
   /// return -1 if variable is not found
   virtual int VariableIndex(const std::string & name) const {
      MATH_ERROR_MSG("Minimizer::VariableIndex","Getting variable index from name not implemented");
      MATH_UNUSED(name);
      return -1;
   }

   /** minimizer configuration parameters **/

   /// set print level
   int PrintLevel() const { return fOptions.PrintLevel(); }

   ///  max number of function calls
   unsigned int MaxFunctionCalls() const { return fOptions.MaxFunctionCalls(); }

   /// max iterations
   unsigned int MaxIterations() const { return fOptions.MaxIterations(); }

   /// absolute tolerance
   double Tolerance() const { return  fOptions.Tolerance(); }

   /// precision of minimizer in the evaluation of the objective function
   /// ( a value <=0 corresponds to the let the minimizer choose its default one)
   double Precision() const { return fOptions.Precision(); }

   /// strategy
   int Strategy() const { return fOptions.Strategy(); }

   /// status code of minimizer
   int Status() const { return fStatus; }

   /// status code of Minos (to be re-implemented by the minimizers supporting Minos)
   virtual int MinosStatus() const { return -1; }

   /// return the statistical scale used for calculate the error
   /// is typically 1 for Chi2 and 0.5 for likelihood minimization
   double ErrorDef() const { return fOptions.ErrorDef(); }

   ///return true if Minimizer has performed a detailed error validation (e.g. run Hesse for Minuit)
   bool IsValidError() const { return fValidError; }

   /// retrieve the minimizer options (implement derived class if needed)
   virtual MinimizerOptions  Options() const {
      return fOptions;
   }

   /// set print level
   void SetPrintLevel(int level) { fOptions.SetPrintLevel(level); }

   ///set maximum of function calls
   void SetMaxFunctionCalls(unsigned int maxfcn) { if (maxfcn > 0) fOptions.SetMaxFunctionCalls(maxfcn); }

   /// set maximum iterations (one iteration can have many function calls)
   void SetMaxIterations(unsigned int maxiter) { if (maxiter > 0) fOptions.SetMaxIterations(maxiter); }

   /// set the tolerance
   void SetTolerance(double tol) { fOptions.SetTolerance(tol); }

   /// set in the minimizer the objective function evaluation precision
   /// ( a value <=0 means the minimizer will choose its optimal value automatically, i.e. default case)
   void SetPrecision(double prec) { fOptions.SetPrecision(prec); }

   ///set the strategy
   void SetStrategy(int strategyLevel) { fOptions.SetStrategy(strategyLevel); }

   /// set scale for calculating the errors
   void SetErrorDef(double up) { fOptions.SetErrorDef(up); }

   /// flag to check if minimizer needs to perform accurate error analysis (e.g. run Hesse for Minuit)
   void SetValidError(bool on) { fValidError = on; }

   /// set all options in one go
   void SetOptions(const MinimizerOptions & opt) {
      fOptions = opt;
   }

   /// reset the defaut options (defined in MinimizerOptions)
   void SetDefaultOptions() {
      fOptions.ResetToDefaultOptions();
   }

protected:


//private:


   // keep protected to be accessible by the derived classes


   bool fValidError;            // flag to control if errors have been validated (Hesse has been run in case of Minuit)
   MinimizerOptions fOptions;   // minimizer options
   int fStatus;                 // status of minimizer
};

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_Minimizer */
