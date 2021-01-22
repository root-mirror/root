/**
 * CMA-ES, Covariance Matrix Evolution Strategy
 * Copyright (c) 2014 INRIA
 * Author: Emmanuel Benazera <emmanuel.benazera@lri.fr>
 *
 * This file is part of libcmaes.
 *
 * libcmaes is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libcmaes is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libcmaes.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ROOT_CMAES_CMAESMinimizer
#define ROOT_CMAES_CMAESMinimizer

#include "Math/Minimizer.h"

#include "Math/IFunctionfwd.h"

#include "libcmaes/cmaes.h"

using namespace libcmaes;

namespace ROOT
{
  namespace cmaes
  {

    class TCMAESMinimizer : public ROOT::Math::Minimizer
    {
    public:
      
      /**
       * Default constructor
       */
      TCMAESMinimizer(); // CMAES_DEFAULT

      /**
       * Constructor with a char (used by PM)
       */
      TCMAESMinimizer(const char *type);

      /**
       * Destructor
       */
      virtual ~TCMAESMinimizer();

    private:
      // blocking copy constructor.
      TCMAESMinimizer(const TCMAESMinimizer &m);
      
      /**
       * assignment operator.
       */
      TCMAESMinimizer& operator = (const TCMAESMinimizer &rhs);

    public:

      // clear resources (parameters) for consecutives minimizations
      virtual void Clear();
      
      /// set the function to minimize
      virtual void SetFunction(const ROOT::Math::IMultiGenFunction & func); 
      
      /// set gradient the function to minimize
      virtual void SetFunction(const ROOT::Math::IMultiGradFunction & func); 
      
      /// set free variable 
      virtual bool SetVariable(unsigned int ivar, const std::string & name, double val, double step); 
      
      /// set lower limit variable  (override if minimizer supports them )
      virtual bool SetLowerLimitedVariable(unsigned int  ivar , const std::string & name , double val , double step , double lower );
      /// set upper limit variable (override if minimizer supports them )
      virtual bool SetUpperLimitedVariable(unsigned int ivar , const std::string & name , double val , double step , double upper ); 
      /// set upper/lower limited variable (override if minimizer supports them )
      virtual bool SetLimitedVariable(unsigned int ivar , const std::string & name , double val , double step , double /* lower */, double /* upper */); 
      /// set fixed variable (override if minimizer supports them )
      virtual bool SetFixedVariable(unsigned int /* ivar */, const std::string & /* name */, double /* val */);  
      /// set variable
      virtual bool SetVariableValue(unsigned int ivar, double val);
      // set variable values
      virtual bool SetVariableValues(const double * val);
      /// set the step size of an already existing variable
      virtual bool SetVariableStepSize(unsigned int ivar, double step );
      /// set the lower-limit of an already existing variable 
      virtual bool SetVariableLowerLimit(unsigned int ivar, double lower);
      /// set the upper-limit of an already existing variable 
      virtual bool SetVariableUpperLimit(unsigned int ivar, double upper);
      /// set the limits of an already existing variable 
      virtual bool SetVariableLimits(unsigned int ivar, double lower, double upper); 
      /// fix an existing variable
      virtual bool FixVariable(unsigned int ivar);
      /// release an existing variable
      //virtual bool ReleaseVariable(unsigned int ivar);
      /// query if an existing variable is fixed (i.e. considered constant in the minimization)
      /// note that by default all variables are not fixed 
      virtual bool IsFixedVariable(unsigned int ivar)  const;
      /// get variable settings in a variable object (like ROOT::Fit::ParamsSettings)
      virtual bool GetVariableSettings(unsigned int ivar, ROOT::Fit::ParameterSettings & varObj) const;
      /// get name of variables (override if minimizer support storing of variable names)
      virtual std::string VariableName(unsigned int ivar) const;
      /// get index of variable given a variable given a name
      /// return -1 if variable is not found
      virtual int VariableIndex(const std::string & name) const;

      // sets libcmaes parameter object based on options.
      template <class TGenoPheno>
	void SetMParameters(CMAParameters<TGenoPheno> &cmaparams,
			    const int &maxiter, const int &maxfevals,
			    const int &noisy, const int &nrestarts,
			    const double &ftarget,
			    const std::string &fplot,
			    const bool &withnumgradient,
			    const bool &mtfeval,
			    const bool &quiet,
			    const int &elitist, // 0: off (default), then 1, 2 or 3
			    const bool &uh);
      
      /** 
	  method to perform the minimization. 
	  Return false in case the minimization did not converge. In this case a 
	  status code different than zero is set 
	  (retrieved by the derived method Minimizer::Status() )" 
	  
	  status = 1    : Covariance was made pos defined -> not yet implemented in libcmaes
	  status = 2    : Hesse is invalid => Not Applicable, kept for correspondence with Minuit2
	  status = 3    : Edm is above max 
	  status = 4    : Reached call limit
	  status = 5    : Any other failure 
      */
      virtual bool Minimize(); 
      
      /// return minimum function value
      virtual double MinValue() const;
      
      /// return expected distance reached from the minimum
      virtual double Edm() const;
      
      /// return  pointer to X values at the minimum 
      virtual const double *  X() const; 
      
      /// return pointer to gradient values at the minimum 
      virtual const double *  MinGradient() const { return 0; } // not available in CMAES 
      
      /// number of function calls to reach the minimum 
      virtual unsigned int NCalls() const;
      
      /// this is <= Function().NDim() which is the total 
      /// number of variables (free+ constrained ones) 
      virtual unsigned int NDim() const { return fDim; }   
      
      /// number of free variables (real dimension of the problem) 
      /// this is <= Function().NDim() which is the total 
      virtual unsigned int NFree() const { return fFreeDim; }  
      
      /// minimizer provides error and error matrix
      virtual bool ProvidesError() const { return true; } 
      
      /// return errors at the minimum
      virtual const double* Errors() const; 
      
      /** 
	  return covariance matrix elements 
	  if the variable is fixed or const the value is zero
	  The ordering of the variables is the same as in errors and parameter value. 
	  This is different from the direct interface of Minuit2 or TMinuit where the 
	  values were obtained only to variable parameters
      */ 
      virtual double CovMatrix(unsigned int i, unsigned int j) const;  
      
      
      /** 
	  Fill the passed array with the  covariance matrix elements 
	  if the variable is fixed or const the value is zero. 
	  The array will be filled as cov[i *ndim + j]
	  The ordering of the variables is the same as in errors and parameter value. 
	  This is different from the direct interface of Minuit2 or TMinuit where the 
	  values were obtained only to variable parameters
      */ 
      virtual bool GetCovMatrix(double * cov) const;  
      
      /** 
	  Fill the passed array with the Hessian matrix elements 
	  The Hessian matrix is the matrix of the second derivatives 
	  and is the inverse of the covariance matrix
	  If the variable is fixed or const the values for that variables are zero. 
	  The array will be filled as h[i *ndim + j]
      */ 
      //virtual bool GetHessianMatrix(double * h) const;  
      
      
      /**
	 return the status of the covariance matrix 
	 status = -1 :  not available (inversion failed or Hesse failed) 
	 status =  0 : available but not positive defined
	 status =  1 : covariance only approximate
	 status =  2 : full matrix but forced pos def 
	 status =  3 : full accurate matrix 
	 
      */
      virtual int CovMatrixStatus() const { return 3; }
      /**
	 return correlation coefficient between variable i and j.
	 If the variable is fixed or const the return value is zero
      */
      virtual double Correlation(unsigned int i, unsigned int j ) const; 
      
      /**
	 get global correlation coefficient for the variable i. This is a number between zero and one which gives 
	 the correlation between the i-th variable  and that linear combination of all other variables which 
	 is most strongly correlated with i.
	 If the variable is fixed or const the return value is zero
      */
      virtual double GlobalCC(unsigned int i) const;
      
      /**
	 get the minos error for parameter i, return false if Minos failed
	 A minimizaiton must be performed befre, return false if no minimization has been done
	 In case of Minos failed the status error is updated as following 
	 status += 10 * minosStatus where the minos status is:
	 status = 1    : maximum number of function calls exceeded when running for lower error
	 status = 2    : maximum number of function calls exceeded when running for upper error
	 status = 3    : new minimum found when running for lower error
	 status = 4    : new minimum found when running for upper error
	 status = 5    : any other failure 
	 
      */
      virtual bool GetMinosError(unsigned int i, double & errLow, double & errUp, int = 0); 
      
      /**
	 scan a parameter i around the minimum. A minimization must have been done before, 
	 return false if it is not the case
      */
      virtual bool Scan(unsigned int i, unsigned int & nstep, double * x, double * y, double xmin = 0, double xmax = 0); 
      
      /**
	 find the contour points (xi,xj) of the function for parameter i and j around the minimum
	 The contour will be find for value of the function = Min + ErrorUp();
      */
      virtual bool Contour(unsigned int i, unsigned int j, unsigned int & npoints, double *xi, double *xj); 
      
      
      /**
	 perform a full calculation of the Hessian matrix for error calculation
	 If a valid minimum exists the calculation is done on the minimum point otherwise is performed 
	 in the current set values of parameters
	 Status code of minimizer is updated according to the following convention (in case Hesse failed)
	 status += 100*hesseStatus where hesse status is: 
	 status = 1 : hesse failed
	 status = 2 : matrix inversion failed
	 status = 3 : matrix is not pos defined 
      */
      //virtual bool Hesse();
      
      
      /// return reference to the objective function
      //virtual const ROOT::Math::IGenFunction & Function() const; 
      
      /// print result of minimization
      virtual void PrintResults(); 
      
      /// set an object to trace operation for each iteration
      /// The object muust implement operator() (unsigned int, MinimumState & state)
      //void SetTraceObject(MnTraceObject & obj);
      
      /// set storage level   = 1 : store all iteration states (default)
      ///                     = 0 : store only first and last state to save memory
      // N/A
      void SetStorageLevel(int level); 
      
    private: 
      unsigned int fDim = 0;       // dimension of the function to be minimized 
      unsigned int fFreeDim = 0; // Number of free dimensions.
      std::string fMinimizer = "cmaes"; // minimizer algo.
      const ROOT::Math::IMultiGenFunction *fObjFunc = nullptr;
      const ROOT::Math::IMultiGradFunction *fObjFuncGrad = nullptr;
      std::vector<double> fLBounds; // Lower bounds of variables
      std::vector<double> fUBounds; // Upper bounds of variables
      std::vector<int> fVariablesType; // 0 for free variable, 1 for fixed variable, 2 for lower bounded, 3 for upper bounded, 4 for lower and upper bounded.
      std::vector<double> fInitialX;
      std::vector<std::string> fNames; // Names of the variables.
      std::vector<double> fInitialSigma; // User-set Initial step-size for each variables.
      std::map<int,double> fFixedVariables; // fixed variables and values.
      CMASolutions fCMAsols;
      CMAParameters<GenoPheno<NoBoundStrategy,NoScalingStrategy>> fCMAparams; // params no bounds.
      CMAParameters<GenoPheno<pwqBoundStrategy,NoScalingStrategy>> fCMAparams_b; // params with bounds.
      CMAParameters<GenoPheno<NoBoundStrategy,linScalingStrategy>> fCMAparams_l; // params no bounds + linear scaling.
      CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>> fCMAparams_lb; // params with bounds + linear scaling.
      mutable std::vector<double> fGlobalCC; // vector of global correlation coefficients.
      mutable std::vector<double> fValues; // X values.
      mutable std::vector<double> fErrors; // X errors.
      bool fWithBounds = false; // whether using box-type constraints as required by parameters.
      bool fWithGradient = false; // whether to use gradient information when available.
      int fWithLinearScaling = 0; // wheter to use linear scaling of objective function parameters. */
    };
    
  }  // end namespace cmaes
} // end namespace ROOT

#endif
