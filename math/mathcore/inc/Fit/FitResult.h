// @(#)root/mathcore:$Id$
// Author: L. Moneta Wed Aug 30 11:05:34 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class FitResult

#ifndef ROOT_Fit_FitResult
#define ROOT_Fit_FitResult

#include "Math/IFunctionfwd.h"
#include "Math/IParamFunctionfwd.h"

#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <memory>

namespace ROOT {

   namespace Math {
      class Minimizer;
   }


   namespace Fit {

      class FitConfig;
      class FitData;
      class BinData;

//___________________________________________________________________________________
/**
   class containg the result of the fit and all the related information
   (fitted parameter values, error, covariance matrix and minimizer result information)
   Contains a pointer also to the fitted (model) function, modified with the fit parameter values.
   When the fit is valid, it is constructed from a  Minimizer and a model function pointer

   @ingroup FitMain
*/
class FitResult {

public:

   typedef  ROOT::Math::IParamMultiFunction IModelFunction;

   /**
      Default constructor for an empty (non valid) fit result
   */
   FitResult ();

   /**
      Constructor from a fit-config for a dummy fit
      (e.g. when only one fcn evaluation is done)
   */
   FitResult (const FitConfig & fconfig);


   /**
      Copy constructor.
   */
   FitResult(const FitResult & rhs);

   /**
      Assignment operator
   */
   FitResult & operator = (const FitResult & rhs);

   /**
      Destructor
   */
   virtual ~FitResult ();


public:

   /**
      Fill the fit result from a Minimizer instance after fitting
      Run also Minos if requested from the configuration
    */
   void FillResult(const std::shared_ptr<ROOT::Math::Minimizer> & min, const FitConfig & fconfig,  const std::shared_ptr<IModelFunction> & f,
              bool isValid, unsigned int sizeOfData = 0, bool binFit = true, const ROOT::Math::IMultiGenFunction * chi2func = 0, unsigned int ncalls = 0);


   /**
      Update the fit result with a new minimization status
      To be run only if same fit is performed with same configuration
      Note that in this case MINOS is not re-run. If one wants to run also MINOS
      a new result must be created
    */
   bool Update(const std::shared_ptr<ROOT::Math::Minimizer> & min, const ROOT::Fit::FitConfig & fconfig, bool isValid, unsigned int ncalls = 0);

   /** minimization quantities **/

   /// minimizer type
   const std::string & MinimizerType() const { return fMinimType; }

   /**
       True if fit successful, otherwise false.
       A fit is considered successful if the minimizer succeded in finding the
       minimum. It could happen that subsequent operations like error analysis (e.g. Minos)
       failed. In that case the status can be still true if the original minimization algorithm
       succeeded in finding the minimum.
       One can query in that case the minimizer return status using Status().
       It is responability to the Minimizer class to tag a found minimum as valid or not
       and to produce also a status code.
   */
   bool IsValid() const { return fValid; }

   /// True if a fit result does not exist (even invalid) with parameter values
   bool IsEmpty() const { return (fParams.size() == 0);  }

   /// Return value of the objective function (chi2 or likelihood) used in the fit
   double MinFcnValue() const { return fVal; }

   ///Number of function calls to find minimum
   unsigned int NCalls() const { return fNCalls; }

   ///Expected distance from minimum
   double Edm() const { return fEdm; }

   ///   get total number of parameters
   unsigned int NTotalParameters() const { return fParams.size(); }
   /// total number of parameters (abbreviation)
   unsigned int NPar() const { return NTotalParameters(); }

   /// get total number of free parameters
   unsigned int NFreeParameters() const { return fNFree; }

   /// minimizer status code
   int Status() const { return fStatus; }

   ///covariance matrix status code
   /// using Minuit convention : =0 not calculated, =1 approximated, =2 made pos def , =3 accurate

   int CovMatrixStatus() const { return fCovStatus; }

   /** fitting quantities **/

   /// Return pointer to model (fit) function with fitted parameter values.
   /// Pointer is managed internally. I must not be deleted
   const IModelFunction * FittedFunction() const {
      return fFitFunc.get();
   }

   /// return BinData used in the fit (return a nullptr in case a different fit is done
   /// or the data are not available
   /// Pointer is managed internally, it must not be deleted
   const BinData * FittedBinData() const;


   /// Chi2 fit value
   /// in case of likelihood must be computed ?
   double Chi2() const { return fChi2; }

   /// Number of degree of freedom
   unsigned int Ndf() const { return fNdf; }

   /// p value of the fit (chi2 probability)
   double Prob() const;

   /// parameter errors (return st::vector)
   const std::vector<double> & Errors() const { return fErrors; }
   /// parameter errors (return const pointer)
   const double * GetErrors() const { return (fErrors.empty()) ? 0 : &fErrors.front(); }

   /// parameter values (return std::vector)
   const std::vector<double> & Parameters() const { return fParams; }
   /// parameter values (return const pointer)
   const double * GetParams() const { return &fParams.front(); }

   /// parameter value by index
   double Value(unsigned int i) const { return fParams[i]; }
   /// parameter value by index
   double Parameter(unsigned int i) const { return fParams[i]; }

   /// parameter error by index
   // (NOTE: this due to conflict with TObject::Error cannot used in derived class which
   // inherits from TObject. Use instead ParError (or Errors()[i] )
   double Error(unsigned int i) const {
      return (i < fErrors.size() ) ? fErrors[i] : 0;
   }
   /// parameter error by index
   double ParError(unsigned int i) const {
      return (i < fErrors.size() ) ? fErrors[i] : 0;
   }

   /// name of the parameter
   std::string ParName(unsigned int i) const;

   /// set the Minos errors for parameter i (called by the Fitter class when running Minos)
   void SetMinosError(unsigned int i, double elow, double eup);

   /// query if parameter i has the Minos error
   bool HasMinosError(unsigned int i) const;

   /// lower Minos error. If Minos has not run for parameter i return the parabolic error
   double LowerError(unsigned int i) const;

   /// upper Minos error. If Minos has not run for parameter i return the parabolic error
   double UpperError(unsigned int i) const;

   /// parameter global correlation coefficient
   double GlobalCC(unsigned int i) const {
      return (i < fGlobalCC.size() ) ? fGlobalCC[i] : -1;
   }


   /// retrieve covariance matrix element
   double CovMatrix (unsigned int i, unsigned int j) const {
      if ( i >= fErrors.size() || j >= fErrors.size() ) return 0;
      if (fCovMatrix.size() == 0) return 0; // no matrix is available in case of non-valid fits
      if ( j < i )
         return fCovMatrix[j + i* (i+1) / 2];
      else
         return fCovMatrix[i + j* (j+1) / 2];
   }

   /// retrieve correlation elements
   double Correlation(unsigned int i, unsigned int j ) const {
      if ( i >= fErrors.size() || j >= fErrors.size() ) return 0;
      if (fCovMatrix.size() == 0) return 0; // no matrix is available in case of non-valid fits
      double tmp = CovMatrix(i,i)*CovMatrix(j,j);
      return ( tmp > 0) ? CovMatrix(i,j)/ std::sqrt(tmp) : 0;
   }

   /// fill covariance matrix elements using a generic matrix class implementing operator(i,j)
   /// the matrix must be previously allocates with right size (npar * npar)
   template<class Matrix>
   void GetCovarianceMatrix(Matrix & mat) const {
      unsigned int npar = fErrors.size();
      if (fCovMatrix.size() != npar*(npar+1)/2 ) return; // do nothing
      for (unsigned int i = 0; i< npar; ++i) {
         for (unsigned int j = 0; j<=i; ++j) {
            mat(i,j) = fCovMatrix[j + i*(i+1)/2 ];
            if (i != j) mat(j,i) = mat(i,j);
         }
      }
   }

   /// fill a correlation matrix elements using a generic symmetric matrix class implementing operator(i,j)
   /// the matrix must be previously allocates with right size (npar * npar)
   template<class Matrix>
   void GetCorrelationMatrix(Matrix & mat) const {
      unsigned int npar = fErrors.size();
      if (fCovMatrix.size() != npar*(npar+1)/2) return; // do nothing
      for (unsigned int i = 0; i< npar; ++i) {
         for (unsigned int j = 0; j<=i; ++j) {
            double tmp = fCovMatrix[i * (i +3)/2 ] * fCovMatrix[ j * (j+3)/2 ];
            mat(i,j) = (tmp > 0) ? fCovMatrix[j + i*(i+1)/2 ] / std::sqrt(tmp) : 0;
            if (i != j) mat(j,i) = mat(i,j);
         }
      }
   }

   /**
      get confidence intervals for an array of n points x.
      stride1 indicates the stride in the coordinate space while stride2 the stride in dimension space.
      For 1-dim points : stride1=1, stride2=1
      for multi-dim points arranged as (x0,x1,...,xN,y0,....yN)          stride1=1      stride2=n
      for multi-dim points arraged  as (x0,y0,..,x1,y1,...,xN,yN,..)     stride1=ndim,  stride2=1

      the confidence interval are returned in the array ci
      cl is the desired confidedence interval value
      norm is a flag to control if the intervals need to be normalized to the chi2/ndf value
      The intervals can be corrected optionally using the chi2/ndf value of the fit if a chi2 fit is performed.
      This has changed since ROOT 6.14, before the interval were corrected by default.
    */
   void GetConfidenceIntervals(unsigned int n, unsigned int stride1, unsigned int stride2, const double * x,  double * ci, double cl=0.95, bool norm = false ) const;

   /**
      evaluate confidence interval for the point specified in the passed data sets
      the confidence interval are returned in the array ci
      cl is the desired confidence interval value.
      This method is mantained for backward compatibility and will be deprecated
   */
   void GetConfidenceIntervals(const BinData & data, double * ci, double cl=0.95, bool norm = false ) const;

   /**
      evaluate confidence interval for the data set used in the last fit
      the confidence interval are returned as a vector of data points
    */
   std::vector<double> GetConfidenceIntervals(double cl=0.95, bool norm = false ) const;

   /**
      scan likelihood value of  parameter and fill the given graph.
    */
   bool Scan(unsigned int ipar, unsigned int &npoints, double *pntsx, double *pntsy, double xmin = 0, double xmax = 0);

   /**
      create contour of two parameters around the minimum
      pass as option confidence level:  default is a value of 0.683
   */
   bool Contour(unsigned int ipar, unsigned int jpar, unsigned int &npoints, double *pntsx, double *pntsy, double confLevel = 0.683);

   /// get index for parameter name (return -1 if not found)
   int Index(const std::string & name) const;

   ///normalize errors using chi2/ndf for chi2 fits
   void NormalizeErrors();

   /// flag to chek if errors are normalized
   bool NormalizedErrors() const { return fNormalized; }

   /// print the result and optionaly covariance matrix and correlations
   void Print(std::ostream & os, bool covmat = false) const;

   ///print error matrix and correlations
   void PrintCovMatrix(std::ostream & os) const;

   /// query if a parameter is bound
   bool IsParameterBound(unsigned int ipar) const;

   /// query if a parameter is fixed
   bool IsParameterFixed(unsigned int ipar) const;

   /// retrieve parameter bounds - return false if parameter is not bound
   bool ParameterBounds(unsigned int ipar, double &lower, double &upper) const;


   /// get name of parameter (deprecated)
   std::string GetParameterName(unsigned int ipar) const {
      return ParName(ipar);
   }


protected:


   /// Return pointer non const pointer to model (fit) function with fitted parameter values.
   /// used by Fitter class
   std::shared_ptr<IModelFunction> ModelFunction()  { return fFitFunc; }
   void SetModelFunction(const std::shared_ptr<IModelFunction> & func) { fFitFunc = func; }

   friend class Fitter;


   bool fValid;             // flag for indicating valid fit
   bool fNormalized;        // flag for indicating is errors are normalized
   unsigned int fNFree;     // number of fit free parameters (total parameters are in size of parameter vector)
   unsigned int fNdf;       // number of degree of freedom
   unsigned int fNCalls;    // number of function calls
   int fStatus;             // minimizer status code
   int fCovStatus;          // covariance matrix status code
   double fVal;             // minimum function value
   double fEdm;             // expected distance from mimimum
   double fChi2;            // fit chi2 value (different than fval in case of chi2 fits)
   std::shared_ptr<ROOT::Math::Minimizer> fMinimizer; //! minimizer object used for fitting
   std::shared_ptr<ROOT::Math::IMultiGenFunction> fObjFunc; //! objective function used for fitting
   std::shared_ptr<IModelFunction> fFitFunc; //! model function resulting  from the fit.
   std::shared_ptr<FitData>    fFitData; //! data set used in the fit
   std::map<unsigned int, bool>           fFixedParams; // list of fixed parameters
   std::map<unsigned int, unsigned int>   fBoundParams; // list of limited parameters
   std::vector<std::pair<double,double> >  fParamBounds; // parameter bounds
   std::vector<double>         fParams;  // parameter values. Size is total number of parameters
   std::vector<double>         fErrors;  // errors
   std::vector<double>         fCovMatrix;  // covariance matrix (size is npar*(npar+1)/2) where npar is total parameters
   std::vector<double>         fGlobalCC;   // global Correlation coefficient
   std::map<unsigned int, std::pair<double,double> > fMinosErrors;   // map contains the two Minos errors
   std::string fMinimType;              // string indicating type of minimizer
   std::vector<std::string> fParNames;  // parameter names (only with FCN only fits, when fFitFunc=0)

};


   } // end namespace Fit

} // end namespace ROOT





#endif /* ROOT_Fit_FitResult */
