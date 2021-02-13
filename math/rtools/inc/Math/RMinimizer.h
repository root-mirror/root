// Author: K. Hermansen and L. Moneta, Aug 2014 

// Implementation file for class RMinimizer

#ifndef ROOT_Math_RMinimizer
#define ROOT_Math_RMinimizer

#include "Math/Functor.h"

#include "Math/IParamFunctionfwd.h"

#include "Math/BasicMinimizer.h"

#include "TMatrixD.h"

#include <vector>
#include <string>

namespace ROOT {
   namespace Math{	

      /*! \brief RMinimizer class.
       *
       *    Minimizer class that uses the ROOT/R interface to pass functions and minimize them in R.
       *    
       *    The class implements the ROOT::Math::Minimizer interface and can be instantiated using the 
       *    ROOT plugin manager (plugin name is "RMinimizer"). The various minimization algorithms 
       *    (BFGS, Nelder-Mead, SANN, etc..) can be passed as an option. 
       *    The default algorithm is BFGS.
       *
       *    The library for this and future R/ROOT classes is currently libRtools.so
       */
      class   RMinimizer  :   public  ROOT::Math::BasicMinimizer    {
         protected:
            std::string fMethod; /*!< minimizer method to be used, must be of a type listed in R optim or optimx descriptions */
         
         private:
         std::vector<double>   fErrors; /*!< vector of parameter errors */
            TMatrixD        fCovMatrix; /*!< covariant matrix */
            TMatrixD       fHessMatrix; /*!< Hessian matrix */
         
         public:
            /*! \brief Default constructor
             *
             * Default constructor with option for the method of minimization, can be any of the following:
            *"Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent" (Brent only for 1D minimization)
            *
            *See R optim or optimx descriptions for more details and options.
            *
            */	
            RMinimizer(Option_t *method);
            ///Destructor
            virtual ~RMinimizer() {}
            ///Function to find the minimum
            virtual bool Minimize();
            ///Returns the number of function calls
            virtual unsigned int NCalls() const;
            ///Returns the ith jth component of the Hessian matrix
            double HessMatrix(unsigned int i, unsigned int j) const;
            /// minimizer provides error and error matrix
            virtual bool ProvidesError() const { return !(fErrors.empty()); }
            /// return errors at the minimum
            virtual const double * Errors() const { return fErrors.data(); }
            /** return covariance matrices element for variables ivar,jvar
            if the variable is fixed the return value is zero
            The ordering of the variables is the same as in the parameter and errors vectors
            */
           virtual double CovMatrix(unsigned int  ivar , unsigned int jvar ) const {
              return fCovMatrix(ivar, jvar);
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
               int ndim = NDim(); 
               if (fCovMatrix.GetNrows() != ndim || fCovMatrix.GetNcols() != ndim ) return false; 
               std::copy(fCovMatrix.GetMatrixArray(), fCovMatrix.GetMatrixArray() + ndim*ndim, covMat);
               return true;
            }
      };

   }
}
#endif
