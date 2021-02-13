// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FumiliFCNBase
#define ROOT_Minuit2_FumiliFCNBase

#include "Minuit2/FCNBase.h"
#include <cassert>
#include <vector>

namespace ROOT {

namespace Minuit2 {

//____________________________________________________________________________________________
/**

Extension of the FCNBase for the Fumili method. Fumili applies only to
minimization problems used for fitting. The method is based on a
linearization of the model function negleting second derivatives.
User needs to provide the model function. The figure-of-merit describing
the difference between the model function and the actual measurements
has to be implemented by the user in a subclass of FumiliFCNBase.
For an example see the FumiliChi2FCN and FumiliStandardChi2FCN classes.


@author  Andras Zsenei and Lorenzo Moneta, Creation date: 23 Aug 2004

@see <A HREF="http://www.cern.ch/winkler/minuit/tutorial/mntutorial.pdf">MINUIT Tutorial</A> on function minimization,
section 5

@see FumiliChi2FCN

@see FumiliStandardChi2FCN

@ingroup Minuit

 */

class FumiliFCNBase : public FCNBase {

public:
   /**
      Default Constructor. Need in this case to create when implementing EvaluateAll the Gradient and Hessian vectors
      with the right size
   */

   FumiliFCNBase() : fNumberOfParameters(0), fValue(0) {}

   /**

      Constructor which initializes the class with the function provided by the
      user for modeling the data.

      @param npar the number of parameters

   */

   FumiliFCNBase(unsigned int npar)
      : fNumberOfParameters(npar), fValue(0), fGradient(std::vector<double>(npar)),
        fHessian(std::vector<double>(static_cast<int>(0.5 * npar * (npar + 1))))
   {
   }

   //   FumiliFCNBase(const ParametricFunction& modelFCN) { fModelFunction = &modelFCN; }

   virtual ~FumiliFCNBase() {}

   /**

      Evaluate function Value, Gradient and Hessian using Fumili approximation, for values of parameters p
      The resul is cached inside and is return from the FumiliFCNBase::Value ,  FumiliFCNBase::Gradient and
      FumiliFCNBase::Hessian methods

      @param par vector of parameters

   **/

   virtual void EvaluateAll(const std::vector<double> &par) = 0;

   /**
      Return cached Value of objective function estimated previously using the  FumiliFCNBase::EvaluateAll method

   **/

   virtual double Value() const { return fValue; }

   /**
      Return cached Value of function Gradient estimated previously using the  FumiliFCNBase::EvaluateAll method
   **/

   virtual const std::vector<double> &Gradient() const { return fGradient; }

   /**
      Return Value of the i-th j-th element of the Hessian matrix estimated previously using the
   FumiliFCNBase::EvaluateAll method
      @param row row Index of the matrix
      @param col col Index of the matrix
   **/

   virtual double Hessian(unsigned int row, unsigned int col) const
   {
      assert(row < fGradient.size() && col < fGradient.size());
      if (row > col)
         return fHessian[col + row * (row + 1) / 2];
      else
         return fHessian[row + col * (col + 1) / 2];
   }

   /**
      return number of function variable (parameters) , i.e. function dimension
   */

   virtual unsigned int Dimension() { return fNumberOfParameters; }

protected:
   /**
      initialize and reset values of gradien and Hessian
   */

   virtual void InitAndReset(unsigned int npar)
   {
      fNumberOfParameters = npar;
      fGradient = std::vector<double>(npar);
      fHessian = std::vector<double>(static_cast<int>(0.5 * npar * (npar + 1)));
   }

   // methods to be used by the derived classes to set the values
   void SetFCNValue(double value) { fValue = value; }

   std::vector<double> &Gradient() { return fGradient; }

   std::vector<double> &Hessian() { return fHessian; }

private:
   unsigned int fNumberOfParameters;
   double fValue;
   std::vector<double> fGradient;
   std::vector<double> fHessian;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FumiliFCNBase
