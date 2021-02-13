// @(#)root/minuit:$Id$
// Author: L. Moneta Wed Oct 25 16:28:55 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for class TMinuitMinimizer

#include "TMinuitMinimizer.h"
#include "Math/IFunction.h"
#include "Fit/ParameterSettings.h"

#include "TMinuit.h"
#include "TROOT.h"
#include "TList.h"

#include "TGraph.h" // needed for scan
#include "TError.h"

#include "TMatrixDSym.h" // needed for inverting the matrix

#include "ThreadLocalStorage.h"

#include <cassert>
#include <algorithm>
#include <functional>
#include <cmath>

//______________________________________________________________________________
//
//  TMinuitMinimizer class implementing the ROOT::Math::Minimizer interface using
//  TMinuit.
//  This class is normally instantiates using the plug-in manager
//  (plug-in with name Minuit or TMinuit)
//  In addition the user can choose the minimizer algorithm: Migrad (the default one), Simplex, or Minimize (combined Migrad + Simplex)
//
//__________________________________________________________________________________________

// initialize the static instances

// Implement a thread local static member
static ROOT::Math::IMultiGenFunction *&GetGlobalFuncPtr() {
   TTHREAD_TLS(ROOT::Math::IMultiGenFunction *) fgFunc = nullptr;
   return fgFunc;
}
TMinuit * TMinuitMinimizer::fgMinuit = 0;
bool TMinuitMinimizer::fgUsed = false;
bool TMinuitMinimizer::fgUseStaticMinuit = true;   // default case use static Minuit instance

ClassImp(TMinuitMinimizer);


TMinuitMinimizer::TMinuitMinimizer(ROOT::Minuit::EMinimizerType type, unsigned int ndim ) :
   fUsed(false),
   fMinosRun(false),
   fDim(ndim),
   fType(type),
   fMinuit(0)
{
   // Constructor for TMinuitMinimier class via an enumeration specifying the minimization
   // algorithm type. Supported types are : kMigrad, kSimplex, kCombined (a combined
   // Migrad + Simplex minimization) and kMigradImproved (a Migrad mininimization folloed by an
   // improved search for global minima). The default type is Migrad (kMigrad).

   // initialize if npar is given
   if (fDim > 0) InitTMinuit(fDim);
}

TMinuitMinimizer::TMinuitMinimizer(const char *  type, unsigned int ndim ) :
   fUsed(false),
   fMinosRun(false),
   fDim(ndim),
   fMinuit(0)
{
   // constructor from a char * for the algorithm type, used by the plug-in manager
   // The names supported (case unsensitive) are:
   //  Migrad (default), Simplex, Minimize (for the combined Migrad+ Simplex) and Migrad_imp

   // select type from the string
   std::string algoname(type);
   std::transform(algoname.begin(), algoname.end(), algoname.begin(), (int(*)(int)) tolower );

   ROOT::Minuit::EMinimizerType algoType = ROOT::Minuit::kMigrad;
   if (algoname == "simplex")   algoType = ROOT::Minuit::kSimplex;
   if (algoname == "minimize" ) algoType = ROOT::Minuit::kCombined;
   if (algoname == "migradimproved" ) algoType = ROOT::Minuit::kMigradImproved;
   if (algoname == "scan" )           algoType = ROOT::Minuit::kScan;
   if (algoname == "seek" )           algoType = ROOT::Minuit::kSeek;

   fType = algoType;

   // initialize if npar is given
   if (fDim > 0) InitTMinuit(fDim);

}

TMinuitMinimizer::~TMinuitMinimizer()
{
   // Destructor implementation.
   if (fMinuit && !fgUseStaticMinuit) {
      delete fMinuit;
      fgMinuit = 0;
   }
}

TMinuitMinimizer::TMinuitMinimizer(const TMinuitMinimizer &) :
   Minimizer()
{
   // Implementation of copy constructor (it is private).
}

TMinuitMinimizer & TMinuitMinimizer::operator = (const TMinuitMinimizer &rhs)
{
   // Implementation of assignment operator (private)
   if (this == &rhs) return *this;  // time saving self-test
   return *this;
}

bool TMinuitMinimizer::UseStaticMinuit(bool on ) {
   // static method to control usage of global TMinuit instance
   bool prev = fgUseStaticMinuit;
   fgUseStaticMinuit = on;
   return prev;
}

void TMinuitMinimizer::InitTMinuit(int dim) {

   // when called a second time check dimension - create only if needed
   // initialize the minuit instance - recreating a new one if needed
   if (fMinuit ==0 ||  dim > fMinuit->fMaxpar) {

      // case not using the global instance - recreate it all the time
      if (fgUseStaticMinuit) {

         // re-use gMinuit as static instance of TMinuit
         // which can be accessed by the user after minimization
         // check if fgMinuit is different than gMinuit
         // case 1: fgMinuit not zero but fgMinuit has been deleted (not in gROOT): set to zero
         // case 2: fgMinuit not zero and exists in global list  : set fgMinuit to gMinuit
         // case 3: fgMinuit zero - and gMinuit not zero: create a new instance locally to avoid conflict
         if (fgMinuit != gMinuit) {
            // if object exists in gROOT remove it to avoid a memory leak
            if (fgMinuit ) {
               if (gROOT->GetListOfSpecials()->FindObject(fgMinuit) == 0) {
                  // case 1: object does not exists in gROOT - means it has been deleted
                  fgMinuit = 0;
               }
               else {
                  // case 2: object exists - but gMinuit points to something else
                  // restore gMinuit to the one used before by TMinuitMinimizer
                  gMinuit = fgMinuit;
               }
            }
            else {
               // case 3: avoid reusing existing one - mantain fgMinuit to zero
               // otherwise we will get a double delete if user deletes externally gMinuit
               // in this case we will loose gMinuit instance
//                fgMinuit = gMinuit;
//                fgUsed = true;  // need to reset in case  other gMinuit instance is later used
            }
         }

         // check if need to create a new TMinuit instance
         if (fgMinuit == 0) {
            fgUsed = false;
            fgMinuit =  new TMinuit(dim);
         }
         else if (fgMinuit->GetNumPars() != int(dim) ) {
            delete fgMinuit;
            fgUsed = false;
            fgMinuit =  new TMinuit(dim);
         }

         fMinuit = fgMinuit;
      }

      else {
         // re- create all the time a new instance of TMinuit (fgUseStaticMinuit is false)
         if (fMinuit) delete fMinuit;
         fMinuit =  new TMinuit(dim);
         fgMinuit = fMinuit;
         fgUsed = false;
      }

   }  // endif fMinuit ==0 || dim > fMaxpar

   fDim = dim;

   R__ASSERT(fMinuit);

   // set print level in TMinuit
   double arglist[1];
   int ierr= 0;
   // TMinuit level is shift by 1 -1 means 0;
   arglist[0] = PrintLevel() - 1;
   fMinuit->mnexcm("SET PRINT",arglist,1,ierr);
   if (PrintLevel() <= 0) SuppressMinuitWarnings();
}


void TMinuitMinimizer::SetFunction(const  ROOT::Math::IMultiGenFunction & func) {
   // Set the objective function to be minimized, by passing a function object implement the
   // basic multi-dim Function interface. In this case the derivatives will be
   // calculated by Minuit
   // Here a TMinuit instance is created since only at this point we know the number of parameters


   fDim = func.NDim();

   // create TMinuit if needed
   InitTMinuit(fDim);

   // assign to the static pointer (NO Thread safety here)
   GetGlobalFuncPtr() = const_cast<ROOT::Math::IMultiGenFunction *>(&func);
   fMinuit->SetFCN(&TMinuitMinimizer::Fcn);

   // switch off gradient calculations
   double arglist[1];
   int ierr = 0;
   fMinuit->mnexcm("SET NOGrad",arglist,0,ierr);
}

void TMinuitMinimizer::SetFunction(const  ROOT::Math::IMultiGradFunction & func) {
   // Set the objective function to be minimized, by passing a function object implement the
   // multi-dim gradient Function interface. In this case the function derivatives are provided
   // by the user via this interface and there not calculated by Minuit.

   fDim = func.NDim();

   // create TMinuit if needed
   InitTMinuit(fDim);

   // assign to the static pointer (NO Thread safety here)
   GetGlobalFuncPtr() = const_cast<ROOT::Math::IMultiGradFunction *>(&func);
   fMinuit->SetFCN(&TMinuitMinimizer::FcnGrad);

   // set gradient
   // by default do not check gradient calculation
   // it cannot be done here, check can be done only after having defined the parameters
   double arglist[1];
   int ierr = 0;
   arglist[0] = 1;
   fMinuit->mnexcm("SET GRAD",arglist,1,ierr);
}

void TMinuitMinimizer::Fcn( int &, double * , double & f, double * x , int /* iflag */) {
   // implementation of FCN static function used internally by TMinuit.
   // Adapt IMultiGenFunction interface to TMinuit FCN static function
   f = GetGlobalFuncPtr()->operator()(x);
}

void TMinuitMinimizer::FcnGrad( int &, double * g, double & f, double * x , int iflag ) {
   // implementation of FCN static function used internally by TMinuit.
   // Adapt IMultiGradFunction interface to TMinuit FCN static function in the case of user
   // provided gradient.
   ROOT::Math::IMultiGradFunction * gFunc = dynamic_cast<ROOT::Math::IMultiGradFunction *> ( GetGlobalFuncPtr());

   assert(gFunc != 0);
   f = gFunc->operator()(x);

   // calculates also derivatives
   if (iflag == 2) gFunc->Gradient(x,g);
}

bool TMinuitMinimizer::SetVariable(unsigned int ivar, const std::string & name, double val, double step) {
   // set a free variable.
   if (!CheckMinuitInstance()) return false;

   if (fgUseStaticMinuit) fUsed = fgUsed;

   // clear after minimization when setting params
   if (fUsed) DoClear();

   // check if parameter was defined and in case it was fixed, release it
   DoReleaseFixParameter(ivar);

   int iret = fMinuit->DefineParameter(ivar , name.c_str(), val, step, 0., 0. );
   return (iret == 0);
}

bool TMinuitMinimizer::SetLimitedVariable(unsigned int ivar, const std::string & name, double val, double step, double lower, double upper) {
   // set a limited variable.
   if (!CheckMinuitInstance()) return false;

   if (fgUseStaticMinuit) fUsed = fgUsed;

   // clear after minimization when setting params
   if (fUsed) DoClear();

   // check if parameter was defined and in case it was fixed, release it
   DoReleaseFixParameter(ivar);

   int iret = fMinuit->DefineParameter(ivar, name.c_str(), val, step, lower, upper );
   return (iret == 0);
}

bool TMinuitMinimizer::SetLowerLimitedVariable(unsigned int  ivar , const std::string & name , double val , double step , double lower ) {
   // set a lower limited variable
   // since is not supported in TMinuit , just use a artificial large value
   Warning("TMinuitMinimizer::SetLowerLimitedVariable","not implemented - use as upper limit 1.E7 instead of +inf");
   return SetLimitedVariable(ivar, name, val , step, lower, lower+ 1.E7);  // use 1.E7 which will make TMinuit happy
}

bool TMinuitMinimizer::SetUpperLimitedVariable(unsigned int  ivar , const std::string & name , double val , double step , double upper ) {
   // set a upper limited variable
   // since is not supported in TMinuit , just use a artificial large negative value
   Warning("TMinuitMinimizer::SetUpperLimitedVariable","not implemented - use as lower limit -1.E7 instead of -inf");
   return SetLimitedVariable(ivar, name, val , step, upper -1.E7, upper);
}


bool TMinuitMinimizer::CheckMinuitInstance() const {
   // check instance of fMinuit
   if (fMinuit == 0) {
      Error("TMinuitMinimizer::CheckMinuitInstance","Invalid TMinuit pointer. Need to call first SetFunction");
      return false;
   }
   return true;
}

bool TMinuitMinimizer::CheckVarIndex(unsigned int ivar) const {
   // check index of Variable (assume fMinuit exists)
   if ((int) ivar >= fMinuit->fNu ) {
      Error("TMinuitMinimizer::CheckVarIndex","Invalid parameter index");
      return false;
   }
   return true;
}


bool TMinuitMinimizer::SetFixedVariable(unsigned int ivar, const std::string & name, double val) {
   // set a fixed variable.
   if (!CheckMinuitInstance()) return false;

   // clear after minimization when setting params
   if (fgUseStaticMinuit) fUsed = fgUsed;

   // clear after minimization when setting params
   if (fUsed) DoClear();

   // put an arbitrary step (0.1*abs(value) otherwise TMinuit consider the parameter as constant
   // constant parameters are treated differently (they are ignored inside TMinuit and not considered in the
   // total list of parameters)
   double step = ( val != 0) ? 0.1 * std::abs(val) : 0.1;
   int iret = fMinuit->DefineParameter(ivar, name.c_str(), val, step, 0., 0. );
   if (iret == 0) iret = fMinuit->FixParameter(ivar);
   return (iret == 0);
}

bool TMinuitMinimizer::SetVariableValue(unsigned int ivar, double val) {
   // set the value of an existing variable
   // parameter must exist or return false
   if (!CheckMinuitInstance()) return false;

   double arglist[2];
   int ierr = 0;

   arglist[0] = ivar+1;  // TMinuit starts from 1
   arglist[1] = val;
   fMinuit->mnexcm("SET PAR",arglist,2,ierr);
   return (ierr==0);
}

bool TMinuitMinimizer::SetVariableStepSize(unsigned int ivar, double step) {
   // set the step-size of an existing variable
   // parameter must exist or return false
   if (!CheckMinuitInstance()) return false;
   // need to re-implement re-calling mnparm
   // get first current parameter values and limits
   double curval,err, lowlim, uplim;
   int iuint;  // internal index
   TString name;
   fMinuit->mnpout(ivar, name, curval, err, lowlim, uplim,iuint);
   if (iuint == -1) return false;
   int iret = fMinuit->DefineParameter(ivar, name, curval, step, lowlim, uplim );
   return (iret == 0);

}

bool TMinuitMinimizer::SetVariableLowerLimit(unsigned int ivar, double lower ) {
   // set the limits of an existing variable
   // parameter must exist or return false
   Warning("TMinuitMinimizer::SetVariableLowerLimit","not implemented - use as upper limit 1.E30 instead of +inf");
   return SetVariableLimits(ivar, lower, 1.E30);
}
bool TMinuitMinimizer::SetVariableUpperLimit(unsigned int ivar, double upper ) {
   // set the limits of an existing variable
   // parameter must exist or return false
   Warning("TMinuitMinimizer::SetVariableUpperLimit","not implemented - - use as lower limit -1.E30 instead of +inf");
   return SetVariableLimits(ivar, -1.E30, upper);
}

bool TMinuitMinimizer::SetVariableLimits(unsigned int ivar, double lower, double upper) {
   // set the limits of an existing variable
   // parameter must exist or return false

   if (!CheckMinuitInstance()) return false;
   // need to re-implement re-calling mnparm
   // get first current parameter values and limits
   double curval,err, lowlim, uplim;
   int iuint;  // internal index
   TString name;
   fMinuit->mnpout(ivar, name, curval, err, lowlim, uplim,iuint);
   if (iuint == -1) return false;
   int iret = fMinuit->DefineParameter(ivar, name, curval, err, lower, upper );
   return (iret == 0);

}

bool TMinuitMinimizer::FixVariable(unsigned int ivar) {
   // Fix an existing variable
   if (!CheckMinuitInstance()) return false;
   if (!CheckVarIndex(ivar)) return false;
   int iret = fMinuit->FixParameter(ivar);
   return (iret == 0);
}

bool TMinuitMinimizer::ReleaseVariable(unsigned int ivar) {
   // Fix an existing variable
   if (!CheckMinuitInstance()) return false;
   if (!CheckVarIndex(ivar)) return false;
   int iret = fMinuit->Release(ivar);
   return (iret == 0);
}

bool TMinuitMinimizer::IsFixedVariable(unsigned int ivar) const {
   // query if variable is fixed
   if (!CheckMinuitInstance()) return false;
   if (!CheckVarIndex(ivar)) return false;
   return (fMinuit->fNiofex[ivar] == 0 );
}

bool TMinuitMinimizer::GetVariableSettings(unsigned int ivar, ROOT::Fit::ParameterSettings & var) const {
   // retrieve variable settings (all set info on the variable)
   if (!CheckMinuitInstance()) return false;
   if (!CheckVarIndex(ivar)) return false;
   double curval,err, lowlim, uplim;
   int iuint;  // internal index
   TString name;
   fMinuit->mnpout(ivar, name, curval, err, lowlim, uplim,iuint);
   if (iuint == -1) return false;
   var.Set(name.Data(), curval, err, lowlim, uplim);
   if (IsFixedVariable(ivar)) var.Fix();
   return true;
}



std::string TMinuitMinimizer::VariableName(unsigned int ivar) const {
   // return the variable name
   if (!CheckMinuitInstance()) return std::string();
   if (!CheckVarIndex(ivar)) return std::string();
   return fMinuit->fCpnam[ivar].Data();
}

int TMinuitMinimizer::VariableIndex(const std::string & ) const {
   // return variable index
   Error("TMinuitMinimizer::VariableIndex"," find index of a variable from its name  is not implemented in TMinuit");
   return -1;
}

bool TMinuitMinimizer::Minimize() {
   // perform the minimization using the algorithm chosen previously by the user
   // By default Migrad is used.
   // Return true if the found minimum is valid and update internal chached values of
   // minimum values, errors and covariance matrix.
   // Status of minimizer is set to:
   // migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult


   if (fMinuit == 0) {
      Error("TMinuitMinimizer::Minimize","invalid TMinuit pointer. Need to call first SetFunction and SetVariable");
      return false;
   }


   // total number of parameter defined in Minuit is fNu
   if (fMinuit->fNu <  static_cast<int>(fDim) ) {
      Error("TMinuitMinimizer::Minimize","The total number of defined parameters is different than the function dimension, npar = %d, dim = %d",fMinuit->fNu, fDim);
      return false;
   }

   int printlevel = PrintLevel();

   // total number of free parameter is 0
   if (fMinuit->fNpar <= 0) {
      // retrieve parameters values  from TMinuit
      RetrieveParams();
      fMinuit->fAmin = (*GetGlobalFuncPtr())(&fParams.front());
      if (printlevel > 0) Info("TMinuitMinimizer::Minimize","There are no free parameter - just compute the function value");
      return true;
   }


   double arglist[10];
   int ierr = 0;


   // set error and print level
   arglist[0] = ErrorDef();
   fMinuit->mnexcm("SET Err",arglist,1,ierr);

   arglist[0] = printlevel - 1;
   fMinuit->mnexcm("SET PRINT",arglist,1,ierr);

   // suppress warning in case Printlevel() == 0
   if (printlevel == 0)    fMinuit->mnexcm("SET NOW",arglist,0,ierr);

   // set precision if needed
   if (Precision() > 0)  {
      arglist[0] = Precision();
      fMinuit->mnexcm("SET EPS",arglist,1,ierr);
   }

   // set strategy
   int strategy = Strategy();
   if (strategy >=0 && strategy <=2 ) {
      arglist[0] = strategy;
      fMinuit->mnexcm("SET STR",arglist,1,ierr);
   }

   arglist[0] = MaxFunctionCalls();
   arglist[1] = Tolerance();

   int nargs = 2;

   switch (fType){
   case ROOT::Minuit::kMigrad:
      // case of Migrad
      fMinuit->mnexcm("MIGRAD",arglist,nargs,ierr);
      break;
   case ROOT::Minuit::kCombined:
      // case of combined (Migrad+ simplex)
      fMinuit->mnexcm("MINIMIZE",arglist,nargs,ierr);
      break;
   case ROOT::Minuit::kSimplex:
      // case of Simlex
      fMinuit->mnexcm("SIMPLEX",arglist,nargs,ierr);
      break;
   case ROOT::Minuit::kScan:
      // case of Scan (scan all parameters with default values)
      nargs = 0;
      fMinuit->mnexcm("SCAN",arglist,nargs,ierr);
      break;
   case ROOT::Minuit::kSeek:
      // case of Seek (random find minimum in a hypercube around current parameter values
      // use Tolerance as measures for standard deviation (if < 1) used default value in Minuit ( supposed to be  3)
      nargs = 1;
      if (arglist[1] >= 1.) nargs = 2;
      fMinuit->mnexcm("SEEK",arglist,nargs,ierr);
      break;
   default:
      // default: use Migrad
      fMinuit->mnexcm("MIGRAD",arglist,nargs,ierr);

   }

   fgUsed = true;
   fUsed = true;

   fStatus = ierr;
   int minErrStatus = ierr;

   if (printlevel>2) Info("TMinuitMinimizer::Minimize","Finished to run MIGRAD - status %d",ierr);

   // run improved if needed
   if (ierr == 0 && fType == ROOT::Minuit::kMigradImproved) {
      fMinuit->mnexcm("IMPROVE",arglist,1,ierr);
      fStatus += 1000*ierr;
      if (printlevel>2) Info("TMinuitMinimizer::Minimize","Finished to run IMPROVE - status %d",ierr);
   }


   // check if Hesse needs to be run
   // Migrad runs inside it automatically for strategy >=1. Do also
   // in case improve or other minimizers are used
   if (minErrStatus == 0  && (IsValidError() || ( strategy >=1 && CovMatrixStatus() < 3) ) ) {
      fMinuit->mnexcm("HESSE",arglist,1,ierr);
      fStatus += 100*ierr;
      if (printlevel>2) Info("TMinuitMinimizer::Minimize","Finished to run HESSE - status %d",ierr);
   }

   // retrieve parameters and errors  from TMinuit
   RetrieveParams();

   if (minErrStatus == 0) {

      // store global min results (only if minimization is OK)
      // ignore cases when Hesse or IMprove return error different than zero
      RetrieveErrorMatrix();

      // need to re-run Minos again if requested
      fMinosRun = false;

      return true;

   }
   return false;

}

void TMinuitMinimizer::RetrieveParams() {
   // retrieve from TMinuit minimum parameter values
   // and errors

   assert(fMinuit != 0);

   // get parameter values
   if (fParams.size() != fDim) fParams.resize( fDim);
   if (fErrors.size() != fDim) fErrors.resize( fDim);
   for (unsigned int i = 0; i < fDim; ++i) {
      fMinuit->GetParameter( i, fParams[i], fErrors[i]);
   }
}

void TMinuitMinimizer::RetrieveErrorMatrix() {
   // get covariance error matrix from TMinuit
   // when some parameters are fixed filled the corresponding rows and column with zero's

   assert(fMinuit != 0);

   unsigned int nfree = NFree();

   unsigned int ndim2 = fDim*fDim;
   if (fCovar.size() != ndim2 )  fCovar.resize(fDim*fDim);
   if (nfree >= fDim) { // no fixed parameters
      fMinuit->mnemat(&fCovar.front(), fDim);
   }
   else {
      // case of fixed params need to take care
      std::vector<double> tmpMat(nfree*nfree);
      fMinuit->mnemat(&tmpMat.front(), nfree);

      unsigned int l = 0;
      for (unsigned int i = 0; i < fDim; ++i) {

         if ( fMinuit->fNiofex[i] > 0 ) {  // not fixed ?
            unsigned int m = 0;
            for (unsigned int j = 0; j <= i; ++j) {
               if ( fMinuit->fNiofex[j] > 0 ) {  //not fixed
                  fCovar[i*fDim + j] = tmpMat[l*nfree + m];
                  fCovar[j*fDim + i] = fCovar[i*fDim + j];
                  m++;
               }
            }
            l++;
         }
      }

   }
}

unsigned int TMinuitMinimizer::NCalls() const {
   // return total number of function calls
   if (fMinuit == 0) return 0;
   return fMinuit->fNfcn;
}

double TMinuitMinimizer::MinValue() const {
   // return minimum function value

   // use part of code from mnstat
   if (!fMinuit) return 0;
   double minval = fMinuit->fAmin;
   if (minval == fMinuit->fUndefi) return 0;
   return minval;
}

double TMinuitMinimizer::Edm() const {
   // return expected distance from the minimum

   // use part of code from mnstat
   if (!fMinuit) return -1;
   if (fMinuit->fAmin == fMinuit->fUndefi || fMinuit->fEDM == fMinuit->fBigedm) return fMinuit->fUp;
   return fMinuit->fEDM;
}

unsigned int TMinuitMinimizer::NFree() const {
    // return number of free parameters
   if (!fMinuit) return 0;
   if (fMinuit->fNpar < 0) return 0;
   return fMinuit->fNpar;
}

bool TMinuitMinimizer::GetCovMatrix(double * cov) const {
   // get covariance matrix
   int covStatus = CovMatrixStatus();
   if ( fCovar.size() != fDim*fDim || covStatus < 2) {
      Error("TMinuitMinimizer::GetHessianMatrix","Hessian matrix has not been computed - status %d",covStatus);
      return false;
   }
   std::copy(fCovar.begin(), fCovar.end(), cov);
   TMatrixDSym cmat(fDim,cov);
   return true;
}

bool TMinuitMinimizer::GetHessianMatrix(double * hes) const {
   // get Hessian - inverse of covariance matrix
   // just invert it
   // but need to get the compact form to avoid the zero for the fixed parameters
   int covStatus = CovMatrixStatus();
   if ( fCovar.size() != fDim*fDim || covStatus < 2) {
      Error("TMinuitMinimizer::GetHessianMatrix","Hessian matrix has not been computed - status %d",covStatus);
      return false;
   }
   // case of fixed params need to take care
   unsigned int nfree = NFree();
   TMatrixDSym mat(nfree);
   fMinuit->mnemat(mat.GetMatrixArray(), nfree);
   // invert the matrix
   // probably need to check if failed. In that case inverse is equal to original
   mat.Invert();

   unsigned int l = 0;
   for (unsigned int i = 0; i < fDim; ++i) {
      if ( fMinuit->fNiofex[i] > 0 ) {  // not fixed ?
         unsigned int m = 0;
         for (unsigned int j = 0; j <= i; ++j) {
            if ( fMinuit->fNiofex[j] > 0 ) {  //not fixed
               hes[i*fDim + j] =  mat(l,m);
               hes[j*fDim + i] = hes[i*fDim + j];
               m++;
            }
         }
         l++;
      }
   }
   return true;
}
//    if  ( fCovar.size() != fDim*fDim ) return false;
//    TMatrixDSym mat(fDim, &fCovar.front() );
//    std::copy(mat.GetMatrixArray(), mat.GetMatrixArray()+ mat.GetNoElements(), hes);
//    return true;
// }

int TMinuitMinimizer::CovMatrixStatus() const {
   // return status of covariance matrix
   //           status:  0= not calculated at all
   //                    1= approximation only, not accurate
   //                    2= full matrix, but forced positive-definite
   //                    3= full accurate covariance matrix

   // use part of code from mnstat
   if (!fMinuit) return 0;
   if (fMinuit->fAmin == fMinuit->fUndefi) return 0;
   return fMinuit->fISW[1];
}

double TMinuitMinimizer::GlobalCC(unsigned int i) const {
   // global correlation coefficient for parameter i
   if (!fMinuit) return 0;
   if (!fMinuit->fGlobcc) return 0;
   if (int(i) >= fMinuit->fNu) return 0;
   // get internal number in Minuit
   int iin = fMinuit->fNiofex[i];
   // index in TMinuit starts from 1
   if (iin < 1) return 0;
   return fMinuit->fGlobcc[iin-1];
}

bool TMinuitMinimizer::GetMinosError(unsigned int i, double & errLow, double & errUp, int ) {
   // Perform Minos analysis for the given parameter  i

   if (fMinuit == 0) {
      Error("TMinuitMinimizer::GetMinosError","invalid TMinuit pointer. Need to call first SetFunction and SetVariable");
      return false;
   }

   // check if parameter is fixed
   if (fMinuit->fNiofex[i] == 0 ) {
      if (PrintLevel() > 0) Info("TMinuitMinimizer::GetMinosError","Parameter %s is fixed. There are no Minos error to calculate. Ignored.",VariableName(i).c_str());
      errLow = 0; errUp = 0;
      return true;
   }

   double arglist[2];
   int ierr = 0;


   // set error, print level, precision and strategy if they have changed
   if (fMinuit->fUp != ErrorDef() ) {
      arglist[0] = ErrorDef();
      fMinuit->mnexcm("SET Err",arglist,1,ierr);
   }

   if (fMinuit->fISW[4] != (PrintLevel()-1) ) {
      arglist[0] = PrintLevel()-1;
      fMinuit->mnexcm("SET PRINT",arglist,1,ierr);
      // suppress warning in case Printlevel() == 0
      if (PrintLevel() == 0)    fMinuit->mnexcm("SET NOW",arglist,0,ierr);
   }
   if (fMinuit->fIstrat != Strategy() ) {
      arglist[0] = Strategy();
      fMinuit->mnexcm("SET STR",arglist,1,ierr);
   }

   if (Precision() > 0 &&  fMinuit->fEpsma2 != Precision() ) {
      arglist[0] = Precision();
      fMinuit->mnexcm("SET EPS",arglist,1,ierr);
   }


   // syntax of MINOS command is:  MINOS [maxcalls] [parno]
   // if parno = 0 all parameters are done
   arglist[0] = MaxFunctionCalls();
   arglist[1] = i+1;  // par number starts from 1 in TMInuit

   int nargs = 2;
   fMinuit->mnexcm("MINOS",arglist,nargs,ierr);
   bool isValid = (ierr == 0);
   // check also the status from fCstatu
   if (isValid && fMinuit->fCstatu != "SUCCESSFUL") {
      if (fMinuit->fCstatu == "FAILURE" ) {
         // in this case MINOS failed on all parameters, so it is not valid !
         ierr = 5;
         isValid = false;
      }
      if (fMinuit->fCstatu == "PROBLEMS") ierr = 6;
      ierr = 7;  // this should be the case UNCHANGED
   }

   fStatus += 10*ierr;
   fMinosStatus = ierr; 

   fMinosRun = true;

   // retrieve parameters in case a new minimum has been found
   if (fMinuit->fCstatu == "SUCCESSFUL")
      RetrieveParams();

   double errParab = 0;
   double gcor = 0;
   // what returns if parameter fixed or constant or at limit ?
   fMinuit->mnerrs(i,errUp,errLow, errParab, gcor);

   // do not flag errors case of PROBLEMS or UNCHANGED (
   return isValid;

}

void TMinuitMinimizer::DoClear() {
   // reset TMinuit

   fMinuit->mncler();

   //reset the internal Minuit random generator to its initial state
   double val = 3;
   int inseed = 12345;
   fMinuit->mnrn15(val,inseed);

   fUsed = false;
   fgUsed = false;

}

void TMinuitMinimizer::DoReleaseFixParameter(int ivar) {
   // check if a parameter is defined and in case it was fixed released
   // TMinuit is not able to release free parameters by redefining them
   // so we need to force the release
   if (fMinuit == 0) return;
   if (fMinuit->GetNumFixedPars() == 0) return;
   // check if parameter has already been defined
   if (int(ivar) >= fMinuit->GetNumPars() ) return;

   // check if parameter is fixed
   for (int i = 0; i < fMinuit->fNpfix; ++i) {
      if (fMinuit->fIpfix[i] == ivar+1 ) {
         // parameter is fixed
         fMinuit->Release(ivar);
         return;
      }
   }

}


void TMinuitMinimizer::PrintResults() {
   // print-out results using classic Minuit format (mnprin)
   if (fMinuit == 0) return;

   // print minimizer result
   if (PrintLevel() > 2)
      fMinuit->mnprin(4,fMinuit->fAmin);
   else
      fMinuit->mnprin(3,fMinuit->fAmin);
}

void TMinuitMinimizer::SuppressMinuitWarnings(bool nowarn) {
   // suppress Minuit2 warnings
   double arglist = 0;
   int ierr = 0;
   if (nowarn)
      fMinuit->mnexcm("SET NOW",&arglist,0,ierr);
   else
      fMinuit->mnexcm("SET WAR",&arglist,0,ierr);
}


bool TMinuitMinimizer::Contour(unsigned int ipar, unsigned int jpar, unsigned int &npoints, double * x, double * y) {
   // contour plot for parameter i and j
   // need a valid FunctionMinimum otherwise exits
   if (fMinuit == 0) {
      Error("TMinuitMinimizer::Contour"," invalid TMinuit instance");
      return false;
   }

   // set error and print level
   double arglist[1];
   int ierr = 0;
   arglist[0] = ErrorDef();
   fMinuit->mnexcm("SET Err",arglist,1,ierr);

   arglist[0] = PrintLevel()-1;
   fMinuit->mnexcm("SET PRINT",arglist,1,ierr);

   // suppress warning in case Printlevel() == 0
   if (PrintLevel() == 0)    fMinuit->mnexcm("SET NOW",arglist,0,ierr);

   // set precision if needed
   if (Precision() > 0)  {
      arglist[0] = Precision();
      fMinuit->mnexcm("SET EPS",arglist,1,ierr);
   }


   if (npoints < 4) {
      Error("TMinuitMinimizer::Contour","Cannot make contour with so few points");
      return false;
   }
   int npfound = 0;
   // parameter numbers in mncont start from zero
   fMinuit->mncont( ipar,jpar,npoints, x, y,npfound);
   if (npfound<4) {
      // mncont did go wrong
      Error("TMinuitMinimizer::Contour","Cannot find more than 4 points");
      return false;
   }
   if (npfound!=(int)npoints) {
      // mncont did go wrong
      Warning("TMinuitMinimizer::Contour","Returning only %d points ",npfound);
      npoints = npfound;
   }
   return true;

}

bool TMinuitMinimizer::Scan(unsigned int ipar, unsigned int & nstep, double * x, double * y, double xmin, double xmax) {
   // scan a parameter (variable) around the minimum value
   // the parameters must have been set before
   // if xmin=0 && xmax == 0  by default scan around 2 sigma of the error
   // if the errors  are also zero then scan from min and max of parameter range
   // (if parameters are limited Minuit scan from min and max instead of 2 sigma by default)
   // (force in that case to use errors)

   // scan is not implemented for TMinuit, the way to return the array is only via the graph
   if (fMinuit == 0) {
      Error("TMinuitMinimizer::Scan"," invalid TMinuit instance");
      return false;
   }

   // case of default xmin and xmax
   if (xmin >= xmax && (int) ipar < fMinuit->GetNumPars() ) {
      double val = 0; double err = 0;
      TString name;
      double xlow = 0; double xup = 0 ;
      int iuint = 0;
      // in mnpout index starts from ze
      fMinuit->mnpout( ipar, name, val, err, xlow, xup, iuint);
      // redefine 2 sigma for all parameters by default (TMinuit does 1 sigma and range if limited)
      if (iuint > 0 && err > 0) {
         xmin = val - 2.*err;
         xmax = val + 2 * err;
      }
   }

   double arglist[4];
   int ierr = 0;

   arglist[0] = PrintLevel()-1;
   fMinuit->mnexcm("SET PRINT",arglist,1,ierr);
   // suppress warning in case Printlevel() == 0
   if (PrintLevel() == 0)    fMinuit->mnexcm("SET NOW",arglist,0,ierr);

   // set precision if needed
   if (Precision() > 0)  {
      arglist[0] = Precision();
      fMinuit->mnexcm("SET EPS",arglist,1,ierr);
   }

   if (nstep == 0) return false;
   arglist[0] = ipar+1;  // TMinuit starts from 1
   arglist[1] = nstep;  //+2; // TMinuit deletes two points
   int nargs = 2;
   if (xmax > xmin ) {
      arglist[2] = xmin;
      arglist[3] = xmax;
      nargs = 4;
   }
   fMinuit->mnexcm("SCAN",arglist,nargs,ierr);
   if (ierr) {
      Error("TMinuitMinimizer::Scan"," Error executing command SCAN");
      return false;
   }
   // get TGraph object
   TGraph * gr = dynamic_cast<TGraph *>(fMinuit->GetPlot() );
   if (!gr) {
      Error("TMinuitMinimizer::Scan"," Error in returned graph object");
      return false;
   }
   nstep = std::min(gr->GetN(), (int) nstep);


   std::copy(gr->GetX(), gr->GetX()+nstep, x);
   std::copy(gr->GetY(), gr->GetY()+nstep, y);
   nstep = gr->GetN();
   return true;
}

bool TMinuitMinimizer::Hesse() {
   // perform calculation of Hessian

   if (fMinuit == 0) {
      Error("TMinuitMinimizer::Hesse","invalid TMinuit pointer. Need to call first SetFunction and SetVariable");
      return false;
   }


   double arglist[10];
   int ierr = 0;

   // set error and print level
   arglist[0] = ErrorDef();
   fMinuit->mnexcm("SET ERR",arglist,1,ierr);

   int printlevel = PrintLevel();
   arglist[0] = printlevel - 1;
   fMinuit->mnexcm("SET PRINT",arglist,1,ierr);

   // suppress warning in case Printlevel() == 0
   if (printlevel == 0)    fMinuit->mnexcm("SET NOW",arglist,0,ierr);

   // set precision if needed
   if (Precision() > 0)  {
      arglist[0] = Precision();
      fMinuit->mnexcm("SET EPS",arglist,1,ierr);
   }

   arglist[0] = MaxFunctionCalls();

   fMinuit->mnexcm("HESSE",arglist,1,ierr);
   fStatus += 100*ierr;

   if (ierr != 0) return false;

   // retrieve results (parameter and error matrix)
   // only if result is OK

   RetrieveParams();
   RetrieveErrorMatrix();

   return true;
}

bool TMinuitMinimizer::SetDebug(bool on) {
   // set debug mode

   if (fMinuit == 0) {
      Error("TMinuitMinimizer::SetDebug","invalid TMinuit pointer. Need to call first SetFunction and SetVariable");
      return false;
   }
   int ierr = 0;
   double arglist[1];
   arglist[0] = 1;
   if (on)
      fMinuit->mnexcm("SET DEBUG",arglist,1,ierr);
   else
      fMinuit->mnexcm("SET NODEBUG",arglist,1,ierr);

   return (ierr == 0);
}
//    } // end namespace Fit

// } // end namespace ROOT
