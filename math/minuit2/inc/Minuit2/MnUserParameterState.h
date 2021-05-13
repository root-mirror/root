// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnUserParameterState
#define ROOT_Minuit2_MnUserParameterState

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserCovariance.h"
#include "Minuit2/MnGlobalCorrelationCoeff.h"

#include <vector>
#include <string>

namespace ROOT {

namespace Minuit2 {

class MinimumState;

//_____________________________________________________________________________
/**
    class which holds the external user and/or internal Minuit representation
    of the parameters and errors;
    transformation internal <-> external on demand;
 */

class MnUserParameterState {

public:
   /// default constructor (invalid state)
   MnUserParameterState()
      : fValid(false), fCovarianceValid(false), fGCCValid(false), fCovStatus(-1), fFVal(0), fEDM(0), fNFcn(0),
        fParameters(MnUserParameters()), fCovariance(MnUserCovariance()), fIntParameters(std::vector<double>()),
        fIntCovariance(MnUserCovariance())
   {
   }

   /// construct from user parameters (before minimization)
   MnUserParameterState(const std::vector<double> &, const std::vector<double> &);

   MnUserParameterState(const MnUserParameters &);

   /// construct from user parameters + covariance (before minimization)
   MnUserParameterState(const std::vector<double> &, const std::vector<double> &, unsigned int);

   MnUserParameterState(const std::vector<double> &, const MnUserCovariance &);

   MnUserParameterState(const MnUserParameters &, const MnUserCovariance &);

   /// construct from internal parameters (after minimization)
   MnUserParameterState(const MinimumState &, double, const MnUserTransformation &);

   // user external representation
   const MnUserParameters &Parameters() const { return fParameters; }
   const MnUserCovariance &Covariance() const { return fCovariance; }
   const MnGlobalCorrelationCoeff &GlobalCC() const { return fGlobalCC; }

   // hessian (inverse of covariance matrix)
   MnUserCovariance Hessian() const;

   // Minuit internal representation
   const std::vector<double> &IntParameters() const { return fIntParameters; }
   const MnUserCovariance &IntCovariance() const { return fIntCovariance; }

   // covariance matrix status (0 = not valid, 1 approximate, 2, full but made pos def, 3 accurate and not pos def
   int CovarianceStatus() const { return fCovStatus; }

   // transformation internal <-> external
   const MnUserTransformation &Trafo() const { return fParameters.Trafo(); }

   bool IsValid() const { return fValid; }
   bool HasCovariance() const { return fCovarianceValid; }
   bool HasGlobalCC() const { return fGCCValid; }

   double Fval() const { return fFVal; }
   double Edm() const { return fEDM; }
   unsigned int NFcn() const { return fNFcn; }

public:
   /** facade: forward interface of MnUserParameters and MnUserTransformation */

   // access to parameters (row-wise)
   const std::vector<ROOT::Minuit2::MinuitParameter> &MinuitParameters() const;
   // access to parameters and errors in column-wise representation
   std::vector<double> Params() const;
   std::vector<double> Errors() const;

   // access to single Parameter
   const MinuitParameter &Parameter(unsigned int i) const;

   // add free Parameter
   void Add(const std::string &name, double val, double err);
   // add limited Parameter
   void Add(const std::string &name, double val, double err, double, double);
   // add const Parameter
   void Add(const std::string &, double);

   // interaction via external number of Parameter
   void Fix(unsigned int);
   void Release(unsigned int);
   void RemoveLimits(unsigned int);
   void SetValue(unsigned int, double);
   void SetError(unsigned int, double);
   void SetLimits(unsigned int, double, double);
   void SetUpperLimit(unsigned int, double);
   void SetLowerLimit(unsigned int, double);
   void SetName(unsigned int iext, const std::string &name) { fParameters.SetName(iext, name); }

   double Value(unsigned int) const;
   double Error(unsigned int) const;

   // interaction via Name of Parameter
   void Fix(const std::string &);
   void Release(const std::string &);
   void SetValue(const std::string &, double);
   void SetError(const std::string &, double);
   void SetLimits(const std::string &, double, double);
   void SetUpperLimit(const std::string &, double);
   void SetLowerLimit(const std::string &, double);
   void RemoveLimits(const std::string &);

   double Value(const std::string &) const;
   double Error(const std::string &) const;

   // convert Name into external number of Parameter
   unsigned int Index(const std::string &) const;
   // convert external number into Name of Parameter
   const std::string &GetName(unsigned int) const;
   // mantain interface with const char * for backward compatibility
   const char *Name(unsigned int) const;

   // transformation internal <-> external
   double Int2ext(unsigned int, double) const;
   double Ext2int(unsigned int, double) const;
   unsigned int IntOfExt(unsigned int) const;
   unsigned int ExtOfInt(unsigned int) const;
   unsigned int VariableParameters() const;
   const MnMachinePrecision &Precision() const;
   void SetPrecision(double eps);

private:
   bool fValid;
   bool fCovarianceValid;
   bool fGCCValid;
   int fCovStatus; // covariance matrix status
   double fFVal;
   double fEDM;
   unsigned int fNFcn;

   MnUserParameters fParameters;
   MnUserCovariance fCovariance;
   MnGlobalCorrelationCoeff fGlobalCC;

   std::vector<double> fIntParameters;
   MnUserCovariance fIntCovariance;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnUserParameterState
