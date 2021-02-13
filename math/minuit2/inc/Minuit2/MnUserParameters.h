// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnUserParameters
#define ROOT_Minuit2_MnUserParameters

#include "Minuit2/MnUserTransformation.h"

#include <vector>
#include <string>

namespace ROOT {

namespace Minuit2 {

class MnMachinePrecision;

/** API class for the user interaction with the parameters;
    serves as input to the minimizer as well as output from it;
    users can interact: Fix/release parameters, set values and errors, etc.;
    parameters can be accessed via their Parameter number (determined
    internally by Minuit and followed the order how the parameters are created)
    or via their user-specified Name (10 character string).
    Minuit has also an internal parameter number which is used during the minimization
    (the fix parameter are skipped). The parameter number used in this class is the external
    one. The class ROOT::Minuit2::MnUserTransformation is used to keep the
    internal <-> external transformation
 */

class MnUserParameters {

public:
   MnUserParameters() : fTransformation(MnUserTransformation()) {}

   MnUserParameters(const std::vector<double> &, const std::vector<double> &);

   ~MnUserParameters() {}

   MnUserParameters(const MnUserParameters &par) : fTransformation(par.fTransformation) {}

   MnUserParameters &operator=(const MnUserParameters &par)
   {
      fTransformation = par.fTransformation;
      return *this;
   }

   const MnUserTransformation &Trafo() const { return fTransformation; }

   unsigned int VariableParameters() const { return fTransformation.VariableParameters(); }

   /// access to parameters (row-wise)
   const std::vector<ROOT::Minuit2::MinuitParameter> &Parameters() const;

   /// access to parameters and errors in column-wise representation
   std::vector<double> Params() const;
   std::vector<double> Errors() const;

   /// access to single Parameter
   const MinuitParameter &Parameter(unsigned int) const;

   /// Add free Parameter Name, Value, Error
   bool Add(const std::string &, double, double);
   /// Add limited Parameter Name, Value, Lower bound, Upper bound
   bool Add(const std::string &, double, double, double, double);
   /// Add const Parameter Name, vale
   bool Add(const std::string &, double);

   /// interaction via external number of Parameter
   void Fix(unsigned int);
   void Release(unsigned int);
   void RemoveLimits(unsigned int);
   void SetValue(unsigned int, double);
   void SetError(unsigned int, double);
   void SetLimits(unsigned int, double, double);
   void SetUpperLimit(unsigned int, double);
   void SetLowerLimit(unsigned int, double);
   void SetName(unsigned int, const std::string &);

   double Value(unsigned int) const;
   double Error(unsigned int) const;

   /// interaction via Name of Parameter
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

   const MnMachinePrecision &Precision() const;
   void SetPrecision(double eps) { fTransformation.SetPrecision(eps); }

private:
   MnUserTransformation fTransformation;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnUserParameters
