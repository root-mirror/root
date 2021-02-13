// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnUserTransformation
#define ROOT_Minuit2_MnUserTransformation

#include "Minuit2/MnConfig.h"
#include "Minuit2/MnMatrix.h"
#include "Minuit2/MinuitParameter.h"
#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/SinParameterTransformation.h"
#include "Minuit2/SqrtLowParameterTransformation.h"
#include "Minuit2/SqrtUpParameterTransformation.h"

#include <vector>
#include <string>
#include <cassert>

namespace ROOT {

namespace Minuit2 {

class MnUserCovariance;

// class MnMachinePrecision;

/**
    class dealing with the transformation  between user specified parameters (external) and
    internal parameters used for minimization
 */

class MnUserTransformation {

public:
   MnUserTransformation()
      : fPrecision(MnMachinePrecision()), fParameters(std::vector<MinuitParameter>()),
        fExtOfInt(std::vector<unsigned int>()), fDoubleLimTrafo(SinParameterTransformation()),
        fUpperLimTrafo(SqrtUpParameterTransformation()), fLowerLimTrafo(SqrtLowParameterTransformation()),
        fCache(std::vector<double>())
   {
   }

   MnUserTransformation(const std::vector<double> &, const std::vector<double> &);

   ~MnUserTransformation() {}

   MnUserTransformation(const MnUserTransformation &trafo)
      : fPrecision(trafo.fPrecision), fParameters(trafo.fParameters), fExtOfInt(trafo.fExtOfInt),
        fDoubleLimTrafo(trafo.fDoubleLimTrafo), fUpperLimTrafo(trafo.fUpperLimTrafo),
        fLowerLimTrafo(trafo.fLowerLimTrafo), fCache(trafo.fCache)
   {
   }

   MnUserTransformation &operator=(const MnUserTransformation &trafo)
   {
      if (this != &trafo) {
         fPrecision = trafo.fPrecision;
         fParameters = trafo.fParameters;
         fExtOfInt = trafo.fExtOfInt;
         fDoubleLimTrafo = trafo.fDoubleLimTrafo;
         fUpperLimTrafo = trafo.fUpperLimTrafo;
         fLowerLimTrafo = trafo.fLowerLimTrafo;
         fCache = trafo.fCache;
      }
      return *this;
   }

   //#ifdef MINUIT2_THREAD_SAFE
   // thread-safe version (do not use cache)
   std::vector<double> operator()(const MnAlgebraicVector &) const;
   //#else // not thread safe
   //   const std::vector<double>& operator()(const MnAlgebraicVector&) const;
   //#endif

   // Index = internal Parameter
   double Int2ext(unsigned int, double) const;

   // Index = internal Parameter
   double Int2extError(unsigned int, double, double) const;

   MnUserCovariance Int2extCovariance(const MnAlgebraicVector &, const MnAlgebraicSymMatrix &) const;

   // Index = external Parameter
   double Ext2int(unsigned int, double) const;

   // Index = internal Parameter
   double DInt2Ext(unsigned int, double) const;

   //   // Index = external Parameter
   //   double dExt2Int(unsigned int, double) const;

   // Index = external Parameter
   unsigned int IntOfExt(unsigned int) const;

   // Index = internal Parameter
   unsigned int ExtOfInt(unsigned int internal) const
   {
      assert(internal < fExtOfInt.size());
      return fExtOfInt[internal];
   }

   const std::vector<MinuitParameter> &Parameters() const { return fParameters; }

   unsigned int VariableParameters() const { return static_cast<unsigned int>(fExtOfInt.size()); }

   // return initial parameter values (useful especially to get fixed parameter values)
   const std::vector<double> &InitialParValues() const { return fCache; }

   /** forwarded interface */

   const MnMachinePrecision &Precision() const { return fPrecision; }
   void SetPrecision(double eps) { fPrecision.SetPrecision(eps); }

   /// access to parameters and errors in column-wise representation

   std::vector<double> Params() const;
   std::vector<double> Errors() const;

   // access to single Parameter
   const MinuitParameter &Parameter(unsigned int) const;

   // add free Parameter
   bool Add(const std::string &, double, double);
   // add limited Parameter
   bool Add(const std::string &, double, double, double, double);
   // add const Parameter
   bool Add(const std::string &, double);

   // interaction via external number of Parameter
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

   // convert Name into external number of Parameter (will assert if parameter  is not found)
   unsigned int Index(const std::string &) const;
   // find parameter index given a name. If it is not found return a -1
   int FindIndex(const std::string &) const;

   // convert external number into Name of Parameter (will assert if index is out of range)
   const std::string &GetName(unsigned int) const;
   // mantain interface with const char * for backward compatibility
   const char *Name(unsigned int) const;

private:
   MnMachinePrecision fPrecision;

   std::vector<MinuitParameter> fParameters;
   std::vector<unsigned int> fExtOfInt;

   SinParameterTransformation fDoubleLimTrafo;
   SqrtUpParameterTransformation fUpperLimTrafo;
   SqrtLowParameterTransformation fLowerLimTrafo;

   mutable std::vector<double> fCache;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnUserTransformation
