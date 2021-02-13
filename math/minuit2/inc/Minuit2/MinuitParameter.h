// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MinuitParameter
#define ROOT_Minuit2_MinuitParameter

#include <algorithm>
#include <memory>
#include <cassert>
#include <string>

namespace ROOT {

namespace Minuit2 {

//____________________________________________________________________________
/**
    class for the individual Minuit Parameter with Name and number;
    contains the input numbers for the minimization or the output result
    from minimization;
    possible interactions: Fix/release, set/remove limits, set Value/error;

    From version 5.20: use string to store the name to avoid limitation of
    name length of 20 characters
 */

class MinuitParameter {

public:
   // default constructor standard with value/error = 0
   MinuitParameter()
      : fNum(0), fValue(0), fError(0.), fConst(false), fFix(false), fLoLimit(0.), fUpLimit(0.), fLoLimValid(false),
        fUpLimValid(false), fName("")
   {
   }

   // constructor for constant Parameter
   MinuitParameter(unsigned int num, const std::string &name, double val)
      : fNum(num), fValue(val), fError(0.), fConst(true), fFix(false), fLoLimit(0.), fUpLimit(0.), fLoLimValid(false),
        fUpLimValid(false), fName(name)
   {
   }

   // constructor for standard Parameter
   MinuitParameter(unsigned int num, const std::string &name, double val, double err)
      : fNum(num), fValue(val), fError(err), fConst(false), fFix(false), fLoLimit(0.), fUpLimit(0.), fLoLimValid(false),
        fUpLimValid(false), fName(name)
   {
   }

   // constructor for limited Parameter
   MinuitParameter(unsigned int num, const std::string &name, double val, double err, double min, double max)
      : fNum(num), fValue(val), fError(err), fConst(false), fFix(false), fLoLimit(min), fUpLimit(max),
        fLoLimValid(true), fUpLimValid(true), fName(name)
   {
      assert(min != max);
      if (min > max) {
         fLoLimit = max;
         fUpLimit = min;
      }
   }

   ~MinuitParameter() {}

   MinuitParameter(const MinuitParameter &par)
      : fNum(par.fNum), fValue(par.fValue), fError(par.fError), fConst(par.fConst), fFix(par.fFix),
        fLoLimit(par.fLoLimit), fUpLimit(par.fUpLimit), fLoLimValid(par.fLoLimValid), fUpLimValid(par.fUpLimValid),
        fName(par.fName)
   {
   }

   MinuitParameter &operator=(const MinuitParameter &par)
   {
      if (this != &par) {
         fNum = par.fNum;
         fName = par.fName;
         fValue = par.fValue;
         fError = par.fError;
         fConst = par.fConst;
         fFix = par.fFix;
         fLoLimit = par.fLoLimit;
         fUpLimit = par.fUpLimit;
         fLoLimValid = par.fLoLimValid;
         fUpLimValid = par.fUpLimValid;
      }
      return *this;
   }

   // access methods
   unsigned int Number() const { return fNum; }
   // new API returning a string
   const std::string &GetName() const { return fName; }
   // return const char * for mantaining backward compatibility
   const char *Name() const { return fName.c_str(); }

   double Value() const { return fValue; }
   double Error() const { return fError; }

   // interaction
   void SetName(const std::string &name) { fName = name; }

   void SetValue(double val) { fValue = val; }
   void SetError(double err) { fError = err; }
   void SetLimits(double low, double up)
   {
      assert(low != up);
      fLoLimit = low;
      fUpLimit = up;
      fLoLimValid = true;
      fUpLimValid = true;
      if (low > up) {
         fLoLimit = up;
         fUpLimit = low;
      }
   }

   void SetUpperLimit(double up)
   {
      fLoLimit = 0.;
      fUpLimit = up;
      fLoLimValid = false;
      fUpLimValid = true;
   }

   void SetLowerLimit(double low)
   {
      fLoLimit = low;
      fUpLimit = 0.;
      fLoLimValid = true;
      fUpLimValid = false;
   }

   void RemoveLimits()
   {
      fLoLimit = 0.;
      fUpLimit = 0.;
      fLoLimValid = false;
      fUpLimValid = false;
   }

   void Fix() { fFix = true; }
   void Release() { fFix = false; }

   // state of Parameter (fixed/const/limited)
   bool IsConst() const { return fConst; }
   bool IsFixed() const { return fFix; }

   bool HasLimits() const { return fLoLimValid || fUpLimValid; }
   bool HasLowerLimit() const { return fLoLimValid; }
   bool HasUpperLimit() const { return fUpLimValid; }
   double LowerLimit() const { return fLoLimit; }
   double UpperLimit() const { return fUpLimit; }

private:
   unsigned int fNum;
   double fValue;
   double fError;
   bool fConst;
   bool fFix;
   double fLoLimit;
   double fUpLimit;
   bool fLoLimValid;
   bool fUpLimValid;
   std::string fName;

private:
   //    void SetName(const std::string & name) {
   //       int l = std::min(int(strlen(name)), 11);
   //       memset(fName, 0, 11*sizeof(char));
   //       memcpy(fName, name, l*sizeof(char));
   //       fName[10] = '\0';
   //    }
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MinuitParameter
