// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnUserTransformation.h"
#include "Minuit2/MnUserCovariance.h"
#include "Minuit2/MnMatrix.h"

#include <algorithm>
#include <cstdio>
#include <string>
#include <sstream>

#include <iostream>

namespace ROOT {

namespace Minuit2 {

class MnParStr {

public:
   MnParStr(const std::string &name) : fName(name) {}

   ~MnParStr() {}

   bool operator()(const MinuitParameter &par) const
   {
      //      return (strcmp(par.Name(), fName) == 0);
      return par.GetName() == fName;
   }

private:
   const std::string &fName;
};

MnUserTransformation::MnUserTransformation(const std::vector<double> &par, const std::vector<double> &err)
   : fPrecision(MnMachinePrecision()), fParameters(std::vector<MinuitParameter>()),
     fExtOfInt(std::vector<unsigned int>()), fDoubleLimTrafo(SinParameterTransformation()),
     fUpperLimTrafo(SqrtUpParameterTransformation()), fLowerLimTrafo(SqrtLowParameterTransformation()),
     fCache(std::vector<double>())
{
   // constructor from a vector of parameter values and a vector of errors (step  sizes)
   // class has as data member the transformation objects (all of the types),
   // the std::vector of MinuitParameter objects and the vector with the index conversions from
   // internal to external (fExtOfInt)

   fParameters.reserve(par.size());
   fExtOfInt.reserve(par.size());
   fCache.reserve(par.size());

   std::string parName;
   for (unsigned int i = 0; i < par.size(); i++) {
      std::ostringstream buf;
      buf << "p" << i;
      parName = buf.str();
      Add(parName, par[i], err[i]);
   }
}

//#ifdef MINUIT2_THREAD_SAFE
//  this if a thread-safe implementation needed if want to share transformation object between the threads
std::vector<double> MnUserTransformation::operator()(const MnAlgebraicVector &pstates) const
{
   // transform an internal  Minuit vector of internal values in a std::vector of external values
   // fixed parameters will have their fixed values
   unsigned int n = pstates.size();
   // need to initialize to the stored (initial values) parameter  values for the fixed ones
   std::vector<double> pcache(fCache);
   for (unsigned int i = 0; i < n; i++) {
      if (fParameters[fExtOfInt[i]].HasLimits()) {
         pcache[fExtOfInt[i]] = Int2ext(i, pstates(i));
      } else {
         pcache[fExtOfInt[i]] = pstates(i);
      }
   }
   return pcache;
}

// #else
// const std::vector<double> & MnUserTransformation::operator()(const MnAlgebraicVector& pstates) const {
//    // transform an internal  Minuit vector of internal values in a std::vector of external values
//    // std::vector<double> Cache(pstates.size() );
//    for(unsigned int i = 0; i < pstates.size(); i++) {
//       if(fParameters[fExtOfInt[i]].HasLimits()) {
//          fCache[fExtOfInt[i]] = Int2ext(i, pstates(i));
//       } else {
//          fCache[fExtOfInt[i]] = pstates(i);
//       }
//    }

//    return fCache;
// }
// #endif

double MnUserTransformation::Int2ext(unsigned int i, double val) const
{
   // return external value from internal value for parameter i
   if (fParameters[fExtOfInt[i]].HasLimits()) {
      if (fParameters[fExtOfInt[i]].HasUpperLimit() && fParameters[fExtOfInt[i]].HasLowerLimit())
         return fDoubleLimTrafo.Int2ext(val, fParameters[fExtOfInt[i]].UpperLimit(),
                                        fParameters[fExtOfInt[i]].LowerLimit());
      else if (fParameters[fExtOfInt[i]].HasUpperLimit() && !fParameters[fExtOfInt[i]].HasLowerLimit())
         return fUpperLimTrafo.Int2ext(val, fParameters[fExtOfInt[i]].UpperLimit());
      else
         return fLowerLimTrafo.Int2ext(val, fParameters[fExtOfInt[i]].LowerLimit());
   }

   return val;
}

double MnUserTransformation::Int2extError(unsigned int i, double val, double err) const
{
   // return external error from internal error for parameter i

   // err = sigma Value == std::sqrt(cov(i,i))
   double dx = err;

   if (fParameters[fExtOfInt[i]].HasLimits()) {
      double ui = Int2ext(i, val);
      double du1 = Int2ext(i, val + dx) - ui;
      double du2 = Int2ext(i, val - dx) - ui;
      if (fParameters[fExtOfInt[i]].HasUpperLimit() && fParameters[fExtOfInt[i]].HasLowerLimit()) {
         //       double al = fParameters[fExtOfInt[i]].Lower();
         //       double ba = fParameters[fExtOfInt[i]].Upper() - al;
         //       double du1 = al + 0.5*(sin(val + dx) + 1.)*ba - ui;
         //       double du2 = al + 0.5*(sin(val - dx) + 1.)*ba - ui;
         //       if(dx > 1.) du1 = ba;
         if (dx > 1.)
            du1 = fParameters[fExtOfInt[i]].UpperLimit() - fParameters[fExtOfInt[i]].LowerLimit();
         dx = 0.5 * (std::fabs(du1) + std::fabs(du2));
      } else {
         dx = 0.5 * (std::fabs(du1) + std::fabs(du2));
      }
   }

   return dx;
}

MnUserCovariance
MnUserTransformation::Int2extCovariance(const MnAlgebraicVector &vec, const MnAlgebraicSymMatrix &cov) const
{
   // return the external covariance matrix from the internal error matrix and the internal parameter value
   // the vector of internal parameter is needed for the derivatives (Jacobian of the transformation)
   // Vext(i,j) = Vint(i,j) * dPext(i)/dPint(i) * dPext(j)/dPint(j)

   MnUserCovariance result(cov.Nrow());
   for (unsigned int i = 0; i < vec.size(); i++) {
      double dxdi = 1.;
      if (fParameters[fExtOfInt[i]].HasLimits()) {
         //       dxdi = 0.5*std::fabs((fParameters[fExtOfInt[i]].Upper() -
         //       fParameters[fExtOfInt[i]].Lower())*cos(vec(i)));
         dxdi = DInt2Ext(i, vec(i));
      }
      for (unsigned int j = i; j < vec.size(); j++) {
         double dxdj = 1.;
         if (fParameters[fExtOfInt[j]].HasLimits()) {
            //   dxdj = 0.5*std::fabs((fParameters[fExtOfInt[j]].Upper() -
            //   fParameters[fExtOfInt[j]].Lower())*cos(vec(j)));
            dxdj = DInt2Ext(j, vec(j));
         }
         result(i, j) = dxdi * cov(i, j) * dxdj;
      }
      //     double diag = Int2extError(i, vec(i), std::sqrt(cov(i,i)));
      //     result(i,i) = diag*diag;
   }

   return result;
}

double MnUserTransformation::Ext2int(unsigned int i, double val) const
{
   // return the internal value for parameter i with external value val

   if (fParameters[i].HasLimits()) {
      if (fParameters[i].HasUpperLimit() && fParameters[i].HasLowerLimit())
         return fDoubleLimTrafo.Ext2int(val, fParameters[i].UpperLimit(), fParameters[i].LowerLimit(), Precision());
      else if (fParameters[i].HasUpperLimit() && !fParameters[i].HasLowerLimit())
         return fUpperLimTrafo.Ext2int(val, fParameters[i].UpperLimit(), Precision());
      else
         return fLowerLimTrafo.Ext2int(val, fParameters[i].LowerLimit(), Precision());
   }

   return val;
}

double MnUserTransformation::DInt2Ext(unsigned int i, double val) const
{
   // return the derivative of the int->ext transformation: dPext(i) / dPint(i)
   // for the parameter i with value val

   double dd = 1.;
   if (fParameters[fExtOfInt[i]].HasLimits()) {
      if (fParameters[fExtOfInt[i]].HasUpperLimit() && fParameters[fExtOfInt[i]].HasLowerLimit())
         //       dd = 0.5*std::fabs((fParameters[fExtOfInt[i]].Upper() -
         //       fParameters[fExtOfInt[i]].Lower())*cos(vec(i)));
         dd = fDoubleLimTrafo.DInt2Ext(val, fParameters[fExtOfInt[i]].UpperLimit(),
                                       fParameters[fExtOfInt[i]].LowerLimit());
      else if (fParameters[fExtOfInt[i]].HasUpperLimit() && !fParameters[fExtOfInt[i]].HasLowerLimit())
         dd = fUpperLimTrafo.DInt2Ext(val, fParameters[fExtOfInt[i]].UpperLimit());
      else
         dd = fLowerLimTrafo.DInt2Ext(val, fParameters[fExtOfInt[i]].LowerLimit());
   }

   return dd;
}

/*
 double MnUserTransformation::dExt2Int(unsigned int, double) const {
    double dd = 1.;

    if(fParameters[fExtOfInt[i]].HasLimits()) {
       if(fParameters[fExtOfInt[i]].HasUpperLimit() && fParameters[fExtOfInt[i]].HasLowerLimit())
          //       dd = 0.5*std::fabs((fParameters[fExtOfInt[i]].Upper() -
 fParameters[fExtOfInt[i]].Lower())*cos(vec(i))); dd = fDoubleLimTrafo.dExt2Int(val,
 fParameters[fExtOfInt[i]].UpperLimit(), fParameters[fExtOfInt[i]].LowerLimit()); else
 if(fParameters[fExtOfInt[i]].HasUpperLimit() && !fParameters[fExtOfInt[i]].HasLowerLimit()) dd =
 fUpperLimTrafo.dExt2Int(val, fParameters[fExtOfInt[i]].UpperLimit()); else dd = fLowerLimTrafo.dExtInt(val,
 fParameters[fExtOfInt[i]].LowerLimit());
    }

    return dd;
 }
 */

unsigned int MnUserTransformation::IntOfExt(unsigned int ext) const
{
   // return internal index given external one ext
   assert(ext < fParameters.size());
   assert(!fParameters[ext].IsFixed());
   assert(!fParameters[ext].IsConst());
   std::vector<unsigned int>::const_iterator iind = std::find(fExtOfInt.begin(), fExtOfInt.end(), ext);
   assert(iind != fExtOfInt.end());

   return (iind - fExtOfInt.begin());
}

std::vector<double> MnUserTransformation::Params() const
{
   // return std::vector of double with parameter values
   unsigned int n = fParameters.size();
   std::vector<double> result(n);
   for (unsigned int i = 0; i < n; ++i)
      result[i] = fParameters[i].Value();

   return result;
}

std::vector<double> MnUserTransformation::Errors() const
{
   // return std::vector of double with parameter errors
   std::vector<double> result;
   result.reserve(fParameters.size());
   for (std::vector<MinuitParameter>::const_iterator ipar = Parameters().begin(); ipar != Parameters().end(); ++ipar)
      result.push_back((*ipar).Error());

   return result;
}

const MinuitParameter &MnUserTransformation::Parameter(unsigned int n) const
{
   // return the MinuitParameter object for index n (external)
   assert(n < fParameters.size());
   return fParameters[n];
}

// bool MnUserTransformation::Remove(const std::string & name) {
//    // remove parameter with name
//    // useful if want to re-define a parameter
//    // return false if parameter does not exist
//    std::vector<MinuitParameter>::iterator itr = std::find_if(fParameters.begin(), fParameters.end(), MnParStr(name)
//    ); if (itr == fParameters.end() ) return false; int n = itr - fParameters.begin(); if (n < 0 || n >=
//    fParameters.size() ) return false; fParameters.erase(itr); fCache.erase( fExtOfInt.begin() + n);
//    std::vector<unsigned int>::iterator iind = std::find(fExtOfInt.begin(), fExtOfInt.end(), n);
//    if (iind != fExtOfInt.end()) fExtOfInt.erase(iind);
// }

bool MnUserTransformation::Add(const std::string &name, double val, double err)
{
   // add a new unlimited parameter giving name, value and err (step size)
   // return false if parameter already exists
   if (std::find_if(fParameters.begin(), fParameters.end(), MnParStr(name)) != fParameters.end())
      return false;
   fExtOfInt.push_back(fParameters.size());
   fCache.push_back(val);
   fParameters.push_back(MinuitParameter(fParameters.size(), name, val, err));
   return true;
}

bool MnUserTransformation::Add(const std::string &name, double val, double err, double low, double up)
{
   // add a new limited parameter giving name, value, err (step size) and lower/upper limits
   // return false if parameter already exists
   if (std::find_if(fParameters.begin(), fParameters.end(), MnParStr(name)) != fParameters.end())
      return false;
   fExtOfInt.push_back(fParameters.size());
   fCache.push_back(val);
   fParameters.push_back(MinuitParameter(fParameters.size(), name, val, err, low, up));
   return true;
}

bool MnUserTransformation::Add(const std::string &name, double val)
{
   // add a new constant parameter giving name and value
   // return false if parameter already exists
   if (std::find_if(fParameters.begin(), fParameters.end(), MnParStr(name)) != fParameters.end())
      return false;
   fCache.push_back(val);
   // costant parameter - do not add in list of internals (fExtOfInt)
   fParameters.push_back(MinuitParameter(fParameters.size(), name, val));
   return true;
}

void MnUserTransformation::Fix(unsigned int n)
{
   // fix parameter n (external index)
   assert(n < fParameters.size());
   std::vector<unsigned int>::iterator iind = std::find(fExtOfInt.begin(), fExtOfInt.end(), n);
   if (iind != fExtOfInt.end())
      fExtOfInt.erase(iind, iind + 1);
   fParameters[n].Fix();
}

void MnUserTransformation::Release(unsigned int n)
{
   // release parameter n (external index)
   assert(n < fParameters.size());
   std::vector<unsigned int>::const_iterator iind = std::find(fExtOfInt.begin(), fExtOfInt.end(), n);
   if (iind == fExtOfInt.end()) {
      fExtOfInt.push_back(n);
      std::sort(fExtOfInt.begin(), fExtOfInt.end());
   }
   fParameters[n].Release();
}

void MnUserTransformation::SetValue(unsigned int n, double val)
{
   // set value for parameter n (external index)
   assert(n < fParameters.size());
   fParameters[n].SetValue(val);
   fCache[n] = val;
}

void MnUserTransformation::SetError(unsigned int n, double err)
{
   // set error for parameter n (external index)
   assert(n < fParameters.size());
   fParameters[n].SetError(err);
}

void MnUserTransformation::SetLimits(unsigned int n, double low, double up)
{
   // set limits (lower/upper) for parameter n (external index)
   assert(n < fParameters.size());
   assert(low != up);
   fParameters[n].SetLimits(low, up);
}

void MnUserTransformation::SetUpperLimit(unsigned int n, double up)
{
   // set upper limit for parameter n (external index)
   assert(n < fParameters.size());
   fParameters[n].SetUpperLimit(up);
}

void MnUserTransformation::SetLowerLimit(unsigned int n, double lo)
{
   // set lower limit for parameter n (external index)
   assert(n < fParameters.size());
   fParameters[n].SetLowerLimit(lo);
}

void MnUserTransformation::RemoveLimits(unsigned int n)
{
   // remove limits for parameter n (external index)
   assert(n < fParameters.size());
   fParameters[n].RemoveLimits();
}

void MnUserTransformation::SetName(unsigned int n, const std::string &name)
{
   // set name for parameter n (external index)
   assert(n < fParameters.size());
   fParameters[n].SetName(name);
}

double MnUserTransformation::Value(unsigned int n) const
{
   // get value for parameter n (external index)
   assert(n < fParameters.size());
   return fParameters[n].Value();
}

double MnUserTransformation::Error(unsigned int n) const
{
   // get error for parameter n (external index)
   assert(n < fParameters.size());
   return fParameters[n].Error();
}

// interface by parameter name

void MnUserTransformation::Fix(const std::string &name)
{
   // fix parameter
   Fix(Index(name));
}

void MnUserTransformation::Release(const std::string &name)
{
   // release parameter
   Release(Index(name));
}

void MnUserTransformation::SetValue(const std::string &name, double val)
{
   // set value for parameter
   SetValue(Index(name), val);
}

void MnUserTransformation::SetError(const std::string &name, double err)
{
   // set error
   SetError(Index(name), err);
}

void MnUserTransformation::SetLimits(const std::string &name, double low, double up)
{
   // set lower/upper limits
   SetLimits(Index(name), low, up);
}

void MnUserTransformation::SetUpperLimit(const std::string &name, double up)
{
   // set upper limit
   SetUpperLimit(Index(name), up);
}

void MnUserTransformation::SetLowerLimit(const std::string &name, double lo)
{
   // set lower limit
   SetLowerLimit(Index(name), lo);
}

void MnUserTransformation::RemoveLimits(const std::string &name)
{
   // remove limits
   RemoveLimits(Index(name));
}

double MnUserTransformation::Value(const std::string &name) const
{
   // get parameter value
   return Value(Index(name));
}

double MnUserTransformation::Error(const std::string &name) const
{
   // get parameter error
   return Error(Index(name));
}

unsigned int MnUserTransformation::Index(const std::string &name) const
{
   // get index (external) corresponding to name
   std::vector<MinuitParameter>::const_iterator ipar =
      std::find_if(fParameters.begin(), fParameters.end(), MnParStr(name));
   assert(ipar != fParameters.end());
   //   return (ipar - fParameters.begin());
   return (*ipar).Number();
}

int MnUserTransformation::FindIndex(const std::string &name) const
{
   // find index (external) corresponding to name - return -1 if not found
   std::vector<MinuitParameter>::const_iterator ipar =
      std::find_if(fParameters.begin(), fParameters.end(), MnParStr(name));
   if (ipar == fParameters.end())
      return -1;
   return (*ipar).Number();
}

const std::string &MnUserTransformation::GetName(unsigned int n) const
{
   // get name corresponding to index (external)
   assert(n < fParameters.size());
   return fParameters[n].GetName();
}

const char *MnUserTransformation::Name(unsigned int n) const
{
   // get name corresponding to index (external)
   return GetName(n).c_str();
}

} // namespace Minuit2

} // namespace ROOT
