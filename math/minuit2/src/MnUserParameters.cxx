// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnUserParameters.h"

namespace ROOT {

namespace Minuit2 {

MnUserParameters::MnUserParameters(const std::vector<double> &par, const std::vector<double> &err)
   : fTransformation(par, err)
{
}

// all implemented forwarding to MnUserTransformation class

const std::vector<MinuitParameter> &MnUserParameters::Parameters() const
{
   // return vector of MinuitParameter objects
   return fTransformation.Parameters();
}

std::vector<double> MnUserParameters::Params() const
{
   // return std::vector of double with parameter values
   return fTransformation.Params();
}

std::vector<double> MnUserParameters::Errors() const
{
   // return std::vector of double with parameter errors
   return fTransformation.Errors();
}

const MinuitParameter &MnUserParameters::Parameter(unsigned int n) const
{
   // return the MinuitParameter object for index n (external)
   return fTransformation.Parameter(n);
}

bool MnUserParameters::Add(const std::string &name, double val, double err)
{
   // add a new unlimited parameter giving name, value and err (step size)
   // return false if parameter already exists
   return fTransformation.Add(name, val, err);
}

bool MnUserParameters::Add(const std::string &name, double val, double err, double low, double up)
{
   // add a new limited parameter giving name, value, err (step size) and lower/upper limits
   // return false if parameter already exists
   return fTransformation.Add(name, val, err, low, up);
}

bool MnUserParameters::Add(const std::string &name, double val)
{
   // add a new unlimited parameter giving name and value
   // return false if parameter already exists
   return fTransformation.Add(name, val);
}

void MnUserParameters::Fix(unsigned int n)
{
   // fix parameter n
   fTransformation.Fix(n);
}

void MnUserParameters::Release(unsigned int n)
{
   // release parameter n
   fTransformation.Release(n);
}

void MnUserParameters::RemoveLimits(unsigned int n)
{
   // remove limits for parameter n
   fTransformation.RemoveLimits(n);
}

void MnUserParameters::SetValue(unsigned int n, double val)
{
   // set value for parameter n
   fTransformation.SetValue(n, val);
}

void MnUserParameters::SetError(unsigned int n, double err)
{
   // set error for parameter n
   fTransformation.SetError(n, err);
}

void MnUserParameters::SetLimits(unsigned int n, double low, double up)
{
   // set limits (lower/upper) for parameter n
   fTransformation.SetLimits(n, low, up);
}

void MnUserParameters::SetUpperLimit(unsigned int n, double up)
{
   // set upper limit for parameter n
   fTransformation.SetUpperLimit(n, up);
}

void MnUserParameters::SetLowerLimit(unsigned int n, double low)
{
   // set lower limit for parameter n
   fTransformation.SetLowerLimit(n, low);
}

void MnUserParameters::SetName(unsigned int n, const std::string &name)
{
   // set name for parameter n
   fTransformation.SetName(n, name);
}

double MnUserParameters::Value(unsigned int n) const
{
   // get value for parameter n
   return fTransformation.Value(n);
}

double MnUserParameters::Error(unsigned int n) const
{
   // get error for parameter n
   return fTransformation.Error(n);
}

// interface using  parameter name

void MnUserParameters::Fix(const std::string &name)
{
   // fix parameter
   Fix(Index(name));
}

void MnUserParameters::Release(const std::string &name)
{
   // release parameter
   Release(Index(name));
}

void MnUserParameters::SetValue(const std::string &name, double val)
{
   // set value for parameter
   SetValue(Index(name), val);
}

void MnUserParameters::SetError(const std::string &name, double err)
{
   // set error
   SetError(Index(name), err);
}

void MnUserParameters::SetLimits(const std::string &name, double low, double up)
{
   // set lower/upper limits
   SetLimits(Index(name), low, up);
}

void MnUserParameters::SetUpperLimit(const std::string &name, double up)
{
   // set upper limit
   fTransformation.SetUpperLimit(Index(name), up);
}

void MnUserParameters::SetLowerLimit(const std::string &name, double low)
{
   // set lower limit
   fTransformation.SetLowerLimit(Index(name), low);
}

void MnUserParameters::RemoveLimits(const std::string &name)
{
   // remove limits
   RemoveLimits(Index(name));
}

double MnUserParameters::Value(const std::string &name) const
{
   // get parameter value
   return Value(Index(name));
}

double MnUserParameters::Error(const std::string &name) const
{
   // get parameter error
   return Error(Index(name));
}

unsigned int MnUserParameters::Index(const std::string &name) const
{
   // get index (external) corresponding to name
   return fTransformation.Index(name);
}

const std::string &MnUserParameters::GetName(unsigned int n) const
{
   // get name corresponding to index (external)
   return fTransformation.GetName(n);
}
const char *MnUserParameters::Name(unsigned int n) const
{
   // get name corresponding to index (external)
   return fTransformation.Name(n);
}

const MnMachinePrecision &MnUserParameters::Precision() const
{
   // get global paramter precision
   return fTransformation.Precision();
}

} // namespace Minuit2

} // namespace ROOT
