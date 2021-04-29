// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveTypes.hxx>
#include <iostream>
using namespace ROOT::Experimental;
namespace REX = ROOT::Experimental;

/** \class REveException
\ingroup REve
Exception class thrown by Eve classes and macros.
*/

////////////////////////////////////////////////////////////////////////////////

bool REX::operator==(const TString& t, const std::string& s)
{ return (s == t.Data()); }

bool REX::operator==(const std::string&  s, const TString& t)
{ return (s == t.Data()); }

// Exc + ops

REveException REX::operator+(const REveException &s1, const std::string &s2)
{ REveException r(s1); r.append(s2); return r; }

REveException REX::operator+(const REveException &s1, const TString &s2)
{ REveException r(s1); r.append(s2.Data()); return r; }

REveException REX::operator+(const REveException &s1,  const char *s2)
{ REveException r(s1); r.append(s2); return r; }

REveException REX::operator+(const REveException &s1, ElementId_t x)
{ REveException r(s1); r.append(std::to_string(x)); return r; }

////////////////////////////////////////////////////////////////////////////////

thread_local REX::REveLog REX::gEveLog;

void REveLog::add(const char* txt) { fLog << txt; std::cout << txt << std::endl;}
void REveLog::add(const std::string& txt) { fLog << txt; std::cout << txt << std::endl;}
bool REveLog::has_contents() { return ! fLog.str().empty(); }
void REveLog::clear() { fLog.clear(); }


REveLog &REveLog::operator<<(std::ostream &(*os)(std::ostream &))
{
   fLog << os;
   std::cout << os;
   return *this;
}