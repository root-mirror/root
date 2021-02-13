// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveMacro.h"

#include "TPRegexp.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TApplication.h"

/** \class TEveMacro
\ingroup TEve
Sub-class of TMacro, overriding Exec to unload the previous version
and cleanup after the execution.
*/

ClassImp(TEveMacro);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TEveMacro::TEveMacro() : TMacro()
{
}

TEveMacro::TEveMacro(const TEveMacro& m) : TMacro(m)
{
   // Copy constructor.
}

TEveMacro::TEveMacro(const char* name) :
   TMacro()
{
   // Constructor with file name.

   if (!name) return;

   fTitle = name;

   TPMERegexp re("([^/]+?)(?:\\.\\w*)?$");
   Int_t nm = re.Match(fTitle);
   if (nm >= 2) {
      fName = re[1];
   } else {
      fName = "<unknown>";
   }
   ReadFile(fTitle);
}


////////////////////////////////////////////////////////////////////////////////
/// Execute the macro.

Long_t TEveMacro::Exec(const char* params, Int_t* error)
{
   Long_t retval = -1;

   if (gROOT->GetGlobalFunction(fName, 0, kTRUE) != 0)
   {
      gROOT->SetExecutingMacro(kTRUE);
      gROOT->SetExecutingMacro(kFALSE);
      retval = gROOT->ProcessLine(Form("%s()", fName.Data()), error);
   }
   else
   {
      // Copy from TMacro::Exec. Difference is that the file is really placed
      // into the /tmp.
      TString fname = "/tmp/";
      {
         //the current implementation uses a file in the current directory.
         //should be replaced by a direct execution from memory by CINT
         fname += GetName();
         fname += ".C";
         SaveSource(fname);
         //disable a possible call to gROOT->Reset from the executed script
         gROOT->SetExecutingMacro(kTRUE);
         //execute script in /tmp
         TString exec = ".x " + fname;
         TString p = params;
         if (p == "") p = fParams;
         if (p != "")
            exec += "(" + p + ")";
         retval = gROOT->ProcessLine(exec, error);
         //enable gROOT->Reset
         gROOT->SetExecutingMacro(kFALSE);
         //delete the temporary file
         gSystem->Unlink(fname);
      }
   }

   //G__unloadfile(fname);

   // In case an exception was thrown (which i do not know how to detect
   // the execution of next macros does not succeed.
   // However strange this might seem, this solves the problem.
   // TTimer::SingleShot(100, "TEveMacro", this, "ResetRoot()");
   //
   // 27.8.07 - ok, this does not work any more. Seems I'll have to fix
   // this real soon now.
   //
   // !!!! FIX MACRO HANDLING !!!!
   //

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Call gROOT->Reset() via interpreter.

void TEveMacro::ResetRoot()
{
   // printf ("TEveMacro::ResetRoot doing 'gROOT->Reset()'.\n");
   gROOT->GetApplication()->ProcessLine("gROOT->Reset()");
}
