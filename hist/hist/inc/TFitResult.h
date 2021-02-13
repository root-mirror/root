// @(#)root/mathcore:$Id$
// Author: David Gonzalez Maline Tue Nov 10 15:01:24 2009

/*************************************************************************
 * Copyright (C) 1995-2009, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFitResult
#define ROOT_TFitResult

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFitResult                                                           //
//                                                                      //
// Provides a way to view the fit result and to store them.             //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

#include "Fit/FitResult.h"

#include "TMatrixDSym.h"

#include <string>

class TGraph;

class TFitResult:public TNamed, public ROOT::Fit::FitResult {

public:

   // Default constructor for I/O
   TFitResult(int status = 0): TNamed("TFitResult","TFitResult"),
                           ROOT::Fit::FitResult() {
      fStatus = status;
   };

   // constructor from name and title
   TFitResult(const char * name, const char * title) :
      TNamed(name,title),
      ROOT::Fit::FitResult()
   {}

   // constructor from an FitResult
   TFitResult(const ROOT::Fit::FitResult& f);

   virtual ~TFitResult() {}


   virtual void  Print(Option_t *option="") const;

   TMatrixDSym GetCovarianceMatrix() const;

   TMatrixDSym GetCorrelationMatrix() const;

   // scan likelihood value of  parameter and fill the given graph.
   using ROOT::Fit::FitResult::Scan;
   bool  Scan(unsigned int ipar, TGraph * gr, double xmin = 0, double xmax = 0);

   // create contour of two parameters around the minimum
   // pass as option confidence level:  default is a value of 0.683
   using ROOT::Fit::FitResult::Contour;
   bool  Contour(unsigned int ipar, unsigned int jpar, TGraph * gr , double confLevel = 0.683);

   using TObject::Error;

   // need to re-implement to solve conflict with TObject::Error
   double Error(unsigned int i) const {
      return ParError(i);
   }

private:
   ClassDef(TFitResult, 0);  // Class holding the result of the fit
};

namespace cling {
   std::string printValue(const TFitResult* val);
}
#endif
