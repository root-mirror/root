// @(#)root/tmva $Id$
// Author: Matt Jachowski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::TActivation                                                     *
 *                                                                                *
 * Description:                                                                   *
 *      Interface for TNeuron activation function classes.                        *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Matt Jachowski  <jachowski@stanford.edu> - Stanford University, USA       *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/


#ifndef ROOT_TMVA_TActivation
#define ROOT_TMVA_TActivation

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TActivation                                                          //
//                                                                      //
// Interface for TNeuron activation function classes                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iosfwd>

#include "Rtypes.h"
#include "TString.h"

namespace TMVA {

   class TActivation {

   public:

      TActivation() {}
      virtual ~TActivation() {}

      // evaluate the activation function
      virtual Double_t Eval(Double_t arg) = 0;

      // evaluate the derivative of the activation function
      virtual Double_t EvalDerivative(Double_t arg) = 0;

      // minimum of the range of activation function
      virtual Double_t GetMin() = 0;

      // maximum of the range of the activation function
      virtual Double_t GetMax() = 0;

      // expression for activation function
      virtual TString GetExpression() = 0;

      // writer of function code
      virtual void MakeFunction(std::ostream& fout, const TString& fncName) = 0;

      ClassDef(TActivation,0); // Interface for TNeuron activation function classes

   };

} // namespace TMVA

#endif
