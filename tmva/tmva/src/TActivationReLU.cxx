// @(#)root/tmva $Id$
// Author: Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TActivationReLU                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Rectified linear unit function  for an ANN.                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Helge Voss                                                                *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::TActivationReLU
\ingroup TMVA
Rectified Linear Unit activation function for TNeuron
*/

#include "TMVA/TActivationReLU.h"

#include "TMVA/TActivation.h"

#include "TString.h"

#include <iostream>

ClassImp(TMVA::TActivationReLU);

////////////////////////////////////////////////////////////////////////////////
/// get expressions for the tanh and its derivative

TString TMVA::TActivationReLU::GetExpression()
{
   TString expr = "max(0,x)";

   return expr;
}

////////////////////////////////////////////////////////////////////////////////
/// writes the Rectified Linear Unit activation function source code

void TMVA::TActivationReLU::MakeFunction( std::ostream& fout, const TString& fncName )
{
   fout << "double " << fncName << "(double x) const {" << std::endl;
   fout << "   // rectified linear unit" << std::endl;
   fout << "   return x>0 ? x : 0; " << std::endl;
   fout << "}" << std::endl;
}
