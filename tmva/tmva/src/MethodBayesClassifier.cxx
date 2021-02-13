// @(#)root/tmva $Id$
// Author: Marcin ....

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodBayesClassifier                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header file for description)                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Abhishek Narain, <narainabhi@gmail.com> - University of Houston           *
 *                                                                                *
 * Copyright (c) 2005-2006:                                                       *
 *      University of Houston,                                                    *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodBayesClassifier
\ingroup TMVA

Description of bayesian classifiers.

*/

#include "TMVA/MethodBayesClassifier.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include <iostream>
#include "TString.h"

REGISTER_METHOD(BayesClassifier)

ClassImp(TMVA::MethodBayesClassifier);

////////////////////////////////////////////////////////////////////////////////
/// standard constructor

   TMVA::MethodBayesClassifier::MethodBayesClassifier( const TString& jobName,
                                                       const TString& methodTitle,
                                                       DataSetInfo& theData,
                                                       const TString& theOption ) :
   TMVA::MethodBase( jobName, Types::kBayesClassifier, methodTitle, theData, theOption)
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor from weight file

TMVA::MethodBayesClassifier::MethodBayesClassifier( DataSetInfo& theData,
                                                    const TString& theWeightFile) :
   TMVA::MethodBase( Types::kBayesClassifier, theData, theWeightFile)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Variable can handle classification with 2 classes

Bool_t TMVA::MethodBayesClassifier::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ )
{
   if( type == Types::kClassification && numberClasses == 2 ) return kTRUE;
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// default initialisation

void TMVA::MethodBayesClassifier::Init( void )
{
}

////////////////////////////////////////////////////////////////////////////////
/// define the options (their key words) that can be set in the option string

void TMVA::MethodBayesClassifier::DeclareOptions()
{
}

////////////////////////////////////////////////////////////////////////////////
/// the option string is decoded, for available options see "DeclareOptions"

void TMVA::MethodBayesClassifier::ProcessOptions()
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodBayesClassifier::~MethodBayesClassifier( void )
{
}

////////////////////////////////////////////////////////////////////////////////
/// some training

void TMVA::MethodBayesClassifier::Train( void )
{
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBayesClassifier::AddWeightsXMLTo( void* /*parent*/ ) const {
   Log() << kFATAL << "Please implement writing of weights as XML" << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// read back the training results from a file (stream)

void  TMVA::MethodBayesClassifier::ReadWeightsFromStream( std::istream & )
{
}

////////////////////////////////////////////////////////////////////////////////
/// returns MVA value for given event

Double_t TMVA::MethodBayesClassifier::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   Double_t myMVA = 0;

   // cannot determine error
   NoErrorCalc(err, errUpper);

   return myMVA;
}

////////////////////////////////////////////////////////////////////////////////
/// write specific classifier response

void TMVA::MethodBayesClassifier::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   fout << "   // not implemented for class: \"" << className << "\"" << std::endl;
   fout << "};" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodBayesClassifier::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "<None>" << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance optimisation:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "<None>" << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "<None>" << Endl;
}
