// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : ResultsClassification                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch>  - CERN, Switzerland          *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2006:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::ResultsClassification
\ingroup TMVA
Class that is the base-class for a vector of result
*/

#include "TMVA/ResultsClassification.h"

#include "TMVA/MsgLogger.h"
#include "TMVA/Results.h"
#include "TMVA/Types.h"

#include "TString.h"

#include <vector>

namespace TMVA {
   class DataSetInfo;
}


////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::ResultsClassification::ResultsClassification( const DataSetInfo* dsi, TString resultsName  )
   : Results( dsi,resultsName  ),
     fRet(1),
     fLogger( new MsgLogger(Form("ResultsClassification%s",resultsName.Data()) , kINFO) )
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::ResultsClassification::~ResultsClassification()
{
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// set MVA response

void TMVA::ResultsClassification::SetValue( Float_t value, Int_t ievt, Bool_t type )
{
   assert(ievt < (int) fMvaValues.size() && ievt < (int) fMvaValuesTypes.size());
   fMvaValues[ievt] = value;
   fMvaValuesTypes[ievt] = type;
}
