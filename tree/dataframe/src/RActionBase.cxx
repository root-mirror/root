// Author: Enrico Guiraud, Danilo Piparo CERN  09/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RDF/RActionBase.hxx"
#include "ROOT/RDF/RLoopManager.hxx"

using namespace ROOT::Internal::RDF;

RActionBase::RActionBase(RLoopManager *lm, const ColumnNames_t &colNames, const RBookedDefines &defines)
   : fLoopManager(lm), fNSlots(lm->GetNSlots()), fColumnNames(colNames), fDefines(defines) { }

// outlined to pin virtual table
RActionBase::~RActionBase() {}
