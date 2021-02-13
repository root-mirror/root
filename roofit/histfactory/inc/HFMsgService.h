// Author: Stephan Hageboeck, CERN  23 Mar 2020
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOFIT_HISTFACTORY_INC_HFMSGSERVICE_H_
#define ROOFIT_HISTFACTORY_INC_HFMSGSERVICE_H_

#include "RooMsgService.h"

// Shortcut definitions to issue HistFactory messages through the RooMsgService.
#define cxcoutDHF oocxcoutD((TObject*)nullptr, HistFactory)
#define cxcoutIHF oocxcoutI((TObject*)nullptr, HistFactory)
#define cxcoutPHF oocxcoutP((TObject*)nullptr, HistFactory)
#define cxcoutWHF oocxcoutW((TObject*)nullptr, HistFactory)
#define cxcoutEHF oocxcoutE((TObject*)nullptr, HistFactory)
#define cxcoutFHF oocxcoutF((TObject*)nullptr, HistFactory)

#endif /* ROOFIT_HISTFACTORY_INC_HFMSGSERVICE_H_ */
