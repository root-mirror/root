// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Fredrik Tegenfeldt, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Config                                                                *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation                                                            *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker    <Andreas.Hocker@cern.ch>     - CERN, Switzerland       *
 *      Joerg Stelzer      <Joerg.Stelzer@cern.ch>      - CERN, Switzerland       *
 *      Fredrik Tegenfeldt <Fredrik.Tegenfeldt@cern.ch> - Iowa State U., USA      *
 *      Helge Voss         <Helge.Voss@cern.ch>         - MPI-K Heidelberg, GER   *
 *                                                                                *
 * Copyright (c) 2006:                                                            *
 *      CERN, Switzerland                                                         *
 *      Iowa State U., USA                                                        *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://mva.sourceforge.net/license.txt)                                       *
 **********************************************************************************/

/*! \class TMVA::Config
\ingroup TMVA

Singleton class for global configuration settings used by TMVA.

*/

#include "TMVA/Config.h"
#include "TMVA/MsgLogger.h"

#include "Rtypes.h"
#include "TString.h"

ClassImp(TMVA::Config);

#if __cplusplus > 199711L
std::atomic<TMVA::Config*> TMVA::Config::fgConfigPtr{ 0 };
#else
TMVA::Config* TMVA::Config::fgConfigPtr = 0;
#endif

TMVA::Config& TMVA::gConfig() { return TMVA::Config::Instance(); }

////////////////////////////////////////////////////////////////////////////////
/// constructor - set defaults

TMVA::Config::Config() :
   fDrawProgressBar      ( kFALSE ),
   fNWorkers             (1),
   fUseColoredConsole    ( kTRUE  ),
   fSilent               ( kFALSE ),
   fWriteOptionsReference( kFALSE ),
   fLogger               (new MsgLogger("Config"))
{
   // plotting
   fVariablePlotting.fTimesRMS = 8.0;
   fVariablePlotting.fNbins1D  = 40;
   fVariablePlotting.fNbins2D  = 300;
   fVariablePlotting.fMaxNumOfAllowedVariables = 200;
   fVariablePlotting.fMaxNumOfAllowedVariablesForScatterPlots = 20;

   fVariablePlotting.fNbinsMVAoutput   = 40;
   fVariablePlotting.fNbinsXOfROCCurve = 100;
   fVariablePlotting.fUsePaperStyle = 0;
   fVariablePlotting.fPlotFormat = VariablePlotting::kPNG;  // format for plotting (use when fUsePaperStyle ==0)

   // IO names
   fIONames.fWeightFileDirPrefix = "";
   fIONames.fWeightFileDir           = "weights";
   fIONames.fWeightFileExtension     = "weights";
   fIONames.fOptionsReferenceFileDir = "optionInfo";

}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::Config::~Config()
{
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// static function: destroy TMVA instance

void TMVA::Config::DestroyInstance()
{
#if __cplusplus > 199711L
   delete fgConfigPtr.exchange(0);
#else
   if (fgConfigPtr != 0) { delete fgConfigPtr; fgConfigPtr = 0;}
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// static function: returns  TMVA instance

TMVA::Config& TMVA::Config::Instance()
{
#if __cplusplus > 199711L
   if(!fgConfigPtr) {
      TMVA::Config* tmp = new Config();
      TMVA::Config* expected = 0;
      if(! fgConfigPtr.compare_exchange_strong(expected,tmp) ) {
         //another thread beat us to the switch
         delete tmp;
      }
   }
   return *fgConfigPtr;
#else
   return fgConfigPtr ? *fgConfigPtr :*(fgConfigPtr = new Config());
#endif
}
