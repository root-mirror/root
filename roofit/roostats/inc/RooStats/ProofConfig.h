// @(#)root/roostats:$Id$
// Author: Kyle Cranmer and Sven Kreiss  July 2010
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_ProofConfig
#define ROOSTATS_ProofConfig




#include "Rtypes.h"

#include "RooWorkspace.h"
#include "RooStudyManager.h"

#include "TROOT.h"

namespace RooStats {

/** \class ProofConfig
    \ingroup Roostats

Holds configuration options for proof and proof-lite.

This class will be expanded in the future to hold more specific configuration
options for the tools in RooStats.

Access to TProof::Mgr for configuration is still possible as usual
(e.g. to set Root Version to be used on workers). You can do:

~~~ {.cpp}
   TProof::Mgr("my.server.url")->ShowROOTVersions()
   TProof::Mgr("my.server.url")->SetROOTVersion("v5-27-06_dbg")
~~~

See doc: http://root.cern/drupal/content/changing-default-root-version
*/


class ProofConfig {

   public:

      // configure proof with number of experiments and host session
      //  in case of Prooflite, it is better to define the number of workers as "worker=n" in the host string

      ProofConfig(RooWorkspace &w, Int_t nExperiments = 0, const char *host = "", Bool_t showGui = kFALSE) :
         fWorkspace(w),
         fNExperiments(nExperiments),
         fHost(host),
         fShowGui(showGui)
      {

            // case of ProofLite
         if (fHost == "" || fHost.Contains("lite") ) {
            fLite = true;


            // get the default value of the machine - use CINT interface until we have a poper PROOF interface that we can call
            int nMaxWorkers = gROOT->ProcessLineFast("TProofLite::GetNumberOfWorkers()");

            if (nExperiments == 0) {
               fNExperiments = nMaxWorkers;
            }

            if (nExperiments > nMaxWorkers)
               std::cout << "ProofConfig - Warning: using a number of workers = " << nExperiments << " which is larger than the number of cores in the machine "
                         << nMaxWorkers << std::endl;

            // set the number of workers in the Host string
            fHost = TString::Format("workers=%d",fNExperiments);
         }
         else {
            fLite = false;
            // have always a default number of experiments
            if (nExperiments == 0) fNExperiments = 8;
         }
      }


      virtual ~ProofConfig() {
         ProofConfig::CloseProof();
      }

      /// close all proof connections
      static void CloseProof(Option_t *option = "s") { RooStudyManager::closeProof(option); }

      // returns fWorkspace
      RooWorkspace& GetWorkspace(void) const { return fWorkspace; }
      // returns fHost
      const char* GetHost(void) const { return fHost; }
      // return fNExperiments
      Int_t GetNExperiments(void) const { return fNExperiments; }
      // return fShowGui
      Bool_t GetShowGui(void) const { return fShowGui; }
      // return true if it is a Lite session (ProofLite)
      Bool_t IsLite() const { return fLite; }

   protected:
      RooWorkspace& fWorkspace;   // workspace that is to be used with the RooStudyManager
      Int_t fNExperiments;        // number of experiments. This is sometimes called "events" in proof; "experiments" in RooStudyManager.
      TString fHost;              // Proof hostname. Use empty string (ie "") for proof-lite. Can also handle options like "workers=2" to run on two nodes.
      Bool_t fShowGui;            // Whether to show the Proof Progress window.
      Bool_t fLite;               // Whether we have a Proof Lite session

   protected:
   ClassDef(ProofConfig,1) // Configuration options for proof.
};
}


#endif
