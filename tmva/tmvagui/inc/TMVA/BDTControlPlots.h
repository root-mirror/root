#ifndef BDTControlPlots__HH
#define BDTControlPlots__HH

#include "TString.h"
#include "tmvaglob.h"

namespace TMVA{

   // input: - Input file (result from TMVA),
   //        - use of TMVA plotting TStyle
   void bdtcontrolplots(TString dataset,TDirectory *);

   void BDTControlPlots(TString dataset, TString fin = "TMVA.root", Bool_t useTMVAStyle = kTRUE );
   void bdtcontrolplots(TString dataset, TDirectory *bdtdir );

}
#endif
