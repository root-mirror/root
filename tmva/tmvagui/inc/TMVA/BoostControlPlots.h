#ifndef BoostControlPlots__HH
#define BoostControlPlots__HH

#include "tmvaglob.h"

namespace TMVA{


   // input: - Input file (result from TMVA),
   //        - use of TMVA plotting TStyle
   // this macro is based on BDTControlPlots.C
   void BoostControlPlots(TString dataset, TString fin = "TMVA.root", Bool_t useTMVAStyle = kTRUE );

   void boostcontrolplots(TString dataset, TDirectory *boostdir ); 

}
#endif
