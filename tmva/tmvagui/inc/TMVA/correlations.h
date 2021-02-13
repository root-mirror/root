#ifndef correlations__HH
#define correlations__HH
#include "tmvaglob.h"
#include "TH1.h"
#include "TH2.h"

namespace TMVA{

   // this macro plots the correlation matrix of the various input
   // variables used in TMVA (e.g. running TMVAnalysis.C).  Signal and
   // Background are plotted separately

   // input: - Input file (result from TMVA),
   //        - use of colors or grey scale
   //        - use of TMVA plotting TStyle
   void correlations(TString dataset, TString fin = "TMVA.root", Bool_t isRegression = kFALSE, 
                     Bool_t greyScale = kFALSE, Bool_t useTMVAStyle = kTRUE );
}
#endif
