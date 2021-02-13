#ifndef MovieMaker__HH
#define MovieMaker__HH
#include "TString.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TKey.h"
#include "TH2F.h"
#include "TPad.h"
#include "TText.h"

#include "network.h"
namespace TMVA{

   void DrawNetworkMovie(TString dataset, TFile* file, const TString& methodType, const TString& methodTitle );

   void DrawMLPoutputMovie(TString dataset, TFile* file, const TString& methodType, const TString& methodTitle );
   // -----------------------------------------------------------------------------

   void MovieMaker(TString dataset, TString methodType = "Method_MLP", TString methodTitle = "MLP" );
}
#endif
