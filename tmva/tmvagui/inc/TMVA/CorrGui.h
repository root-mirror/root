#ifndef CorrGui__HH
#define CorrGui__HH
///////////
//New Gui for easier plotting of scatter corelations
// L. Ancu 04/04/07
////////////

#include "TControlBar.h"
#include "tmvaglob.h"

namespace TMVA{

   //   static TControlBar* CorrGui_Global__cbar = 0;

   void CorrGui(TString dataset,  TString fin = "TMVA.root", TString dirName = "InputVariables_Id", TString title = "TMVA Input Variable",
                Bool_t isRegression = kFALSE );
   void CorrGui_DeleteTBar();

}
#endif
