#ifndef TMVARegGui__HH
#define TMVARegGui__HH

#include "TString.h"
#include "TControlBar.h"
#include "tmvaglob.h"

namespace TMVA{

   TList* RegGuiGetKeyList( const TString& pattern );

   // utility function
   void RegGuiActionButton( TControlBar* cbar,
                            const TString& title, const TString& macro, const TString& comment,
                            const TString& buttonType, TString requiredKey = "" );

   // main GUI
   void TMVARegGui( const char* fName = "TMVAReg.root", TString dataset="");
}
#endif
