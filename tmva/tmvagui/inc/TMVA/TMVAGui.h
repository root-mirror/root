#ifndef TMVAGui__HH
#define TMVAGui__HH


#include "TList.h"
#include "TKey.h"
#include "TString.h"
#include "TControlBar.h"

#include "tmvaglob.h"

namespace TMVA{


   TList* GetKeyList( const TString& pattern );

   // utility function
   void ActionButton( TControlBar* cbar, 
                      const TString& title, const TString& macro, const TString& comment, 
                      const TString& buttonType, TString requiredKey = ""); 

   // main GUI
   void TMVAGui( const char* fName = "TMVA.root",TString dataset = "");

   struct  TMVAGUI {
      TMVAGUI(TString name = "TMVA.root", TString dataset="") {
         TMVA::TMVAGui(name.Data(),dataset);
      }
   };
   
}


#endif
