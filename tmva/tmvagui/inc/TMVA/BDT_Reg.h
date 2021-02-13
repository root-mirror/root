#ifndef BDT_Reg__HH
#define BDT_Reg__HH

#include <vector>

#include "TMVA/tmvaglob.h"

#include "RQ_OBJECT.h"

#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TControlBar.h"

#include "TGWindow.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"

#include "TMVA/DecisionTree.h"
#include "TMVA/Tools.h"
#include "TXMLEngine.h"
#include "TMVA/BDT.h"
// Uncomment this only if the link problem is solved. The include statement tends
// to use the ROOT classes rather than the local TMVA release
// #include "TMVA/DecisionTree.h"
// #include "TMVA/DecisionTreeNode.h"

namespace TMVA{



   class StatDialogBDTReg {

      RQ_OBJECT("StatDialogBDTReg")

         public:

      StatDialogBDTReg(TString dataset, const TGWindow* p, TString wfile,
                       TString methName = "BDT", Int_t itree = 0 );
      virtual ~StatDialogBDTReg() {
         TMVA::DecisionTreeNode::SetIsTraining(false);
         fThis = 0;
         fMain->CloseWindow();
         fMain->Cleanup();
         if(gROOT->GetListOfCanvases()->FindObject(fCanvas))
            delete fCanvas;
      }

      // draw method
      void DrawTree( Int_t itree );

      void RaiseDialog() { if (fMain) { fMain->RaiseWindow(); fMain->Layout(); fMain->MapWindow(); } }

   private:

      TGMainFrame *fMain;
      Int_t        fItree;
      Int_t        fNtrees;
      TCanvas*     fCanvas;
      TString fDataset;

      TGNumberEntry* fInput;

      TGHorizontalFrame* fButtons;
      TGTextButton* fDrawButton;
      TGTextButton* fCloseButton;

      void UpdateCanvases();

      // draw methods
      TMVA::DecisionTree* ReadTree( TString * &vars, Int_t itree );
      void                DrawNode( TMVA::DecisionTreeNode *n,
                                    Double_t x, Double_t y, Double_t xscale,  Double_t yscale, TString* vars );
      void GetNtrees();

      TString fWfile;
      TString fMethName;

   public:

      // static function for external deletion
      static void Delete() { if (fThis != 0) { delete fThis; fThis = 0; } }

      // slots
      void SetItree(); //*SIGNAL*
      void Redraw(); //*SIGNAL*
      void Close(); //*SIGNAL*

   private:

      static StatDialogBDTReg* fThis;

   };

   // ========================================================================================

   extern std::vector<TControlBar*> BDTReg_Global__cbar;

   // intermediate GUI
   void BDT_Reg(TString dataset, const TString& fin = "TMVAReg.root" );
   void BDTReg_DeleteTBar(int i);

   void BDT_Reg(TString dataset, Int_t itree, TString wfile = "", TString methName = "BDT", Bool_t useTMVAStyle = kTRUE );


}
#endif
