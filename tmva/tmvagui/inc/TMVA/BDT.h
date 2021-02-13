#ifndef BDT__HH
#define BDT__HH

#include <vector>

#include "tmvaglob.h"

#include "RQ_OBJECT.h"

#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
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

// Uncomment this only if the link problem is solved. The include statement tends
// to use the ROOT classes rather than the local TMVA release
// #include "TMVA/DecisionTree.h"
// #include "TMVA/DecisionTreeNode.h"
namespace TMVA{

   // this macro displays a decision tree read in from the weight file


   inline Int_t getSigColorF () {return  TColor::GetColor( "#0000FF" );}  // Pure Signal
   inline Int_t getBkgColorF () {return  TColor::GetColor( "#FF0000" );}  // Pure Backgr.
   inline Int_t getIntColorF () {return  TColor::GetColor( "#33aa77" );}  // novel green


   inline Int_t getSigColorT () {return  10;}
   inline Int_t getBkgColorT () {return  10;}
   inline Int_t getIntColorT () {return  10;}



   class StatDialogBDT {

      RQ_OBJECT("StatDialogBDT")

         public:

      StatDialogBDT(TString dataset, const TGWindow* p, TString wfile,
                    TString methName = "BDT", Int_t itree = 0 );
      virtual ~StatDialogBDT() {
         TMVA::DecisionTreeNode::SetIsTraining(false);
         fThis = 0;
         fMain->CloseWindow();
         fMain->Cleanup();
         if(gROOT->GetListOfCanvases()->FindObject(fCanvas))
            delete fCanvas;
      }

      // draw method
      void DrawTree(Int_t itree );

      void RaiseDialog() { if (fMain) { fMain->RaiseWindow(); fMain->Layout(); fMain->MapWindow(); } }

   private:

      TGMainFrame *fMain;
      Int_t        fItree;
      Int_t        fNtrees;
      TCanvas*     fCanvas;


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
      TString fDataset;
      Int_t   fColorOffset;

   public:

      // static function for external deletion
      static void Delete() { if (fThis != 0) { delete fThis; fThis = 0; } }

      // slots
      void SetItree(); //*SIGNAL*
      void Redraw(); //*SIGNAL*
      void Close(); //*SIGNAL*

   private:

      static StatDialogBDT* fThis;

   };
   // ========================================================================================

   extern std::vector<TControlBar*> BDT_Global__cbar;

   // intermediate GUI
   void BDT(TString dataset, const TString& fin = "TMVA.root" );
   void BDT_DeleteTBar(int i);
   void BDT(TString dataset, Int_t itree, TString wfile , TString methName = "BDT", Bool_t useTMVAStyle = kTRUE ) ;

}
#endif
