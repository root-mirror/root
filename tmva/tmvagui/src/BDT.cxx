#include "TMVA/BDT.h"
#include <iostream>
#include <fstream>

#include "RQ_OBJECT.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TFile.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TObjString.h"
#include "TControlBar.h"

#include "TGWindow.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"

#include "TMVA/DecisionTree.h"
#include "TMVA/Tools.h"
#include "TXMLEngine.h"

std::vector<TControlBar*> TMVA::BDT_Global__cbar;

TMVA::StatDialogBDT* TMVA::StatDialogBDT::fThis = 0;

void TMVA::StatDialogBDT::SetItree()
{
   fItree = Int_t(fInput->GetNumber());
}

void TMVA::StatDialogBDT::Redraw()
{
   UpdateCanvases();
}

void TMVA::StatDialogBDT::Close()
{
   delete this;
}

TMVA::StatDialogBDT::StatDialogBDT(TString dataset, const TGWindow* p, TString wfile, TString methName, Int_t itree )
   : fMain( 0 ),
     fItree(itree),
     fNtrees(0),
     fCanvas(0),
     fInput(0),
     fButtons(0),
     fDrawButton(0),
     fCloseButton(0),
     fWfile( wfile ),
     fMethName( methName ),
     fDataset(dataset)
{
   UInt_t totalWidth  = 500;
   UInt_t totalHeight = 200;

   fThis = this;

   TMVA::DecisionTreeNode::SetIsTraining(true);

   // read number of decision trees from weight file
   GetNtrees();

   // main frame
   fMain = new TGMainFrame(p, totalWidth, totalHeight, kMainFrame | kVerticalFrame);

   TGLabel *sigLab = new TGLabel( fMain, Form( "Decision tree [%i-%i]",0,fNtrees-1 ) );
   fMain->AddFrame(sigLab, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,5,5));

   fInput = new TGNumberEntry(fMain, (Double_t) fItree,5,-1,(TGNumberFormat::EStyle) 5);
   fMain->AddFrame(fInput, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,5,5));
   fInput->Resize(100,24);
   fInput->SetLimits(TGNumberFormat::kNELLimitMinMax,0,fNtrees-1);

   fButtons = new TGHorizontalFrame(fMain, totalWidth,30);

   fCloseButton = new TGTextButton(fButtons,"&Close");
   fButtons->AddFrame(fCloseButton, new TGLayoutHints(kLHintsLeft | kLHintsTop));

   fDrawButton = new TGTextButton(fButtons,"&Draw");
   fButtons->AddFrame(fDrawButton, new TGLayoutHints(kLHintsRight | kLHintsTop,15));

   fMain->AddFrame(fButtons,new TGLayoutHints(kLHintsLeft | kLHintsBottom,5,5,5,5));

   fMain->SetWindowName("Decision tree");
   fMain->SetWMPosition(0,0);
   fMain->MapSubwindows();
   fMain->Resize(fMain->GetDefaultSize());
   fMain->MapWindow();

   fInput->Connect("ValueSet(Long_t)","TMVA::StatDialogBDT",this, "SetItree()");

   // doesn't seem to exist .. gives an 'error message' and seems to work just fine without ... :)
   //   fDrawButton->Connect("ValueSet(Long_t)","TGNumberEntry",fInput, "Clicked()");
   fDrawButton->Connect("Clicked()", "TMVA::StatDialogBDT", this, "Redraw()");

   fCloseButton->Connect("Clicked()", "TMVA::StatDialogBDT", this, "Close()");
}

void TMVA::StatDialogBDT::UpdateCanvases()
{
   DrawTree(fItree );
}

void TMVA::StatDialogBDT::GetNtrees()
{
   if(!fWfile.EndsWith(".xml") ){
      std::ifstream fin( fWfile );
      if (!fin.good( )) { // file not found --> Error
         cout << "*** ERROR: Weight file: " << fWfile << " does not exist" << endl;
         return;
      }

      TString dummy = "";

      // read total number of trees, and check whether requested tree is in range
      Int_t nc = 0;
      while (!dummy.Contains("NTrees")) {
         fin >> dummy;
         nc++;
         if (nc > 200) {
            cout << endl;
            cout << "*** Huge problem: could not locate term \"NTrees\" in BDT weight file: "
                 << fWfile << endl;
            cout << "==> panic abort (please contact the TMVA authors)" << endl;
            cout << endl;
            exit(1);
         }
      }
      fin >> dummy;
      fNtrees = dummy.ReplaceAll("\"","").Atoi();
      fin.close();
   }
   else{
      void* doc = TMVA::gTools().xmlengine().ParseFile(fWfile);
      void* rootnode = TMVA::gTools().xmlengine().DocGetRootElement(doc);
      void* ch = TMVA::gTools().xmlengine().GetChild(rootnode);
      while(ch){
         TString nodeName = TString( TMVA::gTools().xmlengine().GetNodeName(ch) );
         if(nodeName=="Weights") {
            TMVA::gTools().ReadAttr( ch, "NTrees", fNtrees );
            break;
         }
         ch = TMVA::gTools().xmlengine().GetNext(ch);
      }
   }
   cout << "--- Found " << fNtrees << " decision trees in weight file" << endl;

}

////////////////////////////////////////////////////////////////////////////////
/// recursively puts an entries in the histogram for the node and its daughters
///

void TMVA::StatDialogBDT::DrawNode( TMVA::DecisionTreeNode *n,
                                    Double_t x, Double_t y,
                                    Double_t xscale,  Double_t yscale, TString * vars)
{
   Float_t xsize=xscale*1.5;
   Float_t ysize=yscale/3;
   if (xsize>0.15) xsize=0.1; //xscale/2;
   if (n->GetLeft() != NULL){
      TLine *a1 = new TLine(x-xscale/4,y-ysize,x-xscale,y-ysize*2);
      a1->SetLineWidth(2);
      a1->Draw();
      DrawNode((TMVA::DecisionTreeNode*) n->GetLeft(), x-xscale, y-yscale, xscale/2, yscale, vars);
   }
   if (n->GetRight() != NULL){
      TLine *a1 = new TLine(x+xscale/4,y-ysize,x+xscale,y-ysize*2);
      a1->SetLineWidth(2);
      a1->Draw();
      DrawNode((TMVA::DecisionTreeNode*) n->GetRight(), x+xscale, y-yscale, xscale/2, yscale, vars  );
   }

   //   TPaveText *t = new TPaveText(x-xscale/2,y-yscale/2,x+xscale/2,y+yscale/2, "NDC");
   TPaveText *t = new TPaveText(x-xsize,y-ysize,x+xsize,y+ysize, "NDC");

   t->SetBorderSize(1);

   t->SetFillStyle(1001);


   Double_t pur=n->GetPurity();
   t->SetFillColor(fColorOffset+Int_t(pur*100));

   char buffer[25];
   sprintf( buffer, "N=%f", n->GetNEvents() );
   if (n->GetNEvents()>0) t->AddText(buffer);
   sprintf( buffer, "S/(S+B)=%4.3f", n->GetPurity() );
   t->AddText(buffer);

   if (n->GetNodeType() == 0){
      if (n->GetCutType()){
         t->AddText(TString(vars[n->GetSelector()]+">"+=::Form("%5.3g",n->GetCutValue())));
      }else{
         t->AddText(TString(vars[n->GetSelector()]+"<"+=::Form("%5.3g",n->GetCutValue())));
      }
   }

   t->Draw();

   return;
}
TMVA::DecisionTree* TMVA::StatDialogBDT::ReadTree( TString* &vars, Int_t itree )
{
   cout << "--- Reading Tree " << itree << " from weight file: " << fWfile << endl;
   TMVA::DecisionTree *d = new TMVA::DecisionTree();
   if(!fWfile.EndsWith(".xml") ){
      std::ifstream fin( fWfile );
      if (!fin.good( )) { // file not found --> Error
         cout << "*** ERROR: Weight file: " << fWfile << " does not exist" << endl;
         delete d;
         d = nullptr;
         return 0;
      }

      TString dummy = "";

      if (itree >= fNtrees) {
         cout << "*** ERROR: requested decision tree: " << itree
              << ", but number of trained trees only: " << fNtrees << endl;
         delete d;
         d = nullptr;
         return 0;
      }

      // file header with name
      while (!dummy.Contains("#VAR")) fin >> dummy;
      fin >> dummy >> dummy >> dummy; // the rest of header line

      // number of variables
      Int_t nVars;
      fin >> dummy >> nVars;

      // variable mins and maxes
      vars = new TString[nVars+1]; // last one is if "fisher cut criterium"
      for (Int_t i = 0; i < nVars; i++) fin >> vars[i] >> dummy >> dummy >> dummy >> dummy;
      vars[nVars]="FisherCrit";

      char buffer[20];
      char line[256];
      sprintf(buffer,"Tree %d",itree);

      while (!dummy.Contains(buffer)) {
         fin.getline(line,256);
         dummy = TString(line);
      }

      d->Read(fin);

      fin.close();
   }
   else{
      if (itree >= fNtrees) {
         cout << "*** ERROR: requested decision tree: " << itree
              << ", but number of trained trees only: " << fNtrees << endl;
         delete d;
         d = nullptr;
         return 0;
      }
      Int_t nVars;
      void* doc = TMVA::gTools().xmlengine().ParseFile(fWfile);
      void* rootnode = TMVA::gTools().xmlengine().DocGetRootElement(doc);
      void* ch = TMVA::gTools().xmlengine().GetChild(rootnode);
      while(ch){
         TString nodeName = TString( TMVA::gTools().xmlengine().GetNodeName(ch) );
         if(nodeName=="Variables"){
            TMVA::gTools().ReadAttr( ch, "NVar", nVars);
            vars = new TString[nVars+1];
            void* varnode =  TMVA::gTools().xmlengine().GetChild(ch);
            for (Int_t i = 0; i < nVars; i++){
               TMVA::gTools().ReadAttr( varnode, "Expression", vars[i]);
               varnode =  TMVA::gTools().xmlengine().GetNext(varnode);
            }
            vars[nVars]="FisherCrit";
         }
         if(nodeName=="Weights") break;
         ch = TMVA::gTools().xmlengine().GetNext(ch);
      }
      ch = TMVA::gTools().xmlengine().GetChild(ch);
      for (int i=0; i<itree; i++) ch = TMVA::gTools().xmlengine().GetNext(ch);
      d->ReadXML(ch);
   }
   return d;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::StatDialogBDT::DrawTree( Int_t itree )
{
   TString *vars;
   TMVA::DecisionTree* d = ReadTree( vars, itree );
   if (d == 0) return;

   UInt_t   depth = d->GetTotalTreeDepth();
   Double_t ystep = 1.0/(depth + 1.0);

   cout << "--- Tree depth: " << depth << endl;

   TStyle* TMVAStyle   = gROOT->GetStyle("Plain"); // our style is based on Plain



   Double_t r[2]    = {1., 0.};
   Double_t g[2]    = {0., 0.};
   Double_t b[2]    = {0., 1.};
   Double_t stop[2] = {0., 1.0};
   fColorOffset = TColor::CreateGradientColorTable(2, stop, r, g, b, 100);

   Int_t MyPalette[100];
   for (int i=0;i<100;i++) MyPalette[i] = fColorOffset+i;
   TMVAStyle->SetPalette(100, MyPalette);



   Int_t   canvasColor = TMVAStyle->GetCanvasColor(); // backup

   TString cbuffer = Form( "Reading weight file: %s", fWfile.Data() );
   TString tbuffer = Form( "Decision Tree no.: %d", itree );
   if (!fCanvas) fCanvas = new TCanvas( "c1", cbuffer, 200, 0, 1000, 600 );
   else          fCanvas->Clear();
   fCanvas->Draw();

   DrawNode( (TMVA::DecisionTreeNode*)d->GetRoot(), 0.5, 1.-0.5*ystep, 0.25, ystep ,vars);

   // make the legend
   Double_t yup=0.99;
   Double_t ydown=yup-ystep/2.5;
   Double_t dy= ystep/2.5 * 0.2;

   TPaveText *whichTree = new TPaveText(0.85,ydown,0.98,yup, "NDC");
   whichTree->SetBorderSize(1);
   whichTree->SetFillStyle(1001);
   whichTree->SetFillColor( TColor::GetColor( "#ffff33" ) );
   whichTree->AddText( tbuffer );
   whichTree->Draw();

   TPaveText *signalleaf = new TPaveText(0.02,ydown ,0.15,yup, "NDC");
   signalleaf->SetBorderSize(1);
   signalleaf->SetFillStyle(1001);
   signalleaf->SetFillColor( getSigColorF() );
   signalleaf->AddText("Pure Signal Nodes");
   signalleaf->SetTextColor( getSigColorT() );
   signalleaf->Draw();

   ydown = ydown - ystep/2.5 -dy;
   yup   = yup - ystep/2.5 -dy;
   TPaveText *backgroundleaf = new TPaveText(0.02,ydown,0.15,yup, "NDC");
   backgroundleaf->SetBorderSize(1);
   backgroundleaf->SetFillStyle(1001);
   backgroundleaf->SetFillColor( getBkgColorF() );

   backgroundleaf->AddText("Pure Backgr. Nodes");
   backgroundleaf->SetTextColor( getBkgColorT() );
   backgroundleaf->Draw();


   fCanvas->Update();
   TString fname = fDataset+Form("/plots/%s_%i", fMethName.Data(), itree );
   cout << "--- Creating image: " << fname << endl;
   TMVAGlob::imgconv( fCanvas, fname );

   TMVAStyle->SetCanvasColor( canvasColor );
}

// ========================================================================================


// intermediate GUI
void TMVA::BDT(TString dataset, const TString& fin  )
{
   // --- read the available BDT weight files

   // destroy all open cavases
   TMVAGlob::DestroyCanvases();

   // checks if file with name "fin" is already open, and if not opens one
   TFile* file = TMVAGlob::OpenFile( fin );

   TDirectory* dir = file->GetDirectory(dataset.Data())->GetDirectory( "Method_BDT" );
   if (!dir) {
      cout << "*** Error in macro \"BDT.C\": cannot find directory \"Method_BDT\" in file: " << fin << endl;
      return;
   }

   // read all directories
   TIter next( dir->GetListOfKeys() );
   TKey *key(0);
   std::vector<TString> methname;
   std::vector<TString> path;
   std::vector<TString> wfile;
   while ((key = (TKey*)next())) {
      TDirectory* mdir = dir->GetDirectory( key->GetName() );
      if (!mdir) {
         cout << "*** Error in macro \"BDT.C\": cannot find sub-directory: " << key->GetName()
              << " in directory: " << dir->GetName() << endl;
         return;
      }

      // retrieve weight file name and path
      TObjString* strPath = (TObjString*)mdir->Get( "TrainingPath" );
      TObjString* strWFile = (TObjString*)mdir->Get( "WeightFileName" );
      if (!strPath || !strWFile) {
         cout << "*** Error in macro \"BDT.C\": could not find TObjStrings \"TrainingPath\" and/or \"WeightFileName\" *** " << endl;
         cout << "*** Maybe you are using TMVA >= 3.8.15 with an older training target file ? *** " << endl;
         return;
      }

      methname.push_back( key->GetName() );
      path    .push_back( strPath->GetString() );
      wfile   .push_back( strWFile->GetString() );
   }

   // create the control bar
   TControlBar* cbar = new TControlBar( "vertical", "Choose weight file:", 50, 50 );
   BDT_Global__cbar.push_back(cbar);

   for (UInt_t im=0; im<path.size(); im++) {
      TString fname = path[im];
      if (fname[fname.Length()-1] != '/') fname += "/";
      fname += wfile[im];
      TString macro = Form( "TMVA::BDT(\"%s\",0,\"%s\",\"%s\")",dataset.Data(), fname.Data(), methname[im].Data() );
      cbar->AddButton( fname, macro, "Plot decision trees from this weight file", "button" );
   }

   // set the style
   cbar->SetTextColor("blue");

   // draw
   cbar->Show();
}

void TMVA::BDT_DeleteTBar(int i)
{
   // destroy all open canvases
   StatDialogBDT::Delete();
   TMVAGlob::DestroyCanvases();

   delete BDT_Global__cbar[i];
   BDT_Global__cbar[i] = 0;
}

// input: - No. of tree
//        - the weight file from which the tree is read
void TMVA::BDT(TString dataset, Int_t itree, TString wfile , TString methName , Bool_t useTMVAStyle  )
{
   // destroy possibly existing dialog windows and/or canvases
   StatDialogBDT::Delete();
   TMVAGlob::DestroyCanvases();
   if(wfile=="")
      wfile = dataset+"/weights/TMVAnalysis_test_BDT.weights.txt";
   // quick check if weight file exist
   if(!wfile.EndsWith(".xml") ){
      std::ifstream fin( wfile );
      if (!fin.good( )) { // file not found --> Error
         cout << "*** ERROR: Weight file: " << wfile << " does not exist" << endl;
         return;
      }
   }
   std::cout << "test1";
   // set style and remove existing canvas'
   TMVAGlob::Initialize( useTMVAStyle );

   StatDialogBDT* gGui = new StatDialogBDT(dataset, gClient->GetRoot(), wfile, methName, itree );

   gGui->DrawTree(itree );

   gGui->RaiseDialog();
}
