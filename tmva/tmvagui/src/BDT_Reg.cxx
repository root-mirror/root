#include "TMVA/BDT_Reg.h"
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

std::vector<TControlBar*> TMVA::BDTReg_Global__cbar;

TMVA::StatDialogBDTReg* TMVA::StatDialogBDTReg::fThis = 0;

void TMVA::StatDialogBDTReg::SetItree() 
{
   fItree = Int_t(fInput->GetNumber());
}

void TMVA::StatDialogBDTReg::Redraw() 
{
   UpdateCanvases();
}

void TMVA::StatDialogBDTReg::Close() 
{
   delete this;
}

TMVA::StatDialogBDTReg::StatDialogBDTReg(TString dataset, const TGWindow* p, TString wfile, TString methName, Int_t itree )
   : fMain( 0 ),
     fItree(itree),
     fNtrees(0),
     fCanvas(0),
     fDataset(dataset),
     fInput(0),
     fButtons(0),
     fDrawButton(0),
     fCloseButton(0),
     fWfile( wfile ),
     fMethName( methName )
{
   UInt_t totalWidth  = 500;
   UInt_t totalHeight = 200;

   fThis = this;

   // read number of decision trees from weight file
   GetNtrees();

   // main frame
   fMain = new TGMainFrame(p, totalWidth, totalHeight, kMainFrame | kVerticalFrame);

   TGLabel *sigLab = new TGLabel( fMain, Form( "Regression tree [%i-%i]",0,fNtrees-1 ) );
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

   fMain->SetWindowName("Regression tree");
   fMain->SetWMPosition(0,0);
   fMain->MapSubwindows();
   fMain->Resize(fMain->GetDefaultSize());
   fMain->MapWindow();

   fInput->Connect("ValueSet(Long_t)","TMVA::StatDialogBDTReg",this, "SetItree()");
   
   // doesn't seem to exist .. gives an 'error message' and seems to work just fine without ... :)
   //   fDrawButton->Connect("Clicked()","TGNumberEntry",fInput, "ValueSet(Long_t)");
   fDrawButton->Connect("Clicked()", "TMVA::StatDialogBDTReg", this, "Redraw()");   

   fCloseButton->Connect("Clicked()", "TMVA::StatDialogBDTReg", this, "Close()");
}

void TMVA::StatDialogBDTReg::UpdateCanvases() 
{
   DrawTree( fItree );
}

void TMVA::StatDialogBDTReg::GetNtrees()
{
   if(!fWfile.EndsWith(".xml") ){
      std::ifstream fin( fWfile );
      if (!fin.good( )) { // file not found --> Error
         std::cout << "*** ERROR: Weight file: " << fWfile << " does not exist" << std::endl;
         return;
      }
   
      TString dummy = "";
      
      // read total number of trees, and check whether requested tree is in range
      Int_t nc = 0;
      while (!dummy.Contains("NTrees")) { 
         fin >> dummy; 
         nc++; 
         if (nc > 200) {
            std::cout << std::endl;
            std::cout << "*** Huge problem: could not locate term \"NTrees\" in BDT weight file: " 
                      << fWfile << std::endl;
            std::cout << "==> panic abort (please contact the TMVA authors)" << std::endl;
            std::cout << std::endl;
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
   std::cout << "--- Found " << fNtrees << " decision trees in weight file" << std::endl;

}

////////////////////////////////////////////////////////////////////////////////
/// recursively puts an entries in the histogram for the node and its daughters
///

void TMVA::StatDialogBDTReg::DrawNode( TMVA::DecisionTreeNode *n, 
                                       Double_t x, Double_t y, 
                                       Double_t xscale,  Double_t yscale, TString * vars) 
{
   Float_t xsize=xscale*1.5;
   Float_t ysize=yscale/3;
   if (xsize>0.15) xsize=0.1;
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
   if      (n->GetNodeType() ==  1) { t->SetFillColor( getSigColorF() ); t->SetTextColor( getSigColorT() ); }
   else if (n->GetNodeType() == -1) { t->SetFillColor( getBkgColorF() ); t->SetTextColor( getBkgColorT() ); }
   else if (n->GetNodeType() ==  0) { t->SetFillColor( getIntColorF() ); t->SetTextColor( getIntColorT() ); }

   char buffer[25];
   //   sprintf( buffer, "N=%f", n->GetNEvents() );
   //   t->AddText(buffer);
   sprintf( buffer, "R=%4.1f +- %4.1f", n->GetResponse(),n->GetRMS() );
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

TMVA::DecisionTree* TMVA::StatDialogBDTReg::ReadTree( TString* &vars, Int_t itree )
{
   std::cout << "--- Reading Tree " << itree << " from weight file: " << fWfile << std::endl;
   TMVA::DecisionTree *d = new TMVA::DecisionTree();


   if(!fWfile.EndsWith(".xml") ){

      std::ifstream fin( fWfile );
      if (!fin.good( )) { // file not found --> Error
         std::cout << "*** ERROR: Weight file: " << fWfile << " does not exist" << std::endl;
         delete d;
         d = nullptr;
         return 0;
      }
      TString dummy = "";
      
      if (itree >= fNtrees) {
         std::cout << "*** ERROR: requested decision tree: " << itree 
                   << ", but number of trained trees only: " << fNtrees << std::endl;
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
      vars = new TString[nVars+1];
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
         std::cout << "*** ERROR: requested decision tree: " << itree 
                   << ", but number of trained trees only: " << fNtrees << std::endl;
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

void TMVA::StatDialogBDTReg::DrawTree( Int_t itree )
{
   TString *vars;   

   TMVA::DecisionTree* d = ReadTree( vars, itree );
   if (d == 0) return;

   UInt_t   depth = d->GetTotalTreeDepth();
   Double_t ystep = 1.0/(depth + 1.0);

   std::cout << "--- Tree depth: " << depth << std::endl;

   TStyle* TMVAStyle   = gROOT->GetStyle("Plain"); // our style is based on Plain
   Int_t   canvasColor = TMVAStyle->GetCanvasColor(); // backup

   TString cbuffer = Form( "Reading weight file: %s", fWfile.Data() );
   TString tbuffer = Form( "Regression Tree no.: %d", itree );
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

   TPaveText *intermediate = new TPaveText(0.02,ydown,0.15,yup, "NDC");
   intermediate->SetBorderSize(1);
   intermediate->SetFillStyle(1001);
   intermediate->SetFillColor( getIntColorF() );
   intermediate->AddText("Intermediate Nodes");
   intermediate->SetTextColor( getIntColorT() );
   intermediate->Draw();

   ydown = ydown - ystep/2.5 -dy;
   yup   = yup - ystep/2.5 -dy;
   TPaveText *signalleaf = new TPaveText(0.02,ydown ,0.15,yup, "NDC");
   signalleaf->SetBorderSize(1);
   signalleaf->SetFillStyle(1001);
   signalleaf->SetFillColor( getSigColorF() );
   signalleaf->AddText("Leaf Nodes");
   signalleaf->SetTextColor( getSigColorT() );
   signalleaf->Draw();
   /*
     ydown = ydown - ystep/2.5 -dy;
     yup   = yup - ystep/2.5 -dy;
     TPaveText *backgroundleaf = new TPaveText(0.02,ydown,0.15,yup, "NDC");
     backgroundleaf->SetBorderSize(1);
     backgroundleaf->SetFillStyle(1001);
     backgroundleaf->SetFillColor( kBkgColorF );

     backgroundleaf->AddText("Backgr. Leaf Nodes");
     backgroundleaf->SetTextColor( kBkgColorT );
     backgroundleaf->Draw();
   */
   fCanvas->Update();
   TString fname = fDataset+Form("/plots/%s_%i", fMethName.Data(), itree );
   std::cout << "--- Creating image: " << fname << std::endl;
   TMVAGlob::imgconv( fCanvas, fname );   

   TMVAStyle->SetCanvasColor( canvasColor );
}   
      
// ========================================================================================

// intermediate GUI
void TMVA::BDT_Reg(TString dataset, const TString& fin )
{
   // --- read the available BDT weight files

   // destroy all open cavases
   TMVAGlob::DestroyCanvases(); 

   // checks if file with name "fin" is already open, and if not opens one
   TFile* file = TMVAGlob::OpenFile( fin );  

   TDirectory* dir = file->GetDirectory(dataset.Data())->GetDirectory( "Method_BDT" );
   if (!dir) {
      std::cout << "*** Error in macro \"BDT_Reg.C\": cannot find directory \"Method_BDT\" in file: " << fin << std::endl;
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
         std::cout << "*** Error in macro \"BDT_Reg.C\": cannot find sub-directory: " << key->GetName() 
                   << " in directory: " << dir->GetName() << std::endl;
         return;
      }

      // retrieve weight file name and path
      TObjString* strPath = (TObjString*)mdir->Get( "TrainingPath" );
      TObjString* strWFile = (TObjString*)mdir->Get( "WeightFileName" );
      if (!strPath || !strWFile) {
         std::cout << "*** Error in macro \"BDT_Reg.C\": could not find TObjStrings \"TrainingPath\" and/or \"WeightFileName\" *** " << std::endl;
         std::cout << "*** Maybe you are using TMVA >= 3.8.15 with an older training target file ? *** " << std::endl;
         return;
      }

      methname.push_back( key->GetName() );
      path    .push_back( strPath->GetString() );
      wfile   .push_back( strWFile->GetString() );
   }

   // create the control bar
   TControlBar* cbar = new TControlBar( "vertical", "Choose weight file:", 50, 50 );
   BDTReg_Global__cbar.push_back(cbar);

   for (UInt_t im=0; im<path.size(); im++) {
      TString fname = path[im];
      if (fname[fname.Length()-1] != '/') fname += "/";
      fname += wfile[im];
      TString macro = Form( "TMVA::BDT_Reg(\"%s\",0,\"%s\",\"%s\")",dataset.Data(), fname.Data(), methname[im].Data() );
      cbar->AddButton( fname, macro, "Plot decision trees from this weight file", "button" );
   }

   // set the style 
   cbar->SetTextColor("blue");

   // draw
   cbar->Show();   
}

void TMVA::BDTReg_DeleteTBar(int i)
{
   // destroy all open canvases
   StatDialogBDTReg::Delete();
   TMVAGlob::DestroyCanvases();

   delete BDTReg_Global__cbar[i];
   BDTReg_Global__cbar[i] = 0;
}

// input: - No. of tree
//        - the weight file from which the tree is read
void TMVA::BDT_Reg(TString dataset, Int_t itree, TString wfile , TString methName, Bool_t useTMVAStyle  ) 
{
   // destroy possibly existing dialog windows and/or canvases
   StatDialogBDTReg::Delete();
   TMVAGlob::DestroyCanvases(); 
   if(wfile=="")
      wfile = dataset+"/weights/TMVARegression_BDT.weights.xml";

   // quick check if weight file exist
   if(!wfile.EndsWith(".xml") ){
      std::ifstream fin( wfile );
      if (!fin.good( )) { // file not found --> Error
         std::cout << "*** ERROR: Weight file: " << wfile << " does not exist" << std::endl;
         return;
      }
   }
   std::cout << "test1";
   // set style and remove existing canvas'
   TMVAGlob::Initialize( useTMVAStyle );

   StatDialogBDTReg* gGui = new StatDialogBDTReg(dataset, gClient->GetRoot(), wfile, methName, itree );

   gGui->DrawTree( itree );

   gGui->RaiseDialog();
}

