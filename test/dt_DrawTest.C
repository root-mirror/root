#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TROOT.h"
#include <TBenchmark.h>

int gHasLibrary = kFALSE;
int gBranchStyle = 1;
TList gSkipped;

void DrawSkippable(TTree* tree, const char* what, const char* where, Bool_t draw = true) {
  //cerr << "Doing " << what << " which is " << skip << endl;
   if (draw) {
      TString cut = what;
      cut.Append(">>");
      cut.Append(where);
      tree->Draw(cut.Data(),"","goff");
      TH1* h = (TH1*)gDirectory->FindObject(where);
      if (h) h->SetTitle(Form("histo made from T->Draw(\"%s\")",what));
   } else {
      gSkipped.Add(new TNamed(where,where));
   }
};

void DrawSkippable(TTree* tree, const char* what, const char* cond,
                   const char* where, Bool_t draw = true) {
   //cerr << "Doing " << what << " which is " << skip << endl;
   if (draw) {
      TString cut = what;
      cut.Append(">>");
      cut.Append(where);
      tree->Draw(cut.Data(),cond,"goff");
      TH1* h = (TH1*)gDirectory->FindObject(where);
      if (h) h->SetTitle(Form("histo made from T->Draw(\"%s\",\"%s\")",what,cond));
   } else {
      gSkipped.Add(new TNamed(where,where));
   }
};

// Rootmarks for fcdflnx1 is 153.4
void DrawMarks() {

  // The base is currently: RunDrawTest.C++("Event.old.split.root",0)
  Float_t rt_base = 2.33/10.0;
  Float_t cp_base = 2.34/10.0;

  Float_t rt = gBenchmark->GetRealTime("DrawTest");
  Float_t ct = gBenchmark->GetCpuTime("DrawTest");

  // gBenchmark->Print("DrawTest");

  Float_t rootmarks = 200*(rt_base + cp_base)/(rt + ct);

  printf("*  ROOTMARKS =%6.1f   *  Root%-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),gROOT->GetVersionDate(),gROOT->GetVersionTime());

}


//_______________________________________________________________
TDirectory* GenerateDrawHist(TTree *tree, int quietLevel = 0, int level = 3)
{
// Test selections via TreeFormula
// tree is a TTree when called by stress9
// tree is a TChain when called from stres11
// This is a quite complex test checking the results of TTree::Draw
// or TChain::Draw with an explicit loop on events.
// Also a good test for the interpreter

   gROOT->cd();
   TDirectory *hfile = gDirectory;

   gBenchmark = new TBenchmark();
   gBenchmark->Start("DrawTest");

   // Each tree->Draw generates an histogram
   DrawSkippable(tree,"GetNtrack()","hGetNtrack",(level>0 && gHasLibrary));

   //gBenchmark->Show("DrawTest");  gBenchmark->Start("DrawTest");

   DrawSkippable(tree,"fNtrack","hNtrack");
   DrawSkippable(tree,"fNseg","hNseg");
   DrawSkippable(tree,"fTemperature","hTemp");

   DrawSkippable(tree,"fH.GetMean()","hHmean");
   if (level>0) DrawSkippable(tree,"fH.fXaxis.fXmax","hHAxisMax");
   if (level>0) DrawSkippable(tree,"fH.fXaxis.GetXmax()","hHAxisGetMax");
   DrawSkippable(tree,"fH.GetXaxis()->GetXmax()","hHGetAxisGetMax",(level>0));
   DrawSkippable(tree,"fH.GetXaxis()->fXmax","hHGetAxisMax",(level>0));
   DrawSkippable(tree,"GetHistogram()->GetXaxis()->GetXmax()","hGetHGetAxisMax",
                 (level>0&&gHasLibrary));
   DrawSkippable(tree,"event.GetHistogram()->GetXaxis()->GetXmax()",
                 "hGetRefHGetAxisMax",(level>0&&gHasLibrary));

   DrawSkippable(tree,"fTracks.fPx","fEvtHdr.fEvtNum%10 == 0","hPx");
   DrawSkippable(tree,"fTracks.fPy","fEvtHdr.fEvtNum%10 == 0","hPy");
   DrawSkippable(tree,"fTracks.fPz","fEvtHdr.fEvtNum%10 == 0","hPz");
   DrawSkippable(tree,"fRandom","3*(fEvtHdr.fEvtNum%10 == 1)","hRandom");
   DrawSkippable(tree,"fMass2",  "fEvtHdr.fEvtNum%10 == 1","hMass2");
   DrawSkippable(tree,"fBx",        "fEvtHdr.fEvtNum%10 == 1","hBx");
   DrawSkippable(tree,"fBy",        "fEvtHdr.fEvtNum%10 == 1","hBy");
   DrawSkippable(tree,"fXfirst","fEvtHdr.fEvtNum%10 == 2","hXfirst");
   DrawSkippable(tree,"fYfirst","fEvtHdr.fEvtNum%10 == 2","hYfirst");
   DrawSkippable(tree,"fZfirst","fEvtHdr.fEvtNum%10 == 2","hZfirst");
   DrawSkippable(tree,"fXlast",  "fEvtHdr.fEvtNum%10 == 3","hXlast");
   DrawSkippable(tree,"fYlast",  "fEvtHdr.fEvtNum%10 == 3","hYlast");
   DrawSkippable(tree,"fZlast",  "fEvtHdr.fEvtNum%10 == 3","hZlast");
   DrawSkippable(tree,"fCharge","fPx < 0","hCharge");
   DrawSkippable(tree,"fNpoint","fPx < 0","hNpoint");
   DrawSkippable(tree,"fValid",  "fPx < 0","hValid");
   DrawSkippable(tree,"fPointValue","hPointValue", gBranchStyle!=0);
   tree->SetAlias("mult","fPx*fPy");
   DrawSkippable(tree,"fEvtHdr.fEvtNum*6+mult", "hAlias", 1);
   tree->SetAlias("track","event.fTracks");
   DrawSkippable(tree,"track.fPx+track.fPy", "hAliasSymbol", 1);
   DrawSkippable(tree,"track.GetPx()+track.GetPy()","hAliasSymbolFunc", gBranchStyle!=0 && gHasLibrary);

   DrawSkippable(tree,"fIsValid","hBool");

   DrawSkippable(tree,"fMatrix","hFullMatrix");
   DrawSkippable(tree,"fMatrix[][0]","hColMatrix");
   DrawSkippable(tree,"fMatrix[1][]","hRowMatrix");
   DrawSkippable(tree,"fMatrix[2][2]","hCellMatrix");

   DrawSkippable(tree,"fMatrix - fVertex","hFullOper");
   DrawSkippable(tree,"fMatrix[2][1] - fVertex[5][1]","hCellOper");
   DrawSkippable(tree,"fMatrix[][1]  - fVertex[5][1]","hColOper");
   DrawSkippable(tree,"fMatrix[2][]  - fVertex[5][2]","hRowOper");
   DrawSkippable(tree,"fMatrix[2][]  - fVertex[5][]","hMatchRowOper");
   DrawSkippable(tree,"fMatrix[][2]  - fVertex[][1]","hMatchColOper");
   DrawSkippable(tree,"fMatrix[][2]  - fVertex[][]","hRowMatOper");
   DrawSkippable(tree,"fMatrix[][2]  - fVertex[5][]","hMatchDiffOper");
   DrawSkippable(tree,"fMatrix[][]   - fVertex[][]","hFullOper2");

   // Test on variable arrays
   DrawSkippable(tree,"fClosestDistance","hClosestDistance",gBranchStyle!=0);
   DrawSkippable(tree,"fClosestDistance[2]","hClosestDistance2",gBranchStyle!=0);
   DrawSkippable(tree,"fClosestDistance[9]","hClosestDistance9",gBranchStyle!=0);

   // Test variable indexing
   DrawSkippable(tree,"fClosestDistance[fNvertex/2]","hClosestDistanceIndex",
                 (level>0)&&gBranchStyle!=0);
   DrawSkippable(tree,"fPx:fPy[fNpoint/6]","fPy[fNpoint/6]>0","hPxInd",(level>0));

   // Test of simple function calls
   DrawSkippable(tree,"sqrt(fNtrack)","hSqrtNtrack",(level>0));

   // Test string operations
   DrawSkippable(tree,"fEvtHdr.fEvtNum","fType==\"type1\" ","hString",(level>0));
   DrawSkippable(tree,"fEvtHdr.fEvtNum","1 && strstr(fType,\"1\") ","+hString",(level>0));
   tree->SetAlias("typ","fType");
   DrawSkippable(tree,"strstr(typ,\"1\") ", "hAliasStr", 1);
   DrawSkippable(tree,"fH.fTitle.fData==\"Event Histogram\"","hStringSpace",1);

   // Test binary operators
   DrawSkippable(tree,"fValid<<4","hShiftValid",(level>0));
   DrawSkippable(tree,"((fValid<<4)>>2)","+hShiftValid",(level>0));
   DrawSkippable(tree,"fValid&0x1","(fNvertex>10) && (fNseg<=6000)"
                 ,"hAndValid",(level>0));

   // Test weight
   DrawSkippable(tree,"fPx","(fBx>.15) || (fBy<=-.15)","hPxBx",(level>0));
   DrawSkippable(tree,"fPx","fBx*fBx*(fBx>.15) + fBy*fBy*(fBy<=-.15)",
                 "hPxBxWeight",(level>0));

   DrawSkippable(tree,"event.fTriggerBits",
                 "hTriggerBits",level>1 && gBranchStyle!=0);
   DrawSkippable(tree,"event.fTriggerBits.fNbits",
                 "event.fTriggerBits.TestBitNumber(28)",
                 "hFiltTriggerBits",level>1 && gBranchStyle!=0);

   DrawSkippable(tree,"event.GetTriggerBits()",
                 "hTriggerBitsFunc",level>1 && gBranchStyle!=0 && gHasLibrary );

   DrawSkippable(tree,"fTracks.fTriggerBits",
                 "hTrackTrigger", level>1 && gBranchStyle!=0);
   DrawSkippable(tree,"fPx",
                 "fTracks.fTriggerBits.TestBitNumber(5)",
                 "hFiltTrackTrigger", level>1 && gBranchStyle!=0);

   DrawSkippable(tree,"TMath::BreitWigner(fPx,3,2)",
                 "",
                 "hBreit", level>1 && gBranchStyle!=0);

   // Test on alternate value
   DrawSkippable(tree,"fMatrix-Alt$(fClosestDistance,0)",
                 "",
                 "hAlt", level>1 && gBranchStyle!=0);

   // Test on the @ notation to access the collection object itself
   DrawSkippable(tree,"event.@fTracks.size()","","hSize",level>2 && !(gBranchStyle==0 && !gHasLibrary));
   DrawSkippable(tree,"event.fTracks@.size()","","+hSize",level>2 && !(gBranchStyle==0 && !gHasLibrary));
   DrawSkippable(tree,"@fTracks.size()","","hSize2",level>2 && gBranchStyle!=0);
   DrawSkippable(tree,"fTracks@.size()","","+hSize2",level>2 && gBranchStyle!=0);
   DrawSkippable(tree,"Sum$(fPx)","","hSumPx",level>2);
   DrawSkippable(tree,"MaxIf$(fPx,fPy>1.0):Max$(fPy)","Sum$(fPy>1.0)>0","hMaxPx",level>2);
   DrawSkippable(tree,"MinIf$(fPx,fPy>1.0):Min$(fPy)","Sum$(fPy>1.0)>0","hMinPx",level>2);

   if (quietLevel<2) gBenchmark->Show("DrawTest");
   else gBenchmark->Stop("DrawTest");
   gBenchmark->Start("DrawTest");

   return hfile;

}

