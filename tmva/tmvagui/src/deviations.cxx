#include "TMVA/deviations.h"
#include "TMVA/Types.h"
#include "TH2.h"
#include "TLine.h"
#include "TObjString.h"


// this macro plots the resulting MVA distributions (Signal and
// Background overlayed) of different MVA methods run in TMVA
// (e.g. running TMVAnalysis.C).


// input: - Input file (result from TMVA)
//        - use of TMVA plotting TStyle
void TMVA::deviations(TString dataset, TString fin, 
                      HistType htype , Bool_t showTarget, Bool_t useTMVAStyle  )
{
   // set style and remove existing canvas'
   TMVAGlob::Initialize( useTMVAStyle );
   gStyle->SetNumberContours(999);

   // checks if file with name "fin" is already open, and if not opens one
   TFile* file = TMVAGlob::OpenFile( fin );  

   // define Canvas layout here!
   const Int_t width = 650;   // size of canvas

   // this defines how many canvases we need
   TCanvas* c[100];

   // counter variables
   Int_t countCanvas = 0;

   // search for the right histograms in full list of keys
   //    TList* methods = new TMap();

   TIter next(file->GetDirectory(dataset.Data())->GetListOfKeys());
   TKey *key(0);   
   while ((key = (TKey*)next())) {

      if (!TString(key->GetName()).BeginsWith("Method_")) continue;
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) continue;

      TString methodName;
      TMVAGlob::GetMethodName(methodName,key);
      cout << "--- Plotting deviation for method: " << methodName << endl;

      TDirectory* mDir = (TDirectory*)key->ReadObj();

      TList* jobNames = new TList();

      TIter keyIt(mDir->GetListOfKeys());
      TKey *titkey;
      while ((titkey = (TKey*)keyIt())) {

         if (!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory")) continue;

         TDirectory *titDir = (TDirectory *)titkey->ReadObj();

         TObjString *jN = new TObjString( titDir->GetName() );
         if (!jobNames->Contains( jN )) jobNames->Add( jN );
         else delete jN;
            
         TString methodTitle;
         TMVAGlob::GetMethodTitle(methodTitle,titDir);

         TString hname = "MVA_" + methodTitle;
         TIter   dirKeyIt( titDir->GetListOfKeys() );
         TKey*   dirKey;

         Int_t countPlots = 0;
         while ((dirKey = (TKey*)dirKeyIt())){
            if (dirKey->ReadObj()->InheritsFrom("TH2F")) {
               TString s(dirKey->ReadObj()->GetName());
               if (s.Contains("_reg_") && 
                   ( (showTarget && s.Contains("_tgt")) || (!showTarget && !s.Contains("_tgt")) ) && 
                   s.Contains( (htype == kCompareType ? "train" : "test" ))) {
                  c[countCanvas] = new TCanvas( Form("canvas%d", countCanvas+1), 
                                                Form( "Regression output deviation versus %s for method: %s",
                                                      (showTarget ? "target" : "input variables"), methodName.Data() ),
                                                countCanvas*50+100, (countCanvas+1)*20, width, (Int_t)width*0.72 ); 
                  c[countCanvas]->SetRightMargin(0.10); // leave space for border
                  TH1* h = (TH1*)dirKey->ReadObj();
                  h->SetTitle( Form("Output deviation for method: %s (%s sample)", 
                                    hname.Data(), (htype == kCompareType ? "training" : "test" )) );
                  //                                    methodName.Data(), (htype == kCompareType ? "training" : "test" )) );
                  h->Draw("colz");
                  TLine* l = new TLine( h->GetXaxis()->GetXmin(), 0, h->GetXaxis()->GetXmax(), 0 );
                  l->SetLineStyle(2);
                  l->Draw();

                  // update and print
                  cout << "plotting logo" << endl;
                  TMVAGlob::plot_logo(1.058);
                  c[countCanvas]->Update();

                  TString fname = Form( "%s/plots/deviation_%s_%s_%s_c%i", 
                                        dataset.Data(),
                                        methodName.Data(), 
                                        (showTarget ? "target" : "vars"),
                                        (htype == kCompareType ? "training" : "test" ), countPlots );
                  TMVAGlob::imgconv( c[countCanvas], fname );

                  countPlots++;
                  countCanvas++;
               }
            }
         }         
      }
   }
}




