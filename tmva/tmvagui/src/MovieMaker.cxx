#include "TMVA/MovieMaker.h"
#include "TString.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TKey.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TText.h"
#include "TObjString.h"

#include <vector>
#include "TMVA/network.h"

void TMVA::DrawNetworkMovie(TString dataset, TFile* file, const TString& methodType, const TString& methodTitle )
{

   TString     dirname  = methodType + "/" + methodTitle + "/" + "EpochMonitoring";
   TDirectory *epochDir = (TDirectory*)file->GetDirectory(dataset.Data())->Get( dirname );
   if (!epochDir) {
      cout << "Big troubles: could not find directory \"" << dirname << "\"" << endl;
      exit(1);
   }
   epochDir->cd();

   // loop over all epoch-wise monitoring histograms
   TIter keyIt(epochDir->GetListOfKeys());
   TKey *key;
   std::vector<TString> epochList;
   Int_t ic = 0;
   while ((key = (TKey*)keyIt())) {
      
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH2F")) continue;
      TString name = key->GetName();
      
      if (!name.BeginsWith("epochmonitoring___")) continue;
      
      // extract epoch
      TObjArray* tokens = name.Tokenize("_");
      TString es = ((TObjString*)tokens->At(2))->GetString();

      // check if done already      
      Bool_t isOld = kFALSE;
      for (std::vector<TString>::const_iterator it = epochList.begin(); it < epochList.end(); ++it) {
         if (*it == es) isOld = kTRUE; 
      }
      if (isOld) continue;
      epochList.push_back( es );

      // create bulk file name
      TString bulkname  = Form( "epochmonitoring___epoch_%s_weights_hist", es.Data() );

      // draw the network
      if (ic <= 60) draw_network(dataset, file, epochDir, bulkname, kTRUE, es );
      ic++;
   }
}


void TMVA::DrawMLPoutputMovie(TString dataset, TFile* file, const TString& methodType, const TString& methodTitle )
{
   gROOT->SetBatch( 1 );

   // define Canvas layout here!
   const Int_t width = 600;   // size of canvas

   // this defines how many canvases we need
   TCanvas* c = 0;

   Float_t nrms = 4;
   Float_t xmin = -1.2;
   Float_t xmax = 1.2;
   Float_t ymin = 0;
   Float_t ymax = 0;
   Float_t maxMult = 6.0;
   Int_t   countCanvas = 0;
   Bool_t  first = kTRUE;
            
   TString     dirname  = methodType + "/" + methodTitle + "/" + "EpochMonitoring";
   TDirectory *epochDir = (TDirectory*)file->GetDirectory(dataset.Data())->Get( dirname );
   if (!epochDir) {
      cout << "Big troubles: could not find directory \"" << dirname << "\"" << endl;
      exit(1);
   }

   // now read all evolution histograms
   TIter keyItTit(epochDir->GetListOfKeys());
   TKey *titkeyTit;
   while ((titkeyTit = (TKey*)keyItTit())) {
      
      if (!gROOT->GetClass(titkeyTit->GetClassName())->InheritsFrom("TH1F")) continue;
      TString name = titkeyTit->GetName();
      
      if (!name.BeginsWith("convergencetest___")) continue;
      if (!name.Contains("_train_"))              continue; // only for training so far
      if (name.EndsWith( "_B"))                   continue;
      
      // must be signal histogram
      if (!name.EndsWith( "_S")) {
         cout << "Big troubles with histogram: " << name << " -> should end with _S" << endl;
         exit(1);
      }
      
      // create canvas
      countCanvas++;
      TString ctitle = Form("TMVA response %s",methodTitle.Data());
      c = new TCanvas( Form("canvas%d", countCanvas), ctitle, 0, 0, width, (Int_t)width*0.78 ); 
      
      TH1F* sig = (TH1F*)titkeyTit->ReadObj();
      sig->SetTitle( Form("TMVA response for classifier: %s", methodTitle.Data()) );
      
      TString dataType = (name.Contains("_train_") ? "(training sample)" : "(test sample)");
      
      // find background
      TString nbn = sig->GetName(); nbn[nbn.Length()-1] = 'B';            
      TH1F* bgd = dynamic_cast<TH1F*>(epochDir->Get( nbn ));
      if (bgd == 0) {
         cout << "Big troubles with histogram: " << bgd << " -> cannot find!" << endl;
         exit(1);
      }
      
      cout << "sig = " << sig->GetName() << endl;
      cout << "bgd = " << bgd->GetName() << endl;
      
      // set the histogram style
      TMVAGlob::SetSignalAndBackgroundStyle( sig, bgd );
      
      // normalise both signal and background
      TMVAGlob::NormalizeHists( sig, bgd );
      
      // set only first time, then same for all plots
      if (first) {
         if (xmin == 0 && xmax == 0) {
            xmin = TMath::Max( TMath::Min(sig->GetMean() - nrms*sig->GetRMS(), 
                                          bgd->GetMean() - nrms*bgd->GetRMS() ),
                               sig->GetXaxis()->GetXmin() );
            xmax = TMath::Min( TMath::Max(sig->GetMean() + nrms*sig->GetRMS(), 
                                          bgd->GetMean() + nrms*bgd->GetRMS() ),
                               sig->GetXaxis()->GetXmax() );
         }
         ymin = 0;
         ymax = TMath::Max( sig->GetMaximum(), bgd->GetMaximum() )*maxMult;
         first = kFALSE;
      }
      
      // build a frame
      Int_t nb = 100;
      TString hFrameName(TString("frame") + methodTitle);
      TObject *o = gROOT->FindObject(hFrameName);
      if(o) delete o;
      TH2F* frame = new TH2F( hFrameName, sig->GetTitle(), 
                              nb, xmin, xmax, nb, ymin, ymax );
      frame->GetXaxis()->SetTitle( methodTitle + " response" );
      frame->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
      TMVAGlob::SetFrameStyle( frame );
      
      // find epoch number (4th token)
      TObjArray* tokens = name.Tokenize("_");
      TString es = ((TObjString*)tokens->At(4))->GetString();
      if (!es.IsFloat()) {
         cout << "Big troubles in epoch parsing: \"" << es << "\" is not float" << endl;
         exit(1);
      }
      Int_t epoch = es.Atoi();
      
      // eventually: draw the frame
      frame->Draw();  
      
      c->GetPad(0)->SetLeftMargin( 0.105 );
      frame->GetYaxis()->SetTitleOffset( 1.2 );
      
      // Draw legend               
      TLegend *legend= new TLegend( c->GetLeftMargin(), 1 - c->GetTopMargin() - 0.12, 
                                    c->GetLeftMargin() + 0.5, 1 - c->GetTopMargin() );
      legend->SetFillStyle( 1 );
      legend->AddEntry(sig,TString("Signal ")     + dataType, "F");
      legend->AddEntry(bgd,TString("Background ") + dataType, "F");
      legend->SetBorderSize(1);
      legend->SetMargin( 0.15 );
      legend->Draw("same");
      
      TText* t = new TText();            
      t->SetTextSize( 0.04 );
      t->SetTextColor( 1 );
      t->SetTextAlign( 31 );
      t->DrawTextNDC( 1 - c->GetRightMargin(), 1 - c->GetTopMargin() + 0.015, Form( "Epoch: %i", epoch) );
      
      // overlay signal and background histograms
      sig->Draw("samehist");
      bgd->Draw("samehist");
      
      // save to file
      TString outdirname  = "movieplots";
      TString foutname = outdirname + "/" + name;
      foutname.Resize( foutname.Length()-2 );
      foutname.ReplaceAll("convergencetest___","");
      foutname += ".gif";
      
      cout << "storing file: " << foutname << endl;
      
      c->Update();
      c->Print(foutname);            
   }
}

// -----------------------------------------------------------------------------

void TMVA::MovieMaker(TString dataset, TString methodType , TString methodTitle )
{
   TString fname = "TMVA.root";
   TFile* file = TMVAGlob::OpenFile( fname );     

   //DrawMLPoutputMovie( file, methodType, methodTitle );
   DrawNetworkMovie(dataset, file, methodType, methodTitle );
}   

