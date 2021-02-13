#include "TMVA/mvasMulticlass.h"
#include "TMVA/Types.h"
#include "TLegend.h"
#include "TText.h"
#include "TH2.h"



// this macro plots the resulting MVA distributions (Signal and
// Background overlayed) of different MVA methods run in TMVA
// (e.g. running TMVAnalysis.C).


// input: - Input file (result from TMVA)
//        - use of TMVA plotting TStyle
void TMVA::mvasMulticlass(TString dataset, TString fin , HistType htype , Bool_t useTMVAStyle  )
{
   // set style and remove existing canvas'
   TMVAGlob::Initialize( useTMVAStyle );

   // switches
   const Bool_t Save_Images = kTRUE;

   // checks if file with name "fin" is already open, and if not opens one
   TFile *file = TMVAGlob::OpenFile(fin);

   // find directory from input variable transformation and extract class names
   TIter keys_dir = file->GetDirectory(dataset.Data())->GetListOfKeys();
   TKey *key;
   TDirectory *tempdir = 0;
   while ((key = (TKey *)keys_dir())) {
      TString name = key->GetName();
      if (name.BeginsWith("InputVariables_")) {
         tempdir = (TDirectory *)(file->GetDirectory(dataset.Data())->Get(name));
         break;
      }
   }
   std::vector<TString> classnames(TMVAGlob::GetClassNames(tempdir));

   // define Canvas layout here!
   // Int_t xPad = 1; // no of plots in x
   // Int_t yPad = 1; // no of plots in y
   //   Int_t noPad = xPad * yPad ;
   const Int_t width = 600;   // size of canvas

   // this defines how many canvases we need
   TCanvas *c = 0;

   // counter variables
   Int_t countCanvas = 0;

   // search for the right histograms in full list of keys
   TIter next(file->GetDirectory(dataset.Data())->GetListOfKeys());
   key = 0;
   while ((key = (TKey*)next())) {

      if (!TString(key->GetName()).BeginsWith("Method_")) continue;
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) continue;

      TString methodName;
      TMVAGlob::GetMethodName(methodName,key);

      TDirectory* mDir = (TDirectory*)key->ReadObj();

      TIter keyIt(mDir->GetListOfKeys());
      TKey *titkey;
      while ((titkey = (TKey*)keyIt())) {

         if (!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory")) continue;

         TDirectory *titDir = (TDirectory *)titkey->ReadObj();
         TString methodTitle;
         TMVAGlob::GetMethodTitle(methodTitle,titDir);

         cout << "--- Found directory for method: " << methodName << "::" << methodTitle << endl;
         TString hname = "MVA_" + methodTitle;
         for(UInt_t icls = 0; icls < classnames.size(); ++icls){
            TObjArray hists;

            std::vector<TString>::iterator classiter = classnames.begin();
            for(; classiter!=classnames.end(); ++classiter){
               TString name(hname+"_Test_"+ classnames.at(icls)
                            + "_prob_for_" + *classiter);
               TH1 *hist = (TH1*)titDir->Get(name);
               if (hist==0){
                  cout << ":\t mva distribution not available (this is normal for Cut classifier)" << endl;
                  continue;
               }
               hists.Add(hist);
            }


            // chop off useless stuff
            ((TH1*)hists.First())->SetTitle( Form("TMVA response for classifier: %s", methodTitle.Data() ));

            // create new canvas
            //cout << "Create canvas..." << endl;
            TString ctitle = ((htype == kMVAType) ?
                              Form("TMVA response for class %s %s", classnames.at(icls).Data(),methodTitle.Data()) :
                              Form("TMVA comparison for class %s %s", classnames.at(icls).Data(),methodTitle.Data())) ;

            c = new TCanvas( Form("canvas%d", countCanvas+1), ctitle,
                             countCanvas*50+200, countCanvas*20, width, (Int_t)width*0.78 );

            // set the histogram style
            //cout << "Set histogram style..." << endl;
            TMVAGlob::SetMultiClassStyle( &hists );

            // normalise all histograms and find maximum
            Float_t histmax = -1;
            for(Int_t i=0; i<hists.GetEntriesFast(); ++i){
               TMVAGlob::NormalizeHist((TH1*)hists[i] );
               if(((TH1*)hists[i])->GetMaximum() > histmax)
                  histmax = ((TH1*)hists[i])->GetMaximum();
            }

            // frame limits (between 0 and 1 per definition)
            Float_t xmin = 0;
            Float_t xmax = 1;
            Float_t ymin = 0;
            Float_t maxMult = (htype == kCompareType) ? 1.3 : 1.2;
            Float_t ymax = histmax*maxMult;
            // build a frame
            Int_t nb = 500;
            TString hFrameName(TString::Format("frame_class%d_",icls) + methodTitle);
            TObject *o = gROOT->FindObject(hFrameName);
            if(o) delete o;
            TH1F* frame = new TH1F( hFrameName, ((TH1*)hists.First())->GetTitle(),
                                    nb, xmin, xmax);
            frame->SetMinimum(ymin);
            frame->SetMaximum(ymax);
            frame->GetXaxis()->SetTitle(methodTitle + " response for " + classnames.at(icls));
            frame->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
            TMVAGlob::SetFrameStyle( frame );

            // eventually: draw the frame
            frame->Draw();

            c->GetPad(0)->SetLeftMargin( 0.105 );
            frame->GetYaxis()->SetTitleOffset( 1.2 );

            // Draw legend
            TLegend *legend= new TLegend( c->GetLeftMargin(), 1 - c->GetTopMargin() - 0.12,
                                          c->GetLeftMargin() + (htype == kCompareType ? 0.40 : 0.3), 1 - c->GetTopMargin() );
            legend->SetFillStyle( 1 );
            classiter = classnames.begin();

            for(Int_t i=0; i<hists.GetEntriesFast(); ++i, ++classiter){
               legend->AddEntry(((TH1*)hists[i]),*classiter,"F");
            }

            legend->SetBorderSize(1);
            legend->SetMargin( 0.3 );
            legend->Draw("same");


            for(Int_t i=0; i<hists.GetEntriesFast(); ++i){

               ((TH1*)hists[i])->Draw("histsame");
               TString ytit = TString("(1/N) ") + ((TH1*)hists[i])->GetYaxis()->GetTitle();
               ((TH1*)hists[i])->GetYaxis()->SetTitle( ytit ); // histograms are normalised

            }


            if (htype == kCompareType) {

               TObjArray othists;
               // if overtraining check, load additional histograms
               classiter = classnames.begin();
               for(; classiter!=classnames.end(); ++classiter){
                  TString name(hname+"_Train_"+ classnames.at(icls)
                               + "_prob_for_" + *classiter);
                  TH1 *hist = (TH1*)titDir->Get(name);
                  if (hist==0){
                     cout << ":\t comparison histogram for overtraining check not available!" << endl;
                     continue;
                  }
                  othists.Add(hist);
               }

               TLegend *legend2= new TLegend( 1 - c->GetRightMargin() - 0.42, 1 - c->GetTopMargin() - 0.12,
                                              1 - c->GetRightMargin(), 1 - c->GetTopMargin() );
               legend2->SetFillStyle( 1 );
               legend2->SetBorderSize(1);

               classiter = classnames.begin();
               for(Int_t i=0; i<othists.GetEntriesFast(); ++i, ++classiter){
                  legend2->AddEntry(((TH1*)othists[i]),*classiter+" (training sample)","P");
               }
               legend2->SetMargin( 0.1 );
               legend2->Draw("same");

               // normalise all histograms and get maximum
               for(Int_t i=0; i<othists.GetEntriesFast(); ++i){
                  TMVAGlob::NormalizeHist((TH1*)othists[i] );
                  if(((TH1*)othists[i])->GetMaximum() > histmax)
                     histmax = ((TH1*)othists[i])->GetMaximum();
               }

               TMVAGlob::SetMultiClassStyle( &othists );
               for(Int_t i=0; i<othists.GetEntriesFast(); ++i){
                  Int_t col = ((TH1*)hists[i])->GetLineColor();
                  ((TH1*)othists[i])->SetMarkerSize( 0.7 );
                  ((TH1*)othists[i])->SetMarkerStyle( 20 );
                  ((TH1*)othists[i])->SetMarkerColor( col );
                  ((TH1*)othists[i])->SetLineWidth( 1 );
                  ((TH1*)othists[i])->Draw("e1same");
               }

               ymax = histmax*maxMult;
               frame->SetMaximum( ymax );

               // for better visibility, plot thinner lines
               TMVAGlob::SetMultiClassStyle( &othists );
               for(Int_t i=0; i<hists.GetEntriesFast(); ++i){
                  ((TH1*)hists[i])->SetLineWidth( 1 );
               }


               // perform K-S test

               cout << "--- Perform Kolmogorov-Smirnov tests" << endl;
               cout << "--- Goodness of consistency for class " << classnames.at(icls)<< endl;
               //TString probatext("Kolmogorov-Smirnov test: ");
               for(Int_t j=0; j<othists.GetEntriesFast(); ++j){
                  Float_t kol = ((TH1*)hists[j])->KolmogorovTest(((TH1*)othists[j]),"X");
                  cout <<  classnames.at(j) << ": " << kol  << endl;
                  //probatext.Append(classnames.at(j)+Form(" %.3f ",kol));
               }



               //TText* tt = new TText( 0.12, 0.74, probatext );
               //tt->SetNDC(); tt->SetTextSize( 0.032 ); tt->AppendPad();

            }


            // redraw axes
            frame->Draw("sameaxis");

            // text for overflows
            //Int_t    nbin = sig->GetNbinsX();
            //Double_t dxu  = sig->GetBinWidth(0);
            //Double_t dxo  = sig->GetBinWidth(nbin+1);
            //TString uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%",
            //                       sig->GetBinContent(0)*dxu*100, bgd->GetBinContent(0)*dxu*100,
            //                      sig->GetBinContent(nbin+1)*dxo*100, bgd->GetBinContent(nbin+1)*dxo*100 );
            //TText* t = new TText( 0.975, 0.115, uoflow );
            //t->SetNDC();
            //t->SetTextSize( 0.030 );
            //t->SetTextAngle( 90 );
            //t->AppendPad();

            // update canvas
            c->Update();

            // save canvas to file

            TMVAGlob::plot_logo(1.058);
            if (Save_Images) {
               if      (htype == kMVAType)     TMVAGlob::imgconv( c, Form("%s/plots/mva_%s_%s",dataset.Data(),classnames.at(icls).Data(), methodTitle.Data()) );
               else if      (htype == kCompareType)     TMVAGlob::imgconv( c, Form("%s/plots/overtrain_%s_%s",dataset.Data(),classnames.at(icls).Data(), methodTitle.Data()) );

            }
            countCanvas++;
         }
      }
      cout << "";
   }
}
