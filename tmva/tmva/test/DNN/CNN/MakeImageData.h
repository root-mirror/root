#include "TH2.h"
#include "TF2.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include <vector>

void makeImages(int n = 10000, int npx = 8, int npy = 8, TString fileName = "imagesData.root") {

   int np = npx*npy;
   auto h1 = new TH2D("h1","h1",npx,0,10,npy,0,10);
   auto h2 = new TH2D("h2","h2",npx,0,10,npy,0,10);

   auto f1 = new TF2("f1","xygaus");
   auto f2 = new TF2("f2","xygaus");
   std::vector<float> x1(np);
   std::vector<float> x2(np);
   TTree sgn("sgn","sgn");
   TTree bkg("bkg","bkg");
   TFile f(fileName,"RECREATE");

   bkg.Branch("ximage",x1.data(),Form("ximage[%d]/F",np));
   sgn.Branch("ximage",x2.data(),Form("ximage[%d]/F",np));

   sgn.SetDirectory(&f);
   bkg.SetDirectory(&f);

   f1->SetParameters(1,5,3,5,3);
   f2->SetParameters(1,5,3.9,5,2.1 );
   gRandom->SetSeed(0);
   for (int i = 0; i < n; ++i) {
      h1->Reset();
      h2->Reset();
      //generate random means in range [3,7] to be not too much on the border
      f1->SetParameter(1,gRandom->Uniform(3,7));
      f1->SetParameter(3,gRandom->Uniform(3,7));
      f2->SetParameter(1,gRandom->Uniform(3,7));
      f2->SetParameter(3,gRandom->Uniform(3,7));

      h1->FillRandom("f1",np*100);
      h2->FillRandom("f2",np*100);


      for (int k = 0; k < npx ; ++k) {
         for (int l = 0; l < npy ; ++l)  {
            int m = k * npy + l;
            x1[m] = h1->GetBinContent(k+1,l+1) + gRandom->Gaus(0,3);
            x2[m] = h2->GetBinContent(k+1,l+1) + gRandom->Gaus(0,3);
         }
      }
      sgn.Fill();
      bkg.Fill();

//       if (n==1) {
//          h1->Draw("COLZ");
//          new TCanvas();
//          h2->Draw("COLZ");
//       }
   }
   sgn.Print();
   bkg.Print();
   sgn.Write();
   bkg.Write();
   f.Close();
}
