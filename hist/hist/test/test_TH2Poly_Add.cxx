// test TH2Poly adding two histograms

#include "gtest/gtest.h"

#include "TH2Poly.h"
#include "TRandom3.h"

TH2Poly *createPoly(Double_t weight = 1) {
   TH2Poly *h2p = new TH2Poly();
   Int_t i,j;
   Int_t nx = 40;
   Int_t ny = 40;
   Double_t xval1,yval1,xval2,yval2;
   Double_t dx=0.2, dy=0.1;
   xval1 = 0.;
   xval2 = dx;
   for (i = 0; i<nx; i++) {
      yval1 = 0.;
      yval2 = dy;
      for (j = 0; j<ny; j++) {
         h2p->AddBin(xval1, yval1, xval2, yval2);
         yval1 = yval2;
         yval2 = yval2+yval2*dy;
      }
      xval1 = xval2;
      xval2 = xval2+xval2*dx;
   }
   TRandom ran;
   for (i = 0; i<30000; i++) {
      h2p->Fill(50*ran.Gaus(2.,1), ran.Gaus(2.,1), weight);
   }

   return h2p;
}

TEST(TH2Poly, Add)
{
   // Create first hist
   TH2Poly *h2p_1 = createPoly(1);

   // Create second hist
   TH2Poly *h2p_2 = createPoly(0.1);

   // Create an added hist
   TH2Poly *h2p_added = (TH2Poly *)(h2p_1->Clone());
   h2p_added->Add(h2p_2, 1);

   EXPECT_LE( 400, h2p_1->GetMaximum());
   EXPECT_LE(  40, h2p_2->GetMaximum());
   EXPECT_LE( 440, h2p_added->GetMaximum());
}
