/// \file
/// \ingroup tutorial_graphics
/// \notebook
/// A “Perceptual” colormap explicitly identifies a fixed value in the data
///
/// On geographical plot this fixed point can, for instance, the "sea level". A perceptual
/// colormap provides a monotonic luminance variations above and below this fixed value.
/// Unlike the rainbow colormap, this colormap provides  a faithful representation of the
/// structures in the data.
///
/// This macro demonstrates how to produce the perceptual colormap shown on the figure 2
/// in [this article](https://root.cern/blog/rainbow-color-map/).
///
/// The function `Perceptual_Colormap` takes two parameters as input:
///   1. `h`, the `TH2D` to be drawn
///   2. `val_cut`, the Z value defining the "sea level"
///
/// Having these parameters this function defines two color maps: one above `val_cut` and one
/// below.
///
/// \macro_image
/// \macro_code
///
/// \author Olivier Couet

void Perceptual_Colormap(TH2D *h, Double_t val_cut) {
   Double_t max     = h->GetMaximum();         // Histogram's maximum
   Double_t min     = h->GetMinimum();         // Histogram's minimum
   Double_t per_cut = (val_cut-min)/(max-min); // normalized value of val_cut
   Double_t eps     = (max-min)*0.00001;       // epsilon

   // Definition of the two palettes below and above val_cut
   const Int_t Number     = 4;
   Double_t Red[Number]   = { 0.11, 0.19 , 0.30, 0.89};
   Double_t Green[Number] = { 0.03, 0.304, 0.60, 0.91};
   Double_t Blue[Number]  = { 0.18, 0.827, 0.50, 0.70};
   Double_t Stops[Number] = { 0., per_cut, per_cut+eps, 1. };

   Int_t nb= 256;
   h->SetContour(nb);

   TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,nb);

   // Histogram drawaing
   h->Draw("colz");
}

void perceptualcolormap() {
   TH2D *h = new TH2D("h","Perceptual Colormap",200,-4,4,200,-4,4);
   h->SetStats(0);

   Double_t a,b;
   for (Int_t i=0;i<1000000;i++) {
      gRandom->Rannor(a,b);
      h->Fill(a-1.5,b-1.5,0.1);
      h->Fill(a+2.,b-3.,0.07);
      h->Fill(a-3.,b+3.,0.05);
      gRandom->Rannor(a,b);
      h->Fill(a+1.5,b+1.5,-0.08);
   }
   Perceptual_Colormap(h, 0.);
}
