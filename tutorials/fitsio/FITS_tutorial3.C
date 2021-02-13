/// \file
/// \ingroup tutorial_FITS
/// \notebook -draw
/// Open a FITS file and retrieve the first plane of the image array
/// as a TImage object.
///
/// \macro_image
/// \macro_code
/// \macro_output
///
/// \author Claudi Martinez

void FITS_tutorial3()
{
   // We open a FITS file that contains several image
   // extensions. The primary HDU contains no data.
   // Data copyright: NASA

   // Open extensions 1 to 5 from file
   TString dir = gROOT->GetTutorialDir();

   auto c = new TCanvas("c1", "FITS tutorial #1", 800, 700);
   c->Divide(2, 3);
   for (auto i : ROOT::TSeqI(1, 6)) {
      TFITSHDU hdu(dir + "/fitsio/sample3.fits", i);

      TImage* im = (TImage *)hdu.ReadAsImage(0);
      c->cd(i);
      im->Draw();
   }
}
