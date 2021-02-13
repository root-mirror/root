/// \file
/// \ingroup tutorial_FITS
/// \notebook -draw
/// Open a FITS file whose primary array represents
/// a spectrum (flux vs wavelength).
///
/// \macro_image
/// \macro_code
/// \macro_output
///
/// \author Claudi Martinez

void FITS_tutorial2()
{
   // We're gonna open a FITS file that contains the primary HDU and a little data table.
   // The primary HDU is an array of 2 rows by 2040 columns, and they represent a radiation
   // spectrum. The first row contains the flux data, whereas the second row the wavelengths.
   // Data copyright: NASA

   TString dir = gROOT->GetTutorialDir();

   // Open primary HDU from file
   TFITSHDU hdu(dir + "/fitsio/sample2.fits");

   // Dump the HDUs within the FITS file
   // and also their metadata
   hdu.Print("F+");

   // We now generate a TGraph from vectors
   std::unique_ptr<TVectorD> Y(hdu.GetArrayRow(0));
   std::unique_ptr<TVectorD> X(hdu.GetArrayRow(1));
   TGraph gr(*X,*Y);

   // Show the graphic
   auto c = new TCanvas("c1", "FITS tutorial #2", 800, 800);
   gr.SetFillColor(kRed);
   gr.DrawClone("BA");
}
