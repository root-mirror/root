/// \file
/// \ingroup tutorial_FITS
/// \notebook -draw
/// Open a FITS file and retrieve the first plane of the image array
/// as a TImage object
///
/// \macro_image
/// \macro_code
/// \macro_output
///
/// \author Claudi Martinez

void FITS_tutorial1()
{
   // Here we open a FITS file that contains only the primary HDU, consisting on an image.
   // The object you will see is a snapshot of the NGC7662 nebula,
   // which was taken by the author on November 2009 in Barcelona.

   TString dir = gROOT->GetTutorialDir();

   // Open primary HDU from file
   TFITSHDU hdu(dir + "/fitsio/sample1.fits");

   // Dump the HDUs within the FITS file
   // and also their metadata
   // printf("Press ENTER to see summary of all data stored in the file:"); getchar();

   hdu.Print("F+");

   // Here we get the exposure time.
   printf("Exposure time = %s\n", hdu.GetKeywordValue("EXPTIME").Data());

   // Read the primary array as a matrix, selecting only layer 0.
   // This function may be useful to do image processing, e.g. custom filtering

   std::unique_ptr<TMatrixD> mat(hdu.ReadAsMatrix(0));
   mat->Print();

   // Read the primary array as an image, selecting only layer 0.
   TImage * im = (TImage *)hdu.ReadAsImage(0);

   // Read the primary array as a histogram. Depending on array dimensions, the returned
   // histogram will be 1D, 2D or 3D.
   TH1 * hist = (TH1 *)hdu.ReadAsHistogram();

   auto c = new TCanvas("c1", "FITS tutorial #1", 1400, 800);
   c->Divide(2, 1);
   c->cd(1);
   im->Draw();
   c->cd(2);
   hist->Draw("COL");
}
