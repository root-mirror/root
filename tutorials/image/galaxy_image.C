/// \file
/// \ingroup tutorial_image
/// \notebook
/// Importing an image and manipulating it.
///
/// \image html image_galaxy_image.png
/// \macro_code
///
/// \author Valeriy Onuchin

void galaxy_image()
{
   TCanvas *canv = new TCanvas("image", "n4254", 40, 40, 812, 700);
   canv->ToggleEventStatus();
   canv->SetRightMargin(0.2);
   canv->SetLeftMargin(0.01);
   canv->SetTopMargin(0.01);
   canv->SetBottomMargin(0.01);

   // read the pixel data from file "galaxy.root"
   // the size of the image is 401 X 401 pixels
   const char *fname = "galaxy.root";
   TFile *gal = 0;
   if (!gSystem->AccessPathName(fname)) {
      gal = TFile::Open(fname);
   } else {
      printf("accessing %s file from http://root.cern/files\n",fname);
      gal = TFile::Open(Form("http://root.cern/files/%s",fname));
   }
   if (!gal) return;
   TImage *img = (TImage*)gal->Get("n4254");
   img->Draw();

   // open the color editor
   img->StartPaletteEditor();

   // zoom the image
   img->Zoom(80, 80, 250, 250);
}
