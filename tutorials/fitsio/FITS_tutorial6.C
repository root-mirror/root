/// \file
/// \ingroup tutorial_FITS
/// \notebook
/// Open a FITS file whose primary array represents
/// a spectrum table (flux vs wavelength) and dump its columns
///
/// \macro_code
/// \macro_output
///
/// \author Claudi Martinez

void FITS_tutorial6()
{
   // We open a table from a FITS file
   // and dump its columns.

   TString dir = gROOT->GetTutorialDir();

   //Open the table
   TFITSHDU hdu(dir + "/fitsio/sample4.fits[1]");

   // Show columns
   const auto nColumns = hdu.GetTabNColumns();
   printf("The table has %d columns:\n", nColumns);
   for (auto i : ROOT::TSeqI(nColumns)) {
      printf(" - Column %d: %s\n", i, hdu.GetColumnName(i).Data());
   }
}
