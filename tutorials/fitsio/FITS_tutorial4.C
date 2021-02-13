/// \file
/// \ingroup tutorial_FITS
/// \notebook
/// Open a FITS file whose primary array represents
/// a spectrum (flux vs wavelength).
///
/// \macro_code
/// \macro_output
///
/// \author Claudi Martinez

void FITS_tutorial4()
{
   // We open a FITS file that contains the primary HDU and a little data table.
   // The data table is extension #1 and it has 2 rows.
   // We want to read only the rows that have the column named DATAMAX greater than 2e-15 (there's only 1
   // matching row Data copyright: NASA

   TString dir = gROOT->GetTutorialDir();

   // Open the table extension number 1)
   TFITSHDU hdu(dir + "/fitsio/sample2.fits[1][DATAMAX > 2e-15]");

   hdu.Print("T");

   hdu.Print("T+");

   std::unique_ptr<TVectorD> vp(hdu.GetTabRealVectorColumn("DATAMAX"));
   const auto &v = *vp;
   std::cout << "v[0] = " << v[0] << std::endl;
   std::cout << "Does the matched row have DATAMAX > 2e-15? :-)" << std::endl;
}
