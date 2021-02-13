/// \file
/// \ingroup tutorial_FITS
/// \notebook
/// Open a FITS file that contains a catalog of astronomical objects
/// and dump some of its columns
///
/// \macro_code
/// \macro_output
///
/// \author Elizabeth Buckley-Geer

void FITS_tutorial7()
{
   // We are going to open a table from a FITS file
   // and print out three columns for some of the objects.
   // This table contains a logical data type so this tutorial tests
   // that we can read it correctly

   TString dir = gROOT->GetTutorialDir();

   // Open the table
   TFITSHDU hdu(dir + "/fitsio/sample5.fits[1]");

   // Read the ra, dec, flux_g and brick_primary columns

   std::unique_ptr<TVectorD> vec1(hdu.GetTabRealVectorColumn("ra"));
   std::unique_ptr<TVectorD> vec2(hdu.GetTabRealVectorColumn("dec"));
   std::unique_ptr<TVectorD> vec3(hdu.GetTabRealVectorColumn("flux_g"));
   std::unique_ptr<TVectorD> vec4(hdu.GetTabRealVectorColumn("brick_primary"));

   for (auto i : ROOT::TSeqI(vec1->GetLwb(), vec1->GetUpb())) {
      const auto bp = (*vec4)[i];
      if (bp) {
         const auto gflux = (*vec3)[i];
         const auto ra = (*vec1)[i];
         const auto dec = (*vec2)[i];
         printf("RA %f DEC %f G-FLUX %f\n", ra, dec, gflux);
      }
   }
}