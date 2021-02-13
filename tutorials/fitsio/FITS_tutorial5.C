/// \file
/// \ingroup tutorial_FITS
/// \notebook
/// Open a FITS file whose primary array represents
/// a spectrum (flux vs wavelength)
///
/// \macro_code
/// \macro_output
///
/// \author Claudi Martinez

using Upvd_t = std::unique_ptr<TVectorD>;

void FITS_tutorial5()
{
   // We open a FITS file that contains a table with 9 rows and 8 columns. Column 4 has name
   // 'mag' and contains a vector of 6 numeric components. The values of vectors in rows 1 and 2 (column 4) are:
   // Row1: (99.0, 24.768, 23.215, 21.68, 21.076, 20.857)
   // Row2: (99.0, 21.689, 20.206, 18.86, 18.32 , 18.128 )
   // WARNING: when coding, row and column indices start from 0

   TString dir = gROOT->GetTutorialDir();

   // Open the table
   TFITSHDU hdu(dir + "/fitsio/sample4.fits[1]");

   // Read vectors at rows 1 and 2 (indices 0 and 1)
   std::array<Upvd_t, 2> vs{Upvd_t(hdu.GetTabRealVectorCell(0, "mag")), Upvd_t(hdu.GetTabRealVectorCell(1, "mag"))};
   for (auto &&v : vs) {
      for(auto i : ROOT::TSeqI(v->GetNoElements())) {
         if (i > 0)
            printf(", ");
         printf("%lg", (*v)[i]); // NOTE: the asterisk is for using the overloaded [] operator of the TVectorD object
      }
      printf(")\n");
   }

   // We now dump all rows using the function GetTabRealVectorCells()
   std::unique_ptr<TObjArray> vectorCollection(hdu.GetTabRealVectorCells("mag"));
   for (auto vObj : *vectorCollection) {
      auto &v = *static_cast<TVectorD*>(vObj);
      for (auto i : ROOT::TSeqI(v.GetNoElements())) {
         if (i > 0) printf(", ");
         printf("%lg", (v[i]));
      }
      printf(")\n");
   }
}
