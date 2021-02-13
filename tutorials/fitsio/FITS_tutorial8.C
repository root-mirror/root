/// \file
/// \ingroup tutorial_FITS
/// \notebook
/// Open a FITS file with columns containing variable-length arrays.
///
/// ABOUT THE DATA: To illustrate such a case we use a Redistribution
/// Matrix File (RMF) conform with the HEASARC OGIP specifications
/// OGIP definition
/// https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html
/// RMF defintion
/// https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#tth_sEc3.1
///
/// Variable-length arrays are used to return the non-null row elements of a
/// sparse matrix containing the probability of detecting a photon of a given
/// energy at a given detector channel.
///
/// The test data file, rmf.fits, is taken from the test data of the sherpa X-ray
/// analysis tools
/// https://cxc.harvard.edu/sherpa/
/// the original file is available under the repository
/// https://github.com/sherpa/sherpa-test-data
///
/// \macro_code
/// \macro_output
///
/// \author Cosimo Nigro

void FITS_tutorial8()
{
   // let us read the first header data unit, opening the file with a FITS viewer
   // we see that the following cells contain variable-length arrays with values
   // row 215 column "F_CHAN": (1, 68)
   // row 215 column "N_CHAN": (66, 130)
   // row 215 column "MATRIX" : (5.8425176E-6, 7.290097E-6, 8.188037E-6, 9.157882E-6, 1.018355E-5, ..)
   TString dir = gROOT->GetTutorialDir();
   TFITSHDU* hdu = new TFITSHDU(dir + "/fitsio/rmf.fits", 1);
   int rownum = 214; // FITS tables are indexed starting from 1
   TString colname1 = "F_CHAN";
   TString colname2 = "N_CHAN";
   TString colname3 = "MATRIX";

   printf("reading in row %d, column %s \n", rownum+1, colname1.Data());
   TArrayD *arr1 = hdu->GetTabVarLengthVectorCell(rownum, colname1);
   printf("(%f, %f) \n", arr1->At(0), arr1->At(1));

   printf("reading in row %d, column %s \n", rownum+1, colname2.Data());
   TArrayD *arr2 = hdu->GetTabVarLengthVectorCell(rownum, colname2);
   printf("(%f, %f) \n", arr2->At(0), arr2->At(1));

   // printing only the first 5 values in the array
   printf("reading in row %d, column %s \n", rownum+1, colname3.Data());
   TArrayD *arr3 = hdu->GetTabVarLengthVectorCell(rownum, colname3);
   printf("(%e, %e, %e, %e, %e, ...) \n", arr3->At(0), arr3->At(1), arr3->At(2), arr3->At(3), arr3->At(4));

}
