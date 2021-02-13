// @(#)root/graf2d:$Id$
// Author:  Claudi Martinez, July 19th 2010

/*************************************************************************
 * Copyright (C) 1995-2010, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/// \defgroup fitsio FITS file
/// \brief Interface to FITS file.
/// \ingroup Graphics2D
///
/// TFITS is an interface that lets you reading Flexible Image Transport System
/// (FITS) files, which are generally used in astronomy. This file format
/// was standardized 1981 and today is still widely used among professional
/// and amateur astronomers. FITS is not only an image file, but also
/// it can contain spectrums, data tables, histograms, and multidimensional
/// data. Furthermore, FITS data can be described itself by containing
/// human-readable information that let us to interpret the data within
/// the FITS file. For example, a FITS could contain a 3D data cube,
/// but an additional description would tell us that we must read it, for
/// example, as a 3-layer image.
///
/// TFITS requires CFITSIO library to be installed on your system. It
/// is currently maintained by NASA/GSFC and can be downloaded from
/// [NASA/GSFC web site](http://fits.gsfc.nasa.gov), as well as documentation.
///
/// Using this interface is easy and straightforward. There is only 1 class
/// called "TFITSHDU" which has several methods to extract data from a
/// FITS file, more specifically, from an HDU within the file. An HDU, or
/// Header Data Unit, is a chunk of data with a header containing several
/// "keyword = value" tokens. The header describes the structure of data
/// within the HDU. An HDU can be of two types: an "image HDU" or a "table
/// HDU". The former can be any kind of multidimensional array of real numbers,
/// by which the name "image" may be confusing: you can store an image, but
/// you can also store a N-dimensional data cube. On the other hand, table
/// HDUs are sets of several rows and columns (a.k.a fields) which contain
/// generic data, as strings, real or complex numbers and even arrays.
///
/// Please have a look to the tutorials ($ROOTSYS/tutorials/fitsio/) to see
/// some examples. IMPORTANT: to run tutorials it is required that
/// you change the current working directory of ROOT (CINT) shell to the
/// tutorials directory. Example:
/// ~~~ {.cpp}
/// root [1] gSystem->ChangeDirectory("tutorials/fitsio")
/// root [1] .x FITS_tutorial1.C
/// ~~~
/// LIST OF TODO
/// - Support for complex values within data tables
/// - Support for reading arrays from table cells
/// - Support for grouping
///
/// IMPLEMENTATION NOTES:
///
/// CFITSIO library uses standard C types ('int', 'long', ...). To avoid
/// confusion, the same types are used internally by the access methods.
/// However, class's fields are ROOT-defined types.

/** \class TFITSHDU
\ingroup fitsio

FITS file interface class

TFITSHDU is a class that allows extracting images and data from FITS files and contains
several methods to manage them.
*/

#include "TFITS.h"
#include "TImage.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "strlcpy.h"

#include "fitsio.h"
#include <cstdlib>

ClassImp(TFITSHDU);

////////////////////////////////////////////////////////////////////////////////
/// Clean path from possible filter and put the result in 'dst'.

void TFITSHDU::CleanFilePath(const char *filepath_with_filter, TString &dst)
{
   dst = filepath_with_filter;

   Ssiz_t ndx = dst.Index("[", 1, 0, TString::kExact);
   if (ndx != kNPOS) {
      dst.Resize(ndx);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// TFITSHDU constructor from file path with HDU selection filter.
/// Please refer to CFITSIO manual for more information about
/// HDU selection filters.
///
/// Examples:
///  - `TFITSHDU("/path/to/myfile.fits")`: just open the PRIMARY HDU
///  - `TFITSHDU("/path/to/myfile.fits[1]")`: open HDU #1
///  - `TFITSHDU("/path/to/myfile.fits[PICS]")`: open HDU called 'PICS'
///  - `TFITSHDU("/path/to/myfile.fits[ACQ][EXPOSURE > 5]")`: open the (table) HDU called 'ACQ' and
///                                                          selects the rows that have column 'EXPOSURE'
///                                                          greater than 5.

TFITSHDU::TFITSHDU(const char *filepath_with_filter)
{
   _initialize_me();

   fFilePath = filepath_with_filter;
   CleanFilePath(filepath_with_filter, fBaseFilePath);

   if (kFALSE == LoadHDU(fFilePath)) {
      _release_resources();
      throw -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TFITSHDU constructor from filepath and extension number.

TFITSHDU::TFITSHDU(const char *filepath, Int_t extension_number)
{
   _initialize_me();
   CleanFilePath(filepath, fBaseFilePath);

   //Add "by extension number" filter
   fFilePath.Form("%s[%d]", fBaseFilePath.Data(), extension_number);

   if (kFALSE == LoadHDU(fFilePath)) {
      _release_resources();
      throw -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TFITSHDU constructor from filepath and extension name.

TFITSHDU::TFITSHDU(const char *filepath, const char *extension_name)
{
   _initialize_me();
   CleanFilePath(filepath, fBaseFilePath);

   //Add "by extension number" filter
   fFilePath.Form("%s[%s]", fBaseFilePath.Data(), extension_name);


   if (kFALSE == LoadHDU(fFilePath)) {
      _release_resources();
      throw -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TFITSHDU destructor.

TFITSHDU::~TFITSHDU()
{
   _release_resources();
}

////////////////////////////////////////////////////////////////////////////////
/// Release internal resources.

void TFITSHDU::_release_resources()
{
   if (fRecords) delete [] fRecords;

   if (fType == kImageHDU) {
      if (fSizes)  delete fSizes;
      if (fPixels) delete fPixels;
   } else {
      if (fColumnsInfo) {
         if (fCells) {
            for (Int_t i = 0; i < fNColumns; i++) {
               if (fColumnsInfo[i].fType == kString) {
                  //Deallocate character arrays allocated for kString columns
                  Int_t offset = i * fNRows;
                  for (Int_t row = 0; row < fNRows; row++) {
                     delete [] fCells[offset+row].fString;
                  }
               } else if (fColumnsInfo[i].fType == kRealArray) {
                  //Deallocate character arrays allocated for kString columns
                  Int_t offset = i * fNRows;
                  for (Int_t row = 0; row < fNRows; row++) {
                     delete [] fCells[offset+row].fRealArray;
                  }
               } else if (fColumnsInfo[i].fType == kRealVector) {
                  // Deallocate character arrays allocated for variable-length array columns
                  Int_t offset = i * fNRows;
                  for (Int_t row = 0; row < fNRows; row++) {
                     delete fCells[offset + row].fRealVector;
                  }
               }
            }

            delete [] fCells;
         }

         delete [] fColumnsInfo;
      }


   }
}

////////////////////////////////////////////////////////////////////////////////
/// Do some initializations.

void TFITSHDU::_initialize_me()
{
   fRecords = 0;
   fPixels = 0;
   fSizes = 0;
   fColumnsInfo = 0;
   fNColumns = fNRows = 0;
   fCells = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Load HDU from fits file satisfying the specified filter.
/// Returns kTRUE if success. Otherwise kFALSE.
/// If filter == "" then the primary array is selected

Bool_t TFITSHDU::LoadHDU(TString& filepath_filter)
{
   fitsfile *fp=0;
   int status = 0;
   char errdescr[FLEN_STATUS+1];

   // Open file with filter
   fits_open_file(&fp, filepath_filter.Data(), READONLY, &status);
   if (status) goto ERR;

   // Read HDU number
   int hdunum;
   fits_get_hdu_num(fp, &hdunum);
   fNumber = Int_t(hdunum);

   // Read HDU type
   int hdutype;
   fits_get_hdu_type(fp, &hdutype, &status);
   if (status) goto ERR;
   fType = (hdutype == IMAGE_HDU) ? kImageHDU : kTableHDU;

   //Read HDU header records
   int nkeys, morekeys;
   char keyname[FLEN_KEYWORD+1];
   char keyvalue[FLEN_VALUE+1];
   char comment[FLEN_COMMENT+1];

   fits_get_hdrspace(fp, &nkeys, &morekeys, &status);
   if (status) goto ERR;

   fRecords = new struct HDURecord[nkeys];

   for (int i = 1; i <= nkeys; i++) {
      fits_read_keyn(fp, i, keyname, keyvalue, comment, &status);
      if (status) goto ERR;
      fRecords[i-1].fKeyword = keyname;
      fRecords[i-1].fValue = keyvalue;
      fRecords[i-1].fComment = comment;
   }

   fNRecords = Int_t(nkeys);

   //Set extension name
   fExtensionName = "PRIMARY"; //Default
   for (int i = 0; i < nkeys; i++) {
      if (fRecords[i].fKeyword == "EXTNAME") {
         fExtensionName = fRecords[i].fValue;
         break;
      }
   }

   //Read HDU's data
   if (fType == kImageHDU) {
      Info("LoadHDU", "The selected HDU contains an Image Extension");

      int param_ndims=0;
      long *param_dimsizes;

      //Read image number of dimensions
      fits_get_img_dim(fp, &param_ndims, &status);
      if (status) goto ERR;
      if (param_ndims > 0) {
         //Read image sizes in each dimension
         param_dimsizes = new long[param_ndims];
         fits_get_img_size(fp, param_ndims, param_dimsizes, &status);
         if (status) {
            delete [] param_dimsizes;
            goto ERR;
         }

         fSizes = new TArrayI(param_ndims);
         fSizes = new TArrayI(param_ndims);
         for (int i = 0; i < param_ndims; i++) { //Use for loop to copy values instead of passing array to constructor, since 'Int_t' size may differ from 'long' size
            fSizes->SetAt(param_dimsizes[i], i);
         }

         delete [] param_dimsizes;

         //Read pixels
         int anynul;
         long *firstpixel = new long[param_ndims];
         double nulval = 0;
         long npixels = 1;

         for (int i = 0; i < param_ndims; i++) {
            npixels *= (long) fSizes->GetAt(i); //Compute total number of pixels
            firstpixel[i] = 1; //Set first pixel to read from.
         }

         double *pixels = new double[npixels];

         fits_read_pix(fp, TDOUBLE, firstpixel, npixels,
                     (void *) &nulval, (void *) pixels, &anynul, &status);

         if (status) {
            delete [] firstpixel;
            delete [] pixels;
            goto ERR;
         }

         fPixels = new TArrayD(npixels, pixels);

         delete [] firstpixel;
         delete [] pixels;

      } else {
         //Null array
         fSizes = new TArrayI();
         fPixels = new TArrayD();
      }
   } else {
      Info("LoadHDU", "The selected HDU contains a Table Extension");

      // Get table's number of rows and columns
      long table_rows;
      int  table_cols;

      fits_get_num_rows(fp, &table_rows, &status);
      if (status) goto ERR;

      fNRows = Int_t(table_rows);

      fits_get_num_cols(fp, &table_cols, &status);
      if (status) goto ERR;

      fNColumns = Int_t(table_cols);

      // Allocate column info array
      fColumnsInfo = new struct Column[table_cols];

      // Read column names
      char colname[80];
      int colnum;

      fits_get_colname(fp, CASEINSEN, (char*) "*", colname, &colnum, &status);
      while (status == COL_NOT_UNIQUE)
      {
         fColumnsInfo[colnum-1].fName = colname;
         fits_get_colname(fp, CASEINSEN, (char*) "*", colname, &colnum, &status);
      }
      if (status != COL_NOT_FOUND) goto ERR;
      status = 0;

      //Allocate cells
      fCells = new union Cell [table_rows * table_cols];

      // Read columns
      int typecode;
      long repeat, width;
      Int_t cellindex;


      for (colnum = 0, cellindex = 0; colnum < fNColumns; colnum++) {
         fits_get_coltype(fp, colnum+1, &typecode, &repeat, &width, &status);

         if (status) goto ERR;

         if ((typecode == TDOUBLE) || (typecode == TSHORT) || (typecode == TLONG)
                                   || (typecode == TFLOAT) || (typecode == TLOGICAL) || (typecode == TBIT)
                                   || (typecode == TBYTE)  || (typecode == TSTRING)) {

            if (typecode == TSTRING) {
               // this column contains strings
               fColumnsInfo[colnum].fType = kString;

               int dispwidth=0;
               fits_get_col_display_width(fp, colnum+1, &dispwidth, &status);
               if (status) goto ERR;


               char *nulval = (char*) "";
               int anynul=0;
               char **array;

               if (dispwidth <= 0) {
                  dispwidth = 1;
               }

               array = new char* [table_rows];
               for (long row = 0; row < table_rows; row++) {
                  array[row] = new char[dispwidth+1]; //also room for end null!
               }

               if (repeat > 0) {
                  fits_read_col(fp, TSTRING, colnum+1, 1, 1, table_rows, nulval, array, &anynul, &status);
                  if (status) {
                     for (long row = 0; row < table_rows; row++) {
                        delete [] array[row];
                     }
                     delete [] array;
                     goto ERR;
                  }

               } else {
                  //No elements: set dummy
                  for (long row = 0; row < table_rows; row++) {
                     strlcpy(array[row], "-",dispwidth+1);
                  }
               }

               //Save values
               for (long row = 0; row < table_rows; row++) {
                  fCells[cellindex++].fString = array[row];
               }

               delete [] array; //Delete temporal array holding pointer to strings, but not delete strings themselves!


            } else {
               // this column contains either a number or a fixed-length array
               double nulval = 0;
               int anynul=0;

               fColumnsInfo[colnum].fDim = (Int_t) repeat;

               double *array = 0;
               char *arrayl = 0;

               if (repeat > 0) {

                  if (typecode == TLOGICAL) {
                     arrayl = new char[table_rows * repeat];
                     fits_read_col(fp, TLOGICAL, colnum + 1, 1, 1, table_rows * repeat, &nulval, arrayl, &anynul,
                                   &status);
                     if (status) {
                        delete[] arrayl;
                        goto ERR;
                     }
                  } else {
                     array = new double[table_rows * repeat]; // Hope you got a big machine! Ask China otherwise :-)
                     fits_read_col(fp, TDOUBLE, colnum + 1, 1, 1, table_rows * repeat, &nulval, array, &anynul,
                                   &status);
                     if (status) {
                        delete[] array;
                        goto ERR;
                     }
                  }

               } else {
                  // No elements: set dummy
                  array = new double[table_rows];
                  for (long row = 0; row < table_rows; row++) {
                     array[row] = 0.0;
                  }
               }

               // Save values
               if (repeat == 1) {
                  // this column contains scalars
                  fColumnsInfo[colnum].fType = kRealNumber;
                  if (typecode == TLOGICAL) {
                     for (long row = 0; row < table_rows; row++) {
                        int temp = (signed char)arrayl[row];
                        fCells[cellindex++].fRealNumber = (double)temp;
                     }
                     delete[] arrayl;
                  } else {
                     for (long row = 0; row < table_rows; row++) {
                        fCells[cellindex++].fRealNumber = array[row];
                     }
                     delete[] array;
                  }
               } else if (repeat > 1) {
                  // this column contains fixed-length arrays
                  fColumnsInfo[colnum].fType = kRealArray;
                  if (typecode == TLOGICAL) {
                     for (long row = 0; row < table_rows; row++) {
                        double *vec = new double[repeat];
                        long offset = row * repeat;
                        for (long component = 0; component < repeat; component++) {
                           int temp = (signed char)arrayl[offset++];
                           vec[component] = (double)temp;
                        }
                        fCells[cellindex++].fRealArray = vec;
                     }
                     delete[] arrayl;
                  } else {
                     for (long row = 0; row < table_rows; row++) {
                        double *vec = new double[repeat];
                        long offset = row * repeat;
                        for (long component = 0; component < repeat; component++) {
                           vec[component] = array[offset++];
                        }
                        fCells[cellindex++].fRealArray = vec;
                     }
                     delete[] array;
                  }
               }
            }
         } else if (typecode < 1) {
            // this column contains variable-length arrays
            fColumnsInfo[colnum].fType = kRealVector;

            // null variables needed by the fits_read_col
            double nulval = 0;
            int anynul=0;
            
            // loop over rows to derive the total memory requirement
            fColumnsInfo[colnum].fRowStart.assign(table_rows+1, 0);
            fColumnsInfo[colnum].fVarLengths.assign(table_rows, 0);
            fColumnsInfo[colnum].fRowStart[0] = 0;

            for (long row = 0; row < table_rows; row++) {
               long offset = 0;
               long repeat_row = 0; 
               
               fits_read_descript(fp, colnum+1, row+1, &repeat_row, &offset, &status);
               // store the starting of each row and the number of elements it contains 
               fColumnsInfo[colnum].fRowStart[row+1] = fColumnsInfo[colnum].fRowStart[row] + repeat_row;
               fColumnsInfo[colnum].fVarLengths[row] = repeat_row;
            }
            
           for (long row = 0; row < table_rows; row++) {
               // number of elements in the cell we want to read, i.e.
               // number of elements in the variable-length array
               int nelements = fColumnsInfo[colnum].fRowStart[row + 1] - fColumnsInfo[colnum].fRowStart[row];
               // size of the variable-length array
               const int size = fColumnsInfo[colnum].fVarLengths[row];
               // vector to store the results
               TArrayD *vec = new TArrayD(size);
               // variable-length array have negative DATATYPE
               int abstype = TMath::Abs(typecode);
             
               // define the array to load the data with the CFITSIO functions
               // a fixed arrays is needed as argument to the fits_read_col function
               // so a vector is defined and then its `.data()` pointer is passed
               // to fits_read_col
               //
               // TODO: add all type cases, for now only short, long, float and double are considered
               //
               if (abstype == 21) {
                  std::vector<short> data; 
                  data.resize(size);
                  fits_read_col(fp, abstype, colnum + 1, row + 1, 1, nelements, &nulval, data.data(), &anynul, &status);
                  for (int i = 0; i < size; i++) 
                     vec->SetAt(data[i], i);
               } else if (abstype == 41) {
                  std::vector<int> data; 
                  data.resize(size);
                  fits_read_col(fp, abstype, colnum + 1, row + 1, 1, nelements, &nulval, data.data(), &anynul, &status);
                  for (int i = 0; i < size; i++) 
                     vec->SetAt(data[i], i);
               } else if (abstype == 42) {
                  std::vector<float> data; 
                  data.resize(size);
                  fits_read_col(fp, abstype, colnum + 1, row + 1, 1, nelements, &nulval, data.data(), &anynul, &status);
                  for (int i = 0; i < size; i++) 
                     vec->SetAt(data[i], i);
              } else if (abstype == 82) {
                  std::vector<double> data; 
                  data.resize(size);
                  fits_read_col(fp, abstype, colnum + 1, row + 1, 1, nelements, &nulval, data.data(), &anynul, &status);
                  for (int i = 0; i < size; i++) 
                     vec->SetAt(data[i], i);
               } else {
                  Error("LoadHDU", "The variable-length array type in column %d is unknown", colnum + 1);
               }
               // place the vector storing the variable-length array in the corresponding cell
               fCells[cellindex++].fRealVector = vec;
            }
         } else {
            Warning("LoadHDU", "error opening FITS file. Column type %d is currently not supported", typecode);
         }
      }

      if (hdutype == ASCII_TBL) {
         // ASCII table

      } else {
         // Binary table
      }
   }

   // Close file
   fits_close_file(fp, &status);
   return kTRUE;

ERR:
   fits_get_errstatus(status, errdescr);
   Warning("LoadHDU", "error opening FITS file. Details: %s", errdescr);
   status = 0;
   if (fp) fits_close_file(fp, &status);
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get record by keyword.

struct TFITSHDU::HDURecord* TFITSHDU::GetRecord(const char *keyword)
{
   for (int i = 0; i < fNRecords; i++) {
      if (fRecords[i].fKeyword == keyword) {
         return &fRecords[i];
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the value of a given keyword. Return "" if not found.

TString& TFITSHDU::GetKeywordValue(const char *keyword)
{
   HDURecord *rec = GetRecord(keyword);
   if (rec) {
      return rec->fValue;
   } else {
      return *(new TString(""));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print records.

void TFITSHDU::PrintHDUMetadata(const Option_t *) const
{
   for (int i = 0; i < fNRecords; i++) {
      if (fRecords[i].fComment.Length() > 0) {
         printf("%-10s = %s / %s\n", fRecords[i].fKeyword.Data(), fRecords[i].fValue.Data(), fRecords[i].fComment.Data());
      } else {
         printf("%-10s = %s\n", fRecords[i].fKeyword.Data(), fRecords[i].fValue.Data());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print HDU's parent file's metadata.

void TFITSHDU::PrintFileMetadata(const Option_t *opt) const
{
   fitsfile *fp=0;
   int status = 0;
   char errdescr[FLEN_STATUS+1];
   int hducount, extnum;
   int hdutype = IMAGE_HDU;
   const char *exttype;
   char extname[FLEN_CARD]="PRIMARY"; //room enough
   int verbose = (opt[0] ? 1 : 0);

   // Open file with no filters: current HDU will be the primary one.
   fits_open_file(&fp, fBaseFilePath.Data(), READONLY, &status);
   if (status) goto ERR;

   // Read HDU count
   fits_get_num_hdus(fp, &hducount, &status);
   if (status) goto ERR;
   printf("Total: %d HDUs\n", hducount);

   extnum = 0;
   while(hducount) {
      // Read HDU type
      fits_get_hdu_type(fp, &hdutype, &status);
      if (status) goto ERR;

      if (hdutype == IMAGE_HDU) {
         exttype="IMAGE";
      } else if (hdutype == ASCII_TBL) {
         exttype="ASCII TABLE";
      } else {
         exttype="BINARY TABLE";
      }

      //Read HDU header records
      int nkeys, morekeys;
      char keyname[FLEN_KEYWORD+1];
      char keyvalue[FLEN_VALUE+1];
      char comment[FLEN_COMMENT+1];

      fits_get_hdrspace(fp, &nkeys, &morekeys, &status);
      if (status) goto ERR;

      struct HDURecord *records = new struct HDURecord[nkeys];

      for (int i = 1; i <= nkeys; i++) {
         fits_read_keyn(fp, i, keyname, keyvalue, comment, &status);
         if (status) {
            delete [] records;
            goto ERR;
         }

         records[i-1].fKeyword = keyname;
         records[i-1].fValue = keyvalue;
         records[i-1].fComment = comment;

         if (strcmp(keyname, "EXTNAME") == 0) {
            //Extension name
            strlcpy(extname, keyvalue,FLEN_CARD);
         }
      }

      //HDU info
      printf("   [%d] %s (%s)\n", extnum, exttype, extname);

      //HDU records
      if (verbose) {
         for (int i = 0; i < nkeys; i++) {
            if (comment[0]) {
               printf("      %-10s = %s / %s\n", records[i].fKeyword.Data(), records[i].fValue.Data(), records[i].fComment.Data());
            } else {
               printf("      %-10s = %s\n", records[i].fKeyword.Data(), records[i].fValue.Data());
            }
         }
      }
      printf("\n");

      delete [] records;

      hducount--;
      extnum++;
      if (hducount){
         fits_movrel_hdu(fp, 1, &hdutype, &status);
         if (status) goto ERR;
      }
   }

   // Close file
   fits_close_file(fp, &status);
   return;

ERR:
   fits_get_errstatus(status, errdescr);
   Warning("PrintFileMetadata", "error opening FITS file. Details: %s", errdescr);
   status = 0;
   if (fp) fits_close_file(fp, &status);
}

////////////////////////////////////////////////////////////////////////////////
/// Print column information

void TFITSHDU::PrintColumnInfo(const Option_t *) const
{
   if (fType != kTableHDU) {
      Warning("PrintColumnInfo", "this is not a table HDU.");
      return;
   }

   for (Int_t i = 0; i < fNColumns; i++) {
      switch (fColumnsInfo[i].fType){
         case kString:
            printf("%-20s : %s\n", fColumnsInfo[i].fName.Data(), "STRING");
            break;
         case kRealNumber:
            printf("%-20s : %s\n", fColumnsInfo[i].fName.Data(), "REAL NUMBER");
            break;
         case kRealArray:
            printf("%-20s : %s\n", fColumnsInfo[i].fName.Data(), "FIXED-LENGTH ARRAY");
            break;
         case kRealVector:
            printf("%-20s : %s\n", fColumnsInfo[i].fName.Data(), "VARIABLE-LENGTH ARRAY");
            break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print full table contents

void TFITSHDU::PrintFullTable(const Option_t *) const
{
   int printed_chars;

   if (fType != kTableHDU) {
      Warning("PrintColumnInfo", "this is not a table HDU.");
      return;
   }

   // check that the table does not contain fixed or variable length arrays
   // in that case is not possible to print a flattened table
   for (Int_t col = 0; col < fNColumns; col++){
      if (fColumnsInfo[col].fType == kRealArray) {
         Warning("PrintColumnInfo", "The table contains column with fixed-length arrays and cannot be flattened for printing.");
         return;
      }
      else if (fColumnsInfo[col].fType == kRealVector) {
         Warning("PrintColumnInfo", "The table contains column with variable-length arrays and cannot be flattened for printing.");
         return;
      }
   }

   // Dump header
   putchar('\n');
   printed_chars = 0;
   for (Int_t col = 0; col < fNColumns; col++) {
      printed_chars += printf("%-10s| ", fColumnsInfo[col].fName.Data());
   }
   putchar('\n');
   while(printed_chars--) {
      putchar('-');
   }
   putchar('\n');

   // Dump rows
   for (Int_t row = 0; row < fNRows; row++) {
      for (Int_t col = 0; col < fNColumns; col++) {
         if (fColumnsInfo[col].fType == kString) {
            printf("%-10s", fCells[col * fNRows + row].fString);
         } else if (fColumnsInfo[col].fType == kRealNumber) {
            printed_chars = printf("%.2lg", fCells[col * fNRows + row].fRealNumber);
            printed_chars -= 10;
            while (printed_chars < 0) {
               putchar(' ');
               printed_chars++;
            }
         }

         if (col <= fNColumns - 1) printf("| ");
      }
      printf("\n");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print metadata.
/// Currently supported options:
///
///  - ""  :  print HDU record data
///  - "F" :  print FITS file's extension names, numbers and types
///  - "F+":  print FITS file's extension names and types and their record data
///  - "T" :  print column information when HDU is a table
///  - "T+" : print full table (columns header and rows)

void TFITSHDU::Print(const Option_t *opt) const
{
   if ((opt[0] == 'F') || (opt[0] == 'f')) {
      PrintFileMetadata((opt[1] == '+') ? "+" : "");
   } else if ((opt[0] == 'T') || (opt[0] == 't')) {
      if (opt[1] == '+') {
         PrintFullTable("");
      } else {
         PrintColumnInfo("");
      }

   } else {
      PrintHDUMetadata("");
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Read image HDU as a displayable image. Return 0 if conversion cannot be done.
/// If the HDU seems to be a multilayer image, 'layer' parameter can be used
/// to retrieve the specified layer (starting from 0)

TImage *TFITSHDU::ReadAsImage(Int_t layer, TImagePalette *pal)
{
   if (fType != kImageHDU) {
      Warning("ReadAsImage", "this is not an image HDU.");
      return 0;
   }

   if (((fSizes->GetSize() != 2) && (fSizes->GetSize() != 3) && (fSizes->GetSize() != 4)) || ((fSizes->GetSize() == 4) && (fSizes->GetAt(3) > 1))) {
      Warning("ReadAsImage", "could not convert image HDU to image because it has %d dimensions.", fSizes->GetSize());
      return 0;
   }

   Int_t width, height;
   UInt_t pixels_per_layer;

   width  = Int_t(fSizes->GetAt(0));
   height = Int_t(fSizes->GetAt(1));

   pixels_per_layer = UInt_t(width) * UInt_t(height);

   if (  ((fSizes->GetSize() == 2) && (layer > 0))
      || (((fSizes->GetSize() == 3) || (fSizes->GetSize() == 4)) && (layer >= fSizes->GetAt(2)))) {

      Warning("ReadAsImage", "layer out of bounds.");
      return 0;
   }

   // Get the maximum and minimum pixel values in the layer to auto-stretch pixels
   Double_t maxval = 0, minval = 0;
   UInt_t i;
   Double_t pixvalue;
   Int_t offset = layer * pixels_per_layer;

   for (i = 0; i < pixels_per_layer; i++) {
      pixvalue = fPixels->GetAt(offset + i);

      if (pixvalue > maxval) {
         maxval = pixvalue;
      }

      if ((i == 0) || (pixvalue < minval)) {
         minval = pixvalue;
      }
   }

   //Build the image stretching pixels into a range from 0.0 to 255.0
   //TImage *im = new TImage(width, height);
   TImage *im = TImage::Create();
   if (!im) return 0;
   TArrayD *layer_pixels = new TArrayD(pixels_per_layer);


   if (maxval == minval) {
      //plain image
      for (i = 0; i < pixels_per_layer; i++) {
         layer_pixels->SetAt(255.0, i);
      }
   } else {
      Double_t factor = 255.0 / (maxval-minval);
      for (i = 0; i < pixels_per_layer; i++) {
         pixvalue = fPixels->GetAt(offset + i);
         layer_pixels->SetAt(factor * (pixvalue-minval), i) ;
      }
   }

   if (pal == 0) {
      // Default palette: grayscale palette
      pal = new TImagePalette(256);
      for (i = 0; i < 256; i++) {
         pal->fPoints[i] = ((Double_t)i)/255.0;
         pal->fColorAlpha[i] = 0xffff;
         pal->fColorBlue[i] = pal->fColorGreen[i] = pal->fColorRed[i] = i << 8;
      }
      pal->fPoints[0] = 0;
      pal->fPoints[255] = 1.0;

      im->SetImage(*layer_pixels, UInt_t(width), pal);

      delete pal;

   } else {
      im->SetImage(*layer_pixels, UInt_t(width), pal);
   }

   delete layer_pixels;

   return im;
}

////////////////////////////////////////////////////////////////////////////////
/// If the HDU is an image, draw the first layer of the primary array
/// To set a title to the canvas, pass it in "opt"

void TFITSHDU::Draw(Option_t *)
{
   if (fType != kImageHDU) {
      Warning("Draw", "cannot draw. This is not an image HDU.");
      return;
   }

   TImage *im = ReadAsImage(0, 0);
   if (im) {
      Int_t width = Int_t(fSizes->GetAt(0));
      Int_t height = Int_t(fSizes->GetAt(1));
      TString cname, ctitle;
      cname.Form("%sHDU", this->GetName());
      ctitle.Form("%d x %d", width, height);
      new TCanvas(cname, ctitle, width, height);
      im->Draw();
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Read image HDU as a matrix. Return 0 if conversion cannot be done
/// If the HDU seems to be a multilayer image, 'layer' parameter can be used
/// to retrieve the specified layer (starting from 0) in matrix form.
/// Options (value of 'opt'):
/// "S": stretch pixel values to a range from 0.0 to 1.0

TMatrixD* TFITSHDU::ReadAsMatrix(Int_t layer, Option_t *opt)
{
   if (fType != kImageHDU) {
      Warning("ReadAsMatrix", "this is not an image HDU.");
      return 0;
   }

   if (((fSizes->GetSize() != 2) && (fSizes->GetSize() != 3) && (fSizes->GetSize() != 4)) || ((fSizes->GetSize() == 4) && (fSizes->GetAt(3) > 1))) {
      Warning("ReadAsMatrix", "could not convert image HDU to image because it has %d dimensions.", fSizes->GetSize());
      return 0;
   }


   if (   ((fSizes->GetSize() == 2) && (layer > 0))
       || (((fSizes->GetSize() == 3) || (fSizes->GetSize() == 4)) && (layer >= fSizes->GetAt(2)))) {
      Warning("ReadAsMatrix", "layer out of bounds.");
      return 0;
   }

   Int_t width, height;
   UInt_t pixels_per_layer;
   Int_t offset;
   UInt_t i;
   TMatrixD *mat=0;

   width  = Int_t(fSizes->GetAt(0));
   height = Int_t(fSizes->GetAt(1));

   pixels_per_layer = UInt_t(width) * UInt_t(height);
   offset = layer * pixels_per_layer;

   double *layer_pixels = new double[pixels_per_layer];

   if ((opt[0] == 'S') || (opt[0] == 's')) {
      //Stretch
      // Get the maximum and minimum pixel values in the layer to auto-stretch pixels
      Double_t factor, maxval=0, minval=0;
      Double_t pixvalue;
      for (i = 0; i < pixels_per_layer; i++) {
         pixvalue = fPixels->GetAt(offset + i);

         if (pixvalue > maxval) {
            maxval = pixvalue;
         }

         if ((i == 0) || (pixvalue < minval)) {
            minval = pixvalue;
         }
      }

      if (maxval == minval) {
         //plain image
         for (i = 0; i < pixels_per_layer; i++) {
            layer_pixels[i] = 1.0;
         }
      } else {
         factor = 1.0 / (maxval-minval);
         mat = new TMatrixD(height, width);

         for (i = 0; i < pixels_per_layer; i++) {
            layer_pixels[i] = factor * (fPixels->GetAt(offset + i) - minval);
         }
      }

   } else {
      //No stretching
      mat = new TMatrixD(height, width);

      for (i = 0; i < pixels_per_layer; i++) {
         layer_pixels[i] = fPixels->GetAt(offset + i);
      }
   }

   if (mat) {
      // mat->Use(height, width, layer_pixels);
      memcpy(mat->GetMatrixArray(), layer_pixels, pixels_per_layer*sizeof(double));
   }

   delete [] layer_pixels;
   return mat;
}


////////////////////////////////////////////////////////////////////////////////
/// Read image HDU as a histogram. Return 0 if conversion cannot be done.
/// The returned object can be TH1D, TH2D or TH3D depending on data dimensionality.
/// Please, check condition (returnedValue->IsA() == TH*D::Class()) to
/// determine the object class.
///
/// NOTE: do not confuse with image histogram! This function interprets
/// the array as a histogram. It does not compute the histogram of pixel
/// values of an image! Here "pixels" are interpreted as number of entries.

TH1 *TFITSHDU::ReadAsHistogram()
{
   if (fType != kImageHDU) {
      Warning("ReadAsHistogram", "this is not an image HDU.");
      return 0;
   }

   TH1 *result=0;

   if ((fSizes->GetSize() != 1) && (fSizes->GetSize() != 2) && (fSizes->GetSize() != 3)) {
      Warning("ReadAsHistogram", "could not convert image HDU to histogram because it has %d dimensions.", fSizes->GetSize());
      return 0;
   }

   if (fSizes->GetSize() == 1) {
      //1D
      UInt_t Nx = UInt_t(fSizes->GetAt(0));
      UInt_t x;

      TH1D *h = new TH1D("", "", Int_t(Nx), 0, Nx-1);

      for (x = 0; x < Nx; x++) {
         Int_t nentries = Int_t(fPixels->GetAt(x));
         if (nentries < 0) nentries = 0; //Crop negative values
         h->Fill(x, nentries);
      }

      result = h;

   } else if (fSizes->GetSize() == 2) {
      //2D
      UInt_t Nx = UInt_t(fSizes->GetAt(0));
      UInt_t Ny = UInt_t(fSizes->GetAt(1));
      UInt_t x,y;

      TH2D *h = new TH2D("", "", Int_t(Nx), 0, Nx-1, Int_t(Ny), 0, Ny-1);

      for (y = 0; y < Ny; y++) {
         UInt_t offset = y * Nx;
         for (x = 0; x < Nx; x++) {
            Int_t nentries = Int_t(fPixels->GetAt(offset + x));
            if (nentries < 0) nentries = 0; //Crop negative values
            h->Fill(x,y, nentries);
         }
      }

      result = h;

   } else if (fSizes->GetSize() == 3) {
      //3D
      UInt_t Nx = UInt_t(fSizes->GetAt(0));
      UInt_t Ny = UInt_t(fSizes->GetAt(1));
      UInt_t Nz = UInt_t(fSizes->GetAt(2));
      UInt_t x,y,z;

      TH3D *h = new TH3D("", "", Int_t(Nx), 0, Nx-1, Int_t(Ny), 0, Ny-1, Int_t(Nz), 0, Nz-1);


      for (z = 0; z < Nz; z++) {
         UInt_t offset1 = z * Nx * Ny;
         for (y = 0; y < Ny; y++) {
            UInt_t offset2 = y * Nx;
            for (x = 0; x < Nx; x++) {
               Int_t nentries = Int_t(fPixels->GetAt(offset1 + offset2 + x));
               if (nentries < 0) nentries = 0; //Crop negative values
               h->Fill(x, y, z, nentries);
            }
         }
      }

      result = h;
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a row from the image HDU when it's a 2D array.

TVectorD* TFITSHDU::GetArrayRow(UInt_t row)
{
   if (fType != kImageHDU) {
      Warning("GetArrayRow", "this is not an image HDU.");
      return 0;
   }

   if (fSizes->GetSize() != 2) {
      Warning("GetArrayRow", "could not get row from HDU because it has %d dimensions.", fSizes->GetSize());
      return 0;
   }

   UInt_t i, offset, W,H;

   W =  UInt_t(fSizes->GetAt(0));
   H =  UInt_t(fSizes->GetAt(1));


   if (row >= H) {
      Warning("GetArrayRow", "index out of bounds.");
      return 0;
   }

   offset = W * row;
   double *v = new double[W];

   for (i = 0; i < W; i++) {
      v[i] = fPixels->GetAt(offset+i);
   }

   TVectorD *vec = new TVectorD(W, v);

   delete [] v;

   return vec;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a column from the image HDU when it's a 2D array.

TVectorD* TFITSHDU::GetArrayColumn(UInt_t col)
{
   if (fType != kImageHDU) {
      Warning("GetArrayColumn", "this is not an image HDU.");
      return 0;
   }

   if (fSizes->GetSize() != 2) {
      Warning("GetArrayColumn", "could not get row from HDU because it has %d dimensions.", fSizes->GetSize());
      return 0;
   }

   UInt_t i, W,H;

   W =  UInt_t(fSizes->GetAt(0));
   H =  UInt_t(fSizes->GetAt(1));


   if (col >= W) {
      Warning("GetArrayColumn", "index out of bounds.");
      return 0;
   }

   double *v = new double[H];

   for (i = 0; i < H; i++) {
      v[i] = fPixels->GetAt(W*i+col);
   }

   TVectorD *vec = new TVectorD(H, v);

   delete [] v;

   return vec;
}


////////////////////////////////////////////////////////////////////////////////
///Get column number given its name

Int_t TFITSHDU::GetColumnNumber(const char *colname)
{
   Int_t colnum;
   for (colnum = 0; colnum < fNColumns; colnum++) {
      if (fColumnsInfo[colnum].fName == colname) {
         return colnum;
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a string-typed column from a table HDU given its column index (>=0).

TObjArray* TFITSHDU::GetTabStringColumn(Int_t colnum)
{
   if (fType != kTableHDU) {
      Warning("GetTabStringColumn", "this is not a table HDU.");
      return 0;
   }

   if ((colnum < 0) || (colnum >= fNColumns)) {
      Warning("GetTabStringColumn", "column index out of bounds.");
      return 0;
   }

   if (fColumnsInfo[colnum].fType != kString) {
      Warning("GetTabStringColumn", "attempting to read a column that is not of type 'kString'.");
      return 0;
   }

   Int_t offset = colnum * fNRows;

   TObjArray *res = new TObjArray();
   for (Int_t row = 0; row < fNRows; row++) {
      res->Add(new TObjString(fCells[offset + row].fString));
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a string-typed column from a table HDU given its name

TObjArray* TFITSHDU::GetTabStringColumn(const char *colname)
{
   if (fType != kTableHDU) {
      Warning("GetTabStringColumn", "this is not a table HDU.");
      return 0;
   }


   Int_t colnum = GetColumnNumber(colname);

   if (colnum == -1) {
      Warning("GetTabStringColumn", "column not found.");
      return 0;
   }

   if (fColumnsInfo[colnum].fType != kString) {
      Warning("GetTabStringColumn", "attempting to read a column that is not of type 'kString'.");
      return 0;
   }

   Int_t offset = colnum * fNRows;

   TObjArray *res = new TObjArray();
   for (Int_t row = 0; row < fNRows; row++) {
      res->Add(new TObjString(fCells[offset + row].fString));
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a real number-typed column from a table HDU given its column index (>=0).

TVectorD* TFITSHDU::GetTabRealVectorColumn(Int_t colnum)
{
   if (fType != kTableHDU) {
      Warning("GetTabRealVectorColumn", "this is not a table HDU.");
      return 0;
   }

   if ((colnum < 0) || (colnum >= fNColumns)) {
      Warning("GetTabRealVectorColumn", "column index out of bounds.");
      return 0;
   }

   if (fColumnsInfo[colnum].fType == kRealArray) {
      Warning("GetTabRealVectorColumn", "attempting to read a column whose cells have embedded fixed-length arrays");
      Info("GetTabRealVectorColumn", "Use GetTabRealVectorCells() or GetTabRealVectorCell() instead.");
      return 0;
   } else if (fColumnsInfo[colnum].fType == kRealVector) {
      Warning("GetTabRealVectorColumn", "attempting to read a column whose cells have embedded variable-length arrays");
      Info("GetTabRealVectorColumn", "Use GetTabVarLengthCell() instead.");
      return 0;
   }

   Int_t offset = colnum * fNRows;

   Double_t *arr = new Double_t[fNRows];

   for (Int_t row = 0; row < fNRows; row++) {
      arr[row] = fCells[offset + row].fRealNumber;
   }

   TVectorD *res = new TVectorD();
   res->Use(fNRows, arr);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a real number-typed column from a table HDU given its name

TVectorD* TFITSHDU::GetTabRealVectorColumn(const char *colname)
{
   if (fType != kTableHDU) {
      Warning("GetTabRealVectorColumn", "this is not a table HDU.");
      return 0;
   }

   Int_t colnum = GetColumnNumber(colname);

   if (colnum == -1) {
      Warning("GetTabRealVectorColumn", "column not found.");
      return 0;
   }

   if (fColumnsInfo[colnum].fType == kRealArray) {
      Warning("GetTabRealVectorColumn", "attempting to read a column whose cells have embedded fixed-length arrays");
      Info("GetTabRealVectorColumn", "Use GetTabRealVectorCells() or GetTabRealVectorCell() instead.");
      return 0;
   } else if (fColumnsInfo[colnum].fType == kRealVector) {
      Warning("GetTabRealVectorColumn", "attempting to read a column whose cells have embedded variable-length arrays");
      Info("GetTabRealVectorColumn", "Use GetTabVarLengthCell() instead.");
      return 0;
   }

   Int_t offset = colnum * fNRows;

   Double_t *arr = new Double_t[fNRows];

   for (Int_t row = 0; row < fNRows; row++) {
      arr[row] = fCells[offset + row].fRealNumber;
   }

   TVectorD *res = new TVectorD();
   res->Use(fNRows, arr);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Change to another HDU given by "filter".
/// The parameter "filter" will be appended to the
/// FITS file's base path. For example:
/// hduObject.Change("[EVENTS][TIME > 5]");
/// Please, see documentation of TFITSHDU(const char *filepath_with_filter) constructor
/// for further information.

Bool_t TFITSHDU::Change(const char *filter)
{
   TString tmppath;
   tmppath.Form("%s%s", fBaseFilePath.Data(), filter);

   _release_resources();
   _initialize_me();

   if (kFALSE == LoadHDU(tmppath)) {
      //Failed! Restore previous hdu
      Warning("Change", "error changing HDU. Restoring the previous one...");

      _release_resources();
      _initialize_me();

      if (kFALSE == LoadHDU(fFilePath)) {
         Warning("Change", "could not restore previous HDU. This object is no longer reliable!!");
      }
      return kFALSE;
   }

   //Set new full path
   fFilePath = tmppath;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Change to another HDU given by extension_number

Bool_t TFITSHDU::Change(Int_t extension_number)
{
   TString tmppath;
   tmppath.Form("[%d]", extension_number);

   return Change(tmppath.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Get a collection of real vectors embedded in cells along a given column from a table HDU. colnum >= 0.

TObjArray *TFITSHDU::GetTabRealVectorCells(Int_t colnum)
{
   if (fType != kTableHDU) {
      Warning("GetTabRealVectorCells", "this is not a table HDU.");
      return 0;
   }

   if ((colnum < 0) || (colnum >= fNColumns)) {
      Warning("GetTabRealVectorCells", "column index out of bounds.");
      return 0;
   }

   if (fColumnsInfo[colnum].fType == kRealVector) {
      Warning("GetTabRealVectorCells", "attempting to read a column whose cells have embedded variable-length arrays");
      Info("GetTabRealVectorCells", "Use GetTabVarLengthCell() instead.");
      return 0;
   }

   Int_t offset = colnum * fNRows;

   TObjArray *res = new TObjArray();
   Int_t dim = fColumnsInfo[colnum].fDim;

   for (Int_t row = 0; row < fNRows; row++) {
      TVectorD *v = new TVectorD();
      v->Use(dim, fCells[offset + row].fRealArray);
      res->Add(v);
   }

   //Make the collection to own the allocated TVectorD objects, so when
   //destroying the collection, the vectors will be destroyed too.
   res->SetOwner(kTRUE);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a collection of real vectors embedded in cells along a given column from a table HDU by name

TObjArray *TFITSHDU::GetTabRealVectorCells(const char *colname)
{
   if (fType != kTableHDU) {
      Warning("GetTabRealVectorCells", "this is not a table HDU.");
      return 0;
   }

   Int_t colnum = GetColumnNumber(colname);

   if (colnum == -1) {
      Warning("GetTabRealVectorCells", "column not found.");
      return 0;
   }

   return GetTabRealVectorCells(colnum);
}

////////////////////////////////////////////////////////////////////////////////
/// Get a real array (with fixed or variable-length) embedded in a cell given by (row>=0, column>=0)

TVectorD *TFITSHDU::GetTabRealVectorCell(Int_t rownum, Int_t colnum)
{
   if (fType != kTableHDU) {
      Warning("GetTabRealVectorCell", "this is not a table HDU.");
      return 0;
   }

   if ((colnum < 0) || (colnum >= fNColumns)) {
      Warning("GetTabRealVectorCell", "column index out of bounds.");
      return 0;
   }

   if ((rownum < 0) || (rownum >= fNRows)) {
      Warning("GetTabRealVectorCell", "row index out of bounds.");
      return 0;
   }
   
   if (fColumnsInfo[colnum].fType == kRealVector) {
      Warning("GetTabRealVectorCells", "attempting to read a column whose cells have embedded variable-length arrays");
      Info("GetTabRealVectorCells", "Use GetTabVarLengthCell() instead.");
      return 0;
   }

   TVectorD *v = new TVectorD();
   v->Use(fColumnsInfo[colnum].fDim, fCells[(colnum * fNRows) + rownum].fRealArray);
   return v;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a real vector embedded in a cell given by (row>=0, column name)

TVectorD *TFITSHDU::GetTabRealVectorCell(Int_t rownum, const char *colname)
{
   if (fType != kTableHDU) {
      Warning("GetTabRealVectorCell", "this is not a table HDU.");
      return 0;
   }

   Int_t colnum = GetColumnNumber(colname);

   if (colnum == -1) {
      Warning("GetTabRealVectorCell", "column not found.");
      return 0;
   }

   return GetTabRealVectorCell(rownum, colnum);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the name of a column given its index (column>=0).
/// In case of error the column name is "".

const TString& TFITSHDU::GetColumnName(Int_t colnum)
{
   static TString noName;
   if (fType != kTableHDU) {
      Error("GetColumnName", "this is not a table HDU.");
      return noName;
   }

   if ((colnum < 0) || (colnum >= fNColumns)) {
      Error("GetColumnName", "column index out of bounds.");
      return noName;
   }
   return fColumnsInfo[colnum].fName;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the variable-length array contained in a cell given by (row>=0, column name)

TArrayD *TFITSHDU::GetTabVarLengthVectorCell(Int_t rownum, Int_t colnum) {

   if (fType != kTableHDU) {
      Warning("GetTabVarLengthVectorCell", "this is not a table HDU.");
      return 0;
   }

   if ((colnum < 0) || (colnum >= fNColumns)) {
      Warning("GetTabVarLengthVectorCell", "column index out of bounds.");
      return 0;
   }

   if ((rownum < 0) || (rownum >= fNRows)) {
      Warning("GetTabVarLengthVectorCell", "row index out of bounds.");
      return 0;
   } 
   
   return fCells[(colnum * fNRows) + rownum].fRealVector;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the variable-length array contained in a cell given by (row>=0, column name)

TArrayD *TFITSHDU::GetTabVarLengthVectorCell(Int_t rownum, const char *colname)
{
   if (fType != kTableHDU) {
      Warning("GetTabVarLengthVectorCell", "this is not a table HDU.");
      return 0;
   }

   Int_t colnum = GetColumnNumber(colname);

   if (colnum == -1) {
      Warning("GetTabVarLengthVectorCell", "column not found.");
      return 0;
   }

   return GetTabVarLengthVectorCell(rownum, colnum);
}
