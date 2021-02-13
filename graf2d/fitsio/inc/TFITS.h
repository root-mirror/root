// @(#)root/graf2d:$Id$
// Author: Claudi Martinez, July 19th 2010

/*************************************************************************
 * Copyright (C) 1995-2010, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFITS
#define ROOT_TFITS

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFITS                                                                //
//                                                                      //
// Interface to FITS astronomical files.                                //
// Please, see TFITS.cxx for info about implementation                  //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TMatrixDfwd.h"
#include "TVectorDfwd.h"
#include <vector>

class TArrayI;
class TArrayD;
class TH1;
class TImage;
class TImagePalette;
class TObjArray;

class TFITSHDU : public TNamed {

private:
   void _release_resources();
   void _initialize_me();

public:
   enum EHDUTypes {         // HDU types
      kImageHDU,
      kTableHDU
   };

   enum EColumnTypes { // Column data types
      kString,         // the column contains a string
      kRealNumber,     // the column contains a scalar
      kRealArray,      // the column contains a fixed-length array
      kRealVector      // the column contains a variable-length array
   };

   struct HDURecord {       // FITS HDU record
      TString fKeyword;
      TString fValue;
      TString fComment;
   };

   struct Column {               // Information of a table column
      TString            fName;       // Column's name
      enum EColumnTypes  fType;       // Column's data type
      Int_t              fDim;        // When cells contain real number arrays, this field indicates
                                      // the dimension of this vector (number of components), being 1 for scalars.
      Bool_t             fVariable;   // kTRUE if the columns has variable-length arrays
      std::vector<int>   fRowStart;   // starting index of each row
      std::vector<int>   fVarLengths; // sizes of the individual varaiable-length arrays
   };

   union Cell {                 //Table cell contents
      Char_t       *fString;
      Double_t      fRealNumber;
      Double_t     *fRealArray;
      TArrayD      *fRealVector;
   };

protected:
   TString             fFilePath;         ///< Path to HDU's file including filter
   TString             fBaseFilePath;     ///< Path to HDU's file excluding filter
   struct HDURecord   *fRecords;          ///< HDU metadata records
   Int_t               fNRecords;         ///< Number of records
   enum EHDUTypes      fType;             ///< HDU type
   TString             fExtensionName;    ///< Extension Name
   Int_t               fNumber;           ///< HDU number (1=PRIMARY)
   TArrayI            *fSizes;            ///< Image sizes in each dimension (when fType == kImageHDU)
   TArrayD            *fPixels;           ///< Image pixels (when fType == kImageHDU)
   struct Column      *fColumnsInfo;      ///< Information about columns (when fType == kTableHDU)
   Int_t               fNColumns;         ///< Number of columns (when fType == kTableHDU)
   Int_t               fNRows;            ///< Number of rows (when fType == kTableHDU)
   union  Cell        *fCells;            ///< Table cells (when fType == kTableHDU). Cells are ordered in the following way:
                                          ///< fCells[0..fNRows-1] -> cells of column 0
                                          ///< fCells[fNRows..2*fNRows-1] -> cells of column 1
                                          ///< fCells[2*fNRows..3*fNRows-1] -> cells of column 2
                                          ///< fCells[(fNColumns-1)*fNRows..fNColumns*fNRows-1] -> cells of column fNColumns-1


   Bool_t            LoadHDU(TString& filepath_filter);
   static void       CleanFilePath(const char *filepath_with_filter, TString &dst);
   void              PrintHDUMetadata(const Option_t *opt="") const;
   void              PrintFileMetadata(const Option_t *opt="") const;
   void              PrintColumnInfo(const Option_t *) const;
   void              PrintFullTable(const Option_t *) const;

public:
   TFITSHDU(const char *filepath_with_filter);
   TFITSHDU(const char *filepath, Int_t extension_number);
   TFITSHDU(const char *filepath, const char *extension_name);
   ~TFITSHDU();

   //Metadata access methods
   Int_t              GetRecordNumber() const { return fNRecords; }
   struct HDURecord  *GetRecord(const char *keyword);
   TString&           GetKeywordValue(const char *keyword);
   void               Print(const Option_t *opt="") const;

   //Image readers
   TH1               *ReadAsHistogram();
   TImage            *ReadAsImage(Int_t layer = 0, TImagePalette *pal = 0);
   TMatrixD          *ReadAsMatrix(Int_t layer = 0, Option_t *opt="");
   TVectorD          *GetArrayRow(UInt_t row);
   TVectorD          *GetArrayColumn(UInt_t col);

   //Table readers
   Int_t              GetTabNColumns() const { return fNColumns; }
   Int_t              GetTabNRows()    const { return fNRows; }
   Int_t              GetColumnNumber(const char *colname);
   const TString&     GetColumnName(Int_t colnum);
   TObjArray         *GetTabStringColumn(Int_t colnum);
   TObjArray         *GetTabStringColumn(const char *colname);
   TVectorD          *GetTabRealVectorColumn(Int_t colnum);
   TVectorD          *GetTabRealVectorColumn(const char *colname);
   TVectorD          *GetTabRealVectorCell(Int_t rownum, Int_t colnum);
   TVectorD          *GetTabRealVectorCell(Int_t rownum, const char *colname);
   TObjArray         *GetTabRealVectorCells(Int_t colnum);
   TObjArray         *GetTabRealVectorCells(const char *colname);
   TArrayD           *GetTabVarLengthVectorCell(Int_t rownum, Int_t colnum);
   TArrayD           *GetTabVarLengthVectorCell(Int_t rownum, const char *colname);

   //Misc
   void               Draw(Option_t *opt="");
   Bool_t             Change(const char *filter);
   Bool_t             Change(Int_t extension_number);


   ClassDef(TFITSHDU,0)  // Class interfacing FITS HDUs
};


#endif
