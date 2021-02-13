// @(#)root/tree:$Id$
// Author: Fons Rademakers   30/11/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TTreeRow
\ingroup tree

Class defining interface to a row of a TTree query result.
Objects of this class are created by TTreeResult methods.

Related classes are TTreeResult.
*/

#include "TTreeRow.h"
#include "TBuffer.h"

ClassImp(TTreeRow);

////////////////////////////////////////////////////////////////////////////////
/// Single row of a query result.

TTreeRow::TTreeRow()
{
   fColumnCount = 0;
   fFields      = 0;
   fOriginal    = 0;
   fRow         = 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Single row of a query result.

TTreeRow::TTreeRow(Int_t nfields)
{
   fColumnCount = nfields;
   fFields      = 0;
   fOriginal    = 0;
   fRow         = 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Single row of a query result.

TTreeRow::TTreeRow(Int_t nfields, const Int_t *fields, const char *row)
{
   fColumnCount = nfields;
   fFields      = 0;
   fOriginal    = 0;
   fRow         = 0;
   SetRow(fields,row);
}

////////////////////////////////////////////////////////////////////////////////
/// This is a shallow copy of a real row, i.e. it only contains
/// a pointer to the original.

TTreeRow::TTreeRow(TSQLRow *original)
{
   fFields      = 0;
   fOriginal    = 0;
   fColumnCount = 0;
   fRow         = 0;

   if (!original) {
      Error("TTreeRow", "original may not be 0");
      return;
   }
   if (original->IsA() != TTreeRow::Class()) {
      Error("TTreeRow", "original must be a TTreeRow");
      return;
   }

   fOriginal = (TTreeRow*) original;
   fColumnCount = fOriginal->fColumnCount;
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy row object.

TTreeRow::~TTreeRow()
{
   if (fFields)
      Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close row.

void TTreeRow::Close(Option_t *)
{
   if (fRow)    delete [] fRow;
   if (fFields) delete [] fFields;
   fColumnCount = 0;
   fOriginal = 0;
   fRow = nullptr;
   fFields = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if row is open and field index within range.

Bool_t TTreeRow::IsValid(Int_t field)
{
   if (!fFields && !fOriginal) {
      Error("IsValid", "row closed");
      return kFALSE;
   }
   if (field < 0 || field >= fColumnCount) {
      Error("IsValid", "field index out of bounds");
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get length in bytes of specified field.

ULong_t TTreeRow::GetFieldLength(Int_t field)
{
   if (!IsValid(field))
      return 0;

   if (fOriginal)
      return fOriginal->GetFieldLength(field);

   if (field > 0) return fFields[field] - fFields[field-1] -1;
   else           return fFields[0] -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get specified field from row (0 <= field < GetFieldCount()).

const char *TTreeRow::GetField(Int_t field)
{
   if (!IsValid(field))
      return 0;

   if (fOriginal)
      return fOriginal->GetField(field);

   if (field > 0) return fRow +fFields[field-1];
   else           return fRow;
}

////////////////////////////////////////////////////////////////////////////////
/// The field and row information.

void TTreeRow::SetRow(const Int_t *fields, const char *row)
{
   if (!fColumnCount) return;
   if (fFields) delete [] fFields;
   Int_t nch    = fields[fColumnCount-1];
   fFields      = new Int_t[fColumnCount];
   fOriginal    = 0;
   if (fRow) delete [] fRow;
   fRow         = new char[nch];
   for (Int_t i=0;i<fColumnCount;i++) fFields[i] = fields[i];
   memcpy(fRow,row,nch);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TTreeRow.

void TTreeRow::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      R__b.ReadVersion(&R__s, &R__c);
      TSQLRow::Streamer(R__b);
      R__b >> fColumnCount;
      fFields = new Int_t[fColumnCount];
      R__b.ReadFastArray(fFields,fColumnCount);
      Int_t nch;
      R__b >> nch;
      fRow = new char[nch];
      R__b.ReadFastArray(fRow,nch);
      R__b.CheckByteCount(R__s, R__c, TTreeRow::IsA());
   } else {
      R__c = R__b.WriteVersion(TTreeRow::Class(),kTRUE);
      TSQLRow::Streamer(R__b);
      R__b << fColumnCount;
      R__b.WriteFastArray(fFields,fColumnCount);
      Int_t nch = fFields ? fFields[fColumnCount-1] : 0;
      R__b << nch;
      R__b.WriteFastArray(fRow,nch);
      R__b.SetByteCount(R__c,kTRUE);
   }
}
