// @(#)root/oracle:$Id$
// Author: Yan Liu and Shaowen Wang   23/11/04

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TOracleRow.h"
#include "TOracleServer.h"
#include "snprintf.h"
#include <cstring>
#include <ctime>

ClassImp(TOracleRow);

using namespace std;
using namespace oracle::occi;


////////////////////////////////////////////////////////////////////////////////
/// Single row of query result.

TOracleRow::TOracleRow(ResultSet *rs, vector<MetaData> *fieldMetaData)
{
   fResult      = rs;
   fFieldInfo   = fieldMetaData;
   fFieldCount  = fFieldInfo->size();

   fFieldsBuffer = nullptr;

   GetRowData();
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy row object.

TOracleRow::~TOracleRow()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close row.

void TOracleRow::Close(Option_t *)
{
   if (fFieldsBuffer) {
      for (int n=0;n<fFieldCount;n++)
        if (fFieldsBuffer[n])
           delete[] fFieldsBuffer[n];
      delete [] fFieldsBuffer;
   }

   fFieldsBuffer = nullptr;
   fFieldInfo   = nullptr;
   fFieldCount  = 0;
   fResult      = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if row is open and field index within range.

Bool_t TOracleRow::IsValid(Int_t field)
{
   if (!fResult) {
      Error("IsValid", "row closed");
      return kFALSE;
   }
   if (field < 0 || field >= (Int_t)fFieldInfo->size()) {
      Error("IsValid", "field index out of bounds");
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get length in bytes of specified field.

ULong_t TOracleRow::GetFieldLength(Int_t field)
{
   if (!IsValid(field) || fFieldInfo->size() <= 0)
      return 0;

   MetaData fieldMD = (*fFieldInfo)[field];

   return fieldMD.getInt(MetaData::ATTR_DATA_SIZE);
}

////////////////////////////////////////////////////////////////////////////////

const char* TOracleRow::GetField(Int_t field)
{
   if ((field<0) || (field>=fFieldCount)) {
      Error("TOracleRow","GetField(): out-of-range or No RowData/ResultSet/MetaData");
      return nullptr;
   }

   return fFieldsBuffer ? fFieldsBuffer[field] : nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void TOracleRow::GetRowData()
{
   if (!fResult || !fFieldInfo || (fFieldCount<=0)) return;

   fFieldsBuffer = new char* [fFieldCount];
   for (int n=0;n<fFieldCount;n++)
     fFieldsBuffer[n] = 0;

   std::string res;

   char str_number[200];

   int fPrecision, fScale, fDataType;
   double double_val;

   try {

   for (int field=0;field<fFieldCount;field++) {
      if (fResult->isNull(field+1)) continue;

      fDataType = (*fFieldInfo)[field].getInt(MetaData::ATTR_DATA_TYPE);

      switch (fDataType) {
        case SQLT_NUM: //NUMBER
           fPrecision = (*fFieldInfo)[field].getInt(MetaData::ATTR_PRECISION);
           fScale = (*fFieldInfo)[field].getInt(MetaData::ATTR_SCALE);

           if ((fScale == 0) || (fPrecision == 0)) {
              res = fResult->getString(field+1);
           } else {
              double_val = fResult->getDouble(field+1);
              snprintf(str_number, sizeof(str_number), TSQLServer::GetFloatFormat(), double_val);
              res = str_number;
           }
           break;

        case SQLT_CHR:  // character string
        case SQLT_VCS:  // variable character string
        case SQLT_AFC: // ansi fixed char
        case SQLT_AVC: // ansi var char
           res = fResult->getString(field+1);
           break;
        case SQLT_DAT:  // Oracle native DATE type
           res = (fResult->getDate(field+1)).toText(TOracleServer::GetDatimeFormat());
           break;
        case SQLT_TIMESTAMP:     // TIMESTAMP
        case SQLT_TIMESTAMP_TZ:  // TIMESTAMP WITH TIMEZONE
        case SQLT_TIMESTAMP_LTZ: // TIMESTAMP WITH LOCAL TIMEZONE
           res = (fResult->getTimestamp(field+1)).toText(TOracleServer::GetDatimeFormat(), 0);
           break;
        case SQLT_IBFLOAT:
        case SQLT_IBDOUBLE:
           res = fResult->getString(field+1);
           break;
        default:
           Error("GetRowData","Oracle type %d was not yet tested - please inform ROOT developers", fDataType);
           continue;
      }

      int len = res.length();
      if (len>0) {
         fFieldsBuffer[field] = new char[len+1];
         strcpy(fFieldsBuffer[field], res.c_str());
      }
   }

   } catch (SQLException &oraex) {
      Error("GetRowData", "%s", (oraex.getMessage()).c_str());
   }
}
