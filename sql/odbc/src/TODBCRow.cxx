// @(#)root/odbc:$Id$
// Author: Sergey Linev   6/02/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TODBCRow.h"

#include <sqlext.h>
#include "strlcpy.h"

ClassImp(TODBCRow);

////////////////////////////////////////////////////////////////////////////////
/// Single row of query result.

TODBCRow::TODBCRow(SQLHSTMT stmt, Int_t fieldcount)
{
   fHstmt = stmt;
   fFieldCount = fieldcount;

   fBuffer = 0;
   fLengths = 0;

   if (fFieldCount>0) {
      fBuffer = new char*[fFieldCount];
      fLengths = new ULong_t[fFieldCount];
      for (Int_t n = 0; n < fFieldCount; n++) {
         fBuffer[n] = 0;
         fLengths[n] = 0;
         CopyFieldValue(n);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy row object.

TODBCRow::~TODBCRow()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close row.

void TODBCRow::Close(Option_t *)
{
   if (fBuffer!=0) {
      for (Int_t n = 0; n < fFieldCount; n++)
         delete[] fBuffer[n];
      delete[] fBuffer;
      fBuffer = 0;
   }

   if (fLengths!=0) {
      delete[] fLengths;
      fLengths = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Extracts field value from statement.
/// First allocates 128 bytes for buffer.
/// If there is not enouth space, bigger buffer is allocated and
/// request is repeated

void TODBCRow::CopyFieldValue(Int_t field)
{
   #define buffer_len 128

   fBuffer[field] = new char[buffer_len];

   SQLLEN ressize;

   SQLRETURN retcode = SQLGetData(fHstmt, field+1, SQL_C_CHAR, fBuffer[field], buffer_len, &ressize);

   if (ressize==SQL_NULL_DATA) {
      delete[] fBuffer[field];
      fBuffer[field] = 0;
      return;
   }

   fLengths[field] = ressize;

   if (retcode==SQL_SUCCESS_WITH_INFO) {
      SQLINTEGER code;
      SQLCHAR state[ 7 ];
      SQLGetDiagRec(SQL_HANDLE_STMT, fHstmt, 1, state, &code, 0, 0, 0);

      if (strcmp((char*)state,"01004")==0) {
//         Info("CopyFieldValue","Before %d %s", ressize, fBuffer[field]);

         char* newbuf = new char[ressize+10];
         strlcpy(newbuf, fBuffer[field], buffer_len);
         delete fBuffer[field];
         fBuffer[field] = newbuf;
         newbuf+=(buffer_len-1); // first data will not be read again
         retcode = SQLGetData(fHstmt, field+1, SQL_C_CHAR, newbuf, ressize+10-buffer_len, &ressize);

//         Info("CopyFieldValue","After %d %s", ressize, fBuffer[field]);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get length in bytes of specified field.

ULong_t TODBCRow::GetFieldLength(Int_t field)
{
   if ((field<0) || (field>=fFieldCount)) return 0;

   return fLengths[field];
}

////////////////////////////////////////////////////////////////////////////////
/// Get specified field from row (0 <= field < GetFieldCount()).

const char *TODBCRow::GetField(Int_t field)
{
   if ((field<0) || (field>=fFieldCount)) return 0;

   return fBuffer[field];
}
