// @(#)root/odbc:$Id$
// Author: Sergey Linev   6/02/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


//________________________________________________________________________
//
//  SQL statement class for ODBC
//
//  See TSQLStatement class documentation for more details
//
//________________________________________________________________________


#include "TODBCStatement.h"
#include "TODBCServer.h"
#include "TDataType.h"
#include "strlcpy.h"
#include "snprintf.h"

#include <sqlext.h>
#include <cstdlib>
#include <iostream>

#define kSqlTime      123781
#define kSqlDate      123782
#define kSqlTimestamp 123783
#define kSqlBinary    123784


ClassImp(TODBCStatement);

////////////////////////////////////////////////////////////////////////////////
///constructor

TODBCStatement::TODBCStatement(SQLHSTMT stmt, Int_t rowarrsize, Bool_t errout) :
   TSQLStatement(errout)
{
   fHstmt = stmt;
   fBufferPreferredSize = rowarrsize;

   fBuffer = 0;
   fStatusBuffer = 0;
   fNumBuffers = 0;
   fBufferLength = 0;
   fBufferCounter = 0;

   fWorkingMode = 0;

   fNumParsProcessed = 0;
   fNumRowsFetched = 0;

   SQLSMALLINT   paramsCount = 0;
   SQLRETURN retcode = SQLNumParams(fHstmt, &paramsCount);
   if (ExtractErrors(retcode,"Constructor"))
      paramsCount = 0;

   if (paramsCount>0) {

      fWorkingMode = 1; // we are now using buffers for parameters
      fNumParsProcessed = 0;

      SQLSetStmtAttr(fHstmt, SQL_ATTR_PARAM_BIND_TYPE, SQL_PARAM_BIND_BY_COLUMN, 0);

      SQLUINTEGER setsize = fBufferPreferredSize;
      retcode = SQLSetStmtAttr(fHstmt, SQL_ATTR_PARAMSET_SIZE, (SQLPOINTER) (long) setsize, 0);
      ExtractErrors(retcode,"Constructor");

      SQLUINTEGER getsize = 0;

      retcode = SQLGetStmtAttr(fHstmt, SQL_ATTR_PARAMSET_SIZE, &getsize, 0, 0);
      ExtractErrors(retcode,"Constructor");

      Int_t bufferlen = fBufferPreferredSize;

      // MySQL is not yet support array of parameters
      if (getsize<=1) bufferlen=1; else
      if (getsize!=setsize) {
         SQLSetStmtAttr(fHstmt, SQL_ATTR_PARAMSET_SIZE, (SQLPOINTER) 1, 0);
         bufferlen = 1;
      }

      SetNumBuffers(paramsCount, bufferlen);

      SQLSetStmtAttr(fHstmt, SQL_ATTR_PARAM_STATUS_PTR, fStatusBuffer, 0);
      SQLSetStmtAttr(fHstmt, SQL_ATTR_PARAMS_PROCESSED_PTR, &fNumParsProcessed, 0);

      // indicates that we are starting
      fBufferCounter = -1;
   }

   fNumRowsFetched = 0;
   fLastResultRow = 0;
}

////////////////////////////////////////////////////////////////////////////////
///destructor

TODBCStatement::~TODBCStatement()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close statement

void TODBCStatement::Close(Option_t *)
{
   FreeBuffers();

   SQLFreeHandle(SQL_HANDLE_STMT, fHstmt);

   fHstmt = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// process statement

Bool_t TODBCStatement::Process()
{
   ClearError();

   SQLRETURN retcode = SQL_SUCCESS;

   if (IsParSettMode()) {

      // check if we start filling buffers, but not complete it
      if (fBufferCounter>=0) {
         // if buffer used not fully, set smaller size of buffer arrays
         if ((fBufferCounter>0) && (fBufferCounter<fBufferLength-1)) {
            SQLUINTEGER setsize = fBufferCounter+1;
            SQLSetStmtAttr(fHstmt, SQL_ATTR_PARAMSET_SIZE, (SQLPOINTER) (long) setsize, 0);
         }
         retcode = SQLExecute(fHstmt);
      }

      // after Process we finish working with parameters data,
      // if necessary, user can try to access resultset of statement
      fWorkingMode = 0;
      FreeBuffers();
      fBufferCounter = -1;
   } else {

      // just execute statement,
      // later one can try to access results of statement
      retcode = SQLExecute(fHstmt);
   }

   return !ExtractErrors(retcode, "Process");
}

////////////////////////////////////////////////////////////////////////////////
///get number of affected rows

Int_t TODBCStatement::GetNumAffectedRows()
{
   ClearError();

   SQLLEN    rowCount;
   SQLRETURN retcode = SQL_SUCCESS;

   retcode = SQLRowCount(fHstmt, &rowCount);

   if (ExtractErrors(retcode, "GetNumAffectedRows")) return -1;

   return rowCount;
}

////////////////////////////////////////////////////////////////////////////////
/// Store result of statement processing.
/// Results set, produced by processing of statement, can be stored, and accessed by
/// TODBCStamenet methoods like NextResultRow(), GetInt(), GetLong() and so on.

Bool_t TODBCStatement::StoreResult()
{
   ClearError();

   if (IsParSettMode()) {
      SetError(-1,"Call Process() method before","StoreResult");
      return kFALSE;
   }

   FreeBuffers();

   SQLSMALLINT columnCount = 0;

   SQLRETURN retcode = SQLNumResultCols(fHstmt, &columnCount);
   if (ExtractErrors(retcode, "StoreResult")) return kFALSE;

   if (columnCount==0) return kFALSE;

   SetNumBuffers(columnCount, fBufferPreferredSize);

   SQLULEN arrsize = fBufferLength;

   SQLSetStmtAttr(fHstmt, SQL_ATTR_ROW_BIND_TYPE, SQL_BIND_BY_COLUMN, 0);
   SQLSetStmtAttr(fHstmt, SQL_ATTR_ROW_ARRAY_SIZE, (SQLPOINTER) arrsize, 0);
   SQLSetStmtAttr(fHstmt, SQL_ATTR_ROW_STATUS_PTR, fStatusBuffer, 0);
   SQLSetStmtAttr(fHstmt, SQL_ATTR_ROWS_FETCHED_PTR, &fNumRowsFetched, 0);

   for (int n=0;n<fNumBuffers;n++) {
      SQLCHAR     columnName[1024];
      SQLSMALLINT nameLength;
      SQLSMALLINT dataType;
      SQLULEN     columnSize;
      SQLSMALLINT decimalDigits;
      SQLSMALLINT nullable;

      retcode = SQLDescribeCol(fHstmt, n+1, columnName, 1024,
                               &nameLength, &dataType,
                               &columnSize, &decimalDigits, &nullable);

      BindColumn(n, dataType, columnSize);

      if (nameLength>0) {
         fBuffer[n].fBnamebuffer = new char[nameLength+1];
         strlcpy(fBuffer[n].fBnamebuffer, (const char*) columnName, nameLength+1);
      }
   }

   fNumRowsFetched = 0;
   fLastResultRow = 0;

   fWorkingMode = 2;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///return number of fields

Int_t TODBCStatement::GetNumFields()
{
   return IsResultSet() ? fNumBuffers : -1;
}

////////////////////////////////////////////////////////////////////////////////
///return field name

const char* TODBCStatement::GetFieldName(Int_t nfield)
{
   ClearError();

   if (!IsResultSet() || (nfield<0) || (nfield>=fNumBuffers)) return 0;

   return fBuffer[nfield].fBnamebuffer;
}


////////////////////////////////////////////////////////////////////////////////
///next result row

Bool_t TODBCStatement::NextResultRow()
{
   ClearError();

   if (!IsResultSet()) return kFALSE;

   if ((fNumRowsFetched==0) ||
       (1.*fBufferCounter >= 1.*(fNumRowsFetched-1))) {

      fBufferCounter = 0;
      fNumRowsFetched = 0;

      SQLRETURN retcode = SQLFetchScroll(fHstmt, SQL_FETCH_NEXT, 0);
      if (retcode==SQL_NO_DATA) fNumRowsFetched=0; else
         ExtractErrors(retcode,"NextResultRow");

      // this is workaround of Oracle Linux ODBC driver
      // it does not returns number of fetched lines, therefore one should
      // calculate it from current row number
      if (!IsError() && (retcode!=SQL_NO_DATA) && (fNumRowsFetched==0)) {
         SQLULEN rownumber = 0;
         SQLRETURN retcode2 = SQLGetStmtAttr(fHstmt, SQL_ATTR_ROW_NUMBER, &rownumber, 0, 0);
         ExtractErrors(retcode2, "NextResultRow");

         if (!IsError()) {
            fNumRowsFetched = rownumber - fLastResultRow;
            fLastResultRow = rownumber;
         }
      }

      if (1.*fNumRowsFetched>fBufferLength)
         SetError(-1, "Missmatch between buffer length and fetched rows number", "NextResultRow");

      if (IsError() || (fNumRowsFetched==0)) {
         fWorkingMode = 0;
         FreeBuffers();
      }

   } else
      fBufferCounter++;

   return IsResultSet();
}

////////////////////////////////////////////////////////////////////////////////
/// Extract errors, produced by last ODBC function call

Bool_t TODBCStatement::ExtractErrors(SQLRETURN retcode, const char* method)
{
   if ((retcode== SQL_SUCCESS) || (retcode == SQL_SUCCESS_WITH_INFO)) return kFALSE;

   SQLINTEGER i = 0;
   SQLINTEGER native;
   SQLCHAR state[ 7 ];
   SQLCHAR text[256];
   SQLSMALLINT len;
   SQLRETURN ret;
   do {
      ret = SQLGetDiagRec(SQL_HANDLE_STMT, fHstmt, ++i, state, &native, text,
                          sizeof(text), &len );
      if (ret == SQL_SUCCESS) SetError(native, (const char*) text, method);
//         Error(method, "%s:%ld:%ld:%s\n", state, i, native, text);
   }
   while( ret == SQL_SUCCESS );
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///run next iteration

Bool_t TODBCStatement::NextIteration()
{
   ClearError();

   if (!IsParSettMode() || (fBuffer==0) || (fBufferLength<=0)) return kFALSE;

   if (fBufferCounter>=fBufferLength-1) {
      SQLRETURN retcode = SQLExecute(fHstmt);
      if (ExtractErrors(retcode,"NextIteration")) return kFALSE;
      fBufferCounter = 0;
   } else
      fBufferCounter++;

   // probably, we do not need it, but anyway
   fStatusBuffer[fBufferCounter] = SQL_ROW_SUCCESS;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///return number of parameters

Int_t TODBCStatement::GetNumParameters()
{
   return IsParSettMode() ? fNumBuffers : 0;
}

////////////////////////////////////////////////////////////////////////////////
///set number of buffers

void TODBCStatement::SetNumBuffers(Int_t isize, Int_t ilen)
{
   FreeBuffers();

   fNumBuffers = isize;
   fBufferLength = ilen;
   fBufferCounter = 0;

   fBuffer = new ODBCBufferRec_t[fNumBuffers];
   for (Int_t n=0;n<fNumBuffers;n++) {
      fBuffer[n].fBroottype = 0;
      fBuffer[n].fBsqltype = 0;
      fBuffer[n].fBsqlctype = 0;
      fBuffer[n].fBbuffer = nullptr;
      fBuffer[n].fBelementsize = 0;
      fBuffer[n].fBlenarray = 0;
      fBuffer[n].fBstrbuffer = 0;
      fBuffer[n].fBnamebuffer = 0;
   }

   fStatusBuffer = new SQLUSMALLINT[fBufferLength];
}

////////////////////////////////////////////////////////////////////////////////
/// Free allocated buffers

void TODBCStatement::FreeBuffers()
{
   if (fBuffer==0) return;
   for (Int_t n=0;n<fNumBuffers;n++) {
      if (fBuffer[n].fBbuffer)
        free(fBuffer[n].fBbuffer);
      delete[] fBuffer[n].fBlenarray;
      delete[] fBuffer[n].fBstrbuffer;
      delete[] fBuffer[n].fBnamebuffer;
   }

   delete[] fStatusBuffer;
   delete[] fBuffer;
   fBuffer = nullptr;
   fNumBuffers = 0;
   fBufferLength = 0;
   fStatusBuffer = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Bind result column to buffer. Allocate buffer of appropriate type

Bool_t TODBCStatement::BindColumn(Int_t ncol, SQLSMALLINT sqltype, SQLUINTEGER size)
{
   ClearError();

   if ((ncol<0) || (ncol>=fNumBuffers)) {
      SetError(-1,"Internal error. Column number invalid","BindColumn");
      return kFALSE;
   }

   if (fBuffer[ncol].fBsqltype!=0) {
      SetError(-1,"Internal error. Bind for column already done","BindColumn");
      return kFALSE;
   }

   SQLSMALLINT sqlctype = 0;
   switch (sqltype) {
      case SQL_CHAR:
      case SQL_VARCHAR:   sqlctype = SQL_C_CHAR; break;
      case SQL_BINARY:
      case SQL_LONGVARBINARY:
      case SQL_VARBINARY: sqlctype = SQL_C_BINARY; break;
      case SQL_LONGVARCHAR: Info("BindColumn","BIG VARCHAR not supported yet"); return kFALSE; break;

      case SQL_DECIMAL:   sqlctype = SQL_C_DOUBLE; break;
      case SQL_NUMERIC:   sqlctype = SQL_C_DOUBLE; break;
      case SQL_SMALLINT:  sqlctype = SQL_C_SLONG; break;
      case SQL_INTEGER:   sqlctype = SQL_C_SLONG; break;
      case SQL_FLOAT:     sqlctype = SQL_C_FLOAT; break;
      case SQL_REAL:
      case SQL_DOUBLE:    sqlctype = SQL_C_DOUBLE; break;
      case SQL_TINYINT:   sqlctype = SQL_C_STINYINT; break;
      case SQL_BIGINT:    sqlctype = SQL_C_SBIGINT; break;
      case SQL_TYPE_DATE: sqlctype = SQL_C_TYPE_DATE; break;
      case SQL_TYPE_TIME: sqlctype = SQL_C_TYPE_TIME; break;
      case SQL_TYPE_TIMESTAMP: sqlctype = SQL_C_TYPE_TIMESTAMP; break;
      default: {
         SetError(-1, Form("SQL type %d not supported",sqltype), "BindColumn");
         return kFALSE;
      }
   }

   int elemsize = 0;

   switch (sqlctype) {
      case SQL_C_ULONG:    elemsize = sizeof(SQLUINTEGER); break;
      case SQL_C_SLONG:    elemsize = sizeof(SQLINTEGER); break;
      case SQL_C_UBIGINT:  elemsize = sizeof(ULong64_t); break; // should be SQLUBIGINT, but it is 64-bit structure on some platforms
      case SQL_C_SBIGINT:  elemsize = sizeof(Long64_t); break; // should be SQLBIGINT, but it is 64-bit structure on some platforms
      case SQL_C_USHORT:   elemsize = sizeof(SQLUSMALLINT); break;
      case SQL_C_SSHORT:   elemsize = sizeof(SQLSMALLINT); break;
      case SQL_C_UTINYINT: elemsize = sizeof(SQLCHAR); break;
      case SQL_C_STINYINT: elemsize = sizeof(SQLSCHAR); break;
      case SQL_C_FLOAT:    elemsize = sizeof(SQLREAL); break;
      case SQL_C_DOUBLE:   elemsize = sizeof(SQLDOUBLE); break;
      case SQL_C_CHAR:     elemsize = size; break;
      case SQL_C_BINARY:   elemsize = size; break;
      case SQL_C_TYPE_DATE: elemsize = sizeof(DATE_STRUCT); break;
      case SQL_C_TYPE_TIME: elemsize = sizeof(TIME_STRUCT); break;
      case SQL_C_TYPE_TIMESTAMP: elemsize = sizeof(TIMESTAMP_STRUCT); break;

      default: {
         SetError(-1, Form("SQL C Type %d is not supported",sqlctype), "BindColumn");
         return kFALSE;
      }
   }

   fBuffer[ncol].fBroottype    = 0;
   fBuffer[ncol].fBsqltype     = sqltype;
   fBuffer[ncol].fBsqlctype    = sqlctype;
   fBuffer[ncol].fBbuffer      = malloc(elemsize * fBufferLength);
   fBuffer[ncol].fBelementsize = elemsize;
   fBuffer[ncol].fBlenarray    = new SQLLEN[fBufferLength];

   SQLRETURN retcode =
      SQLBindCol(fHstmt, ncol+1, sqlctype, fBuffer[ncol].fBbuffer,
                 elemsize,
                 fBuffer[ncol].fBlenarray);

   return !ExtractErrors(retcode, "BindColumn");
}

////////////////////////////////////////////////////////////////////////////////
/// Bind query parameter with buffer. Creates buffer of appropriate type

Bool_t TODBCStatement::BindParam(Int_t npar, Int_t roottype, Int_t size)
{
   ClearError();

   if ((npar<0) || (npar>=fNumBuffers)) return kFALSE;

   if (fBuffer[npar].fBroottype!=0) {
      SetError(-1,Form("ParameterType for par %d already specified", npar),"BindParam");
      return kFALSE;
   }

   SQLSMALLINT sqltype = 0, sqlctype = 0;
   int elemsize = 0;

   switch (roottype) {
      case kUInt_t:     sqltype = SQL_INTEGER; sqlctype = SQL_C_ULONG;    elemsize = sizeof(SQLUINTEGER); break;
      case kInt_t:      sqltype = SQL_INTEGER; sqlctype = SQL_C_SLONG;    elemsize = sizeof(SQLINTEGER); break;
      case kULong_t:    sqltype = SQL_INTEGER; sqlctype = SQL_C_ULONG;    elemsize = sizeof(SQLUINTEGER); break;
      case kLong_t:     sqltype = SQL_INTEGER; sqlctype = SQL_C_SLONG;    elemsize = sizeof(SQLINTEGER); break;

      // here SQLUBIGINT/SQLBIGINT types should be used,
       // but on 32-bit platforms it is structures, which makes its usage inconvinient
      case kULong64_t:  sqltype = SQL_BIGINT;  sqlctype = SQL_C_UBIGINT;  elemsize = sizeof(ULong64_t); break;
      case kLong64_t:   sqltype = SQL_BIGINT;  sqlctype = SQL_C_SBIGINT;  elemsize = sizeof(Long64_t); break;

      case kUShort_t:   sqltype = SQL_SMALLINT;sqlctype = SQL_C_USHORT;   elemsize = sizeof(SQLUSMALLINT); break;
      case kShort_t:    sqltype = SQL_SMALLINT;sqlctype = SQL_C_SSHORT;   elemsize = sizeof(SQLSMALLINT); break;
      case kUChar_t:    sqltype = SQL_TINYINT; sqlctype = SQL_C_UTINYINT; elemsize = sizeof(SQLCHAR); break;
      case kChar_t:     sqltype = SQL_TINYINT; sqlctype = SQL_C_STINYINT; elemsize = sizeof(SQLSCHAR); break;
      case kBool_t:     sqltype = SQL_TINYINT; sqlctype = SQL_C_UTINYINT; elemsize = sizeof(SQLCHAR); break;
      case kFloat_t:    sqltype = SQL_FLOAT;   sqlctype = SQL_C_FLOAT;    elemsize = sizeof(SQLREAL); break;
      case kFloat16_t:  sqltype = SQL_FLOAT;   sqlctype = SQL_C_FLOAT;    elemsize = sizeof(SQLREAL); break;
      case kDouble_t:   sqltype = SQL_DOUBLE;  sqlctype = SQL_C_DOUBLE;   elemsize = sizeof(SQLDOUBLE); break;
      case kDouble32_t: sqltype = SQL_DOUBLE;  sqlctype = SQL_C_DOUBLE;   elemsize = sizeof(SQLDOUBLE); break;
      case kCharStar:   sqltype = SQL_CHAR;    sqlctype = SQL_C_CHAR;     elemsize = size; break;
      case kSqlBinary:  sqltype = SQL_BINARY;  sqlctype = SQL_C_BINARY;   elemsize = size; break;
      case kSqlDate:    sqltype = SQL_TYPE_DATE; sqlctype = SQL_C_TYPE_DATE; elemsize = sizeof(DATE_STRUCT); break;
      case kSqlTime:    sqltype = SQL_TYPE_TIME; sqlctype = SQL_C_TYPE_TIME; elemsize = sizeof(TIME_STRUCT); break;
      case kSqlTimestamp: sqltype = SQL_TYPE_TIMESTAMP; sqlctype = SQL_C_TYPE_TIMESTAMP; elemsize = sizeof(TIMESTAMP_STRUCT); break;
      default: {
         SetError(-1, Form("Root type %d is not supported", roottype), "BindParam");
         return kFALSE;
      }
   }

   void* buffer = malloc(elemsize * fBufferLength);
   SQLLEN* lenarray = new SQLLEN[fBufferLength];
   SQLRETURN retcode =
      SQLBindParameter(fHstmt, npar+1, SQL_PARAM_INPUT,
                       sqlctype, sqltype, 0, 0,
                       buffer, elemsize, lenarray);

   if (ExtractErrors(retcode, "BindParam")) {
      free(buffer);
      delete[] lenarray;
      return kFALSE;
   }

   fBuffer[npar].fBroottype = roottype;
   fBuffer[npar].fBsqlctype = sqlctype;
   fBuffer[npar].fBsqltype = sqltype;
   fBuffer[npar].fBbuffer = buffer;
   fBuffer[npar].fBelementsize = elemsize;
   fBuffer[npar].fBlenarray = lenarray;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get parameter address

void* TODBCStatement::GetParAddr(Int_t npar, Int_t roottype, Int_t length)
{
   ClearError();

   if ((fBuffer==0) || (npar<0) || (npar>=fNumBuffers) || (fBufferCounter<0)) {
      SetError(-1, "Invalid parameter number","GetParAddr");
      return 0;
   }

   if (fBuffer[npar].fBbuffer==0) {
      if (IsParSettMode() && (roottype!=0) && (fBufferCounter==0))
         if (!BindParam(npar, roottype, length)) return 0;

      if (fBuffer[npar].fBbuffer==0) return 0;
   }

   if (roottype!=0)
      if (fBuffer[npar].fBroottype!=roottype) return 0;

   return (char*)fBuffer[npar].fBbuffer + fBufferCounter*fBuffer[npar].fBelementsize;
}

////////////////////////////////////////////////////////////////////////////////
///convert to numeric type

long double TODBCStatement::ConvertToNumeric(Int_t npar)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return 0;

   switch (fBuffer[npar].fBsqlctype) {
      case SQL_C_ULONG:    return *((SQLUINTEGER*) addr); break;
      case SQL_C_SLONG:    return *((SQLINTEGER*) addr); break;
      case SQL_C_UBIGINT:  return *((ULong64_t*) addr); break;
      case SQL_C_SBIGINT:  return *((Long64_t*) addr); break;
      case SQL_C_USHORT:   return *((SQLUSMALLINT*) addr); break;
      case SQL_C_SSHORT:   return *((SQLSMALLINT*) addr); break;
      case SQL_C_UTINYINT: return *((SQLCHAR*) addr); break;
      case SQL_C_STINYINT: return *((SQLSCHAR*) addr); break;
      case SQL_C_FLOAT:    return *((SQLREAL*) addr); break;
      case SQL_C_DOUBLE:   return *((SQLDOUBLE*) addr); break;
      case SQL_C_TYPE_DATE: {
         DATE_STRUCT* dt = (DATE_STRUCT*) addr;
         TDatime rtm(dt->year, dt->month,  dt->day, 0, 0, 0);
         return rtm.GetDate();
         break;
      }
      case SQL_C_TYPE_TIME: {
         TIME_STRUCT* tm = (TIME_STRUCT*) addr;
         TDatime rtm(2000, 1, 1, tm->hour, tm->minute, tm->second);
         return rtm.GetTime();
         break;
      }
      case SQL_C_TYPE_TIMESTAMP: {
         TIMESTAMP_STRUCT* tm = (TIMESTAMP_STRUCT*) addr;
         TDatime rtm(tm->year, tm->month,  tm->day,
                     tm->hour, tm->minute, tm->second);
         return rtm.Get();
         break;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///convert to string

const char* TODBCStatement::ConvertToString(Int_t npar)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return 0;
   if (fBuffer[npar].fBstrbuffer==0)
      fBuffer[npar].fBstrbuffer = new char[100];

   char* buf = fBuffer[npar].fBstrbuffer;

   switch(fBuffer[npar].fBsqlctype) {
#if (SIZEOF_LONG == 8)
      case SQL_C_SLONG:   snprintf(buf, 100, "%d", *((SQLINTEGER*) addr)); break;
      case SQL_C_ULONG:   snprintf(buf, 100, "%u", *((SQLUINTEGER*) addr)); break;
#else
      case SQL_C_SLONG:   snprintf(buf, 100, "%ld", (long)*((SQLINTEGER*) addr)); break;
      case SQL_C_ULONG:   snprintf(buf, 100, "%lu", (unsigned long)*((SQLUINTEGER*) addr)); break;
#endif
      case SQL_C_SBIGINT: snprintf(buf, 100, "%lld", *((Long64_t*) addr)); break;
      case SQL_C_UBIGINT: snprintf(buf, 100, "%llu", *((ULong64_t*) addr)); break;
      case SQL_C_SSHORT:  snprintf(buf, 100, "%hd", *((SQLSMALLINT*) addr)); break;
      case SQL_C_USHORT:  snprintf(buf, 100, "%hu", *((SQLUSMALLINT*) addr)); break;
      case SQL_C_STINYINT:snprintf(buf, 100, "%d", *((SQLSCHAR*) addr)); break;
      case SQL_C_UTINYINT:snprintf(buf, 100, "%u", *((SQLCHAR*) addr)); break;
      case SQL_C_FLOAT:   snprintf(buf, 100, TSQLServer::GetFloatFormat(), *((SQLREAL*) addr)); break;
      case SQL_C_DOUBLE:  snprintf(buf, 100, TSQLServer::GetFloatFormat(), *((SQLDOUBLE*) addr)); break;
      case SQL_C_TYPE_DATE: {
         DATE_STRUCT* dt = (DATE_STRUCT*) addr;
         snprintf(buf,100,"%4.4d-%2.2d-%2.2d",
                  dt->year, dt->month,  dt->day);
         break;
      }
      case SQL_C_TYPE_TIME: {
         TIME_STRUCT* tm = (TIME_STRUCT*) addr;
         snprintf(buf,100,"%2.2d:%2.2d:%2.2d",
                  tm->hour, tm->minute, tm->second);
         break;
      }
      case SQL_C_TYPE_TIMESTAMP: {
         TIMESTAMP_STRUCT* tm = (TIMESTAMP_STRUCT*) addr;
         snprintf(buf,100,"%4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d",
                  tm->year, tm->month,  tm->day,
                  tm->hour, tm->minute, tm->second);
         break;
      }
      default: return 0;
   }

   return buf;
}

////////////////////////////////////////////////////////////////////////////////
/// Verifies if field value is NULL

Bool_t TODBCStatement::IsNull(Int_t npar)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return kTRUE;

   return fBuffer[npar].fBlenarray[fBufferCounter] == SQL_NULL_DATA;
}

////////////////////////////////////////////////////////////////////////////////
///get parameter as integer

Int_t TODBCStatement::GetInt(Int_t npar)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return 0;

   if (fBuffer[npar].fBsqlctype==SQL_C_SLONG)
      return (Int_t) *((SQLINTEGER*) addr);

   return (Int_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
///get parameter as unsigned integer

UInt_t TODBCStatement::GetUInt(Int_t npar)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return 0;

   if (fBuffer[npar].fBsqlctype==SQL_C_ULONG)
      return (UInt_t) *((SQLUINTEGER*) addr);

   return (UInt_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
///get parameter as Long_t

Long_t TODBCStatement::GetLong(Int_t npar)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return 0;

   if (fBuffer[npar].fBsqlctype==SQL_C_SLONG)
     return (Long_t) *((SQLINTEGER*) addr);

   return (Long_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
///get parameter as Long64_t

Long64_t TODBCStatement::GetLong64(Int_t npar)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return 0;

   if (fBuffer[npar].fBsqlctype==SQL_C_SBIGINT)
     return *((Long64_t*) addr);

   return (Long64_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
///get parameter as ULong64_t

ULong64_t TODBCStatement::GetULong64(Int_t npar)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return 0;

   if (fBuffer[npar].fBsqlctype==SQL_C_UBIGINT)
     return *((ULong64_t*) addr);

   return (ULong64_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
///get parameter as Double_t

Double_t TODBCStatement::GetDouble(Int_t npar)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return 0;

   if (fBuffer[npar].fBsqlctype==SQL_C_DOUBLE)
     return *((SQLDOUBLE*) addr);

   return (Double_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
///get parameter as string

const char* TODBCStatement::GetString(Int_t npar)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return 0;

   if (fBuffer[npar].fBsqlctype==SQL_C_CHAR) {
      // first check if string is null

      int len = fBuffer[npar].fBlenarray[fBufferCounter];

      if ((len == SQL_NULL_DATA) || (len==0)) return 0;

      char* res = (char*) addr;
      if (len < fBuffer[npar].fBelementsize) {
         *(res + len) = 0;
         return res;
      }

      if (len > fBuffer[npar].fBelementsize) {
         SetError(-1, Form("Problems with string size %d", len), "GetString");
         return 0;
      }

      if (fBuffer[npar].fBstrbuffer==0)
         fBuffer[npar].fBstrbuffer = new char[len+1];

      strlcpy(fBuffer[npar].fBstrbuffer, res, len+1);

      res = fBuffer[npar].fBstrbuffer;
      *(res + len) = 0;
      return res;
   }

   return ConvertToString(npar);
}

////////////////////////////////////////////////////////////////////////////////
/// return parameter as binary data

Bool_t TODBCStatement::GetBinary(Int_t npar, void* &mem, Long_t& size)
{
   mem = 0;
   size = 0;

   void* addr = GetParAddr(npar);
   if (addr==0) return kFALSE;

   if ((fBuffer[npar].fBsqlctype==SQL_C_BINARY) ||
       (fBuffer[npar].fBsqlctype==SQL_C_CHAR)) {

      // first check if data length is null
      int len = fBuffer[npar].fBlenarray[fBufferCounter];

      if ((len == SQL_NULL_DATA) || (len==0)) return kTRUE;

      size = len;

      if (fBuffer[npar].fBstrbuffer==0)
         fBuffer[npar].fBstrbuffer = new char[size];

      memcpy(fBuffer[npar].fBstrbuffer, addr, size);

      mem = fBuffer[npar].fBstrbuffer;

      return kTRUE;
   }

   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// return field value as date

Bool_t TODBCStatement::GetDate(Int_t npar, Int_t& year, Int_t& month, Int_t& day)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return kFALSE;

   if (fBuffer[npar].fBsqlctype!=SQL_C_TYPE_DATE) return kFALSE;

   DATE_STRUCT* dt = (DATE_STRUCT*) addr;
   year = dt->year;
   month = dt->month;
   day = dt->day;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as time

Bool_t TODBCStatement::GetTime(Int_t npar, Int_t& hour, Int_t& min, Int_t& sec)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return kFALSE;

   if (fBuffer[npar].fBsqlctype!=SQL_C_TYPE_TIME) return kFALSE;

   TIME_STRUCT* tm = (TIME_STRUCT*) addr;
   hour = tm->hour;
   min = tm->minute;
   sec = tm->second;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as date & time

Bool_t TODBCStatement::GetDatime(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return kFALSE;

   if (fBuffer[npar].fBsqlctype!=SQL_C_TYPE_TIMESTAMP) return kFALSE;

   TIMESTAMP_STRUCT* tm = (TIMESTAMP_STRUCT*) addr;

   year = tm->year;
   month = tm->month;
   day = tm->day;
   hour = tm->hour;
   min = tm->minute;
   sec = tm->second;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as time stamp

Bool_t TODBCStatement::GetTimestamp(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec, Int_t& frac)
{
   void* addr = GetParAddr(npar);
   if (addr==0) return kFALSE;

   if (fBuffer[npar].fBsqlctype!=SQL_C_TYPE_TIMESTAMP) return kFALSE;

   TIMESTAMP_STRUCT* tm = (TIMESTAMP_STRUCT*) addr;

   year = tm->year;
   month = tm->month;
   day = tm->day;
   hour = tm->hour;
   min = tm->minute;
   sec = tm->second;
   frac = tm->fraction;
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Set NULL as parameter value
/// If NULL should be set for statement parameter during first iteration,
/// one should call before proper Set... method to identify type of argument for
/// the future. For instance, if one suppose to have double as type of parameter,
/// code should look like:
///    stmt->SetDouble(2, 0.);
///    stmt->SetNull(2);

Bool_t TODBCStatement::SetNull(Int_t npar)
{
   void* addr = GetParAddr(npar, kInt_t);
   if (addr!=0)
      *((SQLINTEGER*) addr) = 0;

   if ((npar>=0) && (npar<fNumBuffers))
      fBuffer[npar].fBlenarray[fBufferCounter] = SQL_NULL_DATA;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///set parameter as Int_t

Bool_t TODBCStatement::SetInt(Int_t npar, Int_t value)
{
   void* addr = GetParAddr(npar, kInt_t);
   if (addr==0) return kFALSE;

   *((SQLINTEGER*) addr) = value;

   fBuffer[npar].fBlenarray[fBufferCounter] = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///set parameter as UInt_t

Bool_t TODBCStatement::SetUInt(Int_t npar, UInt_t value)
{
   void* addr = GetParAddr(npar, kUInt_t);
   if (addr==0) return kFALSE;

   *((SQLUINTEGER*) addr) = value;

   fBuffer[npar].fBlenarray[fBufferCounter] = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///set parameter as Long_t

Bool_t TODBCStatement::SetLong(Int_t npar, Long_t value)
{
   void* addr = GetParAddr(npar, kLong_t);
   if (addr==0) return kFALSE;

   *((SQLINTEGER*) addr) = value;

   fBuffer[npar].fBlenarray[fBufferCounter] = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///set parameter as Long64_t

Bool_t TODBCStatement::SetLong64(Int_t npar, Long64_t value)
{
   void* addr = GetParAddr(npar, kLong64_t);
   if (addr==0) return kFALSE;

   *((Long64_t*) addr) = value;

   fBuffer[npar].fBlenarray[fBufferCounter] = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///set parameter as ULong64_t

Bool_t TODBCStatement::SetULong64(Int_t npar, ULong64_t value)
{
   void* addr = GetParAddr(npar, kULong64_t);
   if (addr==0) return kFALSE;

   *((ULong64_t*) addr) = value;

   fBuffer[npar].fBlenarray[fBufferCounter] = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///set parameter as Double_t

Bool_t TODBCStatement::SetDouble(Int_t npar, Double_t value)
{
   void* addr = GetParAddr(npar, kDouble_t);
   if (addr==0) return kFALSE;

   *((SQLDOUBLE*) addr) = value;

   fBuffer[npar].fBlenarray[fBufferCounter] = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///set parameter as string

Bool_t TODBCStatement::SetString(Int_t npar, const char* value, Int_t maxsize)
{
   void* addr = GetParAddr(npar, kCharStar, maxsize);

   if (addr==0) return kFALSE;

   if (value) {
      int len = strlen(value);

      if (len>=fBuffer[npar].fBelementsize) {
         len = fBuffer[npar].fBelementsize;
         strlcpy((char*) addr, value, len+1);
         fBuffer[npar].fBlenarray[fBufferCounter] = len;

      } else if (len>0) {
         strlcpy((char*) addr, value, maxsize);
         fBuffer[npar].fBlenarray[fBufferCounter] = SQL_NTS;
      } else {
         *((char*) addr) = 0;
         fBuffer[npar].fBlenarray[fBufferCounter] = SQL_NTS;
      }
   } else {
      *((char*) addr) = 0;
      fBuffer[npar].fBlenarray[fBufferCounter] = SQL_NTS;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///set parameter value as binary data

Bool_t TODBCStatement::SetBinary(Int_t npar, void* mem, Long_t size, Long_t maxsize)
{
   void* addr = GetParAddr(npar, kSqlBinary, maxsize);
   if (addr==0) return kFALSE;

   if (size>fBuffer[npar].fBelementsize)
      size = fBuffer[npar].fBelementsize;

   memcpy(addr, mem, size);
   fBuffer[npar].fBlenarray[fBufferCounter] = size;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// set parameter value as date

Bool_t TODBCStatement::SetDate(Int_t npar, Int_t year, Int_t month, Int_t day)
{
   void* addr = GetParAddr(npar, kSqlDate);
   if (addr==0) return kFALSE;

   DATE_STRUCT* dt = (DATE_STRUCT*) addr;
   dt->year = year;
   dt->month = month;
   dt->day = day;

   fBuffer[npar].fBlenarray[fBufferCounter] = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// set parameter value as time

Bool_t TODBCStatement::SetTime(Int_t npar, Int_t hour, Int_t min, Int_t sec)
{
   void* addr = GetParAddr(npar, kSqlTime);
   if (addr==0) return kFALSE;

   TIME_STRUCT* tm = (TIME_STRUCT*) addr;
   tm->hour = hour;
   tm->minute = min;
   tm->second = sec;

   fBuffer[npar].fBlenarray[fBufferCounter] = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// set parameter value as date & time

Bool_t TODBCStatement::SetDatime(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec)
{
   void* addr = GetParAddr(npar, kSqlTimestamp);
   if (addr==0) return kFALSE;

   TIMESTAMP_STRUCT* tm = (TIMESTAMP_STRUCT*) addr;
   tm->year = year;
   tm->month = month;
   tm->day = day;
   tm->hour = hour;
   tm->minute = min;
   tm->second = sec;
   tm->fraction = 0;

   fBuffer[npar].fBlenarray[fBufferCounter] = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// set parameter value as timestamp

Bool_t TODBCStatement::SetTimestamp(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec, Int_t frac)
{
   void* addr = GetParAddr(npar, kSqlTimestamp);
   if (addr==0) return kFALSE;

   TIMESTAMP_STRUCT* tm = (TIMESTAMP_STRUCT*) addr;
   tm->year = year;
   tm->month = month;
   tm->day = day;
   tm->hour = hour;
   tm->minute = min;
   tm->second = sec;
   tm->fraction = frac;

   fBuffer[npar].fBlenarray[fBufferCounter] = 0;

   return kTRUE;
}
