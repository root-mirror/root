// @(#)root/mysql:$Id$
// Author: Sergey Linev   6/02/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  SQL statement class for MySQL                                       //
//                                                                      //
//  See TSQLStatement class documentation for more details.             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMySQLStatement.h"
#include "TMySQLServer.h"
#include "TDataType.h"
#include "TDatime.h"
#include "snprintf.h"
#include <cstdlib>

ClassImp(TMySQLStatement);

ULong64_t TMySQLStatement::fgAllocSizeLimit = 0x8000000; // 128 Mb

#if MYSQL_VERSION_ID >= 40100

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.
/// Checks if statement contains parameters tags.

TMySQLStatement::TMySQLStatement(MYSQL_STMT* stmt, Bool_t errout) :
   TSQLStatement(errout),
   fStmt(stmt)
{
   ULong_t paramcount = mysql_stmt_param_count(fStmt);

   if (paramcount>0) {
      fWorkingMode = 1;
      SetBuffersNumber(paramcount);
      fNeedParBind = kTRUE;
      fIterationCount = -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TMySQLStatement::~TMySQLStatement()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close statement.

void TMySQLStatement::Close(Option_t *)
{
   if (fStmt)
      mysql_stmt_close(fStmt);

   fStmt = nullptr;

   FreeBuffers();
}


// Reset error and check that statement exists
#define CheckStmt(method, res)                          \
   {                                                    \
      ClearError();                                     \
      if (fStmt==0) {                                   \
         SetError(-1,"Statement handle is 0",method);   \
         return res;                                    \
      }                                                 \
   }

// check last mysql statement error code
#define CheckErrNo(method, force, res)                  \
   {                                                    \
      unsigned int stmterrno = mysql_stmt_errno(fStmt);     \
      if ((stmterrno!=0) || force) {                        \
         const char* stmterrmsg = mysql_stmt_error(fStmt);  \
         if (stmterrno==0) { stmterrno = 11111; stmterrmsg = "MySQL statement error"; } \
         SetError(stmterrno, stmterrmsg, method);               \
         return res;                                    \
      }                                                 \
   }


// check last mysql statement error code
#define CheckGetField(method, res)                      \
   {                                                    \
      ClearError();                                     \
      if (!IsResultSetMode()) {                         \
         SetError(-1,"Cannot get statement parameters",method); \
         return res;                                    \
      }                                                 \
      if ((npar<0) || (npar>=fNumBuffers)) {            \
         SetError(-1,Form("Invalid parameter number %d", npar),method); \
         return res;                                    \
      }                                                 \
   }

////////////////////////////////////////////////////////////////////////////////
/// Process statement.

Bool_t TMySQLStatement::Process()
{
   CheckStmt("Process",kFALSE);

   // if parameters was set, processing just means of closing parameters and variables
   if (IsSetParsMode()) {
      if (fIterationCount>=0)
         if (!NextIteration()) return kFALSE;
      fWorkingMode = 0;
      fIterationCount = -1;
      FreeBuffers();
      return kTRUE;
   }

   if (mysql_stmt_execute(fStmt))
      CheckErrNo("Process",kTRUE, kFALSE);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of affected rows after statement is processed.

Int_t TMySQLStatement::GetNumAffectedRows()
{
   CheckStmt("Process", -1);

   my_ulonglong res = mysql_stmt_affected_rows(fStmt);

   if (res == (my_ulonglong) -1)
      CheckErrNo("GetNumAffectedRows", kTRUE, -1);

   return (Int_t) res;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of statement parameters.

Int_t TMySQLStatement::GetNumParameters()
{
   CheckStmt("GetNumParameters", -1);

   Int_t res = mysql_stmt_param_count(fStmt);

   CheckErrNo("GetNumParameters", kFALSE, -1);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Store result of statement processing to access them
/// via GetInt(), GetDouble() and so on methods.

Bool_t TMySQLStatement::StoreResult()
{
   CheckStmt("StoreResult", kFALSE);
   if (fWorkingMode!=0) {
      SetError(-1,"Cannot store result for that statement","StoreResult");
      return kFALSE;
   }

   if (mysql_stmt_store_result(fStmt))
      CheckErrNo("StoreResult",kTRUE, kFALSE);

   // allocate memeory for data reading from query
   MYSQL_RES* meta = mysql_stmt_result_metadata(fStmt);
   if (meta) {
      int count = mysql_num_fields(meta);

      SetBuffersNumber(count);

      MYSQL_FIELD *fields = mysql_fetch_fields(meta);

      for (int n=0;n<count;n++) {
         SetSQLParamType(n, fields[n].type, (fields[n].flags & UNSIGNED_FLAG) == 0, fields[n].length);
         if (fields[n].name)
            fBuffer[n].fFieldName = fields[n].name;
      }

      mysql_free_result(meta);
   }

   if (!fBind) return kFALSE;

   /* Bind the buffers */
   if (mysql_stmt_bind_result(fStmt, fBind))
      CheckErrNo("StoreResult",kTRUE, kFALSE);

   fWorkingMode = 2;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of fields in result set.

Int_t TMySQLStatement::GetNumFields()
{
   return IsResultSetMode() ? fNumBuffers : -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns field name in result set.

const char* TMySQLStatement::GetFieldName(Int_t nfield)
{
   if (!IsResultSetMode() || (nfield<0) || (nfield>=fNumBuffers)) return nullptr;

   return fBuffer[nfield].fFieldName.empty() ? nullptr : fBuffer[nfield].fFieldName.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// Shift cursor to nect row in result set.

Bool_t TMySQLStatement::NextResultRow()
{
   if ((fStmt==0) || !IsResultSetMode()) return kFALSE;

   Bool_t res = !mysql_stmt_fetch(fStmt);

   if (!res) {
      fWorkingMode = 0;
      FreeBuffers();
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment iteration counter for statement, where parameter can be set.
/// Statement with parameters of previous iteration
/// automatically will be applied to database.

Bool_t TMySQLStatement::NextIteration()
{
   ClearError();

   if (!IsSetParsMode() || (fBind==0)) {
      SetError(-1,"Cannot call for that statement","NextIteration");
      return kFALSE;
   }

   fIterationCount++;

   if (fIterationCount==0) return kTRUE;

   if (fNeedParBind) {
      fNeedParBind = kFALSE;
      if (mysql_stmt_bind_param(fStmt, fBind))
         CheckErrNo("NextIteration",kTRUE, kFALSE);
   }

   if (mysql_stmt_execute(fStmt))
      CheckErrNo("NextIteration", kTRUE, kFALSE);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Release all buffers, used by statement.

void TMySQLStatement::FreeBuffers()
{
   if (fBuffer) {
      for (Int_t n=0; n<fNumBuffers;n++) {
         free(fBuffer[n].fMem);
      }
      delete[] fBuffer;
   }

   if (fBind)
      delete[] fBind;

   fBuffer = nullptr;
   fBind = nullptr;
   fNumBuffers = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate buffers for statement parameters/ result fields.

void TMySQLStatement::SetBuffersNumber(Int_t numpars)
{
   FreeBuffers();
   if (numpars<=0) return;

   fNumBuffers = numpars;

   fBind = new MYSQL_BIND[fNumBuffers];
   memset(fBind, 0, sizeof(MYSQL_BIND)*fNumBuffers);

   fBuffer = new TParamData[fNumBuffers];
   for (int n=0;n<fNumBuffers;++n) {
      fBuffer[n].fMem = nullptr;
      fBuffer[n].fSize = 0;
      fBuffer[n].fSqlType = 0;
      fBuffer[n].fSign = kFALSE;
      fBuffer[n].fResLength = 0;
      fBuffer[n].fResNull = false;
      fBuffer[n].fStrBuffer.clear();
      fBuffer[n].fFieldName.clear();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Convert field value to string.

const char *TMySQLStatement::ConvertToString(Int_t npar)
{
   if (fBuffer[npar].fResNull)
      return nullptr;

   void *addr = fBuffer[npar].fMem;
   Bool_t sig = fBuffer[npar].fSign;

   if (!addr)
      return nullptr;

   if ((fBind[npar].buffer_type==MYSQL_TYPE_STRING) ||
      (fBind[npar].buffer_type==MYSQL_TYPE_VAR_STRING))
      return (const char *) addr;

   constexpr int kSize = 100;
   char buf[kSize];
   int len = 0;

   switch(fBind[npar].buffer_type) {
      case MYSQL_TYPE_LONG:
         if (sig) len = snprintf(buf, kSize, "%d",*((int*) addr));
             else len = snprintf(buf, kSize, "%u",*((unsigned int*) addr));
         break;
      case MYSQL_TYPE_LONGLONG:
         if (sig) len = snprintf(buf, kSize, "%lld",*((Long64_t*) addr)); else
                  len = snprintf(buf, kSize, "%llu",*((ULong64_t*) addr));
         break;
      case MYSQL_TYPE_SHORT:
         if (sig) len = snprintf(buf, kSize, "%hd",*((short*) addr)); else
                  len = snprintf(buf, kSize, "%hu",*((unsigned short*) addr));
         break;
      case MYSQL_TYPE_TINY:
         if (sig) len = snprintf(buf, kSize, "%d",*((char*) addr)); else
                  len = snprintf(buf, kSize, "%u",*((unsigned char*) addr));
         break;
      case MYSQL_TYPE_FLOAT:
         len = snprintf(buf, kSize, TSQLServer::GetFloatFormat(), *((float*) addr));
         break;
      case MYSQL_TYPE_DOUBLE:
         len = snprintf(buf, kSize, TSQLServer::GetFloatFormat(), *((double*) addr));
         break;
      case MYSQL_TYPE_DATETIME:
      case MYSQL_TYPE_TIMESTAMP: {
         MYSQL_TIME* tm = (MYSQL_TIME*) addr;
         len = snprintf(buf, kSize, "%4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d",
                            tm->year, tm->month,  tm->day,
                            tm->hour, tm->minute, tm->second);
         break;
      }
      case MYSQL_TYPE_TIME: {
         MYSQL_TIME* tm = (MYSQL_TIME*) addr;
         len = snprintf(buf, kSize, "%2.2d:%2.2d:%2.2d",
                             tm->hour, tm->minute, tm->second);
         break;
      }
      case MYSQL_TYPE_DATE: {
         MYSQL_TIME* tm = (MYSQL_TIME*) addr;
         len = snprintf(buf, kSize, "%4.4d-%2.2d-%2.2d",
                             tm->year, tm->month,  tm->day);
         break;
      }
      default:
         return nullptr;
   }

   if (len >= kSize)
      SetError(-1, Form("Cannot convert param %d into string - buffer too small", npar));

   fBuffer[npar].fStrBuffer = buf;

   return fBuffer[npar].fStrBuffer.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// Convert field to numeric value.

long double TMySQLStatement::ConvertToNumeric(Int_t npar)
{
   if (fBuffer[npar].fResNull) return 0;

   void* addr = fBuffer[npar].fMem;
   Bool_t sig = fBuffer[npar].fSign;

   if (addr==0) return 0;

   switch(fBind[npar].buffer_type) {
      case MYSQL_TYPE_LONG:
         if (sig) return *((int*) addr); else
                  return *((unsigned int*) addr);
         break;
      case MYSQL_TYPE_LONGLONG:
         if (sig) return *((Long64_t*) addr); else
                  return *((ULong64_t*) addr);
         break;
      case MYSQL_TYPE_SHORT:
         if (sig) return *((short*) addr); else
                  return *((unsigned short*) addr);
         break;
      case MYSQL_TYPE_TINY:
         if (sig) return *((char*) addr); else
                  return *((unsigned char*) addr);
         break;
      case MYSQL_TYPE_FLOAT:
         return *((float*) addr);
         break;
      case MYSQL_TYPE_DOUBLE:
         return *((double*) addr);
         break;
#if MYSQL_VERSION_ID >= 50022
      case MYSQL_TYPE_NEWDECIMAL /* new MYSQL_TYPE fixed precision decimal */:
#endif
      case MYSQL_TYPE_STRING:
      case MYSQL_TYPE_VAR_STRING:
      case MYSQL_TYPE_BLOB: {
         char* str = (char*) addr;
         ULong_t len = fBuffer[npar].fResLength;
         if ((str==0) || (*str==0) || (len==0)) return 0;
         Int_t size = fBuffer[npar].fSize;
         if (1.*len<size)
            str[len] = 0;
         else
            str[size-1] = 0;
         long double buf = 0;
         sscanf(str,"%Lf",&buf);
         return buf;
         break;
      }
      case MYSQL_TYPE_DATETIME:
      case MYSQL_TYPE_TIMESTAMP: {
         MYSQL_TIME* tm = (MYSQL_TIME*) addr;
         TDatime rtm(tm->year, tm->month,  tm->day,
                  tm->hour, tm->minute, tm->second);
         return rtm.Get();
         break;
      }
      case MYSQL_TYPE_DATE: {
         MYSQL_TIME* tm = (MYSQL_TIME*) addr;
         TDatime rtm(tm->year, tm->month,  tm->day, 0, 0, 0);
         return rtm.GetDate();
         break;
      }
      case MYSQL_TYPE_TIME: {
         MYSQL_TIME* tm = (MYSQL_TIME*) addr;
         TDatime rtm(2000, 1, 1, tm->hour, tm->minute, tm->second);
         return rtm.GetTime();
         break;
      }

      default:
         return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if field value is null.

Bool_t TMySQLStatement::IsNull(Int_t npar)
{
   CheckGetField("IsNull", kTRUE);

   return fBuffer[npar].fResNull;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as integer.

Int_t TMySQLStatement::GetInt(Int_t npar)
{
   CheckGetField("GetInt", 0);

   if ((fBuffer[npar].fSqlType==MYSQL_TYPE_LONG) && fBuffer[npar].fSign)
     return (Int_t) *((int*) fBuffer[npar].fMem);

   return (Int_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as unsigned integer.

UInt_t TMySQLStatement::GetUInt(Int_t npar)
{
   CheckGetField("GetUInt", 0);

   if ((fBuffer[npar].fSqlType==MYSQL_TYPE_LONG) && !fBuffer[npar].fSign)
     return (UInt_t) *((unsigned int*) fBuffer[npar].fMem);

   return (UInt_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as long integer.

Long_t TMySQLStatement::GetLong(Int_t npar)
{
   CheckGetField("GetLong", 0);

   if ((fBuffer[npar].fSqlType==MYSQL_TYPE_LONG) && fBuffer[npar].fSign)
     return (Long_t) *((int*) fBuffer[npar].fMem);

   return (Long_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as 64-bit integer.

Long64_t TMySQLStatement::GetLong64(Int_t npar)
{
   CheckGetField("GetLong64", 0);

   if ((fBuffer[npar].fSqlType==MYSQL_TYPE_LONGLONG) && fBuffer[npar].fSign)
     return (Long64_t) *((Long64_t*) fBuffer[npar].fMem);

   return (Long64_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as unsigned 64-bit integer.

ULong64_t TMySQLStatement::GetULong64(Int_t npar)
{
   CheckGetField("GetULong64", 0);

   if ((fBuffer[npar].fSqlType==MYSQL_TYPE_LONGLONG) && !fBuffer[npar].fSign)
     return (ULong64_t) *((ULong64_t*) fBuffer[npar].fMem);

   return (ULong64_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as double.

Double_t TMySQLStatement::GetDouble(Int_t npar)
{
   CheckGetField("GetDouble", 0);

   if (fBuffer[npar].fSqlType==MYSQL_TYPE_DOUBLE)
     return (Double_t) *((double*) fBuffer[npar].fMem);

   return (Double_t) ConvertToNumeric(npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as string.

const char *TMySQLStatement::GetString(Int_t npar)
{
   CheckGetField("GetString", 0);

   if ((fBind[npar].buffer_type==MYSQL_TYPE_STRING)
      || (fBind[npar].buffer_type==MYSQL_TYPE_BLOB)
      || (fBind[npar].buffer_type==MYSQL_TYPE_VAR_STRING)
#if MYSQL_VERSION_ID >= 50022
      || (fBuffer[npar].fSqlType==MYSQL_TYPE_NEWDECIMAL)
#endif
       ) {
         if (fBuffer[npar].fResNull) return nullptr;
         char *str = (char *) fBuffer[npar].fMem;
         ULong_t len = fBuffer[npar].fResLength;
         Int_t size = fBuffer[npar].fSize;
         if (1.*len<size) str[len] = 0; else
                          str[size-1] = 0;
         return str;
      }

   return ConvertToString(npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as binary array.

Bool_t TMySQLStatement::GetBinary(Int_t npar, void* &mem, Long_t& size)
{
   mem = 0;
   size = 0;

   CheckGetField("GetBinary", kFALSE);

   if ((fBind[npar].buffer_type==MYSQL_TYPE_STRING) ||
       (fBind[npar].buffer_type==MYSQL_TYPE_VAR_STRING) ||
       (fBind[npar].buffer_type==MYSQL_TYPE_BLOB) ||
       (fBind[npar].buffer_type==MYSQL_TYPE_TINY_BLOB) ||
       (fBind[npar].buffer_type==MYSQL_TYPE_MEDIUM_BLOB) ||
       (fBind[npar].buffer_type==MYSQL_TYPE_LONG_BLOB)) {
         if (fBuffer[npar].fResNull) return kTRUE;
         mem = fBuffer[npar].fMem;
         size = fBuffer[npar].fResLength;
         return kTRUE;
      }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as date.

Bool_t TMySQLStatement::GetDate(Int_t npar, Int_t& year, Int_t& month, Int_t& day)
{
   CheckGetField("GetDate", kFALSE);

   if (fBuffer[npar].fResNull) return kFALSE;

   switch(fBind[npar].buffer_type) {
      case MYSQL_TYPE_DATETIME:
      case MYSQL_TYPE_TIMESTAMP:
      case MYSQL_TYPE_DATE: {
         MYSQL_TIME* tm = (MYSQL_TIME*) fBuffer[npar].fMem;
         if (tm==0) return kFALSE;
         year = tm->year;
         month = tm->month;
         day = tm->day;
         break;
      }
      default:
         return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as time.

Bool_t TMySQLStatement::GetTime(Int_t npar, Int_t& hour, Int_t& min, Int_t& sec)
{
   CheckGetField("GetTime", kFALSE);

   if (fBuffer[npar].fResNull) return kFALSE;

   switch(fBind[npar].buffer_type) {
      case MYSQL_TYPE_DATETIME:
      case MYSQL_TYPE_TIMESTAMP:
      case MYSQL_TYPE_TIME: {
         MYSQL_TIME* tm = (MYSQL_TIME*) fBuffer[npar].fMem;
         if (tm==0) return kFALSE;
         hour = tm->hour;
         min = tm->minute;
         sec = tm->second;
         break;
      }
      default:
         return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as date & time.

Bool_t TMySQLStatement::GetDatime(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec)
{
   CheckGetField("GetDatime", kFALSE);

   if (fBuffer[npar].fResNull) return kFALSE;

   switch(fBind[npar].buffer_type) {
      case MYSQL_TYPE_DATETIME:
      case MYSQL_TYPE_TIMESTAMP: {
         MYSQL_TIME* tm = (MYSQL_TIME*) fBuffer[npar].fMem;
         if (tm==0) return kFALSE;
         year = tm->year;
         month = tm->month;
         day = tm->day;
         hour = tm->hour;
         min = tm->minute;
         sec = tm->second;
         break;
      }
      default:
         return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as time stamp.

Bool_t TMySQLStatement::GetTimestamp(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec, Int_t& frac)
{
   CheckGetField("GetTimstamp", kFALSE);

   if (fBuffer[npar].fResNull) return kFALSE;

   switch(fBind[npar].buffer_type) {
      case MYSQL_TYPE_DATETIME:
      case MYSQL_TYPE_TIMESTAMP: {
         MYSQL_TIME* tm = (MYSQL_TIME*) fBuffer[npar].fMem;
         if (tm==0) return kFALSE;
         year = tm->year;
         month = tm->month;
         day = tm->day;
         hour = tm->hour;
         min = tm->minute;
         sec = tm->second;
         frac = 0;
         break;
      }
      default:
         return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter type to be used as buffer.
/// Used in both setting data to database and retrieving data from data base.
/// Initialize proper MYSQL_BIND structure and allocate required buffers.

Bool_t TMySQLStatement::SetSQLParamType(Int_t npar, int sqltype, Bool_t sig, ULong_t sqlsize)
{
   if ((npar<0) || (npar>=fNumBuffers)) return kFALSE;

   fBuffer[npar].fMem = nullptr;
   fBuffer[npar].fSize = 0;
   fBuffer[npar].fResLength = 0;
   fBuffer[npar].fResNull = false;
   fBuffer[npar].fStrBuffer.clear();

   ULong64_t allocsize = 0;

   Bool_t doreset = false;

   switch (sqltype) {
      case MYSQL_TYPE_LONG:     allocsize = sizeof(int);  break;
      case MYSQL_TYPE_LONGLONG: allocsize = sizeof(Long64_t); break;
      case MYSQL_TYPE_SHORT:    allocsize = sizeof(short); break;
      case MYSQL_TYPE_TINY:     allocsize = sizeof(char); break;
      case MYSQL_TYPE_FLOAT:    allocsize = sizeof(float); break;
      case MYSQL_TYPE_DOUBLE:   allocsize = sizeof(double); break;
#if MYSQL_VERSION_ID >= 50022
      case MYSQL_TYPE_NEWDECIMAL /* new MYSQL_TYPE fixed precision decimal */:
#endif
      case MYSQL_TYPE_STRING:   allocsize = sqlsize > 256 ? sqlsize : 256; break;
      case MYSQL_TYPE_VAR_STRING: allocsize = sqlsize > 256 ? sqlsize : 256; break;
      case MYSQL_TYPE_MEDIUM_BLOB:
      case MYSQL_TYPE_LONG_BLOB:
      case MYSQL_TYPE_BLOB:     allocsize = sqlsize >= 65525 ? sqlsize : 65535; break;
      case MYSQL_TYPE_TINY_BLOB:   allocsize = sqlsize > 255 ? sqlsize : 255; break;
      case MYSQL_TYPE_TIME:
      case MYSQL_TYPE_DATE:
      case MYSQL_TYPE_TIMESTAMP:
      case MYSQL_TYPE_DATETIME: allocsize = sizeof(MYSQL_TIME); doreset = true; break;
      default: SetError(-1,"Nonsupported SQL type","SetSQLParamType"); return kFALSE;
   }

   if (allocsize > fgAllocSizeLimit) allocsize = fgAllocSizeLimit;

   fBuffer[npar].fMem = malloc(allocsize);
   fBuffer[npar].fSize = allocsize;
   fBuffer[npar].fSqlType = sqltype;
   fBuffer[npar].fSign = sig;

   if ((allocsize>0) && fBuffer[npar].fMem && doreset)
      memset(fBuffer[npar].fMem, 0, allocsize);

   fBind[npar].buffer_type = enum_field_types(sqltype);
   fBind[npar].buffer = fBuffer[npar].fMem;
   fBind[npar].buffer_length = allocsize;
   fBind[npar].is_null= &(fBuffer[npar].fResNull);
   fBind[npar].length = &(fBuffer[npar].fResLength);
   fBind[npar].is_unsigned = !sig;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check boundary condition before setting value of parameter.
/// Return address of parameter buffer.

void *TMySQLStatement::BeforeSet(const char* method, Int_t npar, Int_t sqltype, Bool_t sig, ULong_t size)
{
   ClearError();

   if (!IsSetParsMode()) {
      SetError(-1,"Cannot set parameter for statement", method);
      return 0;
   }

   if ((npar<0) || (npar>=fNumBuffers)) {
      SetError(-1,Form("Invalid parameter number %d",npar), method);
      return 0;
   }

   if ((fIterationCount==0) && (fBuffer[npar].fSqlType==0))
      if (!SetSQLParamType(npar, sqltype, sig, size)) {
         SetError(-1,"Cannot initialize parameter buffer", method);
         return 0;
      }

   if ((fBuffer[npar].fSqlType!=sqltype) ||
      (fBuffer[npar].fSign != sig)) return 0;

   fBuffer[npar].fResNull = false;

   return fBuffer[npar].fMem;
}

////////////////////////////////////////////////////////////////////////////////
/// Set NULL as parameter value.
/// If NULL should be set for statement parameter during first iteration,
/// one should call before proper Set... method to identify type of argument for
/// the future. For instance, if one suppose to have double as type of parameter,
/// code should look like:
///    stmt->SetDouble(2, 0.);
///    stmt->SetNull(2);

Bool_t TMySQLStatement::SetNull(Int_t npar)
{
   void* addr = BeforeSet("SetNull", npar, MYSQL_TYPE_LONG);

   if (addr!=0)
      *((int*) addr) = 0;

   if ((npar>=0) && (npar<fNumBuffers))
      fBuffer[npar].fResNull = true;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as integer.

Bool_t TMySQLStatement::SetInt(Int_t npar, Int_t value)
{
   void* addr = BeforeSet("SetInt", npar, MYSQL_TYPE_LONG);

   if (addr!=0)
      *((int*) addr) = value;

   return (addr!=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as unsigned integer.

Bool_t TMySQLStatement::SetUInt(Int_t npar, UInt_t value)
{
   void* addr = BeforeSet("SetUInt", npar, MYSQL_TYPE_LONG, kFALSE);

   if (addr!=0)
      *((unsigned int*) addr) = value;

   return (addr!=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as long integer.

Bool_t TMySQLStatement::SetLong(Int_t npar, Long_t value)
{
   void* addr = BeforeSet("SetLong", npar, MYSQL_TYPE_LONG);

   if (addr!=0)
      *((int*) addr) = value;

   return (addr!=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as 64-bit integer.

Bool_t TMySQLStatement::SetLong64(Int_t npar, Long64_t value)
{
   void* addr = BeforeSet("SetLong64", npar, MYSQL_TYPE_LONGLONG);

   if (addr!=0)
      *((Long64_t*) addr) = value;

   return (addr!=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as unsigned 64-bit integer.

Bool_t TMySQLStatement::SetULong64(Int_t npar, ULong64_t value)
{
   void* addr = BeforeSet("SetULong64", npar, MYSQL_TYPE_LONGLONG, kFALSE);

   if (addr!=0)
      *((ULong64_t*) addr) = value;

   return (addr!=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as double.

Bool_t TMySQLStatement::SetDouble(Int_t npar, Double_t value)
{
   void* addr = BeforeSet("SetDouble", npar, MYSQL_TYPE_DOUBLE, kFALSE);

   if (addr!=0)
      *((double*) addr) = value;

   return (addr!=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as string.

Bool_t TMySQLStatement::SetString(Int_t npar, const char* value, Int_t maxsize)
{
   Int_t len = value ? strlen(value) : 0;

   void* addr = BeforeSet("SetString", npar, MYSQL_TYPE_STRING, true, maxsize);

   if (addr==0) return kFALSE;

   if (len >= fBuffer[npar].fSize) {
      free(fBuffer[npar].fMem);

      fBuffer[npar].fMem = malloc(len+1);
      fBuffer[npar].fSize = len + 1;

      fBind[npar].buffer = fBuffer[npar].fMem;
      fBind[npar].buffer_length = fBuffer[npar].fSize;

      addr = fBuffer[npar].fMem;
      fNeedParBind = kTRUE;
   }

   if (value) strcpy((char*) addr, value);
   else ((char*)addr)[0]='\0';

   fBuffer[npar].fResLength = len;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as binary data.

Bool_t TMySQLStatement::SetBinary(Int_t npar, void* mem, Long_t size, Long_t maxsize)
{
   if (size>=maxsize) maxsize = size + 1;

   int bin_type = MYSQL_TYPE_BLOB;
   if (maxsize > 65525) bin_type = MYSQL_TYPE_MEDIUM_BLOB;
   if (maxsize > 16777205) bin_type = MYSQL_TYPE_LONG_BLOB;

   void* addr = BeforeSet("SetBinary", npar, bin_type, true, maxsize);

   if (addr==0) return kFALSE;

   if (size >= fBuffer[npar].fSize) {
      free(fBuffer[npar].fMem);

      fBuffer[npar].fMem = malloc(size+1);
      fBuffer[npar].fSize = size + 1;

      fBind[npar].buffer = fBuffer[npar].fMem;
      fBind[npar].buffer_length = fBuffer[npar].fSize;

      addr = fBuffer[npar].fMem;
      fNeedParBind = kTRUE;
   }

   memcpy(addr, mem, size);

   fBuffer[npar].fResLength = size;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as date.

Bool_t TMySQLStatement::SetDate(Int_t npar, Int_t year, Int_t month, Int_t day)
{
   MYSQL_TIME* addr = (MYSQL_TIME*) BeforeSet("SetDate", npar, MYSQL_TYPE_DATE);

   if (addr!=0) {
      addr->year = year;
      addr->month = month;
      addr->day = day;
   }

   return (addr!=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as time.

Bool_t TMySQLStatement::SetTime(Int_t npar, Int_t hour, Int_t min, Int_t sec)
{
   MYSQL_TIME* addr = (MYSQL_TIME*) BeforeSet("SetTime", npar, MYSQL_TYPE_TIME);

   if (addr!=0) {
      addr->hour = hour;
      addr->minute = min;
      addr->second = sec;
   }

   return (addr!=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as date & time.

Bool_t TMySQLStatement::SetDatime(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec)
{
   MYSQL_TIME* addr = (MYSQL_TIME*) BeforeSet("SetDatime", npar, MYSQL_TYPE_DATETIME);

   if (addr!=0) {
      addr->year = year;
      addr->month = month;
      addr->day = day;
      addr->hour = hour;
      addr->minute = min;
      addr->second = sec;
   }

   return (addr!=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as timestamp.

Bool_t TMySQLStatement::SetTimestamp(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec, Int_t)
{
   MYSQL_TIME* addr = (MYSQL_TIME*) BeforeSet("SetTimestamp", npar, MYSQL_TYPE_TIMESTAMP);

   if (addr!=0) {
      addr->year = year;
      addr->month = month;
      addr->day = day;
      addr->hour = hour;
      addr->minute = min;
      addr->second = sec;
   }

   return (addr!=0);
}

#else

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.
/// For MySQL version < 4.1 no statement is supported

TMySQLStatement::TMySQLStatement(MYSQL_STMT*, Bool_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TMySQLStatement::~TMySQLStatement()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Close statement

void TMySQLStatement::Close(Option_t *)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Process statement.

Bool_t TMySQLStatement::Process()
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of affected rows after statement is processed.

Int_t TMySQLStatement::GetNumAffectedRows()
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of statement parameters.

Int_t TMySQLStatement::GetNumParameters()
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Store result of statement processing to access them
/// via GetInt(), GetDouble() and so on methods.

Bool_t TMySQLStatement::StoreResult()
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of fields in result set.

Int_t TMySQLStatement::GetNumFields()
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns field name in result set.

const char* TMySQLStatement::GetFieldName(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Shift cursor to nect row in result set.

Bool_t TMySQLStatement::NextResultRow()
{
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment iteration counter for statement, where parameter can be set.
/// Statement with parameters of previous iteration
/// automatically will be applied to database.

Bool_t TMySQLStatement::NextIteration()
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Release all buffers, used by statement.

void TMySQLStatement::FreeBuffers()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate buffers for statement parameters/ result fields.

void TMySQLStatement::SetBuffersNumber(Int_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Convert field value to string.

const char* TMySQLStatement::ConvertToString(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert field to numeric value.

long double TMySQLStatement::ConvertToNumeric(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if field value is null.

Bool_t TMySQLStatement::IsNull(Int_t)
{
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as integer.

Int_t TMySQLStatement::GetInt(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as unsigned integer.

UInt_t TMySQLStatement::GetUInt(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as long integer.

Long_t TMySQLStatement::GetLong(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as 64-bit integer.

Long64_t TMySQLStatement::GetLong64(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as unsigned 64-bit integer.

ULong64_t TMySQLStatement::GetULong64(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as double.

Double_t TMySQLStatement::GetDouble(Int_t)
{
   return 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as string.

const char *TMySQLStatement::GetString(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as binary array.

Bool_t TMySQLStatement::GetBinary(Int_t, void* &, Long_t&)
{
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Return field value as date.

Bool_t TMySQLStatement::GetDate(Int_t, Int_t&, Int_t&, Int_t&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as time.

Bool_t TMySQLStatement::GetTime(Int_t, Int_t&, Int_t&, Int_t&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as date & time.

Bool_t TMySQLStatement::GetDatime(Int_t, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as time stamp.

Bool_t TMySQLStatement::GetTimestamp(Int_t, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter type to be used as buffer.
/// Used in both setting data to database and retriving data from data base.
/// Initialize proper MYSQL_BIND structure and allocate required buffers.

Bool_t TMySQLStatement::SetSQLParamType(Int_t, int, Bool_t, ULong_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check boundary condition before setting value of parameter.
/// Return address of parameter buffer.

void *TMySQLStatement::BeforeSet(const char*, Int_t, Int_t, Bool_t, ULong_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set NULL as parameter value.
/// If NULL should be set for statement parameter during first iteration,
/// one should call before proper Set... method to identify type of argument for
/// the future. For instance, if one suppose to have double as type of parameter,
/// code should look like:
///    stmt->SetDouble(2, 0.);
///    stmt->SetNull(2);

Bool_t TMySQLStatement::SetNull(Int_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as integer.

Bool_t TMySQLStatement::SetInt(Int_t, Int_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as unsigned integer.

Bool_t TMySQLStatement::SetUInt(Int_t, UInt_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as long integer.

Bool_t TMySQLStatement::SetLong(Int_t, Long_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as 64-bit integer.

Bool_t TMySQLStatement::SetLong64(Int_t, Long64_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as unsigned 64-bit integer.

Bool_t TMySQLStatement::SetULong64(Int_t, ULong64_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as double.

Bool_t TMySQLStatement::SetDouble(Int_t, Double_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as string.

Bool_t TMySQLStatement::SetString(Int_t, const char*, Int_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as binary data.

Bool_t TMySQLStatement::SetBinary(Int_t, void*, Long_t, Long_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as date.

Bool_t TMySQLStatement::SetDate(Int_t, Int_t, Int_t, Int_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as time.

Bool_t TMySQLStatement::SetTime(Int_t, Int_t, Int_t, Int_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as date & time.

Bool_t TMySQLStatement::SetDatime(Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as timestamp.

Bool_t TMySQLStatement::SetTimestamp(Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t)
{
   return kFALSE;
}

#endif // MYSQL_VERSION_ID > 40100
