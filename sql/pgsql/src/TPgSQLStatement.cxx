   // @(#)root/pgsql:$Id$
// Author: Dennis Box (dbox@fnal.gov)  3/12/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  SQL statement class for PgSQL                                       //
//                                                                      //
//  See TSQLStatement class documentation for more details.             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TPgSQLStatement.h"
#include "TDataType.h"
#include "TDatime.h"
#include "TTimeStamp.h"

#include <stdlib.h>

ClassImp(TPgSQLStatement)

#ifdef PG_VERSION_NUM

#include "libpq/libpq-fs.h"

static const Int_t kBindStringSize = 25;

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.
/// Checks if statement contains parameters tags.

TPgSQLStatement::TPgSQLStatement(PgSQL_Stmt_t* stmt, Bool_t errout):
   TSQLStatement(errout),
   fStmt(stmt),
   fNumBuffers(0),
   fBind(0),
   fFieldName(0),
   fWorkingMode(0),
   fIterationCount(0),
   fParamLengths(0),
   fParamFormats(0),
   fNumResultRows(0),
   fNumResultCols(0)
{
   // Given fRes not used, we retrieve the statement using the connection.
   if (fStmt->fRes != NULL) {
      PQclear(fStmt->fRes);
   }

   fStmt->fRes = PQdescribePrepared(fStmt->fConn,"preparedstmt");
   unsigned long paramcount = PQnparams(fStmt->fRes);
   fNumResultCols = PQnfields(fStmt->fRes);
   fIterationCount = -1;

   if (paramcount>0) {
      fWorkingMode = 1;
      SetBuffersNumber(paramcount);
   } else {
      fWorkingMode = 2;
      SetBuffersNumber(fNumResultCols);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TPgSQLStatement::~TPgSQLStatement()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close statement.

void TPgSQLStatement::Close(Option_t *)
{
   if (fStmt->fRes)
      PQclear(fStmt->fRes);

   fStmt->fRes = 0;

   PGresult *res=PQexec(fStmt->fConn,"DEALLOCATE preparedstmt;");
   PQclear(res);

   FreeBuffers();
   //TPgSQLServers responsibility to free connection
   fStmt->fConn=0;
   delete fStmt;
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

#define CheckErrNo(method, force, wtf)                  \
   {                                                    \
      int stmterrno = PQresultStatus(fStmt->fRes);      \
      if ((stmterrno!=0) || force) {                        \
         const char* stmterrmsg = PQresultErrorMessage(fStmt->fRes);  \
         if (stmterrno==0) { stmterrno = -1; stmterrmsg = "PgSQL statement error"; } \
         SetError(stmterrno, stmterrmsg, method);               \
         return wtf;                                    \
      }                                                 \
   }

#define CheckErrResult(method, pqresult, retVal)         \
   {                                                     \
      ExecStatusType stmterrno=PQresultStatus(pqresult); \
      if (!pgsql_success(stmterrno)) {                   \
       const char* stmterrmsg = PQresultErrorMessage(fStmt->fRes);  \
       SetError(stmterrno, stmterrmsg, method);                     \
       PQclear(res);                                     \
       return retVal;                                    \
     }                                                   \
   }

#define RollBackTransaction(method)                          \
   {                                                         \
      PGresult *resnum=PQexec(fStmt->fConn,"COMMIT");        \
      CheckErrResult("RollBackTransaction", resnum, kFALSE); \
      PQclear(res);                                          \
   }

// check last pgsql statement error code
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

Bool_t TPgSQLStatement::Process()
{
   CheckStmt("Process",kFALSE);

   // We create the prepared statement below, MUST delete the old one
   // from our constructor first!
   if (fStmt->fRes != NULL) {
      PQclear(fStmt->fRes);
   }

   if (IsSetParsMode()) {
      fStmt->fRes= PQexecPrepared(fStmt->fConn,"preparedstmt",fNumBuffers,
                                 (const char* const*)fBind,
                                 0,0,0);

   } else { //result set mode

      fStmt->fRes= PQexecPrepared(fStmt->fConn,"preparedstmt",0,(const char* const*) 0,0,0,0);
   }
   ExecStatusType stat = PQresultStatus(fStmt->fRes);
   if (!pgsql_success(stat))
      CheckErrNo("Process",kTRUE, kFALSE);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of affected rows after statement is processed.

Int_t TPgSQLStatement::GetNumAffectedRows()
{
   CheckStmt("GetNumAffectedRows", -1);

   return (Int_t) atoi(PQcmdTuples(fStmt->fRes));
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of statement parameters.

Int_t TPgSQLStatement::GetNumParameters()
{
   CheckStmt("GetNumParameters", -1);

   if (IsSetParsMode()) {
      return fNumBuffers;
   } else {
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Store result of statement processing to access them
/// via GetInt(), GetDouble() and so on methods.

Bool_t TPgSQLStatement::StoreResult()
{
   int i;
   for (i=0;i<fNumResultCols;i++){
      fFieldName[i] = PQfname(fStmt->fRes,i);
      fParamFormats[i]=PQftype(fStmt->fRes,i);
      fParamLengths[i]=PQfsize(fStmt->fRes,i);

   }
   fNumResultRows=PQntuples(fStmt->fRes);
   ExecStatusType stat = PQresultStatus(fStmt->fRes);
   fWorkingMode = 2;
   if (!pgsql_success(stat))
      CheckErrNo("StoreResult",kTRUE, kFALSE);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of fields in result set.

Int_t TPgSQLStatement::GetNumFields()
{
   if (fWorkingMode==1)
      return fNumBuffers;
   if (fWorkingMode==2)
      return fNumResultCols;
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns field name in result set.

const char* TPgSQLStatement::GetFieldName(Int_t nfield)
{
   if (!IsResultSetMode() || (nfield<0) || (nfield>=fNumBuffers)) return 0;

   return fFieldName[nfield];
}

////////////////////////////////////////////////////////////////////////////////
/// Shift cursor to nect row in result set.

Bool_t TPgSQLStatement::NextResultRow()
{
   if ((fStmt==0) || !IsResultSetMode()) return kFALSE;

   Bool_t res=kTRUE;

   fIterationCount++;
   if (fIterationCount>=fNumResultRows)
     res=kFALSE;
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment iteration counter for statement, where parameter can be set.
/// Statement with parameters of previous iteration
/// automatically will be applied to database.

Bool_t TPgSQLStatement::NextIteration()
{
   ClearError();

   if (!IsSetParsMode() || (fBind==0)) {
      SetError(-1,"Cannot call for that statement","NextIteration");
      return kFALSE;
   }

   fIterationCount++;

   if (fIterationCount==0) return kTRUE;

   fStmt->fRes= PQexecPrepared(fStmt->fConn,"preparedstmt",fNumBuffers,
                               (const char* const*)fBind,
                               0,//fParamLengths,
                               0,//fParamFormats,
                               0);
   ExecStatusType stat = PQresultStatus(fStmt->fRes);
   if (!pgsql_success(stat) ){
      CheckErrNo("NextIteration", kTRUE, kFALSE) ;
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Release all buffers, used by statement.

void TPgSQLStatement::FreeBuffers()
{
  //individual field names free()'ed by PQclear of fStmt->fRes
   if (fFieldName)
      delete[] fFieldName;

   if (fBind){
      for (Int_t i=0;i<fNumBuffers;i++)
         delete [] fBind[i];
      delete[] fBind;
   }

   if (fParamLengths)
      delete [] fParamLengths;

   if (fParamFormats)
      delete [] fParamFormats;

   fFieldName = 0;
   fBind = 0;
   fNumBuffers = 0;
   fParamLengths = 0;
   fParamFormats = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate buffers for statement parameters/ result fields.

void TPgSQLStatement::SetBuffersNumber(Int_t numpars)
{
   FreeBuffers();
   if (numpars<=0) return;

   fNumBuffers = numpars;

   fBind = new char*[fNumBuffers];
   for(int i=0; i<fNumBuffers; ++i){
      fBind[i] = new char[kBindStringSize]; //big enough to handle text rep. of 64 bit number
   }
   fFieldName = new char*[fNumBuffers];

   fParamLengths = new int[fNumBuffers];
   memset(fParamLengths, 0, sizeof(int)*fNumBuffers);

   fParamFormats = new int[fNumBuffers];
   memset(fParamFormats, 0, sizeof(int)*fNumBuffers);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert field value to string.

const char* TPgSQLStatement::ConvertToString(Int_t npar)
{
   const char *buf = PQgetvalue(fStmt->fRes, fIterationCount, npar);
   return buf;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert field to numeric.

long double TPgSQLStatement::ConvertToNumeric(Int_t npar)
{
   if (PQgetisnull(fStmt->fRes,fIterationCount,npar))
      return (long double)0;

   return (long double) atof(PQgetvalue(fStmt->fRes,fIterationCount,npar));
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if field value is null.

Bool_t TPgSQLStatement::IsNull(Int_t npar)
{
   CheckGetField("IsNull", kTRUE);

   return PQgetisnull(fStmt->fRes,fIterationCount,npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Get integer.

Int_t TPgSQLStatement::GetInt(Int_t npar)
{
   if (PQgetisnull(fStmt->fRes,fIterationCount,npar))
      return (Int_t)0;

   return (Int_t) atoi(PQgetvalue(fStmt->fRes,fIterationCount,npar));
}

////////////////////////////////////////////////////////////////////////////////
/// Get unsigned integer.

UInt_t TPgSQLStatement::GetUInt(Int_t npar)
{
   if (PQgetisnull(fStmt->fRes,fIterationCount,npar))
      return (UInt_t)0;

   return (UInt_t) atoi(PQgetvalue(fStmt->fRes,fIterationCount,npar));
}

////////////////////////////////////////////////////////////////////////////////
/// Get long.

Long_t TPgSQLStatement::GetLong(Int_t npar)
{
   if (PQgetisnull(fStmt->fRes,fIterationCount,npar))
      return (Long_t)0;

   return (Long_t) atol(PQgetvalue(fStmt->fRes,fIterationCount,npar));
}

////////////////////////////////////////////////////////////////////////////////
/// Get long64.

Long64_t TPgSQLStatement::GetLong64(Int_t npar)
{
   if (PQgetisnull(fStmt->fRes,fIterationCount,npar))
      return (Long64_t)0;

#ifndef R__WIN32
   return (Long64_t) atoll(PQgetvalue(fStmt->fRes,fIterationCount,npar));
#else
   return (Long64_t) _atoi64(PQgetvalue(fStmt->fRes,fIterationCount,npar));
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as unsigned 64-bit integer

ULong64_t TPgSQLStatement::GetULong64(Int_t npar)
{
   if (PQgetisnull(fStmt->fRes,fIterationCount,npar))
      return (ULong64_t)0;

#ifndef R__WIN32
   return (ULong64_t) atoll(PQgetvalue(fStmt->fRes,fIterationCount,npar));
#else
   return (ULong64_t) _atoi64(PQgetvalue(fStmt->fRes,fIterationCount,npar));
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as double.

Double_t TPgSQLStatement::GetDouble(Int_t npar)
{
   if (PQgetisnull(fStmt->fRes,fIterationCount,npar))
      return (Double_t)0;
   return (Double_t) atof(PQgetvalue(fStmt->fRes,fIterationCount,npar));
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as string.

const char *TPgSQLStatement::GetString(Int_t npar)
{
   return PQgetvalue(fStmt->fRes,fIterationCount,npar);
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as binary array.
/// Note PQgetvalue mallocs/frees and ROOT classes expect new/delete.

Bool_t TPgSQLStatement::GetBinary(Int_t npar, void* &mem, Long_t& size)
{
   size_t sz;
   char *cptr = PQgetvalue(fStmt->fRes,fIterationCount,npar);
   unsigned char * mptr = PQunescapeBytea((const unsigned char*)cptr,&sz);
   if ((Long_t)sz>size) {
      delete [] (unsigned char*) mem;
      mem = (void*) new unsigned char[sz];
   }
   size=sz;
   memcpy(mem,mptr,sz);
   PQfreemem(mptr);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return large object whose oid is in the given field.

Bool_t TPgSQLStatement::GetLargeObject(Int_t npar, void* &mem, Long_t& size)
{
   Int_t objID = atoi(PQgetvalue(fStmt->fRes,fIterationCount,npar));

   // All this needs to happen inside a transaction, or it will NOT work.
   PGresult *res=PQexec(fStmt->fConn,"BEGIN");

   CheckErrResult("GetLargeObject", res, kFALSE);
   PQclear(res);

   Int_t lObjFD = lo_open(fStmt->fConn, objID, INV_READ);

   if (lObjFD<0) {
      Error("GetLargeObject", "SQL Error on lo_open: %s", PQerrorMessage(fStmt->fConn));
      RollBackTransaction("GetLargeObject");
      return kFALSE;
   }
   // Object size is not known beforehand.
   // Possible fast ways to get it are:
   // (1) Create a function that does fopen, fseek, ftell on server
   // (2) Query large object table with size()
   // Both can not be expected to work in general,
   // as  (1) needs permissions and changes DB,
   // and (2) needs permission.
   // So we use
   // (3) fopen, fseek and ftell locally.

   lo_lseek(fStmt->fConn, lObjFD, 0, SEEK_END);
   Long_t sz = lo_tell(fStmt->fConn, lObjFD);
   lo_lseek(fStmt->fConn, lObjFD, 0, SEEK_SET);

   if ((Long_t)sz>size) {
      delete [] (unsigned char*) mem;
      mem = (void*) new unsigned char[sz];
      size=sz;
   }

   Int_t readBytes = lo_read(fStmt->fConn, lObjFD, (char*)mem, size);

   if (readBytes != sz) {
      Error("GetLargeObject", "SQL Error on lo_read: %s", PQerrorMessage(fStmt->fConn));
      RollBackTransaction("GetLargeObject");
      return kFALSE;
   }

   if (lo_close(fStmt->fConn, lObjFD) != 0) {
      Error("GetLargeObject", "SQL Error on lo_close: %s", PQerrorMessage(fStmt->fConn));
      RollBackTransaction("GetLargeObject");
      return kFALSE;
   }

   res=PQexec(fStmt->fConn,"COMMIT");

   ExecStatusType stat = PQresultStatus(res);
   if (!pgsql_success(stat)) {
      Error("GetLargeObject", "SQL Error on COMMIT: %s", PQerrorMessage(fStmt->fConn));
      RollBackTransaction("GetLargeObject");
      return kFALSE;
   }
   PQclear(res);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as date, in UTC.

Bool_t TPgSQLStatement::GetDate(Int_t npar, Int_t& year, Int_t& month, Int_t& day)
{
   TString val=PQgetvalue(fStmt->fRes,fIterationCount,npar);
   TDatime d = TDatime(val.Data());
   year = d.GetYear();
   month = d.GetMonth();
   day= d.GetDay();
   Int_t hour = d.GetHour();
   Int_t min = d.GetMinute();
   Int_t sec = d.GetSecond();
   ConvertTimeToUTC(val, year, month, day, hour, min, sec);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field as time, in UTC.

Bool_t TPgSQLStatement::GetTime(Int_t npar, Int_t& hour, Int_t& min, Int_t& sec)
{
   TString val=PQgetvalue(fStmt->fRes,fIterationCount,npar);
   TDatime d = TDatime(val.Data());
   hour = d.GetHour();
   min = d.GetMinute();
   sec= d.GetSecond();
   Int_t year = d.GetYear();
   Int_t month = d.GetMonth();
   Int_t day = d.GetDay();
   ConvertTimeToUTC(val, day, month, year, hour, min, sec);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as date & time, in UTC.

Bool_t TPgSQLStatement::GetDatime(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec)
{
   TString val=PQgetvalue(fStmt->fRes,fIterationCount,npar);
   TDatime d = TDatime(val.Data());
   year = d.GetYear();
   month = d.GetMonth();
   day= d.GetDay();
   hour = d.GetHour();
   min = d.GetMinute();
   sec= d.GetSecond();
   ConvertTimeToUTC(val, year, month, day, hour, min, sec);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert timestamp value to UTC if a zone is included.

void TPgSQLStatement::ConvertTimeToUTC(const TString &PQvalue, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec)
{
   Ssiz_t p = PQvalue.Last('.');
   // Check if timestamp has timezone
   TSubString *s_zone = nullptr;
   Bool_t hasZone = kFALSE;
   Ssiz_t tzP = PQvalue.Last('+');
   if ((tzP != kNPOS) && (tzP > p) ) {
      s_zone = new TSubString(PQvalue(tzP+1,PQvalue.Length()-tzP));
      hasZone=kTRUE;
   } else {
      Ssiz_t tzM = PQvalue.Last('-');
      if ((tzM != kNPOS) && (tzM > p) ) {
         s_zone = new TSubString(PQvalue(tzM+1,PQvalue.Length()-tzM));
         hasZone = kTRUE;
      }
   }
   if (hasZone == kTRUE) {
      // Parse timezone, might look like e.g. +00 or -00:00
      Int_t hourOffset, minuteOffset = 0;
      Int_t conversions=sscanf(s_zone->Data(), "%2d:%2d", &hourOffset, &minuteOffset);
      Int_t secondOffset = hourOffset*3600;
      if (conversions>1) {
         // Use sign from hour also for minute
         secondOffset += (TMath::Sign(minuteOffset, hourOffset))*60;
      }
      // Use TTimeStamp so we do not have to take care of over-/underflows
      TTimeStamp ts(year, month, day, hour, min, sec, 0, kTRUE, -secondOffset);
      UInt_t uyear, umonth, uday, uhour, umin, usec;
      ts.GetDate(kTRUE, 0, &uyear, &umonth, &uday);
      ts.GetTime(kTRUE, 0, &uhour, &umin, &usec);
      year=uyear;
      month=umonth;
      day=uday;
      hour=uhour;
      min=umin;
      sec=usec;
      delete s_zone;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return field as timestamp, in UTC.
/// Second fraction is to be interpreted as in the following example:
/// 2013-01-12 12:10:23.093854+02
/// Fraction is '93854', precision is fixed in this method to 6 decimal places.
/// This means the returned frac-value is always in microseconds.

Bool_t TPgSQLStatement::GetTimestamp(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec, Int_t& frac)
{
   TString val=PQgetvalue(fStmt->fRes,fIterationCount,npar);
   TDatime d(val.Data());
   year = d.GetYear();
   month = d.GetMonth();
   day= d.GetDay();
   hour = d.GetHour();
   min = d.GetMinute();
   sec= d.GetSecond();

   ConvertTimeToUTC(val, year, month, day, hour, min, sec);

   Ssiz_t p = val.Last('.');
   TSubString s_frac = val(p,val.Length()-p+1);

   // atoi ignores timezone part.
   // We MUST use atof here to correctly convert the fraction of
   // "12:23:01.093854" and put a limitation on precision,
   // as we can only return an Int_t.
   frac=(Int_t) (atof(s_frac.Data())*1.E6);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set NULL as parameter value.
/// If NULL should be set for statement parameter during first iteration,
/// one should call before proper Set... method to identify type of argument for
/// the future. For instance, if one suppose to have double as type of parameter,
/// code should look like:
///    stmt->SetDouble(2, 0.);
///    stmt->SetNull(2);

Bool_t TPgSQLStatement::SetNull(Int_t npar)
{
   fBind[npar] = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as integer.

Bool_t TPgSQLStatement::SetInt(Int_t npar, Int_t value)
{
   snprintf(fBind[npar],kBindStringSize,"%d",value);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as unsinged integer.

Bool_t TPgSQLStatement::SetUInt(Int_t npar, UInt_t value)
{
   snprintf(fBind[npar],kBindStringSize,"%u",value);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as long.

Bool_t TPgSQLStatement::SetLong(Int_t npar, Long_t value)
{
   snprintf(fBind[npar],kBindStringSize,"%ld",value);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as 64-bit integer.

Bool_t TPgSQLStatement::SetLong64(Int_t npar, Long64_t value)
{
   snprintf(fBind[npar],kBindStringSize,"%lld",(Long64_t)value);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as unsinged 64-bit integer.

Bool_t TPgSQLStatement::SetULong64(Int_t npar, ULong64_t value)
{
   snprintf(fBind[npar],kBindStringSize,"%llu",(ULong64_t)value);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as double value.

Bool_t TPgSQLStatement::SetDouble(Int_t npar, Double_t value)
{
   snprintf(fBind[npar],kBindStringSize,"%lf",value);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as string.

Bool_t TPgSQLStatement::SetString(Int_t npar, const char* value, Int_t maxsize)
{
   if(sizeof(fBind[npar])<(unsigned)maxsize){
      delete [] fBind[npar];
      fBind[npar] = new char[maxsize];
   }
   strlcpy(fBind[npar],value,maxsize);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as binary data.

Bool_t TPgSQLStatement::SetBinary(Int_t npar, void* mem, Long_t size, Long_t maxsize)
{
   // Set parameter value as binary data.

   size_t sz = size, mxsz = maxsize;
   unsigned char* escape_ptr = PQescapeBytea((const unsigned char*)mem, sz, &mxsz);
   unsigned char* binary_ptr = PQunescapeBytea((const unsigned char*)escape_ptr, &mxsz);
   PQfreemem(escape_ptr);

   delete [] fBind[npar];
   fBind[npar] = new char[mxsz+1];
   fBind[npar][mxsz] = '\0';
   memcpy(fBind[npar], binary_ptr, mxsz);

   PQfreemem(binary_ptr);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value to large object and immediately insert the large object into DB.

Bool_t TPgSQLStatement::SetLargeObject(Int_t npar, void* mem, Long_t size, Long_t /*maxsize*/)
{
   // All this needs to happen inside a transaction, or it will NOT work.
   PGresult *res=PQexec(fStmt->fConn,"BEGIN");

   CheckErrResult("GetLargeObject", res, kFALSE);
   PQclear(res);

   Int_t lObjID = lo_creat(fStmt->fConn, INV_READ | INV_WRITE);
   if (lObjID<0) {
      Error("SetLargeObject", "Error in SetLargeObject: %s", PQerrorMessage(fStmt->fConn));
      RollBackTransaction("GetLargeObject");
      return kFALSE;
   }

   Int_t lObjFD = lo_open(fStmt->fConn, lObjID, INV_READ | INV_WRITE);
   if (lObjFD<0) {
      Error("SetLargeObject", "Error in SetLargeObject: %s", PQerrorMessage(fStmt->fConn));
      RollBackTransaction("GetLargeObject");
      return kFALSE;
   }

   Int_t writtenBytes = lo_write(fStmt->fConn, lObjFD, (char*)mem, size);

   if (writtenBytes != size) {
      Error("SetLargeObject", "SQL Error on lo_write: %s", PQerrorMessage(fStmt->fConn));
      RollBackTransaction("GetLargeObject");
      return kFALSE;
   }

   if (lo_close(fStmt->fConn, lObjFD) != 0) {
      Error("SetLargeObject", "SQL Error on lo_close: %s", PQerrorMessage(fStmt->fConn));
      RollBackTransaction("GetLargeObject");
      return kFALSE;
   }

   res=PQexec(fStmt->fConn,"COMMIT");
   ExecStatusType stat = PQresultStatus(res);
   if (!pgsql_success(stat)) {
      Error("SetLargeObject", "SQL Error on COMMIT: %s", PQerrorMessage(fStmt->fConn));
      PQclear(res);
      return kFALSE;
   }
   PQclear(res);

   snprintf(fBind[npar],kBindStringSize,"%d",lObjID);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as date.

Bool_t TPgSQLStatement::SetDate(Int_t npar, Int_t year, Int_t month, Int_t day)
{
   TDatime d =TDatime(year,month,day,0,0,0);
   snprintf(fBind[npar],kBindStringSize,"%s",(char*)d.AsSQLString());

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as time.

Bool_t TPgSQLStatement::SetTime(Int_t npar, Int_t hour, Int_t min, Int_t sec)
{
   TDatime d=TDatime(2000,1,1,hour,min,sec);
   snprintf(fBind[npar],kBindStringSize,"%s",(char*)d.AsSQLString());
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as date & time.

Bool_t TPgSQLStatement::SetDatime(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec)
{
   TDatime d=TDatime(year,month,day,hour,min,sec);
   snprintf(fBind[npar],kBindStringSize,"%s",(char*)d.AsSQLString());
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as timestamp.
/// Second fraction is assumed as value in microseconds,
/// i.e. as a fraction with six decimal places.
/// See GetTimestamp() for an example.

Bool_t TPgSQLStatement::SetTimestamp(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec, Int_t frac)
{
   TDatime d(year,month,day,hour,min,sec);
   snprintf(fBind[npar],kBindStringSize,"%s.%06d",(char*)d.AsSQLString(),frac);
   return kTRUE;
}

#else

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.
/// For PgSQL version < 8.2 no statement is supported.

TPgSQLStatement::TPgSQLStatement(PgSQL_Stmt_t*, Bool_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TPgSQLStatement::~TPgSQLStatement()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Close statement.

void TPgSQLStatement::Close(Option_t *)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Process statement.

Bool_t TPgSQLStatement::Process()
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of affected rows after statement is processed.

Int_t TPgSQLStatement::GetNumAffectedRows()
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of statement parameters.

Int_t TPgSQLStatement::GetNumParameters()
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Store result of statement processing to access them
/// via GetInt(), GetDouble() and so on methods.

Bool_t TPgSQLStatement::StoreResult()
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of fields in result set.

Int_t TPgSQLStatement::GetNumFields()
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns field name in result set.

const char* TPgSQLStatement::GetFieldName(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Shift cursor to nect row in result set.

Bool_t TPgSQLStatement::NextResultRow()
{
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Increment iteration counter for statement, where parameter can be set.
/// Statement with parameters of previous iteration
/// automatically will be applied to database.

Bool_t TPgSQLStatement::NextIteration()
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Release all buffers, used by statement.

void TPgSQLStatement::FreeBuffers()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate buffers for statement parameters/ result fields.

void TPgSQLStatement::SetBuffersNumber(Int_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Convert field value to string.

const char* TPgSQLStatement::ConvertToString(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert field to numeric value.

long double TPgSQLStatement::ConvertToNumeric(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if field value is null.

Bool_t TPgSQLStatement::IsNull(Int_t)
{
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as integer.

Int_t TPgSQLStatement::GetInt(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as unsigned integer.

UInt_t TPgSQLStatement::GetUInt(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as long integer.

Long_t TPgSQLStatement::GetLong(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as 64-bit integer.

Long64_t TPgSQLStatement::GetLong64(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as unsigned 64-bit integer.

ULong64_t TPgSQLStatement::GetULong64(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as double.

Double_t TPgSQLStatement::GetDouble(Int_t)
{
   return 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as string.

const char *TPgSQLStatement::GetString(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as binary array.

Bool_t TPgSQLStatement::GetBinary(Int_t, void* &, Long_t&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return large object whose oid is in the given field.

Bool_t TPgSQLStatement::GetLargeObject(Int_t, void* &, Long_t&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as date.

Bool_t TPgSQLStatement::GetDate(Int_t, Int_t&, Int_t&, Int_t&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as time.

Bool_t TPgSQLStatement::GetTime(Int_t, Int_t&, Int_t&, Int_t&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as date & time.

Bool_t TPgSQLStatement::GetDatime(Int_t, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as time stamp.

Bool_t TPgSQLStatement::GetTimestamp(Int_t, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter type to be used as buffer.
/// Used in both setting data to database and retriving data from data base.
/// Initialize proper PGSQL_BIND structure and allocate required buffers.

Bool_t TPgSQLStatement::SetSQLParamType(Int_t, int, bool, int)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set NULL as parameter value.
/// If NULL should be set for statement parameter during first iteration,
/// one should call before proper Set... method to identify type of argument for
/// the future. For instance, if one suppose to have double as type of parameter,
/// code should look like:
///    stmt->SetDouble(2, 0.);
///    stmt->SetNull(2);

Bool_t TPgSQLStatement::SetNull(Int_t npar)
{
   if ((npar >= 0) && (npar < fNumBuffers))
      fBind[npar] = 0;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as integer.

Bool_t TPgSQLStatement::SetInt(Int_t, Int_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as unsigned integer.

Bool_t TPgSQLStatement::SetUInt(Int_t, UInt_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as long integer.

Bool_t TPgSQLStatement::SetLong(Int_t, Long_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as 64-bit integer.

Bool_t TPgSQLStatement::SetLong64(Int_t, Long64_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as unsigned 64-bit integer.

Bool_t TPgSQLStatement::SetULong64(Int_t, ULong64_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as double.

Bool_t TPgSQLStatement::SetDouble(Int_t, Double_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as string.

Bool_t TPgSQLStatement::SetString(Int_t, const char*, Int_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as binary data.

Bool_t TPgSQLStatement::SetBinary(Int_t, void*, Long_t, Long_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value to large object and immediately insert the large object into DB.

Bool_t TPgSQLStatement::SetLargeObject(Int_t, void*, Long_t, Long_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as date.

Bool_t TPgSQLStatement::SetDate(Int_t, Int_t, Int_t, Int_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as time.

Bool_t TPgSQLStatement::SetTime(Int_t, Int_t, Int_t, Int_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as date & time.

Bool_t TPgSQLStatement::SetDatime(Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameter value as timestamp.

Bool_t TPgSQLStatement::SetTimestamp(Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t)
{
   return kFALSE;
}

#endif  //PG_VERSION_NUM
