// @(#)root/oracle:$Id$
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
//  SQL statement class for Oracle                                      //
//                                                                      //
//  See TSQLStatement class documentation for more details.             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TOracleStatement.h"
#include "TOracleServer.h"
#include "TDataType.h"
#include "snprintf.h"
#include <cstdlib>

ClassImp(TOracleStatement);

using namespace std;
using namespace oracle::occi;


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor of TOracleStatement class
/// On creation time specifies buffer length, which should be
/// used in data fetching or data inserting

TOracleStatement::TOracleStatement(Environment* env, Connection* conn, Statement* stmt, Int_t niter, Bool_t errout) :
   TSQLStatement(errout),
   fEnv(env),
   fConn(conn),
   fStmt(stmt),
   fNumIterations(niter),
   fTimeFmt(TOracleServer::GetDatimeFormat())
{
   if (fStmt) {
      fStmt->setPrefetchMemorySize(1000000);
      fStmt->setPrefetchRowCount(niter);
      fStmt->setMaxIterations(niter);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor of TOracleStatement clas

TOracleStatement::~TOracleStatement()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close Oracle statement
/// Removes and destroys all buffers and metainfo

void TOracleStatement::Close(Option_t *)
{

   if (fFieldInfo)
      delete fFieldInfo;

   if (fResult && fStmt)
      fStmt->closeResultSet(fResult);

   if (fConn && fStmt)
      fConn->terminateStatement(fStmt);

   CloseBuffer();

   fConn = nullptr;
   fStmt =  nullptr;
   fResult = nullptr;
   fFieldInfo = nullptr;
   fIterCounter = 0;
}

// Check that statement is ready for use
#define CheckStatement(method, res)                     \
   {                                                    \
      ClearError();                                     \
      if (fStmt==0) {                                   \
         SetError(-1,"Statement is not correctly initialized",method); \
         return res;                                    \
      }                                                 \
   }

// Check that parameter can be set for statement
#define CheckSetPar(method)                             \
   {                                                    \
      CheckStatement(method, kFALSE);                   \
      if (!IsParSettMode()) {                           \
         SetError(-1,"Parameters cannot be set for this statement", method); \
         return kFALSE;                                 \
      }                                                 \
      if (npar<0) {                                     \
         TString errmsg("Invalid parameter number ");   \
         errmsg+= npar;                                 \
         SetError(-1,errmsg.Data(),method);             \
         return kFALSE;                                 \
      }                                                 \
   }

#define CheckGetField(method, defres)                   \
   {                                                    \
      ClearError();                                     \
      if (!IsResultSet()) {                             \
         SetError(-1,"There is no result set for statement", method); \
         return defres;                                 \
      }                                                 \
      if ((npar<0) || (npar>=fBufferSize)) {                                     \
         TString errmsg("Invalid parameter number ");   \
         errmsg+= npar;                                 \
         SetError(-1,errmsg.Data(),method);             \
         return defres;                                 \
      }                                                 \
   }

////////////////////////////////////////////////////////////////////////////////
/// Set buffer size, which is used to keep string values of
/// currently fetched column.

void TOracleStatement::SetBufferSize(Int_t size)
{
    CloseBuffer();
    if (size<=0) return;
    fBufferSize = size;
    fBuffer = new TBufferRec[size];
    for (Int_t n=0;n<fBufferSize;n++) {
       fBuffer[n].membuf = nullptr;
       fBuffer[n].bufsize = -1;
    }
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy buffers, used in data fetching

void TOracleStatement::CloseBuffer()
{
   if (fBuffer) {
      for (Int_t n=0;n<fBufferSize;n++) {
         if (fBuffer[n].membuf)
            free(fBuffer[n].membuf);
      }

      delete[] fBuffer;
   }
   fBuffer = nullptr;
   fBufferSize = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Process SQL statement

Bool_t TOracleStatement::Process()
{
   CheckStatement("Process", kFALSE);

   try {

      if (IsParSettMode()) {
         fStmt->executeUpdate();
         fWorkingMode = 0;
      } else {
         fStmt->execute();
      }

      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "Process");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of affected rows after statement Process() was called
/// Make sense for queries like SELECT, INSERT, UPDATE

Int_t TOracleStatement::GetNumAffectedRows()
{
   CheckStatement("GetNumAffectedRows", -1);

   try {
      return fStmt->getUpdateCount();
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetNumAffectedRows");
   }
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Return number of parameters in statement
/// Not yet implemented for Oracle

Int_t TOracleStatement::GetNumParameters()
{
   CheckStatement("GetNumParameters", -1);

   Info("GetParametersNumber","Not implemented");

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set NULL as value of parameter npar

Bool_t TOracleStatement::SetNull(Int_t npar)
{
   CheckSetPar("SetNull");

   try {
      fStmt->setNull(npar+1, OCCIINT);

      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetNull");
   }

   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Set integer value for parameter npar

Bool_t TOracleStatement::SetInt(Int_t npar, Int_t value)
{
   CheckSetPar("SetInt");

   try {
      fStmt->setInt(npar+1, value);

      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetInt");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set unsigned integer value for parameter npar

Bool_t TOracleStatement::SetUInt(Int_t npar, UInt_t value)
{
   CheckSetPar("SetUInt");

   try {
      fStmt->setUInt(npar+1, value);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetUInt");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set long integer value for parameter npar

Bool_t TOracleStatement::SetLong(Int_t npar, Long_t value)
{
   CheckSetPar("SetLong");

   try {
      fStmt->setNumber(npar+1, Number(value));
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetLong");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set 64-bit integer value for parameter npar

Bool_t TOracleStatement::SetLong64(Int_t npar, Long64_t value)
{
   CheckSetPar("SetLong64");

   try {
      fStmt->setNumber(npar+1, Number((long double)value));
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetLong64");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set unsigned 64-bit integer value for parameter npar

Bool_t TOracleStatement::SetULong64(Int_t npar, ULong64_t value)
{
   CheckSetPar("SetULong64");

   try {
      fStmt->setNumber(npar+1, Number((long double)value));
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetULong64");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set double value for parameter npar

Bool_t TOracleStatement::SetDouble(Int_t npar, Double_t value)
{
   CheckSetPar("SetDouble");

   try {
      fStmt->setDouble(npar+1, value);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetDouble");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set string value for parameter npar

Bool_t TOracleStatement::SetString(Int_t npar, const char* value, Int_t maxsize)
{
   CheckSetPar("SetString");

   try {

   // this is when NextIteration is called first time
      if (fIterCounter==1) {
         fStmt->setDatabaseNCHARParam(npar+1, true);
         fStmt->setMaxParamSize(npar+1, maxsize);
      }

      fStmt->setString(npar+1, value);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetString");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// set parameter value as binary data

Bool_t TOracleStatement::SetBinary(Int_t npar, void* mem, Long_t size, Long_t maxsize)
{
   CheckSetPar("SetBinary");

   try {

      // this is when NextIteration is called first time
      if (fIterCounter==1)
         fStmt->setMaxParamSize(npar+1, maxsize);

      Bytes buf((unsigned char*) mem, size);

      fStmt->setBytes(npar+1, buf);

      return kTRUE;

   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetBinary");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set date value for parameter npar

Bool_t TOracleStatement::SetDate(Int_t npar, Int_t year, Int_t month, Int_t day)
{
   CheckSetPar("SetDate");

   try {
      Date tm = fStmt->getDate(npar+1);
      int o_year;
      unsigned int o_month, o_day, o_hour, o_minute, o_second;
      tm.getDate(o_year, o_month, o_day, o_hour, o_minute, o_second);
      tm.setDate(year, month, day, o_hour, o_minute, o_second);
      fStmt->setDate(npar+1, tm);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetDate");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set time value for parameter npar

Bool_t TOracleStatement::SetTime(Int_t npar, Int_t hour, Int_t min, Int_t sec)
{
   CheckSetPar("SetTime");

   try {
      Date tm = fStmt->getDate(npar+1);
      int o_year;
      unsigned int o_month, o_day, o_hour, o_minute, o_second;
      tm.getDate(o_year, o_month, o_day, o_hour, o_minute, o_second);
      tm.setDate(o_year, o_month, o_day, hour, min, sec);
      fStmt->setDate(npar+1, tm);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetTime");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set date & time value for parameter npar

Bool_t TOracleStatement::SetDatime(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec)
{
   CheckSetPar("SetDatime");

   try {
      Date tm(fEnv, year, month, day, hour, min, sec);
      fStmt->setDate(npar+1, tm);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetDatime");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set date & time value for parameter npar

Bool_t TOracleStatement::SetTimestamp(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec, Int_t frac)
{
   CheckSetPar("SetTimestamp");

   try {
      Timestamp tm(fEnv, year, month, day, hour, min, sec, frac);
      fStmt->setTimestamp(npar+1, tm);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetTimestamp");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set vector of integer values for parameter npar

Bool_t TOracleStatement::SetVInt(Int_t npar, const std::vector<Int_t> value, const char* schemaName, const char* typeName)
{
   CheckSetPar("SetVInt");

   try {
      setVector(fStmt, npar+1, value, schemaName, typeName);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetVInt");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set vector of unsigned integer values for parameter npar

Bool_t TOracleStatement::SetVUInt(Int_t npar, const std::vector<UInt_t> value, const char* schemaName, const char* typeName)
{
   CheckSetPar("SetVUInt");

   try {
      setVector(fStmt, npar+1, value, schemaName, typeName);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetVUInt");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set vector of long integer values for parameter npar

Bool_t TOracleStatement::SetVLong(Int_t npar, const std::vector<Long_t> value, const char* schemaName, const char* typeName)
{
   CheckSetPar("SetVLong");

   try {
      std::vector<Number> nvec;
      for (std::vector<Long_t>::const_iterator it = value.begin();
           it != value.end();
           ++it) {
         nvec.push_back(Number(*it));
      }
      setVector(fStmt, npar+1, nvec, schemaName, typeName);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetVLong");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set vector of 64-bit integer values for parameter npar

Bool_t TOracleStatement::SetVLong64(Int_t npar, const std::vector<Long64_t> value, const char* schemaName, const char* typeName)
{
   CheckSetPar("SetVLong64");

   try {
      std::vector<Number> nvec;
      for (std::vector<Long64_t>::const_iterator it = value.begin();
           it != value.end();
           ++it) {
        nvec.push_back(Number((long double)*it));
      }
      setVector(fStmt, npar+1, nvec, schemaName, typeName);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetVLong64");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set vector of unsigned 64-bit integer values for parameter npar

Bool_t TOracleStatement::SetVULong64(Int_t npar, std::vector<ULong64_t> value, const char* schemaName, const char* typeName)
{
   CheckSetPar("SetVULong64");

   try {
      std::vector<Number> nvec;
      for (std::vector<ULong64_t>::const_iterator it = value.begin();
           it != value.end();
           ++it) {
        nvec.push_back(Number((long double)*it));
      }
      setVector(fStmt, npar+1, nvec, schemaName, typeName);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetVULong64");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set vector of double values for parameter npar

Bool_t TOracleStatement::SetVDouble(Int_t npar, const std::vector<Double_t> value, const char* schemaName, const char* typeName)
{
   CheckSetPar("SetVDouble");

   try {
      setVector(fStmt, npar+1, value, schemaName, typeName);
      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetVDouble");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add next iteration for statement with parameters

Bool_t TOracleStatement::NextIteration()
{
   CheckStatement("NextIteration", kFALSE);

   try {
      fWorkingMode=1;
      // if number of iterations achievs limit, execute it and continue to fill
      if ((fIterCounter % fNumIterations == 0) && (fIterCounter>0)) {
         fStmt->executeUpdate();
      }

      if (fIterCounter % fNumIterations != 0) {
         fStmt->addIteration();
      }

      fIterCounter++;

      return kTRUE;
   } catch (SQLException &oraex)  {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "NextIteration");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Store result of statement processing.
/// Required to access results of SELECT queries

Bool_t TOracleStatement::StoreResult()
{
   CheckStatement("StoreResult", kFALSE);

   try {
      if (fStmt->status() == Statement::RESULT_SET_AVAILABLE) {
         fResult      = fStmt->getResultSet();
         fFieldInfo   = (fResult==0) ? 0 : new std::vector<MetaData>(fResult->getColumnListMetaData());
         Int_t count  = (fFieldInfo==0) ? 0 : fFieldInfo->size();
         SetBufferSize(count);
         if ((fResult!=0) && (count>0)) fWorkingMode = 2;

         return IsResultSet();
      }
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "StoreResult");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Defines maximum size for field which must be used for read or write operation
/// Some Oracle types as LONG (long binary continer) requires this call
/// before any data can be read from database. Call it once before first call to NextResultRow()

Bool_t TOracleStatement::SetMaxFieldSize(Int_t nfield, Long_t maxsize)
{
   CheckStatement("SetMaxFieldSize", kFALSE);

   try {
      if (fResult)
         fResult->setMaxColumnSize(nfield+1, maxsize);
      else
         fStmt->setMaxParamSize(nfield+1, maxsize);
      return kTRUE;
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "SetMaxFieldSize");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns number of fields in result set

Int_t TOracleStatement::GetNumFields()
{
   return IsResultSet() ?  fBufferSize : -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field name in result set

const char* TOracleStatement::GetFieldName(Int_t npar)
{
   CheckGetField("GetFieldName", 0);

   if (!IsResultSet() || (npar<0) || (npar>=fBufferSize)) return nullptr;

   if (fBuffer[npar].namebuf.empty())
      fBuffer[npar].namebuf = (*fFieldInfo)[npar].getString(MetaData::ATTR_NAME);

   return fBuffer[npar].namebuf.empty() ? nullptr : fBuffer[npar].namebuf.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// Move cursor to next row in result set.
/// For Oracle it may lead to additional request to database

Bool_t TOracleStatement::NextResultRow()
{
   ClearError();

   if (fResult==0) {
      SetError(-1,"There is no result set for statement", "NextResultRow");
      return kFALSE;
   }

   if (fResult==0) return kFALSE;

   try {
      for (int n=0;n<fBufferSize;n++) {
        if (fBuffer[n].membuf) {
           free(fBuffer[n].membuf);
           fBuffer[n].membuf = nullptr;
        }
        fBuffer[n].bufsize = -1;
      }
      if (fResult->next() == oracle::occi::ResultSet::END_OF_FETCH) {
         fWorkingMode = 0;
         CloseBuffer();
         return kFALSE;
      }
      return kTRUE;
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "NextResultRow");

      if (oraex.getErrorCode()==32108)
         Info("NextResultRow", "Use TSQLStatement::SetMaxFieldSize() to solve a problem");

   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if fieled value in result set is NULL

Bool_t TOracleStatement::IsNull(Int_t npar)
{
   CheckGetField("IsNull", kFALSE);

   try {
      return fResult->isNull(npar+1);
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "IsNull");
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as integer

Int_t TOracleStatement::GetInt(Int_t npar)
{
   CheckGetField("GetInt", 0);

   Int_t res = 0;

   try {
      if (!fResult->isNull(npar+1))
        res = fResult->getInt(npar+1);
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetInt");
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as unsigned integer

UInt_t TOracleStatement::GetUInt(Int_t npar)
{
   CheckGetField("GetUInt", 0);

   UInt_t res = 0;

   try {
      if (!fResult->isNull(npar+1))
        res = fResult->getUInt(npar+1);
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetUInt");
   }

   return res;
}


////////////////////////////////////////////////////////////////////////////////
/// return field value as long integer

Long_t TOracleStatement::GetLong(Int_t npar)
{
   CheckGetField("GetLong", 0);

   Long_t res = 0;

   try {
      if (!fResult->isNull(npar+1))
        res = (Long_t) fResult->getNumber(npar+1);
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetLong");
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as 64-bit integer

Long64_t TOracleStatement::GetLong64(Int_t npar)
{
   CheckGetField("GetLong64", 0);

   Long64_t res = 0;

   try {
      if (!fResult->isNull(npar+1))
        res = (Long64_t) (long double) fResult->getNumber(npar+1);
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetLong64");
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as unsigned 64-bit integer

ULong64_t TOracleStatement::GetULong64(Int_t npar)
{
   CheckGetField("GetULong64", 0);

   ULong64_t res = 0;

   try {
      if (!fResult->isNull(npar+1))
        res = (ULong64_t) (long double) fResult->getNumber(npar+1);
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetULong64");
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as double

Double_t TOracleStatement::GetDouble(Int_t npar)
{
   CheckGetField("GetDouble", 0.);

   Double_t res = 0;

   try {
      if (!fResult->isNull(npar+1))
        res = fResult->getDouble(npar+1);
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetDouble");
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as string

const char* TOracleStatement::GetString(Int_t npar)
{
   CheckGetField("GetString", 0);

   if (fBuffer[npar].membuf)
      return (const char *) fBuffer[npar].membuf;

   try {
      if (fResult->isNull(npar+1)) return 0;

      int datatype = (*fFieldInfo)[npar].getInt(MetaData::ATTR_DATA_TYPE);

      std::string res;

      switch (datatype) {
        case SQLT_NUM: { // oracle numeric NUMBER
           int prec = (*fFieldInfo)[npar].getInt(MetaData::ATTR_PRECISION);
           int scale = (*fFieldInfo)[npar].getInt(MetaData::ATTR_SCALE);

           if ((scale == 0) || (prec == 0)) {
              res = fResult->getString(npar+1);
           } else {
              double double_val = fResult->getDouble(npar+1);
              char str_number[50];
              snprintf(str_number, sizeof(str_number), TSQLServer::GetFloatFormat(), double_val);
              res = str_number;
           }
           break;
        }
        case SQLT_CHR:  // character string
        case SQLT_VCS:  // variable character string
        case SQLT_AFC: // ansi fixed char
        case SQLT_AVC: // ansi var char
           res = fResult->getString(npar+1);
           break;
        case SQLT_DAT:  // Oracle native DATE type
           res = (fResult->getDate(npar+1)).toText(fTimeFmt.Data());
           break;
        case SQLT_TIMESTAMP:     // TIMESTAMP
        case SQLT_TIMESTAMP_TZ:  // TIMESTAMP WITH TIMEZONE
        case SQLT_TIMESTAMP_LTZ: // TIMESTAMP WITH LOCAL TIMEZONE
           res = (fResult->getTimestamp(npar+1)).toText(fTimeFmt.Data(), 0);
           break;
        default:
           res = fResult->getString(npar+1);
           Info("getString","Type %d may not be supported", datatype);
      }

      int len = res.length();

      if (len > 0) {
          fBuffer[npar].membuf = malloc(len+1);
          fBuffer[npar].bufsize = len+1;
          strncpy((char *) fBuffer[npar].membuf, res.c_str(), len+1);
      }

      return (const char *)fBuffer[npar].membuf;

   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetString");
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return field value as binary array
/// Supports LONG, BLOB, CLOB, BFILE, CFILE types of columns
/// Reads complete content of the column, therefore not suitable for
/// big structures

Bool_t TOracleStatement::GetBinary(Int_t npar, void* &mem, Long_t& size)
{
   mem = 0;
   size = 0;

   CheckGetField("GetBinary", kFALSE);

   if (fBuffer[npar].bufsize >= 0) {
      mem = fBuffer[npar].membuf;
      size = fBuffer[npar].bufsize;
      return kTRUE;
   }

   try {
      if (fResult->isNull(npar+1)) return kTRUE;

      int datatype = (*fFieldInfo)[npar].getInt(MetaData::ATTR_DATA_TYPE);

      switch (datatype) {
         case SQLT_LNG: {
            Bytes parbytes = fResult->getBytes(npar+1);

            size = parbytes.length();

            fBuffer[npar].bufsize = size;

            if (size > 0) {
               mem = malloc(size);

               fBuffer[npar].membuf = mem;

               parbytes.getBytes((unsigned char *) mem, size);
            }

            break;
         }

         case SQLT_BLOB: {
            Blob parblob = fResult->getBlob(npar+1);

            size = parblob.length();

            fBuffer[npar].bufsize = size;

            if (size > 0) {
               mem = malloc(size);

               fBuffer[npar].membuf = mem;

               parblob.read(size, (unsigned char *) mem, size);
            }

            break;
         }

         case SQLT_CLOB: {
            Clob parclob = fResult->getClob(npar+1);

            size = parclob.length();

            fBuffer[npar].bufsize = size;

            if (size > 0) {
               mem = malloc(size);

               fBuffer[npar].membuf = mem;

               parclob.read(size, (unsigned char *) mem, size);
            }

            break;
         }

         case SQLT_BFILEE:
         case SQLT_CFILEE: {

            Bfile parbfile = fResult->getBfile(npar+1);

            size = parbfile.length();

            fBuffer[npar].bufsize = size;

            if (size>0) {
               mem = malloc(size);

               fBuffer[npar].membuf = mem;

               parbfile.read(size, (unsigned char *) mem, size);
            }

            break;
         }

         default:
           Error("GetBinary", "Oracle data type %d not supported", datatype);
           SetError(-1, "Unsupported type for binary convertion", "GetBinary");
           return false;
      }

      return kTRUE;

   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetBinary");
   }

   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// return field value as date

Bool_t TOracleStatement::GetDate(Int_t npar, Int_t& year, Int_t& month, Int_t& day)
{
   Int_t hour, min, sec;

   return GetDatime(npar, year, month, day, hour, min, sec);
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as time

Bool_t TOracleStatement::GetTime(Int_t npar, Int_t& hour, Int_t& min, Int_t& sec)
{
   Int_t year, month, day;

   return GetDatime(npar, year, month, day, hour, min, sec);
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as date & time

Bool_t TOracleStatement::GetDatime(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec)
{
   CheckGetField("GetDatime", kFALSE);

   try {
      if (!fResult->isNull(npar+1)) {
         int datatype = (*fFieldInfo)[npar].getInt(MetaData::ATTR_DATA_TYPE);

         if (datatype!=SQLT_DAT) return kFALSE;

         Date tm = fResult->getDate(npar+1);
         int o_year;
         unsigned int o_month, o_day, o_hour, o_minute, o_second;
         tm.getDate(o_year, o_month, o_day, o_hour, o_minute, o_second);
         year = (Int_t) o_year;
         month = (Int_t) o_month;
         day = (Int_t) o_day;
         hour = (Int_t) o_hour;
         min = (Int_t) o_minute;
         sec = (Int_t) o_second;
         return kTRUE;
      }
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetDatime");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as date & time

Bool_t TOracleStatement::GetTimestamp(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec, Int_t& frac)
{
   CheckGetField("GetTimestamp", kFALSE);

   try {
      if (!fResult->isNull(npar+1)) {
         int datatype = (*fFieldInfo)[npar].getInt(MetaData::ATTR_DATA_TYPE);

         if ((datatype!=SQLT_TIMESTAMP) &&
             (datatype!=SQLT_TIMESTAMP_TZ) &&
             (datatype!=SQLT_TIMESTAMP_LTZ)) return kFALSE;

         Timestamp tm = fResult->getTimestamp(npar+1);
         int o_year;
         unsigned int o_month, o_day, o_hour, o_minute, o_second, o_frac;
         tm.getDate(o_year, o_month, o_day);
         tm.getTime(o_hour, o_minute, o_second, o_frac);
         year = (Int_t) o_year;
         month = (Int_t) o_month;
         day = (Int_t) o_day;
         hour = (Int_t) o_hour;
         min = (Int_t) o_minute;
         sec = (Int_t) o_second;
         frac = (Int_t) o_frac;
         return kTRUE;
      }
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetTimestamp");
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as vector of integers

Bool_t TOracleStatement::GetVInt(Int_t npar, std::vector<Int_t> &value)
{
   CheckGetField("GetVInt", kFALSE);
   try {
      if (!fResult->isNull(npar+1))
         getVector(fResult, npar+1, value);
      return kTRUE;
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetVInt");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as vector of unsigned integers

Bool_t TOracleStatement::GetVUInt(Int_t npar, std::vector<UInt_t> &value)
{
   CheckGetField("GetVUInt", kFALSE);
   try {
      if (!fResult->isNull(npar+1))
         getVector(fResult, npar+1, value);
      return kTRUE;
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetVUInt");
   }
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// return field value as vector of long integers

Bool_t TOracleStatement::GetVLong(Int_t npar, std::vector<Long_t> &value)
{
   CheckGetField("GetVLong", kFALSE);
   try {
      std::vector<Number> res;
      if (!fResult->isNull(npar+1))
         getVector(fResult, npar+1, res);
      for (std::vector<Number>::const_iterator it = res.begin();
           it != res.end();
           ++it ) {
         value.push_back((Long_t)*it);
      }
      return kTRUE;
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetVLong");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as vector of 64-bit integers

Bool_t TOracleStatement::GetVLong64(Int_t npar, std::vector<Long64_t> &value)
{
   CheckGetField("GetVLong64", kFALSE);
   try {
      std::vector<Number> res;
      if (!fResult->isNull(npar+1))
         getVector(fResult, npar+1, res);
      for (std::vector<Number>::const_iterator it = res.begin();
           it != res.end();
           ++it ) {
         value.push_back((Long_t)*it);
      }
      return kTRUE;
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetVLong64");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as vector of unsigned 64-bit integers

Bool_t TOracleStatement::GetVULong64(Int_t npar, std::vector<ULong64_t> &value)
{
   CheckGetField("GetVULong64", kFALSE);
   try {
      std::vector<Number> res;
      if (!fResult->isNull(npar+1))
         getVector(fResult, npar+1, res);
      for (std::vector<Number>::const_iterator it = res.begin();
           it != res.end();
           ++it ) {
        value.push_back((Long_t)(long double)*it);
      }
      return kTRUE;
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetVULong64");
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// return field value as vector of doubles

Bool_t TOracleStatement::GetVDouble(Int_t npar, std::vector<Double_t> &value)
{
   CheckGetField("GetVDouble", kFALSE);
   try {
      if (!fResult->isNull(npar+1))
         getVector(fResult, npar+1, value);
      return kTRUE;
   } catch (SQLException &oraex) {
      SetError(oraex.getErrorCode(), oraex.getMessage().c_str(), "GetVDouble");
   }
   return kFALSE;
}

