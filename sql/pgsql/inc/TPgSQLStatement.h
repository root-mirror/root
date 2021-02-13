// @(#)root/mysql:$Id$
// Author: Dennis Box (dbox@fnal.gov)  3/12/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPgSQLStatement
#define ROOT_TPgSQLStatement

#include "TSQLStatement.h"

#include <libpq-fe.h>
#include <pg_config.h> // to get PG_VERSION_NUM

struct PgSQL_Stmt_t {
   PGconn   *fConn;
   PGresult *fRes;
};

class TPgSQLStatement : public TSQLStatement {

private:


   PgSQL_Stmt_t         *fStmt{nullptr};          //! executed statement
   Int_t                 fNumBuffers{0};          //! number of statement parameters
   char                **fBind{nullptr};          //! array of data for input
   char                **fFieldName{nullptr};     //! array of column names
   Int_t                 fWorkingMode{0};         //! 1 - setting parameters, 2 - retrieving results
   Int_t                 fIterationCount{0};      //! number of iteration
   int                  *fParamLengths{nullptr};  //! length of column
   int                  *fParamFormats{nullptr};  //! data type (OID)
   Int_t                 fNumResultRows{0};
   Int_t                 fNumResultCols{0};

   Bool_t      IsSetParsMode() const { return fWorkingMode==1; }
   Bool_t      IsResultSetMode() const { return fWorkingMode==2; }

   Bool_t      SetSQLParamType(Int_t npar, Bool_t isbinary = kFALSE, Int_t param_len = 0, Int_t maxsize = 0);

   long double ConvertToNumeric(Int_t npar);
   const char *ConvertToString(Int_t npar);

   void        FreeBuffers();
   void        SetBuffersNumber(Int_t n);

   void        ConvertTimeToUTC(const TString &PQvalue, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec);

public:
   TPgSQLStatement(PgSQL_Stmt_t* stmt, Bool_t errout = kTRUE);
   virtual ~TPgSQLStatement();

   void        Close(Option_t * = "") final;

   Int_t       GetBufferLength() const final { return 1; }
   Int_t       GetNumParameters() final;

   Bool_t      SetNull(Int_t npar) final;
   Bool_t      SetInt(Int_t npar, Int_t value) final;
   Bool_t      SetUInt(Int_t npar, UInt_t value) final;
   Bool_t      SetLong(Int_t npar, Long_t value) final;
   Bool_t      SetLong64(Int_t npar, Long64_t value) final;
   Bool_t      SetULong64(Int_t npar, ULong64_t value) final;
   Bool_t      SetDouble(Int_t npar, Double_t value) final;
   Bool_t      SetString(Int_t npar, const char* value, Int_t maxsize = 256) final;
   Bool_t      SetBinary(Int_t npar, void* mem, Long_t size, Long_t maxsize = 0x1000) final;
   Bool_t      SetLargeObject(Int_t npar, void* mem, Long_t size, Long_t maxsize = 0x1000) final;
   Bool_t      SetDate(Int_t npar, Int_t year, Int_t month, Int_t day) final;
   Bool_t      SetTime(Int_t npar, Int_t hour, Int_t min, Int_t sec) final;
   Bool_t      SetDatime(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec) final;
   Bool_t      SetTimestamp(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec, Int_t frac = 0) final;
   Bool_t      SetTimestamp(Int_t npar, const TTimeStamp& tm) final;

   Bool_t      NextIteration() final;

   Bool_t      Process() final;
   Int_t       GetNumAffectedRows() final;

   Bool_t      StoreResult() final;
   Int_t       GetNumFields() final;
   const char *GetFieldName(Int_t nfield) final;
   Bool_t      NextResultRow() final;

   Bool_t      IsNull(Int_t npar) final;
   Int_t       GetInt(Int_t npar) final;
   UInt_t      GetUInt(Int_t npar) final;
   Long_t      GetLong(Int_t npar) final;
   Long64_t    GetLong64(Int_t npar) final;
   ULong64_t   GetULong64(Int_t npar) final;
   Double_t    GetDouble(Int_t npar) final;
   const char *GetString(Int_t npar) final;
   Bool_t      GetBinary(Int_t npar, void* &mem, Long_t& size) final;
   Bool_t      GetLargeObject(Int_t npar, void* &mem, Long_t& size) final;
   Bool_t      GetDate(Int_t npar, Int_t& year, Int_t& month, Int_t& day) final;
   Bool_t      GetTime(Int_t npar, Int_t& hour, Int_t& min, Int_t& sec) final;
   Bool_t      GetDatime(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec) final;
   Bool_t      GetTimestamp(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec, Int_t&) final;
   Bool_t      GetTimestamp(Int_t npar, TTimeStamp& tm) final;

   ClassDefOverride(TPgSQLStatement, 0);  // SQL statement class for PgSQL DB
};

#endif
