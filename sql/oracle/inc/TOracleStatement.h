// @(#)root/oracle:$Id$
// Author: Sergey Linev   6/02/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TOracleStatement
#define ROOT_TOracleStatement

#include "TSQLStatement.h"

#include <occi.h>

#ifdef CONST
#undef CONST
#endif

class TOracleStatement : public TSQLStatement {

protected:

   struct TBufferRec {
      void *membuf{nullptr};
      Long_t bufsize{-1};
      std::string namebuf;
   };

   oracle::occi::Environment *fEnv{nullptr};                 // environment
   oracle::occi::Connection  *fConn{nullptr};                // connection to Oracle
   oracle::occi::Statement   *fStmt{nullptr};                // executed statement
   oracle::occi::ResultSet   *fResult{nullptr};              // query result (rows)
   std::vector<oracle::occi::MetaData> *fFieldInfo{nullptr}; // info for each field in the row
   TBufferRec            *fBuffer{nullptr};                  // buffer of values and field names
   Int_t                  fBufferSize{0};                    // size of fBuffer
   Int_t                  fNumIterations{0};                 // size of internal statement buffer
   Int_t                  fIterCounter{0};                   //counts nextiteration calls and process iterations, if required
   Int_t                  fWorkingMode{0};                   // 1 - settingpars, 2 - getting results
   TString                fTimeFmt;                          // format for date to string conversion, default "MM/DD/YYYY, HH24:MI:SS"

   Bool_t      IsParSettMode() const { return fWorkingMode==1; }
   Bool_t      IsResultSet() const { return (fWorkingMode==2) && (fResult!=0); }

   void        SetBufferSize(Int_t size);
   void        CloseBuffer();

public:
   TOracleStatement(oracle::occi::Environment* env,
                    oracle::occi::Connection* conn,
                    oracle::occi::Statement* stmt,
                    Int_t niter, Bool_t errout = kTRUE);
   virtual ~TOracleStatement();

   TOracleStatement(const TOracleStatement &) = delete;
   TOracleStatement& operator=(const TOracleStatement &) = delete;

   void        Close(Option_t * = "") final;

   Int_t       GetBufferLength() const final { return fNumIterations; }
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
   Bool_t      SetDate(Int_t npar, Int_t year, Int_t month, Int_t day) final;
   Bool_t      SetTime(Int_t npar, Int_t hour, Int_t min, Int_t sec) final;
   Bool_t      SetDatime(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec) final;
   using TSQLStatement::SetTimestamp;
   Bool_t      SetTimestamp(Int_t npar, Int_t year, Int_t month, Int_t day, Int_t hour, Int_t min, Int_t sec, Int_t frac = 0) final;
   void        SetTimeFormating(const char *fmt) final { fTimeFmt = fmt; }
   Bool_t      SetVInt(Int_t npar, const std::vector<Int_t> value, const char* schemaName, const char* typeName) final;
   Bool_t      SetVUInt(Int_t npar, const std::vector<UInt_t> value, const char* schemaName, const char* typeName) final;
   Bool_t      SetVLong(Int_t npar, const std::vector<Long_t> value, const char* schemaName, const char* typeName) final;
   Bool_t      SetVLong64(Int_t npar, const std::vector<Long64_t> value, const char* schemaName, const char* typeName) final;
   Bool_t      SetVULong64(Int_t npar, const std::vector<ULong64_t> value, const char* schemaName, const char* typeName) final;
   Bool_t      SetVDouble(Int_t npar, const std::vector<Double_t> value, const char* schemaName, const char* typeName) final;

   Bool_t      NextIteration() final;

   Bool_t      Process() final;
   Int_t       GetNumAffectedRows() final;

   Bool_t      StoreResult() final;
   Int_t       GetNumFields() final;
   const char *GetFieldName(Int_t nfield) final;
   Bool_t      SetMaxFieldSize(Int_t nfield, Long_t maxsize) final;
   Bool_t      NextResultRow() final;

   Bool_t      IsNull(Int_t) final;
   Int_t       GetInt(Int_t npar) final;
   UInt_t      GetUInt(Int_t npar) final;
   Long_t      GetLong(Int_t npar) final;
   Long64_t    GetLong64(Int_t npar) final;
   ULong64_t   GetULong64(Int_t npar) final;
   Double_t    GetDouble(Int_t npar) final;
   const char *GetString(Int_t npar) final;
   Bool_t      GetBinary(Int_t npar, void* &mem, Long_t& size) final;
   Bool_t      GetDate(Int_t npar, Int_t& year, Int_t& month, Int_t& day) final;
   Bool_t      GetTime(Int_t npar, Int_t& hour, Int_t& min, Int_t& sec) final;
   Bool_t      GetDatime(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec) final;
   using TSQLStatement::GetTimestamp;
   Bool_t      GetTimestamp(Int_t npar, Int_t& year, Int_t& month, Int_t& day, Int_t& hour, Int_t& min, Int_t& sec, Int_t& frac) final;
   Bool_t      GetVInt(Int_t npar, std::vector<Int_t> &value) final;
   Bool_t      GetVUInt(Int_t npar, std::vector<UInt_t> &value) final;
   Bool_t      GetVLong(Int_t npar, std::vector<Long_t> &value) final;
   Bool_t      GetVLong64(Int_t npar, std::vector<Long64_t> &value) final;
   Bool_t      GetVULong64(Int_t npar, std::vector<ULong64_t> &value) final;
   Bool_t      GetVDouble(Int_t npar, std::vector<Double_t> &value) final;

   ClassDefOverride(TOracleStatement, 0); // SQL statement class for Oracle
};

#endif
