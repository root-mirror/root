// @(#)root/odbc:$Id$
// Author: Sergey Linev   6/02/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TODBCServer.h"

#include "TODBCRow.h"
#include "TODBCResult.h"
#include "TODBCStatement.h"
#include "TSQLColumnInfo.h"
#include "TSQLTableInfo.h"
#include "TUrl.h"
#include "TString.h"
#include "TObjString.h"
#include "TList.h"
#include "strlcpy.h"

#include <iostream>


#include <sqlext.h>


ClassImp(TODBCServer);

////////////////////////////////////////////////////////////////////////////////
/// Open a connection to a ODBC server. The db arguments can be:
/// 1. Form "odbc://[user[:passwd]@]<host>[:<port>][/<database>][?Driver]",
///    e.g.: "odbc://pcroot.cern.ch:3306/test?MySQL".
///    Driver argument specifies ODBC driver, which should be used for
///    connection. By default, MyODBC driver name is used.
///    The uid is the username and pw the password that should be used
///    for the connection.
///    If uid and pw are not specified (==0), user and passwd arguments from
///    URL will be used. Works only with MySQL ODBC, probably with PostrSQL
///    ODBC.
/// 2. Form "odbcd://DRIVER={MyODBC};SERVER=pcroot.cern.ch;DATABASE=test;USER=user;PASSWORD=pass;OPTION=3;PORT=3306;"
///    This is a form, which is accepted by SQLDriverConnect function of ODBC.
///    Here some other arguments can be specified, which are not included
///    in standard URL format.
/// 3. Form "odbcn://MySpecialConfig", where MySpecialConfig is entry,
///    defined in user DSN (user data source). Here uid and pw should be
///    always specified.
///
///    Configuring unixODBC under Linux: http://www.unixodbc.org/odbcinst.html
///    Remarks: for variants 1 & 2 it is enough to create/configure
///    odbcinst.ini file. For variant 3 file odbc.ini should be created.
///    Path to this files can be specified in environmental variables like
///      export ODBCINI=/home/my/unixODBC/etc/odbc.ini
///      export ODBCSYSINI=/home/my/unixODBC/etc
///
///    Configuring MySQL ODBC under Windows.
///    Installing ODBC driver for MySQL is enough to use it under Windows.
///    Afer odbcd:// variant can be used with DRIVER={MySQL ODBC 3.51 Driver};
///    To configure User DSN, go into Start menu -> Settings ->
///    Control panel -> Administrative tools-> Data Sources (ODBC).
///
///    To install Oracle ODBC driver for Windows, one should download
///    and install either complete Oracle client (~500 MB), or so-called
///    Instant Client Basic and Instant Client ODBC (~20 MB together).
///    Some remark about Instant Client:
///       1) Two additional DLLs are required: mfc71.dll & msver71.dll
///          They can be found either in MS VC++ 7.1 Free Toolkit or
///          downloaded from other Internet sites
///       2) ORACLE_HOME environment variable should be specified and point to
///           location, where Instant Client files are extracted
///       3) Run odbc_install.exe from account with administrative rights
///       3) In $ORACLE_HOME/network/admin/ directory appropriate *.ora files
///          like ldap.ora, sqlnet.ora, tnsnames.ora should be installed.
///          Contact your Oracle administrator to get these files.
///    After Oracle ODBC driver is installed, appropriate entry in ODBC drivers
///    list like "Oracle in instantclient10_2" should appiar. Connection
///    string example:
///     "odbcd://DRIVER={Oracle in instantclient10_2};DBQ=db-test;UID=user_name;PWD=user_pass;";

TODBCServer::TODBCServer(const char *db, const char *uid, const char *pw) :
   TSQLServer()
{
   TString connstr;
   Bool_t simpleconnect = kTRUE;

   SQLRETURN retcode;
   SQLHWND  hwnd;

   fPort = 1; // indicate that we are connected

   if ((strncmp(db, "odbc", 4)!=0) || (strlen(db)<8)) {
      SetError(-1, "db argument should be started from odbc...","TODBCServer");
      goto zombie;
   }

   if (strncmp(db, "odbc://", 7)==0) {
      TUrl url(db);
      if (!url.IsValid()) {
         SetError(-1, Form("not valid URL: %s", db), "TODBCServer");
         goto zombie;
      }
      const char* driver = "MyODBC";
      const char* dbase = url.GetFile();
      if (dbase!=0)
         if (*dbase=='/') dbase++; //skip leading "/" if appears

      if (((uid==0) || (*uid==0)) && (strlen(url.GetUser())>0)) {
         uid = url.GetUser();
         pw = url.GetPasswd();
      }

      if (strlen(url.GetOptions())!=0) driver = url.GetOptions();

      connstr.Form("DRIVER={%s};"
                   "SERVER=%s;"
                   "DATABASE=%s;"
                   "USER=%s;"
                   "PASSWORD=%s;"
                   "OPTION=3;",
                    driver, url.GetHost(), dbase, uid, pw);
      if (url.GetPort()>0)
          connstr += Form("PORT=%d;", url.GetPort());

      fHost = url.GetHost();
      fPort = url.GetPort()>0 ? url.GetPort() : 1;
      fDB = dbase;
      simpleconnect = kFALSE;
   } else
   if (strncmp(db, "odbcd://", 8)==0) {
      connstr = db+8;
      simpleconnect = kFALSE;
   } else
   if (strncmp(db, "odbcn://", 8)==0) {
      connstr = db+8;
      simpleconnect = kTRUE;
   } else {
      SetError(-1, "db argument is invalid", "TODBCServer");
      goto zombie;
   }

   retcode = SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &fHenv);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;

   /* Set the ODBC version environment attribute */
   retcode = SQLSetEnvAttr(fHenv, SQL_ATTR_ODBC_VERSION, (void*)SQL_OV_ODBC3, 0);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;

   /* Allocate connection handle */
   retcode = SQLAllocHandle(SQL_HANDLE_DBC, fHenv, &fHdbc);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;

   /* Set login timeout to 5 seconds. */
   retcode = SQLSetConnectAttr(fHdbc, SQL_LOGIN_TIMEOUT, (SQLPOINTER) 5, 0);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;

   char sbuf[2048];

   SQLSMALLINT reslen;
   SQLINTEGER reslen1;

   hwnd = 0;

   if (simpleconnect)
      retcode = SQLConnect(fHdbc, (SQLCHAR*) connstr.Data(), SQL_NTS,
                                  (SQLCHAR*) uid, SQL_NTS,
                                  (SQLCHAR*) pw, SQL_NTS);
   else
      retcode = SQLDriverConnect(fHdbc, hwnd,
                                 (SQLCHAR*) connstr.Data(), SQL_NTS,
                                 (SQLCHAR*) sbuf, sizeof(sbuf), &reslen, SQL_DRIVER_NOPROMPT);

   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;

   fType = "ODBC";

   retcode = SQLGetInfo(fHdbc, SQL_USER_NAME, sbuf, sizeof(sbuf), &reslen);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;
   fUserId = sbuf;

   retcode = SQLGetInfo(fHdbc, SQL_DBMS_NAME, sbuf, sizeof(sbuf), &reslen);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;
   fServerInfo = sbuf;
   fType = sbuf;

   retcode = SQLGetInfo(fHdbc, SQL_DBMS_VER, sbuf, sizeof(sbuf), &reslen);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;
   fServerInfo += " ";
   fServerInfo += sbuf;

   // take current catalog - database name
   retcode = SQLGetConnectAttr(fHdbc, SQL_ATTR_CURRENT_CATALOG, sbuf, sizeof(sbuf), &reslen1);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;
   if (fDB.Length()==0) fDB = sbuf;

   retcode = SQLGetInfo(fHdbc, SQL_SERVER_NAME, sbuf, sizeof(sbuf), &reslen);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;
   if (fHost.Length()==0) fHost = sbuf;

/*

   SQLUINTEGER iinfo;
   retcode = SQLGetInfo(fHdbc, SQL_PARAM_ARRAY_ROW_COUNTS, &iinfo, sizeof(iinfo), 0);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;
   Info("Constr", "SQL_PARAM_ARRAY_ROW_COUNTS = %u", iinfo);

   retcode = SQLGetInfo(fHdbc, SQL_PARAM_ARRAY_SELECTS, &iinfo, sizeof(iinfo), 0);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;
   Info("Constr", "SQL_PARAM_ARRAY_SELECTS = %u", iinfo);

   retcode = SQLGetInfo(fHdbc, SQL_BATCH_ROW_COUNT, &iinfo, sizeof(iinfo), 0);
   if (ExtractErrors(retcode, "TODBCServer")) goto zombie;
   Info("Constr", "SQL_BATCH_ROW_COUNT = %u", iinfo);
*/

   return;

zombie:
   fPort = -1;
   fHost = "";
   MakeZombie();
}

////////////////////////////////////////////////////////////////////////////////
/// Close connection to MySQL DB server.

TODBCServer::~TODBCServer()
{
   if (IsConnected())
      Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Produce TList object with list of available
/// ODBC drivers (isdrivers = kTRUE) or data sources (isdrivers = kFALSE)

TList* TODBCServer::ListData(Bool_t isdrivers)
{
   SQLHENV   henv;
   SQLRETURN retcode;

   retcode = SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &henv);
   if ((retcode!=SQL_SUCCESS) && (retcode!=SQL_SUCCESS_WITH_INFO)) return 0;

   retcode = SQLSetEnvAttr(henv, SQL_ATTR_ODBC_VERSION, (void*)SQL_OV_ODBC3, 0);
   if ((retcode!=SQL_SUCCESS) && (retcode!=SQL_SUCCESS_WITH_INFO)) return 0;

   TList* lst = 0;

   char namebuf[2048], optbuf[2048];
   SQLSMALLINT reslen1, reslen2;

   do {
      strlcpy(namebuf, "",2048);
      strlcpy(optbuf, "",2048);
      if (isdrivers)
         retcode = SQLDrivers(henv, (lst==0 ? SQL_FETCH_FIRST : SQL_FETCH_NEXT),
                     (SQLCHAR*) namebuf, sizeof(namebuf), &reslen1,
                     (SQLCHAR*) optbuf, sizeof(optbuf), &reslen2);
      else
         retcode = SQLDataSources(henv, (lst==0 ? SQL_FETCH_FIRST : SQL_FETCH_NEXT),
                     (SQLCHAR*) namebuf, sizeof(namebuf), &reslen1,
                     (SQLCHAR*) optbuf, sizeof(optbuf), &reslen2);

      if (retcode==SQL_NO_DATA) break;
      if ((retcode==SQL_SUCCESS) || (retcode==SQL_SUCCESS_WITH_INFO)) {
         if (lst==0) {
            lst = new TList;
            lst->SetOwner(kTRUE);
         }
         for (int n=0;n<reslen2-1;n++)
            if (optbuf[n]=='\0') optbuf[n] = ';';

         lst->Add(new TNamed(namebuf, optbuf));
      }
   } while ((retcode==SQL_SUCCESS) || (retcode==SQL_SUCCESS_WITH_INFO));

   SQLFreeHandle(SQL_HANDLE_ENV, henv);

   return lst;

}


////////////////////////////////////////////////////////////////////////////////
/// Produce TList object with list of available ODBC drivers
/// User must delete TList object aftewards
/// Name of driver can be used in connecting to data base in form
/// TSQLServer::Connect("odbcd://DRIVER={<drivername>};DBQ=<dbname>;UID=user;PWD=pass;", 0, 0);

TList* TODBCServer::GetDrivers()
{
   return ListData(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Print list of ODBC drivers in form:
///   <name> : <options list>

void TODBCServer::PrintDrivers()
{
   TList* lst = GetDrivers();
   std::cout << "List of ODBC drivers:" << std::endl;
   TIter iter(lst);
   TNamed* n = 0;
   while ((n = (TNamed*) iter()) != 0)
      std::cout << "  " << n->GetName() << " : " << n->GetTitle() << std::endl;
   delete lst;
}

////////////////////////////////////////////////////////////////////////////////
/// Produce TList object with list of available ODBC data sources
/// User must delete TList object aftewards
/// Name of data source can be used later for connection:
/// TSQLServer::Connect("odbcn://<data_source_name>", "user", "pass");

TList* TODBCServer::GetDataSources()
{
   return ListData(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Print list of ODBC data sources in form:
///   <name> : <options list>

void TODBCServer::PrintDataSources()
{
   TList* lst = GetDataSources();
   std::cout << "List of ODBC data sources:" << std::endl;
   TIter iter(lst);
   TNamed* n = 0;
   while ((n = (TNamed*) iter()) != 0)
      std::cout << "  " << n->GetName() << " : " << n->GetTitle() << std::endl;
   delete lst;
}

////////////////////////////////////////////////////////////////////////////////
/// Extract errors, produced by last ODBC function call

Bool_t TODBCServer::ExtractErrors(SQLRETURN retcode, const char* method)
{
   if ((retcode==SQL_SUCCESS) || (retcode==SQL_SUCCESS_WITH_INFO)) return kFALSE;

   SQLINTEGER i = 0;
   SQLINTEGER native;
   SQLCHAR state[7];
   SQLCHAR text[256];
   SQLSMALLINT len;

   while (SQLGetDiagRec(SQL_HANDLE_ENV, fHenv, ++i, state, &native, text,
                          sizeof(text), &len ) == SQL_SUCCESS)
      //Error(method, "%s:%ld:%ld:%s\n", state, i, native, text);
      SetError(native, (const char*) text, method);

   i = 0;

   while (SQLGetDiagRec(SQL_HANDLE_DBC, fHdbc, ++i, state, &native, text,
                          sizeof(text), &len ) == SQL_SUCCESS)
//      Error(method, "%s:%ld:%ld:%s\n", state, i, native, text);
      SetError(native, (const char*) text, method);

   return kTRUE;
}

// Reset error and check that server connected
#define CheckConnect(method, res)                       \
   {                                                    \
      ClearError();                                     \
      if (!IsConnected()) {                             \
         SetError(-1,"ODBC driver is not connected",method); \
         return res;                                    \
      }                                                 \
   }

////////////////////////////////////////////////////////////////////////////////
/// Close connection to MySQL DB server.

void TODBCServer::Close(Option_t *)
{
   SQLDisconnect(fHdbc);
   SQLFreeHandle(SQL_HANDLE_DBC, fHdbc);
   SQLFreeHandle(SQL_HANDLE_ENV, fHenv);
   fPort = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute SQL command. Result object must be deleted by the user.
/// Returns a pointer to a TSQLResult object if successful, 0 otherwise.
/// The result object must be deleted by the user.

TSQLResult *TODBCServer::Query(const char *sql)
{
   CheckConnect("Query", 0);

   SQLRETURN    retcode;
   SQLHSTMT     hstmt;

   SQLAllocHandle(SQL_HANDLE_STMT, fHdbc, &hstmt);

   retcode = SQLExecDirect(hstmt, (SQLCHAR*) sql, SQL_NTS);
   if (ExtractErrors(retcode, "Query")) {
      SQLFreeHandle(SQL_HANDLE_STMT, hstmt);
      return 0;
   }

   return new TODBCResult(hstmt);
}

////////////////////////////////////////////////////////////////////////////////
/// Executes query which does not produce any results set
/// Return kTRUE if successfull

Bool_t TODBCServer::Exec(const char* sql)
{
   CheckConnect("Exec", 0);

   SQLRETURN    retcode;
   SQLHSTMT     hstmt;

   SQLAllocHandle(SQL_HANDLE_STMT, fHdbc, &hstmt);

   retcode = SQLExecDirect(hstmt, (SQLCHAR*) sql, SQL_NTS);

   Bool_t res = !ExtractErrors(retcode, "Exec");

   SQLFreeHandle(SQL_HANDLE_STMT, hstmt);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Select a database. Returns 0 if successful, non-zero otherwise.
/// Not all RDBMS support selecting of database (catalog) after connecting
/// Normally user should specify database name at time of connection

Int_t TODBCServer::SelectDataBase(const char *db)
{
   CheckConnect("SelectDataBase", -1);

   SQLRETURN retcode = SQLSetConnectAttr(fHdbc, SQL_ATTR_CURRENT_CATALOG, (SQLCHAR*) db, SQL_NTS);
   if (ExtractErrors(retcode, "SelectDataBase")) return -1;

   fDB = db;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// List all available databases. Wild is for wildcarding "t%" list all
/// databases starting with "t".
/// Returns a pointer to a TSQLResult object if successful, 0 otherwise.
/// The result object must be deleted by the user.

TSQLResult *TODBCServer::GetDataBases(const char *)
{
   CheckConnect("GetDataBases", 0);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// List all tables in the specified database. Wild is for wildcarding
/// "t%" list all tables starting with "t".
/// Returns a pointer to a TSQLResult object if successful, 0 otherwise.
/// The result object must be deleted by the user.

TSQLResult *TODBCServer::GetTables(const char*, const char* wild)
{
   CheckConnect("GetTables", 0);

   SQLRETURN    retcode;
   SQLHSTMT     hstmt;

   SQLAllocHandle(SQL_HANDLE_STMT, fHdbc, &hstmt);

   SQLCHAR* schemaName = 0;
   SQLSMALLINT schemaNameLength = 0;

/*
   TString schemabuf;
   // schema is used by Oracle to specify to which user belong table
   // therefore, to see correct tables, schema name is set to user name
   if ((fUserId.Length()>0) && (fServerInfo.Contains("Oracle"))) {
      schemabuf = fUserId;
      schemabuf.ToUpper();
      schemaName = (SQLCHAR*) schemabuf.Data();
      schemaNameLength = schemabuf.Length();
   }
*/

   SQLCHAR* tableName = 0;
   SQLSMALLINT tableNameLength = 0;

   if ((wild!=0) && (strlen(wild)!=0)) {
      tableName = (SQLCHAR*) wild;
      tableNameLength = strlen(wild);
      SQLSetStmtAttr(hstmt, SQL_ATTR_METADATA_ID, (SQLPOINTER) SQL_FALSE, 0);
   }

   retcode = SQLTables(hstmt, NULL, 0, schemaName, schemaNameLength, tableName, tableNameLength, (SQLCHAR*) "TABLE", 5);
   if (ExtractErrors(retcode, "GetTables")) {
      SQLFreeHandle(SQL_HANDLE_STMT, hstmt);
      return 0;
   }

   return new TODBCResult(hstmt);
}

////////////////////////////////////////////////////////////////////////////////
/// Return list of tables in database
/// See TSQLServer::GetTablesList() for details.

TList* TODBCServer::GetTablesList(const char* wild)
{
   CheckConnect("GetTablesList", 0);

   TSQLResult* res = GetTables(0, wild);
   if (res==0) return 0;

   TList* lst = 0;

   TSQLRow* row = 0;

   while ((row = res->Next())!=0) {
      const char* tablename = row->GetField(2);
      if (tablename!=0) {
//         Info("List","%s %s %s %s %s", tablename, row->GetField(0), row->GetField(1), row->GetField(3), row->GetField(4));
         if (lst==0) {
            lst = new TList;
            lst->SetOwner(kTRUE);
         }
         lst->Add(new TObjString(tablename));
      }
      delete row;
   }

   delete res;

   return lst;
}


////////////////////////////////////////////////////////////////////////////////
/// Produces SQL table info
/// Object must be deleted by user

TSQLTableInfo* TODBCServer::GetTableInfo(const char* tablename)
{
   CheckConnect("GetTableInfo", 0);

   #define STR_LEN 128+1
   #define REM_LEN 254+1

   /* Declare buffers for result set data */

   SQLCHAR       szCatalog[STR_LEN], szSchema[STR_LEN];
   SQLCHAR       szTableName[STR_LEN], szColumnName[STR_LEN];
   SQLCHAR       szTypeName[STR_LEN], szRemarks[REM_LEN];
   SQLCHAR       szColumnDefault[STR_LEN], szIsNullable[STR_LEN];
   SQLLEN        columnSize, bufferLength, charOctetLength, ordinalPosition;
   SQLSMALLINT   dataType, decimalDigits, numPrecRadix, nullable;
   SQLSMALLINT   sqlDataType, datetimeSubtypeCode;
   SQLRETURN     retcode;
   SQLHSTMT      hstmt;

   /* Declare buffers for bytes available to return */

   SQLLEN cbCatalog, cbSchema, cbTableName, cbColumnName;
   SQLLEN cbDataType, cbTypeName, cbColumnSize, cbBufferLength;
   SQLLEN cbDecimalDigits, cbNumPrecRadix, cbNullable, cbRemarks;
   SQLLEN cbColumnDefault, cbSQLDataType, cbDatetimeSubtypeCode, cbCharOctetLength;
   SQLLEN cbOrdinalPosition, cbIsNullable;


   SQLAllocHandle(SQL_HANDLE_STMT, fHdbc, &hstmt);

   retcode = SQLColumns(hstmt, NULL, 0, NULL, 0, (SQLCHAR*) tablename, SQL_NTS, NULL, 0);
   if (ExtractErrors(retcode, "GetTableInfo")) {
      SQLFreeHandle(SQL_HANDLE_STMT, hstmt);
      return 0;
   }

   TList* lst = 0;

   /* Bind columns in result set to buffers */

   SQLBindCol(hstmt, 1, SQL_C_CHAR, szCatalog, STR_LEN,&cbCatalog);
   SQLBindCol(hstmt, 2, SQL_C_CHAR, szSchema, STR_LEN, &cbSchema);
   SQLBindCol(hstmt, 3, SQL_C_CHAR, szTableName, STR_LEN,&cbTableName);
   SQLBindCol(hstmt, 4, SQL_C_CHAR, szColumnName, STR_LEN, &cbColumnName);
   SQLBindCol(hstmt, 5, SQL_C_SSHORT, &dataType, 0, &cbDataType);
   SQLBindCol(hstmt, 6, SQL_C_CHAR, szTypeName, STR_LEN, &cbTypeName);
   SQLBindCol(hstmt, 7, SQL_C_SLONG, &columnSize, 0, &cbColumnSize);
   SQLBindCol(hstmt, 8, SQL_C_SLONG, &bufferLength, 0, &cbBufferLength);
   SQLBindCol(hstmt, 9, SQL_C_SSHORT, &decimalDigits, 0, &cbDecimalDigits);
   SQLBindCol(hstmt, 10, SQL_C_SSHORT, &numPrecRadix, 0, &cbNumPrecRadix);
   SQLBindCol(hstmt, 11, SQL_C_SSHORT, &nullable, 0, &cbNullable);
   SQLBindCol(hstmt, 12, SQL_C_CHAR, szRemarks, REM_LEN, &cbRemarks);
   SQLBindCol(hstmt, 13, SQL_C_CHAR, szColumnDefault, STR_LEN, &cbColumnDefault);
   SQLBindCol(hstmt, 14, SQL_C_SSHORT, &sqlDataType, 0, &cbSQLDataType);
   SQLBindCol(hstmt, 15, SQL_C_SSHORT, &datetimeSubtypeCode, 0, &cbDatetimeSubtypeCode);
   SQLBindCol(hstmt, 16, SQL_C_SLONG, &charOctetLength, 0, &cbCharOctetLength);
   SQLBindCol(hstmt, 17, SQL_C_SLONG, &ordinalPosition, 0, &cbOrdinalPosition);
   SQLBindCol(hstmt, 18, SQL_C_CHAR, szIsNullable, STR_LEN, &cbIsNullable);

   retcode = SQLFetch(hstmt);

   while ((retcode==SQL_SUCCESS) || (retcode==SQL_SUCCESS_WITH_INFO)) {

      Int_t sqltype = kSQL_NONE;

      Int_t data_size = -1;    // size in bytes
      Int_t data_length = -1;  // declaration like VARCHAR(n) or NUMERIC(n)
      Int_t data_scale = -1;   // second argument in declaration
      Int_t data_sign = -1; // no info about sign

      switch (dataType) {
         case SQL_CHAR:
            sqltype = kSQL_CHAR;
            data_size = columnSize;
            data_length = charOctetLength;
            break;
         case SQL_VARCHAR:
         case SQL_LONGVARCHAR:
            sqltype = kSQL_VARCHAR;
            data_size = columnSize;
            data_length = charOctetLength;
            break;
         case SQL_DECIMAL:
         case SQL_NUMERIC:
            sqltype = kSQL_NUMERIC;
            data_size = columnSize; // size of column in database
            data_length = columnSize;
            data_scale = decimalDigits;
            break;
         case SQL_INTEGER:
         case SQL_TINYINT:
         case SQL_BIGINT:
            sqltype = kSQL_INTEGER;
            data_size = columnSize;
            break;
         case SQL_REAL:
         case SQL_FLOAT:
            sqltype = kSQL_FLOAT;
            data_size = columnSize;
            data_sign = 1;
            break;
         case SQL_DOUBLE:
            sqltype = kSQL_DOUBLE;
            data_size = columnSize;
            data_sign = 1;
            break;
         case SQL_BINARY:
         case SQL_VARBINARY:
         case SQL_LONGVARBINARY:
            sqltype = kSQL_BINARY;
            data_size = columnSize;
            break;
         case SQL_TYPE_TIMESTAMP:
            sqltype = kSQL_TIMESTAMP;
            data_size = columnSize;
            break;
      }

      if (lst==0) lst = new TList;

      lst->Add(new TSQLColumnInfo((const char*) szColumnName,
                                  (const char*) szTypeName,
                                  nullable!=0,
                                  sqltype,
                                  data_size,
                                  data_length,
                                  data_scale,
                                  data_sign));

      retcode = SQLFetch(hstmt);
   }

   SQLFreeHandle(SQL_HANDLE_STMT, hstmt);

   return new TSQLTableInfo(tablename, lst);
}

////////////////////////////////////////////////////////////////////////////////
/// List all columns in specified table in the specified database.
/// Wild is for wildcarding "t%" list all columns starting with "t".
/// Returns a pointer to a TSQLResult object if successful, 0 otherwise.
/// The result object must be deleted by the user.

TSQLResult *TODBCServer::GetColumns(const char*, const char *table, const char*)
{
   CheckConnect("GetColumns", 0);

   SQLRETURN    retcode;
   SQLHSTMT     hstmt;

   SQLAllocHandle(SQL_HANDLE_STMT, fHdbc, &hstmt);

   retcode = SQLColumns(hstmt, NULL, 0, NULL, 0, (SQLCHAR*) table, SQL_NTS, NULL, 0);
   if (ExtractErrors(retcode, "GetColumns")) {
      SQLFreeHandle(SQL_HANDLE_STMT, hstmt);
      return 0;
   }

   return new TODBCResult(hstmt);
}

////////////////////////////////////////////////////////////////////////////////
/// returns maximum allowed length of identifier (table name, column name, index name)

Int_t TODBCServer::GetMaxIdentifierLength()
{
   CheckConnect("GetMaxIdentifierLength", 20);

   SQLUINTEGER info = 0;
   SQLRETURN retcode;

   retcode = SQLGetInfo(fHdbc, SQL_MAX_IDENTIFIER_LEN, (SQLPOINTER)&info, sizeof(info), NULL);

   if (ExtractErrors(retcode, "GetMaxIdentifierLength")) return 20;

   return info;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a database. Returns 0 if successful, non-zero otherwise.

Int_t TODBCServer::CreateDataBase(const char*)
{
   CheckConnect("CreateDataBase", -1);

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Drop (i.e. delete) a database. Returns 0 if successful, non-zero
/// otherwise.

Int_t TODBCServer::DropDataBase(const char*)
{
   CheckConnect("DropDataBase", -1);

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Reload permission tables. Returns 0 if successful, non-zero
/// otherwise. User must have reload permissions.

Int_t TODBCServer::Reload()
{
   CheckConnect("Reload", -1);

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Shutdown the database server. Returns 0 if successful, non-zero
/// otherwise. User must have shutdown permissions.

Int_t TODBCServer::Shutdown()
{
   CheckConnect("Shutdown", -1);

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return server info.

const char *TODBCServer::ServerInfo()
{
   CheckConnect("ServerInfo", 0);

   return fServerInfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates ODBC statement for provided query.
/// See TSQLStatement class for more details.

TSQLStatement *TODBCServer::Statement(const char *sql, Int_t bufsize)
{
   CheckConnect("Statement", 0);

   if (!sql || !*sql) {
      SetError(-1, "no query string specified", "Statement");
      return 0;
   }

//   SQLUINTEGER info = 0;
//   SQLGetInfo(fHdbc, SQL_PARAM_ARRAY_ROW_COUNTS, (SQLPOINTER)&info, sizeof(info), NULL);
//   if (info==SQL_PARC_BATCH) Info("Statement","info==SQL_PARC_BATCH"); else
//   if (info==SQL_PARC_NO_BATCH) Info("Statement","info==SQL_PARC_NO_BATCH"); else
//    Info("Statement","info==%u", info);


   SQLRETURN    retcode;
   SQLHSTMT     hstmt;

   retcode = SQLAllocHandle(SQL_HANDLE_STMT, fHdbc, &hstmt);
   if (ExtractErrors(retcode, "Statement")) return 0;

   retcode = SQLPrepare(hstmt, (SQLCHAR*) sql, SQL_NTS);
   if (ExtractErrors(retcode, "Statement")) {
      SQLFreeHandle(SQL_HANDLE_STMT, hstmt);
      return 0;
   }

   return new TODBCStatement(hstmt, bufsize, fErrorOut);
}

////////////////////////////////////////////////////////////////////////////////
/// Starts transaction.
/// Check for transaction support.
/// Switch off autocommitment mode.

Bool_t TODBCServer::StartTransaction()
{
   CheckConnect("StartTransaction", kFALSE);

   SQLUINTEGER info = 0;
   SQLRETURN retcode;

   retcode = SQLGetInfo(fHdbc, SQL_TXN_CAPABLE, (SQLPOINTER)&info, sizeof(info), NULL);
   if (ExtractErrors(retcode, "StartTransaction")) return kFALSE;

   if (info==0) {
      SetError(-1,"Transactions not supported","StartTransaction");
      return kFALSE;
   }

   if (!Commit()) return kFALSE;

   retcode = SQLSetConnectAttr(fHdbc, SQL_ATTR_AUTOCOMMIT, (SQLPOINTER) SQL_AUTOCOMMIT_OFF, 0);
   if (ExtractErrors(retcode, "StartTransaction")) return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Complete current transaction (commit = kTRUE) or rollback
/// Switches on autocommit mode of ODBC driver

Bool_t TODBCServer::EndTransaction(Bool_t commit)
{
   const char* method = commit ? "Commit" : "Rollback";

   CheckConnect(method, kFALSE);

   SQLRETURN retcode = SQLEndTran(SQL_HANDLE_DBC, fHdbc, commit ? SQL_COMMIT : SQL_ROLLBACK);
   if (ExtractErrors(retcode, method)) return kFALSE;

   retcode = SQLSetConnectAttr(fHdbc, SQL_ATTR_AUTOCOMMIT, (SQLPOINTER) SQL_AUTOCOMMIT_ON, 0);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Commit transaction

Bool_t TODBCServer::Commit()
{
   return EndTransaction(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Rollback transaction

Bool_t TODBCServer::Rollback()
{
   return EndTransaction(kFALSE);
}
