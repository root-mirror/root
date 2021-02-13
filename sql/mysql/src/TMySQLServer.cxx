// @(#)root/mysql:$Id$
// Author: Fons Rademakers   15/02/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMySQLServer                                                         //
//                                                                      //
// MySQL server plugin implementing the TSQLServer interface.           //
//                                                                      //
// To open a connection to a server use the static method Connect().    //
// The db argument of Connect() is of the form:                         //
//    mysql://<host>[:<port>][/<database>], e.g.                        //
// mysql://pcroot.cern.ch:3456/test                                     //
//                                                                      //
// As an example of connecting to mysql we assume that the server is    //
// running on the local host and that you have access to a database     //
// named "test" by connecting using an account that has a username and  //
// password of "tuser" and "tpass". You can set up this account         //
// by using the "mysql" program to connect to the server as the MySQL   //
// root user and issuing the following statement:                       //
//                                                                      //
// mysql> GRANT ALL ON test.* TO 'tuser'@'localhost' IDENTIFIED BY 'tpass';
//                                                                      //
// If the test database does not exist, create it with this statement:  //
//                                                                      //
// mysql> CREATE DATABASE test;                                         //
//                                                                      //
// If you want to use a different server host, username, password,      //
// or database name, just substitute the appropriate values.            //
// To connect do:                                                       //
//                                                                      //
// TSQLServer *db = TSQLServer::Connect("mysql://localhost/test", "tuser", "tpass");
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMySQLServer.h"
#include "TMySQLResult.h"
#include "TMySQLStatement.h"
#include "TSQLColumnInfo.h"
#include "TSQLTableInfo.h"
#include "TSQLRow.h"
#include "TUrl.h"
#include "TList.h"
#include "TObjString.h"
#include "TObjArray.h"

ClassImp(TMySQLServer);

////////////////////////////////////////////////////////////////////////////////
/// Open a connection to a MySQL DB server. The db arguments should be
/// of the form "mysql://<host>[:<port>][/<database>]", e.g.:
/// "mysql://pcroot.cern.ch:3456/test". The uid is the username and pw
/// the password that should be used for the connection.
///
/// In addition, several parameters can be specified in url after "?" symbol:
///    timeout=N           n is connect timeout is seconds
///    socket=socketname   socketname should be name of Unix socket, used
///                        for connection
///    multi_statements    tell the server that the client may send multiple
///                        statements in a single string (separated by ;);
///    multi_results       tell the server that the client can handle multiple
///                        result sets from multiple-statement executions or
///                        stored procedures
///    reconnect=0|1       enable or disable automatic reconnection to the server
///                        if the connection is found to have been lost
///    compress            use the compressed client/server protocol
///    cnf_file=filename   Read options from the named option file instead of
///                        from my.cnf
///    cnf_group=groupname Read options from the named group from my.cnf or the
///                        file specified with cnf_file option
/// If several parameters are specified, they should be separated by "&" symbol
/// Example of connection argument:
///    TSQLServer::Connect("mysql://host.domain/test?timeout=10&multi_statements");

TMySQLServer::TMySQLServer(const char *db, const char *uid, const char *pw)
{
   fMySQL = nullptr;
   fInfo = "MySQL";

   TUrl url(db);

   if (!url.IsValid()) {
      TString errmsg("malformed db argument ");
      errmsg += db;
      SetError(-1, errmsg.Data(), "TMySQLServer");
      MakeZombie();
      return;
   }

   if (strncmp(url.GetProtocol(), "mysql", 5)) {
      SetError(-1, "protocol in db argument should be mysql://", "TMySQLServer");
      MakeZombie();
      return;
   }

   const char* dbase = url.GetFile();
   if (dbase)
      if (*dbase=='/') dbase++; //skip leading "/" if appears

   fMySQL = new MYSQL;
   mysql_init(fMySQL);

   ULong_t client_flag = 0;
   TString socket;

   TString optstr = url.GetOptions();
   TObjArray* optarr = optstr.Tokenize("&");
   if (optarr!=0) {
      TIter next(optarr);
      TObject *obj = 0;
      while ((obj = next()) != 0) {
         TString opt = obj->GetName();
         opt.ToLower();
         opt.ReplaceAll(" ","");
         if (opt.Contains("timeout=")) {
            opt.Remove(0, 8);
            Int_t timeout = opt.Atoi();
            if (timeout > 0) {
               UInt_t mysqltimeout = (UInt_t) timeout;
               mysql_options(fMySQL, MYSQL_OPT_CONNECT_TIMEOUT, (const char*) &mysqltimeout);
               if (gDebug) Info("TMySQLServer","Set timeout %d",timeout);
            }
         } else
         if (opt.Contains("read_timeout=")) {
           #if MYSQL_VERSION_ID >= 40101
            opt.Remove(0, 13);
            Int_t timeout = opt.Atoi();
            if (timeout > 0) {
               UInt_t mysqltimeout = (UInt_t) timeout;
               mysql_options(fMySQL, MYSQL_OPT_READ_TIMEOUT, (const char*) &mysqltimeout);
               if (gDebug) Info("TMySQLServer","Set read timeout %d", timeout);
            }
           #else
            Warning("TMySQLServer","MYSQL_OPT_READ_TIMEOUT option not supported by this version of MySql");
           #endif

         } else
         if (opt.Contains("write_timeout=")) {
           #if MYSQL_VERSION_ID >= 40101
            opt.Remove(0, 14);
            Int_t timeout = opt.Atoi();
            if (timeout > 0) {
               UInt_t mysqltimeout = (UInt_t) timeout;
               mysql_options(fMySQL, MYSQL_OPT_WRITE_TIMEOUT, (const char*) &mysqltimeout);
               if (gDebug) Info("TMySQLServer","Set write timeout %d", timeout);
            }
           #else
            Warning("TMySQLServer","MYSQL_OPT_WRITE_TIMEOUT option not supported by this version of MySql");
           #endif
         } else
         if (opt.Contains("reconnect=")) {
           #if MYSQL_VERSION_ID >= 50013
            opt.Remove(0, 10);
            bool reconnect_on = (opt=="1") || (opt=="true");
            mysql_options(fMySQL, MYSQL_OPT_RECONNECT, (const char*) &reconnect_on);
            if (gDebug) Info("TMySQLServer","Set reconnect options %s", (reconnect_on ? "ON" : "OFF"));
           #else
            Warning("TMySQLServer","MYSQL_OPT_RECONNECT option not supported by this version of MySql");
           #endif
         } else
         if (opt.Contains("socket=")) {
            socket = (obj->GetName()+7);
            if (gDebug) Info("TMySQLServer","Use socket %s", socket.Data());
         } else
         if (opt.Contains("multi_statements")) {
           #if MYSQL_VERSION_ID >= 40100
            client_flag = client_flag | CLIENT_MULTI_STATEMENTS;
            if (gDebug) Info("TMySQLServer","Use CLIENT_MULTI_STATEMENTS");
           #else
            Warning("TMySQLServer","CLIENT_MULTI_STATEMENTS not supported by this version of MySql");
           #endif
         } else
         if (opt.Contains("multi_results")) {
           #if MYSQL_VERSION_ID >= 40100
            client_flag = client_flag | CLIENT_MULTI_RESULTS;
            if (gDebug) Info("TMySQLServer","Use CLIENT_MULTI_RESULTS");
           #else
            Warning("TMySQLServer","CLIENT_MULTI_RESULTS not supported by this version of MySql");
           #endif
         } else
         if (opt.Contains("compress")) {
            mysql_options(fMySQL, MYSQL_OPT_COMPRESS, 0);
            if (gDebug) Info("TMySQLServer","Use compressed client/server protocol");
         } else
         if (opt.Contains("cnf_file=")) {
            const char* filename = (obj->GetName()+9);
            mysql_options(fMySQL, MYSQL_READ_DEFAULT_FILE, filename);
            if (gDebug) Info("TMySQLServer","Read mysql options from %s file", filename);
         } else
         if (opt.Contains("cnf_group=")) {
            const char* groupname = (obj->GetName()+10);
            mysql_options(fMySQL, MYSQL_READ_DEFAULT_GROUP, groupname);
            if (gDebug) Info("TMySQLServer","Read mysql options from %s group of my.cnf file", groupname);
         }
      }
      optarr->Delete();
      delete optarr;
   }

   Int_t port = 3306;
   if (url.GetPort()>0) port = url.GetPort();

   if (mysql_real_connect(fMySQL, url.GetHost(), uid, pw, dbase, port,
                         (socket.Length()>0) ? socket.Data() : 0 , client_flag)) {
      fType = "MySQL";
      fHost = url.GetHost();
      fDB   = dbase;
      fPort = port;
   } else {
      SetError(mysql_errno(fMySQL), mysql_error(fMySQL), "TMySQLServer");
      MakeZombie();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Close connection to MySQL DB server.

TMySQLServer::~TMySQLServer()
{
   if (IsConnected())
      Close();
   delete fMySQL;
}

// Reset error and check that server connected
#define CheckConnect(method, res)                       \
   {                                                    \
      ClearError();                                     \
      if (!IsConnected()) {                             \
         SetError(-1,"MySQL server is not connected",method); \
         return res;                                    \
      }                                                 \
   }


// check last mysql error code
#define CheckErrNo(method, force, res)                  \
   {                                                    \
      unsigned int sqlerrno = mysql_errno(fMySQL);         \
      if ((sqlerrno!=0) || force) {                        \
         const char* sqlerrmsg = mysql_error(fMySQL);      \
         if (sqlerrno==0) { sqlerrno = 11111; sqlerrmsg = "MySQL error"; } \
         SetError(sqlerrno, sqlerrmsg, method);               \
         return res;                                    \
      }                                                 \
   }


////////////////////////////////////////////////////////////////////////////////
/// Close connection to MySQL DB server.

void TMySQLServer::Close(Option_t *)
{
   ClearError();

   if (!fMySQL)
      return;

   mysql_close(fMySQL);
   fPort = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute SQL command. Result object must be deleted by the user.
/// Returns a pointer to a TSQLResult object if successful, 0 otherwise.
/// The result object must be deleted by the user.

TSQLResult *TMySQLServer::Query(const char *sql)
{
   CheckConnect("Query", 0);

   if (mysql_query(fMySQL, sql))
      CheckErrNo("Query",kTRUE,0);

   MYSQL_RES *res = mysql_store_result(fMySQL);
   CheckErrNo("Query", kFALSE, 0);

   return new TMySQLResult(res);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute SQL command which does not produce any result sets.
/// Returns kTRUE if successful.

Bool_t TMySQLServer::Exec(const char* sql)
{
   CheckConnect("Exec", kFALSE);

   if (mysql_query(fMySQL, sql))
      CheckErrNo("Exec",kTRUE,kFALSE);

   return !IsError();
}

////////////////////////////////////////////////////////////////////////////////
/// Select a database. Returns 0 if successful, non-zero otherwise.

Int_t TMySQLServer::SelectDataBase(const char *dbname)
{
   CheckConnect("SelectDataBase", -1);

   Int_t res = mysql_select_db(fMySQL, dbname);
   if (res==0) fDB = dbname;
          else CheckErrNo("SelectDataBase", kTRUE, res);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// List all available databases. Wild is for wildcarding "t%" list all
/// databases starting with "t".
/// Returns a pointer to a TSQLResult object if successful, 0 otherwise.
/// The result object must be deleted by the user.

TSQLResult *TMySQLServer::GetDataBases(const char *wild)
{
   CheckConnect("GetDataBases", 0);

   MYSQL_RES *res = mysql_list_dbs(fMySQL, wild);

   CheckErrNo("GetDataBases", kFALSE, 0);

   return new TMySQLResult(res);
}

////////////////////////////////////////////////////////////////////////////////
/// List all tables in the specified database. Wild is for wildcarding
/// "t%" list all tables starting with "t".
/// Returns a pointer to a TSQLResult object if successful, 0 otherwise.
/// The result object must be deleted by the user.

TSQLResult *TMySQLServer::GetTables(const char *dbname, const char *wild)
{
   CheckConnect("GetTables", 0);

   if (SelectDataBase(dbname) != 0) return nullptr;

   MYSQL_RES *res = mysql_list_tables(fMySQL, wild);

   CheckErrNo("GetTables", kFALSE, 0);

   return new TMySQLResult(res);
}


////////////////////////////////////////////////////////////////////////////////
/// Return list of tables with specified wildcard.

TList* TMySQLServer::GetTablesList(const char* wild)
{
   CheckConnect("GetTablesList", 0);

   MYSQL_RES *res = mysql_list_tables(fMySQL, wild);

   CheckErrNo("GetTablesList", kFALSE, 0);

   MYSQL_ROW row = mysql_fetch_row(res);

   TList *lst = nullptr;

   while (row!=0) {
      CheckErrNo("GetTablesList", kFALSE, lst);

      const char *tablename = row[0];

      if (tablename) {
         if (!lst) {
            lst = new TList();
            lst->SetOwner(kTRUE);
         }
         lst->Add(new TObjString(tablename));
      }

      row = mysql_fetch_row(res);
   }

   mysql_free_result(res);

   return lst;
}

////////////////////////////////////////////////////////////////////////////////
/// Produces SQL table info.
/// Object must be deleted by user.

TSQLTableInfo *TMySQLServer::GetTableInfo(const char* tablename)
{
   CheckConnect("GetTableInfo", 0);

   if (!tablename || (*tablename==0)) return nullptr;

   TString sql;
   sql.Form("SELECT * FROM `%s` LIMIT 1", tablename);

   if (mysql_query(fMySQL, sql.Data()) != 0)
      CheckErrNo("GetTableInfo", kTRUE, 0);

   MYSQL_RES *res = mysql_store_result(fMySQL);
   CheckErrNo("GetTableInfo", kFALSE, 0);

   unsigned int numfields = mysql_num_fields(res);

   MYSQL_FIELD* fields = mysql_fetch_fields(res);

   sql.Form("SHOW COLUMNS FROM `%s`", tablename);
   TSQLResult* showres = Query(sql.Data());

   if (!showres) {
      mysql_free_result(res);
      return nullptr;
   }

   TList *lst = nullptr;

   unsigned int nfield = 0;

   TSQLRow* row = nullptr;

   while ((row = showres->Next()) != 0) {
      const char* column_name = row->GetField(0);
      const char* type_name = row->GetField(1);

      if ((nfield>=numfields) ||
          (strcmp(column_name, fields[nfield].name)!=0))
      {
         SetError(-1,"mismatch in column names","GetTableInfo");
         break;
      }

      Int_t sqltype = kSQL_NONE;

      Int_t data_size = -1;    // size in bytes
      Int_t data_length = -1;  // declaration like VARCHAR(n) or NUMERIC(n)
      Int_t data_scale = -1;   // second argument in declaration
      Int_t data_sign = -1; // signed type or not

      if (IS_NUM(fields[nfield].type)) {
         if (fields[nfield].flags & UNSIGNED_FLAG)
            data_sign = 0;
         else
            data_sign = 1;
      }

      Bool_t nullable = (fields[nfield].flags & NOT_NULL_FLAG) == 0;

      data_length = fields[nfield].length;
      if (data_length==0) data_length = -1;

#if MYSQL_VERSION_ID >= 40100

      switch (fields[nfield].type) {
         case MYSQL_TYPE_TINY:
         case MYSQL_TYPE_SHORT:
         case MYSQL_TYPE_LONG:
         case MYSQL_TYPE_INT24:
         case MYSQL_TYPE_LONGLONG:
            sqltype = kSQL_INTEGER;
            break;
         case MYSQL_TYPE_DECIMAL:
            sqltype = kSQL_NUMERIC;
            data_scale = fields[nfield].decimals;
            break;
         case MYSQL_TYPE_FLOAT:
            sqltype = kSQL_FLOAT;
            break;
         case MYSQL_TYPE_DOUBLE:
            sqltype = kSQL_DOUBLE;
            break;
         case MYSQL_TYPE_TIMESTAMP:
            sqltype = kSQL_TIMESTAMP;
            break;
         case MYSQL_TYPE_DATE:
         case MYSQL_TYPE_TIME:
         case MYSQL_TYPE_DATETIME:
         case MYSQL_TYPE_YEAR:
            break;
         case MYSQL_TYPE_STRING:
            if (fields[nfield].charsetnr==63)
               sqltype = kSQL_BINARY;
            else
               sqltype = kSQL_CHAR;
            data_size = data_length;
            break;
         case MYSQL_TYPE_VAR_STRING:
            if (fields[nfield].charsetnr==63)
               sqltype = kSQL_BINARY;
            else
               sqltype = kSQL_VARCHAR;
            data_size = data_length;
            break;
         case MYSQL_TYPE_BLOB:
            if (fields[nfield].charsetnr==63)
               sqltype = kSQL_BINARY;
            else
               sqltype = kSQL_VARCHAR;
            data_size = data_length;
            break;
         case MYSQL_TYPE_SET:
         case MYSQL_TYPE_ENUM:
         case MYSQL_TYPE_GEOMETRY:
         case MYSQL_TYPE_NULL:
            break;
         default:
            if (IS_NUM(fields[nfield].type))
               sqltype = kSQL_NUMERIC;
      }

#endif

      if (!lst)
         lst = new TList;
      lst->Add(new TSQLColumnInfo(column_name,
                                  type_name,
                                  nullable,
                                  sqltype,
                                  data_size,
                                  data_length,
                                  data_scale,
                                  data_sign));

      nfield++;
      delete row;
   }

   mysql_free_result(res);
   delete showres;

   sql.Form("SHOW TABLE STATUS LIKE '%s'", tablename);

   TSQLTableInfo* info = 0;

   TSQLResult* stats = Query(sql.Data());

   if (stats!=0) {
      row = 0;

      while ((row = stats->Next()) != 0) {
         if (strcmp(row->GetField(0), tablename)!=0) {
            delete row;
            continue;
         }
         const char* comments = 0;
         const char* engine = 0;
         const char* create_time = 0;
         const char* update_time = 0;

         for (int n=1;n<stats->GetFieldCount();n++) {
            TString fname = stats->GetFieldName(n);
            fname.ToLower();
            if (fname=="engine") engine = row->GetField(n); else
            if (fname=="comment") comments = row->GetField(n); else
            if (fname=="create_time") create_time = row->GetField(n); else
            if (fname=="update_time") update_time = row->GetField(n);
         }

         info = new TSQLTableInfo(tablename,
                                  lst,
                                  comments,
                                  engine,
                                  create_time,
                                  update_time);

         delete row;
         break;
      }
      delete stats;
   }

   if (info==0)
      info = new TSQLTableInfo(tablename, lst);

   return info;
}

////////////////////////////////////////////////////////////////////////////////
/// List all columns in specified table in the specified database.
/// Wild is for wildcarding "t%" list all columns starting with "t".
/// Returns a pointer to a TSQLResult object if successful, 0 otherwise.
/// The result object must be deleted by the user.

TSQLResult *TMySQLServer::GetColumns(const char *dbname, const char *table,
                                     const char *wild)
{
   CheckConnect("GetColumns", 0);

   if (SelectDataBase(dbname) != 0) return nullptr;

   TString sql;
   if (wild)
      sql.Form("SHOW COLUMNS FROM %s LIKE '%s'", table, wild);
   else
      sql.Form("SHOW COLUMNS FROM %s", table);

   return Query(sql.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Create a database. Returns 0 if successful, non-zero otherwise.

Int_t TMySQLServer::CreateDataBase(const char *dbname)
{
   CheckConnect("CreateDataBase", -1);

   Int_t res = mysql_query(fMySQL, Form("CREATE DATABASE %s",dbname));

   CheckErrNo("CreateDataBase", kFALSE, res);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Drop (i.e. delete) a database. Returns 0 if successful, non-zero
/// otherwise.

Int_t TMySQLServer::DropDataBase(const char *dbname)
{
   CheckConnect("DropDataBase", -1);

   Int_t res = mysql_query(fMySQL, Form("DROP DATABASE %s",dbname));

   CheckErrNo("DropDataBase", kFALSE, res);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Reload permission tables. Returns 0 if successful, non-zero
/// otherwise. User must have reload permissions.

Int_t TMySQLServer::Reload()
{
   CheckConnect("Reload", -1);

   Int_t res = mysql_reload(fMySQL);

   CheckErrNo("Reload", kFALSE, res);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Shutdown the database server. Returns 0 if successful, non-zero
/// otherwise. User must have shutdown permissions.

Int_t TMySQLServer::Shutdown()
{
   CheckConnect("Shutdown", -1);

   Int_t res;

#if MYSQL_VERSION_ID >= 50001 || \
    (MYSQL_VERSION_ID < 50000 && MYSQL_VERSION_ID >= 40103)
   res = mysql_shutdown(fMySQL, SHUTDOWN_DEFAULT);
#else
   res = mysql_shutdown(fMySQL);
#endif

   CheckErrNo("Shutdown", kFALSE, res);

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Return server info in form "MySQL <vesrion>".

const char *TMySQLServer::ServerInfo()
{
   CheckConnect("ServerInfo", 0);

   const char* res = mysql_get_server_info(fMySQL);

   CheckErrNo("ServerInfo", kFALSE, res);

   fInfo = "MySQL ";
   fInfo += res;

   return fInfo.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if TSQLStatement class is supported.
/// Starts from MySQL 4.1.

Bool_t TMySQLServer::HasStatement() const
{
#if MYSQL_VERSION_ID < 40100
   return kFALSE;
#else
   return kTRUE;
#endif
}


////////////////////////////////////////////////////////////////////////////////
/// Produce TMySQLStatement.

TSQLStatement *TMySQLServer::Statement(const char *sql, Int_t)
{
#if MYSQL_VERSION_ID < 40100
   ClearError();
   SetError(-1, "Statement class does not supported by MySQL version < 4.1", "Statement");
   return nullptr;
#else

   CheckConnect("Statement", 0);

   if (!sql || !*sql) {
      SetError(-1, "no query string specified","Statement");
      return nullptr;
   }

   MYSQL_STMT *stmt = mysql_stmt_init(fMySQL);
   if (!stmt)
      CheckErrNo("Statement", kTRUE, 0);

   if (mysql_stmt_prepare(stmt, sql, strlen(sql))) {
      SetError(mysql_errno(fMySQL), mysql_error(fMySQL), "Statement");
      mysql_stmt_close(stmt);
      return nullptr;
   }

   return new TMySQLStatement(stmt, fErrorOut);

#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Start transaction

Bool_t TMySQLServer::StartTransaction()
{
   CheckConnect("StartTransaction", kFALSE);

   return TSQLServer::StartTransaction();
}

////////////////////////////////////////////////////////////////////////////////
/// Commit changes

Bool_t TMySQLServer::Commit()
{
   CheckConnect("Commit", kFALSE);

#if MYSQL_VERSION_ID >= 40100

   if (mysql_commit(fMySQL))
      CheckErrNo("Commit", kTRUE, kFALSE);

   return kTRUE;

#else

   return TSQLServer::Commit();

#endif

}

////////////////////////////////////////////////////////////////////////////////
/// Rollback changes

Bool_t TMySQLServer::Rollback()
{
   CheckConnect("Rollback", kFALSE);

#if MYSQL_VERSION_ID >= 40100

   if (mysql_rollback(fMySQL))
      CheckErrNo("Rollback", kTRUE, kFALSE);

   return kTRUE;

#else

   return TSQLServer::Rollback();

#endif

}

////////////////////////////////////////////////////////////////////////////////
/// Execute Ping to SQL Connection.
/// Since mysql_ping tries to reconnect by itself,
/// a double call to the mysql function is implemented.
/// Returns kTRUE if successful

Bool_t TMySQLServer::PingVerify()
{
   CheckConnect("Ping", kFALSE);

   if (mysql_ping(fMySQL)) {
      if (mysql_ping(fMySQL)) {
         Error("PingVerify", "not able to automatically reconnect a second time");
         CheckErrNo("Ping", kTRUE, kFALSE);
      } else
         Info("PingVerify", "connection was lost, but could automatically reconnect");
   }

   return !IsError();
}

////////////////////////////////////////////////////////////////////////////////
/// Execute Ping to SQL Connection using the mysql_ping function.
/// Returns 0 if successful, non-zero in case an error occured.

Int_t TMySQLServer::Ping()
{
   CheckConnect("PingInt", kFALSE);

   return mysql_ping(fMySQL);
}
