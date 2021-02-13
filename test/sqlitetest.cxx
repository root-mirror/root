#include <iostream>
#include <iomanip>
#include <unistd.h>

#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "TSQLStatement.h"
#include "TSQLTableInfo.h"
#include "TTimeStamp.h"
#include "TList.h"

int main() {

   // Create a new DB called testdb:
   TSQLServer *serv=TSQLServer::Connect("sqlite://testdb.sqlite", "", "");
   if (serv == NULL) {
      std::cerr << "Connection failed!" << std::endl;
      _exit(1);
   }

   // First, some debug-checks:
   std::cout << "DB: " << serv->GetDB() << std::endl;
   std::cout << "DBMS: " << serv->GetDBMS() << std::endl;
   std::cout << "HOST: " << serv->GetHost() << std::endl;
   std::cout << "PORT: " << serv->GetPort() << std::endl;
   std::cout << "Info: " << serv->ServerInfo() << std::endl;

   // Create table:
   if ((serv!=0) && serv->IsConnected()) {
      // create statement instance
      TSQLStatement* stmt = serv->Statement("CREATE TABLE TESTTABLE (ID1 FOO, ID2 FOO, ID3 FOO, ID4 FOO, ID5 FOO, ID6 FOO, ID7 FOO, ID8 FOO, ID9 FOO, ID10 FOO)");
      // process statement
      stmt->Process();
      // destroy object
      delete stmt;
   }

   serv->StartTransaction();//Exec("BEGIN TRANSACTION;");

   // Fill with data:
   TSQLStatement* stmt = serv->Statement("INSERT INTO TESTTABLE (ID1, ID2, ID3, ID4, ID5, ID6, ID7, ID8, ID9, ID10) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", 100);
   std::cout << "statement pars: " << stmt->GetNumParameters() << std::endl;


   for (int n=0;n<10;n++) {
      if (stmt->NextIteration()) {
         stmt->SetInt(0, n);
         stmt->SetUInt(1, n);
         stmt->SetLong(2, n);
         stmt->SetULong64(3, n);
         stmt->SetDate(4, 2013,1,n);
         stmt->SetTime(5, 1,1,n);
         stmt->SetDatime(6, 2013,1,n,1,1,n);
         stmt->SetTimestamp(7, 2013,1,1,1,1,n,102+n);
         TString foo;
         foo.Form("testbinary %d", n);
         void *binary=const_cast<char*>(foo.Data());
         stmt->SetBinary(8, binary, foo.Length());
         stmt->SetString(9, Form("%d", n), 200);
      }
   }

   stmt->Process();
   delete stmt;

   serv->Commit();

   delete serv;

   std::cout << "Testing tool created DB with sample data. Now reopening and selecting all." << std::endl;

   serv=TSQLServer::Connect("sqlite://testdb.sqlite", "", "");
   if (serv == NULL) {
      std::cerr << "Connection failed!" << std::endl;
      _exit(1);
   }

   std::cout << "Select table names:" << std::endl;
   TSQLResult *res=serv->GetTables(NULL, NULL);

   int fields=res->GetFieldCount();
   std::cout << "Table list has field count: " << fields << std::endl;

   int rowcount=res->GetRowCount();
   std::cout << "Table list has row count: " << rowcount << std::endl;

   TSQLRow *row=NULL;
   while ((row=res->Next()) != NULL) {
      for (int i=0; i<fields; i++) {
         std::cout << row->GetField(i);
         TSQLResult *res2=serv->GetColumns(NULL, row->GetField(i), NULL);
         if (res2 == NULL) continue;
         std::cout << "|";
         int fields2=res2->GetFieldCount();
         std::cout << "Cols: " << fields2 << " ";
         TSQLRow *row2=NULL;
         while ((row2=res2->Next()) != NULL) {
            std::cout << "(";
            for (int ii=0; ii<fields2; ii++) {
               if (row2->GetField(ii) == NULL) continue;
               std::cout << row2->GetField(ii) << "|";
            }
            std::cout << ")";
            delete row2;
         }
      }
      std::cout << std::endl;
      delete row;
   }
   delete res;
   std::cout << std::endl;

   std::cout << "Alternate way using GetTablesList:" << std::endl;
   TList *tables = serv->GetTablesList();
   tables->Print();
   delete tables;

   std::cout << "Completed listing tables. Now selecting * from testtable, first using Query() and string output:" << std::endl;
   row=NULL;
   res=serv->Query("SELECT * from TESTTABLE;");
   fields=res->GetFieldCount();
   for (int i=0; i<fields; i++) {
      std::cout << "|" << std::setw(19) << res->GetFieldName(i) << "|";
   }
   std::cout << std::endl;
   while ((row=res->Next()) != NULL) {
      for (int i=0; i<fields; i++) {
         std::cout << "|" << std::setw(19) << row->GetField(i) << "|";
      }
      std::cout << std::endl;
      delete row;
   }
   delete res;
   std::cout << std::endl;

   std::cout << "Now using TSQLStatement-methods with appropriate types:" << std::endl;

   stmt = serv->Statement("SELECT * FROM TESTTABLE;", 100);
   // process statement
   if (stmt->Process()) {
      std::cout << "iteration..." << std::endl;
      // store result of statement in buffer
      stmt->StoreResult();

      // display info about selected field
      std::cout << "NumFields = " << stmt->GetNumFields() << std::endl;
      for (int n=0;n<stmt->GetNumFields();n++) {
         std::cout << "|" << std::setw(19) << stmt->GetFieldName(n) << "|";
      }
      std::cout << std::endl;

      // extract rows one after another
      while (stmt->NextResultRow()) {
         Int_t id1 = stmt->GetInt(0);
         std::cout << "|" << std::setw(19) << id1 << "|";

         UInt_t id2 = stmt->GetUInt(1);
         std::cout << "|" << std::setw(19) << id2 << "|";

         Long_t id3 = stmt->GetLong(2);
         std::cout << "|" << std::setw(19) << id3 << "|";

         ULong64_t id4 = stmt->GetULong64(3);
         std::cout << "|" << std::setw(19) << id4 << "|";

         Int_t year=0, month=0, day=0, hour=0, minute=0, second=0, frac=0;

         stmt->GetDate(4, year, month, day);
         TString id5;
         id5.Form("%04d-%02d-%02d %02d:%02d:%02d", year, month, day, 0, 0, 0);
         std::cout << "|" << std::setw(19) << id5 << "|";

         stmt->GetTime(5, hour, minute, second);
         TString id6;
         id6.Form("%04d-%02d-%02d %02d:%02d:%02d", 2000, 01, 01, hour, minute, second);
         std::cout << "|" << std::setw(19) << id6 << "|";

         stmt->GetDatime(6, year, month, day, hour, minute, second);
         TString id7;
         id7.Form("%04d-%02d-%02d %02d:%02d:%02d", year, month, day, hour, minute, second);
         std::cout << "|" << std::setw(19) << id7 << "|";

         stmt->GetTimestamp(7, year, month, day, hour, minute, second, frac);
         TTimeStamp ts(year,month,day,hour,minute,second);
         TString id8;
         // Frac is in milliseconds for SQLite, thus use %3d here.
         id8.Form("%04d-%02d-%02d %02d:%02d:%02d.%3d", year, month, day, hour, minute, second, frac);
         std::cout << "|" << std::setw(23) << id8 << "|";

         void* id9 = new unsigned char[1];
         Long_t binSize=1;
         stmt->GetBinary(8, id9, binSize);
         char* id9str = new char[binSize+1];
         memcpy(id9str, id9, binSize);
         id9str[binSize]='\0';
         std::cout << "|" << std::setw(19) << id9str << "|";
         delete [] (unsigned char*) id9;
         delete [] id9str;

         TString id10 = stmt->GetString(9);
         std::cout << "|" << std::setw(19) << id10 << "|";

         std::cout << std::endl;
      }
   }
   delete stmt;

   // Test tableinfo:
   std::cout << "Tableinfo:" << std::endl;
   TSQLTableInfo *ti=serv->GetTableInfo("TESTTABLE");
   ti->Print();

   delete serv;

   return 0;
}
