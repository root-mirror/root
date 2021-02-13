// @(#)root/tree:$Id$
// Author: Philippe Canal and al. 08/2004

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TTreeSQL
\ingroup tree

Implement TTree for a SQL backend
*/

#include <vector>
#include <map>
#include <cstdlib>

#include "TString.h"
#include "TError.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TList.h"

#include "TSQLRow.h"
#include "TSQLResult.h"
#include "TSQLServer.h"
#include "TSQLTableInfo.h"
#include "TSQLColumnInfo.h"

#include "TTreeSQL.h"
#include "TBasketSQL.h"

ClassImp(TTreeSQL);

////////////////////////////////////////////////////////////////////////////////
/// Constructor with an explicit TSQLServer

TTreeSQL::TTreeSQL(TSQLServer *server, TString DB, const TString& table) :
   TTree(table.Data(), "Database read from table: " + table, 0), fDB(DB),
   fTable(table.Data()),
   fResult(0), fRow(0),
   fServer(server),
   fBranchChecked(kFALSE),
   fTableInfo(0)
{
   fCurrentEntry = -1;
   fQuery = TString("Select * from " + fTable);
   fEntries = 0;

   if (fServer==0) {
      Error("TTreeSQL","No TSQLServer specified");
      return;
   }
   if (CheckTable(fTable.Data())) {
      Init();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented yet

TBranch* TTreeSQL::BranchImp(const char *, const char *,
                             TClass *, void *, Int_t ,
                             Int_t )
{
   Fatal("BranchImp","Not implemented yet");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented yet

TBranch* TTreeSQL::BranchImp(const char *, TClass *,
                             void *, Int_t , Int_t )
{
   Fatal("BranchImp","Not implemented yet");
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
/// Not implemented yet

Int_t TTreeSQL::Branch(TCollection *, Int_t,
                       Int_t, const char *)
{
   Fatal("Branch","Not implemented yet");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented yet

Int_t TTreeSQL::Branch(TList *, Int_t, Int_t)
{
   Fatal("Branch","Not implemented yet");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented yet

Int_t TTreeSQL::Branch(const char *, Int_t ,
                       Int_t)
{
   Fatal("Branch","Not implemented yet");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented yet

TBranch* TTreeSQL::Bronch(const char *, const char *, void *,
                          Int_t, Int_t)
{
   Fatal("Bronch","Not implemented yet");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented yet

TBranch* TTreeSQL::BranchOld(const char *, const char *,
                             void *, Int_t, Int_t)
{
   Fatal("BranchOld","Not implemented yet");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented yet

TBranch *TTreeSQL::Branch(const char *, const char *, void *,
                          Int_t, Int_t)
{
   Fatal("Branch","Not implemented yet");
   return 0;
}

///////////////////////////////////////////////////////////////////////////////
/// Create a branch

TBranch * TTreeSQL::Branch(const char *name, void *address,
                           const char *leaflist, Int_t bufsize)
{
   Int_t nb = fBranches.GetEntriesFast();
   TBranch *branch;
   TString brName;

   for (int i=0;i<nb;i++) {
      branch = (TBranch*)fBranches.UncheckedAt(i);
      brName = branch->GetName();
      if (brName.CompareTo(name) == 0) {
         // Now if the branch exists in db, root gives a warning and exit
         // Dealing with duplicate branch has been done, but not tested yet.
         // So if you want to allow duplicate branch, just comment Fatal() line and uncomment commented
         // below Fatal() line

         Fatal("Branch()", "Duplicate branch!!!");

         /* Commented. If uncommented, should comment Fatal line.
         // this is a duplicate branch. So reset data structure memory address and return.
         branch->SetAddress(address);
         return branch;
         */
      }
   }

   return TTree::Branch(name, address, leaflist, bufsize);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the basket is properly setup

void TTreeSQL::CheckBasket(TBranch *branch)
{
   TBasketSQL* basket = (TBasketSQL *)branch->GetBasket(0);

   if (basket==0) {
      basket = (TBasketSQL*)CreateBasket(branch);
      if (basket==0) return;
      //++(branch->fNBaskets);
      branch->GetListOfBaskets()->AddAtAndExpand(basket,0);
   }
   TBuffer * buffer = basket->GetBufferRef();

   if(buffer == 0){
      std::vector<Int_t> *columns = GetColumnIndice(branch);
      if (columns) basket->CreateBuffer(branch->GetName(),"A", columns, branch, &fResult);
   }

   Int_t nb = branch->GetListOfBranches()->GetEntriesFast();
   for (int i=0;i<nb;i++) {
      TBranch * subbranch = (TBranch*)branch->GetListOfBranches()->UncheckedAt(i);
      if(subbranch) CheckBasket(subbranch);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the table has a column corresponding the branch
/// and that the resultset are properly setup

Bool_t TTreeSQL::CheckBranch(TBranch * tb)
{
   if (fServer==0) {
      return kFALSE;
   }
   TString leafName;
   TLeaf *leaf;
   Int_t nl;
   TString str = "";
   TString typeName = "";

   if (!tb) return kFALSE;

   TBasketSQL *basket = (TBasketSQL *)tb->GetBasket(0);
   if (!basket) return kFALSE;

   TSQLResult *rs = basket->GetResultSet();
   if (!rs) {
      Error("CheckBranch","%s has basket but no resultset yet",tb->GetName());
      return kFALSE;
   }

   nl = tb->GetNleaves();

   for(int j=0;j<nl;j++) {
      leaf = (TLeaf*)tb->GetListOfLeaves()->UncheckedAt(j);
      typeName = leaf->GetTypeName();
      typeName = ConvertTypeName(leaf->GetTypeName());
      leafName = leaf->GetName();
      str = "";
      str = tb->GetName();
      str += "__";
      str += leafName;

      for (int i=0; i< rs->GetFieldCount(); ++i) {
         if (str.CompareTo(rs->GetFieldName(i),TString::kIgnoreCase) == 0) return kTRUE;
      }
      // We assume that if ONE of the leaf is in the table, then ALL the leaf are in
      // the table.
      // TODO: this assumption is harmful if user changes branch structure while keep its name
      CreateBranch(str, typeName);
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the table exist in the database

Bool_t TTreeSQL::CheckTable(const TString &table) const
{
   if (fServer==0) return kFALSE;
   TSQLResult * tables = fServer->GetTables(fDB.Data(),table);
   if (!tables) return kFALSE;
   TSQLRow * row = 0;
   while( (row = tables->Next()) ) {
      if(table.CompareTo(row->GetField(0),TString::kIgnoreCase)==0){
         return kTRUE;
      }
   }
   // The table is a not a permanent table, let's see if it is a 'temporary' table
   Int_t before = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kFatal;
   TSQLResult *res = fServer->GetColumns(fDB.Data(),table);
   if (res) {
      delete res;
      return kTRUE;
   }
   gErrorIgnoreLevel = before;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert from ROOT typename to SQL typename

TString TTreeSQL::ConvertTypeName(const TString& typeName )
{
   TString tn = "";

   if(typeName == "Char_t"){
      tn = "TEXT";
   }
   else if(typeName == "Int_t") {
      tn = "INTEGER";
   }
   else if(typeName == "Short_t") {
      tn = "SMALLINT";
   }
   else if( typeName == "UShort_t") {
      tn = "SMALLINT UNSIGNED";
   }
   else if(typeName == "Float_t"){
      tn = "FLOAT";
   }
   else if(typeName == "Float16_t"){
      tn = "FLOAT";
   }
   else if(typeName == "Double_t"){
      tn = "DOUBLE";
   }
   else if(typeName == "Double32_t"){
      tn = "FLOAT";
   }
   else if(typeName == "UInt_t") {
      tn = "INT UNSIGNED";
   }
   else if( typeName == "Long_t") {
      tn = "INTEGER";
   }
   else if( typeName == "ULong_t") {
      tn = "INTEGER UNSIGNED";
   }
   else if( typeName == "Long64_t") {
      tn = "BIGINT";
   }
   else if( typeName == "ULong64_t") {
      tn = "BIGINT UNSIGNED";
   }
   else if( typeName == "Bool_t") {
      tn = "BOOL";
   }
   else if( typeName == "TString") {
      tn = "TEXT";
   }

   else {
      Error("ConvertTypeName","TypeName (%s) not found",typeName.Data());
      return "";
   }

   return tn;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TBasketSQL

TBasket * TTreeSQL::CreateBasket(TBranch * tb)
{
   if (fServer==0) {
      Error("CreateBasket","No TSQLServer specified");
      return 0;
   }
   std::vector<Int_t> *columnVec = GetColumnIndice(tb);
   if (columnVec) {
      return new TBasketSQL(tb->GetName(), tb->GetName(), tb,
                            &fResult, &fInsertQuery, columnVec, &fRow);
   } else {
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create the column(s) in the database that correspond to the branch/

void TTreeSQL::CreateBranch(const TString &branchName, const TString &typeName)
{
   if (fServer==0) {
      Error("CreateBranch","No TSQLServer specified");
      return;
   }
   TString alterSQL = "";
   alterSQL = "";
   alterSQL = "ALTER TABLE ";
   alterSQL += fTable.Data();
   alterSQL += " ADD ";
   alterSQL += branchName.Data();;
   alterSQL += " ";
   alterSQL += typeName;
   alterSQL += " ";

   fServer->Query(alterSQL);
}

////////////////////////////////////////////////////////////////////////////////
/// determine leaf description string

void TTreeSQL::CreateBranches()
{
   TList * columns = fTableInfo->GetColumns();
   if(!columns) return;

   TIter next(columns);

   TString branchName;
   TString type;
   TString leafName;
   TBranch * br = 0;
   TSQLColumnInfo * info;
   while ( (info = ((TSQLColumnInfo*) next()) ))
   {
      type = info->GetTypeName();
      branchName = info->GetName();


      Int_t pos;
      if ((pos=branchName.Index("__"))!=kNPOS) {
          leafName = branchName(pos+2,branchName.Length());
          branchName.Remove(pos);
      } else {
          leafName = branchName;
      }

      TString str;
      int i;
      unsigned ui;
      double d;
      float f;

      br = 0;

      if(type.CompareTo("varchar",TString::kIgnoreCase)==0 ||
         type.CompareTo("varchar2",TString::kIgnoreCase)==0 ||
         type.CompareTo("char",TString::kIgnoreCase)==0 ||
         type.CompareTo("longvarchar",TString::kIgnoreCase)==0 ||
         type.CompareTo("longvarbinary",TString::kIgnoreCase)==0 ||
         type.CompareTo("varbinary",TString::kIgnoreCase)==0 ||
         type.CompareTo("text",TString::kIgnoreCase )==0 ) {
         br = TTree::Branch(leafName,&str);

      }
      else if(type.CompareTo("int",TString::kIgnoreCase)==0 ){
         br = TTree::Branch(leafName,&i);
      }

      //Somehow it should be possible to special-case the time classes
      //but I think we'd need to create a new TSQLTime or something like that...
      else if( type.CompareTo("date",TString::kIgnoreCase)==0 ||
               type.CompareTo("time",TString::kIgnoreCase)==0 ||
               type.CompareTo("timestamp",TString::kIgnoreCase)==0 ||
               type.CompareTo("datetime",TString::kIgnoreCase)==0 ) {
         br = TTree::Branch(leafName,&str);

      }

      else if(type.CompareTo("bit",TString::kIgnoreCase)==0 ||
              type.CompareTo("tinyint",TString::kIgnoreCase)==0 ||
              type.CompareTo("smallint",TString::kIgnoreCase)==0 ) {
         br = TTree::Branch(leafName,&ui);
      }

      else if( type.CompareTo("decimal",TString::kIgnoreCase)==0 ||
               type.CompareTo("numeric",TString::kIgnoreCase)==0 ||
               type.CompareTo("double",TString::kIgnoreCase)==0 ||
               type.CompareTo("float",TString::kIgnoreCase)==0 )
      {
         br = TTree::Branch(leafName,&f);
      }
      else if( type.CompareTo("bigint",TString::kIgnoreCase)==0 ||
               type.CompareTo("real",TString::kIgnoreCase) == 0)
      {
         br = TTree::Branch(leafName,&d);
      }

      if (br == 0)
      {
         Error("CreateBranches", "Skipped %s", branchName.Data());
         continue;
      }

      br->ResetAddress();

      (br->GetBasketEntry())[0] = 0;
      (br->GetBasketEntry())[1] = fEntries;
      br->SetEntries(fEntries);

      //++(br->fNBaskets);
      br->GetListOfBaskets()->AddAtAndExpand(CreateBasket(br),0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create the database table corresponding to this TTree.

Bool_t TTreeSQL::CreateTable(const TString &table)
{
   if (fServer==0) {
      Error("CreateTable","No TSQLServer specified");
      return false;
   }
   Int_t i, j;
   TString branchName, leafName, typeName;
   TString createSQL, alterSQL, str;
   Int_t nb = fBranches.GetEntriesFast();
   Int_t nl = 0;

   TBranch *branch;
   TLeaf *leaf;

   for (i=0;i<nb;i++) {
      branch = (TBranch*)fBranches.UncheckedAt(i);
      branchName = branch->GetName();
      nl = branch->GetNleaves();
      for(j=0;j<nl;j++) {
         leaf = (TLeaf*)branch->GetListOfLeaves()->UncheckedAt(j);
         leafName = leaf->GetName();
         typeName = ConvertTypeName(leaf->GetTypeName());
         // length = leaf->GetLenStatic();

         if(i == 0 && j == 0) {
            createSQL = "";
            createSQL += "CREATE TABLE ";
            createSQL += table;
            createSQL += " (";
            createSQL += branchName;
            createSQL += "__";
            createSQL += leafName;
            createSQL += " ";
            createSQL += typeName;
            createSQL += " ";
            createSQL += ")";

            TSQLResult *sres = fServer->Query(createSQL.Data());
            if (!sres) {
               Error("CreateTable","May have failed");
               return false;
            }
         }
         else {
            str = "";
            str = branchName;
            str += "__";
            str += leafName;
            CreateBranch(str, typeName);
         } //else
      }  // inner for loop
   } // outer for loop
   // retrieve table to initialize fResult
   delete fResult;
   fResult = fServer->Query(fQuery.Data());
   return (fResult!=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Initializeation routine

void TTreeSQL::Init()
{
   fCurrentEntry = -1;

   GetEntries();

   delete fResult;
   fResult = fServer->Query(fQuery.Data());
   if(!fResult) return;

   if (fDB != "") {
      fServer->SelectDataBase(fDB);
   }
   fTableInfo = fServer->GetTableInfo(fTable);
   CreateBranches();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the information from the user object to the TTree

Int_t TTreeSQL::Fill()
{
   Int_t nb = fBranches.GetEntriesFast();
   TString typeName;
   TBranch *branch;

   if (fServer==0) return 0;

   if(!CheckTable(fTable.Data())) {
      if (!CreateTable(fTable.Data())) {
         return -1;
      }
   }

   PrepEntry(fEntries);

   for (int i=0;i<nb;i++) {
      branch = (TBranch*)fBranches.UncheckedAt(i);
      CheckBasket(branch);
   }

   if (!fBranchChecked) {
      for(int i=0;i<nb;i++) {
         branch = (TBranch*)fBranches.UncheckedAt(i);
         if (!CheckBranch(branch)) {
            Error("Fill","CheckBranch for %s failed",branch->GetName());
         }
      }
      fBranchChecked = kTRUE;
   }
   ResetQuery();

   TTree::Fill();

   if (fInsertQuery[fInsertQuery.Length()-1]!='(') {
      fInsertQuery.Remove(fInsertQuery.Length()-1);
      fInsertQuery += ")";
      TSQLResult *res = fServer?fServer->Query(fInsertQuery):0;

      if (res) {
         return res->GetRowCount();
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a vector of columns index corresponding to the
/// current SQL table and the branch given as argument
/// Returns 0 if no columns indices is found
/// Otherwise returns a pointer to a vector to be deleted by the caller

std::vector<Int_t> *TTreeSQL::GetColumnIndice(TBranch *branch)
{
   if (!CheckTable(fTable)) return 0;

   std::vector<Int_t> *columns = new std::vector<Int_t>;

   Int_t nl = branch->GetNleaves();

   std::vector<TString> names;

   TList *col_list = fTableInfo->GetColumns();
   if (col_list==0) {
      delete columns;
      return 0;
   }

   std::pair<TString,Int_t> value;

   TIter next(col_list);
   TSQLColumnInfo * cinfo;
   int rows = 0;
   while ((cinfo = (TSQLColumnInfo*) next())) {
      names.push_back( cinfo->GetName() );
      rows++;
   }

   for(int j=0;j<nl;j++) {

      Int_t col = -1;
      TLeaf *leaf = (TLeaf*)branch->GetListOfLeaves()->UncheckedAt(j);
      TString leafName = leaf->GetName();
      TString str;

      str = "";
      str = branch->GetName();
      str += "__";
      str += leafName;
      for (Int_t i=0;i<rows;++i) {
         if (str.CompareTo(names[i],TString::kIgnoreCase)==0) {
            col = i;
            break;
         }
      }
      if (col<0) {
         str = leafName;
         for (Int_t i=0;i<rows;++i) {
            if (str.CompareTo(names[i],TString::kIgnoreCase)==0) {
               col = i;
               break;
            }
         }
      }
      if(col>=0){
         columns->push_back(col);
      } else Error("GetColumnIndice","Error finding column %d %s",j,str.Data());
   }
   if (columns->empty()) {
      delete columns;
      return 0;
   } else
      return columns;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the number of rows in the database

Long64_t  TTreeSQL::GetEntries() const
{
   if (fServer==0) return GetEntriesFast();
   if (!CheckTable(fTable.Data())) return 0;

   TTreeSQL* thisvar = const_cast<TTreeSQL*>(this);

   // What if the user already started to call GetEntry
   // What about the initial value of fEntries is it really 0?

   TString counting = "select count(*) from " + fTable;
   TSQLResult *count = fServer->Query(counting);

   if (count==0) {
      thisvar->fEntries = 0;
   } else {
      TSQLRow * row = count->Next();
      if (row) {
         TString val = row->GetField(0);
         Long_t ret;
         sscanf(val.Data(), "%ld",&(ret) );
         thisvar->fEntries = ret;
      } else {
         thisvar->fEntries = 0;
      }
   }
   return fEntries;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of entries as of the last check.
/// Use GetEntries for a more accurate count.

Long64_t  TTreeSQL::GetEntriesFast()    const
{
   return fEntries;
}

////////////////////////////////////////////////////////////////////////////////
/// Load the data for the entry from the database.

Int_t TTreeSQL::GetEntry(Long64_t entry, Int_t getall)
{
   if (PrepEntry(entry)>=0) return TTree::GetEntry(entry,getall);
   else return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the tree to the load the specified entry.

Long64_t TTreeSQL::LoadTree(Long64_t entry)
{
   fReadEntry = entry;
   return PrepEntry(entry);
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure the server and result set are setup for the requested entry

Long64_t TTreeSQL::PrepEntry(Long64_t entry)
{
   if (entry < 0 || entry >= fEntries || fServer==0) return 0;
   fReadEntry = entry;

   if(entry == fCurrentEntry) return entry;

   if(entry < fCurrentEntry || fResult==0){
      delete fResult;
      fResult = fServer->Query(fQuery.Data());
      fCurrentEntry = -1;
   }

   Bool_t reset = false;
   while ( fResult && fCurrentEntry < entry ) {
      ++fCurrentEntry;
      delete fRow;
      fRow = fResult->Next();
      if (fRow==0 && !reset) {
         delete fResult;
         fResult = fServer->Query(fQuery.Data());
         fCurrentEntry = -1;
         reset = true;
      }
   }
   if (fRow==0) return -1;
   return entry;
}

//______________________________________________________________________________
// void TTreeSQL::LoadNumberEntries()
// {
//    R__ASSERT(0);

//    fResult =    fServer->Query(fQuery.Data());
//    fEntries=0;

//    while(fResult->Next()){
//       fEntries++;
//    }
//    fResult =    fServer->Query(fQuery.Data());
// }

////////////////////////////////////////////////////////////////////////////////
///  Refresh contents of this Tree and its branches from the current
///  Tree status in the database
///  One can call this function in case the Tree on its file is being
///  updated by another process

void TTreeSQL::Refresh()
{
   // Note : something to be done?
   GetEntries(); // Re-load the number of entries
   fCurrentEntry = -1;
   delete fResult; fResult = 0;
   delete fRow; fRow = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the internal query

void TTreeSQL::ResetQuery()
{
   fInsertQuery = "INSERT INTO " + fTable + " VALUES (";
}


////////////////////////////////////////////////////////////////////////////////
// Destructor

TTreeSQL::~TTreeSQL()
{
   delete fTableInfo;
   delete fResult;
   delete fRow;
}
