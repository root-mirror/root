// Author: Jakob Blomer CERN  07/2018

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RSQLITEDS
#define ROOT_RSQLITEDS

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDataSource.hxx"
#include "ROOT/RStringView.hxx"

#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace ROOT {

namespace RDF {

namespace Internal {
// Members are defined in RSqliteDS.cxx in order to not pullute this header file with sqlite3.h
struct RSqliteDSDataSet;
}

// clang-format off
/**
\class ROOT::RDF::RSqliteDS
\ingroup dataframe
\brief RSqliteDS is an RDF data source implementation for SQL result sets from sqlite3 files.

The RSqliteDS is able to feed an RDataFrame with data from a SQlite SELECT query. One can use it like

    auto rdf = ROOT::RDF::MakeSqliteDataFrame("/path/to/file.sqlite", "select name from table");
    auto h = rdf.Define("lName", "name.length()").Histo1D("lName");

The data source has to provide column types for all the columns. Determining column types in SQlite is tricky
as it is dynamically typed and in principle each row can have different column types. The following heuristics
is used:

  - If a table column is queried as is ("SELECT colname FROM table"), the default/declared column type is taken.
  - For expressions ("SELECT 1+1 FROM table"), the type of the first row of the result set determines the column type.
    That can result in a column to be of thought of type NULL where subsequent rows actually have meaningful values.
    The provided SELECT query can be used to avoid such ambiguities.
*/
class RSqliteDS final : public ROOT::RDF::RDataSource {
private:
   // clang-format off
   /// All the types known to SQlite. Changes require changing fgTypeNames, too.
   enum class ETypes {
      kInteger,
      kReal,
      kText,
      kBlob,
      kNull
   };
   // clang-format on

   /// Used to hold a single "cell" of the SELECT query's result table. Can be changed to std::variant once available.
   struct Value_t {
      explicit Value_t(ETypes type);

      ETypes fType;
      bool fIsActive; ///< Not all columns of the query are necessarily used by the RDF. Allows for skipping them.
      Long64_t fInteger;
      double fReal;
      std::string fText;
      std::vector<unsigned char> fBlob;
      void *fNull;
      void *fPtr; ///< Points to one of the values; an address to this pointer is returned by GetColumnReadersImpl.
   };

   void SqliteError(int errcode);

   std::unique_ptr<Internal::RSqliteDSDataSet> fDataSet;
   unsigned int fNSlots;
   ULong64_t fNRow;
   std::vector<std::string> fColumnNames;
   std::vector<ETypes> fColumnTypes;
   /// The data source is inherently single-threaded and returns only one row at a time. This vector holds the results.
   std::vector<Value_t> fValues;

   // clang-format off
   /// Corresponds to the types defined in ETypes.
   static constexpr char const *fgTypeNames[] = {
      "Long64_t",
      "double",
      "std::string",
      "std::vector<unsigned char>",
      "void *"
   };
   // clang-format on

public:
   RSqliteDS(const std::string &fileName, const std::string &query);
   ~RSqliteDS();
   void SetNSlots(unsigned int nSlots) final;
   const std::vector<std::string> &GetColumnNames() const final;
   bool HasColumn(std::string_view colName) const final;
   std::string GetTypeName(std::string_view colName) const final;
   std::vector<std::pair<ULong64_t, ULong64_t>> GetEntryRanges() final;
   bool SetEntry(unsigned int slot, ULong64_t entry) final;
   void Initialise() final;
   std::string GetLabel() final;

protected:
   Record_t GetColumnReadersImpl(std::string_view name, const std::type_info &) final;
};

RDataFrame MakeSqliteDataFrame(std::string_view fileName, std::string_view query);

} // namespace RDF

} // namespace ROOT

#endif
