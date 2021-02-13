// Author: Enrico Guiraud, Danilo Piparo CERN  9/2017

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RTRIVIALTDS
#define ROOT_RTRIVIALTDS

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataSource.hxx"

namespace ROOT {

namespace RDF {

class RTrivialDS final : public ROOT::RDF::RDataSource {
private:
   unsigned int fNSlots = 0U;
   ULong64_t fSize = 0ULL;
   bool fSkipEvenEntries = false;
   std::vector<std::pair<ULong64_t, ULong64_t>> fEntryRanges;
   std::vector<std::string> fColNames{"col0"};
   std::vector<ULong64_t> fCounter;
   std::vector<ULong64_t *> fCounterAddr;
   std::vector<void *> GetColumnReadersImpl(std::string_view name, const std::type_info &);

protected:
   std::string AsString() { return "trivial data source"; };

public:
   RTrivialDS(ULong64_t size, bool skipEvenEntries = false);
   /// This ctor produces a data-source that returns infinite entries
   RTrivialDS();
   ~RTrivialDS();
   const std::vector<std::string> &GetColumnNames() const;
   bool HasColumn(std::string_view colName) const;
   std::string GetTypeName(std::string_view) const;
   std::vector<std::pair<ULong64_t, ULong64_t>> GetEntryRanges();
   bool SetEntry(unsigned int slot, ULong64_t entry);
   void SetNSlots(unsigned int nSlots);
   void Initialise();
   std::string GetLabel();
};

// Make a RDF wrapping a RTrivialDS with the specified amount of entries
RInterface<RDFDetail::RLoopManager, RTrivialDS> MakeTrivialDataFrame(ULong64_t size, bool skipEvenEntries = false);
// Make a RDF wrapping a RTrivialDS with infinite entries
RInterface<RDFDetail::RLoopManager, RTrivialDS> MakeTrivialDataFrame();

} // ns RDF

} // ns ROOT

#endif
