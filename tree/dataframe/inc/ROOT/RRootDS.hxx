// Author: Enrico Guiraud, Danilo Piparo CERN  9/2017

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RROOTTDS
#define ROOT_RROOTTDS

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDataSource.hxx"
#include <TChain.h>

#include <memory>

namespace ROOT {

namespace Internal {

namespace RDF {

class RRootDS final : public ROOT::RDF::RDataSource {
private:
   unsigned int fNSlots = 0U;
   std::string fTreeName;
   std::string fFileNameGlob;
   mutable TChain fModelChain; // Mutable needed for getting the column type name
   std::vector<double *> fAddressesToFree;
   std::vector<std::string> fListOfBranches;
   std::vector<std::pair<ULong64_t, ULong64_t>> fEntryRanges;
   std::vector<std::vector<void *>> fBranchAddresses; // first container-> slot, second -> column;
   std::vector<std::unique_ptr<TChain>> fChains;

   std::vector<void *> GetColumnReadersImpl(std::string_view, const std::type_info &);

protected:
   std::string AsString() { return "ROOT data source"; };

public:
   RRootDS(std::string_view treeName, std::string_view fileNameGlob);
   ~RRootDS();
   std::string GetTypeName(std::string_view colName) const;
   const std::vector<std::string> &GetColumnNames() const;
   bool HasColumn(std::string_view colName) const;
   void InitSlot(unsigned int slot, ULong64_t firstEntry);
   void FinaliseSlot(unsigned int slot);
   std::vector<std::pair<ULong64_t, ULong64_t>> GetEntryRanges();
   bool SetEntry(unsigned int slot, ULong64_t entry);
   void SetNSlots(unsigned int nSlots);
   void Initialise();
   std::string GetLabel();
};

RDataFrame MakeRootDataFrame(std::string_view treeName, std::string_view fileNameGlob);

} // ns RDF

} // ns Internal

} // ns ROOT

#endif
