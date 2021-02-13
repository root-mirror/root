// Author: Enric Tejedor CERN  10/2017

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RCSVTDS
#define ROOT_RCSVTDS

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDataSource.hxx"

#include <cstdint>
#include <deque>
#include <list>
#include <map>
#include <memory>
#include <vector>

#include <TRegexp.h>

namespace ROOT {

namespace Internal {
class RRawFile;
}

namespace RDF {

class RCsvDS final : public ROOT::RDF::RDataSource {

private:
   // Possible values are d, b, l, s. This is possible only because we treat double, bool, Long64_t and string
   using ColType_t = char;
   static const std::map<ColType_t, std::string> fgColTypeMap;

   // Regular expressions for type inference
   static const TRegexp fgIntRegex, fgDoubleRegex1, fgDoubleRegex2, fgDoubleRegex3, fgTrueRegex, fgFalseRegex;

   std::uint64_t fDataPos = 0;
   bool fReadHeaders = false;
   unsigned int fNSlots = 0U;
   std::unique_ptr<ROOT::Internal::RRawFile> fCsvFile;
   const char fDelimiter;
   const Long64_t fLinesChunkSize;
   ULong64_t fEntryRangesRequested = 0ULL;
   ULong64_t fProcessedLines = 0ULL; // marks the progress of the consumption of the csv lines
   std::vector<std::string> fHeaders;
   std::map<std::string, ColType_t> fColTypes;
   std::list<ColType_t> fColTypesList;
   std::vector<std::vector<void *>> fColAddresses;         // fColAddresses[column][slot]
   std::vector<Record_t> fRecords;                         // fRecords[entry][column]
   std::vector<std::vector<double>> fDoubleEvtValues;      // one per column per slot
   std::vector<std::vector<Long64_t>> fLong64EvtValues;    // one per column per slot
   std::vector<std::vector<std::string>> fStringEvtValues; // one per column per slot
   // This must be a deque to avoid the specialisation vector<bool>. This would not
   // work given that the pointer to the boolean in that case cannot be taken
   std::vector<std::deque<bool>> fBoolEvtValues; // one per column per slot

   void FillHeaders(const std::string &);
   void FillRecord(const std::string &, Record_t &);
   void GenerateHeaders(size_t);
   std::vector<void *> GetColumnReadersImpl(std::string_view, const std::type_info &);
   void InferColTypes(std::vector<std::string> &);
   void InferType(const std::string &, unsigned int);
   std::vector<std::string> ParseColumns(const std::string &);
   size_t ParseValue(const std::string &, std::vector<std::string> &, size_t);
   ColType_t GetType(std::string_view colName) const;

protected:
   std::string AsString();

public:
   RCsvDS(std::string_view fileName, bool readHeaders = true, char delimiter = ',', Long64_t linesChunkSize = -1LL);
   void Finalise();
   void FreeRecords();
   ~RCsvDS();
   const std::vector<std::string> &GetColumnNames() const;
   std::vector<std::pair<ULong64_t, ULong64_t>> GetEntryRanges();
   std::string GetTypeName(std::string_view colName) const;
   bool HasColumn(std::string_view colName) const;
   bool SetEntry(unsigned int slot, ULong64_t entry);
   void SetNSlots(unsigned int nSlots);
   std::string GetLabel();
};

////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Factory method to create a CSV RDataFrame.
/// \param[in] fileName Path of the CSV file.
/// \param[in] readHeaders `true` if the CSV file contains headers as first row, `false` otherwise
///                        (default `true`).
/// \param[in] delimiter Delimiter character (default ',').
RDataFrame MakeCsvDataFrame(std::string_view fileName, bool readHeaders = true, char delimiter = ',',
                            Long64_t linesChunkSize = -1LL);

} // ns RDF

} // ns ROOT

#endif
