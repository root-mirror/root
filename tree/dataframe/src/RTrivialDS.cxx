#include <ROOT/RDF/Utils.hxx>
#include <ROOT/TSeq.hxx>
#include <ROOT/RTrivialDS.hxx>
#include <ROOT/RMakeUnique.hxx>
#include <TError.h>

#include <limits>

namespace ROOT {

namespace RDF {

std::vector<void *> RTrivialDS::GetColumnReadersImpl(std::string_view, const std::type_info &ti)
{
   // We know we have only one column and that it's holding ULong64_t's
   if (ti != typeid(ULong64_t)) {
      throw std::runtime_error("The type specified for the column \"col0\" is not ULong64_t.");
   }
   std::vector<void *> ret;
   for (auto i : ROOT::TSeqU(fNSlots)) {
      fCounterAddr[i] = &fCounter[i];
      ret.emplace_back((void *)(&fCounterAddr[i]));
   }
   return ret;
}

RTrivialDS::RTrivialDS(ULong64_t size, bool skipEvenEntries) : fSize(size), fSkipEvenEntries(skipEvenEntries)
{
}

RTrivialDS::RTrivialDS() : fSize(std::numeric_limits<ULong64_t>::max()), fSkipEvenEntries(false)
{
}

RTrivialDS::~RTrivialDS()
{
}

const std::vector<std::string> &RTrivialDS::GetColumnNames() const
{
   return fColNames;
}

bool RTrivialDS::HasColumn(std::string_view colName) const
{
   return colName == fColNames[0];
}

std::string RTrivialDS::GetTypeName(std::string_view) const
{
   return "ULong64_t";
}

std::vector<std::pair<ULong64_t, ULong64_t>> RTrivialDS::GetEntryRanges()
{
   if (fSize == std::numeric_limits<ULong64_t>::max()) {
      auto currentEntry = *std::max_element(fCounter.begin(), fCounter.end());
      // infinite source, just make some ranges up
      std::vector<std::pair<ULong64_t, ULong64_t>> ranges(fNSlots);
      for (auto &range : ranges) {
         range = std::make_pair(currentEntry, currentEntry + 10);
         currentEntry += 10;
      }
      return ranges;
   }

   // empty fEntryRanges so we'll return an empty vector on subsequent calls
   auto ranges = std::move(fEntryRanges);
   return ranges;
}

bool RTrivialDS::SetEntry(unsigned int slot, ULong64_t entry)
{
   if (fSkipEvenEntries && 0 == entry % 2) {
      return false;
   }
   fCounter[slot] = entry;
   return true;
}

void RTrivialDS::SetNSlots(unsigned int nSlots)
{
   R__ASSERT(0U == fNSlots && "Setting the number of slots even if the number of slots is different from zero.");

   fNSlots = nSlots;
   fCounter.resize(fNSlots);
   fCounterAddr.resize(fNSlots);
}

void RTrivialDS::Initialise()
{
   if (fSize == std::numeric_limits<ULong64_t>::max()) {
      // infinite source, nothing to do here
      return;
   }

   // initialize fEntryRanges
   const auto chunkSize = fSize / fNSlots;
   auto start = 0UL;
   auto end = 0UL;
   for (auto i : ROOT::TSeqUL(fNSlots)) {
      start = end;
      end += chunkSize;
      fEntryRanges.emplace_back(start, end);
      (void)i;
   }
   // TODO: redistribute reminder to all slots
   fEntryRanges.back().second += fSize % fNSlots;
}

std::string RTrivialDS::GetLabel()
{
   return "TrivialDS";
}

RInterface<RDFDetail::RLoopManager, RTrivialDS> MakeTrivialDataFrame(ULong64_t size, bool skipEvenEntries)
{
   auto lm = std::make_unique<RDFDetail::RLoopManager>(std::make_unique<RTrivialDS>(size, skipEvenEntries),
                                                       RDFInternal::ColumnNames_t{});
   return RInterface<RDFDetail::RLoopManager, RTrivialDS>(std::move(lm));
}

RInterface<RDFDetail::RLoopManager, RTrivialDS> MakeTrivialDataFrame()
{
   auto lm = std::make_unique<RDFDetail::RLoopManager>(std::make_unique<RTrivialDS>(), RDFInternal::ColumnNames_t{});
   return RInterface<RDFDetail::RLoopManager, RTrivialDS>(std::move(lm));
}

} // ns RDF

} // ns ROOT
