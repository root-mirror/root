#ifndef ROOT_RNONCOPIABLECOLUMNDS
#define ROOT_RNONCOPIABLECOLUMNDS

#include "ROOT/RDataSource.hxx"
#include <string>
#include <vector>

class RNonCopiable {
public:
   using type = unsigned int;
   unsigned int fValue = 42;
   RNonCopiable(const RNonCopiable &) = delete;
   RNonCopiable() = default;
};

class RNonCopiableColumnDS final : public ROOT::RDF::RDataSource {
private:
   std::vector<std::pair<ULong64_t, ULong64_t>> fEntryRanges = {{0UL, 1UL}};
   std::vector<std::string> fColNames{fgColumnName};
   RNonCopiable fNonCopiable;
   RNonCopiable *fCounterAddr = &fNonCopiable;
   std::vector<void *> GetColumnReadersImpl(std::string_view, const std::type_info &)
   {
      std::vector<void *> ret{(void *)(&fCounterAddr)};
      return ret;
   }

public:
   using NonCopiable_t = RNonCopiable;
   constexpr const static auto fgColumnName = "nonCopiable";
   RNonCopiableColumnDS(){};
   ~RNonCopiableColumnDS(){};
   const std::vector<std::string> &GetColumnNames() const { return fColNames; };
   bool HasColumn(std::string_view colName) const { return colName == fColNames[0]; };
   std::string GetTypeName(std::string_view) const { return "RNonCopiable"; };
   std::vector<std::pair<ULong64_t, ULong64_t>> GetEntryRanges()
   {
      auto entryRanges(std::move(fEntryRanges)); // empty fEntryRanges
      return entryRanges;
   };
   bool SetEntry(unsigned int, ULong64_t){ return true;};
   void SetNSlots(unsigned int){};
   std::string GetLabel(){
      return "NonCopiableColumnDS";
   }
};

#endif
