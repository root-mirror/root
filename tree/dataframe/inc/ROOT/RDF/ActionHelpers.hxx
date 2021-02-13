/**
 \file ROOT/RDF/ActionHelpers.hxx
 \ingroup dataframe
 \author Enrico Guiraud, CERN
 \author Danilo Piparo, CERN
 \date 2016-12
 \author Vincenzo Eduardo Padulano
 \date 2020-06
*/

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RDFOPERATIONS
#define ROOT_RDFOPERATIONS

#include "Compression.h"
#include "ROOT/RIntegerSequence.hxx"
#include "ROOT/RStringView.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/TBufferMerger.hxx" // for SnapshotHelper
#include "ROOT/RDF/RCutFlowReport.hxx"
#include "ROOT/RDF/Utils.hxx"
#include "ROOT/RMakeUnique.hxx"
#include "ROOT/RSnapshotOptions.hxx"
#include "ROOT/TypeTraits.hxx"
#include "ROOT/RDF/RDisplay.hxx"
#include "RtypesCore.h"
#include "TBranch.h"
#include "TClassEdit.h"
#include "TClassRef.h"
#include "TDirectory.h"
#include "TFile.h" // for SnapshotHelper
#include "TH1.h"
#include "TGraph.h"
#include "TLeaf.h"
#include "TObject.h"
#include "TTree.h"
#include "TTreeReader.h" // for SnapshotHelper
#include "ROOT/RDF/RMergeableValue.hxx"

#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <iomanip>
#include <numeric> // std::accumulate in MeanHelper

/// \cond HIDDEN_SYMBOLS

namespace ROOT {
namespace Detail {
namespace RDF {
template <typename Helper>
class RActionImpl {
public:
   virtual ~RActionImpl() = default;
   // call Helper::FinalizeTask if present, do nothing otherwise
   template <typename T = Helper>
   auto CallFinalizeTask(unsigned int slot) -> decltype(&T::FinalizeTask, void())
   {
      static_cast<Helper *>(this)->FinalizeTask(slot);
   }

   template <typename... Args>
   void CallFinalizeTask(unsigned int, Args...) {}

   // Helper functions for RMergeableValue
   virtual std::unique_ptr<RMergeableValueBase> GetMergeableValue() const
   {
      throw std::logic_error("`GetMergeableValue` is not implemented for this type of action.");
   }
};

} // namespace RDF
} // namespace Detail

namespace Internal {
namespace RDF {
using namespace ROOT::TypeTraits;
using namespace ROOT::VecOps;
using namespace ROOT::RDF;
using namespace ROOT::Detail::RDF;

using Hist_t = ::TH1D;

/// The container type for each thread's partial result in an action helper
// We have to avoid to instantiate std::vector<bool> as that makes it impossible to return a reference to one of
// the thread-local results. In addition, a common definition for the type of the container makes it easy to swap
// the type of the underlying container if e.g. we see problems with false sharing of the thread-local results..
template <typename T>
using Results = typename std::conditional<std::is_same<T, bool>::value, std::deque<T>, std::vector<T>>::type;

template <typename F>
class ForeachSlotHelper : public RActionImpl<ForeachSlotHelper<F>> {
   F fCallable;

public:
   using ColumnTypes_t = RemoveFirstParameter_t<typename CallableTraits<F>::arg_types>;
   ForeachSlotHelper(F &&f) : fCallable(f) {}
   ForeachSlotHelper(ForeachSlotHelper &&) = default;
   ForeachSlotHelper(const ForeachSlotHelper &) = delete;

   void InitTask(TTreeReader *, unsigned int) {}

   template <typename... Args>
   void Exec(unsigned int slot, Args &&... args)
   {
      // check that the decayed types of Args are the same as the branch types
      static_assert(std::is_same<TypeList<typename std::decay<Args>::type...>, ColumnTypes_t>::value, "");
      fCallable(slot, std::forward<Args>(args)...);
   }

   void Initialize() { /* noop */}

   void Finalize() { /* noop */}

   std::string GetActionName() { return "ForeachSlot"; }
};

class CountHelper : public RActionImpl<CountHelper> {
   const std::shared_ptr<ULong64_t> fResultCount;
   Results<ULong64_t> fCounts;

public:
   using ColumnTypes_t = TypeList<>;
   CountHelper(const std::shared_ptr<ULong64_t> &resultCount, const unsigned int nSlots);
   CountHelper(CountHelper &&) = default;
   CountHelper(const CountHelper &) = delete;
   void InitTask(TTreeReader *, unsigned int) {}
   void Exec(unsigned int slot);
   void Initialize() { /* noop */}
   void Finalize();

   // Helper functions for RMergeableValue
   std::unique_ptr<RMergeableValueBase> GetMergeableValue() const final
   {
      return std::make_unique<RMergeableCount>(*fResultCount);
   }

   ULong64_t &PartialUpdate(unsigned int slot);

   std::string GetActionName() { return "Count"; }
};

template <typename ProxiedVal_t>
class ReportHelper : public RActionImpl<ReportHelper<ProxiedVal_t>> {
   const std::shared_ptr<RCutFlowReport> fReport;
   // Here we have a weak pointer since we need to keep track of the validity
   // of the proxied node. It can happen that the user does not trigger the
   // event loop by looking into the RResultPtr and the chain goes out of scope
   // before the Finalize method is invoked.
   std::weak_ptr<ProxiedVal_t> fProxiedWPtr;
   bool fReturnEmptyReport;

public:
   using ColumnTypes_t = TypeList<>;
   ReportHelper(const std::shared_ptr<RCutFlowReport> &report, const std::shared_ptr<ProxiedVal_t> &pp, bool emptyRep)
      : fReport(report), fProxiedWPtr(pp), fReturnEmptyReport(emptyRep){};
   ReportHelper(ReportHelper &&) = default;
   ReportHelper(const ReportHelper &) = delete;
   void InitTask(TTreeReader *, unsigned int) {}
   void Exec(unsigned int /* slot */) {}
   void Initialize() { /* noop */}
   void Finalize()
   {
      // We need the weak_ptr in order to avoid crashes at tear down
      if (!fReturnEmptyReport && !fProxiedWPtr.expired())
         fProxiedWPtr.lock()->Report(*fReport);
   }

   std::string GetActionName() { return "Report"; }
};

class FillHelper : public RActionImpl<FillHelper> {
   // this sets a total initial size of 16 MB for the buffers (can increase)
   static constexpr unsigned int fgTotalBufSize = 2097152;
   using BufEl_t = double;
   using Buf_t = std::vector<BufEl_t>;

   std::vector<Buf_t> fBuffers;
   std::vector<Buf_t> fWBuffers;
   const std::shared_ptr<Hist_t> fResultHist;
   unsigned int fNSlots;
   unsigned int fBufSize;
   /// Histograms containing "snapshots" of partial results. Non-null only if a registered callback requires it.
   Results<std::unique_ptr<Hist_t>> fPartialHists;
   Buf_t fMin;
   Buf_t fMax;

   void UpdateMinMax(unsigned int slot, double v);

public:
   FillHelper(const std::shared_ptr<Hist_t> &h, const unsigned int nSlots);
   FillHelper(FillHelper &&) = default;
   FillHelper(const FillHelper &) = delete;
   void InitTask(TTreeReader *, unsigned int) {}
   void Exec(unsigned int slot, double v);
   void Exec(unsigned int slot, double v, double w);

   template <typename T, typename std::enable_if<IsDataContainer<T>::value || std::is_same<T, std::string>::value, int>::type = 0>
   void Exec(unsigned int slot, const T &vs)
   {
      auto &thisBuf = fBuffers[slot];
      // range-based for results in warnings on some compilers due to vector<bool>'s custom reference type
      for (auto v = vs.begin(); v != vs.end(); ++v) {
         UpdateMinMax(slot, *v);
         thisBuf.emplace_back(*v); // TODO: Can be optimised in case T == BufEl_t
      }
   }

   template <typename T, typename W,
             typename std::enable_if<IsDataContainer<T>::value && IsDataContainer<W>::value, int>::type = 0>
   void Exec(unsigned int slot, const T &vs, const W &ws)
   {
      auto &thisBuf = fBuffers[slot];

      for (auto &v : vs) {
         UpdateMinMax(slot, v);
         thisBuf.emplace_back(v);
      }

      auto &thisWBuf = fWBuffers[slot];
      for (auto &w : ws) {
         thisWBuf.emplace_back(w); // TODO: Can be optimised in case T == BufEl_t
      }
   }

   template <typename T, typename W,
             typename std::enable_if<IsDataContainer<T>::value && !IsDataContainer<W>::value, int>::type = 0>
   void Exec(unsigned int slot, const T &vs, const W w)
   {
      auto &thisBuf = fBuffers[slot];
      for (auto &v : vs) {
         UpdateMinMax(slot, v);
         thisBuf.emplace_back(v); // TODO: Can be optimised in case T == BufEl_t
      }

      auto &thisWBuf = fWBuffers[slot];
      thisWBuf.insert(thisWBuf.end(), vs.size(), w);
   }

   // ROOT-10092: Filling with a scalar as first column and a collection as second is not supported
   template <typename T, typename W,
             typename std::enable_if<IsDataContainer<W>::value && !IsDataContainer<T>::value, int>::type = 0>
   void Exec(unsigned int, const T &, const W &)
   {
      throw std::runtime_error(
        "Cannot fill object if the type of the first column is a scalar and the one of the second a container.");
   }

   Hist_t &PartialUpdate(unsigned int);

   void Initialize() { /* noop */}

   void Finalize();

   // Helper functions for RMergeableValue
   std::unique_ptr<RMergeableValueBase> GetMergeableValue() const final
   {
      return std::make_unique<RMergeableFill<Hist_t>>(*fResultHist);
   }

   std::string GetActionName() { return "Fill"; }
};

extern template void FillHelper::Exec(unsigned int, const std::vector<float> &);
extern template void FillHelper::Exec(unsigned int, const std::vector<double> &);
extern template void FillHelper::Exec(unsigned int, const std::vector<char> &);
extern template void FillHelper::Exec(unsigned int, const std::vector<int> &);
extern template void FillHelper::Exec(unsigned int, const std::vector<unsigned int> &);
extern template void FillHelper::Exec(unsigned int, const std::vector<float> &, const std::vector<float> &);
extern template void FillHelper::Exec(unsigned int, const std::vector<double> &, const std::vector<double> &);
extern template void FillHelper::Exec(unsigned int, const std::vector<char> &, const std::vector<char> &);
extern template void FillHelper::Exec(unsigned int, const std::vector<int> &, const std::vector<int> &);
extern template void
FillHelper::Exec(unsigned int, const std::vector<unsigned int> &, const std::vector<unsigned int> &);

template <typename HIST = Hist_t>
class FillParHelper : public RActionImpl<FillParHelper<HIST>> {
   std::vector<HIST *> fObjects;

public:
   FillParHelper(FillParHelper &&) = default;
   FillParHelper(const FillParHelper &) = delete;

   FillParHelper(const std::shared_ptr<HIST> &h, const unsigned int nSlots) : fObjects(nSlots, nullptr)
   {
      fObjects[0] = h.get();
      // Initialise all other slots
      for (unsigned int i = 1; i < nSlots; ++i) {
         fObjects[i] = new HIST(*fObjects[0]);
         if (auto objAsHist = dynamic_cast<TH1*>(fObjects[i])) {
            objAsHist->SetDirectory(nullptr);
         }
      }
   }

   void InitTask(TTreeReader *, unsigned int) {}

   void Exec(unsigned int slot, double x0) // 1D histos
   {
      fObjects[slot]->Fill(x0);
   }

   void Exec(unsigned int slot, double x0, double x1) // 1D weighted and 2D histos
   {
      fObjects[slot]->Fill(x0, x1);
   }

   void Exec(unsigned int slot, double x0, double x1, double x2) // 2D weighted and 3D histos
   {
      fObjects[slot]->Fill(x0, x1, x2);
   }

   void Exec(unsigned int slot, double x0, double x1, double x2, double x3) // 3D weighted histos
   {
      fObjects[slot]->Fill(x0, x1, x2, x3);
   }

   template <typename X0, typename std::enable_if<IsDataContainer<X0>::value || std::is_same<X0, std::string>::value, int>::type = 0>
   void Exec(unsigned int slot, const X0 &x0s)
   {
      auto thisSlotH = fObjects[slot];
      for (auto x0 = x0s.begin(); x0 != x0s.end(); x0++) {
         thisSlotH->Fill(*x0); // TODO: Can be optimised in case T == vector<double>
      }
   }

   // ROOT-10092: Filling with a scalar as first column and a collection as second is not supported
   template <typename X0, typename X1,
             typename std::enable_if<IsDataContainer<X1>::value && !IsDataContainer<X0>::value, int>::type = 0>
   void Exec(unsigned int , const X0 &, const X1 &)
   {
      throw std::runtime_error(
        "Cannot fill object if the type of the first column is a scalar and the one of the second a container.");
   }

   template <typename X0, typename X1,
             typename std::enable_if<IsDataContainer<X0>::value && IsDataContainer<X1>::value, int>::type = 0>
   void Exec(unsigned int slot, const X0 &x0s, const X1 &x1s)
   {
      auto thisSlotH = fObjects[slot];
      if (x0s.size() != x1s.size()) {
         throw std::runtime_error("Cannot fill histogram with values in containers of different sizes.");
      }
      auto x0sIt = std::begin(x0s);
      const auto x0sEnd = std::end(x0s);
      auto x1sIt = std::begin(x1s);
      for (; x0sIt != x0sEnd; x0sIt++, x1sIt++) {
         thisSlotH->Fill(*x0sIt, *x1sIt); // TODO: Can be optimised in case T == vector<double>
      }
   }

   template <typename X0, typename W,
             typename std::enable_if<IsDataContainer<X0>::value && !IsDataContainer<W>::value, int>::type = 0>
   void Exec(unsigned int slot, const X0 &x0s, const W w)
   {
      auto thisSlotH = fObjects[slot];
      for (auto &&x : x0s) {
         thisSlotH->Fill(x, w);
      }
   }

   template <typename X0, typename X1, typename X2,
             typename std::enable_if<IsDataContainer<X0>::value && IsDataContainer<X1>::value && IsDataContainer<X2>::value,
                                     int>::type = 0>
   void Exec(unsigned int slot, const X0 &x0s, const X1 &x1s, const X2 &x2s)
   {
      auto thisSlotH = fObjects[slot];
      if (!(x0s.size() == x1s.size() && x1s.size() == x2s.size())) {
         throw std::runtime_error("Cannot fill histogram with values in containers of different sizes.");
      }
      auto x0sIt = std::begin(x0s);
      const auto x0sEnd = std::end(x0s);
      auto x1sIt = std::begin(x1s);
      auto x2sIt = std::begin(x2s);
      for (; x0sIt != x0sEnd; x0sIt++, x1sIt++, x2sIt++) {
         thisSlotH->Fill(*x0sIt, *x1sIt, *x2sIt); // TODO: Can be optimised in case T == vector<double>
      }
   }

   template <typename X0, typename X1, typename W,
             typename std::enable_if<IsDataContainer<X0>::value && IsDataContainer<X1>::value && !IsDataContainer<W>::value,
                                     int>::type = 0>
   void Exec(unsigned int slot, const X0 &x0s, const X1 &x1s, const W w)
   {
      auto thisSlotH = fObjects[slot];
      if (x0s.size() != x1s.size()) {
         throw std::runtime_error("Cannot fill histogram with values in containers of different sizes.");
      }
      auto x0sIt = std::begin(x0s);
      const auto x0sEnd = std::end(x0s);
      auto x1sIt = std::begin(x1s);
      for (; x0sIt != x0sEnd; x0sIt++, x1sIt++) {
         thisSlotH->Fill(*x0sIt, *x1sIt, w); // TODO: Can be optimised in case T == vector<double>
      }
   }

   template <typename X0, typename X1, typename X2, typename X3,
             typename std::enable_if<IsDataContainer<X0>::value && IsDataContainer<X1>::value && IsDataContainer<X2>::value &&
                                        IsDataContainer<X3>::value,
                                     int>::type = 0>
   void Exec(unsigned int slot, const X0 &x0s, const X1 &x1s, const X2 &x2s, const X3 &x3s)
   {
      auto thisSlotH = fObjects[slot];
      if (!(x0s.size() == x1s.size() && x1s.size() == x2s.size() && x1s.size() == x3s.size())) {
         throw std::runtime_error("Cannot fill histogram with values in containers of different sizes.");
      }
      auto x0sIt = std::begin(x0s);
      const auto x0sEnd = std::end(x0s);
      auto x1sIt = std::begin(x1s);
      auto x2sIt = std::begin(x2s);
      auto x3sIt = std::begin(x3s);
      for (; x0sIt != x0sEnd; x0sIt++, x1sIt++, x2sIt++, x3sIt++) {
         thisSlotH->Fill(*x0sIt, *x1sIt, *x2sIt, *x3sIt); // TODO: Can be optimised in case T == vector<double>
      }
   }

   template <typename X0, typename X1, typename X2, typename W,
             typename std::enable_if<IsDataContainer<X0>::value && IsDataContainer<X1>::value && IsDataContainer<X2>::value &&
                                        !IsDataContainer<W>::value,
                                     int>::type = 0>
   void Exec(unsigned int slot, const X0 &x0s, const X1 &x1s, const X2 &x2s, const W w)
   {
      auto thisSlotH = fObjects[slot];
      if (!(x0s.size() == x1s.size() && x1s.size() == x2s.size())) {
         throw std::runtime_error("Cannot fill histogram with values in containers of different sizes.");
      }
      auto x0sIt = std::begin(x0s);
      const auto x0sEnd = std::end(x0s);
      auto x1sIt = std::begin(x1s);
      auto x2sIt = std::begin(x2s);
      for (; x0sIt != x0sEnd; x0sIt++, x1sIt++, x2sIt++) {
         thisSlotH->Fill(*x0sIt, *x1sIt, *x2sIt, w);
      }
   }

   void Initialize() { /* noop */}

   void Finalize()
   {
      auto resObj = fObjects[0];
      const auto nSlots = fObjects.size();
      TList l;
      l.SetOwner(); // The list will free the memory associated to its elements upon destruction
      for (unsigned int slot = 1; slot < nSlots; ++slot) {
         l.Add(fObjects[slot]);
      }

      resObj->Merge(&l);
   }

   HIST &PartialUpdate(unsigned int slot) { return *fObjects[slot]; }

   // Helper functions for RMergeableValue
   std::unique_ptr<RMergeableValueBase> GetMergeableValue() const final
   {
      return std::make_unique<RMergeableFill<HIST>>(*fObjects[0]);
   }

   std::string GetActionName() { return "FillPar"; }
};

class FillTGraphHelper : public ROOT::Detail::RDF::RActionImpl<FillTGraphHelper> {
public:
   using Result_t = ::TGraph;

private:
   std::vector<::TGraph *> fGraphs;

public:
   FillTGraphHelper(FillTGraphHelper &&) = default;
   FillTGraphHelper(const FillTGraphHelper &) = delete;

   // The last parameter is always false, as at the moment there is no way to propagate the parameter from the user to
   // this method
   FillTGraphHelper(const std::shared_ptr<::TGraph> &g, const unsigned int nSlots) : fGraphs(nSlots, nullptr)
   {
      fGraphs[0] = g.get();
      // Initialise all other slots
      for (unsigned int i = 1; i < nSlots; ++i) {
         fGraphs[i] = new TGraph(*fGraphs[0]);
      }
   }

   void Initialize() {}
   void InitTask(TTreeReader *, unsigned int) {}

   template <typename X0, typename X1,
             typename std::enable_if<IsDataContainer<X0>::value && IsDataContainer<X1>::value, int>::type = 0>
   void Exec(unsigned int slot, const X0 &x0s, const X1 &x1s)
   {
      if (x0s.size() != x1s.size()) {
         throw std::runtime_error("Cannot fill Graph with values in containers of different sizes.");
      }
      auto thisSlotG = fGraphs[slot];
      auto x0sIt = std::begin(x0s);
      const auto x0sEnd = std::end(x0s);
      auto x1sIt = std::begin(x1s);
      for (; x0sIt != x0sEnd; x0sIt++, x1sIt++) {
         thisSlotG->SetPoint(thisSlotG->GetN(), *x0sIt, *x1sIt);
      }
   }

   template <typename X0, typename X1>
   void Exec(unsigned int slot, X0 x0, X1 x1)
   {
      auto thisSlotG = fGraphs[slot];
      thisSlotG->SetPoint(thisSlotG->GetN(), x0, x1);
   }

   void Finalize()
   {
      const auto nSlots = fGraphs.size();
      auto resGraph = fGraphs[0];
      TList l;
      l.SetOwner(); // The list will free the memory associated to its elements upon destruction
      for (unsigned int slot = 1; slot < nSlots; ++slot) {
         l.Add(fGraphs[slot]);
      }
      resGraph->Merge(&l);
   }

   // Helper functions for RMergeableValue
   std::unique_ptr<RMergeableValueBase> GetMergeableValue() const final
   {
      return std::make_unique<RMergeableFill<Result_t>>(*fGraphs[0]);
   }

   std::string GetActionName() { return "Graph"; }

   Result_t &PartialUpdate(unsigned int slot) { return *fGraphs[slot]; }
};

// In case of the take helper we have 4 cases:
// 1. The column is not an RVec, the collection is not a vector
// 2. The column is not an RVec, the collection is a vector
// 3. The column is an RVec, the collection is not a vector
// 4. The column is an RVec, the collection is a vector

template <typename V, typename COLL>
void FillColl(V&& v, COLL& c) {
   c.emplace_back(v);
}

// Use push_back for bool since some compilers do not support emplace_back.
template <typename COLL>
void FillColl(bool v, COLL& c) {
   c.push_back(v);
}

// Case 1.: The column is not an RVec, the collection is not a vector
// No optimisations, no transformations: just copies.
template <typename RealT_t, typename T, typename COLL>
class TakeHelper : public RActionImpl<TakeHelper<RealT_t, T, COLL>> {
   Results<std::shared_ptr<COLL>> fColls;

public:
   using ColumnTypes_t = TypeList<T>;
   TakeHelper(const std::shared_ptr<COLL> &resultColl, const unsigned int nSlots)
   {
      fColls.emplace_back(resultColl);
      for (unsigned int i = 1; i < nSlots; ++i)
         fColls.emplace_back(std::make_shared<COLL>());
   }
   TakeHelper(TakeHelper &&);
   TakeHelper(const TakeHelper &) = delete;

   void InitTask(TTreeReader *, unsigned int) {}

   void Exec(unsigned int slot, T &v) { FillColl(v, *fColls[slot]); }

   void Initialize() { /* noop */}

   void Finalize()
   {
      auto rColl = fColls[0];
      for (unsigned int i = 1; i < fColls.size(); ++i) {
         const auto &coll = fColls[i];
         const auto end = coll->end();
         // Use an explicit loop here to prevent compiler warnings introduced by
         // clang's range-based loop analysis and vector<bool> references.
         for (auto j = coll->begin(); j != end; j++) {
            FillColl(*j, *rColl);
         }
      }
   }

   COLL &PartialUpdate(unsigned int slot) { return *fColls[slot].get(); }

   std::string GetActionName() { return "Take"; }
};

// Case 2.: The column is not an RVec, the collection is a vector
// Optimisations, no transformations: just copies.
template <typename RealT_t, typename T>
class TakeHelper<RealT_t, T, std::vector<T>> : public RActionImpl<TakeHelper<RealT_t, T, std::vector<T>>> {
   Results<std::shared_ptr<std::vector<T>>> fColls;

public:
   using ColumnTypes_t = TypeList<T>;
   TakeHelper(const std::shared_ptr<std::vector<T>> &resultColl, const unsigned int nSlots)
   {
      fColls.emplace_back(resultColl);
      for (unsigned int i = 1; i < nSlots; ++i) {
         auto v = std::make_shared<std::vector<T>>();
         v->reserve(1024);
         fColls.emplace_back(v);
      }
   }
   TakeHelper(TakeHelper &&);
   TakeHelper(const TakeHelper &) = delete;

   void InitTask(TTreeReader *, unsigned int) {}

   void Exec(unsigned int slot, T &v) { FillColl(v, *fColls[slot]); }

   void Initialize() { /* noop */}

   // This is optimised to treat vectors
   void Finalize()
   {
      ULong64_t totSize = 0;
      for (auto &coll : fColls)
         totSize += coll->size();
      auto rColl = fColls[0];
      rColl->reserve(totSize);
      for (unsigned int i = 1; i < fColls.size(); ++i) {
         auto &coll = fColls[i];
         rColl->insert(rColl->end(), coll->begin(), coll->end());
      }
   }

   std::vector<T> &PartialUpdate(unsigned int slot) { return *fColls[slot]; }

   std::string GetActionName() { return "Take"; }
};

// Case 3.: The column is a RVec, the collection is not a vector
// No optimisations, transformations from RVecs to vectors
template <typename RealT_t, typename COLL>
class TakeHelper<RealT_t, RVec<RealT_t>, COLL> : public RActionImpl<TakeHelper<RealT_t, RVec<RealT_t>, COLL>> {
   Results<std::shared_ptr<COLL>> fColls;

public:
   using ColumnTypes_t = TypeList<RVec<RealT_t>>;
   TakeHelper(const std::shared_ptr<COLL> &resultColl, const unsigned int nSlots)
   {
      fColls.emplace_back(resultColl);
      for (unsigned int i = 1; i < nSlots; ++i)
         fColls.emplace_back(std::make_shared<COLL>());
   }
   TakeHelper(TakeHelper &&);
   TakeHelper(const TakeHelper &) = delete;

   void InitTask(TTreeReader *, unsigned int) {}

   void Exec(unsigned int slot, RVec<RealT_t> av) { fColls[slot]->emplace_back(av.begin(), av.end()); }

   void Initialize() { /* noop */}

   void Finalize()
   {
      auto rColl = fColls[0];
      for (unsigned int i = 1; i < fColls.size(); ++i) {
         auto &coll = fColls[i];
         for (auto &v : *coll) {
            rColl->emplace_back(v);
         }
      }
   }

   std::string GetActionName() { return "Take"; }
};

// Case 4.: The column is an RVec, the collection is a vector
// Optimisations, transformations from RVecs to vectors
template <typename RealT_t>
class TakeHelper<RealT_t, RVec<RealT_t>, std::vector<RealT_t>>
   : public RActionImpl<TakeHelper<RealT_t, RVec<RealT_t>, std::vector<RealT_t>>> {

   Results<std::shared_ptr<std::vector<std::vector<RealT_t>>>> fColls;

public:
   using ColumnTypes_t = TypeList<RVec<RealT_t>>;
   TakeHelper(const std::shared_ptr<std::vector<std::vector<RealT_t>>> &resultColl, const unsigned int nSlots)
   {
      fColls.emplace_back(resultColl);
      for (unsigned int i = 1; i < nSlots; ++i) {
         auto v = std::make_shared<std::vector<RealT_t>>();
         v->reserve(1024);
         fColls.emplace_back(v);
      }
   }
   TakeHelper(TakeHelper &&);
   TakeHelper(const TakeHelper &) = delete;

   void InitTask(TTreeReader *, unsigned int) {}

   void Exec(unsigned int slot, RVec<RealT_t> av) { fColls[slot]->emplace_back(av.begin(), av.end()); }

   void Initialize() { /* noop */}

   // This is optimised to treat vectors
   void Finalize()
   {
      ULong64_t totSize = 0;
      for (auto &coll : fColls)
         totSize += coll->size();
      auto rColl = fColls[0];
      rColl->reserve(totSize);
      for (unsigned int i = 1; i < fColls.size(); ++i) {
         auto &coll = fColls[i];
         rColl->insert(rColl->end(), coll->begin(), coll->end());
      }
   }

   std::string GetActionName() { return "Take"; }
};

// Extern templates for TakeHelper
// NOTE: The move-constructor of specializations declared as extern templates
// must be defined out of line, otherwise cling fails to find its symbol.
template <typename RealT_t, typename T, typename COLL>
TakeHelper<RealT_t, T, COLL>::TakeHelper(TakeHelper<RealT_t, T, COLL> &&) = default;
template <typename RealT_t, typename T>
TakeHelper<RealT_t, T, std::vector<T>>::TakeHelper(TakeHelper<RealT_t, T, std::vector<T>> &&) = default;
template <typename RealT_t, typename COLL>
TakeHelper<RealT_t, RVec<RealT_t>, COLL>::TakeHelper(TakeHelper<RealT_t, RVec<RealT_t>, COLL> &&) = default;
template <typename RealT_t>
TakeHelper<RealT_t, RVec<RealT_t>, std::vector<RealT_t>>::TakeHelper(TakeHelper<RealT_t, RVec<RealT_t>, std::vector<RealT_t>> &&) = default;

// External templates are disabled for gcc5 since this version wrongly omits the C++11 ABI attribute
#if __GNUC__ > 5
extern template class TakeHelper<bool, bool, std::vector<bool>>;
extern template class TakeHelper<unsigned int, unsigned int, std::vector<unsigned int>>;
extern template class TakeHelper<unsigned long, unsigned long, std::vector<unsigned long>>;
extern template class TakeHelper<unsigned long long, unsigned long long, std::vector<unsigned long long>>;
extern template class TakeHelper<int, int, std::vector<int>>;
extern template class TakeHelper<long, long, std::vector<long>>;
extern template class TakeHelper<long long, long long, std::vector<long long>>;
extern template class TakeHelper<float, float, std::vector<float>>;
extern template class TakeHelper<double, double, std::vector<double>>;
#endif


template <typename ResultType>
class MinHelper : public RActionImpl<MinHelper<ResultType>> {
   const std::shared_ptr<ResultType> fResultMin;
   Results<ResultType> fMins;

public:
   MinHelper(MinHelper &&) = default;
   MinHelper(const std::shared_ptr<ResultType> &minVPtr, const unsigned int nSlots)
      : fResultMin(minVPtr), fMins(nSlots, std::numeric_limits<ResultType>::max())
   {
   }

   void Exec(unsigned int slot, ResultType v) { fMins[slot] = std::min(v, fMins[slot]); }

   void InitTask(TTreeReader *, unsigned int) {}

   template <typename T, typename std::enable_if<IsDataContainer<T>::value, int>::type = 0>
   void Exec(unsigned int slot, const T &vs)
   {
      for (auto &&v : vs)
         fMins[slot] = std::min(static_cast<ResultType>(v), fMins[slot]);
   }

   void Initialize() { /* noop */}

   void Finalize()
   {
      *fResultMin = std::numeric_limits<ResultType>::max();
      for (auto &m : fMins)
         *fResultMin = std::min(m, *fResultMin);
   }

   // Helper functions for RMergeableValue
   std::unique_ptr<RMergeableValueBase> GetMergeableValue() const final
   {
      return std::make_unique<RMergeableMin<ResultType>>(*fResultMin);
   }

   ResultType &PartialUpdate(unsigned int slot) { return fMins[slot]; }

   std::string GetActionName() { return "Min"; }
};

// TODO
// extern template void MinHelper::Exec(unsigned int, const std::vector<float> &);
// extern template void MinHelper::Exec(unsigned int, const std::vector<double> &);
// extern template void MinHelper::Exec(unsigned int, const std::vector<char> &);
// extern template void MinHelper::Exec(unsigned int, const std::vector<int> &);
// extern template void MinHelper::Exec(unsigned int, const std::vector<unsigned int> &);

template <typename ResultType>
class MaxHelper : public RActionImpl<MaxHelper<ResultType>> {
   const std::shared_ptr<ResultType> fResultMax;
   Results<ResultType> fMaxs;

public:
   MaxHelper(MaxHelper &&) = default;
   MaxHelper(const MaxHelper &) = delete;
   MaxHelper(const std::shared_ptr<ResultType> &maxVPtr, const unsigned int nSlots)
      : fResultMax(maxVPtr), fMaxs(nSlots, std::numeric_limits<ResultType>::lowest())
   {
   }

   void InitTask(TTreeReader *, unsigned int) {}
   void Exec(unsigned int slot, ResultType v) { fMaxs[slot] = std::max(v, fMaxs[slot]); }

   template <typename T, typename std::enable_if<IsDataContainer<T>::value, int>::type = 0>
   void Exec(unsigned int slot, const T &vs)
   {
      for (auto &&v : vs)
         fMaxs[slot] = std::max(static_cast<ResultType>(v), fMaxs[slot]);
   }

   void Initialize() { /* noop */}

   void Finalize()
   {
      *fResultMax = std::numeric_limits<ResultType>::lowest();
      for (auto &m : fMaxs) {
         *fResultMax = std::max(m, *fResultMax);
      }
   }

   // Helper functions for RMergeableValue
   std::unique_ptr<RMergeableValueBase> GetMergeableValue() const final
   {
      return std::make_unique<RMergeableMax<ResultType>>(*fResultMax);
   }

   ResultType &PartialUpdate(unsigned int slot) { return fMaxs[slot]; }

   std::string GetActionName() { return "Max"; }
};

// TODO
// extern template void MaxHelper::Exec(unsigned int, const std::vector<float> &);
// extern template void MaxHelper::Exec(unsigned int, const std::vector<double> &);
// extern template void MaxHelper::Exec(unsigned int, const std::vector<char> &);
// extern template void MaxHelper::Exec(unsigned int, const std::vector<int> &);
// extern template void MaxHelper::Exec(unsigned int, const std::vector<unsigned int> &);

template <typename ResultType>
class SumHelper : public RActionImpl<SumHelper<ResultType>> {
   const std::shared_ptr<ResultType> fResultSum;
   Results<ResultType> fSums;

   /// Evaluate neutral element for this type and the sum operation.
   /// This is assumed to be any_value - any_value if operator- is defined
   /// for the type, otherwise a default-constructed ResultType{} is used.
   template <typename T = ResultType>
   auto NeutralElement(const T &v, int /*overloadresolver*/) -> decltype(v - v)
   {
      return v - v;
   }

   template <typename T = ResultType, typename Dummy = int>
   ResultType NeutralElement(const T &, Dummy) // this overload has lower priority thanks to the template arg
   {
      return ResultType{};
   }

public:
   SumHelper(SumHelper &&) = default;
   SumHelper(const SumHelper &) = delete;
   SumHelper(const std::shared_ptr<ResultType> &sumVPtr, const unsigned int nSlots)
      : fResultSum(sumVPtr), fSums(nSlots, NeutralElement(*sumVPtr, -1))
   {
   }

   void InitTask(TTreeReader *, unsigned int) {}
   void Exec(unsigned int slot, ResultType v) { fSums[slot] += v; }

   template <typename T, typename std::enable_if<IsDataContainer<T>::value, int>::type = 0>
   void Exec(unsigned int slot, const T &vs)
   {
      for (auto &&v : vs)
         fSums[slot] += static_cast<ResultType>(v);
   }

   void Initialize() { /* noop */}

   void Finalize()
   {
      for (auto &m : fSums)
         *fResultSum += m;
   }

   // Helper functions for RMergeableValue
   std::unique_ptr<RMergeableValueBase> GetMergeableValue() const final
   {
      return std::make_unique<RMergeableSum<ResultType>>(*fResultSum);
   }

   ResultType &PartialUpdate(unsigned int slot) { return fSums[slot]; }

   std::string GetActionName() { return "Sum"; }
};

class MeanHelper : public RActionImpl<MeanHelper> {
   const std::shared_ptr<double> fResultMean;
   std::vector<ULong64_t> fCounts;
   std::vector<double> fSums;
   std::vector<double> fPartialMeans;

public:
   MeanHelper(const std::shared_ptr<double> &meanVPtr, const unsigned int nSlots);
   MeanHelper(MeanHelper &&) = default;
   MeanHelper(const MeanHelper &) = delete;
   void InitTask(TTreeReader *, unsigned int) {}
   void Exec(unsigned int slot, double v);

   template <typename T, typename std::enable_if<IsDataContainer<T>::value, int>::type = 0>
   void Exec(unsigned int slot, const T &vs)
   {
      for (auto &&v : vs) {
         fSums[slot] += v;
         fCounts[slot]++;
      }
   }

   void Initialize() { /* noop */}

   void Finalize();

   // Helper functions for RMergeableValue
   std::unique_ptr<RMergeableValueBase> GetMergeableValue() const final
   {
      const ULong64_t counts = std::accumulate(fCounts.begin(), fCounts.end(), 0ull);
      return std::make_unique<RMergeableMean>(*fResultMean, counts);
   }

   double &PartialUpdate(unsigned int slot);

   std::string GetActionName() { return "Mean"; }
};

extern template void MeanHelper::Exec(unsigned int, const std::vector<float> &);
extern template void MeanHelper::Exec(unsigned int, const std::vector<double> &);
extern template void MeanHelper::Exec(unsigned int, const std::vector<char> &);
extern template void MeanHelper::Exec(unsigned int, const std::vector<int> &);
extern template void MeanHelper::Exec(unsigned int, const std::vector<unsigned int> &);

class StdDevHelper : public RActionImpl<StdDevHelper> {
   // Number of subsets of data
   const unsigned int fNSlots;
   const std::shared_ptr<double> fResultStdDev;
   // Number of element for each slot
   std::vector<ULong64_t> fCounts;
   // Mean of each slot
   std::vector<double> fMeans;
   // Squared distance from the mean
   std::vector<double> fDistancesfromMean;

public:
   StdDevHelper(const std::shared_ptr<double> &meanVPtr, const unsigned int nSlots);
   StdDevHelper(StdDevHelper &&) = default;
   StdDevHelper(const StdDevHelper &) = delete;
   void InitTask(TTreeReader *, unsigned int) {}
   void Exec(unsigned int slot, double v);

   template <typename T, typename std::enable_if<IsDataContainer<T>::value, int>::type = 0>
   void Exec(unsigned int slot, const T &vs)
   {
      for (auto &&v : vs) {
         Exec(slot, v);
      }
   }

   void Initialize() { /* noop */}

   void Finalize();

   // Helper functions for RMergeableValue
   std::unique_ptr<RMergeableValueBase> GetMergeableValue() const final
   {
      const ULong64_t counts = std::accumulate(fCounts.begin(), fCounts.end(), 0ull);
      const Double_t mean =
         std::inner_product(fMeans.begin(), fMeans.end(), fCounts.begin(), 0.) / static_cast<Double_t>(counts);
      return std::make_unique<RMergeableStdDev>(*fResultStdDev, counts, mean);
   }

   std::string GetActionName() { return "StdDev"; }
};

extern template void StdDevHelper::Exec(unsigned int, const std::vector<float> &);
extern template void StdDevHelper::Exec(unsigned int, const std::vector<double> &);
extern template void StdDevHelper::Exec(unsigned int, const std::vector<char> &);
extern template void StdDevHelper::Exec(unsigned int, const std::vector<int> &);
extern template void StdDevHelper::Exec(unsigned int, const std::vector<unsigned int> &);

template <typename PrevNodeType>
class DisplayHelper : public RActionImpl<DisplayHelper<PrevNodeType>> {
private:
   using Display_t = ROOT::RDF::RDisplay;
   const std::shared_ptr<Display_t> fDisplayerHelper;
   const std::shared_ptr<PrevNodeType> fPrevNode;

public:
   DisplayHelper(const std::shared_ptr<Display_t> &d, const std::shared_ptr<PrevNodeType> &prevNode)
      : fDisplayerHelper(d), fPrevNode(prevNode)
   {
   }
   DisplayHelper(DisplayHelper &&) = default;
   DisplayHelper(const DisplayHelper &) = delete;
   void InitTask(TTreeReader *, unsigned int) {}

   template <typename... Columns>
   void Exec(unsigned int, Columns &... columns)
   {
      fDisplayerHelper->AddRow(columns...);
      if (!fDisplayerHelper->HasNext()) {
         fPrevNode->StopProcessing();
      }
   }

   void Initialize() {}

   void Finalize() {}

   std::string GetActionName() { return "Display"; }
};

// std::vector<bool> is special, and not in a good way. As a consequence Snapshot of RVec<bool> needs to be treated
// specially. In particular, if RVec<bool> is filled with a (fixed or variable size) boolean array coming from
// a ROOT file, when writing out the correspinding branch from a Snapshot we do not have an address to set for the
// TTree branch (std::vector<bool> and, consequently, RVec<bool> do not provide a `data()` method).
// Bools is a lightweight wrapper around a C array of booleans that is meant to provide a stable address for the
// output TTree to read the contents of the snapshotted branches at Fill time.
class BoolArray {
   std::size_t fSize = 0;
   bool *fBools = nullptr;

   bool *CopyVector(const RVec<bool> &v)
   {
      auto b = new bool[fSize];
      std::copy(v.begin(), v.end(), b);
      return b;
   }

   bool *CopyArray(bool *o, std::size_t size)
   {
      auto b = new bool[size];
      for (auto i = 0u; i < size; ++i)
         b[i] = o[i];
      return b;
   }

public:
   // this generic constructor could be replaced with a constexpr if in SetBranchesHelper
   BoolArray() = default;
   template <typename T>
   BoolArray(const T &) { throw std::runtime_error("This constructor should never be called"); }
   BoolArray(const RVec<bool> &v) : fSize(v.size()), fBools(CopyVector(v)) {}
   BoolArray(const BoolArray &b)
   {
      CopyArray(b.fBools, b.fSize);
   }
   BoolArray &operator=(const BoolArray &b)
   {
      delete[] fBools;
      CopyArray(b.fBools, b.fSize);
      return *this;
   }
   BoolArray(BoolArray &&b)
   {
      fSize = b.fSize;
      fBools = b.fBools;
      b.fSize = 0;
      b.fBools = nullptr;
   }
   BoolArray &operator=(BoolArray &&b)
   {
      delete[] fBools;
      fSize = b.fSize;
      fBools = b.fBools;
      b.fSize = 0;
      b.fBools = nullptr;
      return *this;
   }
   ~BoolArray() { delete[] fBools; }
   std::size_t Size() const { return fSize; }
   bool *Data() { return fBools; }
};
using BoolArrayMap = std::map<std::string, BoolArray>;

inline bool *UpdateBoolArrayIfBool(BoolArrayMap &boolArrays, RVec<bool> &v, const std::string &outName)
{
   // create a boolArrays entry
   boolArrays[outName] = BoolArray(v);
   return boolArrays[outName].Data();
}

template <typename T>
T *UpdateBoolArrayIfBool(BoolArrayMap &, RVec<T> &v, const std::string &)
{
   return v.data();
}

// Helper which gets the return value of the data() method if the type is an
// RVec (of anything but a bool), nullptr otherwise.
inline void *GetData(ROOT::VecOps::RVec<bool> & /*v*/)
{
   return nullptr;
}

template <typename T>
void *GetData(ROOT::VecOps::RVec<T> &v)
{
   return v.data();
}

template <typename T>
void *GetData(T & /*v*/)
{
   return nullptr;
}

template <typename T>
void SetBranchesHelper(BoolArrayMap &, TTree *inputTree, TTree &outputTree, const std::string &inName,
                       const std::string &name, TBranch *&branch, void *&branchAddress, T *address)
{
   auto *inputBranch = inputTree ? inputTree->GetBranch(inName.c_str()) : nullptr;
   if (inputBranch) {
      // Respect the original bufsize and splitlevel arguments
      // In particular, by keeping splitlevel equal to 0 if this was the case for `inputBranch`, we avoid
      // writing garbage when unsplit objects cannot be written as split objects (e.g. in case of a polymorphic
      // TObject branch, see https://bit.ly/2EjLMId ).
      const auto bufSize = inputBranch->GetBasketSize();
      const auto splitLevel = inputBranch->GetSplitLevel();

      static TClassRef tbo_cl("TBranchObject");
      if (inputBranch->IsA() == tbo_cl) {
         // Need to pass a pointer to pointer
         outputTree.Branch(name.c_str(), (T **)inputBranch->GetAddress(), bufSize, splitLevel);
      } else {
         outputTree.Branch(name.c_str(), address, bufSize, splitLevel);
      }
   } else {
      outputTree.Branch(name.c_str(), address);
   }
   // This is not an array branch, so we don't need to register the address of the input branch.
   branch = nullptr;
   branchAddress = nullptr;
}

/// Helper function for SnapshotHelper and SnapshotHelperMT. It creates new branches for the output TTree of a Snapshot.
/// This overload is called for columns of type `RVec<T>`. For RDF, these can represent:
/// 1. c-style arrays in ROOT files, so we are sure that there are input trees to which we can ask the correct branch title
/// 2. RVecs coming from a custom column or a source
/// 3. vectors coming from ROOT files
/// 4. TClonesArray
///
/// In case of 1., we keep aside the pointer to the branch and the pointer to the input value (in `branch` and
/// `branchAddress`) so we can intercept changes in the address of the input branch and tell the output branch.
template <typename T>
void SetBranchesHelper(BoolArrayMap &boolArrays, TTree *inputTree, TTree &outputTree, const std::string &inName,
                       const std::string &outName, TBranch *&branch, void *&branchAddress, RVec<T> *ab)
{
   TBranch *inputBranch = nullptr;
   if (inputTree) {
      inputBranch = inputTree->GetBranch(inName.c_str());
      if (!inputBranch) {
         // try harder
         inputBranch = inputTree->FindBranch(inName.c_str());
      }
   }
   const bool isTClonesArray = inputBranch != nullptr && std::string(inputBranch->GetClassName()) == "TClonesArray";
   const auto mustWriteStdVec = !inputBranch || isTClonesArray ||
                                ROOT::ESTLType::kSTLvector == TClassEdit::IsSTLCont(inputBranch->GetClassName());

   if (mustWriteStdVec) {
      // Treat:
      // 2. RVec coming from a custom column or a source
      // 3. RVec coming from a column on disk of type vector (the RVec is adopting the data of that vector)
      // 4. TClonesArray.
      // In all cases, we write out a std::vector<T> when the column is RVec<T>
      if (isTClonesArray) {
         Warning("Snapshot",
                 "Branch \"%s\" contains TClonesArrays but the type specified to Snapshot was RVec<T>. The branch will "
                 "be written out as a std::vector instead of a TClonesArray. Specify that the type of the branch is "
                 "TClonesArray as a Snapshot template parameter to write out a TClonesArray instead.", inName.c_str());
      }
      outputTree.Branch(outName.c_str(), &ab->AsVector());
      return;
   }

   // Treat 1, the C-array case
   auto *const leaf = static_cast<TLeaf *>(inputBranch->GetListOfLeaves()->UncheckedAt(0));
   const auto bname = leaf->GetName();
   const auto counterStr =
      leaf->GetLeafCount() ? std::string(leaf->GetLeafCount()->GetName()) : std::to_string(leaf->GetLenStatic());
   const auto btype = leaf->GetTypeName();
   const auto rootbtype = TypeName2ROOTTypeName(btype);
   const auto leaflist = std::string(bname) + "[" + counterStr + "]/" + rootbtype;

   /// RVec<bool> is special because std::vector<bool> is special. In particular, it has no `data()`,
   /// so we need to explicitly manage storage of the data that the tree needs to Fill branches with.
   auto dataPtr = UpdateBoolArrayIfBool(boolArrays, *ab, outName);

   auto *const outputBranch = outputTree.Branch(outName.c_str(), dataPtr, leaflist.c_str());
   outputBranch->SetTitle(inputBranch->GetTitle());

   // Record the branch ptr and the address associated to it if this is not a bool array
   if (!std::is_same<bool, T>::value) {
      branch = outputBranch;
      branchAddress = GetData(*ab);
   }
}

// generic version, no-op
template <typename T>
void UpdateBoolArray(BoolArrayMap &, T&, const std::string &, TTree &) {}

// RVec<bool> overload, update boolArrays if needed
inline void UpdateBoolArray(BoolArrayMap &boolArrays, RVec<bool> &v, const std::string &outName, TTree &t)
{
   // in case the RVec<bool> does not correspond to a bool C-array
   if (boolArrays.find(outName) == boolArrays.end())
      return;

   if (v.size() > boolArrays[outName].Size()) {
      boolArrays[outName] = BoolArray(v); // resize and copy
      t.SetBranchAddress(outName.c_str(), boolArrays[outName].Data());
   }
   else {
      std::copy(v.begin(), v.end(), boolArrays[outName].Data()); // just copy
   }
}

void ValidateSnapshotOutput(const RSnapshotOptions &opts, const std::string &treeName, const std::string &fileName);

/// Helper object for a single-thread Snapshot action
template <typename... ColTypes>
class SnapshotHelper : public RActionImpl<SnapshotHelper<ColTypes...>> {
   const std::string fFileName;
   const std::string fDirName;
   const std::string fTreeName;
   const RSnapshotOptions fOptions;
   std::unique_ptr<TFile> fOutputFile;
   std::unique_ptr<TTree> fOutputTree; // must be a ptr because TTrees are not copy/move constructible
   bool fIsFirstEvent{true};
   const ColumnNames_t fInputBranchNames; // This contains the resolved aliases
   const ColumnNames_t fOutputBranchNames;
   TTree *fInputTree = nullptr; // Current input tree. Set at initialization time (`InitTask`)
   BoolArrayMap fBoolArrays; // Storage for C arrays of bools to be written out
   std::vector<TBranch *> fBranches;     // Addresses of branches in output, non-null only for the ones holding C arrays
   std::vector<void *> fBranchAddresses; // Addresses associated to output branches, non-null only for the ones holding C arrays

public:
   using ColumnTypes_t = TypeList<ColTypes...>;
   SnapshotHelper(std::string_view filename, std::string_view dirname, std::string_view treename,
                  const ColumnNames_t &vbnames, const ColumnNames_t &bnames, const RSnapshotOptions &options)
      : fFileName(filename), fDirName(dirname), fTreeName(treename), fOptions(options), fInputBranchNames(vbnames),
        fOutputBranchNames(ReplaceDotWithUnderscore(bnames)), fBranches(vbnames.size(), nullptr),
        fBranchAddresses(vbnames.size(), nullptr)
   {
      ValidateSnapshotOutput(fOptions, fTreeName, fFileName);
   }

   SnapshotHelper(const SnapshotHelper &) = delete;
   SnapshotHelper(SnapshotHelper &&) = default;

   void InitTask(TTreeReader *r, unsigned int /* slot */)
   {
      if (r)
         fInputTree = r->GetTree();
   }

   void Exec(unsigned int /* slot */, ColTypes &... values)
   {
      using ind_t = std::index_sequence_for<ColTypes...>;
      if (! fIsFirstEvent) {
         UpdateCArraysPtrs(values..., ind_t{});
      } else {
         SetBranches(values..., ind_t{});
         // AddClone guarantees that if the input file changes the branches of the output tree are updated with the new
         // addresses of the branch values
         if (fInputTree != nullptr)
            fInputTree->AddClone(fOutputTree.get());
         fIsFirstEvent = false;
      }
      UpdateBoolArrays(values..., ind_t{});
      fOutputTree->Fill();
   }

   template <std::size_t... S>
   void UpdateCArraysPtrs(ColTypes &... values, std::index_sequence<S...> /*dummy*/)
   {
      // This code deals with branches which hold C arrays of variable size. It can happen that the buffers
      // associated to those is re-allocated. As a result the value of the pointer can change therewith
      // leaving associated to the branch of the output tree an invalid pointer.
      // With this code, we set the value of the pointer in the output branch anew when needed.
      // Nota bene: the extra ",0" after the invocation of SetAddress, is because that method returns void and 
      // we need an int for the expander list.
      int expander[] = {(fBranches[S] && fBranchAddresses[S] != GetData(values)
                         ? fBranches[S]->SetAddress(GetData(values)),
                         fBranchAddresses[S] = GetData(values), 0 : 0, 0)...,
                        0};
      (void)expander; // avoid unused variable warnings for older compilers such as gcc 4.9
   }

   template <std::size_t... S>
   void SetBranches(ColTypes &... values, std::index_sequence<S...> /*dummy*/)
   {
      // create branches in output tree (and fill fBoolArrays for RVec<bool> columns)
      int expander[] = {(SetBranchesHelper(fBoolArrays, fInputTree, *fOutputTree, fInputBranchNames[S],
                                           fOutputBranchNames[S], fBranches[S], fBranchAddresses[S], &values),
                         0)...,
                        0};
      (void)expander; // avoid unused variable warnings for older compilers such as gcc 4.9
   }

   template <std::size_t... S>
   void UpdateBoolArrays(ColTypes &...values, std::index_sequence<S...> /*dummy*/)
   {
      int expander[] = {(UpdateBoolArray(fBoolArrays, values, fOutputBranchNames[S], *fOutputTree), 0)..., 0};
      (void)expander; // avoid unused variable warnings for older compilers such as gcc 4.9
   }

   void Initialize()
   {
      fOutputFile.reset(
         TFile::Open(fFileName.c_str(), fOptions.fMode.c_str(), /*ftitle=*/"",
                     ROOT::CompressionSettings(fOptions.fCompressionAlgorithm, fOptions.fCompressionLevel)));
      if(!fOutputFile)
         throw std::runtime_error("Snapshot: could not create output file " + fFileName);

      TDirectory *outputDir = fOutputFile.get();
      if (!fDirName.empty()) {
         TString checkupdate = fOptions.fMode;
         checkupdate.ToLower();
         if (checkupdate == "update")
            outputDir = fOutputFile->mkdir(fDirName.c_str(), "", true);  // do not overwrite existing directory
         else
            outputDir = fOutputFile->mkdir(fDirName.c_str());
      }

      fOutputTree =
         std::make_unique<TTree>(fTreeName.c_str(), fTreeName.c_str(), fOptions.fSplitLevel, /*dir=*/outputDir);

      if (fOptions.fAutoFlush)
         fOutputTree->SetAutoFlush(fOptions.fAutoFlush);
   }

   void Finalize()
   {
      if (fOutputFile && fOutputTree) {
         // because TTree::Write writes in gDirectory, not in fDirectory
         fOutputTree->AutoSave("flushbaskets");
         // must destroy the TTree first, otherwise TFile will delete it too leading to a double delete
         fOutputTree.reset();
         fOutputFile->Close();
      } else {
         Warning("Snapshot", "A lazy Snapshot action was booked but never triggered.");
      }
   }

   std::string GetActionName() { return "Snapshot"; }
};

/// Helper object for a multi-thread Snapshot action
template <typename... ColTypes>
class SnapshotHelperMT : public RActionImpl<SnapshotHelperMT<ColTypes...>> {
   const unsigned int fNSlots;
   std::unique_ptr<ROOT::Experimental::TBufferMerger> fMerger; // must use a ptr because TBufferMerger is not movable
   std::vector<std::shared_ptr<ROOT::Experimental::TBufferMergerFile>> fOutputFiles;
   std::vector<std::unique_ptr<TTree>> fOutputTrees;
   std::vector<int> fIsFirstEvent;        // vector<bool> does not allow concurrent writing of different elements
   const std::string fFileName;           // name of the output file name
   const std::string fDirName;            // name of TFile subdirectory in which output must be written (possibly empty)
   const std::string fTreeName;           // name of output tree
   const RSnapshotOptions fOptions;       // struct holding options to pass down to TFile and TTree in this action
   const ColumnNames_t fInputBranchNames; // This contains the resolved aliases
   const ColumnNames_t fOutputBranchNames;
   std::vector<TTree *> fInputTrees; // Current input trees. Set at initialization time (`InitTask`)
   std::vector<BoolArrayMap> fBoolArrays; // Per-thread storage for C arrays of bools to be written out
   // Addresses of branches in output per slot, non-null only for the ones holding C arrays
   std::vector<std::vector<TBranch *>> fBranches;
   // Addresses associated to output branches per slot, non-null only for the ones holding C arrays
   std::vector<std::vector<void *>> fBranchAddresses; 

public:
   using ColumnTypes_t = TypeList<ColTypes...>;
   SnapshotHelperMT(const unsigned int nSlots, std::string_view filename, std::string_view dirname,
                    std::string_view treename, const ColumnNames_t &vbnames, const ColumnNames_t &bnames,
                    const RSnapshotOptions &options)
      : fNSlots(nSlots), fOutputFiles(fNSlots), fOutputTrees(fNSlots), fIsFirstEvent(fNSlots, 1), fFileName(filename),
        fDirName(dirname), fTreeName(treename), fOptions(options), fInputBranchNames(vbnames),
        fOutputBranchNames(ReplaceDotWithUnderscore(bnames)), fInputTrees(fNSlots), fBoolArrays(fNSlots),
        fBranches(fNSlots, std::vector<TBranch *>(vbnames.size(), nullptr)), 
        fBranchAddresses(fNSlots, std::vector<void *>(vbnames.size(), nullptr))
   {
      ValidateSnapshotOutput(fOptions, fTreeName, fFileName);
   }
   SnapshotHelperMT(const SnapshotHelperMT &) = delete;
   SnapshotHelperMT(SnapshotHelperMT &&) = default;

   void InitTask(TTreeReader *r, unsigned int slot)
   {
      ::TDirectory::TContext c; // do not let tasks change the thread-local gDirectory
      if (!fOutputFiles[slot]) {
         // first time this thread executes something, let's create a TBufferMerger output directory
         fOutputFiles[slot] = fMerger->GetFile();
      }
      TDirectory *treeDirectory = fOutputFiles[slot].get();
      if (!fDirName.empty()) {
         // call returnExistingDirectory=true since MT can end up making this call multiple times
         treeDirectory = fOutputFiles[slot]->mkdir(fDirName.c_str(), "", true);
      }
      // re-create output tree as we need to create its branches again, with new input variables
      // TODO we could instead create the output tree and its branches, change addresses of input variables in each task
      fOutputTrees[slot] =
         std::make_unique<TTree>(fTreeName.c_str(), fTreeName.c_str(), fOptions.fSplitLevel, /*dir=*/treeDirectory);
      fOutputTrees[slot]->SetBit(TTree::kEntriesReshuffled);
      // TODO can be removed when RDF supports interleaved TBB task execution properly, see ROOT-10269
      fOutputTrees[slot]->SetImplicitMT(false);
      if (fOptions.fAutoFlush)
         fOutputTrees[slot]->SetAutoFlush(fOptions.fAutoFlush);
      if (r) {
         // not an empty-source RDF
         fInputTrees[slot] = r->GetTree();
         // AddClone guarantees that if the input file changes the branches of the output tree are updated with the new
         // addresses of the branch values. We need this in case of friend trees with different cluster granularity
         // than the main tree.
         // FIXME: AddClone might result in many many (safe) warnings printed by TTree::CopyAddresses, see ROOT-9487.
         const auto friendsListPtr = fInputTrees[slot]->GetListOfFriends();
         if (friendsListPtr && friendsListPtr->GetEntries() > 0)
            fInputTrees[slot]->AddClone(fOutputTrees[slot].get());
      }
      fIsFirstEvent[slot] = 1; // reset first event flag for this slot
   }

   void FinalizeTask(unsigned int slot)
   {
      if (fOutputTrees[slot]->GetEntries() > 0)
         fOutputFiles[slot]->Write();
      // clear now to avoid concurrent destruction of output trees and input tree (which has them listed as fClones)
      fOutputTrees[slot].reset(nullptr);
   }

   void Exec(unsigned int slot, ColTypes &... values)
   {
      using ind_t = std::index_sequence_for<ColTypes...>;
      if (!fIsFirstEvent[slot]) {
         UpdateCArraysPtrs(slot, values..., ind_t{});
      } else {
         SetBranches(slot, values..., ind_t{});
         fIsFirstEvent[slot] = 0;
      }
      UpdateBoolArrays(slot, values..., ind_t{});
      fOutputTrees[slot]->Fill();
      auto entries = fOutputTrees[slot]->GetEntries();
      auto autoFlush = fOutputTrees[slot]->GetAutoFlush();
      if ((autoFlush > 0) && (entries % autoFlush == 0))
         fOutputFiles[slot]->Write();
   }

   template <std::size_t... S>
   void UpdateCArraysPtrs(unsigned int slot, ColTypes &... values, std::index_sequence<S...> /*dummy*/)
   {
      // This code deals with branches which hold C arrays of variable size. It can happen that the buffers
      // associated to those is re-allocated. As a result the value of the pointer can change therewith
      // leaving associated to the branch of the output tree an invalid pointer.
      // With this code, we set the value of the pointer in the output branch anew when needed.
      // Nota bene: the extra ",0" after the invocation of SetAddress, is because that method returns void and
      // we need an int for the expander list.
      (void)slot; // avoid bogus 'unused parameter' warning
      int expander[] = {(fBranches[slot][S] && fBranchAddresses[slot][S] != GetData(values)
                         ? fBranches[slot][S]->SetAddress(GetData(values)),
                         fBranchAddresses[slot][S] = GetData(values), 0 : 0, 0)...,
                        0};
      (void)expander; // avoid unused variable warnings for older compilers such as gcc 4.9
   }

   template <std::size_t... S>
   void SetBranches(unsigned int slot, ColTypes &... values, std::index_sequence<S...> /*dummy*/)
   {
         // hack to call TTree::Branch on all variadic template arguments
         int expander[] = {
            (SetBranchesHelper(fBoolArrays[slot], fInputTrees[slot], *fOutputTrees[slot], fInputBranchNames[S],
                               fOutputBranchNames[S], fBranches[slot][S], fBranchAddresses[slot][S], &values),
             0)...,
            0};
      (void)expander; // avoid unused variable warnings for older compilers such as gcc 4.9
      (void)slot;     // avoid unused variable warnings in gcc6.2
   }

   template <std::size_t... S>
   void UpdateBoolArrays(unsigned int slot, ColTypes &... values, std::index_sequence<S...> /*dummy*/)
   {
      (void)slot; // avoid bogus 'unused parameter' warning
      int expander[] = {
         (UpdateBoolArray(fBoolArrays[slot], values, fOutputBranchNames[S], *fOutputTrees[slot]), 0)..., 0};
      (void)expander; // avoid unused variable warnings for older compilers such as gcc 4.9
   }

   void Initialize()
   {
      const auto cs = ROOT::CompressionSettings(fOptions.fCompressionAlgorithm, fOptions.fCompressionLevel);
      auto out_file = TFile::Open(fFileName.c_str(), fOptions.fMode.c_str(), /*ftitle=*/fFileName.c_str(), cs);
      if(!out_file)
         throw std::runtime_error("Snapshot: could not create output file " + fFileName);
      fMerger = std::make_unique<ROOT::Experimental::TBufferMerger>(std::unique_ptr<TFile>(out_file));
   }

   void Finalize()
   {
      auto fileWritten = false;
      for (auto &file : fOutputFiles) {
         if (file) {
            file->Write();
            file->Close();
            fileWritten = true;
         }
      }

      if (!fileWritten) {
         if (std::none_of(fOutputFiles.begin(), fOutputFiles.end(), [] (const std::shared_ptr<ROOT::Experimental::TBufferMergerFile> ptr) { return bool(ptr); })) {
            Warning("Snapshot",
                    "No input entries (input TTree was empty or no entry passed the Filters). Output TTree is empty.");
         } else {
            Warning("Snapshot", "A lazy Snapshot action was booked but never triggered.");
         }
      }

      // flush all buffers to disk by destroying the TBufferMerger
      fOutputFiles.clear();
      fMerger.reset();
   }

   std::string GetActionName() { return "Snapshot"; }
};

template <typename Acc, typename Merge, typename R, typename T, typename U,
          bool MustCopyAssign = std::is_same<R, U>::value>
class AggregateHelper : public RActionImpl<AggregateHelper<Acc, Merge, R, T, U, MustCopyAssign>> {
   Acc fAggregate;
   Merge fMerge;
   const std::shared_ptr<U> fResult;
   Results<U> fAggregators;

public:
   using ColumnTypes_t = TypeList<T>;
   AggregateHelper(Acc &&f, Merge &&m, const std::shared_ptr<U> &result, const unsigned int nSlots)
      : fAggregate(std::move(f)), fMerge(std::move(m)), fResult(result), fAggregators(nSlots, *result)
   {
   }
   AggregateHelper(AggregateHelper &&) = default;
   AggregateHelper(const AggregateHelper &) = delete;

   void InitTask(TTreeReader *, unsigned int) {}

   template <bool MustCopyAssign_ = MustCopyAssign, typename std::enable_if<MustCopyAssign_, int>::type = 0>
   void Exec(unsigned int slot, const T &value)
   {
      fAggregators[slot] = fAggregate(fAggregators[slot], value);
   }

   template <bool MustCopyAssign_ = MustCopyAssign, typename std::enable_if<!MustCopyAssign_, int>::type = 0>
   void Exec(unsigned int slot, const T &value)
   {
      fAggregate(fAggregators[slot], value);
   }

   void Initialize() { /* noop */}

   template <typename MergeRet = typename CallableTraits<Merge>::ret_type,
             bool MergeAll = std::is_same<void, MergeRet>::value>
   typename std::enable_if<MergeAll, void>::type Finalize()
   {
      fMerge(fAggregators);
      *fResult = fAggregators[0];
   }

   template <typename MergeRet = typename CallableTraits<Merge>::ret_type,
             bool MergeTwoByTwo = std::is_same<U, MergeRet>::value>
   typename std::enable_if<MergeTwoByTwo, void>::type Finalize(...) // ... needed to let compiler distinguish overloads
   {
      for (const auto &acc : fAggregators)
         *fResult = fMerge(*fResult, acc);
   }

   U &PartialUpdate(unsigned int slot) { return fAggregators[slot]; }

   std::string GetActionName() { return "Aggregate"; }
};

} // end of NS RDF
} // end of NS Internal
} // end of NS ROOT

/// \endcond

#endif
