// Author: Enrico Guiraud, Danilo Piparo CERN  12/2016

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
  \defgroup dataframe Data Frame
The ROOT Data Frame allows to analyse data stored in TTrees with a high level interface.
*/


#ifndef ROOT_TDATAFRAME
#define ROOT_TDATAFRAME

#include "ROOT/RArrayView.hxx"
#include "ROOT/TDFOperations.hxx"
#include "ROOT/TDFTraitsUtils.hxx"
#include "TChain.h"
#include "TH1F.h" // For Histo actions
#include "TH2F.h" // For Histo actions
#include "TH3F.h" // For Histo actions
#include "TInterpreter.h"
#include "TProfile.h" // For Histo actions
#include "TProfile2D.h" // For Histo actions
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <algorithm> // std::find
#include <array>
#include <map>
#include <memory>
#include <string>
#include <typeinfo>
#include <vector>

namespace ROOT {

using BranchNames_t = std::vector<std::string>;

// Fwd declarations
namespace Detail {
class TDataFrameImpl;
}

namespace Experimental {

/// Smart pointer for the return type of actions
/**
\class ROOT::Experimental::TActionResultProxy
\ingroup dataframe
\brief A wrapper around the result of TDataFrame actions able to trigger calculations lazily.
\tparam T Type of the action result

A smart pointer which allows to access the result of a TDataFrame action. The
methods of the encapsulated object can be accessed via the arrow operator.
Upon invocation of the arrow operator or dereferencing (`operator*`), the
loop on the events and calculations of all scheduled actions are executed
if needed.
It is possible to iterate on the result proxy if the proxied object is a collection.
~~~{.cpp}
for (auto& myItem : myResultProxy) { ... };
~~~
If iteration is not supported by the type of the proxied object, a compilation error is thrown.

*/
template <typename T>
class TActionResultProxy {
/// \cond HIDDEN_SYMBOLS
   template<typename V, bool isCont = ROOT::Internal::TDFTraitsUtils::TIsContainer<V>::fgValue>
   struct TIterationHelper{
      using Iterator_t = void;
      void GetBegin(const V& ){static_assert(sizeof(V) == 0, "It does not make sense to ask begin for this class.");}
      void GetEnd(const V& ){static_assert(sizeof(V) == 0, "It does not make sense to ask end for this class.");}
   };

   template<typename V>
   struct TIterationHelper<V,true>{
      using Iterator_t = decltype(std::begin(std::declval<V>()));
      static Iterator_t GetBegin(const V &v) {return std::begin(v);};
      static Iterator_t GetEnd(const V &v) {return std::end(v);};
   };
/// \endcond
   using SPT_t = std::shared_ptr<T>;
   using SPTDFI_t = std::shared_ptr<ROOT::Detail::TDataFrameImpl>;
   using WPTDFI_t = std::weak_ptr<ROOT::Detail::TDataFrameImpl>;
   using ShrdPtrBool_t = std::shared_ptr<bool>;
   friend class ROOT::Detail::TDataFrameImpl;

   ShrdPtrBool_t fReadiness = std::make_shared<bool>(false); ///< State registered also in the TDataFrameImpl until the event loop is executed
   WPTDFI_t fImplWeakPtr;                                    ///< Points to the TDataFrameImpl at the root of the functional graph
   SPT_t fObjPtr;                                            ///< Shared pointer encapsulating the wrapped result

   /// Triggers the event loop in the TDataFrameImpl instance to which it's associated via the fImplWeakPtr
   void TriggerRun();

   /// Get the pointer to the encapsulated result.
   /// Ownership is not transferred to the caller.
   /// Triggers event loop and execution of all actions booked in the associated TDataFrameImpl.
   T *Get()
   {
      if (!*fReadiness) TriggerRun();
      return fObjPtr.get();
   }

   TActionResultProxy(const SPT_t &objPtr, const ShrdPtrBool_t &readiness, const SPTDFI_t &firstData)
      : fReadiness(readiness), fImplWeakPtr(firstData), fObjPtr(objPtr) { }

   /// Factory to allow to keep the constructor private
   static TActionResultProxy<T>
   MakeActionResultProxy(const SPT_t &objPtr, const ShrdPtrBool_t &readiness, const SPTDFI_t &firstData)
   {
      return TActionResultProxy(objPtr, readiness, firstData);
   }
public:
   TActionResultProxy() = delete;

   /// Get a reference to the encapsulated object.
   /// Triggers event loop and execution of all actions booked in the associated TDataFrameImpl.
   T &operator*() { return *Get(); }

   /// Get a pointer to the encapsulated object.
   /// Ownership is not transferred to the caller.
   /// Triggers event loop and execution of all actions booked in the associated TDataFrameImpl.
   T *operator->() { return Get(); }

   /// Return an iterator to the beginning of the contained object if this makes
   /// sense, throw a compilation error otherwise
   typename TIterationHelper<T>::Iterator_t
   begin()
   {
      if (!*fReadiness) TriggerRun();
      return TIterationHelper<T>::GetBegin(*fObjPtr);
   }

   /// Return an iterator to the end of the contained object if this makes
   /// sense, throw a compilation error otherwise
   typename TIterationHelper<T>::Iterator_t
   end()
   {
      if (!*fReadiness) TriggerRun();
      return TIterationHelper<T>::GetEnd(*fObjPtr);
   }
};

} // end NS Experimental

} // end NS ROOT

// Internal classes

namespace ROOT {

namespace Detail {
class TDataFrameImpl;
}

namespace Internal {

const char *ToConstCharPtr(const char *s);
const char *ToConstCharPtr(const std::string& s);
std::string ColumnName2ColumnTypeName(const std::string &colName, ROOT::Detail::TDataFrameImpl &df);
unsigned int GetNSlots();

template <typename TDFNode, typename ActionType, typename BranchType, typename ActionResultType>
ROOT::Experimental::TActionResultProxy<ActionResultType>
CallCreateAction(TDFNode* node, const BranchNames_t &bl, const std::shared_ptr<ActionResultType> &r,
                 BranchType*)
{
   return node->template CreateAction<ActionType,BranchType,ActionResultType>(bl, r, nullptr);
}

using TVBPtr_t = std::shared_ptr<TTreeReaderValueBase>;
using TVBVec_t = std::vector<TVBPtr_t>;

template<typename BranchType>
std::shared_ptr<ROOT::Internal::TTreeReaderValueBase>
ReaderValueOrArray(TTreeReader &r, const std::string &branch, TDFTraitsUtils::TTypeList<BranchType>) {
   return std::make_shared<TTreeReaderValue<BranchType>>(r, branch.c_str());
}


template<typename BranchType>
std::shared_ptr<ROOT::Internal::TTreeReaderValueBase>
ReaderValueOrArray(TTreeReader &r, const std::string &branch, TDFTraitsUtils::TTypeList<std::array_view<BranchType>>) {
   return std::make_shared<TTreeReaderArray<BranchType>>(r, branch.c_str());
}



template <int... S, typename... BranchTypes>
TVBVec_t BuildReaderValues(TTreeReader &r, const BranchNames_t &bl, const BranchNames_t &tmpbl,
                           TDFTraitsUtils::TTypeList<BranchTypes...>,
                           TDFTraitsUtils::TStaticSeq<S...>)
{
   // isTmpBranch has length bl.size(). Elements are true if the corresponding
   // branch is a temporary branch created with AddColumn, false if they are
   // actual branches present in the TTree.
   std::array<bool, sizeof...(S)> isTmpBranch;
   for (unsigned int i = 0; i < isTmpBranch.size(); ++i)
      isTmpBranch[i] = std::find(tmpbl.begin(), tmpbl.end(), bl.at(i)) != tmpbl.end();

   // Build vector of pointers to TTreeReaderValueBase.
   // tvb[i] points to a TTreeReader{Value,Array} specialized for the i-th BranchType,
   // corresponding to the i-th branch in bl
   // For temporary branches (declared with AddColumn) a nullptr is created instead
   // S is expected to be a sequence of sizeof...(BranchTypes) integers
   // Note that here TTypeList only contains one single type
   TVBVec_t tvb{isTmpBranch[S] ? nullptr : ReaderValueOrArray(r, bl.at(S), TDFTraitsUtils::TTypeList<BranchTypes>())
                ...}; // "..." expands BranchTypes and S simultaneously

   return tvb;
}

template <typename Filter>
void CheckFilter(Filter&)
{
   using FilterRet_t = typename TDFTraitsUtils::TFunctionTraits<Filter>::Ret_t;
   static_assert(std::is_same<FilterRet_t, bool>::value, "filter functions must return a bool");
}

void CheckTmpBranch(const std::string &branchName, TTree *treePtr);

///////////////////////////////////////////////////////////////////////////////
/// Check that the callable passed to TDataFrameInterface::Reduce:
/// - takes exactly two arguments of the same type
/// - has a return value of the same type as the arguments
template<typename F, typename T>
void CheckReduce(F&, ROOT::Internal::TDFTraitsUtils::TTypeList<T,T>)
{
   using Ret_t = typename ROOT::Internal::TDFTraitsUtils::TFunctionTraits<F>::Ret_t;
   static_assert(std::is_same<Ret_t, T>::value,
      "reduce function must have return type equal to argument type");
   return;
}

///////////////////////////////////////////////////////////////////////////////
/// This overload of CheckReduce is called if T is not a TTypeList<T,T>
template<typename F, typename T>
void CheckReduce(F&, T)
{
   static_assert(sizeof(F) == 0,
      "reduce function must take exactly two arguments of the same type");
}

/// Returns local BranchNames or default BranchNames according to which one should be used
const BranchNames_t &PickBranchNames(unsigned int nArgs, const BranchNames_t &bl, const BranchNames_t &defBl);

class TDataFrameActionBase {
protected:
   ROOT::Detail::TDataFrameImpl *fImplPtr; ///< A raw pointer to the TDataFrameImpl at the root of this functional graph. It is only guaranteed to contain a valid address during an event loop.
   const BranchNames_t fTmpBranches;
   std::vector<TVBVec_t> fReaderValues;
public:
   TDataFrameActionBase(ROOT::Detail::TDataFrameImpl *implPtr, const BranchNames_t &tmpBranches);
   virtual ~TDataFrameActionBase() {}
   virtual void Run(unsigned int slot, Long64_t entry) = 0;
   virtual void BuildReaderValues(TTreeReader &r, unsigned int slot) = 0;
   void CreateSlots(unsigned int nSlots);
};

using ActionBasePtr_t = std::shared_ptr<TDataFrameActionBase>;
using ActionBaseVec_t = std::vector<ActionBasePtr_t>;

// Forward declarations
template<typename T>
T &GetBranchValue(TVBPtr_t &readerValues, unsigned int slot, Long64_t entry, const std::string &branch,
                  ROOT::Detail::TDataFrameImpl *df, TDFTraitsUtils::TTypeList<T>);
template<typename T>
std::array_view<T> GetBranchValue(TVBPtr_t &readerValues, unsigned int slot, Long64_t entry, const std::string &branch,
                                  ROOT::Detail::TDataFrameImpl *df, TDFTraitsUtils::TTypeList<std::array_view<T>>);


template <typename Helper, typename PrevDataFrame,
          typename BranchTypes_t = typename Helper::BranchTypes_t>
class TDataFrameAction final : public TDataFrameActionBase {
   using TypeInd_t = typename TDFTraitsUtils::TGenStaticSeq<BranchTypes_t::fgSize>::Type_t;

   Helper fHelper;
   const BranchNames_t fBranches;
   PrevDataFrame &fPrevData;

public:
   TDataFrameAction(Helper &&h, const BranchNames_t &bl, PrevDataFrame &pd)
      : TDataFrameActionBase(pd.GetImplPtr(), pd.GetTmpBranches()),
        fHelper(std::move(h)), fBranches(bl), fPrevData(pd) { }

   TDataFrameAction(const TDataFrameAction &) = delete;

   void BuildReaderValues(TTreeReader &r, unsigned int slot) final
   {
      fReaderValues[slot] = ROOT::Internal::BuildReaderValues(r, fBranches, fTmpBranches, BranchTypes_t(), TypeInd_t());
   }

   void Run(unsigned int slot, Long64_t entry) final
   {
      // check if entry passes all filters
      if (fPrevData.CheckFilters(slot, entry))
         Exec(slot, entry, TypeInd_t(), BranchTypes_t());
   }

   template <int... S, typename... BranchTypes>
   void Exec(unsigned int slot, Long64_t entry,
             TDFTraitsUtils::TStaticSeq<S...>,
             TDFTraitsUtils::TTypeList<BranchTypes...>)
   {
      (void) entry; // avoid bogus unused-but-set-parameter warning by gcc
      // Take each pointer in tvb, cast it to a pointer to the
      // correct specialization of TTreeReaderValue, and get its content.
      // S expands to a sequence of integers 0 to sizeof...(types)-1
      // S and BranchTypes are expanded simultaneously by "..."
      fHelper.Exec(slot, GetBranchValue(fReaderValues[slot][S], slot, entry,
                                        fBranches[S], fImplPtr,
                                        TDFTraitsUtils::TTypeList<BranchTypes>())
                   ...);
   }

   ~TDataFrameAction() {
      fHelper.Finalize();
   }
};

namespace ActionTypes {
   struct Histo1D {};
   struct Min {};
   struct Max {};
   struct Mean {};
}

// Utilities to accommodate v7
namespace TDFV7Utils {

template<typename T, bool ISV7HISTO = !std::is_base_of<TH1, T>::value>
struct TIsV7Histo {
   const static bool fgValue = ISV7HISTO;
};

template<typename T, bool ISV7HISTO = TIsV7Histo<T>::fgValue>
struct Histo {
   static void SetCanExtendAllAxes(T &h)
   {
      h.SetCanExtend(::TH1::kAllAxes);
   }
   static bool HasAxisLimits(T &h)
   {
      auto xaxis = h.GetXaxis();
      return !(xaxis->GetXmin() == 0. && xaxis->GetXmax() == 0.);
   }
};

template<typename T>
struct Histo<T, true> {
   static void SetCanExtendAllAxes(T&) { }
   static bool HasAxisLimits(T&) {return true;}
};

} // end NS TDFV7Utils

} // end NS Internal

namespace Detail {

// forward declarations for TDataFrameInterface
class TDataFrameFilterBase;
template <typename F, typename PrevData>
class TDataFrameFilter;
class TDataFrameBranchBase;
template <typename F, typename PrevData>
class TDataFrameBranch;
class TDataFrameImpl;
class TDataFrameGuessedType{};
}

namespace Experimental {

class TDataFrame;

} // end Experimental
} // end ROOT

namespace cling {
   std::string printValue(ROOT::Experimental::TDataFrame *tdf); // For a nice printing at the prompt
}

namespace ROOT {
namespace Experimental {

/**
* \class ROOT::Experimental::TDataFrameInterface
* \ingroup dataframe
* \brief The public interface to the TDataFrame federation of classes: TDataFrameImpl, TDataFrameFilter, TDataFrameBranch
* \tparam T One of the TDataFrameImpl, TDataFrameFilter, TDataFrameBranch classes. The user never specifies this type manually.
*/
template <typename Proxied>
class TDataFrameInterface {
   friend std::string cling::printValue(ROOT::Experimental::TDataFrame *tdf); // For a nice printing at the prompt
   template<typename T> friend class TDataFrameInterface;
   template <typename TDFNode, typename ActionType, typename BranchType, typename ActionResultType>
   friend TActionResultProxy<ActionResultType> ROOT::Internal::CallCreateAction(TDFNode*, const BranchNames_t&, const std::shared_ptr<ActionResultType>&, BranchType*);
public:

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Append a filter to the call graph.
   /// \param[in] f Function, lambda expression, functor class or any other callable object. It must return a `bool` signalling whether the event has passed the selection (true) or not (false).
   /// \param[in] bn Names of the branches in input to the filter function.
   /// \param[in] name Optional name of this filter. See `Report`.
   ///
   /// Append a filter node at the point of the call graph corresponding to the
   /// object this method is called on.
   /// The callable `f` should not have side-effects (e.g. modification of an
   /// external or static variable) to ensure correct results when implicit
   /// multi-threading is active.
   ///
   /// TDataFrame only evaluates filters when necessary: if multiple filters
   /// are chained one after another, they are executed in order and the first
   /// one returning false causes the event to be discarded.
   /// Even if multiple actions or transformations depend on the same filter,
   /// it is executed once per entry. If its result is requested more than
   /// once, the cached result is served.
   template <typename F>
   TDataFrameInterface<ROOT::Detail::TDataFrameFilterBase>
   Filter(F f, const BranchNames_t &bn = {}, const std::string &name = "")
   {
      ROOT::Internal::CheckFilter(f);
      auto df = GetDataFrameChecked();
      const BranchNames_t &defBl = df->GetDefaultBranches();
      auto nArgs = ROOT::Internal::TDFTraitsUtils::TFunctionTraits<F>::Args_t::fgSize;
      const BranchNames_t &actualBl = ROOT::Internal::PickBranchNames(nArgs, bn, defBl);
      using DFF_t = ROOT::Detail::TDataFrameFilter<F, Proxied>;
      auto FilterPtr = std::make_shared<DFF_t> (std::move(f), actualBl, *fProxiedPtr, name);
      df->Book(FilterPtr);
      TDataFrameInterface<ROOT::Detail::TDataFrameFilterBase> tdf_f(FilterPtr, fImplWeakPtr);
      return tdf_f;
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Append a filter to the call graph.
   /// \param[in] f Function, lambda expression, functor class or any other callable object. It must return a `bool` signalling whether the event has passed the selection (true) or not (false).
   /// \param[in] name Optional name of this filter. See `Report`.
   ///
   /// Refer to the first overload of this method for the full documentation.
   template <typename F>
   TDataFrameInterface<ROOT::Detail::TDataFrameFilterBase>
   Filter(F f, const std::string &name)
   {
      return Filter(f, {}, name);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Append a filter to the call graph.
   /// \param[in] f Function, lambda expression, functor class or any other callable object. It must return a `bool` signalling whether the event has passed the selection (true) or not (false).
   /// \param[in] bn Names of the branches in input to the filter function.
   ///
   /// Refer to the first overload of this method for the full documentation.
   template <typename F>
   TDataFrameInterface<ROOT::Detail::TDataFrameFilterBase>
   Filter(F f, const std::initializer_list<std::string> &bn)
   {
      return Filter(f, BranchNames_t{bn});
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Creates a temporary branch
   /// \param[in] name The name of the temporary branch.
   /// \param[in] expression Function, lambda expression, functor class or any other callable object producing the temporary value. Returns the value that will be assigned to the temporary branch.
   /// \param[in] bl Names of the branches in input to the producer function.
   ///
   /// Create a temporary branch that will be visible from all subsequent nodes
   /// of the functional chain. The `expression` is only evaluated for entries that pass
   /// all the preceding filters.
   /// A new variable is created called `name`, accessible as if it was contained
   /// in the dataset from subsequent transformations/actions.
   ///
   /// Use cases include:
   ///
   /// * caching the results of complex calculations for easy and efficient multiple access
   /// * extraction of quantities of interest from complex objects
   /// * branch aliasing, i.e. changing the name of a branch
   ///
   /// An exception is thrown if the name of the new branch is already in use
   /// for another branch in the TTree.
   template <typename F>
   TDataFrameInterface<ROOT::Detail::TDataFrameBranchBase>
   AddColumn(const std::string &name, F expression, const BranchNames_t &bl = {})
   {
      auto df = GetDataFrameChecked();
      ROOT::Internal::CheckTmpBranch(name, df->GetTree());
      const BranchNames_t &defBl = df->GetDefaultBranches();
      auto nArgs = ROOT::Internal::TDFTraitsUtils::TFunctionTraits<F>::Args_t::fgSize;
      const BranchNames_t &actualBl = ROOT::Internal::PickBranchNames(nArgs, bl, defBl);
      using DFB_t = ROOT::Detail::TDataFrameBranch<F, Proxied>;
      auto BranchPtr = std::make_shared<DFB_t>(name, std::move(expression), actualBl, *fProxiedPtr);
      df->Book(BranchPtr);
      TDataFrameInterface<ROOT::Detail::TDataFrameBranchBase> tdf_b(BranchPtr, fImplWeakPtr);
      return tdf_b;
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Execute a user-defined function on each entry (*instant action*)
   /// \param[in] f Function, lambda expression, functor class or any other callable object performing user defined calculations.
   /// \param[in] bl Names of the branches in input to the user function.
   ///
   /// The callable `f` is invoked once per entry. This is an *instant action*:
   /// upon invocation, an event loop as well as execution of all scheduled actions
   /// is triggered.
   /// Users are responsible for the thread-safety of this callable when executing
   /// with implicit multi-threading enabled (i.e. ROOT::EnableImplicitMT).
   template <typename F>
   void Foreach(F f, const BranchNames_t &bl = {})
   {
      namespace IU = ROOT::Internal::TDFTraitsUtils;
      using Args_t = typename IU::TFunctionTraits<decltype(f)>::ArgsNoDecay_t;
      using Ret_t = typename IU::TFunctionTraits<decltype(f)>::Ret_t;
      ForeachSlot(IU::AddSlotParameter<Ret_t>(f, Args_t()), bl);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Execute a user-defined function requiring a processing slot index on each entry (*instant action*)
   /// \param[in] f Function, lambda expression, functor class or any other callable object performing user defined calculations.
   /// \param[in] bl Names of the branches in input to the user function.
   ///
   /// Same as `Foreach`, but the user-defined function takes an extra
   /// `unsigned int` as its first parameter, the *processing slot index*.
   /// This *slot index* will be assigned a different value, `0` to `poolSize - 1`,
   /// for each thread of execution.
   /// This is meant as a helper in writing thread-safe `Foreach`
   /// actions when using `TDataFrame` after `ROOT::EnableImplicitMT()`.
   /// The user-defined processing callable is able to follow different
   /// *streams of processing* indexed by the first parameter.
   /// `ForeachSlot` works just as well with single-thread execution: in that
   /// case `slot` will always be `0`.
   template<typename F>
   void ForeachSlot(F f, const BranchNames_t &bl = {})
   {
      auto df = GetDataFrameChecked();
      const BranchNames_t &defBl= df->GetDefaultBranches();
      auto nArgs = ROOT::Internal::TDFTraitsUtils::TFunctionTraits<F>::Args_t::fgSize;
      const BranchNames_t &actualBl = ROOT::Internal::PickBranchNames(nArgs-1, bl, defBl);
      using Op_t = ROOT::Internal::Operations::ForeachSlotOperation<F>;
      using DFA_t  = ROOT::Internal::TDataFrameAction<Op_t, Proxied>;
      df->Book(std::make_shared<DFA_t>(Op_t(std::move(f)), actualBl, *fProxiedPtr));
      df->Run();
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Execute a user-defined reduce operation on the values of a branch
   /// \tparam F The type of the reduce callable. Automatically deduced.
   /// \tparam T The type of the branch to apply the reduction to. Automatically deduced.
   /// \param[in] f A callable with signature `T(T,T)`
   /// \param[in] branchName The branch to be reduced. If omitted, the default branch is used instead.
   ///
   /// A reduction takes two values of a branch and merges them into one (e.g.
   /// by summing them, taking the maximum, etc). This action performs the
   /// specified reduction operation on all branch values, returning
   /// a single value of the same type. The callable f must satisfy the general
   /// requirements of a *processing function* besides having signature `T(T,T)`
   /// where `T` is the type of branch.
   ///
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultPtr documentation.
   template<typename F, typename T = typename ROOT::Internal::TDFTraitsUtils::TFunctionTraits<F>::Ret_t>
   TActionResultProxy<T>
   Reduce(F f, const std::string &branchName = {})
   {
      static_assert(std::is_default_constructible<T>::value,
         "reduce object cannot be default-constructed. Please provide an initialisation value (initValue)");
      return Reduce(std::move(f), branchName, T());
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Execute a user-defined reduce operation on the values of a branch
   /// \tparam F The type of the reduce callable. Automatically deduced.
   /// \tparam T The type of the branch to apply the reduction to. Automatically deduced.
   /// \param[in] f A callable with signature `T(T,T)`
   /// \param[in] branchName The branch to be reduced. If omitted, the default branch is used instead.
   /// \param[in] initValue The reduced object is initialised to this value rather than being default-constructed
   ///
   /// See the description of the other Reduce overload for more information.
   template<typename F, typename T = typename ROOT::Internal::TDFTraitsUtils::TFunctionTraits<F>::Ret_t>
   TActionResultProxy<T>
   Reduce(F f, const std::string &branchName, const T &initValue)
   {
      using Args_t = typename ROOT::Internal::TDFTraitsUtils::TFunctionTraits<F>::Args_t;
      ROOT::Internal::CheckReduce(f, Args_t());
      auto df = GetDataFrameChecked();
      unsigned int nSlots = df->GetNSlots();
      auto bl = GetBranchNames<T>({branchName}, "reduce branch values");
      auto redObjPtr = std::make_shared<T>(initValue);
      using Op_t = ROOT::Internal::Operations::ReduceOperation<F,T>;
      using DFA_t = typename ROOT::Internal::TDataFrameAction<Op_t, Proxied>;
      df->Book(std::make_shared<DFA_t>(Op_t(std::move(f), redObjPtr, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(redObjPtr);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Return the number of entries processed (*lazy action*)
   ///
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   TActionResultProxy<unsigned int> Count()
   {
      auto df = GetDataFrameChecked();
      unsigned int nSlots = df->GetNSlots();
      auto cSPtr = std::make_shared<unsigned int>(0);
      using Op_t = ROOT::Internal::Operations::CountOperation;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied>;
      df->Book(std::make_shared<DFA_t>(Op_t(cSPtr, nSlots), BranchNames_t({}), *fProxiedPtr));
      return df->MakeActionResultProxy(cSPtr);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Return a collection of values of a branch (*lazy action*)
   /// \tparam T The type of the branch.
   /// \tparam COLL The type of collection used to store the values.
   /// \param[in] branchName The name of the branch of which the values are to be collected
   ///
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   template <typename T, typename COLL = std::vector<T>>
   TActionResultProxy<COLL>
   Take(const std::string &branchName = "")
   {
      auto df = GetDataFrameChecked();
      unsigned int nSlots = df->GetNSlots();
      auto bl = GetBranchNames<T>({branchName}, "get the values of the branch");
      auto valuesPtr = std::make_shared<COLL>();
      using Op_t = ROOT::Internal::Operations::TakeOperation<T,COLL>;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied>;
      df->Book(std::make_shared<DFA_t>(Op_t(valuesPtr, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(valuesPtr);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return a one-dimensional histogram with the values of a branch (*lazy action*)
   /// \tparam T The type of the branch the values of which are used to fill the histogram.
   /// \param[in] valBranchName The name of the branch of which the values are to be collected.
   /// \param[in] weightBranchName The name of the branch of which the weights are to be collected.
   /// \param[in] model The model to be considered to build the new return value.
   ///
   /// If no branch type is specified and no weight branch is specified, the implementation will try to guess one.
   /// The returned histogram is independent of the input one.
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   /// The user renounces to the ownership of the model. The value to be used is the
   /// returned one.
   template <typename T = ROOT::Detail::TDataFrameGuessedType, typename W = void>
   TActionResultProxy<::TH1F> Histo1D(::TH1F &&model, const std::string &valBranchName = "", const std::string &weightBranchName ="")
   {
      auto bl = GetBranchNames<T,W>({valBranchName, weightBranchName},"fill the histogram");
      auto h = std::make_shared<::TH1F>(model);
      return Histo1DImpl<T,W>((W*)nullptr, bl,h);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return a one-dimensional histogram with the values of a branch (*lazy action*)
   /// \tparam T The type of the branch the values of which are used to fill the histogram.
   /// \param[in] valBranchName The name of the branch of which the values are to be collected.
   /// \param[in] weightBranchName The name of the branch of which the weights are to be collected.
   /// \param[in] nbins The number of bins.
   /// \param[in] minVal The lower value of the xaxis.
   /// \param[in] maxVal The upper value of the xaxis.
   ///
   /// If no branch type is specified, the implementation will try to guess one.
   ///
   /// If no axes boundaries are specified, all entries are buffered: at the end of
   /// the loop on the entries, the histogram is filled. If the axis boundaries are
   /// specified, the histogram (or histograms in the parallel case) are filled. This
   /// latter mode may result in a reduced memory footprint.
   ///
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   template <typename T = ROOT::Detail::TDataFrameGuessedType, typename W = void>
   TActionResultProxy<::TH1F> Histo1D(const std::string &valBranchName = "", int nBins = 128, double minVal = 0., double maxVal = 0., const std::string &weightBranchName ="")
   {
      auto bl = GetBranchNames<T,W>({valBranchName, weightBranchName},"fill the histogram");
      auto blSize = bl.size();
      ::TH1F h("", "", nBins, minVal, maxVal);
      if (minVal == maxVal) {
         ROOT::Internal::TDFV7Utils::Histo<::TH1F>::SetCanExtendAllAxes(h);
      }

      // A weighted histogram
      return Histo1D<T,W>(std::move(h), bl[0], blSize == 1 ? "" : bl[1]);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return a one-dimensional histogram with the values of a branch (*lazy action*)
   /// \tparam T The type of the branch the values of which are used to fill the histogram.
   /// \param[in] valBranchName The name of the branch of which the values are to be collected.
   /// \param[in] weightBranchName The name of the branch of which the weights are to be collected.
   ///
   /// If no branch type is specified, the implementation will try to guess one.
   ///
   /// If no axes boundaries are specified, all entries are buffered: at the end of
   /// the loop on the entries, the histogram is filled. If the axis boundaries are
   /// specified, the histogram (or histograms in the parallel case) are filled. This
   /// latter mode may result in a reduced memory footprint.
   ///
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   template <typename T = ROOT::Detail::TDataFrameGuessedType, typename W = void>
   TActionResultProxy<::TH1F> Histo1D(const std::string &valBranchName, const std::string &weightBranchName)
   {
      return Histo1D<T,W>(valBranchName, 128, 0., 0., weightBranchName);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return a two-dimensional histogram with the values of a branch (*lazy action*)
   /// \tparam B0 The type of the branch the values of which are used to fill the histogram.
   /// \tparam B1 The type of the branch the values of which are used to fill the histogram.
   /// \param[in] model The model to be considered to build the new return value.
   /// \param[in] b0BranchName The name of the branch of which the x values are to be collected.
   /// \param[in] b1BranchName The name of the branch of which the y values are to be collected.
   ///
   /// The returned histogram is independent of the input one.
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   /// The user renounces to the ownership of the model. The value to be used is the
   /// returned one.
   template <typename B0, typename B1>
   TActionResultProxy<::TH2F> Histo2D(::TH2F &&model, const std::string &b0BranchName = "", const std::string &b1BranchName = "")
   {
      auto h = std::make_shared<::TH2F>(model);
      if (!ROOT::Internal::TDFV7Utils::Histo<::TH2F>::HasAxisLimits(*h)) {
         throw std::runtime_error("2D histograms with no axes limits are not supported yet.");
      }
      auto bl = GetBranchNames<B0,B1>({b0BranchName, b1BranchName},"fill the histogram");
      using Op_t = ROOT::Internal::Operations::FillTOOperation<::TH2F>;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<B0,B1>>;
      auto df = GetDataFrameChecked();
      auto nSlots = df->GetNSlots();
      df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(h);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return a two-dimensional histogram with the values of a branch (*lazy action*)
   /// \tparam B0 The type of the branch the values of which are used to fill the histogram.
   /// \tparam B1 The type of the branch the values of which are used to fill the histogram.
   /// \tparam W The type of the branch the weights of which are used to fill the histogram.
   /// \param[in] model The model to be considered to build the new return value.
   /// \param[in] b0BranchName The name of the branch of which the x values are to be collected.
   /// \param[in] b1BranchName The name of the branch of which the y values are to be collected.
   /// \param[in] wBranchName The name of the branch of which the weight values are to be collected.
   ///
   /// The returned histogram is independent of the input one.
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   /// The user renounces to the ownership of the model. The value to be used is the
   /// returned one.
   template <typename B0, typename B1, typename W>
   TActionResultProxy<::TH2F> Histo2D(::TH2F &&model, const std::string &b0BranchName = "", const std::string &b1BranchName = "", const std::string &wBranchName = "")
   {
      auto h = std::make_shared<::TH2F>(model);
      if (!ROOT::Internal::TDFV7Utils::Histo<::TH2F>::HasAxisLimits(*h)) {
         throw std::runtime_error("2D histograms with no axes limits are not supported yet.");
      }
      auto bl = GetBranchNames<B0,B1,W>({b0BranchName, b1BranchName, wBranchName},"fill the histogram");
      using Op_t = ROOT::Internal::Operations::FillTOOperation<::TH2F>;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<B0,B1,W>>;
      auto df = GetDataFrameChecked();
      auto nSlots = df->GetNSlots();
      df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(h);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return a three-dimensional histogram with the values of a branch (*lazy action*)
   /// \tparam B0 The type of the branch the values of which are used to fill the histogram.
   /// \tparam B1 The type of the branch the values of which are used to fill the histogram.
   /// \tparam B2 The type of the branch the values of which are used to fill the histogram.
   /// \param[in] model The model to be considered to build the new return value.
   /// \param[in] b0BranchName The name of the branch of which the x values are to be collected.
   /// \param[in] b1BranchName The name of the branch of which the y values are to be collected.
   /// \param[in] b2BranchName The name of the branch of which the z values are to be collected.
   ///
   /// The returned histogram is independent of the input one.
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   /// The user renounces to the ownership of the model. The value to be used is the
   /// returned one.
   template <typename B0, typename B1, typename B2>
   TActionResultProxy<::TH3F> Histo3D(::TH3F &&model, const std::string &b0BranchName = "", const std::string &b1BranchName = "", const std::string &b2BranchName = "", const std::string &wBranchName = "")
   {
      auto h = std::make_shared<::TH3F>(model);
      if (!ROOT::Internal::TDFV7Utils::Histo<::TH3F>::HasAxisLimits(*h)) {
         throw std::runtime_error("2D histograms with no axes limits are not supported yet.");
      }
      auto bl = GetBranchNames<B0,B1,B2>({b0BranchName, b1BranchName, b2BranchName, wBranchName},"fill the histogram");
      using Op_t = ROOT::Internal::Operations::FillTOOperation<::TH3F>;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<B0,B1,B2>>;
      auto df = GetDataFrameChecked();
      auto nSlots = df->GetNSlots();
      df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(h);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return a three-dimensional histogram with the values of a branch (*lazy action*)
   /// \tparam B0 The type of the branch the values of which are used to fill the histogram.
   /// \tparam B1 The type of the branch the values of which are used to fill the histogram.
   /// \tparam B2 The type of the branch the values of which are used to fill the histogram.
   /// \tparam W The type of the branch the weights of which are used to fill the histogram.
   /// \param[in] model The model to be considered to build the new return value.
   /// \param[in] b0BranchName The name of the branch of which the x values are to be collected.
   /// \param[in] b1BranchName The name of the branch of which the y values are to be collected.
   /// \param[in] b2BranchName The name of the branch of which the z values are to be collected.
   /// \param[in] wBranchName The name of the branch of which the weight values are to be collected.
   ///
   /// The returned histogram is independent of the input one.
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   /// The user renounces to the ownership of the model. The value to be used is the
   /// returned one.
   template <typename B0, typename B1, typename B2, typename W>
   TActionResultProxy<::TH3F> Histo3D(::TH3F &&model, const std::string &b0BranchName = "", const std::string &b1BranchName = "", const std::string &b2BranchName = "", const std::string &wBranchName = "")
   {
      auto h = std::make_shared<::TH3F>(model);
      if (!ROOT::Internal::TDFV7Utils::Histo<::TH3F>::HasAxisLimits(*h)) {
         throw std::runtime_error("2D histograms with no axes limits are not supported yet.");
      }
      auto bl = GetBranchNames<B0,B1,B2,W>({b0BranchName, b1BranchName, b2BranchName, wBranchName},"fill the histogram");
      using Op_t = ROOT::Internal::Operations::FillTOOperation<::TH3F>;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<B0,B1,B2,W>>;
      auto df = GetDataFrameChecked();
      auto nSlots = df->GetNSlots();
      df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(h);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return a profile (*lazy action*)
   /// \tparam B0 The type of the branch the values of which are used to fill the profile.
   /// \tparam B1 The type of the branch the values of which are used to fill the profile.
   /// \param[in] model The model to be considered to build the new return value.
   /// \param[in] b0BranchName The name of the branch of which the x values are to be collected.
   /// \param[in] b1BranchName The name of the branch of which the y values are to be collected.
   ///
   /// The returned profile is independent of the input one.
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   /// The user renounces to the ownership of the model. The value to be used is the
   /// returned one.
   template <typename B0, typename B1>
   TActionResultProxy<::TProfile> Profile1D(::TProfile &&model, const std::string &b0BranchName = "", const std::string &b1BranchName = "")
   {
      auto h = std::make_shared<::TProfile>(model);
      if (!ROOT::Internal::TDFV7Utils::Histo<::TProfile>::HasAxisLimits(*h)) {
         throw std::runtime_error("Profiles with no axes limits are not supported yet.");
      }
      auto bl = GetBranchNames<B0,B1>({b0BranchName, b1BranchName},"fill the profile");
      using Op_t = ROOT::Internal::Operations::FillTOOperation<::TProfile>;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<B0,B1>>;
      auto df = GetDataFrameChecked();
      auto nSlots = df->GetNSlots();
      df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(h);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return a profile (*lazy action*)
   /// \tparam B0 The type of the branch the values of which are used to fill the profile.
   /// \tparam B1 The type of the branch the values of which are used to fill the profile.
   /// \tparam W The type of the branch the weights of which are used to fill the profile.
   /// \param[in] model The model to be considered to build the new return value.
   /// \param[in] b0BranchName The name of the branch of which the x values are to be collected.
   /// \param[in] b1BranchName The name of the branch of which the y values are to be collected.
   /// \param[in] wBranchName The name of the branch of which the weight values are to be collected.
   ///
   /// The returned profile is independent of the input one.
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   /// The user renounces to the ownership of the model. The value to be used is the
   /// returned one.
   template <typename B0, typename B1, typename W>
   TActionResultProxy<::TProfile> Profile1D(::TProfile &&model, const std::string &b0BranchName = "", const std::string &b1BranchName = "", const std::string &wBranchName = "")
   {
      auto h = std::make_shared<::TProfile>(model);
      if (!ROOT::Internal::TDFV7Utils::Histo<::TProfile>::HasAxisLimits(*h)) {
         throw std::runtime_error("Profiles with no axes limits are not supported yet.");
      }
      auto bl = GetBranchNames<B0,B1,W>({b0BranchName, b1BranchName, wBranchName},"fill the profile");
      using Op_t = ROOT::Internal::Operations::FillTOOperation<::TProfile>;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<B0,B1,W>>;
      auto df = GetDataFrameChecked();
      auto nSlots = df->GetNSlots();
      df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(h);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return a two-dimensional profile (*lazy action*)
   /// \tparam B0 The type of the branch the values of which are used to fill the profile.
   /// \tparam B1 The type of the branch the values of which are used to fill the profile.
   /// \param[in] model The model to be considered to build the new return value.
   /// \param[in] b0BranchName The name of the branch of which the x values are to be collected.
   /// \param[in] b1BranchName The name of the branch of which the y values are to be collected.
   ///
   /// The returned profile is independent of the input one.
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   /// The user renounces to the ownership of the model. The value to be used is the
   /// returned one.
   template <typename B0, typename B1>
   TActionResultProxy<::TProfile2D> Profile2D(::TProfile2D &&model, const std::string &b0BranchName = "", const std::string &b1BranchName = "")
   {
      auto h = std::make_shared<::TProfile2D>(model);
      if (!ROOT::Internal::TDFV7Utils::Histo<::TProfile2D>::HasAxisLimits(*h)) {
         throw std::runtime_error("2D profiles with no axes limits are not supported yet.");
      }
      auto bl = GetBranchNames<B0,B1>({b0BranchName, b1BranchName},"fill the profile");
      using Op_t = ROOT::Internal::Operations::FillTOOperation<::TProfile2D>;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<B0,B1>>;
      auto df = GetDataFrameChecked();
      auto nSlots = df->GetNSlots();
      df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(h);
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return a two-dimensional profile (*lazy action*)
   /// \tparam B0 The type of the branch the values of which are used to fill the profile.
   /// \tparam B1 The type of the branch the values of which are used to fill the profile.
   /// \tparam W The type of the branch the weights of which are used to fill the profile.
   /// \param[in] model The model to be considered to build the new return value.
   /// \param[in] b0BranchName The name of the branch of which the x values are to be collected.
   /// \param[in] b1BranchName The name of the branch of which the y values are to be collected.
   /// \param[in] wBranchName The name of the branch of which the weight values are to be collected.
   ///
   /// The returned profile is independent of the input one.
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   /// The user renounces to the ownership of the model. The value to be used is the
   /// returned one.
   template <typename B0, typename B1, typename W>
   TActionResultProxy<::TProfile2D> Profile2D(::TProfile2D &&model, const std::string &b0BranchName = "", const std::string &b1BranchName = "", const std::string &wBranchName = "")
   {
      auto h = std::make_shared<::TProfile2D>(model);
      if (!ROOT::Internal::TDFV7Utils::Histo<::TProfile2D>::HasAxisLimits(*h)) {
         throw std::runtime_error("2D profiles with no axes limits are not supported yet.");
      }
      auto bl = GetBranchNames<B0,B1,W>({b0BranchName, b1BranchName, wBranchName},"fill the profile");
      using Op_t = ROOT::Internal::Operations::FillTOOperation<::TProfile2D>;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<B0,B1,W>>;
      auto df = GetDataFrameChecked();
      auto nSlots = df->GetNSlots();
      df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(h);
   }


   ////////////////////////////////////////////////////////////////////////////
   /// \brief Fill and return any entity with a Fill method (*lazy action*)
   /// \tparam BRANCHTYPES The types of the branches the values of which are used to fill the object.
   /// \param[in] model The model to be considered to build the new return value.
   /// \param[in] bl The name of the branches read to fill the object.
   ///
   /// The returned object is independent of the input one.
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   /// The user renounces to the ownership of the model. The value to be used is the
   /// returned one.
   /// It is compulsory to express the branches to be considered.
   template <typename... BRANCHTYPES, typename T>
   TActionResultProxy<T> Fill(T &&model, const BranchNames_t &bl)
   {
      auto h = std::make_shared<T>(model);
      if (!ROOT::Internal::TDFV7Utils::Histo<T>::HasAxisLimits(*h)) {
         throw std::runtime_error("The absence of axes limits is not supported yet.");
      }
      using Op_t = ROOT::Internal::Operations::FillTOOperation<T>;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<BRANCHTYPES...>>;
      auto df = GetDataFrameChecked();
      auto nSlots = df->GetNSlots();
      df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(h);
   }


   ////////////////////////////////////////////////////////////////////////////
   /// \brief Return the minimum of processed branch values (*lazy action*)
   /// \tparam T The type of the branch.
   /// \param[in] branchName The name of the branch to be treated.
   ///
   /// If no branch type is specified, the implementation will try to guess one.
   ///
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   template <typename T = ROOT::Detail::TDataFrameGuessedType>
   TActionResultProxy<double> Min(const std::string &branchName = "")
   {
      auto bl = GetBranchNames<T>({branchName}, "calculate the minimum");
      auto minV = std::make_shared<double>(std::numeric_limits<double>::max());
      return CreateAction<ROOT::Internal::ActionTypes::Min>(bl, minV, (T*)(nullptr));
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Return the maximum of processed branch values (*lazy action*)
   /// \tparam T The type of the branch.
   /// \param[in] branchName The name of the branch to be treated.
   ///
   /// If no branch type is specified, the implementation will try to guess one.
   ///
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   template <typename T = ROOT::Detail::TDataFrameGuessedType>
   TActionResultProxy<double> Max(const std::string &branchName = "")
   {
      auto bl = GetBranchNames<T>({branchName}, "calculate the maximum");
      auto maxV = std::make_shared<double>(std::numeric_limits<double>::min());
      return CreateAction<ROOT::Internal::ActionTypes::Max>(bl, maxV, (T*)(nullptr));
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Return the mean of processed branch values (*lazy action*)
   /// \tparam T The type of the branch.
   /// \param[in] branchName The name of the branch to be treated.
   ///
   /// If no branch type is specified, the implementation will try to guess one.
   ///
   /// This action is *lazy*: upon invocation of this method the calculation is
   /// booked but not executed. See TActionResultProxy documentation.
   template <typename T = ROOT::Detail::TDataFrameGuessedType>
   TActionResultProxy<double> Mean(const std::string &branchName = "")
   {
      auto bl = GetBranchNames<T>({branchName}, "calculate the mean");
      auto meanV = std::make_shared<double>(0);
      return CreateAction<ROOT::Internal::ActionTypes::Mean>(bl, meanV, (T*)(nullptr));
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Print filtering statistics on screen
   ///
   /// Calling `Report` on the main `TDataFrame` object prints stats for
   /// all named filters in the call graph. Calling this method on a
   /// stored chain state (i.e. a graph node different from the first) prints
   /// the stats for all named filters in the chain section between the original
   /// `TDataFrame` and that node (included). Stats are printed in the same
   /// order as the named filters have been added to the graph.
   void Report() {
      auto df = GetDataFrameChecked();
      if (!df->HasRunAtLeastOnce())
         df->Run();
      fProxiedPtr->Report();
   }

private:

   inline const char* GetNodeTypeName() {return "";};

   /// Returns the default branches if needed, takes care of the error handling.
   template<typename T1, typename T2 = void, typename T3 = void, typename T4 = void>
   BranchNames_t GetBranchNames(BranchNames_t bl, const std::string &actionNameForErr)
   {
      constexpr auto isT2Void = std::is_same<T2, void>::value;
      constexpr auto isT3Void = std::is_same<T3, void>::value;
      constexpr auto isT4Void = std::is_same<T4, void>::value;

      unsigned int neededBranches = 1 + !isT2Void + !isT3Void + !isT4Void;

      unsigned int providedBranches = 0;
      std::for_each(bl.begin(), bl.end(), [&providedBranches](const std::string &s) {if (!s.empty()) providedBranches++;} );

      if (neededBranches == providedBranches) return bl;

      return GetDefaultBranchNames(neededBranches, actionNameForErr);
   }

   // Two overloaded template methods which allow to avoid branching in the Histo1D method.
   // W == void: histogram w/o weights
   template<typename X, typename W>
   TActionResultProxy<::TH1F>
   Histo1DImpl(void*, const BranchNames_t &bl, const std::shared_ptr<::TH1F> &h)
   {
      // perform type guessing if needed and build the action
      return CreateAction<ROOT::Internal::ActionTypes::Histo1D>(bl, h, (X*)(nullptr));
   }

   // W != void: histogram w/ weights
   template<typename X, typename W>
   TActionResultProxy<::TH1F>
   Histo1DImpl(W*, const BranchNames_t &bl, const std::shared_ptr<::TH1F> &h)
   {
      // weighted histograms never need to do type guessing, we can build
      // the action here
      auto df = GetDataFrameChecked();
      auto hasAxisLimits = ROOT::Internal::TDFV7Utils::Histo<::TH1F>::HasAxisLimits(*h);
      auto nSlots = df->GetNSlots();
      if (hasAxisLimits) {
         using Op_t = ROOT::Internal::Operations::FillTOOperation<::TH1F>;
         using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<X,W>>;
         df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      } else {
         using Op_t = ROOT::Internal::Operations::FillOperation;
         using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<X,W>>;
         df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      }
      return df->MakeActionResultProxy(h);
   }

   /// \cond HIDDEN_SYMBOLS
   template <typename BranchType>
   TActionResultProxy<::TH1F>
   BuildAndBook(const BranchNames_t &bl, const std::shared_ptr<::TH1F> &h,
                unsigned int nSlots, ROOT::Internal::ActionTypes::Histo1D*)
   {
      auto df = GetDataFrameChecked();
      auto hasAxisLimits = ROOT::Internal::TDFV7Utils::Histo<::TH1F>::HasAxisLimits(*h);

      if (hasAxisLimits) {
         using Op_t = ROOT::Internal::Operations::FillTOOperation<::TH1F>;
         using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<BranchType>>;
         df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      } else {
         using Op_t = ROOT::Internal::Operations::FillOperation;
         using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<BranchType>>;
         df->Book(std::make_shared<DFA_t>(Op_t(h, nSlots), bl, *fProxiedPtr));
      }
      return df->MakeActionResultProxy(h);
   }

   template <typename BranchType>
   TActionResultProxy<double>
   BuildAndBook(const BranchNames_t &bl, const std::shared_ptr<double> &minV,
                unsigned int nSlots, ROOT::Internal::ActionTypes::Min*)
   {
      using Op_t = ROOT::Internal::Operations::MinOperation;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<BranchType>>;
      auto df = GetDataFrameChecked();
      df->Book(std::make_shared<DFA_t>(Op_t(minV, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(minV);
   }

   template <typename BranchType>
   TActionResultProxy<double>
   BuildAndBook(const BranchNames_t &bl, const std::shared_ptr<double> &maxV,
                unsigned int nSlots, ROOT::Internal::ActionTypes::Max*)
   {
      using Op_t = ROOT::Internal::Operations::MaxOperation;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<BranchType>>;
      auto df = GetDataFrameChecked();
      df->Book(std::make_shared<DFA_t>(Op_t(maxV, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(maxV);
   }


   template <typename BranchType>
   TActionResultProxy<double>
   BuildAndBook(const BranchNames_t &bl, const std::shared_ptr<double> &meanV,
                unsigned int nSlots, ROOT::Internal::ActionTypes::Mean*)
   {
      using Op_t = ROOT::Internal::Operations::MeanOperation;
      using DFA_t = ROOT::Internal::TDataFrameAction<Op_t, Proxied, ROOT::Internal::TDFTraitsUtils::TTypeList<BranchType>>;
      auto df = GetDataFrameChecked();
      df->Book(std::make_shared<DFA_t>(Op_t(meanV, nSlots), bl, *fProxiedPtr));
      return df->MakeActionResultProxy(meanV);
   }
   /// \endcond

   // Type was specified by the user, no need to guess it
   template <typename ActionType, typename BranchType, typename ActionResultType>
   TActionResultProxy<ActionResultType>
   CreateAction(const BranchNames_t &bl, const std::shared_ptr<ActionResultType> &r,
                BranchType*)
   {
      auto df = GetDataFrameChecked();
      unsigned int nSlots = df->GetNSlots();
      return BuildAndBook<BranchType>(bl, r, nSlots, (ActionType*)nullptr);
   }

   // User did not specify type, do type guessing
   template <typename ActionType, typename ActionResultType>
   TActionResultProxy<ActionResultType>
   CreateAction(const BranchNames_t &bl, const std::shared_ptr<ActionResultType> &r,
                ROOT::Detail::TDataFrameGuessedType*)
   {
      gInterpreter->ProcessLine("#include \"ROOT/TDataFrame.hxx\"");
      auto df = GetDataFrameChecked();
      const auto& theBranchName = bl[0];
      const auto theBranchTypeName = ROOT::Internal::ColumnName2ColumnTypeName(theBranchName, *df);
      const auto actionResultTypeName = TClass::GetClass(typeid(std::shared_ptr<ActionResultType>))->GetName();
      const auto actionTypeName = TClass::GetClass(typeid(ActionType))->GetName();
      if (theBranchTypeName.empty()) {
         std::string exceptionText = "The type of branch ";
         exceptionText += theBranchName;
         exceptionText += " could not be guessed. Please specify one.";
         throw std::runtime_error(exceptionText.c_str());
      }
      std::stringstream createAction_str;

      createAction_str << "ROOT::Internal::CallCreateAction<"
                       << GetNodeTypeName() << ", "
                       << actionTypeName << ", "
                       << theBranchTypeName << ", "
                       << actionResultTypeName << "::element_type>("
                       << "(" << GetNodeTypeName() << "*)" << this << ", "
                       << "*(ROOT::BranchNames_t*)" << &bl << ", "
                       << "*(" << actionResultTypeName << "*)" << &r << ", "
                       << "nullptr);";
      return *(TActionResultProxy<ActionResultType>*) gInterpreter->ProcessLine(createAction_str.str().c_str());
   }

protected:
   /// Get the TDataFrameImpl if reachable. If not, throw.
   std::shared_ptr<ROOT::Detail::TDataFrameImpl>
   GetDataFrameChecked()
   {
      auto df = fImplWeakPtr.lock();
      if (!df) {
         throw std::runtime_error("The main TDataFrame is not reachable: did it go out of scope?");
      }
      return df;
   }

   const BranchNames_t GetDefaultBranchNames(unsigned int nExpectedBranches, const std::string &actionNameForErr)
   {
      auto df = GetDataFrameChecked();
      const BranchNames_t &defaultBranches = df->GetDefaultBranches();
      const auto dBSize = defaultBranches.size();
      if (nExpectedBranches > dBSize) {
         std::string msg("Trying to deduce the branches from the default list in order to ");
         msg += actionNameForErr;
         msg += ". A set of branches of size ";
         msg += std::to_string(dBSize);
         msg += " was found. ";
         msg += std::to_string(nExpectedBranches);
         msg += 1 != nExpectedBranches ? " are" : " is";
         msg += " needed. Please specify the branches explicitly.";
         throw std::runtime_error(msg);
      }
      auto bnBegin = defaultBranches.begin();
      return BranchNames_t(bnBegin, bnBegin + nExpectedBranches);
   }

   TDataFrameInterface(const std::shared_ptr<Proxied> &proxied,
                       const std::weak_ptr<ROOT::Detail::TDataFrameImpl> &impl)
      : fProxiedPtr(proxied), fImplWeakPtr(impl) {}

   /// Only enabled when building a TDataFrameInterface<TDataFrameImpl>
   template<typename T = Proxied, typename std::enable_if<std::is_same<T,ROOT::Detail::TDataFrameImpl>::value, int>::type = 0>
   TDataFrameInterface(const std::shared_ptr<Proxied> &proxied)
      : fProxiedPtr(proxied), fImplWeakPtr(proxied->GetSharedPtr()) {}

   std::shared_ptr<Proxied> fProxiedPtr;
   std::weak_ptr<ROOT::Detail::TDataFrameImpl> fImplWeakPtr;
};

class TDataFrame : public TDataFrameInterface<ROOT::Detail::TDataFrameImpl> {
private:
   std::shared_ptr<TTree> fTree;
   void InitTree(TTree &tree, bool ownsTree);
public:
   TDataFrame(const std::string &treeName, const std::string &filenameglob, const BranchNames_t &defaultBranches = {});
   ////////////////////////////////////////////////////////////////////////////
   /// \brief Build the dataframe
   /// \tparam FILENAMESCOLL The type of the file collection: only requirement: must have begin and end.
   /// \param[in] treeName Name of the tree contained in the directory
   /// \param[in] filenamescoll Collection of file names, for example a list of strings.
   /// \param[in] defaultBranches Collection of default branches.
   ///
   /// The default branches are looked at in case no branch is specified in the
   /// booking of actions or transformations.
   /// See ROOT::Experimental::TDataFrameInterface for the documentation of the
   /// methods available.
   template<typename FILENAMESCOLL, typename std::enable_if<ROOT::Internal::TDFTraitsUtils::TIsContainer<FILENAMESCOLL>::fgValue, int>::type = 0>
   TDataFrame(const std::string &treeName, const FILENAMESCOLL &filenamescoll, const BranchNames_t &defaultBranches = {});
   TDataFrame(const std::string &treeName, ::TDirectory *dirPtr, const BranchNames_t &defaultBranches = {});
   TDataFrame(TTree &tree, const BranchNames_t &defaultBranches = {});
};

} // end NS Experimental

namespace Detail {

class TDataFrameBranchBase {
protected:
   TDataFrameImpl *fImplPtr; ///< A raw pointer to the TDataFrameImpl at the root of this functional graph. It is only guaranteed to contain a valid address during an event loop.
   BranchNames_t fTmpBranches;
   const std::string fName;
public:
   TDataFrameBranchBase(TDataFrameImpl *df, const BranchNames_t &tmpBranches, const std::string &name);
   virtual ~TDataFrameBranchBase() {}
   virtual void BuildReaderValues(TTreeReader &r, unsigned int slot) = 0;
   virtual void CreateSlots(unsigned int nSlots) = 0;
   virtual void *GetValue(unsigned int slot, Long64_t entry) = 0;
   virtual const std::type_info &GetTypeId() const = 0;
   virtual bool CheckFilters(unsigned int slot, Long64_t entry) = 0;
   TDataFrameImpl *GetImplPtr() const;
   virtual void Report() const = 0;
   virtual void PartialReport() const = 0;
   std::string GetName() const;
   BranchNames_t GetTmpBranches() const;
};
using TmpBranchBasePtr_t = std::shared_ptr<TDataFrameBranchBase>;

template <typename F, typename PrevData>
class TDataFrameBranch final : public TDataFrameBranchBase {
   using BranchTypes_t = typename ROOT::Internal::TDFTraitsUtils::TFunctionTraits<F>::Args_t;
   using TypeInd_t = typename ROOT::Internal::TDFTraitsUtils::TGenStaticSeq<BranchTypes_t::fgSize>::Type_t;
   using Ret_t = typename ROOT::Internal::TDFTraitsUtils::TFunctionTraits<F>::Ret_t;

   F fExpression;
   const BranchNames_t fBranches;

   std::vector<ROOT::Internal::TVBVec_t> fReaderValues;
   std::vector<std::shared_ptr<Ret_t>> fLastResultPtr;
   PrevData &fPrevData;
   std::vector<Long64_t> fLastCheckedEntry = {-1};

public:
   TDataFrameBranch(const std::string &name, F &&expression, const BranchNames_t &bl, PrevData &pd)
      : TDataFrameBranchBase(pd.GetImplPtr(), pd.GetTmpBranches(), name),
        fExpression(std::move(expression)), fBranches(bl), fPrevData(pd)
   {
      fTmpBranches.emplace_back(name);
   }

   TDataFrameBranch(const TDataFrameBranch &) = delete;

   void BuildReaderValues(TTreeReader &r, unsigned int slot) final
   {
      fReaderValues[slot] = ROOT::Internal::BuildReaderValues(r, fBranches, fTmpBranches, BranchTypes_t(), TypeInd_t());
   }

   void *GetValue(unsigned int slot, Long64_t entry) final
   {
      if (entry != fLastCheckedEntry[slot]) {
         // evaluate this filter, cache the result
         auto newValuePtr = GetValueHelper(BranchTypes_t(), TypeInd_t(), slot, entry);
         fLastResultPtr[slot] = newValuePtr;
         fLastCheckedEntry[slot] = entry;
      }
      return static_cast<void *>(fLastResultPtr[slot].get());
   }

   const std::type_info &GetTypeId() const { return typeid(Ret_t); }

   void CreateSlots(unsigned int nSlots) final
   {
      fReaderValues.resize(nSlots);
      fLastCheckedEntry.resize(nSlots, -1);
      fLastResultPtr.resize(nSlots);
   }

   bool CheckFilters(unsigned int slot, Long64_t entry) final
   {
      // dummy call: it just forwards to the previous object in the chain
      return fPrevData.CheckFilters(slot, entry);
   }

   template <int... S, typename... BranchTypes>
   std::shared_ptr<Ret_t>
   GetValueHelper(ROOT::Internal::TDFTraitsUtils::TTypeList<BranchTypes...>,
                  ROOT::Internal::TDFTraitsUtils::TStaticSeq<S...>,
                  unsigned int slot, Long64_t entry)
   {
      auto valuePtr = std::make_shared<Ret_t>(fExpression(
         ROOT::Internal::GetBranchValue(fReaderValues[slot][S], slot, entry, fBranches[S],
                                        fImplPtr, ROOT::Internal::TDFTraitsUtils::TTypeList<BranchTypes>())...));
      return valuePtr;
   }

   // recursive chain of `Report`s
   // TDataFrameBranch simply forwards the call to the previous node
   void Report() const final {
      fPrevData.PartialReport();
   }

   void PartialReport() const final {
      fPrevData.PartialReport();
   }

};

class TDataFrameFilterBase {
protected:
   TDataFrameImpl *fImplPtr; ///< A raw pointer to the TDataFrameImpl at the root of this functional graph. It is only guaranteed to contain a valid address during an event loop.
   const BranchNames_t fTmpBranches;
   std::vector<ROOT::Internal::TVBVec_t> fReaderValues = {};
   std::vector<Long64_t> fLastCheckedEntry = {-1};
   std::vector<int> fLastResult = {true}; // std::vector<bool> cannot be used in a MT context safely
   std::vector<ULong64_t> fAccepted = {0};
   std::vector<ULong64_t> fRejected = {0};
   const std::string fName;

public:
   TDataFrameFilterBase(TDataFrameImpl *df, const BranchNames_t &tmpBranches, const std::string &name);
   virtual ~TDataFrameFilterBase() {}
   virtual void BuildReaderValues(TTreeReader &r, unsigned int slot) = 0;
   virtual bool CheckFilters(unsigned int slot, Long64_t entry) = 0;
   virtual void Report() const = 0;
   virtual void PartialReport() const = 0;
   TDataFrameImpl *GetImplPtr() const;
   BranchNames_t GetTmpBranches() const;
   bool HasName() const;
   void CreateSlots(unsigned int nSlots);
   void PrintReport() const;
};
using FilterBasePtr_t = std::shared_ptr<TDataFrameFilterBase>;
using FilterBaseVec_t = std::vector<FilterBasePtr_t>;

template <typename FilterF, typename PrevDataFrame>
class TDataFrameFilter final : public TDataFrameFilterBase {
   using BranchTypes_t = typename ROOT::Internal::TDFTraitsUtils::TFunctionTraits<FilterF>::Args_t;
   using TypeInd_t = typename ROOT::Internal::TDFTraitsUtils::TGenStaticSeq<BranchTypes_t::fgSize>::Type_t;

   FilterF fFilter;
   const BranchNames_t fBranches;
   PrevDataFrame &fPrevData;

public:
   TDataFrameFilter(FilterF &&f, const BranchNames_t &bl,
                    PrevDataFrame &pd, const std::string &name = "")
      : TDataFrameFilterBase(pd.GetImplPtr(), pd.GetTmpBranches(), name),
        fFilter(std::move(f)), fBranches(bl), fPrevData(pd) { }

   TDataFrameFilter(const TDataFrameFilter &) = delete;

   bool CheckFilters(unsigned int slot, Long64_t entry) final
   {
      if (entry != fLastCheckedEntry[slot]) {
         if (!fPrevData.CheckFilters(slot, entry)) {
            // a filter upstream returned false, cache the result
            fLastResult[slot] = false;
         } else {
            // evaluate this filter, cache the result
            auto passed = CheckFilterHelper(BranchTypes_t(), TypeInd_t(), slot, entry);
            passed ? ++fAccepted[slot] : ++fRejected[slot];
            fLastResult[slot] = passed;
         }
         fLastCheckedEntry[slot] = entry;
      }
      return fLastResult[slot];
   }

   template <int... S, typename... BranchTypes>
   bool CheckFilterHelper(ROOT::Internal::TDFTraitsUtils::TTypeList<BranchTypes...>,
                          ROOT::Internal::TDFTraitsUtils::TStaticSeq<S...>,
                          unsigned int slot, Long64_t entry)
   {
      // Take each pointer in tvb, cast it to a pointer to the
      // correct specialization of TTreeReaderValue, and get its content.
      // S expands to a sequence of integers 0 to `sizeof...(types)-1
      // S and types are expanded simultaneously by "..."
      (void) slot; // avoid bogus unused-but-set-parameter warning by gcc
      (void) entry; // avoid bogus unused-but-set-parameter warning by gcc
      return fFilter(ROOT::Internal::GetBranchValue(fReaderValues[slot][S], slot, entry, fBranches[S],
                     fImplPtr, ROOT::Internal::TDFTraitsUtils::TTypeList<BranchTypes>())...);
   }

   void BuildReaderValues(TTreeReader &r, unsigned int slot) final
   {
      fReaderValues[slot] = ROOT::Internal::BuildReaderValues(r, fBranches, fTmpBranches, BranchTypes_t(), TypeInd_t());
   }


   // recursive chain of `Report`s
   void Report() const final {
      PartialReport();
   }

   void PartialReport() const final {
      fPrevData.PartialReport();
      PrintReport();
   }
};

class TDataFrameImpl : public std::enable_shared_from_this<TDataFrameImpl> {

   ROOT::Internal::ActionBaseVec_t fBookedActions;
   ROOT::Detail::FilterBaseVec_t fBookedFilters;
   ROOT::Detail::FilterBaseVec_t fBookedNamedFilters;
   std::map<std::string, TmpBranchBasePtr_t> fBookedBranches;
   std::vector<std::shared_ptr<bool>> fResProxyReadiness;
   ::TDirectory *fDirPtr{nullptr};
   TTree *fTree{nullptr};
   const BranchNames_t fDefaultBranches;
   const unsigned int fNSlots{0};
   bool fHasRunAtLeastOnce{false};

public:
   TDataFrameImpl(TTree *tree, const BranchNames_t &defaultBranches);
   TDataFrameImpl(const TDataFrameImpl &) = delete;
   ~TDataFrameImpl(){};
   void Run();
   void BuildAllReaderValues(TTreeReader &r, unsigned int slot);
   void CreateSlots(unsigned int nSlots);
   TDataFrameImpl *GetImplPtr();
   std::shared_ptr<TDataFrameImpl> GetSharedPtr() { return shared_from_this(); }
   const BranchNames_t &GetDefaultBranches() const;
   const BranchNames_t GetTmpBranches() const { return {}; };
   TTree *GetTree() const;
   const TDataFrameBranchBase &GetBookedBranch(const std::string &name) const;
   void *GetTmpBranchValue(const std::string &branch, unsigned int slot, Long64_t entry);
   ::TDirectory *GetDirectory() const;
   std::string GetTreeName() const;
   void Book(const ROOT::Internal::ActionBasePtr_t &actionPtr);
   void Book(const ROOT::Detail::FilterBasePtr_t &filterPtr);
   void Book(const ROOT::Detail::TmpBranchBasePtr_t &branchPtr);
   bool CheckFilters(int, unsigned int);
   unsigned int GetNSlots() const;
   template<typename T>
   Experimental::TActionResultProxy<T> MakeActionResultProxy(const std::shared_ptr<T> &r)
   {
      auto readiness = std::make_shared<bool>(false);
      const auto &df = shared_from_this();
      auto resPtr = Experimental::TActionResultProxy<T>::MakeActionResultProxy(r, readiness, df);
      fResProxyReadiness.emplace_back(readiness);
      return resPtr;
   }
   bool HasRunAtLeastOnce() const { return fHasRunAtLeastOnce; }
   void Report() const;
   /// End of recursive chain of calls, does nothing
   void PartialReport() const {}
   void SetTree(TTree *tree) {fTree = tree;}
};

} // end NS ROOT::Detail

} // end NS ROOT

// Functions and method implementations
namespace ROOT {

namespace Experimental {

template<typename T>
void Experimental::TActionResultProxy<T>::TriggerRun()
{
   auto df = fImplWeakPtr.lock();
   if (!df) {
      throw std::runtime_error("The main TDataFrame is not reachable: did it go out of scope?");
   }
   df->Run();
}

template<typename FILENAMESCOLL, typename std::enable_if<ROOT::Internal::TDFTraitsUtils::TIsContainer<FILENAMESCOLL>::fgValue, int>::type>
TDataFrame::TDataFrame(const std::string &treeName, const FILENAMESCOLL &filenamescoll, const BranchNames_t &defaultBranches)
   : TDataFrameInterface<ROOT::Detail::TDataFrameImpl>(std::make_shared<ROOT::Detail::TDataFrameImpl>(nullptr, defaultBranches))
{
   auto chain = new TChain(treeName.c_str());
   for (auto &fileName : filenamescoll)
      chain->Add(ROOT::Internal::ToConstCharPtr(fileName));
   fTree = std::make_shared<TTree>(static_cast<TTree*>(chain));
   fProxiedPtr->SetTree(chain);
}

template<>
inline const char* TDataFrameInterface<ROOT::Detail::TDataFrameFilterBase>::GetNodeTypeName() { return "ROOT::Experimental::TDataFrameInterface<ROOT::Detail::TDataFrameFilterBase>";}
template<>
inline const char* TDataFrameInterface<ROOT::Detail::TDataFrameBranchBase>::GetNodeTypeName() { return "ROOT::Experimental::TDataFrameInterface<ROOT::Detail::TDataFrameBranchBase>";}
template<>
inline const char* TDataFrameInterface<ROOT::Detail::TDataFrameImpl>::GetNodeTypeName() { return "ROOT::Experimental::TDataFrameInterface<ROOT::Detail::TDataFrameImpl>";}

// Before we had to specialise the GetNodeTypeName method
extern template class TDataFrameInterface<ROOT::Detail::TDataFrameFilterBase>;
extern template class TDataFrameInterface<ROOT::Detail::TDataFrameBranchBase>;

} // end NS Experimental

namespace Internal {
template <typename T>
T &GetBranchValue(TVBPtr_t &readerValue, unsigned int slot, Long64_t entry, const std::string &branch,
                  ROOT::Detail::TDataFrameImpl *df, TDFTraitsUtils::TTypeList<T>)
{
   if (readerValue == nullptr) {
      // temporary branch
      void *tmpBranchVal = df->GetTmpBranchValue(branch, slot, entry);
      return *static_cast<T *>(tmpBranchVal);
   } else {
      // real branch
      return **std::static_pointer_cast<TTreeReaderValue<T>>(readerValue);
   }
}

template<typename T>
std::array_view<T> GetBranchValue(TVBPtr_t &readerValue, unsigned int slot,
                                  Long64_t entry, const std::string &branch,
                                  ROOT::Detail::TDataFrameImpl *df,
                                  TDFTraitsUtils::TTypeList<std::array_view<T>>)
{
   if(readerValue == nullptr) {
      // temporary branch
      void *tmpBranchVal = df->GetTmpBranchValue(branch, slot, entry);
      auto &tra = *static_cast<TTreeReaderArray<T> *>(tmpBranchVal);
      return std::array_view<T>(tra.begin(), tra.end());
   } else {
      // real branch
      auto &tra = *std::static_pointer_cast<TTreeReaderArray<T>>(readerValue);
      if (tra.GetSize() > 1 &&
          1 != (&tra[1] - &tra[0])) {
         std::string exceptionText = "Branch ";
         exceptionText += branch;
         exceptionText += " hangs from a non-split branch. For this reason, it cannot be accessed via an array_view. Please read the top level branch instead.";
         throw std::runtime_error(exceptionText.c_str());
      }
      return std::array_view<T>(tra.begin(), tra.end());
   }
}

} // end NS Internal

} // end NS ROOT

////////////////////////////////////////////////////////////////////////////////
/// Print a TDataFrame at the prompt:
namespace cling {
inline std::string printValue(ROOT::Experimental::TDataFrame *tdf)
{
   auto df = tdf->GetDataFrameChecked();
   auto treeName = df->GetTreeName();
   auto defBranches = df->GetDefaultBranches();
   auto tmpBranches = df->GetTmpBranches();

   std::ostringstream ret;
   ret << "A data frame built on top of the " << treeName << " dataset.";
   if (!defBranches.empty()) {
      if(defBranches.size() == 1) ret << "\nDefault branch: " << defBranches[0];
      else {
         ret << "\nDefault branches:\n";
         for (auto &&branch : defBranches) {
            ret << " - " << branch << "\n";
         }
      }
   }

   return ret.str();
}
}
#endif // ROOT_TDATAFRAME
