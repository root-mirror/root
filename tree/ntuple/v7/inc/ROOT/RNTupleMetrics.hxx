/// \file ROOT/RNTupleMetrics.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2019-08-27
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RNTupleMetrics
#define ROOT7_RNTupleMetrics

#include <ROOT/RConfig.hxx>
#include <ROOT/RStringView.hxx>

#include <TError.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <ctime> // for CPU time measurement with clock()
#include <functional>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <vector>
#include <utility>

namespace ROOT {
namespace Experimental {
namespace Detail {

class RNTupleMetrics;

// clang-format off
/**
\class ROOT::Experimental::Detail::RNTuplePerfCounter
\ingroup NTuple
\brief A performance counter with a name and a unit, which can be activated on demand

Derived classes decide on the counter type and implement printing of the value.
*/
// clang-format on
class RNTuplePerfCounter {
private:
   /// Symbol to split name, unit, description, and value when printing
   static constexpr char kFieldSeperator = '|';

   std::string fName;
   std::string fUnit;
   std::string fDescription;
   bool fIsEnabled = false;

public:
   RNTuplePerfCounter(const std::string &name, const std::string &unit, const std::string &desc)
      : fName(name), fUnit(unit), fDescription(desc) {}
   virtual ~RNTuplePerfCounter();
   void Enable() { fIsEnabled = true; }
   bool IsEnabled() const { return fIsEnabled; }
   std::string GetName() const { return fName; }
   std::string GetDescription() const { return fDescription; }
   std::string GetUnit() const { return fUnit; }

   virtual std::int64_t GetValueAsInt() const = 0;
   virtual std::string GetValueAsString() const = 0;
   std::string ToString() const;
};


// clang-format off
/**
\class ROOT::Experimental::Detail::RNTuplePlainCounter
\ingroup NTuple
\brief A non thread-safe integral performance counter
*/
// clang-format on
class RNTuplePlainCounter : public RNTuplePerfCounter {
private:
   std::int64_t fCounter = 0;

public:
   RNTuplePlainCounter(const std::string &name, const std::string &unit, const std::string &desc)
      : RNTuplePerfCounter(name, unit, desc)
   {
   }

   R__ALWAYS_INLINE void Inc() { ++fCounter; }
   R__ALWAYS_INLINE void Dec() { --fCounter; }
   R__ALWAYS_INLINE void Add(int64_t delta) { fCounter += delta; }
   R__ALWAYS_INLINE int64_t GetValue() const { return fCounter; }
   R__ALWAYS_INLINE void SetValue(int64_t val) { fCounter = val; }
   std::int64_t GetValueAsInt() const override { return fCounter; }
   std::string GetValueAsString() const override { return std::to_string(fCounter); }
};


// clang-format off
/**
\class ROOT::Experimental::Detail::RNTupleAtomicCounter
\ingroup NTuple
\brief A thread-safe integral performance counter
*/
// clang-format on
class RNTupleAtomicCounter : public RNTuplePerfCounter {
private:
   std::atomic<std::int64_t> fCounter{0};

public:
   RNTupleAtomicCounter(const std::string &name, const std::string &unit, const std::string &desc)
      : RNTuplePerfCounter(name, unit, desc) { }

   R__ALWAYS_INLINE
   void Inc() {
      if (R__unlikely(IsEnabled()))
         ++fCounter;
   }

   R__ALWAYS_INLINE
   void Dec() {
      if (R__unlikely(IsEnabled()))
         --fCounter;
   }

   R__ALWAYS_INLINE
   void Add(int64_t delta) {
      if (R__unlikely(IsEnabled()))
         fCounter += delta;
   }

   R__ALWAYS_INLINE
   int64_t XAdd(int64_t delta) {
      if (R__unlikely(IsEnabled()))
         return fCounter.fetch_add(delta);
      return 0;
   }

   R__ALWAYS_INLINE
   int64_t GetValue() const {
      if (R__unlikely(IsEnabled()))
         return fCounter.load();
      return 0;
   }

   R__ALWAYS_INLINE
   void SetValue(int64_t val) {
      if (R__unlikely(IsEnabled()))
         fCounter.store(val);
   }

   std::int64_t GetValueAsInt() const override { return GetValue(); }
   std::string GetValueAsString() const override { return std::to_string(GetValue()); }
};


// clang-format off
/**
\class ROOT::Experimental::Detail::RNTupleCalcPerf
\ingroup NTuple
\brief A metric element that computes its floating point value from other counters.
*/
// clang-format on
class RNTupleCalcPerf : public RNTuplePerfCounter {
public:
   using MetricFunc_t = std::function<std::pair<bool, double>(const RNTupleMetrics &)>;

private:
   RNTupleMetrics &fMetrics;
   const MetricFunc_t fFunc;

public:
   RNTupleCalcPerf(const std::string &name, const std::string &unit, const std::string &desc,
                   RNTupleMetrics &metrics, MetricFunc_t &&func)
      : RNTuplePerfCounter(name, unit, desc), fMetrics(metrics), fFunc(std::move(func))
   {
   }

   double GetValue() const {
      auto ret = fFunc(fMetrics);
      if (ret.first)
         return ret.second;
      return std::numeric_limits<double>::quiet_NaN();
   }

   std::int64_t GetValueAsInt() const override {
      return static_cast<std::int64_t>(GetValue());
   }

   std::string GetValueAsString() const override {
      return std::to_string(GetValue());
   }
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RNTupleTickCounter
\ingroup NTuple
\brief An either thread-safe or non thread safe counter for CPU ticks

On print, the value is converted from ticks to ns.
*/
// clang-format on
template <typename BaseCounterT>
class RNTupleTickCounter : public BaseCounterT {
public:
   RNTupleTickCounter(const std::string &name, const std::string &unit, const std::string &desc)
      : BaseCounterT(name, unit, desc)
   {
      R__ASSERT(unit == "ns");
   }

   std::int64_t GetValueAsInt() const final {
      auto ticks = BaseCounterT::GetValue();
      return std::uint64_t((double(ticks) / double(CLOCKS_PER_SEC)) * (1000. * 1000. * 1000.));
   }

   std::string GetValueAsString() const final {
      return std::to_string(GetValueAsInt());
   }
};


// clang-format off
/**
\class ROOT::Experimental::Detail::RNTupleTimer
\ingroup NTuple
\brief Record wall time and CPU time between construction and destruction

Uses RAII as a stop watch. Only the wall time counter is used to determine whether the timer is active.
*/
// clang-format on
template <typename WallTimeT, typename CpuTimeT>
class RNTupleTimer {
private:
   using Clock_t = std::chrono::steady_clock;

   WallTimeT &fCtrWallTime;
   CpuTimeT &fCtrCpuTicks;
   /// Wall clock time
   Clock_t::time_point fStartTime;
   /// CPU time
   clock_t fStartTicks = 0;

public:
   RNTupleTimer(WallTimeT &ctrWallTime, CpuTimeT &ctrCpuTicks)
      : fCtrWallTime(ctrWallTime), fCtrCpuTicks(ctrCpuTicks)
   {
      if (!fCtrWallTime.IsEnabled())
         return;
      fStartTime = Clock_t::now();
      fStartTicks = clock();
   }

   ~RNTupleTimer() {
      if (!fCtrWallTime.IsEnabled())
         return;
      auto wallTimeNs = std::chrono::duration_cast<std::chrono::nanoseconds>(Clock_t::now() - fStartTime);
      fCtrWallTime.Add(wallTimeNs.count());
      fCtrCpuTicks.Add(clock() - fStartTicks);
   }

   RNTupleTimer(const RNTupleTimer &other) = delete;
   RNTupleTimer &operator =(const RNTupleTimer &other) = delete;
};

using RNTuplePlainTimer = RNTupleTimer<RNTuplePlainCounter, RNTupleTickCounter<RNTuplePlainCounter>>;
using RNTupleAtomicTimer = RNTupleTimer<RNTupleAtomicCounter, RNTupleTickCounter<RNTupleAtomicCounter>>;


// clang-format off
/**
\class ROOT::Experimental::Detail::RNTupleMetrics
\ingroup NTuple
\brief A collection of Counter objects with a name, a unit, and a description.

The class owns the counters; on registration of a new
*/
// clang-format on
class RNTupleMetrics {
private:
   /// Symbol to split metrics name from counter / sub metrics name
   static constexpr char kNamespaceSeperator = '.';

   std::vector<std::unique_ptr<RNTuplePerfCounter>> fCounters;
   std::vector<RNTupleMetrics *> fObservedMetrics;
   std::string fName;
   bool fIsEnabled = false;

   bool Contains(const std::string &name) const;

public:
   explicit RNTupleMetrics(const std::string &name) : fName(name) {}
   RNTupleMetrics(const RNTupleMetrics &other) = delete;
   RNTupleMetrics & operator=(const RNTupleMetrics &other) = delete;
   RNTupleMetrics(RNTupleMetrics &&other) = default;
   RNTupleMetrics & operator=(RNTupleMetrics &&other) = default;
   ~RNTupleMetrics() = default;

   // TODO(jblomer): return a reference
   template <typename CounterPtrT, class... Args>
   CounterPtrT MakeCounter(const std::string &name, Args&&... args)
   {
      R__ASSERT(!Contains(name));
      auto counter = std::make_unique<std::remove_pointer_t<CounterPtrT>>(name, std::forward<Args>(args)...);
      auto ptrCounter = counter.get();
      fCounters.emplace_back(std::move(counter));
      return ptrCounter;
   }

   /// Searches this object and all the observed sub metrics. Returns nullptr if name is not found.
   const RNTuplePerfCounter *GetCounter(std::string_view name) const;

   void ObserveMetrics(RNTupleMetrics &observee);

   void Print(std::ostream &output, const std::string &prefix = "") const;
   void Enable();
   bool IsEnabled() const { return fIsEnabled; }
};

} // namespace Detail
} // namespace Experimental
} // namespace ROOT

#endif
