/**
 \file ROOT/RDF/RMergeableValue.hxx
 \ingroup dataframe
 \author Vincenzo Eduardo Padulano
 \author Enrico Guiraud
 \date 2020-06
*/

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RDF_RMERGEABLEVALUE
#define ROOT_RDF_RMERGEABLEVALUE

#include <memory>
#include <stdexcept>
#include <algorithm> // std::min, std::max

#include "RtypesCore.h"
#include "TList.h" // RMergeableFill::Merge

namespace ROOT {
namespace Detail {
namespace RDF {

// Fwd declarations for RMergeableValue
template <typename T>
class RMergeableValue;

template <typename T, typename... Ts>
std::unique_ptr<RMergeableValue<T>> MergeValues(std::unique_ptr<RMergeableValue<T>> OutputMergeable,
                                                std::unique_ptr<RMergeableValue<Ts>>... InputMergeables);

template <typename T, typename... Ts>
void MergeValues(RMergeableValue<T> &OutputMergeable, const RMergeableValue<Ts> &... InputMergeables);

/**
\class ROOT::Detail::RDF::RMergeableValueBase
\brief Base class of RMergeableValue.
\ingroup dataframe
Base class of the mergeable RDataFrame results family of classes. Provides a
non-templated custom type to allow passing a `std::unique_ptr` to the mergeable
object along the call chain. This class is never used in the public API and has
no meaning for the final user.
*/
class RMergeableValueBase {
public:
   virtual ~RMergeableValueBase() = default;
   /**
      Default constructor. Needed to allow serialization of ROOT objects. See
      [TBufferFile::WriteObjectClass]
      (classTBufferFile.html#a209078a4cb58373b627390790bf0c9c1)
   */
   RMergeableValueBase() = default;
};

/**
\class ROOT::Detail::RDF::RMergeableValue
\ingroup dataframe
\brief A result of an RDataFrame execution, that knows how to merge with other
results of the same type.
\tparam T Type of the action result.

Results of the execution of an RDataFrame computation graph do not natively
know how to merge with other results of the same type. In a distributed
environment it is often needed to have a merging mechanism for partial results
coming from the execution of an analysis on different chunks of the same dataset
that has happened on different executors. In order to achieve this,
RMergeableValue stores the result of the RDataFrame action and has a `Merge`
method to allow the aggregation of information coming from another similar
result into the current.

A mergeable value can be retrieved from an RResultPtr through the
[GetMergeableValue]
(namespaceROOT_1_1Detail_1_1RDF.html#a8b3a9c7b416826acc952d78a56d14ecb) free
function and a sequence of mergeables can be merged together with the helper
function [MergeValues]
(namespaceROOT_1_1Detail_1_1RDF.html#af16fefbe2d120983123ddf8a1e137277).
All the classes and functions involved are inside the `ROOT::Detail::RDF`
namespace.

In a nutshell:
~~~{.cpp}
using namespace ROOT::Detail::RDF;
ROOT::RDataFrame d("myTree", "file_*.root");
auto h1 = d.Histo1D("Branch_A");
auto h2 = d.Histo1D("Branch_A");

// Retrieve mergeables from the `RResultPtr`s
auto mergeableh1 = GetMergeableValue(h1);
auto mergeableh2 = GetMergeableValue(h2);

// Merge the values and get another mergeable back
auto mergedptr = MergeValues(std::move(mergeableh1), std::move(mergeableh2));

// Retrieve the merged TH1D object
const auto &mergedhisto = mergedptr->GetValue();
~~~

Though this snippet can run on a single thread of a single machine, it is
straightforward to generalize it to a distributed case, e.g. where `mergeableh1`
and `mergeableh2` are created on separate machines and sent to a `reduce`
process where the `MergeValues` function is called. The final user would then
just be given the final merged result coming from `mergedptr->GetValue`.

RMergeableValue is the base class for all the different specializations that may
be needed according to the peculiarities of the result types. The following
subclasses, their names hinting at the action operation of the result, are
currently available:

- RMergeableCount
- RMergeableFill, responsible for the following actions:
   - Graph
   - Histo{1,2,3}D
   - Profile{1,2}D
   - Stats
- RMergeableMax
- RMergeableMean
- RMergeableMin
- RMergeableStdDev
- RMergeableSum
*/
template <typename T>
class RMergeableValue : public RMergeableValueBase {
   // Friend function declarations
   template <typename T1, typename... Ts>
   friend std::unique_ptr<RMergeableValue<T1>> MergeValues(std::unique_ptr<RMergeableValue<T1>> OutputMergeable,
                                                           std::unique_ptr<RMergeableValue<Ts>>... InputMergeables);
   template <typename T1, typename... Ts>
   friend void MergeValues(RMergeableValue<T1> &OutputMergeable, const RMergeableValue<Ts> &... InputMergeables);

   /////////////////////////////////////////////////////////////////////////////
   /// \brief Aggregate the information contained in another RMergeableValue
   ///        into this.
   ///
   /// Virtual function reimplemented in all the subclasses.
   ///
   /// \note All the `Merge` methods in the RMergeableValue family are private.
   /// To merge multiple RMergeableValue objects please use [MergeValues]
   /// (namespaceROOT_1_1Detail_1_1RDF.html#af16fefbe2d120983123ddf8a1e137277).
   virtual void Merge(const RMergeableValue<T> &) = 0;

protected:
   T fValue;

public:
   /**
      Constructor taking the action result by const reference. This involves a
      copy of the result into the data member, but gives full ownership of data
      to the mergeable.
   */
   RMergeableValue(const T &value) : fValue{value} {}
   /**
      Default constructor. Needed to allow serialization of ROOT objects. See
      [TBufferFile::WriteObjectClass]
      (classTBufferFile.html#a209078a4cb58373b627390790bf0c9c1)
   */
   RMergeableValue() = default;
   /////////////////////////////////////////////////////////////////////////////
   /// \brief Retrieve the result wrapped by this mergeable.
   const T &GetValue() const { return fValue; }
};

/**
\class ROOT::Detail::RDF::RMergeableCount
\ingroup dataframe
\brief Specialization of RMergeableValue for the
[Count](classROOT_1_1RDF_1_1RInterface.html#a9678150c9c18cddd7b599690ba854734)
action.
*/
class RMergeableCount final : public RMergeableValue<ULong64_t> {
   /////////////////////////////////////////////////////////////////////////////
   /// \brief Aggregate the information contained in another RMergeableValue
   ///        into this.
   /// \param[in] other Another RMergeableValue object.
   /// \throws std::invalid_argument If the cast of the other object to the same
   ///         type as this one fails.
   ///
   /// The other RMergeableValue object is cast to the same type as this object.
   /// This is needed to make sure that only results of the same type of action
   /// are merged together. Then the two results are added together to update
   /// the value held by the current object.
   ///
   /// \note All the `Merge` methods in the RMergeableValue family are private.
   /// To merge multiple RMergeableValue objects please use [MergeValues]
   /// (namespaceROOT_1_1Detail_1_1RDF.html#af16fefbe2d120983123ddf8a1e137277).
   void Merge(const RMergeableValue<ULong64_t> &other) final
   {
      try {
         const auto &othercast = dynamic_cast<const RMergeableCount &>(other);
         this->fValue += othercast.fValue;
      } catch (const std::bad_cast &) {
         throw std::invalid_argument("Results from different actions cannot be merged together.");
      }
   }

public:
   /////////////////////////////////////////////////////////////////////////////
   /// \brief Constructor that initializes data members.
   /// \param[in] value The action result.
   RMergeableCount(ULong64_t value) : RMergeableValue<ULong64_t>(value) {}
   /**
      Default constructor. Needed to allow serialization of ROOT objects. See
      [TBufferFile::WriteObjectClass]
      (classTBufferFile.html#a209078a4cb58373b627390790bf0c9c1)
   */
   RMergeableCount() = default;
   RMergeableCount(RMergeableCount &&) = default;
   RMergeableCount(const RMergeableCount &) = delete;
};

/**
\class ROOT::Detail::RDF::RMergeableFill
\ingroup dataframe
\brief Specialization of RMergeableValue for histograms and statistics.

This subclass is responsible for merging results coming from the following
actions:
- [Graph](classROOT_1_1RDF_1_1RInterface.html#a804b466ebdbddef5c7e3400cc6b89301)
- [Histo{1D,2D,3D}]
  (classROOT_1_1RDF_1_1RInterface.html#a247ca3aeb7ce5b95015b7fae72983055)
- [Profile{1D,2D}]
  (classROOT_1_1RDF_1_1RInterface.html#a8ef7dc16b0e9f7bc9cfbe2d9e5de0cef)
- [Stats](classROOT_1_1RDF_1_1RInterface.html#abc68922c464e472f5f856e8981955af6)

*/
template <typename T>
class RMergeableFill final : public RMergeableValue<T> {
   /////////////////////////////////////////////////////////////////////////////
   /// \brief Aggregate the information contained in another RMergeableValue
   ///        into this.
   /// \param[in] other Another RMergeableValue object.
   /// \throws std::invalid_argument If the cast of the other object to the same
   ///         type as this one fails.
   ///
   /// The other RMergeableValue object is cast to the same type as this object.
   /// This is needed to make sure that only results of the same type of action
   /// are merged together. The function then calls the right `Merge` method
   /// according to the class of the fValue data member.
   ///
   /// \note All the `Merge` methods in the RMergeableValue family are private.
   /// To merge multiple RMergeableValue objects please use [MergeValues]
   /// (namespaceROOT_1_1Detail_1_1RDF.html#af16fefbe2d120983123ddf8a1e137277).
   void Merge(const RMergeableValue<T> &other) final
   {
      try {
         const auto &othercast = dynamic_cast<const RMergeableFill<T> &>(other);
         TList l;                                     // The `Merge` method accepts a TList
         l.Add(const_cast<T *>(&(othercast.fValue))); // Ugly but needed because of the signature of TList::Add
         this->fValue.Merge(&l); // if `T == TH1D` Eventually calls TH1::ExtendAxis that creates new instances of TH1D
      } catch (const std::bad_cast &) {
         throw std::invalid_argument("Results from different actions cannot be merged together.");
      }
   }

public:
   /////////////////////////////////////////////////////////////////////////////
   /// \brief Constructor that initializes data members.
   /// \param[in] value The action result.
   RMergeableFill(const T &value) : RMergeableValue<T>(value) {}
   /**
      Default constructor. Needed to allow serialization of ROOT objects. See
      [TBufferFile::WriteObjectClass]
      (classTBufferFile.html#a209078a4cb58373b627390790bf0c9c1)
   */
   RMergeableFill() = default;
   RMergeableFill(RMergeableFill &&) = default;
   RMergeableFill(const RMergeableFill &) = delete;
};

template <typename T>
class RMergeableMax final : public RMergeableValue<T> {

   void Merge(const RMergeableValue<T> &other) final
   {
      try {
         const auto &othercast = dynamic_cast<const RMergeableMax<T> &>(other);
         this->fValue = std::max(this->fValue, othercast.fValue);
      } catch (const std::bad_cast &) {
         throw std::invalid_argument("Results from different actions cannot be merged together.");
      }
   }

public:
   /////////////////////////////////////////////////////////////////////////////
   /// \brief Constructor that initializes data members.
   /// \param[in] value The action result.
   RMergeableMax(const T &value) : RMergeableValue<T>(value) {}
   /**
      Default constructor. Needed to allow serialization of ROOT objects. See
      [TBufferFile::WriteObjectClass]
      (classTBufferFile.html#a209078a4cb58373b627390790bf0c9c1)
   */
   RMergeableMax() = default;
   RMergeableMax(RMergeableMax &&) = default;
   RMergeableMax(const RMergeableMax &) = delete;
};

/**
\class ROOT::Detail::RDF::RMergeableMean
\ingroup dataframe
\brief Specialization of RMergeableValue for the
[Mean](classROOT_1_1RDF_1_1RInterface.html#ade6b020284f2f4fe9d3b09246b5f376a)
action.

This subclass is responsible for merging results coming from Mean actions. Other
than the result itself, the number of entries that were used to compute that
mean is also stored in the object.
*/
class RMergeableMean final : public RMergeableValue<Double_t> {
   ULong64_t fCounts; ///< The number of entries used to compute the mean.

   /////////////////////////////////////////////////////////////////////////////
   /// \brief Aggregate the information contained in another RMergeableValue
   ///        into this.
   /// \param[in] other Another RMergeableValue object.
   /// \throws std::invalid_argument If the cast of the other object to the same
   ///         type as this one fails.
   ///
   /// The other RMergeableValue object is cast to the same type as this object.
   /// This is needed to make sure that only results of the same type of action
   /// are merged together. The function then computes the weighted mean of the
   /// two means held by the mergeables.
   ///
   /// \note All the `Merge` methods in the RMergeableValue family are private.
   /// To merge multiple RMergeableValue objects please use [MergeValues]
   /// (namespaceROOT_1_1Detail_1_1RDF.html#af16fefbe2d120983123ddf8a1e137277).
   void Merge(const RMergeableValue<Double_t> &other) final
   {
      try {
         const auto &othercast = dynamic_cast<const RMergeableMean &>(other);
         const auto &othervalue = othercast.fValue;
         const auto &othercounts = othercast.fCounts;

         // Compute numerator and denumerator of the weighted mean
         const auto num = this->fValue * fCounts + othervalue * othercounts;
         const auto denum = static_cast<Double_t>(fCounts + othercounts);

         // Update data members
         this->fValue = num / denum;
         fCounts += othercounts;
      } catch (const std::bad_cast &) {
         throw std::invalid_argument("Results from different actions cannot be merged together.");
      }
   }

public:
   /////////////////////////////////////////////////////////////////////////////
   /// \brief Constructor that initializes data members.
   /// \param[in] value The action result.
   /// \param[in] counts The number of entries used to compute that result.
   RMergeableMean(Double_t value, ULong64_t counts) : RMergeableValue<Double_t>(value), fCounts{counts} {}
   /**
      Default constructor. Needed to allow serialization of ROOT objects. See
      [TBufferFile::WriteObjectClass]
      (classTBufferFile.html#a209078a4cb58373b627390790bf0c9c1)
   */
   RMergeableMean() = default;
   RMergeableMean(RMergeableMean &&) = default;
   RMergeableMean(const RMergeableMean &) = delete;
};

template <typename T>
class RMergeableMin final : public RMergeableValue<T> {

   void Merge(const RMergeableValue<T> &other) final
   {
      try {
         const auto &othercast = dynamic_cast<const RMergeableMin<T> &>(other);
         this->fValue = std::min(this->fValue, othercast.fValue);
      } catch (const std::bad_cast &) {
         throw std::invalid_argument("Results from different actions cannot be merged together.");
      }
   }

public:
   /////////////////////////////////////////////////////////////////////////////
   /// \brief Constructor that initializes data members.
   /// \param[in] value The action result.
   RMergeableMin(const T &value) : RMergeableValue<T>(value) {}
   /**
      Default constructor. Needed to allow serialization of ROOT objects. See
      [TBufferFile::WriteObjectClass]
      (classTBufferFile.html#a209078a4cb58373b627390790bf0c9c1)
   */
   RMergeableMin() = default;
   RMergeableMin(RMergeableMin &&) = default;
   RMergeableMin(const RMergeableMin &) = delete;
};

/**
\class ROOT::Detail::RDF::RMergeableStdDev
\ingroup dataframe
\brief Specialization of RMergeableValue for the
[StdDev](classROOT_1_1RDF_1_1RInterface.html#a482c4e4f81fe1e421c016f89cd281572)
action.

This class also stores information about the number of entries and the average
used to compute the standard deviation.
*/
class RMergeableStdDev final : public RMergeableValue<Double_t> {
   ULong64_t fCounts; ///< Number of entries of the set.
   Double_t fMean;    ///< Average of the set.

   /////////////////////////////////////////////////////////////////////////////
   /// \brief Aggregate the information contained in another RMergeableValue
   ///        into this.
   /// \param[in] other Another RMergeableValue object.
   /// \throws std::invalid_argument If the cast of the other object to the same
   ///         type as this one fails.
   ///
   /// The other RMergeableValue object is cast to the same type as this object.
   /// This is needed to make sure that only results of the same type of action
   /// are merged together. The function then computes the aggregated standard
   /// deviation of the two samples using an algorithm by
   /// [Chan et al. (1979)]
   /// (http://i.stanford.edu/pub/cstr/reports/cs/tr/79/773/CS-TR-79-773.pdf)
   ///
   /// \note All the `Merge` methods in the RMergeableValue family are private.
   /// To merge multiple RMergeableValue objects please use [MergeValues]
   /// (namespaceROOT_1_1Detail_1_1RDF.html#af16fefbe2d120983123ddf8a1e137277).
   void Merge(const RMergeableValue<Double_t> &other) final
   {
      try {
         const auto &othercast = dynamic_cast<const RMergeableStdDev &>(other);
         const auto &othercounts = othercast.fCounts;
         const auto &othermean = othercast.fMean;

         // Compute the aggregated variance using an algorithm by Chan et al.
         // See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Parallel_algorithm
         const auto thisvariance = std::pow(this->fValue, 2);
         const auto othervariance = std::pow(othercast.fValue, 2);

         const auto delta = othermean - fMean;

         const auto m_a = thisvariance * (fCounts - 1);
         const auto m_b = othervariance * (othercounts - 1);

         const auto sumcounts = static_cast<Double_t>(fCounts + othercounts);

         const auto M2 = m_a + m_b + std::pow(delta, 2) * fCounts * othercounts / sumcounts;

         const auto meannum = fMean * fCounts + othermean * othercounts;

         // Update the data members
         this->fValue = std::sqrt(M2 / (sumcounts - 1));
         fMean = meannum / sumcounts;
         fCounts += othercounts;
      } catch (const std::bad_cast &) {
         throw std::invalid_argument("Results from different actions cannot be merged together.");
      }
   }

public:
   /////////////////////////////////////////////////////////////////////////////
   /// \brief Constructor that initializes data members.
   /// \param[in] value The action result.
   /// \param[in] counts The number of entries of the set.
   /// \param[in] mean The average of the set.
   RMergeableStdDev(Double_t value, ULong64_t counts, Double_t mean)
      : RMergeableValue<Double_t>(value), fCounts{counts}, fMean{mean}
   {
   }
   /**
      Default constructor. Needed to allow serialization of ROOT objects. See
      [TBufferFile::WriteObjectClass]
      (classTBufferFile.html#a209078a4cb58373b627390790bf0c9c1)
   */
   RMergeableStdDev() = default;
   RMergeableStdDev(RMergeableStdDev &&) = default;
   RMergeableStdDev(const RMergeableStdDev &) = delete;
};

template <typename T>
class RMergeableSum final : public RMergeableValue<T> {

   void Merge(const RMergeableValue<T> &other) final
   {
      try {
         const auto &othercast = dynamic_cast<const RMergeableSum<T> &>(other);
         this->fValue += othercast.fValue;
      } catch (const std::bad_cast &) {
         throw std::invalid_argument("Results from different actions cannot be merged together.");
      }
   }

public:
   /////////////////////////////////////////////////////////////////////////////
   /// \brief Constructor that initializes data members.
   /// \param[in] value The action result.
   RMergeableSum(const T &value) : RMergeableValue<T>(value) {}
   /**
      Default constructor. Needed to allow serialization of ROOT objects. See
      [TBufferFile::WriteObjectClass]
      (classTBufferFile.html#a209078a4cb58373b627390790bf0c9c1)
   */
   RMergeableSum() = default;
   RMergeableSum(RMergeableSum &&) = default;
   RMergeableSum(const RMergeableSum &) = delete;
};

/// \cond HIDDEN_SYMBOLS
// What follows mimics C++17 std::conjunction without using recursive template instantiations.
// Used in `MergeValues` to check that all the mergeables hold values of the same type.
template <bool...>
struct bool_pack {
};
template <class... Ts>
using conjunction = std::is_same<bool_pack<true, Ts::value...>, bool_pack<Ts::value..., true>>;
/// \endcond

////////////////////////////////////////////////////////////////////////////////
/// \brief Merge multiple RMergeableValue objects into one.
/// \param[in] OutputMergeable The mergeable object where all the information
///            will be aggregated.
/// \param[in] InputMergeables Other mergeables containing the partial results.
/// \returns An RMergeableValue holding the aggregated value wrapped in an
///          `std::unique_ptr`.
///
/// This is the recommended way of merging multiple RMergeableValue objects.
/// This overload takes ownership of the mergeables and gives back to the user
/// a mergeable with the aggregated information. All the mergeables with the
/// partial results get destroyed in the process.
///
/// Example usage:
/// ~~~{.cpp}
/// using namespace ROOT::Detail::RDF;
/// // mh1, mh2, mh3 are std::unique_ptr<RMergeableValue<TH1D>>
/// auto mergedptr = MergeValues(std::move(mh1), std::move(mh2), std::move(mh3));
/// const auto &mergedhisto = mergedptr->GetValue(); // Final merged histogram
/// // Do stuff with it
/// mergedhisto.Draw();
/// ~~~
template <typename T, typename... Ts>
std::unique_ptr<RMergeableValue<T>> MergeValues(std::unique_ptr<RMergeableValue<T>> OutputMergeable,
                                                std::unique_ptr<RMergeableValue<Ts>>... InputMergeables)
{
   // Check all mergeables have the same template type
   static_assert(conjunction<std::is_same<Ts, T>...>::value, "Values must all be of the same type.");

   // Using dummy array initialization inspired by https://stackoverflow.com/a/25683817
   using expander = int[];
   // Cast to void to suppress unused-value warning in Clang
   (void)expander{0, (OutputMergeable->Merge(*InputMergeables), 0)...};

   return OutputMergeable;
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Merge multiple RMergeableValue objects into one.
/// \param[in,out] OutputMergeable The mergeable object where all the
///                information will be aggregated.
/// \param[in] InputMergeables Other mergeables containing the partial results.
///
/// This overload modifies the mergeable objects in-place. The ownership is left
/// to the caller. The first argument to the function will get all the
/// values contained in the other arguments merged into itself. This is a
/// convenience overload introduced for the ROOT Python API.
///
/// Example usage:
/// ~~~{.cpp}
/// // mh1, mh2, mh3 are std::unique_ptr<RMergeableValue<TH1D>>
/// ROOT::Detail::RDF::MergeValues(*mh1, *mh2, *mh3);
/// const auto &mergedhisto = mh1->GetValue(); // Final merged histogram
/// // Do stuff with it
/// mergedhisto.Draw();
/// ~~~
template <typename T, typename... Ts>
void MergeValues(RMergeableValue<T> &OutputMergeable, const RMergeableValue<Ts> &... InputMergeables)
{
   // Check all mergeables are of the same type
   static_assert(conjunction<std::is_same<Ts, T>...>::value, "Values must all be of the same type.");

   // Using dummy array initialization inspired by https://stackoverflow.com/a/25683817
   using expander = int[];
   // Cast to void to suppress unused-value warning in Clang
   (void)expander{0, (OutputMergeable.Merge(InputMergeables), 0)...};
}
} // namespace RDF
} // namespace Detail
} // namespace ROOT

#endif // ROOT_RDF_RMERGEABLEVALUE