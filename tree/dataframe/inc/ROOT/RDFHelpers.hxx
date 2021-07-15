// Author: Enrico Guiraud, Danilo Piparo CERN  02/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// This header contains helper free functions that slim down RDataFrame's programming model

#ifndef ROOT_RDF_HELPERS
#define ROOT_RDF_HELPERS

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultHandle.hxx>
#include <ROOT/RDF/GraphUtils.hxx>
#include <ROOT/TypeTraits.hxx>

#include <algorithm> // std::transform
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <type_traits>
#include <utility> // std::index_sequence
#include <vector>
#include <mutex>
#include <chrono>
#include <array>

namespace ROOT {
namespace Internal {
namespace RDF {
template <typename... ArgTypes, typename F>
std::function<bool(ArgTypes...)> NotHelper(ROOT::TypeTraits::TypeList<ArgTypes...>, F &&f)
{
   return std::function<bool(ArgTypes...)>([=](ArgTypes... args) mutable { return !f(args...); });
}

template <typename... ArgTypes, typename Ret, typename... Args>
std::function<bool(ArgTypes...)> NotHelper(ROOT::TypeTraits::TypeList<ArgTypes...>, Ret (*f)(Args...))
{
   return std::function<bool(ArgTypes...)>([=](ArgTypes... args) mutable { return !f(args...); });
}

template <typename I, typename T, typename F>
class PassAsVecHelper;

template <std::size_t... N, typename T, typename F>
class PassAsVecHelper<std::index_sequence<N...>, T, F> {
   template <std::size_t Idx>
   using AlwaysT = T;
   std::decay_t<F> fFunc;

public:
   PassAsVecHelper(F &&f) : fFunc(std::forward<F>(f)) {}
   auto operator()(AlwaysT<N>... args) -> decltype(fFunc({args...})) { return fFunc({args...}); }
};

template <std::size_t N, typename T, typename F>
auto PassAsVec(F &&f) -> PassAsVecHelper<std::make_index_sequence<N>, T, F>
{
   return PassAsVecHelper<std::make_index_sequence<N>, T, F>(std::forward<F>(f));
}

} // namespace RDF
} // namespace Internal

namespace RDF {
namespace RDFInternal = ROOT::Internal::RDF;


// clag-format off
/// Given a callable with signature bool(T1, T2, ...) return a callable with same signature that returns the negated result
///
/// The callable must have one single non-template definition of operator(). This is a limitation with respect to
/// std::not_fn, required for interoperability with RDataFrame.
// clang-format on
template <typename F,
          typename Args = typename ROOT::TypeTraits::CallableTraits<std::decay_t<F>>::arg_types_nodecay,
          typename Ret = typename ROOT::TypeTraits::CallableTraits<std::decay_t<F>>::ret_type>
auto Not(F &&f) -> decltype(RDFInternal::NotHelper(Args(), std::forward<F>(f)))
{
   static_assert(std::is_same<Ret, bool>::value, "RDF::Not requires a callable that returns a bool.");
   return RDFInternal::NotHelper(Args(), std::forward<F>(f));
}

// clang-format off
/// PassAsVec is a callable generator that allows passing N variables of type T to a function as a single collection.
///
/// PassAsVec<N, T>(func) returns a callable that takes N arguments of type T, passes them down to function `func` as
/// an initializer list `{t1, t2, t3,..., tN}` and returns whatever f({t1, t2, t3, ..., tN}) returns.
///
/// Note that for this to work with RDataFrame the type of all columns that the callable is applied to must be exactly T.
/// Example usage together with RDataFrame ("varX" columns must all be `float` variables):
/// \code
/// bool myVecFunc(std::vector<float> args);
/// df.Filter(PassAsVec<3, float>(myVecFunc), {"var1", "var2", "var3"});
/// \endcode
// clang-format on
template <std::size_t N, typename T, typename F>
auto PassAsVec(F &&f) -> RDFInternal::PassAsVecHelper<std::make_index_sequence<N>, T, F>
{
   return RDFInternal::PassAsVecHelper<std::make_index_sequence<N>, T, F>(std::forward<F>(f));
}

// clang-format off
/// Create a graphviz representation of the dataframe computation graph, return it as a string.
/// \param[in] node any node of the graph. Called on the head (first) node, it prints the entire graph. Otherwise, only the branch the node belongs to.
///
/// The output can be displayed with a command akin to `dot -Tpng output.dot > output.png && open output.png`.
///
/// Note that "hanging" Defines, i.e. Defines without downstream nodes, will not be displayed by SaveGraph as they are
/// effectively optimized away from the computation graph.
///
/// Note that SaveGraph is not thread-safe and must not be called concurrently from different threads.
// clang-format on
template <typename NodeType>
std::string SaveGraph(NodeType node)
{
   ROOT::Internal::RDF::GraphDrawing::GraphCreatorHelper helper;
   return helper(node);
}

// clang-format off
/// Create a graphviz representation of the dataframe computation graph, write it to the specified file.
/// \param[in] node any node of the graph. Called on the head (first) node, it prints the entire graph. Otherwise, only the branch the node belongs to.
/// \param[in] outputFile file where to save the representation.
///
/// The output can be displayed with a command akin to `dot -Tpng output.dot > output.png && open output.png`.
///
/// Note that "hanging" Defines, i.e. Defines without downstream nodes, will not be displayed by SaveGraph as they are
/// effectively optimized away from the computation graph.
///
/// Note that SaveGraph is not thread-safe and must not be called concurrently from different threads.
// clang-format on
template <typename NodeType>
void SaveGraph(NodeType node, const std::string &outputFile)
{
   ROOT::Internal::RDF::GraphDrawing::GraphCreatorHelper helper;
   std::string dotGraph = helper(node);

   std::ofstream out(outputFile);
   if (!out.is_open()) {
      throw std::runtime_error("Could not open output file \"" + outputFile  + "\"for reading");
   }

   out << dotGraph;
   out.close();
}

// clang-format off
/// Cast a RDataFrame node to the common type ROOT::RDF::RNode
/// \param[in] node Any node of a RDataFrame graph
// clang-format on
template <typename NodeType>
RNode AsRNode(NodeType node)
{
   return node;
}

// clang-format off
/// Trigger the event loop of multiple RDataFrames concurrently
/// \param[in] handles A vector of RResultHandles
///
/// This function triggers the event loop of all computation graphs which relate to the
/// given RResultHandles. The advantage compared to running the event loop implicitly by accessing the
/// RResultPtr is that the event loops will run concurrently. Therefore, the overall
/// computation of all results is generally more efficient.
/// It should be noted that user-defined operations (e.g., Filters and Defines) of the different RDataFrame graphs are assumed to be safe to call concurrently.
///
/// ~~~{.cpp}
/// ROOT::RDataFrame df1("tree1", "file1.root");
/// auto r1 = df1.Histo1D("var1");
///
/// ROOT::RDataFrame df2("tree2", "file2.root");
/// auto r2 = df2.Sum("var2");
///
/// // RResultPtr -> RResultHandle conversion is automatic
/// ROOT::RDF::RunGraphs({r1, r2});
/// ~~~
// clang-format on
void RunGraphs(std::vector<RResultHandle> handles);


// clang-format off
/// RDF progress helper.
/// This class provides callback functions to RDataFrame, which record event statistics
/// and print them to the terminal every second.
/// ProgressHelper::operator()(unsigned int, T&) is thread safe, and can be used as callback in MT mode.
// clang-format on
class ProgressHelper {
public:

  /// Create a progress helper.
  /// \param increment RDF callbacks are called every `n` events. Pass this `n` here.
  /// \param maxEvent If known, pass the number of events in the dataset. This will allow for computing completion
  /// statistics.
  /// \param progressBarWidth Number of characters the progress bar will occupy.
  /// \param useShellColours Use shell colour codes to colour the output.
  ProgressHelper(std::size_t increment,
      std::size_t maxEvent,
      unsigned int progressBarWidth = 40,
      bool useShellColours = true) :
  fMaxEvents{maxEvent},
  fIncrement{increment},
  fBarWidth{progressBarWidth},
  fUseShellColours{useShellColours}
  { }

  ~ProgressHelper() = default;

  // clang-format off
  /// Thread-safe callback for RDataFrame.
  /// It will record elapsed times and event statistics, and print a progress bar every second.
  /// \param slot Ignored.
  /// \param value Ignored.
  // clang-format on
  template<typename T>
  void operator()(unsigned int /*slot*/, T& value) { operator()(value); }

  // clang-format off
  /// Thread-safe callback for RDataFrame.
  /// It will record elapsed times and event statistics, and print a progress bar every second.
  /// \param value Ignored.
  // clang-format on
  template<typename T>
  void operator()(T& /*value*/) {
    using namespace std::chrono;
    // ***************************************************
    // Warning: Here, everything needs to be thread safe:
    // ***************************************************
    fProcessedEvents += fIncrement;

    // We only print every 1s, and only one thread does it. Try to take the lock:
    if (duration_cast<seconds>(system_clock::now() - fLastPrintTime).count() == 0 || !fPrintMutex.try_lock())
      return;

    // ***************************************************
    // Protected by lock from here:
    // ***************************************************
    std::lock_guard<std::mutex> lockGuard(fPrintMutex, std::adopt_lock);

    std::size_t eventCount;
    seconds elapsedSeconds;
    std::tie(eventCount, elapsedSeconds) = RecordEvtCountAndTime();

    std::cout << "\r";
    PrintProgressbar(std::cout, eventCount);
    PrintStats(std::cout, eventCount, elapsedSeconds);
    std::cout << std::flush;
  }


private:
  double EvtPerSec() const;
  std::pair<std::size_t, std::chrono::seconds> RecordEvtCountAndTime();
  void PrintStats(std::ostream& stream, std::size_t currentEventCount, std::chrono::seconds totalElapsedSeconds) const;
  void PrintProgressbar(std::ostream& stream, std::size_t currentEventCount) const;

  const std::chrono::time_point<std::chrono::system_clock> fBeginTime = std::chrono::system_clock::now();
  std::chrono::time_point<std::chrono::system_clock> fLastPrintTime = fBeginTime;

  std::atomic<std::size_t> fProcessedEvents{ 0 };
  std::size_t fLastProcessedEvents = 0;
  const std::size_t fMaxEvents;
  const std::size_t fIncrement;

  std::array<double, 20> fEventsPerSecondStatistics;
  std::size_t fEventsPerSecondStatisticsIndex = 0;

  const unsigned int fBarWidth;

  std::mutex fPrintMutex;
  const bool fUseShellColours;
};


// clang-format off
/// Count events in a tree in the file described by `fileUrl`.
/// \note This will open and close the file, so don't use this in loops.
///
/// \param treename Name of the tree.
/// \param fileUrl Path or URL to a ROOT file.
/// \return Event counts. Returns 0 in case of errors.
// clang-format on
std::size_t CountEvents(const char* treename, const char* fileUrl);

} // namespace RDF
} // namespace ROOT
#endif
