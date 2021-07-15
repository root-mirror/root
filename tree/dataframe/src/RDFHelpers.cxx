// Author: Stefan Wunsch, Enrico Guiraud CERN  09/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RDFHelpers.hxx"
#include "TROOT.h"      // IsImplicitMTEnabled
#include "TError.h"     // Warning
#include "RConfigure.h" // R__USE_IMT
#include "TFile.h"
#include "TTree.h"
#ifdef R__USE_IMT
#include "ROOT/TThreadExecutor.hxx"
#endif // R__USE_IMT

#include <set>
#include <numeric>

using ROOT::RDF::RResultHandle;

void ROOT::RDF::RunGraphs(std::vector<RResultHandle> handles)
{
   if (handles.empty()) {
      Warning("RunGraphs", "Got an empty list of handles");
      return;
   }

   // Check that there are results which have not yet been run
   unsigned int nNotRun = 0;
   for (const auto &h : handles) {
      if (!h.IsReady())
         nNotRun++;
   }
   if (nNotRun < handles.size()) {
      Warning("RunGraphs", "Got %lu handles from which %lu link to results which are already ready.", handles.size(),
              handles.size() - nNotRun);
      return;
   }
   if (nNotRun == 0)
      return;

   // Find the unique event loops
   auto sameGraph = [](const RResultHandle &a, const RResultHandle &b) { return a.fLoopManager < b.fLoopManager; };
   std::set<RResultHandle, decltype(sameGraph)> s(handles.begin(), handles.end(), sameGraph);
   std::vector<RResultHandle> uniqueLoops(s.begin(), s.end());

   // Trigger the unique event loops
   auto run = [](RResultHandle &h) { h.fLoopManager->Run(); };
#ifdef R__USE_IMT
   if (ROOT::IsImplicitMTEnabled()) {
      ROOT::TThreadExecutor{}.Foreach(run, uniqueLoops);
      return;
   }
#endif // R__USE_IMT
   for (auto &h : uniqueLoops)
      run(h);
}


/// Compute a running mean of events/s.
double ROOT::RDF::ProgressHelper::EvtPerSec() const {
  if (fEventsPerSecondStatisticsIndex < fEventsPerSecondStatistics.size())
    return std::accumulate(fEventsPerSecondStatistics.begin(), fEventsPerSecondStatistics.begin() + fEventsPerSecondStatisticsIndex, 0.) / fEventsPerSecondStatisticsIndex;
  else
    return std::accumulate(fEventsPerSecondStatistics.begin(), fEventsPerSecondStatistics.end(), 0.) / fEventsPerSecondStatistics.size();
}

/// Record current event counts and time stamp, populate evts/s statistics array.
std::pair<std::size_t, std::chrono::seconds>
ROOT::RDF::ProgressHelper::RecordEvtCountAndTime() {
  using namespace std::chrono;

  const auto currentEventCount = fProcessedEvents.load();
  const auto eventsPerTimeInterval = currentEventCount - fLastProcessedEvents;
  fLastProcessedEvents = currentEventCount;

  const auto oldPrintTime = fLastPrintTime;
  fLastPrintTime = system_clock::now();

  const duration<double> secondsCurrentInterval = fLastPrintTime - oldPrintTime;
  fEventsPerSecondStatistics[fEventsPerSecondStatisticsIndex++ % fEventsPerSecondStatistics.size()] = eventsPerTimeInterval / secondsCurrentInterval.count();

  return {currentEventCount, duration_cast<seconds>(fLastPrintTime - fBeginTime)};
}

namespace {
  struct PacksOfThree {
    std::array<unsigned int, 7> packsOfThree; // More not supported by std::size_t
    int packCounter = 0;

    PacksOfThree(std::size_t count) {
      for (; count > 0; ++packCounter) {
        assert(packCounter < static_cast<int>(packsOfThree.size()));
        packsOfThree[packCounter] = count % 1000;
        count /= 1000;
      }
      --packCounter;
    }
  };

  /// Format event counts as `6.346.362k`.
  std::ostream& operator<<(std::ostream& stream, const PacksOfThree& packs) {
    for (int i = static_cast<int>(packs.packCounter); i >= 1; --i) {
      if (i == packs.packCounter) stream << packs.packsOfThree[i];
      else stream << std::setw(3) << std::setfill('0') << std::right << packs.packsOfThree[i];

      stream << (i > 1 ? '.' : 'k');
    }
    return stream << std::setfill(' ');
  }

  /// Format std::chrono::seconds as `1:30m`.
  std::ostream& operator<<(std::ostream& stream, std::chrono::seconds elapsedSeconds) {
    stream << std::chrono::duration_cast<std::chrono::minutes>(elapsedSeconds).count()
        << ":" << std::setw(2) << std::right << std::setfill('0') << elapsedSeconds.count() % 60;
    return stream << 'm' << std::setfill(' ');
  }

  struct RestoreStreamState {
    RestoreStreamState(std::ostream& stream) :
      fStream(stream),
      fFlags(stream.flags()),
      fFillChar(stream.fill()) { }
    ~RestoreStreamState() {
      fStream.setf(fFlags);
      fStream.fill(fFillChar);
    }

    std::ostream& fStream;
    std::ios_base::fmtflags fFlags;
    std::ostream::char_type fFillChar;
  };
}

/// Print event and time statistics.
void ROOT::RDF::ProgressHelper::PrintStats(std::ostream& stream, std::size_t currentEventCount, std::chrono::seconds elapsedSeconds) const {
  const auto evtpersec = EvtPerSec();
  RestoreStreamState restore(stream);

  stream << "[" << elapsedSeconds << "  ";

  // Event counts:
  if (fUseShellColours) stream << "\e[32m";

  stream << PacksOfThree(currentEventCount);
  if (fMaxEvents != 0) {
    stream << "/" << PacksOfThree(fMaxEvents);
  }
  stream << " evt  ";

  if (fUseShellColours) stream << "\e[0m";


  // events/s
  stream << std::scientific << std::setprecision(2) << evtpersec << " evt/s";


  // Time statistics:
  if (fMaxEvents != 0) {
    if (fUseShellColours) stream << "\e[35m" << "  ";
    const std::chrono::seconds remainingSeconds( static_cast<long long>((fMaxEvents - currentEventCount) / evtpersec) );
    stream << remainingSeconds << " remaining";
    if (fUseShellColours) stream << "\e[0m";
  }

  stream << "]   ";
}

/// Print a progress bar of width `ProgressHelper::fBarWidth` if `fMaxEvents` is known.
void ROOT::RDF::ProgressHelper::PrintProgressbar(std::ostream& stream, std::size_t currentEventCount) const {
  if (fMaxEvents == 0)
    return;
  RestoreStreamState restore(stream);

  const double completion = double(currentEventCount) / fMaxEvents;
  const unsigned int nBar = completion * fBarWidth;

  std::string bars(std::max(nBar, 1u), '=');
  bars.back() = (nBar == fBarWidth) ? '=' : '>';

  if (fUseShellColours) stream << "\e[33m";
  stream << '|' << std::setfill(' ') << std::setw(fBarWidth) << std::left << bars << "|   ";
  if (fUseShellColours) stream << "\e[0m";
}

std::size_t ROOT::RDF::CountEvents(const char* treename, const char* fileUrl) {
  std::unique_ptr<TFile> file( TFile::Open(fileUrl, "READ") );
  if (!file)
    return 0u;

  std::unique_ptr<TTree> tree( file->Get<TTree>(treename) );
  if (tree) return tree->GetEntries();

  return 0u;
}
