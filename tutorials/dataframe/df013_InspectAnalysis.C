/// \file
/// \ingroup tutorial_dataframe
/// \notebook -draw
/// Use callbacks to update a plot and a progress bar during the event loop.
///
/// Showcase registration of callback functions that act on partial results while
/// the event-loop is running using `OnPartialResult` and `OnPartialResultSlot`.
/// This tutorial is not meant to run in batch mode.
///
/// \macro_code
///
/// \date September 2017
/// \author Enrico Guiraud (CERN)

using namespace ROOT; // RDataFrame lives in here

void df013_InspectAnalysis()
{
   ROOT::EnableImplicitMT();
   const auto poolSize = ROOT::GetThreadPoolSize();
   const auto nSlots = 0 == poolSize ? 1 : poolSize;

   // ## Setup a simple RDataFrame
   // We start by creating a RDataFrame with a good number of empty events
   const auto nEvents = nSlots * 10000ull;
   RDataFrame d(nEvents);

   // `heavyWork` is a lambda that fakes some interesting computation and just returns a normally distributed double
   TRandom r;
   auto heavyWork = [&r]() {
      for (volatile int i = 0; i < 1000000; ++i)
         ;
      return r.Gaus();
   };

   // Let's define a column "x" produced by invoking `heavyWork` for each event
   // `df` stores a modified data-frame that contains "x"
   auto df = d.Define("x", heavyWork);

   // Now we register a histogram-filling action with the RDataFrame.
   // `h` can be used just like a pointer to TH1D but it is actually a TResultProxy<TH1D>, a smart object that triggers
   // an event-loop to fill the pointee histogram if needed.
   auto h = df.Histo1D<double>({"browserHisto", "", 100, -2., 2.}, "x");

   // ## Use the callback mechanism to draw the histogram on a TBrowser while it is being filled
   // So far we have registered a column "x" to a data-frame with `nEvents` events and we registered the filling of a
   // histogram with the values of column "x".
   // In the following we will register three functions for execution during the event-loop:
   // - one is to be executed once just before the loop and adds a partially-filled histogram to a TBrowser
   // - the next is executed every 50 events and draws the partial histogram on the TBrowser's TPad
   // - another callback is responsible of updating a simple progress bar from multiple threads

   // First off we create a TBrowser that contains a "RDFResults" directory
   auto dfDirectory = new TMemFile("RDFResults", "RECREATE");
   auto browser = new TBrowser("b", dfDirectory);
   // The global pad should now be set to the TBrowser's canvas, let's store its value in a local variable
   auto browserPad = gPad;

   // A useful feature of `TResultProxy` is its `OnPartialResult` method: it allows us to register a callback that is
   // executed once per specified number of events during the event-loop, on "partial" versions of the result objects
   // contained in the `TResultProxy`. In this case, the partial result is going to be a histogram filled with an
   // increasing number of events.
   // Instead of requesting the callback to be executed every N entries, this time we use the special value `kOnce` to
   // request that it is executed once right before starting the event-loop.
   // The callback is a C++11 lambda that registers the partial result object in `dfDirectory`.
   h.OnPartialResult(h.kOnce, [dfDirectory](TH1D &h_) { dfDirectory->Add(&h_); });
   // Note that we called `OnPartialResult` with a dot, `.`, since this is a method of `TResultProxy` itself.
   // We do not want to call `OnPartialResult` on the pointee histogram!)

   // Multiple callbacks can be registered on the same `TResultProxy` (they are executed one after the other in the
   // same order as they were registered). We now request that the partial result is drawn and the TBrowser's TPad is
   // updated every 50 events.
   h.OnPartialResult(50, [&browserPad](TH1D &hist) {
      if (!browserPad)
         return; // in case root -b was invoked
      browserPad->cd();
      hist.Draw();
      browserPad->Update();
      // This call tells ROOT to process all pending GUI events
      // It allows users to use the TBrowser as usual while the event-loop is running
      gSystem->ProcessEvents();
   });

   // Finally, we would like to print a progress bar on the terminal to show how the event-loop is progressing.
   // To take into account _all_ events we use `OnPartialResultSlot`: when Implicit Multi-Threading is enabled, in fact,
   // `OnPartialResult` invokes the callback only in one of the worker threads, and always returns that worker threads'
   // partial result. This is useful because it means we don't have to worry about concurrent execution and
   // thread-safety of the callbacks if we are happy with just one threads' partial result.
   // `OnPartialResultSlot`, on the other hand, invokes the callback in each one of the worker threads, every time a
   // thread finishes processing a batch of `everyN` events. This is what we want for the progress bar, but we need to
   // take care that two threads will not print to terminal at the same time: we need a std::mutex for synchronization.
   std::string progressBar;
   std::mutex barMutex; // Only one thread at a time can lock a mutex. Let's use this to avoid concurrent printing.
   // Magic numbers that yield good progress bars for nSlots = 1,2,4,8
   const auto everyN = nSlots == 8 ? 1000 : 100ull * nSlots;
   const auto barWidth = nEvents / everyN;
   h.OnPartialResultSlot(everyN, [&barWidth, &progressBar, &barMutex](unsigned int /*slot*/, TH1D & /*partialHist*/) {
      std::lock_guard<std::mutex> l(barMutex); // lock_guard locks the mutex at construction, releases it at destruction
      progressBar.push_back('#');
      // re-print the line with the progress bar
      std::cout << "\r[" << std::left << std::setw(barWidth) << progressBar << ']' << std::flush;
   });

   // ## Running the analysis
   // So far we told RDataFrame what we want to happen during the event-loop, but we have not actually run any of those
   // actions: the TBrowser is still empty, the progress bar has not been printed even once, and we haven't produced
   // a single data-point!
   // As usual with RDataFrame, the event-loop is triggered by accessing the contents of a TResultProxy for the first
   // time. Let's run!
   std::cout << "Analysis running..." << std::endl;
   h->Draw(); // the final, complete result will be drawn after the event-loop has completed.
   std::cout << "\nDone!" << std::endl;

   // Finally, some book-keeping: in the TMemFile that we are using as TBrowser directory, we substitute the partial
   // result with a clone of the final result (the "original" final result will be deleted at the end of the macro).
   dfDirectory->Clear();
   auto clone = static_cast<TH1D *>(h->Clone());
   clone->SetDirectory(nullptr);
   dfDirectory->Add(clone);
   if (!browserPad)
      return; // in case root -b was invoked
   browserPad->cd();
   clone->Draw();
}
