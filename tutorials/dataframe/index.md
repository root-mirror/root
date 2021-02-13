\defgroup tutorial_dataframe Data Frame tutorials
\ingroup Tutorials
\brief These examples show various features of [RDataFrame](classROOT_1_1RDataFrame.html): ROOT's declarative analysis interface.

[RDataFrame](classROOT_1_1RDataFrame.html) offers a high level interface for the analysis of data stored in [TTree](classTTree.html)s, [CSV files](classROOT_1_1RDF_1_1RCsvDS.html) and [other data formats](classROOT_1_1RDF_1_1RDataSource.html).

In addition, multi-threading and other low-level optimisations allow users to exploit all the resources available on their machines transparently.

In a nutshell:
~~~{.cpp}
ROOT::EnableImplicitMT(); // Enable ROOT's implicit multi-threading
ROOT::RDataFrame d("myTree", "file_*.root"); // Interface to TTree and TChain
auto histoA = d.Histo1D("Branch_A");  // Book the filling of a histogram
auto histoB = d.Histo1D("Branch_B");  // Book the filling of another histogram
// Data processing is triggered by the next line, which accesses a booked result for the first time
// All booked results are evaluated during the same parallel event loop.
histoA->Draw(); // <-- event loop runs here!
histoB->Draw(); // HistoB has already been filled, no event loop is run here
~~~

Explore the examples below or go to [RDataFrame's user guide](classROOT_1_1RDataFrame.html).
