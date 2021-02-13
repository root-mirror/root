/// \file
/// \ingroup tutorial_dataframe
/// \notebook
/// Display cut/Filter efficiencies with RDataFrame.
///
/// This tutorial shows how to get information about the efficiency of the filters
/// applied
///
/// \macro_code
/// \macro_output
///
/// \date December 2016
/// \author Danilo Piparo (CERN)

using FourVector = ROOT::Math::XYZTVector;
using FourVectors = std::vector<FourVector>;
using CylFourVector = ROOT::Math::RhoEtaPhiVector;

// A simple helper function to fill a test tree: this makes the example
// stand-alone.
void fill_tree(const char *treeName, const char *fileName)
{
   ROOT::RDataFrame d(50);
   int i(0);
   d.Define("b1", [&i]() { return (double)i; })
      .Define("b2",
              [&i]() {
                 auto j = i * i;
                 ++i;
                 return j;
              })
      .Snapshot(treeName, fileName);
}

void df004_cutFlowReport()
{

   // We prepare an input tree to run on
   auto fileName = "df004_cutFlowReport.root";
   auto treeName = "myTree";
   fill_tree(treeName, fileName);

   // We read the tree from the file and create a RDataFrame
   ROOT::RDataFrame d(treeName, fileName, {"b1", "b2"});

   // ## Define cuts and create the report
   // Here we define two simple cuts
   auto cut1 = [](double b1) { return b1 > 25.; };
   auto cut2 = [](int b2) { return 0 == b2 % 2; };

   // An optional string parameter name can be passed to the Filter method to create a named filter.
   // Named filters work as usual, but also keep track of how many entries they accept and reject.
   auto filtered1 = d.Filter(cut1, {"b1"}, "Cut1");
   auto filtered2 = d.Filter(cut2, {"b2"}, "Cut2");

   auto augmented1 = filtered2.Define("b3", [](double b1, int b2) { return b1 / b2; });
   auto cut3 = [](double x) { return x < .5; };
   auto filtered3 = augmented1.Filter(cut3, {"b3"}, "Cut3");

   // Statistics are retrieved through a call to the Report method:
   // when Report is called on the main RDataFrame object, it retrieves stats
   // for all named filters declared up to that point.
   // When called on a stored chain state (i.e. a chain/graph node), it
   // retrieves stats for all named filters in the section of the chain between
   // the main RDataFrame and that node (included).
   // Stats are printed in the same order as named filters have been added to
   // the graph, and refer to the latest event-loop that has been run using the
   // relevant RDataFrame.
   std::cout << "Cut3 stats:" << std::endl;
   filtered3.Report()->Print();

   // It is not only possible to print the information about cuts, but also to
   // retrieve it to then use it programmatically.
   std::cout << "All stats:" << std::endl;
   auto allCutsReport = d.Report();
   allCutsReport->Print();

   // We can now loop on the cuts
   std::cout << "Name\tAll\tPass\tEfficiency" << std::endl;
   for (auto &&cutInfo : allCutsReport) {
      std::cout << cutInfo.GetName() << "\t" << cutInfo.GetAll() << "\t" << cutInfo.GetPass() << "\t"
                << cutInfo.GetEff() << " %" << std::endl;
   }

   // Or get information about them individually
   auto cutName = "Cut1";
   auto cut = allCutsReport->At("Cut1");
   std::cout << cutName << " efficiency is " << cut.GetEff() << " %" << std::endl;
}
