/// \file
/// \ingroup tutorial_dataframe
/// \notebook -nodraw
/// Basic RDataFrame usage.
///
/// This tutorial illustrates the basic features of the RDataFrame class,
/// a utility which allows to interact with data stored in TTrees following
/// a functional-chain like approach.
///
/// \macro_code
/// \macro_output
///
/// \date December 2016
/// \author Enrico Guiraud (CERN)

// ## Preparation

// A simple helper function to fill a test tree: this makes the example
// stand-alone.
void fill_tree(const char *treeName, const char *fileName)
{
   ROOT::RDataFrame d(10);
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

int df001_introduction()
{

   // We prepare an input tree to run on
   auto fileName = "df001_introduction.root";
   auto treeName = "myTree";
   fill_tree(treeName, fileName);

   // We read the tree from the file and create a RDataFrame, a class that
   // allows us to interact with the data contained in the tree.
   // We select a default column, a *branch* to adopt ROOT jargon, which will
   // be looked at if none is specified by the user when dealing with filters
   // and actions.
   ROOT::RDataFrame d(treeName, fileName, {"b1"});

   // ## Operations on the dataframe
   // We now review some *actions* which can be performed on the data frame.
   // All actions but ForEach return a TActionResultPtr<T>. The series of
   // operations on the data frame is not executed until one of those pointers
   // is accessed. If the Foreach action is invoked, the execution is immediate.
   // But first of all, let us we define now our cut-flow with two lambda
   // functions. We can use free functions too.
   auto cutb1 = [](double b1) { return b1 < 5.; };
   auto cutb1b2 = [](int b2, double b1) { return b2 % 2 && b1 < 4.; };

   // ### `Count` action
   // The `Count` allows to retrieve the number of the entries that passed the
   // filters. Here we show how the automatic selection of the column kicks
   // in in case the user specifies none.
   auto entries1 = d.Filter(cutb1) // <- no column name specified here!
                      .Filter(cutb1b2, {"b2", "b1"})
                      .Count();

   std::cout << *entries1 << " entries passed all filters" << std::endl;

   // Filters can be expressed as strings. The content must be C++ code. The
   // name of the variables must be the name of the branches. The code is
   // just in time compiled.
   auto entries2 = d.Filter("b1 < 5.").Count();
   std::cout << *entries2 << " entries passed the string filter" << std::endl;

   // ### `Min`, `Max` and `Mean` actions
   // These actions allow to retrieve statistical information about the entries
   // passing the cuts, if any.
   auto b1b2_cut = d.Filter(cutb1b2, {"b2", "b1"});
   auto minVal = b1b2_cut.Min();
   auto maxVal = b1b2_cut.Max();
   auto meanVal = b1b2_cut.Mean();
   auto nonDefmeanVal = b1b2_cut.Mean("b2"); // <- Column is not the default
   std::cout << "The mean is always included between the min and the max: " << *minVal << " <= " << *meanVal
             << " <= " << *maxVal << std::endl;

   // ### `Take` action
   // The `Take` action allows to retrieve all values of the variable stored in a
   // particular column that passed filters we specified. The values are stored
   // in a list by default, but other collections can be chosen.
   auto b1_cut = d.Filter(cutb1);
   auto b1Vec = b1_cut.Take<double>();
   auto b1List = b1_cut.Take<double, std::list<double>>();

   std::cout << "Selected b1 entries" << std::endl;
   for (auto b1_entry : *b1List)
      std::cout << b1_entry << " ";
   std::cout << std::endl;
   auto b1VecCl = ROOT::GetClass(b1Vec.GetPtr());
   std::cout << "The type of b1Vec is " << b1VecCl->GetName() << std::endl;

   // ### `Histo1D` action
   // The `Histo1D` action allows to fill an histogram. It returns a TH1D filled
   // with values of the column that passed the filters. For the most common
   // types, the type of the values stored in the column is automatically
   // guessed.
   auto hist = d.Filter(cutb1).Histo1D();
   std::cout << "Filled h " << hist->GetEntries() << " times, mean: " << hist->GetMean() << std::endl;

   // ### `Foreach` action
   // The most generic action of all: an operation is applied to all entries.
   // In this case we fill a histogram. In some sense this is a violation of a
   // purely functional paradigm - C++ allows to do that.
   TH1F h("h", "h", 12, -1, 11);
   d.Filter([](int b2) { return b2 % 2 == 0; }, {"b2"}).Foreach([&h](double b1) { h.Fill(b1); });

   std::cout << "Filled h with " << h.GetEntries() << " entries" << std::endl;

   // ## Express your chain of operations with clarity!
   // We are discussing an example here but it is not hard to imagine much more
   // complex pipelines of actions acting on data. Those might require code
   // which is well organised, for example allowing to conditionally add filters
   // or again to clearly separate filters and actions without the need of
   // writing the entire pipeline on one line. This can be easily achieved.
   // We'll show this re-working the `Count` example:
   auto cutb1_result = d.Filter(cutb1);
   auto cutb1b2_result = d.Filter(cutb1b2, {"b2", "b1"});
   auto cutb1_cutb1b2_result = cutb1_result.Filter(cutb1b2, {"b2", "b1"});
   // Now we want to count:
   auto evts_cutb1_result = cutb1_result.Count();
   auto evts_cutb1b2_result = cutb1b2_result.Count();
   auto evts_cutb1_cutb1b2_result = cutb1_cutb1b2_result.Count();

   std::cout << "Events passing cutb1: " << *evts_cutb1_result << std::endl
             << "Events passing cutb1b2: " << *evts_cutb1b2_result << std::endl
             << "Events passing both: " << *evts_cutb1_cutb1b2_result << std::endl;

   // ## Calculating quantities starting from existing columns
   // Often, operations need to be carried out on quantities calculated starting
   // from the ones present in the columns. We'll create in this example a third
   // column the values of which are the sum of the *b1* and *b2* ones, entry by
   // entry. The way in which the new quantity is defined is via a callable.
   // It is important to note two aspects at this point:
   // - The value is created on the fly only if the entry passed the existing
   // filters.
   // - The newly created column behaves as the one present on the file on disk.
   // - The operation creates a new value, without modifying anything. De facto,
   // this is like having a general container at disposal able to accommodate
   // any value of any type.
   // Let's dive in an example:
   auto entries_sum = d.Define("sum", [](double b1, int b2) { return b2 + b1; }, {"b1", "b2"})
                         .Filter([](double sum) { return sum > 4.2; }, {"sum"})
                         .Count();
   std::cout << *entries_sum << std::endl;

   // Additional columns can be expressed as strings. The content must be C++
   // code. The name of the variables must be the name of the branches. The code
   // is just in time compiled.
   auto entries_sum2 = d.Define("sum2", "b1 + b2").Filter("sum2 > 4.2").Count();
   std::cout << *entries_sum2 << std::endl;

   // It is possible at any moment to read the entry number and the processing
   // slot number. The latter may change when implicit multithreading is active.
   // The special columns which provide the entry number and the slot index are
   // called "rdfentry_" and "rdfslot_" respectively. Their types are an unsigned
   // 64 bit integer and an unsigned integer.
   auto printEntrySlot = [](ULong64_t iEntry, unsigned int slot) {
      std::cout << "Entry: " << iEntry << " Slot: " << slot << std::endl;
   };
   d.Foreach(printEntrySlot, {"rdfentry_", "rdfslot_"});

   return 0;
}
