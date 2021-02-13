/// \file
/// \ingroup tutorial_dataframe
/// \notebook -draw
/// Use Range to limit the amount of data processed.
///
/// This tutorial shows how to express the concept of ranges when working with the RDataFrame.
///
/// \macro_code
/// \macro_output
///
/// \date March 2017
/// \author Danilo Piparo (CERN)

// A simple helper function to fill a test tree: this makes the example
// stand-alone.
void fill_tree(const char *treeName, const char *fileName)
{
   ROOT::RDataFrame d(100);
   int i(0);
   d.Define("b1", [&i]() { return i; })
      .Define("b2",
              [&i]() {
                 float j = i * i;
                 ++i;
                 return j;
              })
      .Snapshot(treeName, fileName);
}

int df006_ranges()
{

   // We prepare an input tree to run on
   auto fileName = "df006_ranges.root";
   auto treeName = "myTree";
   fill_tree(treeName, fileName);

   // We read the tree from the file and create a RDataFrame.
   ROOT::RDataFrame d(treeName, fileName);

   // ## Usage of ranges
   // Now we'll count some entries using ranges
   auto c_all = d.Count();

   // This is how you can express a range of the first 30 entries
   auto d_0_30 = d.Range(0, 30);
   auto c_0_30 = d_0_30.Count();

   // This is how you pick all entries from 15 onwards
   auto d_15_end = d.Range(15, 0);
   auto c_15_end = d_15_end.Count();

   // We can use a stride too, in this case we pick an event every 3
   auto d_15_end_3 = d.Range(15, 0, 3);
   auto c_15_end_3 = d_15_end_3.Count();

   // The Range is a 1st class citizen in the RDataFrame graph:
   // not only actions (like Count) but also filters and new columns can be added to it.
   auto d_0_50 = d.Range(0, 50);
   auto c_0_50_odd_b1 = d_0_50.Filter("1 == b1 % 2").Count();

   // An important thing to notice is that the counts of a filter are relative to the
   // number of entries a filter "sees". Therefore, if a Range depends on a filter,
   // the Range will act on the entries passing the filter only.
   auto c_0_3_after_even_b1 = d.Filter("0 == b1 % 2").Range(0, 3).Count();

   // Ok, time to wrap up: let's print all counts!
   cout << "Usage of ranges:\n"
        << " - All entries: " << *c_all << endl
        << " - Entries from 0 to 30: " << *c_0_30 << endl
        << " - Entries from 15 onwards: " << *c_15_end << endl
        << " - Entries from 15 onwards in steps of 3: " << *c_15_end_3 << endl
        << " - Entries from 0 to 50, odd only: " << *c_0_50_odd_b1 << endl
        << " - First three entries of all even entries: " << *c_0_3_after_even_b1 << endl;

   return 0;
}
