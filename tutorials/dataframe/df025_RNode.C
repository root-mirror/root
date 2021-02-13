/// \file
/// \ingroup tutorial_dataframe
/// \notebook
/// Manipulate RDF objects in functions, loops and conditional branches.
///
/// Each RDataFrame object has its own type. It helps with performance,
/// but sometimes it gets in the way of writing simple code that manages RDF objects.
/// Luckily, every RDF object can be converted to the generic RNode type.
/// This tutorial shows how to take advantage of RNode to easily manipulate RDataFrames.
///
/// \macro_code
/// \macro_output
///
/// \date June 2020
/// \authors Danilo Piparo, Enrico Guiraud (CERN)

/// A generic function that takes an RDF object and applies a string filter
ROOT::RDF::RNode AddFilter(ROOT::RDF::RNode node, string_view filterStr)
{
   return node.Filter(filterStr);
}

void df025_RNode()
{
   ROOT::RDataFrame df(8);

   // Using the generic AddFilter helper function defined above: RNode in, RNode out
   auto f1 = AddFilter(df, "rdfentry_ > 0");
   auto f2 = f1.Filter([](ULong64_t e) { return e > 1; }, {"rdfentry_"});

   // Conditionally applying a filter is simple with ROOT::RDF::RNode
   bool someCondition = true;
   auto maybe_filtered = ROOT::RDF::RNode(f2);
   if (someCondition)
      maybe_filtered = maybe_filtered.Filter("rdfentry_ > 3");

   // Adding new columns with Define in a loop is simple thanks to ROOT::RDF::RNode
   auto with_columns = ROOT::RDF::RNode(maybe_filtered);
   for (auto i = 0; i < 3; ++i)
      with_columns = with_columns.Define("x" + std::to_string(i), "42");

   // RNodes can be used exactly like any other RDF object
   std::cout << "Entries passing the selection: " << with_columns.Count().GetValue() << std::endl;
}
