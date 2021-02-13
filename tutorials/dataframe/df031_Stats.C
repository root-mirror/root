/// \file
/// \ingroup tutorial_dataframe
/// \notebook
/// Use the Stats action to extract the statistics of a column.
///
/// Extract the statistics relative to RDataFrame columns and store them
/// in TStatistic instances.
///
/// \macro_code
/// \macro_output
///
/// \date April 2019
/// \author Danilo Piparo (CERN)

void df031_Stats() {

    // Create a data frame and add two columns: one for the values and one for the weight.
    ROOT::RDataFrame r(256);
    auto rr = r.Define("v", [](ULong64_t e){return e;}, {"rdfentry_"})
               .Define("w", [](ULong64_t e){return 1./(e+1);}, {"v"});

    // Now extract the statistics, weighted, unweighted - with and without explicit types.
    auto stats_eu = rr.Stats<ULong64_t>("v");
    auto stats_ew = rr.Stats<ULong64_t, double>("v", "w");
    auto stats_iu = rr.Stats("v");
    auto stats_iw = rr.Stats("v", "w");

    // Now print them: they are all identical of course!
    stats_eu->Print();
    stats_ew->Print();
    stats_iu->Print();
    stats_iw->Print();
}
