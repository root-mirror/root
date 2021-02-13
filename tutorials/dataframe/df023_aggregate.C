/// \file
/// \ingroup tutorial_dataframe
/// \notebook
/// Use the Aggregate action to specify arbitrary data aggregations.
///
/// This tutorial shows how to use the Aggregate action to evaluate the product of all the elements of a column.
/// This operation may be performed using a Reduce action, however aggregate is used for the sake of the tutorial
///
/// \macro_code
/// \macro_output
///
/// \date July 2018
/// \authors Enrico Guiraud, Danilo Piparo (CERN), Massimo Tumolo (Politecnico di Torino)

void df023_aggregate()
{

   // Column to be aggregated
   const std::string columnName = "x";

   ROOT::EnableImplicitMT(2);
   auto rdf = ROOT::RDataFrame(5);
   auto d = rdf.Define(columnName, "rdfentry_ + 1.");

   // Aggregator function. It receives an accumulator (acc) and a column value (x). The variable acc is shared among the
   // calls, so the function has to specify how the value has to be aggregated in the accumulator.
   auto aggregator = [](double acc, double x) { return acc * x; };

   // If multithread is enabled, the aggregator function will be called by more threads and will produce a vector of
   // partial accumulators. The merger function performs the final aggregation of these partial results.
   auto merger = [](std::vector<double> &accumulators) {
      auto size = accumulators.size();
      for (int i = 1; i < size; ++i) {
         accumulators[0] *= accumulators[i];
      }
   };

   // The accumulator is initialized at this value by every thread.
   double initValue = 1.;

   // Multiplies all elements of the column "x"
   auto result = d.Aggregate(aggregator, merger, columnName, initValue);

   std::cout << *result << std::endl;
}
