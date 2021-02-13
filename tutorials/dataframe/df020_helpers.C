/// \file
/// \ingroup tutorial_dataframe
/// \notebook
/// Show usage of RDataFrame's helper tools, contained in ROOT/RDFHelpers.hxx.
///
/// \macro_code
///
/// \date July 2018
/// \author Enrico Guiraud (CERN)

void df020_helpers()
{
   // First of all, we create a dataframe with 3 entries and define two simple columns
   const auto nEntries = 3;
   ROOT::RDataFrame _df(nEntries);
   auto df = _df.Define("one", [] { return 1; }).Define("two", [] { return 2; });

   // *** Not ***
   // This helper takes a callable `f` (which must return a `bool`) and produces a new callable which takes the same
   // arguments as `f` but returns its negated result. `Not` is useful to invert the check performed by a given Filter.
   // Here we define a simple lambda that checks whether a value is equal to 1, and invert it with Not:
   auto isOne = [] (int a) { return a == 1; };
   auto isNotOne = ROOT::RDF::Not(isOne);

   // Both `isOne` and `isNotOne` are callables that we can use in `Filters`:
   auto c1 = df.Filter(isOne, {"one"}).Count();
   auto c2 = df.Filter(isNotOne, {"two"}).Count();
   // Both counts are equal to the total number of entries, as both Filters always pass.
   R__ASSERT(*c1 == nEntries);
   R__ASSERT(*c2 == nEntries);

   // *** PassAsVec ***
   // Consider the following function, which checks if a vector consists of two elements equal to 1 and 2:
   auto checkOneTwo = [] (const std::vector<int> &v) { return v.size() == 2 && v[0] == 1 && v[1] == 2; };
   // The following line, although it looks reasonable, would _not_ run correctly:
   // df.Filter(checkOneTwo, {"one", "two"});
   // The reason is that `Filter(..., {"one", "two"})` expects a callable that takes exactly two integers, while
   // `checkOneTwo` actually takes a vector of integers (i.e. it does not have the right signature).
   // PassAsVec helps passing down the single values "one", "two" to `checkOneTwo` as a collection: it takes a callable
   // `f` that expects a collection as argument and returns a new callable that takes single arguments instead, passes
   // them down to `f` and returns what `f` returns.
   // PassAsVec requires that number of arguments and their type is specified as template argument.
   // Here's an example usage (remember, PassAsVec(f) returns a new callable!):
   auto c3 = df.Filter(ROOT::RDF::PassAsVec<2, int>(checkOneTwo), {"one", "two"}).Count();
   R__ASSERT(*c3 == nEntries);
}
