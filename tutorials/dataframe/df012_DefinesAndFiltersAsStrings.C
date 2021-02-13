/// \file
/// \ingroup tutorial_dataframe
/// \notebook -nodraw
/// Use just-in-time-compiled Filters and Defines for quick prototyping.
///
/// This tutorial illustrates how to save some typing when using RDataFrame
/// by invoking functions that perform jit-compiling at runtime.
///
/// \macro_code
/// \macro_output
///
/// \date October 2017
/// \author Guilherme Amadio (CERN)

void df012_DefinesAndFiltersAsStrings()
{
   // We will inefficiently calculate an approximation of pi by generating
   // some data and doing very simple filtering and analysis on it

   // We start by creating an empty dataframe where we will insert 10 million
   // random points in a square of side 2.0 (that is, with an inscribed circle
   // of radius 1.0)

   size_t npoints = 10000000;
   ROOT::RDataFrame df(npoints);

   // Define what we want inside the dataframe. We do not need to define p as an array,
   // but we do it here to demonstrate how to use jitting with RDataFrame

   // NOTE: Although it's possible to use "for (auto&& x : p)" below, it will
   // shadow the name of the data column "x", and may cause compilation failures
   // if the local variable and the data column are of different types or the
   // local x variable is declared in the global scope of the lambda function

   auto pidf = df.Define("x", "gRandom->Uniform(-1.0, 1.0)")
                 .Define("y", "gRandom->Uniform(-1.0, 1.0)")
                 .Define("p", "std::array<double, 2> v{x, y}; return v;")
                 .Define("r", "double r2 = 0.0; for (auto&& x : p) r2 += x*x; return sqrt(r2);");

   // Now we have a dataframe with columns x, y, p (which is a point based on x
   // and y), and the radius r = sqrt(x*x + y*y). In order to approximate pi, we
   // need to know how many of our data points fall inside the unit circle compared
   // with the total number of points. The ratio of the areas is
   //
   //     A_circle / A_square = pi r*r / l * l, where r = 1.0, and l = 2.0
   //
   // Therefore, we can approximate pi with 4 times the number of points inside the
   // unit circle over the total number of points in our dataframe:

   auto incircle = *(pidf.Filter("r <= 1.0").Count());

   double pi_approx = 4.0 * incircle / npoints;

   std::cout << "pi is approximately equal to " << pi_approx << std::endl;
}
