/// \file
/// \ingroup tutorial_dataframe
/// \notebook
/// Use the Display action to inspect entry values.
///
/// \macro_code
/// \macro_output
///
/// \date August 2018
/// \authors Enrico Guiraud, Danilo Piparo, Enric Tejedor Saavedra (CERN), Massimo Tumolo (Politecnico di Torino)

void df024_Display()
{
   // Setting up a Dataframe with some data
   unsigned long long y = 1;
   int x = 1;
   double w = 1;
   double z = 1;
   ROOT::RDataFrame df(10);
   auto d = df.Define("y", [&y]() { return y *= 100; }) // A column with ulongs
              .Define("x",
                      [&x]() {
                         return std::vector<int>({x++, x++, x++, x++});
                      })                                // A column with four-elements collection
              .Define("w", [&w]() { return w *= 1.8; }) // A column with doubles
              .Define("z", [&z]() {
                 z *= 1.1;
                 return std::vector<std::vector<double>>({{z, ++z}, {z, ++z}, {z, ++z}});
              }); // A column of matrices

   // Preparing the RResultPtr<RDisplay> object with all columns and default number of entries
   auto d1 = d.Display("");
   // Preparing the RResultPtr<RDisplay> object with two columns and default number of entries
   auto d2 = d.Display({"x", "y"});

   // Printing the short representations, the event loop will run
   std::cout << "The following is the representation of all columns with the default nr of entries" << std::endl;
   d1->Print();
   std::cout << "\n\nThe following is the representation of two columns with the default nr of entries" << std::endl;
   d2->Print();
}
