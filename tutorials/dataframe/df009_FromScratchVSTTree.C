/// \file
/// \ingroup tutorial_dataframe
/// \notebook -nodraw
/// Compare creation of a ROOT dataset with RDataFrame and TTree.
///
/// This tutorial illustrates how simpler it can be to use a
/// RDataFrame to create a dataset with respect to the usage
/// of the TTree interfaces.
///
/// \macro_code
///
/// \date August 2017
/// \author Danilo Piparo (CERN)

// ##This is the classic way of creating a ROOT dataset
// The steps are:
// - Create a file
// - Create a tree associated to the file
// - Define the variables to write in the entries
// - Define the branches associated to those variables
// - Write the event loop to set the right value to the variables
//   - Call TTree::Fill to save the value of the variables
// - Write the TTree
// - Close the file
void classicWay()
{
   TFile f("df009_FromScratchVSTTree_classic.root", "RECREATE");
   TTree t("treeName", "treeName");
   double b1;
   int b2;
   t.Branch("b1", &b1);
   t.Branch("b2", &b2);
   for (int i = 0; i < 10; ++i) {
      b1 = i;
      b2 = i * i;
      t.Fill();
   }
   t.Write();
   f.Close();
}

// ##This is the RDF way of creating a ROOT dataset
// Few lines are needed to achieve the same result.
// Parallel creation of the TTree is not supported in the
// classic method.
// In this case the steps are:
// - Create an empty RDataFrame
// - If needed, define variables for the functions used to fill the branches
// - Create new columns expressing their content with lambdas, functors, functions or strings
// - Invoke the Snapshot action
//
// Parallelism is not the only advantage. Starting from an existing dataset and
// filter it, enrich it with new columns, leave aside some other columns and
// write a new dataset becomes very easy to do.
void RDFWay()
{
   ROOT::RDataFrame df(10);
   auto b = 0.;
   df.Define("b1", [&b]() { return b++; })
     .Define("b2", "(int) b1 * b1") // This can even be a string
     .Snapshot("treeName", "df009_FromScratchVSTTree_df.root");
}

void df009_FromScratchVSTTree()
{

   classicWay();
   RDFWay();
}
