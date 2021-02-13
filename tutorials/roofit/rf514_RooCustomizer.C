/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Using the RooCustomizer to create multiple PDFs that share a lot of properties, but have unique parameters for each category.
/// As an extra complication, some of the new parameters need to be functions
/// of a mass parameter.
///
/// \macro_output
/// \macro_code
///
/// \author Stephan Hageboeck, CERN


#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooCustomizer.h"
#include "RooCategory.h"
#include "RooFormulaVar.h"
#include <iostream>

void rf514_RooCustomizer() {

  // Define a proto model that will be used as the template for each category
  // ---------------------------------------------------------------------------

  RooRealVar E("Energy","Energy",0,3000);

  RooRealVar meanG("meanG","meanG", 100., 0., 3000.);
  RooRealVar sigmaG("sigmaG","sigmaG", 3.);
  RooGaussian gauss("gauss", "gauss", E, meanG, sigmaG);

  RooRealVar pol1("pol1", "Constant of the polynomial", 1, -10, 10);
  RooPolynomial linear("linear", "linear", E, pol1);

  RooRealVar yieldSig("yieldSig", "yieldSig", 1, 0, 1.E4);
  RooRealVar yieldBkg("yieldBkg", "yieldBkg", 1, 0, 1.E4);

  RooAddPdf model("model", "S + B model",
      RooArgList(gauss,linear),
      RooArgList(yieldSig, yieldBkg));

  std::cout << "The proto model before customisation:" << std::endl;
  model.Print("T"); // "T" prints the model as a tree


  // Build the categories
  RooCategory sample("sample","sample");
  sample["Sample1"] = 1;
  sample["Sample2"] = 2;
  sample["Sample3"] = 3;


  // Start to customise the proto model that was defined above.
  // ---------------------------------------------------------------------------

  // We need two sets for bookkeeping of PDF nodes:
  RooArgSet newLeafs;           // This set collects leafs that are created in the process.
  RooArgSet allCustomiserNodes; // This set lists leafs that have been used in a replacement operation.


  // 1. Each sample should have its own mean for the gaussian
  // The customiser will make copies of `meanG` for each category.
  // These will all appear in the set `newLeafs`, which will own the new nodes.
  RooCustomizer cust(model, sample, newLeafs, &allCustomiserNodes);
  cust.splitArg(meanG, sample);


  // 2. Each sample should have its own signal yield, but there is an extra complication:
  // We need the yields 1 and 2 to be a function of the variable "mass".
  // For this, we pre-define nodes with exacly the names that the customiser would have created automatically,
  // that is, "<nodeName>_<categoryName>", and we register them in the set of customiser nodes.
  // The customiser will pick them up instead of creating new ones.
  // If we don't provide one (e.g. for "yieldSig_Sample3"), it will be created automatically by cloning `yieldSig`.
  RooRealVar mass("M", "M", 1, 0, 12000);
  RooFormulaVar yield1("yieldSig_Sample1", "Signal yield in the first sample", "M/3.360779", mass);
  RooFormulaVar yield2("yieldSig_Sample2", "Signal yield in the second sample", "M/2", mass);
  allCustomiserNodes.add(yield1);
  allCustomiserNodes.add(yield2);

  // Instruct the customiser to replace all yieldSig nodes for each sample:
  cust.splitArg(yieldSig, sample);


  // Now we can start building the PDFs for all categories:
  auto pdf1 = cust.build("Sample1");
  auto pdf2 = cust.build("Sample2");
  auto pdf3 = cust.build("Sample3");

  // And we inspect the two PDFs
  std::cout << "\nPDF 1 with a yield depending on M:" << std::endl;
  pdf1->Print("T");
  std::cout << "\nPDF 2 with a yield depending on M:" << std::endl;
  pdf2->Print("T");
  std::cout << "\nPDF 3 with a free yield:" << std::endl;
  pdf3->Print("T");

  std::cout << "\nThe following leafs have been created automatically while customising:" << std::endl;
  newLeafs.Print("V");


  // If we needed to set reasonable values for the means of the gaussians, this could be done as follows:
  auto& meanG1 = static_cast<RooRealVar&>(allCustomiserNodes["meanG_Sample1"]);
  meanG1.setVal(200);
  auto& meanG2 = static_cast<RooRealVar&>(allCustomiserNodes["meanG_Sample2"]);
  meanG2.setVal(300);

  std::cout << "\nThe following leafs have been used while customising"
    << "\n\t(partial overlap with the set of automatically created leaves."
    << "\n\ta new customiser for a different PDF could reuse them if necessary.):" << std::endl;
  allCustomiserNodes.Print("V");


}
