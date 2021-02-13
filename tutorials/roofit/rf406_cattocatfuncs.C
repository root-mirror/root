/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Data and categories: demonstration of discrete-->discrete (invertible) functions
///
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPolynomial.h"
#include "RooCategory.h"
#include "RooMappedCategory.h"
#include "RooMultiCategory.h"
#include "RooSuperCategory.h"
#include "Roo1DTable.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf406_cattocatfuncs()
{
   // C o n s t r u c t  t w o   c a t e g o r i e s
   // ----------------------------------------------

   // Define a category with labels only
   RooCategory tagCat("tagCat", "Tagging category");
   tagCat.defineType("Lepton");
   tagCat.defineType("Kaon");
   tagCat.defineType("NetTagger-1");
   tagCat.defineType("NetTagger-2");
   tagCat.Print();

   // Define a category with explicitly numbered states
   RooCategory b0flav("b0flav", "B0 flavour eigenstate");
   b0flav.defineType("B0", -1);
   b0flav.defineType("B0bar", 1);
   b0flav.Print();

   // Construct a dummy dataset with random values of tagCat and b0flav
   RooRealVar x("x", "x", 0, 10);
   RooPolynomial p("p", "p", x);
   RooDataSet *data = p.generate(RooArgSet(x, b0flav, tagCat), 10000);

   // C r e a t e   a   c a t - > c a t   m  a p p i n g   c a t e g o r y
   // ---------------------------------------------------------------------

   // A RooMappedCategory is category->category mapping function based on string expression
   // The constructor takes an input category an a default state name to which unassigned
   // states are mapped
   RooMappedCategory tcatType("tcatType", "tagCat type", tagCat, "Cut based");

   // Enter fully specified state mappings
   tcatType.map("Lepton", "Cut based");
   tcatType.map("Kaon", "Cut based");

   // Enter a wilcard expression mapping
   tcatType.map("NetTagger*", "Neural Network");

   // Make a table of the mapped category state multiplicity in data
   Roo1DTable *mtable = data->table(tcatType);
   mtable->Print("v");

   // C r e a t e   a   c a t   X   c a t   p r o d u c t   c a t e g o r y
   // ----------------------------------------------------------------------

   // A SUPER-category is 'product' of _lvalue_ categories. The state names of a super
   // category is a composite of the state labels of the input categories
   RooSuperCategory b0Xtcat("b0Xtcat", "b0flav X tagCat", RooArgSet(b0flav, tagCat));

   // Make a table of the product category state multiplicity in data
   Roo1DTable *stable = data->table(b0Xtcat);
   stable->Print("v");

   // Since the super category is an lvalue, assignment is explicitly possible
   b0Xtcat.setLabel("{B0bar;Lepton}");

   // A MULTI-category is a 'product' of any category (function). The state names of a super
   // category is a composite of the state labels of the input categories
   RooMultiCategory b0Xttype("b0Xttype", "b0flav X tagType", RooArgSet(b0flav, tcatType));

   // Make a table of the product category state multiplicity in data
   Roo1DTable *xtable = data->table(b0Xttype);
   xtable->Print("v");
}
