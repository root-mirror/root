/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Organization and simultaneous fits: RooArgSet and RooArgList tools and tricks
///
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooCategory.h"
using namespace RooFit;

void rf508_listsetmanip()
{

   // C r e a t e   d u m m y   o b j e c t s
   // ---------------------------------------

   // Create some variables
   RooRealVar a("a", "a", 1, -10, 10);
   RooRealVar b("b", "b", 2, -10, 10);
   RooRealVar c("c", "c", 3, -10, 10);
   RooRealVar d("d", "d", 4, -10, 10);
   RooRealVar x("x", "x", 0, -10, 10);
   c.setError(0.5);
   a.setConstant();
   b.setConstant();

   // Create a category
   RooCategory e("e", "e");
   e.defineType("sig");
   e.defineType("bkg");

   // Create a pdf
   RooGaussian g("g", "g", x, a, b);

   // C r e a t i n g ,   f i l l i n g   R o o A r g S e t s
   // -------------------------------------------------------

   // A RooArgSet is a set of RooAbsArg objects. Each object in the set must have
   // a unique name

   // Set constructors exists with up to 9 initial arguments
   RooArgSet s(a, b);

   // At any time objects can be added with add()
   s.add(e);

   // Add up to 9 additional arguments in one call
   s.add(RooArgSet(c, d));

   // Sets can contain any type of RooAbsArg, also pdf and functions
   s.add(g);

   // Remove element d
   s.remove(d);

   // A c c e s s i n g   R o o A r g S e t   c o n t e n t s
   // -------------------------------------------------------

   // You can look up objects by name
   RooAbsArg *aptr = s.find("a");

   // Construct a subset by name
   RooArgSet *subset1 = (RooArgSet *)s.selectByName("a,b,c");

   // Construct asubset by attribute
   RooArgSet *subset2 = (RooArgSet *)s.selectByAttrib("Constant", kTRUE);

   // Construct the subset of overlapping contents with another set
   RooArgSet s1(a, b, c);
   RooArgSet s2(c, d, e);
   RooArgSet *subset3 = (RooArgSet *)s1.selectCommon(s2);

   // O w n i n g   R o o A r g S e t s
   // ---------------------------------

   // Create a RooArgSet that owns its components
   // A set either owns all of its components or none,
   // so once addOwned() is used, add() can no longer be
   // used and will result in an error message

   RooRealVar *ac = (RooRealVar *)a.clone("a");
   RooRealVar *bc = (RooRealVar *)b.clone("b");
   RooRealVar *cc = (RooRealVar *)c.clone("c");

   RooArgSet s3;
   s3.addOwned(RooArgSet(*ac, *bc, *cc));

   // Another possibility is to add an owned clone
   // of an object instead of the original
   s3.addClone(RooArgSet(d, e, g));

   // A clone of a owning set is non-owning and its
   // contents is owned by the originating owning set
   RooArgSet *sclone = (RooArgSet *)s3.Clone("sclone");

   // To make a clone of a set and its contents use
   // the snapshot method
   RooArgSet *sclone2 = (RooArgSet *)s3.snapshot();

   // If a set contains function objects, only the head node
   // is cloned in a snapshot. To make a snapshot of all
   // servers of a function object do as follows. The result
   // of a RooArgSet snapshot with deepCloning option is a set
   // of cloned objects, and all their clone (recursive) server
   // dependencies, that together form a self-consistent
   // set that is free of external dependencies

   RooArgSet *sclone3 = (RooArgSet *)s3.snapshot(kTRUE);

   // S e t   p r i n t i n g
   // ------------------------

   // Inline printing only show list of names of contained objects
   cout << "sclone = " << (*sclone) << endl;

   // Plain print shows the same, prefixed by name of the set
   sclone->Print();

   // Standard printing shows one line for each item with the items name, class name and value
   sclone->Print("s");

   // Verbose printing adds each items arguments, address and 'extras' as defined by the object
   sclone->Print("v");

   // U s i n g   R o o A r g L i s t s
   // ---------------------------------

   // List constructors exists with up to 9 initial arguments
   RooArgList l(a, b, c, d);

   // Lists have an explicit order and allow multiple arguments with the same name
   l.add(RooArgList(a, b, c, d));

   // Access by index is provided
   RooAbsArg *arg4 = l.at(4);
}
