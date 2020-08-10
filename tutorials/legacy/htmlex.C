/// \file
/// \ingroup tutorial_legacy
/// This file demonstrates how THtml can document sources.
///
/// See the [Users Guide](https://root.cern/root/htmldoc/guides/users-guide/ROOTUsersGuide.html)
/// chapter [Automatic HTML Documentation](https://root.cern/root/htmldoc/guides/users-guide/ROOTUsersGuideChapters/HTMLDoc.pdf), and
/// [THtml's class documentation](https://root.cern/doc/master/classTHtml.html).
///
/// To see this demo script in action start up ROOT and run
/// ~~~{.cpp}
///   root [0] .x $(ROOTSYS)/tutorials/htmlex.C+
/// ~~~
/// and check the output in `./htmldoc`.
///
/// and of course we can put HTML code into comments, too.
///
/// \macro_code
///
/// \author Axel Naumann

#include "THtml.h"

class THtmlDemo: public TObject {
public:
   THtmlDemo(): fHtml(0)
   {
      printf("This class is for demonstration purposes only!\n");
   }
   ~THtmlDemo() { if (fHtml) delete fHtml; }

   // inline methods can have their documentation in front
   // of the declaration. DontDoMuch is so short - where
   // else would one put it?
   void DontDoMuch() {}

   void Convert()
   {
      // Create a "beautified" version of this source file.
      // It will be called htmldoc/htmlex.C.html.

      GetHtml()->SetSourceDir("$(ROOTSYS)/tutorials");
      GetHtml()->Convert("htmlex.C", "Example of THtml", "./htmldoc/", "./");
   }

   void ReferenceDoc()
   {
      // This function documents THtmlDemo.
      // It will create THtmlDemo.html and src/THtmlDemo.cxx.html
      // - the beautified version of the source file

      GetHtml()->SetSourceDir("$(ROOTSYS)/tutorials");
      GetHtml()->SetOutputDir("./htmldoc");
      GetHtml()->MakeIndex("THtmlDemo"); // create ClassIndex.html and the javascript and CSS files
      GetHtml()->MakeClass("THtmlDemo"); // update the class doc
   }

   void MakeDocForAllClasses(Bool_t evenForROOT = kFALSE)
   {
      // Creates the documentation pages for all classes that have
      // been loaded, and that are accessible from "./".
      // If evenForROOT is set, we'll try to document ROOT's classes,
      // too - you will end up with a copy of ROOT's class reference.
      // The documentation will end up in the subdirectory htmldoc/.

      if (evenForROOT)
         GetHtml()->SetSourceDir(".:$(ROOTSYS)");
      else
         GetHtml()->SetSourceDir(".");
      GetHtml()->SetOutputDir("./htmldoc");
      GetHtml()->MakeAll();
   }

   void RunAll() {
      // Show off a bit - do everything we can.
      MakeDocForAllClasses();
      ReferenceDoc();
      Convert();
   }

protected:
   THtml* GetHtml()
   {
      // Return out THtml object, and create it if it doesn't exist.
      if (!fHtml) fHtml = new THtml();
      return fHtml;
   }

private:
   Int_t fVeryUselessMember; // This is a very useless member.
   THtml* fHtml; // our local THtml instance.
   ClassDef(THtmlDemo, 0); // A demo of THtml.
};

void htmlex() {
   THtmlDemo htmldemo;
   htmldemo.RunAll();
}
