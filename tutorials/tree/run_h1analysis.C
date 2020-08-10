/// \file
/// \ingroup tutorial_tree
/// \notebook -nodraw
/// Macro driving the analysis can specify file name and type
///
/// .- type == 0 : normal
///  - type =  1 : use AClic to compile selector
///  - type =  2 : use a fill list and then process the fill list///
///
/// \macro_code
///
/// \author

//----------------------------------------

void run_h1analysis(int type = 0, const char * h1dir = 0) {

   std::cout << "Run h1 analysis " << std::endl;

   // create first the chain with all the files

   TChain chain("h42");

   if (h1dir) {
      gSystem->Setenv("H1",h1dir);
   }
   else
      gSystem->Setenv("H1","http://root.cern/files/h1/");


   std::cout << "Creating the chain" << std::endl;

   chain.SetCacheSize(20*1024*1024);
   chain.Add("$H1/dstarmb.root");
   chain.Add("$H1/dstarp1a.root");
   chain.Add("$H1/dstarp1b.root");
   chain.Add("$H1/dstarp2.root");

   TString selectionMacro = gSystem->GetDirName(__FILE__) + "/h1analysis.C";

   if (type == 0)
      chain.Process(selectionMacro);
   else if (type == 1)   {
      // use AClic ( add a + at the end
      selectionMacro += "+";
      chain.Process(selectionMacro);
   }
   else if (type == 2) {
      chain.Process(selectionMacro,"fillList");
      chain.Process(selectionMacro,"useList");
   }
}



