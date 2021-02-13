/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example on how to use the trained multiclass
/// classifiers within an analysis module
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TMVAMulticlassApplication
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker


#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TH1F.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace TMVA;

void TMVAMulticlassApplication( TString myMethodList = "" )
{

   TMVA::Tools::Instance();

   //---------------------------------------------------------------
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;
   Use["MLP"]             = 1;
   Use["BDTG"]            = 1;
   Use["DL_CPU"]          = 1;
   Use["DL_GPU"]          = 1;
   Use["FDA_GA"]          = 1;
   Use["PDEFoam"]         = 1;
   //---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAMulticlassApp" << std::endl;
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " " << std::endl;
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }


   // create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // create a set of variables and declare them to the reader
   // - the variable names must corresponds in name and type to
   // those given in the weight file(s) that you use
   Float_t var1, var2, var3, var4;
   reader->AddVariable( "var1", &var1 );
   reader->AddVariable( "var2", &var2 );
   reader->AddVariable( "var3", &var3 );
   reader->AddVariable( "var4", &var4 );

   // book the MVA methods
   TString dir    = "dataset/weights/";
   TString prefix = "TMVAMulticlass";

   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
        TString methodName = TString(it->first) + TString(" method");
        TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
        // check if file existing (i.e. method has been trained)
        if (!gSystem->AccessPathName( weightfile ))
           // file exists
           reader->BookMVA( methodName, weightfile );
        else {
           std::cout << "TMVAMultiClassApplication: Skip " << methodName << " since it has not been trained !" << std::endl;
           it->second = 0;
        }
      }
   }

   // book output histograms
   UInt_t nbin = 100;
   TH1F *histMLP_signal(0), *histBDTG_signal(0), *histFDAGA_signal(0), *histPDEFoam_signal(0);
   TH1F *histDLCPU_signal(0), *histDLGPU_signal(0); 
   if (Use["MLP"])
      histMLP_signal    = new TH1F( "MVA_MLP_signal",    "MVA_MLP_signal",    nbin, 0., 1.1 );
   if (Use["BDTG"])
      histBDTG_signal  = new TH1F( "MVA_BDTG_signal",   "MVA_BDTG_signal",   nbin, 0., 1.1 );
   if (Use["DL_CPU"])
      histDLCPU_signal = new TH1F("MVA_DLCPU_signal", "MVA_DLCPU_signal", nbin, 0., 1.1);
   if (Use["DL_GPU"])
      histDLGPU_signal = new TH1F("MVA_DLGPU_signal", "MVA_DLGPU_signal", nbin, 0., 1.1);
   if (Use["FDA_GA"])
      histFDAGA_signal = new TH1F( "MVA_FDA_GA_signal", "MVA_FDA_GA_signal", nbin, 0., 1.1 );
   if (Use["PDEFoam"])
      histPDEFoam_signal = new TH1F( "MVA_PDEFoam_signal", "MVA_PDEFoam_signal", nbin, 0., 1.1 );


   TFile *input(0);
   TString fname = "./tmva_example_multiclass.root";
   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   }
   else {
      TFile::SetCacheFileDir(".");
      input = TFile::Open("http://root.cern.ch/files/tmva_multiclass_example.root", "CACHEREAD");
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAMulticlassApp : Using input file: " << input->GetName() << std::endl;

   // prepare the tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop

   TTree* theTree = (TTree*)input->Get("TreeS");
   std::cout << "--- Select signal sample" << std::endl;
   theTree->SetBranchAddress( "var1", &var1 );
   theTree->SetBranchAddress( "var2", &var2 );
   theTree->SetBranchAddress( "var3", &var3 );
   theTree->SetBranchAddress( "var4", &var4 );

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
      if (ievt%1000 == 0){
         std::cout << "--- ... Processing event: " << ievt << std::endl;
      }
      theTree->GetEntry(ievt);

      if (Use["MLP"])
         histMLP_signal->Fill((reader->EvaluateMulticlass( "MLP method" ))[0]);
      if (Use["BDTG"])
         histBDTG_signal->Fill((reader->EvaluateMulticlass( "BDTG method" ))[0]);
      if (Use["DL_CPU"])
         histDLCPU_signal->Fill((reader->EvaluateMulticlass("DL_CPU method"))[0]);
      if (Use["DL_GPU"])
         histDLGPU_signal->Fill((reader->EvaluateMulticlass("DL_GPU method"))[0]);
      if (Use["FDA_GA"])
         histFDAGA_signal->Fill((reader->EvaluateMulticlass( "FDA_GA method" ))[0]);
      if (Use["PDEFoam"])
         histPDEFoam_signal->Fill((reader->EvaluateMulticlass( "PDEFoam method" ))[0]);

   }

   // get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   TFile *target  = new TFile( "TMVAMulticlassApp.root","RECREATE" );
   if (Use["MLP"])
      histMLP_signal->Write();
   if (Use["BDTG"])
      histBDTG_signal->Write();
   if (Use["DL_CPU"])
      histDLCPU_signal->Write();
   if (Use["DL_GPU"])
      histDLGPU_signal->Write();
   if (Use["FDA_GA"])
      histFDAGA_signal->Write();
   if (Use["PDEFoam"])
      histPDEFoam_signal->Write();

   target->Close();
   std::cout << "--- Created root file: \"TMVMulticlassApp.root\" containing the MVA output histograms" << std::endl;

   delete reader;

   std::cout << "==> TMVAMulticlassApp is done!" << std::endl << std::endl;
}

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   TMVAMulticlassApplication(methodList);
   return 0;
}
