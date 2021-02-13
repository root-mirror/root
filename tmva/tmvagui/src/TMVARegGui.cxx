#include "TMVA/TMVARegGui.h"
#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TControlBar.h"
#include "TObjString.h"


// some global lists
static TList*               TMVARegGui_keyContent;
static std::vector<TString> TMVARegGui_inactiveButtons;


TList* TMVA::RegGuiGetKeyList( const TString& pattern )
{
   TList* list = new TList();

   TIter next( TMVARegGui_keyContent );
   TKey* key(0);
   while ((key = (TKey*)next())) {         
      if (TString(key->GetName()).Contains( pattern )) { list->Add( new TObjString( key->GetName() ) ); }
   }
   return list;
}

// utility function
void TMVA::RegGuiActionButton( TControlBar* cbar, 
                               const TString& title, const TString& macro, const TString& comment, 
                               const TString& buttonType, TString requiredKey) 
{
   cbar->AddButton( title, macro, comment, buttonType );

   // search    
   if (requiredKey != "") {
      Bool_t found = kFALSE;
      TIter next( TMVARegGui_keyContent );
      TKey* key(0);
      while ((key = (TKey*)next())) {         
         if (TString(key->GetName()).Contains( requiredKey )) { found = kTRUE; break; }
      }
      if (!found) TMVARegGui_inactiveButtons.push_back( title );
   }
}

// main GUI
void TMVA::TMVARegGui( const char* fName ,TString dataset) 
{   
   // Use this script in order to run the various individual macros
   // that plot the output of TMVA (e.g. running TMVARegression.cxx),
   // stored in the file "TMVA.Regroot"

   TString curMacroPath(gROOT->GetMacroPath());
   // uncomment next line for macros submitted to next root version
   gROOT->SetMacroPath(curMacroPath+":./:$ROOTSYS/tmva/test/:");

   // for the sourceforge version, including $ROOTSYS/tmva/test in the
   // macro path is a mistake, especially if "./" was not part of path
   // add ../macros to the path (comment out next line for the ROOT version of TMVA)
   // gROOT->SetMacroPath(curMacroPath+":../macros:");

   TString curIncludePath=gSystem->GetIncludePath();
   //std::cout <<"inc path="<<curIncludePath<<std::endl;
   TString newIncludePath=TString("-I../ ")+curIncludePath;
   gSystem->SetIncludePath(newIncludePath);
  
   cout << "--- Launch TMVA GUI to view input file: " << fName << endl;

   // init
   TMVARegGui_inactiveButtons.clear();

   // check if file exist
   TFile* file = TFile::Open( fName );
   if (!file) {
      cout << "==> Abort TMVARegGui, please verify filename" << endl;
      return;
   }
   //
   if(file->GetListOfKeys()->GetEntries()<=0)
      {
         cout << "==> Abort TMVARegGui, please verify if dataset exist" << endl;
         return;
      }
   if( (dataset==""||dataset.IsWhitespace()) && (file->GetListOfKeys()->GetEntries()==1))
      {
         TKey *key=(TKey*)file->GetListOfKeys()->At(0);
         dataset=key->GetName();
      }else if((dataset==""||dataset.IsWhitespace()) && (file->GetListOfKeys()->GetEntries()>=1))
      {
         gROOT->Reset();
         gStyle->SetScreenFactor(2); // if you have a large screen, select 1,2 or 1.4
       
         TControlBar *bar=new TControlBar("vertical","Select dataset", 0, 0);
         bar->SetButtonWidth(300);
         for(int i=0;i<file->GetListOfKeys()->GetEntries();i++)
            {
               TKey *key=(TKey*)file->GetListOfKeys()->At(i);
               dataset=key->GetName();
               bar->AddButton(dataset.Data(),Form("TMVA::TMVARegGui(\"%s\",\"%s\")",fName,dataset.Data()),dataset.Data());
            }
       
         bar->AddSeparator();
         bar->AddButton( "Quit",   ".q", "Quit", "button");

         // set the style 
         bar->SetTextColor("black");
         bar->Show();
         gROOT->SaveContext();
         return ;
      }
   // find all references   
   TMVARegGui_keyContent = (TList*)file->GetDirectory(dataset.Data())->GetListOfKeys()->Clone();

   // close file
   file->Close();

   TString defaultRequiredClassifier = "";

   //   gROOT->Reset();
   //   gStyle->SetScreenFactor(2); // if you have a large screen, select 1,2 or 1.4

   // create the control bar
   TControlBar* cbar = new TControlBar( "vertical", "TMVA Plotting Macros for Regression", 0, 0 );

   const TString buttonType( "button" );

   // configure buttons   
   Int_t ic = 1;

   // find all input variables types
   TList* keylist = RegGuiGetKeyList( "InputVariables" );
   TListIter it( keylist );
   TObjString* str = 0;
   char ch = 'a';
   while ( (str = (TObjString*)it()) ) {
      TString tmp   = str->GetString();
      TString title = Form( "Input variables and target(s) '%s'-transformed (training sample)", 
                            tmp.ReplaceAll("InputVariables_","").Data() );
      if (tmp.Contains( "Id" )) title = "Input variables and target(s) (training sample)";
      RegGuiActionButton( cbar, 
                          Form( "    (%i%c) %s    ", ic, ch++, title.Data() ),
                          Form( "TMVA::variables(\"%s\",\"%s\",\"%s\",\"%s\",kTRUE)",dataset.Data() , fName, str->GetString().Data(), title.Data() ),
                          Form( "Plots all '%s'-transformed input variables and target(s) (macro variables(...))", 
                                str->GetString().Data() ),
                          buttonType, str->GetString() );
   }      
   ic++;

   // correlation scatter plots 
   it.Reset(); ch = 'a';
   while ( (str = (TObjString*)it()) ) {
      TString tmp   = str->GetString();
      TString title = Form( "Input variable correlations '%s'-transformed (scatter profiles)", 
                            tmp.ReplaceAll("InputVariables_","").Data() );
      if (tmp.Contains( "Id" )) title = "Input variable correlations (scatter profiles)";
      RegGuiActionButton( cbar, 
                          Form( "(%i%c) %s", ic, ch++, title.Data() ),
                          Form( "TMVA::CorrGui(\"%s\",\"%s\",\"%s\",\"%s\",kTRUE)",dataset.Data() , fName, str->GetString().Data(), title.Data() ),
                          Form( "Plots all correlation profiles between '%s'-transformed input variables (macro CorrGui(...))", 
                                str->GetString().Data() ),
                          buttonType, str->GetString() );
   }      
   
   // coefficients
   RegGuiActionButton( cbar,  
                       Form( "(%i) Input Variable Linear Correlation Coefficients", ++ic ),
                       Form( "TMVA::correlations(\"%s\",\"%s\",kTRUE)",dataset.Data(), fName ),
                       "Plots signal and background correlation summaries for all input variables (macro correlations.cxx)", 
                       buttonType );

   RegGuiActionButton( cbar,  
                       Form( "(%ia) Regression Output Deviation versus Target (test sample)", ++ic ),
                       Form( "TMVA::deviations(\"%s\",\"%s\",TMVA::kMVAType,kTRUE)",dataset.Data(), fName ),
                       "Plots the deviation between regression output and target versus target on test data (macro deviations(...,0))",
                       buttonType, defaultRequiredClassifier );

   RegGuiActionButton( cbar,  
                       Form( "(%ib) Regression Output Deviation versus Target (training sample)", ic ),
                       Form( "TMVA::deviations(\"%s\",\"%s\",TMVA::kCompareType,kTRUE)",dataset.Data() , fName ),
                       "Plots the deviation between regression output and target versus target on test data (macro deviations(...,0))",
                       buttonType, defaultRequiredClassifier );

   RegGuiActionButton( cbar,  
                       Form( "(%ic) Regression Output Deviation versus Input Variables (test sample)", ic ),
                       Form( "TMVA::deviations(\"%s\",\"%s\",TMVA::kMVAType,kFALSE)",dataset.Data(), fName ),
                       "Plots the deviation between regression output and target versus target on test data (macro deviations(...,0))",
                       buttonType, defaultRequiredClassifier );

   RegGuiActionButton( cbar,  
                       Form( "   (%id) Regression Output Deviation versus Input Variables (training sample)   ", ic ),
                       Form( "TMVA::deviations(\"%s\",\"%s\",TMVA::kCompareType,kFALSE)",dataset.Data() , fName ),
                       "Plots the deviation between regression output and target versus target on test data (macro deviations(...,0))",
                       buttonType, defaultRequiredClassifier );

   RegGuiActionButton( cbar,  
                       Form( "(%i) Summary of Average Regression Deviations ", ++ic ),
                       Form( "TMVA::regression_averagedevs(\"%s\",\"%s\")",dataset.Data() , fName ),
                       "Plot Summary of average deviations: MVAvalue - target (macro regression_averagedevs.cxx)",
                       buttonType );

   RegGuiActionButton( cbar,  
                       Form( "(%ia) Network Architecture", ++ic ),
                       Form( "TMVA::network(\"%s\",\"%s\")",dataset.Data(), fName ), 
                       "Plots the MLP weights (macro network.cxx)",
                       buttonType, "MLP" );

   RegGuiActionButton( cbar,  
                       Form( "(%ib) Network Convergence Test", ic ),
                       Form( "TMVA::annconvergencetest(\"%s\",\"%s\")",dataset.Data() , fName ), 
                       "Plots error estimator versus training epoch for training and test samples (macro annconvergencetest.cxx)",
                       buttonType, "MLP" );

   RegGuiActionButton( cbar,  
                       Form( "(%i) Plot Foams", ++ic ),
                       Form("TMVA::PlotFoams(\"%s/weights/TMVARegression_PDEFoam.weights_foams.root\")",dataset.Data()),
                       "Plot Foams (macro PlotFoams.cxx)",
                       buttonType, "PDEFoam" );

   RegGuiActionButton( cbar,  
                       Form( "(%i) Regression Trees (BDT)", ++ic ),
                       Form( "TMVA::BDT_Reg(\"%s\",\"%s\")",dataset.Data() , fName ),
                       "Plots the Regression Trees trained by BDT algorithms (macro BDT_Reg(itree,...))",
                       buttonType, "BDT" );

   RegGuiActionButton( cbar,  
                       Form( "(%i) Regression Tree Control Plots (BDT)", ++ic ),
                       Form( "TMVA::BDTControlPlots(\"%s\",\"%s\")",dataset.Data(), fName ),
                       "Plots to monitor boosting and pruning of regression trees (macro BDTControlPlots.cxx)",
                       buttonType, "BDT" );

   cbar->AddSeparator();

   cbar->AddButton( Form( "(%i) Quit", ++ic ),   ".q", "Quit", buttonType );

   // set the style 
   cbar->SetTextColor("black");

   // there seems to be a bug in ROOT: font jumps back to default after pressing on >2 different buttons
   // cbar->SetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
   
   // draw
   cbar->Show();

   // indicate inactive buttons
   for (UInt_t i=0; i<TMVARegGui_inactiveButtons.size(); i++) cbar->SetButtonState( TMVARegGui_inactiveButtons[i], 3 );
   if (TMVARegGui_inactiveButtons.size() > 0) {
      cout << "=== Note: inactive buttons indicate that the corresponding methods were not trained ===" << endl;
   }

   gROOT->SaveContext();
}
