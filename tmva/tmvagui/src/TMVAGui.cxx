#include "TMVA/TMVAGui.h"
#include "TMVA/Types.h"
#include <iostream>
#include <vector>

#include "TList.h"
#include "TROOT.h"
#include "TKey.h"
#include "TString.h"
#include "TControlBar.h"
#include "TObjString.h"


// some global lists
static TList*               TMVAGui_keyContent;
static std::vector<TString> TMVAGui_inactiveButtons;

TList* TMVA::GetKeyList( const TString& pattern )
{
   TList* list = new TList();

   TIter next( TMVAGui_keyContent );
   TKey* key(0);
   while ((key = (TKey*)next())) {
      if (TString(key->GetName()).Contains( pattern )) { list->Add( new TObjString( key->GetName() ) ); }
   }
   return list;
}

// utility function
void TMVA::ActionButton( TControlBar* cbar,
                         const TString& title, const TString& macro, const TString& comment,
                         const TString& buttonType, TString requiredKey )
{
   cbar->AddButton( title, macro, comment, buttonType );

   // search
   if (requiredKey != "") {
      Bool_t found = kFALSE;
      TIter next( TMVAGui_keyContent );
      TKey* key(0);
      while ((key = (TKey*)next())) {
         if (TString(key->GetName()).Contains( requiredKey )) { found = kTRUE; break; }
      }
      if (!found) TMVAGui_inactiveButtons.push_back( title );
   }
}

// main GUI
void TMVA::TMVAGui( const char* fName  , TString dataset)
{
   // Use this script in order to run the various individual macros
   // that plot the output of TMVA (e.g. running TMVAClassification.C),
   // stored in the file "TMVA.root"

   TString curMacroPath(gROOT->GetMacroPath());
   // uncomment next line for macros submitted to next root version
   gROOT->SetMacroPath(curMacroPath+":./:$ROOTSYS/tmva/test/:");

   // for the sourceforge version, including $ROOTSYS/tmva/test in the
   // macro path is a mistake, especially if "./" was not part of path
   // add ../macros to the path (comment out next line for the ROOT version of TMVA)
   // gROOT->SetMacroPath(curMacroPath+":../macros:");

   TString curIncludePath=gSystem->GetIncludePath();
   TString newIncludePath=TString("-I../ ")+curIncludePath;
   gSystem->SetIncludePath(newIncludePath);

   cout << "--- Launch TMVA GUI to view input file: " << fName << endl;

   // init
   TMVAGui_inactiveButtons.clear();

   // check if file exist
   TFile* file = TFile::Open( fName );
   if (!file) {
      cout << "==> Abort TMVAGui, please verify filename" << endl;
      return;
   }

   if(file->GetListOfKeys()->GetEntries()<=0)
      {
         cout << "==> Abort TMVAGui, please verify if dataset exist" << endl;
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
               bar->AddButton(dataset.Data(),Form("TMVA::TMVAGui(\"%s\",\"%s\")",fName,dataset.Data()),dataset.Data());
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
   TMVAGui_keyContent = (TList*)file->GetDirectory(dataset.Data())->GetListOfKeys()->Clone();

   // close file
   file->Close();

   TString defaultRequiredClassifier = "";

   //   gROOT->Reset();
   //   gStyle->SetScreenFactor(2); // if you have a large screen, select 1,2 or 1.4

   // create the control bar
   TControlBar* cbar = new TControlBar( "vertical", "TMVA Plotting Macros for Classification", 0, 0 );

   const TString buttonType( "button" );

   // configure buttons
   Int_t ic = 1;

   // find all input variables types
   TList* keylist = GetKeyList( "InputVariables" );
   TListIter it( keylist );
   TObjString* str = 0;
   char ch = 'a';
   while ((str = (TObjString*)it())) {
      TString tmp   = str->GetString();
      TString title = Form( "Input variables '%s'-transformed (training sample)",
                            tmp.ReplaceAll("InputVariables_","").Data() );
      if (tmp.Contains( "Id" )) title = "Input variables (training sample)";
      ActionButton( cbar,
                    Form( "(%i%c) %s", ic, ch++, title.Data() ),
                    Form( "TMVA::variables(\"%s\",\"%s\",\"%s\",\"%s\")",dataset.Data(), fName, str->GetString().Data(), title.Data() ),
                    Form( "Plots all '%s'-transformed input variables (macro variables(...))", str->GetString().Data() ),
                    buttonType, str->GetString() );
   }
   ic++;

   // correlation scatter plots
   it.Reset(); ch = 'a';
   while ((str = (TObjString*)it())) {
      TString tmp   = str->GetString();
      TString title = Form( "Input variable correlations '%s'-transformed (scatter profiles)",
                            tmp.ReplaceAll("InputVariables_","").Data() );
      if (tmp.Contains( "Id" )) title = "Input variable correlations (scatter profiles)";
      ActionButton( cbar,
                    Form( "(%i%c) %s", ic, ch++, title.Data() ),
                    Form( "TMVA::CorrGui(\"%s\",\"%s\",\"%s\",\"%s\")",dataset.Data(), fName, str->GetString().Data(), title.Data() ),
                    Form( "Plots all correlation profiles between '%s'-transformed input variables (macro CorrGui(...))",
                          str->GetString().Data() ),
                    buttonType, str->GetString() );
   }

   TString title;
   // coefficients
   title =Form( "(%i) Input Variable Linear Correlation Coefficients", ++ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::correlations(\"%s\",\"%s\")",dataset.Data(), fName ),
                 "Plots signal and background correlation summaries for all input variables (macro correlations.C)",
                 buttonType );

   title =Form( "(%ia) Classifier Output Distributions (test sample)", ++ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::mvas(\"%s\",\"%s\",  TMVA::kMVAType)",dataset.Data(), fName ),
                 "Plots the output of each classifier for the test data (macro mvas(...,0))",
                 buttonType, defaultRequiredClassifier );

   title =Form( "(%ib) Classifier Output Distributions (test and training samples superimposed)", ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::mvas(\"%s\",\"%s\",  TMVA::kCompareType )",dataset.Data(), fName),
                 "Plots the output of each classifier for the test (histograms) and training (dots) data (macro mvas(...,3))",
                 buttonType, defaultRequiredClassifier );

   title = Form( "(%ic) Classifier Probability Distributions (test sample)", ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::mvas(\"%s\",\"%s\", TMVA::kProbaType)",dataset.Data(), fName ),
                 "Plots the probability of each classifier for the test data (macro mvas(...,1))",
                 buttonType, defaultRequiredClassifier );

   title =Form( "(%id) Classifier Rarity Distributions (test sample)", ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::mvas(\"%s\",\"%s\", TMVA::kRarityType)",dataset.Data(), fName ),
                 "Plots the Rarity of each classifier for the test data (macro mvas(...,2)) - background distribution should be uniform",
                 buttonType, defaultRequiredClassifier );

   title =Form( "(%ia) Classifier Cut Efficiencies", ++ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::mvaeffs(\"%s\",\"%s\")",dataset.Data(), fName ),
                 "Plots signal and background efficiencies versus cut on classifier output (macro mvaeffs.cxx)",
                 buttonType, defaultRequiredClassifier );

   title = Form( "(%ib) Classifier Background Rejection vs Signal Efficiency (ROC curve)", ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::efficiencies(\"%s\",\"%s\")",dataset.Data(), fName ),
                 "Plots background rejection vs signal efficiencies (macro efficiencies.cxx) [\"ROC\" stands for \"Receiver Operation Characteristics\"]",
                 buttonType, defaultRequiredClassifier );

   title = Form( "(%ib) Classifier 1/(Backgr. Efficiency) vs Signal Efficiency (ROC curve)", ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::efficiencies(\"%s\",\"%s\",%d)",dataset.Data(), fName, 3 ),
                 "Plots 1/(background eff.)  vs signal efficiencies (macro efficiencies.cxx) [\"ROC\" stands for \"Receiver Operation Characteristics\"]",
                 buttonType, defaultRequiredClassifier );

   title = Form( "(%i) Parallel Coordinates (requires ROOT-version >= 5.17)", ++ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::paracoor(\"%s\",\"%s\")",dataset.Data(), fName ),
                 "Plots parallel coordinates for classifiers and input variables (macro paracoor.cxx, requires ROOT >= 5.17)",
                 buttonType, defaultRequiredClassifier );

   title =Form( "(%i) PDFs of Classifiers (requires \"CreateMVAPdfs\" option set)", ++ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::probas(\"%s\",\"%s\")",dataset.Data(), fName ),
                 "Plots the PDFs of the classifier output distributions for signal and background - if requested (macro probas.cxx)",
                 buttonType, defaultRequiredClassifier );

   title = Form( "(%i) Training History", ++ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::training_history(\"%s\",\"%s\")",dataset.Data(), fName ),
                 "Plot training history of classifiers with multiple passed (eg Neural Networks) ",
                 buttonType, defaultRequiredClassifier );

   title = Form( "(%i) Likelihood Reference Distributiuons", ++ic);
   ActionButton( cbar,
                 title,
                 Form( "TMVA::likelihoodrefs(\"%s\",\"%s\")",dataset.Data(), fName ),
                 "Plots to verify the likelihood reference distributions (macro likelihoodrefs.cxx)",
                 buttonType, "Likelihood" );

   title = Form( "(%ia) Network Architecture (MLP)", ++ic );
   TString call = Form( "TMVA::network(\"%s\",\"%s\")", dataset.Data(),fName );
   ActionButton( cbar,
                 title,
                 call,
                 "Plots the MLP weights (macro network.cxx)",
                 buttonType, "MLP" );

   title = Form( "(%ib) Network Convergence Test (MLP)", ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::annconvergencetest(\"%s\",\"%s\")",dataset.Data(), fName ),
                 "Plots error estimator versus training epoch for training and test samples (macro annconvergencetest.C)",
                 buttonType, "MLP" );

   title = Form( "(%i) Decision Trees (BDT)", ++ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::BDT(\"%s\",\"%s\")",dataset.Data() , fName ),
                 "Plots the Decision Trees trained by BDT algorithms (macro BDT(itree,...))",
                 buttonType, "BDT" );

   title = Form( "(%i) Decision Tree Control Plots (BDT)", ++ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::BDTControlPlots(\"%s\",\"%s\")",dataset.Data() , fName ),
                 "Plots to monitor boosting and pruning of decision trees (macro BDTControlPlots.cxx)",
                 buttonType, "BDT" );
   //    ActionButton( cbar,
   //                  Form( "(%i) Rule Ensemble Importance Plots (RuleFit)", ++ic ),
   //                  Form( "TMVA::rulevis(\"%s\",0)", fName ),
   //                  "Plots all input variables with rule ensemble weights, including linear terms (macro rulevis.cxx)",
   //                  buttonType, "RuleFit" );

   title = Form( "(%i) Plot Foams (PDEFoam)", ++ic );
   ActionButton( cbar,
                 title,
                 Form("TMVA::PlotFoams(\"%s/weights/TMVAClassification_PDEFoam.weights_foams.root\")",dataset.Data()),
                 "Plot Foams (macro PlotFoams.cxx)",
                 buttonType, "PDEFoam" );

   title = Form( "(%i) General Boost Control Plots", ++ic );
   ActionButton( cbar,
                 title,
                 Form( "TMVA::BoostControlPlots(\"%s\",\"%s\")",dataset.Data() , fName ),
                 "Plots to monitor boosting of general classifiers (macro BoostControlPlots)",
                 buttonType, "Boost" );

   cbar->AddSeparator();

   cbar->AddButton( Form( "(%i) Quit", ++ic ),   ".q", "Quit", buttonType );

   // set the style
   cbar->SetTextColor("black");

   // there seems to be a bug in ROOT: font jumps back to default after pressing on >2 different buttons
   // cbar->SetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");

   // draw
   cbar->Show();

   // indicate inactive buttons
   for (UInt_t i=0; i<TMVAGui_inactiveButtons.size(); i++) cbar->SetButtonState(TMVAGui_inactiveButtons[i], 3 );
   if (TMVAGui_inactiveButtons.size() > 0) {
      cout << "=== Note: inactive buttons indicate classifiers that were not trained, ===" << endl;
      cout << "===       or functionalities that were not invoked during the training ===" << endl;
   }

   gROOT->SaveContext();
}
