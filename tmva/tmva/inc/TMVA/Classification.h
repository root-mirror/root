// @(#)root/tmva $Id$ 2017
// Authors:  Omar Zapata, Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss, Kai Voss, Eckhard von Toerne,
// Jan Therhaag

#ifndef ROOT_TMVA_Classification
#define ROOT_TMVA_Classification

#include <TString.h>
#include <TMultiGraph.h>
#include <vector>
#include <map>

#include <TMVA/IMethod.h>
#include <TMVA/MethodBase.h>
#include <TMVA/Configurable.h>
#include <TMVA/Types.h>
#include <TMVA/DataSet.h>
#include <TMVA/Event.h>
#include <TMVA/Results.h>
#include <TMVA/ResultsClassification.h>
#include <TMVA/ResultsMulticlass.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/OptionMap.h>
#include <TMVA/Envelope.h>

/*! \class TMVA::ClassificationResult
 * Class to save the results of the classifier.
 * Every machine learning method booked have an object for the results
 * in the classification process, in this class is stored the mvas,
 * data loader name and ml method name and title.
 * You can to display the resutls calling the method Show, get the ROC-integral with the
 * method GetROCIntegral or get the TMVA::ROCCurve object calling GetROC.
\ingroup TMVA
*/

/*! \class TMVA::Classification
 * Class to perform two class classification.
 * The first step before any analysis is to preperate the data,
 * to do that you need to create an object of TMVA::DataLoader,
 * in this object you need to configure the variables and the number of events
 * to train/test.
 * The class TMVA::Experimental::Classification needs a TMVA::DataLoader object,
 * optional a TFile object to save the results and some extra options in a string
 * like "V:Color:Transformations=I;D;P;U;G:Silent:DrawProgressBar:ModelPersistence:Jobs=2" where:
 * V                = verbose output
 * Color            = coloured screen output
 * Silent           = batch mode: boolean silent flag inhibiting any output from TMVA
 * Transformations  = list of transformations to test.
 * DrawProgressBar  = draw progress bar to display training and testing.
 * ModelPersistence = to save the trained model in xml or serialized files.
 * Jobs             = number of ml methods to test/train in parallel using MultiProc, requires to call Evaluate method.
 * Basic example.
 * \code
void classification(UInt_t jobs = 2)
{
   TMVA::Tools::Instance();

   TFile *input(0);
   TString fname = "./tmva_class_example.root";
   if (!gSystem->AccessPathName(fname)) {
      input = TFile::Open(fname); // check if file in local directory exists
   } else {
      TFile::SetCacheFileDir(".");
      input = TFile::Open("http://root.cern.ch/files/tmva_class_example.root", "CACHEREAD");
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }

   // Register the training and test trees

   TTree *signalTree = (TTree *)input->Get("TreeS");
   TTree *background = (TTree *)input->Get("TreeB");

   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

   dataloader->AddVariable("myvar1 := var1+var2", 'F');
   dataloader->AddVariable("myvar2 := var1-var2", "Expression 2", "", 'F');
   dataloader->AddVariable("var3", "Variable 3", "units", 'F');
   dataloader->AddVariable("var4", "Variable 4", "units", 'F');

   dataloader->AddSpectator("spec1 := var1*2", "Spectator 1", "units", 'F');
   dataloader->AddSpectator("spec2 := var1*3", "Spectator 2", "units", 'F');

   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight = 1.0;
   Double_t backgroundWeight = 1.0;

   dataloader->SetBackgroundWeightExpression("weight");

   TMVA::Experimental::Classification *cl = new TMVA::Experimental::Classification(dataloader, Form("Jobs=%d", jobs));

   cl->BookMethod(TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=2000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:"
                                             "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");
   cl->BookMethod(TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm");

   cl->Evaluate(); // Train and Test all methods

   auto &results = cl->GetResults();

   TCanvas *c = new TCanvas(Form("ROC"));
   c->SetTitle("ROC-Integral Curve");

   auto mg = new TMultiGraph();
   for (UInt_t i = 0; i < results.size(); i++) {
      auto roc = results[i].GetROCGraph();
      roc->SetLineColorAlpha(i + 1, 0.1);
      mg->Add(roc);
   }
   mg->Draw("AL");
   mg->GetXaxis()->SetTitle(" Signal Efficiency ");
   mg->GetYaxis()->SetTitle(" Background Rejection ");
   c->BuildLegend(0.15, 0.15, 0.3, 0.3);
   c->Draw();

   delete cl;
}
 * \endcode
 *
\ingroup TMVA
*/

namespace TMVA {
class ResultsClassification;
namespace Experimental {
class ClassificationResult : public TObject {
   friend class Classification;

private:
   OptionMap fMethod;                                                             //
   TString fDataLoaderName;                                                       //
   std::map<UInt_t, std::vector<std::tuple<Float_t, Float_t, Bool_t>>> fMvaTrain; // Mvas for two-class classification
   std::map<UInt_t, std::vector<std::tuple<Float_t, Float_t, Bool_t>>>
      fMvaTest;                      // Mvas for two-class and multiclass classification
   std::vector<TString> fClassNames; //

   Bool_t IsMethod(TString methodname, TString methodtitle);
   Bool_t fIsCuts;        // if it is a method cuts need special output
   Double_t fROCIntegral; //

public:
   ClassificationResult();
   ClassificationResult(const ClassificationResult &cr);
   ~ClassificationResult() {}

   const TString GetMethodName() const { return fMethod.GetValue<TString>("MethodName"); }
   const TString GetMethodTitle() const { return fMethod.GetValue<TString>("MethodTitle"); }
   ROCCurve *GetROC(UInt_t iClass = 0, TMVA::Types::ETreeType type = TMVA::Types::kTesting);
   Double_t GetROCIntegral(UInt_t iClass = 0, TMVA::Types::ETreeType type = TMVA::Types::kTesting);
   TString GetDataLoaderName() { return fDataLoaderName; }
   Bool_t IsCutsMethod() { return fIsCuts; }

   void Show();

   TGraph *GetROCGraph(UInt_t iClass = 0, TMVA::Types::ETreeType type = TMVA::Types::kTesting);
   ClassificationResult &operator=(const ClassificationResult &r);

   ClassDef(ClassificationResult, 3);
};

class Classification : public Envelope {
   std::vector<ClassificationResult> fResults; //!
   std::vector<IMethod *> fIMethods;           //! vector of objects with booked methods
   Types::EAnalysisType fAnalysisType;         //!
   Bool_t fCorrelations;                       //!
   Bool_t fROC;                                //!
public:
   explicit Classification(DataLoader *loader, TFile *file, TString options);
   explicit Classification(DataLoader *loader, TString options);
   ~Classification();

   virtual void Train();
   virtual void TrainMethod(TString methodname, TString methodtitle);
   virtual void TrainMethod(Types::EMVA method, TString methodtitle);

   virtual void Test();
   virtual void TestMethod(TString methodname, TString methodtitle);
   virtual void TestMethod(Types::EMVA method, TString methodtitle);

   virtual void Evaluate();

   std::vector<ClassificationResult> &GetResults();

   MethodBase *GetMethod(TString methodname, TString methodtitle);

protected:
   TString GetMethodOptions(TString methodname, TString methodtitle);
   Bool_t HasMethodObject(TString methodname, TString methodtitle, Int_t &index);
   Bool_t IsCutsMethod(TMVA::MethodBase *method);
   TMVA::ROCCurve *
   GetROC(TMVA::MethodBase *method, UInt_t iClass = 0, TMVA::Types::ETreeType type = TMVA::Types::kTesting);
   TMVA::ROCCurve *GetROC(TString methodname, TString methodtitle, UInt_t iClass = 0,
                          TMVA::Types::ETreeType type = TMVA::Types::kTesting);

   Double_t GetROCIntegral(TString methodname, TString methodtitle, UInt_t iClass = 0);

   ClassificationResult &GetResults(TString methodname, TString methodtitle);
   void CopyFrom(TDirectory *src, TFile *file);
   void MergeFiles();

   ClassDef(Classification, 0);
};
} // namespace Experimental
} // namespace TMVA

#endif // ROOT_TMVA_Classification
