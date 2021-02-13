// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Eckhard von Toerne

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodCompositeBase                                                   *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Virtual base class for all MVA method                                     *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker   <Andreas.Hocker@cern.ch>   - CERN, Switzerland          *
 *      Nadim Sah         <Nadim.Sah@cern.ch>        - Berlin, Germany            *
 *      Peter Speckmayer  <Peter.Speckmazer@cern.ch> - CERN, Switzerland          *
 *      Joerg Stelzer     <Joerg.Stelzer@cern.ch>    - MSU East Lansing, USA      *
 *      Helge Voss        <Helge.Voss@cern.ch>       - MPI-K Heidelberg, Germany  *
 *      Jan Therhaag      <Jan.Therhaag@cern.ch>     - U of Bonn, Germany         *
 *      Eckhard v. Toerne <evt@uni-bonn.de>          - U of Bonn, Germany         *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      MSU East Lansing, USA                                                     *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodCategory
\ingroup TMVA

Class for categorizing the phase space

This class is meant to allow categorisation of the data. For different
categories, different classifiers may be booked and different variables
may be considered. The aim is to account for the difference that
is due to different locations/angles.
*/


#include "TMVA/MethodCategory.h"

#include <algorithm>
#include <vector>

#include "TRandom3.h"
#include "TH1F.h"
#include "TSpline.h"
#include "TDirectory.h"
#include "TTreeFormula.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Config.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/DataSetManager.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodCompositeBase.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/PDF.h"
#include "TMVA/Ranking.h"
#include "TMVA/Timer.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/VariableInfo.h"
#include "TMVA/VariableRearrangeTransform.h"

REGISTER_METHOD(Category)

ClassImp(TMVA::MethodCategory);

////////////////////////////////////////////////////////////////////////////////
/// standard constructor

   TMVA::MethodCategory::MethodCategory( const TString& jobName,
                                         const TString& methodTitle,
                                         DataSetInfo& theData,
                                         const TString& theOption )
   : TMVA::MethodCompositeBase( jobName, Types::kCategory, methodTitle, theData, theOption),
   fCatTree(0),
   fDataSetManager(NULL)
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor from weight file

TMVA::MethodCategory::MethodCategory( DataSetInfo& dsi,
                                      const TString& theWeightFile)
   : TMVA::MethodCompositeBase( Types::kCategory, dsi, theWeightFile),
     fCatTree(0),
     fDataSetManager(NULL)
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodCategory::~MethodCategory( void )
{
   std::vector<TTreeFormula*>::iterator formIt = fCatFormulas.begin();
   std::vector<TTreeFormula*>::iterator lastF = fCatFormulas.end();
   for(;formIt!=lastF; ++formIt) delete *formIt;
   delete fCatTree;
}

////////////////////////////////////////////////////////////////////////////////
/// check whether method category has analysis type
/// the method type has to be the same for all sub-methods

Bool_t TMVA::MethodCategory::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets )
{
   std::vector<IMethod*>::iterator itrMethod = fMethods.begin();

   // iterate over methods and check whether they have the analysis type
   for(; itrMethod != fMethods.end(); ++itrMethod ) {
      if ( !(*itrMethod)->HasAnalysisType(type, numberClasses, numberTargets) )
         return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// options for this method

void TMVA::MethodCategory::DeclareOptions()
{
}

////////////////////////////////////////////////////////////////////////////////
/// adds sub-classifier for a category

TMVA::IMethod* TMVA::MethodCategory::AddMethod( const TCut& theCut,
                                                const TString& theVariables,
                                                Types::EMVA theMethod ,
                                                const TString& theTitle,
                                                const TString& theOptions )
{
   std::string addedMethodName(Types::Instance().GetMethodName(theMethod).Data());

   Log() << kINFO << "Adding sub-classifier: " << addedMethodName << "::" << theTitle << Endl;

   DataSetInfo& dsi = CreateCategoryDSI(theCut, theVariables, theTitle);

   IMethod* addedMethod = ClassifierFactory::Instance().Create(addedMethodName,GetJobName(),theTitle,dsi,theOptions);

   MethodBase *method = (dynamic_cast<MethodBase*>(addedMethod));
   if(method==0) return 0;

   if(fModelPersistence) method->SetWeightFileDir(fFileDir);
   method->SetModelPersistence(fModelPersistence);
   method->SetAnalysisType( fAnalysisType );
   method->SetupMethod();
   method->ParseOptions();
   method->ProcessSetup();
   method->SetFile(fFile);
   method->SetSilentFile(IsSilentFile());


   // set or create correct method base dir for added method
   const TString dirName(Form("Method_%s",method->GetMethodTypeName().Data()));
   TDirectory * dir = BaseDir()->GetDirectory(dirName);
   if (dir != 0) method->SetMethodBaseDir( dir );
   else method->SetMethodBaseDir( BaseDir()->mkdir(dirName,Form("Directory for all %s methods", method->GetMethodTypeName().Data())) );

   // method->SetBaseDir(eigenes base dir, gucken ob Fisher dir existiert, sonst erzeugen )

   // check-for-unused-options is performed; may be overridden by derived
   // classes
   method->CheckSetup();

   // disable writing of XML files and standalone classes for sub methods
   method->DisableWriting( kTRUE );

   // store method, cut and variable names and create cut formula
   fMethods.push_back(method);
   fCategoryCuts.push_back(theCut);
   fVars.push_back(theVariables);

   DataSetInfo& primaryDSI = DataInfo();

   UInt_t newSpectatorIndex = primaryDSI.GetSpectatorInfos().size();
   fCategorySpecIdx.push_back(newSpectatorIndex);

   primaryDSI.AddSpectator( Form("%s_cat%i:=%s", GetName(),(int)fMethods.size(),theCut.GetTitle()),
                            Form("%s:%s",GetName(),method->GetName()),
                            "pass", 0, 0, 'C' );

   return method;
}

////////////////////////////////////////////////////////////////////////////////
/// create a DataSetInfo object for a sub-classifier

TMVA::DataSetInfo& TMVA::MethodCategory::CreateCategoryDSI(const TCut& theCut,
                                                           const TString& theVariables,
                                                           const TString& theTitle)
{
   // create a new dsi with name: theTitle+"_dsi"
   TString dsiName=theTitle+"_dsi";
   DataSetInfo& oldDSI = DataInfo();
   DataSetInfo* dsi = new DataSetInfo(dsiName);

   // register the new dsi
   //   DataSetManager::Instance().AddDataSetInfo(*dsi); // DSMTEST replaced by following line
   fDataSetManager->AddDataSetInfo(*dsi);

   // copy the targets and spectators from the old dsi to the new dsi
   std::vector<VariableInfo>::iterator itrVarInfo;

   for (itrVarInfo = oldDSI.GetTargetInfos().begin(); itrVarInfo != oldDSI.GetTargetInfos().end(); ++itrVarInfo)
      dsi->AddTarget(*itrVarInfo);

   for (itrVarInfo = oldDSI.GetSpectatorInfos().begin(); itrVarInfo != oldDSI.GetSpectatorInfos().end(); ++itrVarInfo)
      dsi->AddSpectator(*itrVarInfo);

   // split string that contains the variables into tiny little pieces
   std::vector<TString> variables = gTools().SplitString(theVariables,':' );

   // prepare to create varMap
   std::vector<UInt_t> varMap;
   Int_t counter=0;

   // add the variables that were specified in theVariables
   std::vector<TString>::iterator itrVariables;
   Bool_t found = kFALSE;

   // iterate over all variables in 'variables' and add them
   for (itrVariables = variables.begin(); itrVariables != variables.end(); ++itrVariables) {
      counter=0;

      // check the variables of the old dsi for the variable that we want to add
      for (itrVarInfo = oldDSI.GetVariableInfos().begin(); itrVarInfo != oldDSI.GetVariableInfos().end(); ++itrVarInfo) {
         if((*itrVariables==itrVarInfo->GetLabel()) ) { // || (*itrVariables==itrVarInfo->GetExpression())) {
            // don't compare the expression, since the user might take two times the same expression, but with different labels
            // and apply different transformations to the variables.
            dsi->AddVariable(*itrVarInfo);
            varMap.push_back(counter);
            found = kTRUE;
         }
         counter++;
      }

      // check the spectators of the old dsi for the variable that we want to add
      for (itrVarInfo = oldDSI.GetSpectatorInfos().begin(); itrVarInfo != oldDSI.GetSpectatorInfos().end(); ++itrVarInfo) {
         if((*itrVariables==itrVarInfo->GetLabel()) ) { // || (*itrVariables==itrVarInfo->GetExpression())) {
            // don't compare the expression, since the user might take two times the same expression, but with different labels
            // and apply different transformations to the variables.
            dsi->AddVariable(*itrVarInfo);
            varMap.push_back(counter);
            found = kTRUE;
         }
         counter++;
      }

      // if the variable is neither in the variables nor in the spectators, we abort
      if (!found) {
         Log() << kFATAL <<"The variable " << itrVariables->Data() << " was not found and could not be added " << Endl;
      }
      found = kFALSE;
   }

   // in the case that no variables are specified, add the default-variables from the original dsi
   if (theVariables=="") {
      for (UInt_t i=0; i<oldDSI.GetVariableInfos().size(); i++) {
         dsi->AddVariable(oldDSI.GetVariableInfos()[i]);
         varMap.push_back(i);
      }
   }

   // add the variable map 'varMap' to the vector of varMaps
   fVarMaps.push_back(varMap);

   // set classes and cuts
   UInt_t nClasses=oldDSI.GetNClasses();
   TString className;

   for (UInt_t i=0; i<nClasses; i++) {
      className = oldDSI.GetClassInfo(i)->GetName();
      dsi->AddClass(className);
      dsi->SetCut(oldDSI.GetCut(i),className);
      dsi->AddCut(theCut,className);
      dsi->SetWeightExpression(oldDSI.GetWeightExpression(i),className);
   }

   // set split options, root dir and normalization for the new dsi
   dsi->SetSplitOptions(oldDSI.GetSplitOptions());
   dsi->SetRootDir(oldDSI.GetRootDir());
   TString norm(oldDSI.GetNormalization().Data());
   dsi->SetNormalization(norm);

   DataSetInfo& dsiReference= (*dsi);

   return dsiReference;
}

////////////////////////////////////////////////////////////////////////////////
/// initialize the method

void TMVA::MethodCategory::Init()
{
}

////////////////////////////////////////////////////////////////////////////////
/// initialize the circular tree

void TMVA::MethodCategory::InitCircularTree(const DataSetInfo& dsi)
{
   delete fCatTree;

   std::vector<VariableInfo>::const_iterator viIt;
   const std::vector<VariableInfo>& vars  = dsi.GetVariableInfos();
   const std::vector<VariableInfo>& specs = dsi.GetSpectatorInfos();

   Bool_t hasAllExternalLinks = kTRUE;
   for (viIt = vars.begin(); viIt != vars.end(); ++viIt)
      if( viIt->GetExternalLink() == 0 ) {
         hasAllExternalLinks = kFALSE;
         break;
      }
   for (viIt = specs.begin(); viIt != specs.end(); ++viIt)
      if( viIt->GetExternalLink() == 0 ) {
         hasAllExternalLinks = kFALSE;
         break;
      }

   if(!hasAllExternalLinks) return;

   {
      // Rather than having TTree::TTree add to the current directory and then remove it, let
      // make sure to not add it in the first place.
      // The add-then-remove can lead to  a problem if gDirectory points to the same directory (for example
      // gROOT) in the current thread and another one (and both try to add to the directory at the same time).
      TDirectory::TContext ctxt(nullptr);
      fCatTree = new TTree(Form("Circ%s",GetMethodName().Data()),"Circular Tree for categorization");
      fCatTree->SetCircular(1);
   }

   for (viIt = vars.begin(); viIt != vars.end(); ++viIt) {
      const VariableInfo& vi = *viIt;
      fCatTree->Branch(vi.GetExpression(),(Float_t*)vi.GetExternalLink(), TString(vi.GetExpression())+TString("/F"));
   }
   for (viIt = specs.begin(); viIt != specs.end(); ++viIt) {
      const VariableInfo& vi = *viIt;
      if(vi.GetVarType()=='C') continue;
      fCatTree->Branch(vi.GetExpression(),(Float_t*)vi.GetExternalLink(), TString(vi.GetExpression())+TString("/F"));
   }

   for(UInt_t cat=0; cat!=fCategoryCuts.size(); ++cat) {
      fCatFormulas.push_back(new TTreeFormula(Form("Category_%i",cat), fCategoryCuts[cat].GetTitle(), fCatTree));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// train all sub-classifiers

void TMVA::MethodCategory::Train()
{
   // specify the minimum # of training events and set 'classification'
   const Int_t  MinNoTrainingEvents = 10;

   Types::EAnalysisType analysisType = GetAnalysisType();

   // start the training
   Log() << kINFO << "Train all sub-classifiers for "
         << (analysisType == Types::kRegression ? "Regression" : "Classification") << " ..." << Endl;

   // don't do anything if no sub-classifier booked
   if (fMethods.empty()) {
      Log() << kINFO << "...nothing found to train" << Endl;
      return;
   }

   std::vector<IMethod*>::iterator itrMethod;

   // iterate over all booked sub-classifiers  and train them
   for (itrMethod = fMethods.begin(); itrMethod != fMethods.end(); ++itrMethod ) {

      MethodBase* mva = dynamic_cast<MethodBase*>(*itrMethod);
      if(!mva) continue;
      mva->SetAnalysisType( analysisType );
      if (!mva->HasAnalysisType( analysisType,
                                 mva->DataInfo().GetNClasses(),
                                 mva->DataInfo().GetNTargets() ) ) {
         Log() << kWARNING << "Method " << mva->GetMethodTypeName() << " is not capable of handling " ;
         if (analysisType == Types::kRegression)
            Log() << "regression with " << mva->DataInfo().GetNTargets() << " targets." << Endl;
         else
            Log() << "classification with " << mva->DataInfo().GetNClasses() << " classes." << Endl;
         itrMethod = fMethods.erase( itrMethod );
         continue;
      }
      if (mva->Data()->GetNTrainingEvents() >= MinNoTrainingEvents) {

         Log() << kINFO << "Train method: " << mva->GetMethodName() << " for "
               << (analysisType == Types::kRegression ? "Regression" : "Classification") << Endl;
         mva->TrainMethod();
         Log() << kINFO << "Training finished" << Endl;

      } else {

         Log() << kWARNING << "Method " << mva->GetMethodName()
               << " not trained (training tree has less entries ["
               << mva->Data()->GetNTrainingEvents()
               << "] than required [" << MinNoTrainingEvents << "]" << Endl;

         Log() << kERROR << " w/o training/test events for that category, I better stop here and let you fix " << Endl;
         Log() << kFATAL << "that one first, otherwise things get too messy later ... " << Endl;

      }
   }

   if (analysisType != Types::kRegression) {

      // variable ranking
      Log() << kINFO << "Begin ranking of input variables..." << Endl;
      for (itrMethod = fMethods.begin(); itrMethod != fMethods.end(); ++itrMethod) {
         MethodBase* mva = dynamic_cast<MethodBase*>(*itrMethod);
         if (mva && mva->Data()->GetNTrainingEvents() >= MinNoTrainingEvents) {
            const Ranking* ranking = (*itrMethod)->CreateRanking();
            if (ranking != 0)
               ranking->Print();
            else
               Log() << kINFO << "No variable ranking supplied by classifier: "
                     << dynamic_cast<MethodBase*>(*itrMethod)->GetMethodName() << Endl;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create XML description of Category classifier

void TMVA::MethodCategory::AddWeightsXMLTo( void* parent ) const
{
   void* wght = gTools().AddChild(parent, "Weights");
   gTools().AddAttr( wght, "NSubMethods", fMethods.size() );
   void* submethod(0);

   // iterate over methods and write them to XML file
   for (UInt_t i=0; i<fMethods.size(); i++) {
      MethodBase* method = dynamic_cast<MethodBase*>(fMethods[i]);
      submethod = gTools().AddChild(wght, "SubMethod");
      gTools().AddAttr(submethod, "Index", i);
      gTools().AddAttr(submethod, "Method", method->GetMethodTypeName() + "::" + method->GetMethodName());
      gTools().AddAttr(submethod, "Cut", fCategoryCuts[i]);
      gTools().AddAttr(submethod, "Variables", fVars[i]);
      method->WriteStateToXML( submethod );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// read weights of sub-classifiers of MethodCategory from xml weight file

void TMVA::MethodCategory::ReadWeightsFromXML( void* wghtnode )
{
   UInt_t nSubMethods;
   TString fullMethodName;
   TString methodType;
   TString methodTitle;
   TString theCutString;
   TString theVariables;
   Int_t titleLength;
   gTools().ReadAttr( wghtnode, "NSubMethods",  nSubMethods );
   void* subMethodNode = gTools().GetChild(wghtnode);

   Log() << kINFO << "Recreating sub-classifiers from XML-file " << Endl;

   // recreate all sub-methods from weight file
   for (UInt_t i=0; i<nSubMethods; i++) {
      gTools().ReadAttr( subMethodNode, "Method",    fullMethodName );
      gTools().ReadAttr( subMethodNode, "Cut",       theCutString   );
      gTools().ReadAttr( subMethodNode, "Variables", theVariables   );

      // determine sub-method type
      methodType = fullMethodName(0,fullMethodName.Index("::"));
      if (methodType.Contains(" ")) methodType = methodType(methodType.Last(' ')+1,methodType.Length());

      // determine sub-method title
      titleLength = fullMethodName.Length()-fullMethodName.Index("::")-2;
      methodTitle = fullMethodName(fullMethodName.Index("::")+2,titleLength);

      // reconstruct dsi for sub-method
      DataSetInfo& dsi = CreateCategoryDSI(TCut(theCutString), theVariables, methodTitle);

      // recreate sub-method from weights and add to fMethods
      MethodBase* method = dynamic_cast<MethodBase*>( ClassifierFactory::Instance().Create( methodType.Data(),
                                                                                            dsi, "none" ) );
      if(method==0)
         Log() << kFATAL << "Could not create sub-method " << method << " from XML." << Endl;

      method->SetupMethod();
      method->ReadStateFromXML(subMethodNode);

      fMethods.push_back(method);
      fCategoryCuts.push_back(TCut(theCutString));
      fVars.push_back(theVariables);

      DataSetInfo& primaryDSI = DataInfo();

      UInt_t spectatorIdx = 10000;
      UInt_t counter=0;

      // find the spectator index
      std::vector<VariableInfo>& spectators=primaryDSI.GetSpectatorInfos();
      std::vector<VariableInfo>::iterator itrVarInfo;
      TString specName= Form("%s_cat%i", GetName(),(int)fCategorySpecIdx.size()+1);

      for (itrVarInfo = spectators.begin(); itrVarInfo != spectators.end(); ++itrVarInfo, ++counter) {
         if((specName==itrVarInfo->GetLabel()) || (specName==itrVarInfo->GetExpression())) {
            spectatorIdx=counter;
            fCategorySpecIdx.push_back(spectatorIdx);
            break;
         }
      }

      subMethodNode = gTools().GetNextChild(subMethodNode);
   }

   InitCircularTree(DataInfo());

}

////////////////////////////////////////////////////////////////////////////////
/// process user options

void TMVA::MethodCategory::ProcessOptions()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodCategory::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "This method allows to define different categories of events. The" <<Endl;
   Log() << "categories are defined via cuts on the variables. For each" << Endl;
   Log() << "category, a different classifier and set of variables can be" <<Endl;
   Log() << "specified. The categories which are defined for this method must" << Endl;
   Log() << "be disjoint." << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// no ranking

const TMVA::Ranking* TMVA::MethodCategory::CreateRanking()
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TMVA::MethodCategory::PassesCut( const Event* ev, UInt_t methodIdx )
{
   // if it's not a simple 'spectator' variable (0 or 1) that the categories are defined by
   // (but rather some 'formula' (i.e. eta>0), then this formulas are stored in fCatTree and that
   // one will be evaluated.. (the formulae return 'true' or 'false'
   if (fCatTree) {
      if (methodIdx>=fCatFormulas.size()) {
         Log() << kFATAL << "Large method index " << methodIdx << ", number of category formulas = "
               << fCatFormulas.size() << Endl;
      }
      TTreeFormula* f = fCatFormulas[methodIdx];
      return f->EvalInstance(0) > 0.5;
   }
   // otherwise, it simply looks if "variable == true"  ("greater 0.5 to be "sure" )
   else {

      // checks whether an event lies within a cut
      if (methodIdx>=fCategorySpecIdx.size()) {
         Log() << kFATAL << "Unknown method index " << methodIdx << " maximum allowed index="
               << fCategorySpecIdx.size() << Endl;
      }
      UInt_t spectatorIdx = fCategorySpecIdx[methodIdx];
      Float_t specVal = ev->GetSpectator(spectatorIdx);
      Bool_t pass = (specVal>0.5);
      return pass;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// returns the mva value of the right sub-classifier

Double_t TMVA::MethodCategory::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   if (fMethods.empty()) return 0;

   UInt_t methodToUse = 0;
   const Event* ev = GetEvent();

   // determine which sub-classifier to use for this event
   Int_t suitableCutsN = 0;

   for (UInt_t i=0; i<fMethods.size(); ++i) {
      if (PassesCut(ev, i)) {
         ++suitableCutsN;
         methodToUse=i;
      }
   }

   if (suitableCutsN == 0) {
      Log() << kWARNING << "Event does not lie within the cut of any sub-classifier." << Endl;
      return 0;
   }

   if (suitableCutsN > 1) {
      Log() << kFATAL << "The defined categories are not disjoint." << Endl;
      return 0;
   }

   // get mva value from the suitable sub-classifier
   ev->SetVariableArrangement(&fVarMaps[methodToUse]);
   Double_t mvaValue = dynamic_cast<MethodBase*>(fMethods[methodToUse])->GetMvaValue(ev,err,errUpper);
   ev->SetVariableArrangement(0);

   return mvaValue;
}

////////////////////////////////////////////////////////////////////////////////
/// returns the mva values of the multi-class right sub-classifier
///
const std::vector<Float_t> &TMVA::MethodCategory::GetMulticlassValues()
{
   if (fMethods.empty())
      return MethodBase::GetMulticlassValues();

   UInt_t methodToUse = 0;
   const Event *ev = GetEvent();

   // determine which sub-classifier to use for this event
   Int_t suitableCutsN = 0;

   for (UInt_t i = 0; i < fMethods.size(); ++i) {
      if (PassesCut(ev, i)) {
         ++suitableCutsN;
         methodToUse = i;
      }
   }

   if (suitableCutsN == 0) {
      Log() << kWARNING << "Event does not lie within the cut of any sub-classifier." << Endl;
      return MethodBase::GetMulticlassValues();
   }

   if (suitableCutsN > 1) {
      Log() << kFATAL << "The defined categories are not disjoint." << Endl;
      return MethodBase::GetMulticlassValues();
   }
   MethodBase *meth = dynamic_cast<MethodBase *>(fMethods[methodToUse]);
   if (!meth) {
      Log() << kFATAL << "method not found in Category Regression method" << Endl;
      return MethodBase::GetMulticlassValues();
   }
   // get mva value from the suitable sub-classifier
   ev->SetVariableArrangement(&fVarMaps[methodToUse]);
   auto & result =  meth->GetMulticlassValues();
   ev->SetVariableArrangement(nullptr);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// returns the mva value of the right sub-classifier

const std::vector<Float_t> &TMVA::MethodCategory::GetRegressionValues()
{
   if (fMethods.empty()) return MethodBase::GetRegressionValues();

   UInt_t methodToUse = 0;
   const Event* ev = GetEvent();

   // determine which sub-classifier to use for this event
   Int_t suitableCutsN = 0;

   for (UInt_t i=0; i<fMethods.size(); ++i) {
      if (PassesCut(ev, i)) {
         ++suitableCutsN;
         methodToUse=i;
      }
   }

   if (suitableCutsN == 0) {
      Log() << kWARNING << "Event does not lie within the cut of any sub-classifier." << Endl;
      return MethodBase::GetRegressionValues();
   }

   if (suitableCutsN > 1) {
      Log() << kFATAL << "The defined categories are not disjoint." << Endl;
      return MethodBase::GetRegressionValues();
   }
   MethodBase* meth = dynamic_cast<MethodBase*>(fMethods[methodToUse]);
   if (!meth){
      Log() << kFATAL << "method not found in Category Regression method" << Endl;
      return MethodBase::GetRegressionValues();
   }
   // get mva value from the suitable sub-classifier
   ev->SetVariableArrangement(&fVarMaps[methodToUse]);
   auto & result =  meth->GetRegressionValues(ev);
   return result;
}
