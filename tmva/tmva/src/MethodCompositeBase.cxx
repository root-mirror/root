// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss,Or Cohen

/*****************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis  *
 * Package: TMVA                                                             *
 * Class  : MethodCompositeBase                                              *
 * Web    : http://tmva.sourceforge.net                                      *
 *                                                                           *
 * Description:                                                              *
 *      Virtual base class for all MVA method                                *
 *                                                                           *
 * Authors (alphabetical):                                                   *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland         *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - MSU, USA                  *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada    *
 *      Or Cohen        <orcohenor@gmail.com>    - Weizmann Inst., Israel    *
 *                                                                           *
 * Copyright (c) 2005:                                                       *
 *      CERN, Switzerland                                                    *
 *      U. of Victoria, Canada                                               *
 *      MPI-K Heidelberg, Germany                                            *
 *      LAPP, Annecy, France                                                 *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted according to the terms listed in LICENSE      *
 * (http://tmva.sourceforge.net/LICENSE)                                     *
 *****************************************************************************/

/*! \class TMVA::MethodCompositeBase
\ingroup TMVA

Virtual base class for combining several TMVA method.

This class is virtual class meant to combine more than one classifier
together. The training of the classifiers is done by classes that are
derived from this one, while the saving and loading of weights file
and the evaluation is done here.
*/

#include "TMVA/MethodCompositeBase.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Factory.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodBoost.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/Config.h"

#include "TRandom3.h"

#include <iostream>
#include <algorithm>
#include <vector>


using std::vector;

ClassImp(TMVA::MethodCompositeBase);

////////////////////////////////////////////////////////////////////////////////

TMVA::MethodCompositeBase::MethodCompositeBase( const TString& jobName,
                                                Types::EMVA methodType,
                                                const TString& methodTitle,
                                                DataSetInfo& theData,
                                                const TString& theOption )
: TMVA::MethodBase( jobName, methodType, methodTitle, theData, theOption),
   fCurrentMethodIdx(0), fCurrentMethod(0)
{}

////////////////////////////////////////////////////////////////////////////////

TMVA::MethodCompositeBase::MethodCompositeBase( Types::EMVA methodType,
                                                DataSetInfo& dsi,
                                                const TString& weightFile)
   : TMVA::MethodBase( methodType, dsi, weightFile),
     fCurrentMethodIdx(0), fCurrentMethod(0)
{}

////////////////////////////////////////////////////////////////////////////////
/// returns pointer to MVA that corresponds to given method title

TMVA::IMethod* TMVA::MethodCompositeBase::GetMethod( const TString &methodTitle ) const
{
   std::vector<IMethod*>::const_iterator itrMethod    = fMethods.begin();
   std::vector<IMethod*>::const_iterator itrMethodEnd = fMethods.end();

   for (; itrMethod != itrMethodEnd; ++itrMethod) {
      MethodBase* mva = dynamic_cast<MethodBase*>(*itrMethod);
      if ( (mva->GetMethodName())==methodTitle ) return mva;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// returns pointer to MVA that corresponds to given method index

TMVA::IMethod* TMVA::MethodCompositeBase::GetMethod( const Int_t index ) const
{
   std::vector<IMethod*>::const_iterator itrMethod = fMethods.begin()+index;
   if (itrMethod<fMethods.end()) return *itrMethod;
   else                          return 0;
}


////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodCompositeBase::AddWeightsXMLTo( void* parent ) const
{
   void* wght = gTools().AddChild(parent, "Weights");
   gTools().AddAttr( wght, "NMethods",   fMethods.size()   );
   for (UInt_t i=0; i< fMethods.size(); i++)
      {
         void* methxml = gTools().AddChild( wght, "Method" );
         MethodBase* method = dynamic_cast<MethodBase*>(fMethods[i]);
         gTools().AddAttr(methxml,"Index",          i );
         gTools().AddAttr(methxml,"Weight",         fMethodWeight[i]);
         gTools().AddAttr(methxml,"MethodSigCut",   method->GetSignalReferenceCut());
         gTools().AddAttr(methxml,"MethodSigCutOrientation", method->GetSignalReferenceCutOrientation());
         gTools().AddAttr(methxml,"MethodTypeName", method->GetMethodTypeName());
         gTools().AddAttr(methxml,"MethodName",     method->GetMethodName()   );
         gTools().AddAttr(methxml,"JobName",        method->GetJobName());
         gTools().AddAttr(methxml,"Options",        method->GetOptions());
         if (method->fTransformationPointer)
            gTools().AddAttr(methxml,"UseMainMethodTransformation", TString("true"));
         else
            gTools().AddAttr(methxml,"UseMainMethodTransformation", TString("false"));
         method->AddWeightsXMLTo(methxml);
      }
}

////////////////////////////////////////////////////////////////////////////////
/// delete methods

TMVA::MethodCompositeBase::~MethodCompositeBase( void )
{
   std::vector<IMethod*>::iterator itrMethod = fMethods.begin();
   for (; itrMethod != fMethods.end(); ++itrMethod) {
      Log() << kVERBOSE << "Delete method: " << (*itrMethod)->GetName() << Endl;
      delete (*itrMethod);
   }
   fMethods.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// XML streamer

void TMVA::MethodCompositeBase::ReadWeightsFromXML( void* wghtnode )
{
   UInt_t nMethods;
   TString methodName, methodTypeName, jobName, optionString;

   for (UInt_t i=0;i<fMethods.size();i++) delete fMethods[i];
   fMethods.clear();
   fMethodWeight.clear();
   gTools().ReadAttr( wghtnode, "NMethods",  nMethods );
   void* ch = gTools().GetChild(wghtnode);
   for (UInt_t i=0; i< nMethods; i++) {
      Double_t methodWeight, methodSigCut, methodSigCutOrientation;
      gTools().ReadAttr( ch, "Weight",   methodWeight   );
      gTools().ReadAttr( ch, "MethodSigCut", methodSigCut);
      gTools().ReadAttr( ch, "MethodSigCutOrientation", methodSigCutOrientation);
      gTools().ReadAttr( ch, "MethodTypeName",  methodTypeName );
      gTools().ReadAttr( ch, "MethodName",  methodName );
      gTools().ReadAttr( ch, "JobName",  jobName );
      gTools().ReadAttr( ch, "Options",  optionString );

      // Bool_t rerouteTransformation = kFALSE;
      if (gTools().HasAttr( ch, "UseMainMethodTransformation")) {
         TString rerouteString("");
         gTools().ReadAttr( ch, "UseMainMethodTransformation", rerouteString );
         rerouteString.ToLower();
         // if (rerouteString=="true")
         //    rerouteTransformation=kTRUE;
      }

      //remove trailing "~" to signal that options have to be reused
      optionString.ReplaceAll("~","");
      //ignore meta-options for method Boost
      optionString.ReplaceAll("Boost_","~Boost_");
      optionString.ReplaceAll("!~","~!");

      if (i==0){
         // the cast on MethodBoost is ugly, but a similar line is also in ReadWeightsFromFile --> needs to be fixed later
         ((TMVA::MethodBoost*)this)->BookMethod( Types::Instance().GetMethodType( methodTypeName), methodName,  optionString );
      }
      fMethods.push_back(
         ClassifierFactory::Instance().Create(methodTypeName.Data(), jobName, methodName, DataInfo(), optionString));

      fMethodWeight.push_back(methodWeight);
      MethodBase* meth = dynamic_cast<MethodBase*>(fMethods.back());

      if(meth==0)
         Log() << kFATAL << "Could not read method from XML" << Endl;

      void* methXML = gTools().GetChild(ch);

      TString _fFileDir= meth->DataInfo().GetName();
      _fFileDir+="/"+gConfig().GetIONames().fWeightFileDir;
      meth->SetWeightFileDir(_fFileDir);
      meth->SetModelPersistence(IsModelPersistence());
      meth->SetSilentFile(IsSilentFile());
      meth->SetupMethod();
      meth->SetMsgType(kWARNING);
      meth->ParseOptions();
      meth->ProcessSetup();
      meth->CheckSetup();
      meth->ReadWeightsFromXML(methXML);
      meth->SetSignalReferenceCut(methodSigCut);
      meth->SetSignalReferenceCutOrientation(methodSigCutOrientation);

      meth->RerouteTransformationHandler (&(this->GetTransformationHandler()));

      ch = gTools().GetNextChild(ch);
   }
   //Log() << kINFO << "Reading methods from XML done " << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// text streamer

void  TMVA::MethodCompositeBase::ReadWeightsFromStream( std::istream& istr )
{
   TString var, dummy;
   TString methodName, methodTitle=GetMethodName(),
      jobName=GetJobName(),optionString=GetOptions();
   UInt_t methodNum; Double_t methodWeight;
   // and read the Weights (BDT coefficients)
   // coverity[tainted_data_argument]
   istr >> dummy >> methodNum;
   Log() << kINFO << "Read " << methodNum << " Classifiers" << Endl;
   for (UInt_t i=0;i<fMethods.size();i++) delete fMethods[i];
   fMethods.clear();
   fMethodWeight.clear();
   for (UInt_t i=0; i<methodNum; i++) {
      istr >> dummy >> methodName >>  dummy >> fCurrentMethodIdx >> dummy >> methodWeight;
      if ((UInt_t)fCurrentMethodIdx != i) {
         Log() << kFATAL << "Error while reading weight file; mismatch MethodIndex="
               << fCurrentMethodIdx << " i=" << i
               << " MethodName " << methodName
               << " dummy " << dummy
               << " MethodWeight= " << methodWeight
               << Endl;
      }
      if (GetMethodType() != Types::kBoost || i==0) {
         istr >> dummy >> jobName;
         istr >> dummy >> methodTitle;
         istr >> dummy >> optionString;
         if (GetMethodType() == Types::kBoost)
            ((TMVA::MethodBoost*)this)->BookMethod( Types::Instance().GetMethodType( methodName), methodTitle,  optionString );
      }
      else methodTitle=Form("%s (%04i)",GetMethodName().Data(),fCurrentMethodIdx);
      fMethods.push_back(
         ClassifierFactory::Instance().Create(methodName.Data(), jobName, methodTitle, DataInfo(), optionString));
      fMethodWeight.push_back( methodWeight );
      if(MethodBase* m = dynamic_cast<MethodBase*>(fMethods.back()) )
         m->ReadWeightsFromStream(istr);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// return composite MVA response

Double_t TMVA::MethodCompositeBase::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   Double_t mvaValue = 0;
   for (UInt_t i=0;i< fMethods.size(); i++) mvaValue+=fMethods[i]->GetMvaValue()*fMethodWeight[i];

   // cannot determine error
   NoErrorCalc(err, errUpper);

   return mvaValue;
}
