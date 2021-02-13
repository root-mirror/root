// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss, Omar Zapata

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VariableTransformBase                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch>  - CERN, Switzerland          *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#include "TMVA/VariableTransformBase.h"
#include "TMVA/VariableIdentityTransform.h"
#include "TMVA/VariableDecorrTransform.h"
#include "TMVA/VariableInfo.h"
#include "TMVA/VariablePCATransform.h"
#include "TMVA/VariableGaussTransform.h"
#include "TMVA/VariableNormalizeTransform.h"

#include "TMVA/Config.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Ranking.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/Version.h"
#include "TMVA/TransformationHandler.h"

#include "THashTable.h"
#include "TList.h"
#include "TObjString.h"

#include <algorithm>
#include <cassert>
#include <exception>
#include <stdexcept>
#include <set>

////////////////////////////////////////////////////////////////////////////////
/// create variable transformations

namespace TMVA {
void CreateVariableTransforms(const TString& trafoDefinitionIn,
                              TMVA::DataSetInfo& dataInfo,
                              TMVA::TransformationHandler& transformationHandler,
                              TMVA::MsgLogger& log)
{
   TString trafoDefinition(trafoDefinitionIn);
   if (trafoDefinition == "None") return; // no transformations

   // workaround for transformations to complicated to be handled by makeclass
   // count number of transformations with incomplete set of variables
   TString trafoDefinitionCheck(trafoDefinitionIn);
   int npartial = 0, ntrafo = 0;
   for (Int_t pos = 0, siz = trafoDefinition.Sizeof(); pos < siz; ++pos) {
      TString ch = trafoDefinition(pos,1);
      if ( ch == "(" ) npartial++;
      if ( ch == "+" || ch == ",") ntrafo++;
   }
   if (npartial>1) {
      log << kWARNING
          << "The use of multiple partial variable transformations during the "
             "application phase can be properly invoked via the \"Reader\", but "
             "it is not yet implemented in \"MakeClass\", the creation mechanism "
             "for standalone C++ application classes. The standalone C++ class "
             "produced by this training job is thus INCOMPLETE AND MUST NOT BE USED! "
             "The transformation in question is: " << trafoDefinitionIn << Endl;
      // ToDo make info and do not write the standalone class
      //
      // this does not work since this function is static
      // fDisableWriting=true; // disable creation of stand-alone class
      // ToDo we need to tell the transformation that it cannot write itself
   }
   // workaround end

   Int_t parenthesisCount = 0;
   for (Int_t position = 0, size = trafoDefinition.Sizeof(); position < size; ++position) {
      TString ch = trafoDefinition(position,1);
      if      (ch == "(") ++parenthesisCount;
      else if (ch == ")") --parenthesisCount;
      else if (ch == "," && parenthesisCount == 0) trafoDefinition.Replace(position,1,'+');
   }

   TList* trList = gTools().ParseFormatLine( trafoDefinition, "+" );
   TListIter trIt(trList);
   while (TObjString* os = (TObjString*)trIt()) {
      TString tdef = os->GetString();
      Int_t idxCls = -1;

      TString variables = "";
      if (tdef.Contains("(")) { // contains selection of variables
         Ssiz_t parStart = tdef.Index( "(" );
         Ssiz_t parLen   = tdef.Index( ")", parStart )-parStart+1;

         variables = tdef(parStart,parLen);
         tdef.Remove(parStart,parLen);
         variables.Remove(parLen-1,1);
         variables.Remove(0,1);
      }

      TList* trClsList = gTools().ParseFormatLine( tdef, "_" ); // split entry to get trf-name and class-name
      TListIter trClsIt(trClsList);
      if (trClsList->GetSize() < 1)
          log << kFATAL <<Form("Dataset[%s] : ",dataInfo.GetName())<< "Incorrect transformation string provided." << Endl;
      const TString& trName = ((TObjString*)trClsList->At(0))->GetString();

      if (trClsList->GetEntries() > 1) {
         TString trCls = "AllClasses";
         ClassInfo *ci = NULL;
         trCls  = ((TObjString*)trClsList->At(1))->GetString();
         if (trCls != "AllClasses") {
            ci = dataInfo.GetClassInfo( trCls );
            if (ci == NULL)
               log << kFATAL <<Form("Dataset[%s] : ",dataInfo.GetName())<< "Class " << trCls << " not known for variable transformation "
                   << trName << ", please check." << Endl;
            else
               idxCls = ci->GetNumber();
         }
      }

      VariableTransformBase* transformation = NULL;
      if (trName == "I" || trName == "Ident" || trName == "Identity") {
         if (variables.Length() == 0) variables = "_V_";
         transformation = new VariableIdentityTransform(dataInfo);
      }
      else if (trName == "D" || trName == "Deco" || trName == "Decorrelate") {
         if (variables.Length() == 0) variables = "_V_";
         transformation = new VariableDecorrTransform(dataInfo);
      }
      else if (trName == "P" || trName == "PCA") {
         if (variables.Length() == 0) variables = "_V_";
         transformation = new VariablePCATransform(dataInfo);
      }
      else if (trName == "U" || trName == "Uniform") {
         if (variables.Length() == 0) variables = "_V_,_T_";
         transformation = new VariableGaussTransform(dataInfo, "Uniform" );
      }
      else if (trName == "G" || trName == "Gauss") {
         if (variables.Length() == 0) variables = "_V_";
         transformation = new VariableGaussTransform(dataInfo);
      }
      else if (trName == "N" || trName == "Norm" || trName == "Normalise" || trName == "Normalize") {
         if (variables.Length() == 0) variables = "_V_,_T_";
         transformation = new VariableNormalizeTransform(dataInfo);
      }
      else if (trName == "S" || trName == "Scale" || trName == "ScaleNorm" ) {
         if (variables.Length() == 0) variables = "_V_,_T_";
         transformation = new VariableNormalizeTransform(dataInfo,"Scale");
      }
      else
         log << kFATAL << Form("Dataset[%s] : ",dataInfo.GetName())
             << "<ProcessOptions> Variable transform '"
             << trName << "' unknown." << Endl;


      if (transformation) {
         ClassInfo* clsInfo = dataInfo.GetClassInfo(idxCls);
         if (clsInfo)
            log << kHEADER << Form("[%s] : ",dataInfo.GetName())
                << "Create Transformation \"" << trName << "\" with reference class "
                << clsInfo->GetName() << "=("<< idxCls <<")" << Endl << Endl;
         else
            log << kHEADER << Form("[%s] : ",dataInfo.GetName())
                << "Create Transformation \"" << trName << "\" with events from all classes."
                << Endl << Endl;

         transformation->SelectInput(variables);
         transformationHandler.AddTransformation(transformation, idxCls);
      }
   }
}

}
