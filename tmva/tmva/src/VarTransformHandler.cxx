/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VarTransformHandler                                                   *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *       Implementation of unsupervised variable transformation methods           *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *       Abhinav Moudgil <abhinav.moudgil@research.iiit.ac.in> - IIIT-H, India    *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *       CERN, Switzerland                                                        *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#include "TMVA/VarTransformHandler.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Event.h"
#include "TMVA/DataInputHandler.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodDNN.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/VariableInfo.h"

#include "TMath.h"
#include "TVectorD.h"
#include "TMatrix.h"
#include "TMatrixTSparse.h"
#include "TMatrixDSparsefwd.h"

#include <algorithm>
#include <iomanip>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::VarTransformHandler::VarTransformHandler( DataLoader* dl )
   : fLogger     ( new MsgLogger(TString("VarTransformHandler").Data(), kINFO) ),
     fDataSetInfo(dl->GetDataSetInfo()),
     fDataLoader (dl),
     fEvents (fDataSetInfo.GetDataSet()->GetEventCollection())
{
   Log() << kINFO << "Number of events - " << fEvents.size() << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::VarTransformHandler::~VarTransformHandler()
{
    // do something
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// Computes variance of all the variables and
/// returns a new DataLoader with the selected variables whose variance is above a specific threshold.
/// Threshold can be provided by user otherwise default value is 0 i.e. remove the variables which have same value in all
/// the events.
///
/// \param[in] threshold value (Double)
///
/// Transformation Definition String Format: "VT(optional float value)"
///
/// Usage examples:
///
/// String    | Description
/// -------   |----------------------------------------
/// "VT"      | Select variables whose variance is above threshold value = 0 (Default)
/// "VT(1.5)" | Select variables whose variance is above threshold value = 1.5

TMVA::DataLoader* TMVA::VarTransformHandler::VarianceThreshold(Double_t threshold)
{
   CalcNorm();
   const UInt_t nvars = fDataSetInfo.GetNVariables();
   Log() << kINFO << "Number of variables before transformation: " << nvars << Endl;
   std::vector<VariableInfo>& vars = fDataSetInfo.GetVariableInfos();

   // return a new dataloader
   // iterate over all variables, ignore the ones whose variance is below specific threshold
   // DataLoader *transformedLoader=(DataLoader *)fDataLoader->Clone("vt_transformed_dataset");
   // TMVA::DataLoader *transformedLoader = new TMVA::DataLoader(fDataSetInfo.GetName());
   TMVA::DataLoader *transformedLoader = new TMVA::DataLoader("vt_transformed_dataset");
   Log() << kINFO << "Selecting variables whose variance is above threshold value = " << threshold << Endl;
   Int_t maxL = fDataSetInfo.GetVariableNameMaxLength();
   maxL = maxL + 16;
   Log() << kINFO << "----------------------------------------------------------------" << Endl;
   Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(maxL) << "Selected Variables";
   Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(10) << "Variance" << Endl;
   Log() << kINFO << "----------------------------------------------------------------" << Endl;
   for (UInt_t ivar=0; ivar<nvars; ivar++) {
      Double_t variance =  vars[ivar].GetVariance();
      if (variance > threshold)
      {
         Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(maxL) << vars[ivar].GetExpression();
         Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(maxL) << variance << Endl;
         transformedLoader->AddVariable(vars[ivar].GetExpression(), vars[ivar].GetVarType());
      }
   }
   CopyDataLoader(transformedLoader,fDataLoader);
   Log() << kINFO << "----------------------------------------------------------------" << Endl;
   // CopyDataLoader(transformedLoader, fDataLoader);
   // DataLoader *transformedLoader=(DataLoader *)fDataLoader->Clone(fDataSetInfo.GetName());
   transformedLoader->PrepareTrainingAndTestTree(fDataLoader->GetDataSetInfo().GetCut("Signal"), fDataLoader->GetDataSetInfo().GetCut("Background"), fDataLoader->GetDataSetInfo().GetSplitOptions());
   Log() << kINFO << "Number of variables after transformation: " << transformedLoader->GetDataSetInfo().GetNVariables() << Endl;

   return transformedLoader;
}

///////////////////////////////////////////////////////////////////////////////
////////////////////////////// Utility methods ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Updates maximum and minimum value of a variable or target

void TMVA::VarTransformHandler::UpdateNorm (Int_t ivar, Double_t x)
{
   Int_t nvars = fDataSetInfo.GetNVariables();
   std::vector<VariableInfo>& vars = fDataSetInfo.GetVariableInfos();
   std::vector<VariableInfo>& tars = fDataSetInfo.GetTargetInfos();
   if( ivar < nvars ){
      if (x < vars[ivar].GetMin()) vars[ivar].SetMin(x);
      if (x > vars[ivar].GetMax()) vars[ivar].SetMax(x);
   }
   else{
      if (x < tars[ivar-nvars].GetMin()) tars[ivar-nvars].SetMin(x);
      if (x > tars[ivar-nvars].GetMax()) tars[ivar-nvars].SetMax(x);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Computes maximum, minimum, mean, RMS and variance for all
/// variables and targets

void TMVA::VarTransformHandler::CalcNorm()
{
   const std::vector<TMVA::Event*>& events = fDataSetInfo.GetDataSet()->GetEventCollection();

   const UInt_t nvars = fDataSetInfo.GetNVariables();
   const UInt_t ntgts = fDataSetInfo.GetNTargets();
   std::vector<VariableInfo>& vars = fDataSetInfo.GetVariableInfos();
   std::vector<VariableInfo>& tars = fDataSetInfo.GetTargetInfos();

   UInt_t nevts = events.size();

   TVectorD x2( nvars+ntgts ); x2 *= 0;
   TVectorD x0( nvars+ntgts ); x0 *= 0;
   TVectorD v0( nvars+ntgts ); v0 *= 0;

   Double_t sumOfWeights = 0;
   for (UInt_t ievt=0; ievt<nevts; ievt++) {
      const Event* ev = events[ievt];

      Double_t weight = ev->GetWeight();
      sumOfWeights += weight;
      for (UInt_t ivar=0; ivar<nvars; ivar++) {
         Double_t x = ev->GetValue(ivar);
         if (ievt==0) {
            vars[ivar].SetMin(x);
            vars[ivar].SetMax(x);
         }
         else {
            UpdateNorm(ivar,  x );
         }
         x0(ivar) += x*weight;
         x2(ivar) += x*x*weight;
      }
      for (UInt_t itgt=0; itgt<ntgts; itgt++) {
         Double_t x = ev->GetTarget(itgt);
         if (ievt==0) {
            tars[itgt].SetMin(x);
            tars[itgt].SetMax(x);
         }
         else {
            UpdateNorm( nvars+itgt,  x );
         }
         x0(nvars+itgt) += x*weight;
         x2(nvars+itgt) += x*x*weight;
      }
   }

   if (sumOfWeights <= 0) {
      Log() << kFATAL << " the sum of event weights calculated for your input is == 0"
            << " or exactly: " << sumOfWeights << " there is obviously some problem..."<< Endl;
   }

   // set Mean and RMS
   for (UInt_t ivar=0; ivar<nvars; ivar++) {
      Double_t mean = x0(ivar)/sumOfWeights;

      vars[ivar].SetMean( mean );
      if (x2(ivar)/sumOfWeights - mean*mean < 0) {
         Log() << kFATAL << " the RMS of your input variable " << ivar
               << " evaluates to an imaginary number: sqrt("<< x2(ivar)/sumOfWeights - mean*mean
               <<") .. sometimes related to a problem with outliers and negative event weights"
               << Endl;
      }
      vars[ivar].SetRMS( TMath::Sqrt( x2(ivar)/sumOfWeights - mean*mean) );
   }
   for (UInt_t itgt=0; itgt<ntgts; itgt++) {
      Double_t mean = x0(nvars+itgt)/sumOfWeights;
      tars[itgt].SetMean( mean );
      if (x2(nvars+itgt)/sumOfWeights - mean*mean < 0) {
         Log() << kFATAL << " the RMS of your target variable " << itgt
               << " evaluates to an imaginary number: sqrt(" << x2(nvars+itgt)/sumOfWeights - mean*mean
               <<") .. sometimes related to a problem with outliers and negative event weights"
               << Endl;
      }
      tars[itgt].SetRMS( TMath::Sqrt( x2(nvars+itgt)/sumOfWeights - mean*mean) );
   }

   // calculate variance
   for (UInt_t ievt=0; ievt<nevts; ievt++) {
      const Event* ev = events[ievt];
      Double_t weight = ev->GetWeight();

      for (UInt_t ivar=0; ivar<nvars; ivar++) {
         Double_t x = ev->GetValue(ivar);
         Double_t mean = vars[ivar].GetMean();
         v0(ivar) += weight*(x-mean)*(x-mean);
      }

      for (UInt_t itgt=0; itgt<ntgts; itgt++) {
         Double_t x = ev->GetTarget(itgt);
         Double_t mean = tars[itgt].GetMean();
         v0(nvars+itgt) += weight*(x-mean)*(x-mean);
      }
   }

   Int_t maxL = fDataSetInfo.GetVariableNameMaxLength();
   maxL = maxL + 8;
   Log() << kINFO << "----------------------------------------------------------------" << Endl;
   Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(maxL) << "Variables";
   Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(10) << "Variance" << Endl;
   Log() << kINFO << "----------------------------------------------------------------" << Endl;

   // set variance
   Log() << std::setprecision(5);
   for (UInt_t ivar=0; ivar<nvars; ivar++) {
      Double_t variance = v0(ivar)/sumOfWeights;
      vars[ivar].SetVariance( variance );
      Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(maxL) << vars[ivar].GetExpression();
      Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(maxL) << variance << Endl;
   }

   maxL = fDataSetInfo.GetTargetNameMaxLength();
   maxL = maxL + 8;
   Log() << kINFO << "----------------------------------------------------------------" << Endl;
   Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(maxL) << "Targets";
   Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(10) << "Variance" << Endl;
   Log() << kINFO << "----------------------------------------------------------------" << Endl;

   for (UInt_t itgt=0; itgt<ntgts; itgt++) {
      Double_t variance = v0(nvars+itgt)/sumOfWeights;
      tars[itgt].SetVariance( variance );
      Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(maxL) << tars[itgt].GetExpression();
      Log() << kINFO << std::setiosflags(std::ios::left) << std::setw(maxL) << variance << Endl;
   }

   Log() << kINFO << "Set minNorm/maxNorm for variables to: " << Endl;
   Log() << std::setprecision(3);
   for (UInt_t ivar=0; ivar<nvars; ivar++)
      Log() << "    " << vars[ivar].GetExpression()
            << "\t: [" << vars[ivar].GetMin() << "\t, " << vars[ivar].GetMax() << "\t] " << Endl;
   Log() << kINFO << "Set minNorm/maxNorm for targets to: " << Endl;
   Log() << std::setprecision(3);
   for (UInt_t itgt=0; itgt<ntgts; itgt++)
      Log() << "    " << tars[itgt].GetExpression()
            << "\t: [" << tars[itgt].GetMin() << "\t, " << tars[itgt].GetMax() << "\t] " << Endl;
   Log() << std::setprecision(5); // reset to better value
}

////////////////////////////////////////////////////////////////////////////////
void TMVA::VarTransformHandler::CopyDataLoader(TMVA::DataLoader* des, TMVA::DataLoader* src)
{
   for( std::vector<TreeInfo>::const_iterator treeinfo=src->DataInput().Sbegin();treeinfo!=src->DataInput().Send();++treeinfo)
   {
      des->AddSignalTree( (*treeinfo).GetTree(), (*treeinfo).GetWeight(),(*treeinfo).GetTreeType());
   }

   for( std::vector<TreeInfo>::const_iterator treeinfo=src->DataInput().Bbegin();treeinfo!=src->DataInput().Bend();++treeinfo)
   {
      des->AddBackgroundTree( (*treeinfo).GetTree(), (*treeinfo).GetWeight(),(*treeinfo).GetTreeType());
   }
}
