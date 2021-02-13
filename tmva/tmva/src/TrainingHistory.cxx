/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TrainingHistory                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Author:                                                                        *
 *      Joseph McKenna        <Joseph.McKenna@cern.ch> - Aarhus, Denmark          *
 *                                                                                *
 * Copyright (c) 2019:                                                            *
 *      Aarhus, Denmark                                                           *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::TrainingHistory
\ingroup TMVA

Tracking data from training. Eg, From deep learning record loss for each Epoch

*/


#include "TMVA/TrainingHistory.h"

#include "TH1D.h"

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::TrainingHistory::TrainingHistory()
{
}

TMVA::TrainingHistory::~TrainingHistory()
{
   for (auto p : fHistoryData) {
     delete p;
   } 
}

void TMVA::TrainingHistory::AddValue(TString Property,Int_t stage, Double_t value)
{
   if (!fHistoryMap.count(Property))
   {
      fHistoryMap[Property]=fHistoryData.size();
      IterationRecord* data=new IterationRecord();
      fHistoryData.push_back(data);
   }
   int iHistory=fHistoryMap.at(Property);
   //std::cout<<"HISTORY!"<<"Adding value ("<<Property<<"):"<<stage<<"\t"<<value<<std::endl;
   fHistoryData.at(iHistory)->push_back({stage,value});
}

void TMVA::TrainingHistory::SaveHistory(TString Name)
{
   //if (fHistoryData.empty()) return;
   for ( const auto &element : fHistoryMap ) {
      TString property = element.first;
      Int_t iHistory = element.second;
      Int_t nBins=fHistoryData.at(iHistory)->size();
      Double_t xMin=fHistoryData.at(iHistory)->front().first;
      Double_t xMax=fHistoryData.at(iHistory)->back().first;
      Double_t BinSize=(xMax-xMin)/(Double_t)(nBins-1);
      TH1D* h=new TH1D("TrainingHistory_"+Name+"_"+property,"TrainingHistory_"+Name+"_"+property,nBins,xMin-0.5*BinSize,xMax+0.5*BinSize);
      for (int i=0; i<nBins; i++) {
         h->AddBinContent(i+1,fHistoryData.at(iHistory)->at(i).second);
      }
      h->Print();
      if (h!=0) {
         h->Write();
         delete h;
      }

   }
}
