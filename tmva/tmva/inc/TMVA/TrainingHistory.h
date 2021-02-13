/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MsgLogger                                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Author:                                                                        *
 *      Joseph McKenna        <Joseph.McKenna@cern.ch> - Aarhus, Denmark          *
 *                                                                                *
* Copyright (c) 2019:                                                             *
 *      Aarhus, Denmark                                                           *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/


#ifndef ROOT_TMVA_TrainingHistory
#define ROOT_TMVA_TrainingHistory

#include <vector>
#include "TString.h"
#include <map>

namespace TMVA {

   class TrainingHistory {

   public:
      typedef std::vector<std::pair<Int_t,Double_t>> IterationRecord;
      TrainingHistory();
      virtual ~TrainingHistory();

      void AddValue(TString Property, Int_t stage, Double_t value);
      void SaveHistory(TString Name);
   private:
      std::map<TString,int> fHistoryMap;
      std::vector<IterationRecord*> fHistoryData;
   
   };

} // namespace TMVA

#endif
