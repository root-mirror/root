#ifndef ROOT_TMVA_QUICKMVAPROBESTIMATOR 
#define ROOT_TMVA_QUICKMVAPROBESTIMATOR 

#include <vector>
#include <algorithm>

#include "TMVA/MsgLogger.h"

namespace TMVA {

   class QuickMVAProbEstimator {
   public:
    
      struct EventInfo{
         Double_t eventValue;
         Double_t eventWeight; 
         Int_t    eventType;  //signal or background
      };
      static bool compare(EventInfo e1, EventInfo e2){return e1.eventValue < e2.eventValue;}
    
   QuickMVAProbEstimator(Int_t nMin=40, Int_t nMax=5000):fIsSorted(false),fNMin(nMin),fNMax(nMax){ fLogger = new MsgLogger("QuickMVAProbEstimator");}
    
    
      virtual ~QuickMVAProbEstimator(){delete fLogger;}
      void AddEvent(Double_t val, Double_t weight, Int_t type);
    
   
      Double_t GetMVAProbAt(Double_t value);
    
    
   private:
      std::vector<EventInfo> fEvtVector;
      Bool_t                 fIsSorted;
      UInt_t                 fNMin;
      UInt_t                 fNMax;
    
      mutable MsgLogger*    fLogger;
      MsgLogger& Log() const { return *fLogger; }
    
      ClassDef(QuickMVAProbEstimator,0); // Interface to different separation critiera used in training algorithms
      
      
   };
}


#endif
