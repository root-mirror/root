/*! \class TMVA::QuickMVAProbEstimator
\ingroup TMVA

*/

#include "TMVA/QuickMVAProbEstimator.h"

#include "TMVA/MsgLogger.h"
#include "TMVA/Types.h"

#include "TMath.h"


void TMVA::QuickMVAProbEstimator::AddEvent(Double_t val, Double_t weight, Int_t type){
   EventInfo ev;
   ev.eventValue=val; ev.eventWeight=weight; ev.eventType=type;

   fEvtVector.push_back(ev);
   if (fIsSorted) fIsSorted=false;

}


Double_t TMVA::QuickMVAProbEstimator::GetMVAProbAt(Double_t value){
   // Well.. if it's fast is actually another question all together, merely
   // it's a quick and dirty simple kNN approach to the 1-Dim signal/backgr. MVA
   // distributions.


   if (!fIsSorted) {
      std::sort(fEvtVector.begin(),fEvtVector.end(),TMVA::QuickMVAProbEstimator::compare), fIsSorted=true;
   }

   Double_t     percentage = 0.1;
   UInt_t  nRange = TMath::Max(fNMin,(UInt_t) (fEvtVector.size() * percentage));
   nRange = TMath::Min(fNMax,nRange);
   // just make sure that nRange > you total number of events
   if (nRange > fEvtVector.size()) {
      nRange = fEvtVector.size()/3.;
      Log() << kWARNING  << " !!  you have only " << fEvtVector.size() << " of events.. . I choose "
            << nRange << " for the quick and dirty kNN MVAProb estimate" << Endl;
   }

   EventInfo tmp; tmp.eventValue=value;
   std::vector<EventInfo>::iterator it = std::upper_bound(fEvtVector.begin(),fEvtVector.end(),tmp,TMVA::QuickMVAProbEstimator::compare);

   UInt_t iLeft=0, iRight=0;
   Double_t nSignal=0;
   Double_t nBackgr=0;

   while ( (iLeft+iRight) < nRange){
      if ( fEvtVector.end() > it+iRight+1){
         iRight++;
         if ( ((it+iRight))->eventType == 0) nSignal+=((it+iRight))->eventWeight;
         else                                nBackgr+=((it+iRight))->eventWeight;
      }
      if ( fEvtVector.begin() <= it-iLeft-1){
         iLeft++;
         if ( ((it-iLeft))->eventType == 0) nSignal+=((it-iLeft))->eventWeight;
         else                               nBackgr+=((it-iLeft))->eventWeight;
      }
   }

   Double_t mvaProb = (nSignal+nBackgr) ? nSignal/(nSignal+nBackgr) : -1 ;
   return mvaProb;

}
