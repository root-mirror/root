// @(#)root/tmva $Id$
// Author: Andrzej Zemla

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : SVWorkingSet                                                          *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation                                                            *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Marcin Wolter  <Marcin.Wolter@cern.ch> - IFJ PAN, Krakow, Poland          *
 *      Andrzej Zemla  <azemla@cern.ch>        - IFJ PAN, Krakow, Poland          *
 *      (IFJ PAN: Henryk Niewodniczanski Inst. Nucl. Physics, Krakow, Poland)     *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *      PAN, Krakow, Poland                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::SVWorkingSet
\ingroup TMVA
Working class for Support Vector Machine
*/

#include "TMVA/SVWorkingSet.h"

#include "TMVA/MsgLogger.h"
#include "TMVA/SVEvent.h"
#include "TMVA/SVKernelFunction.h"
#include "TMVA/SVKernelMatrix.h"
#include "TMVA/Types.h"


#include "TMath.h"
#include "TRandom3.h"

#include <vector>

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::SVWorkingSet::SVWorkingSet()
   : fdoRegression(kFALSE),
     fInputData(0),
     fSupVec(0),
     fKFunction(0),
     fKMatrix(0),
     fTEventUp(0),
     fTEventLow(0),
     fB_low(1.),
     fB_up(-1.),
     fTolerance(0.01),
     fLogger( new MsgLogger( "SVWorkingSet", kINFO ) )
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::SVWorkingSet::SVWorkingSet(std::vector<TMVA::SVEvent*>*inputVectors, SVKernelFunction* kernelFunction,
                                 Float_t tol, Bool_t doreg)
   : fdoRegression(doreg),
     fInputData(inputVectors),
     fSupVec(0),
     fKFunction(kernelFunction),
     fTEventUp(0),
     fTEventLow(0),
     fB_low(1.),
     fB_up(-1.),
     fTolerance(tol),
     fLogger( new MsgLogger( "SVWorkingSet", kINFO ) )
{
   fKMatrix = new TMVA::SVKernelMatrix(inputVectors, kernelFunction);
   Float_t *pt;
   for( UInt_t i = 0; i < fInputData->size(); i++){
      pt = fKMatrix->GetLine(i);
      fInputData->at(i)->SetLine(pt);
      fInputData->at(i)->SetNs(i);
      if(fdoRegression) fInputData->at(i)->SetErrorCache(fInputData->at(i)->GetTarget());
   }
   TRandom3 rand;
   UInt_t kk = rand.Integer(fInputData->size());
   if(fdoRegression) {
      fTEventLow = fTEventUp =fInputData->at(0);
      fB_low = fTEventUp ->GetTarget() - fTolerance;
      fB_up  = fTEventLow->GetTarget() + fTolerance;
   }
   else{
      while(1){
         if(fInputData->at(kk)->GetTypeFlag()==-1){
            fTEventLow = fInputData->at(kk);
            break;
         }
         kk = rand.Integer(fInputData->size());
      }

      while (1){
         if (fInputData->at(kk)->GetTypeFlag()==1) {
            fTEventUp = fInputData->at(kk);
            break;
         }
         kk = rand.Integer(fInputData->size());
      }
   }
   fTEventUp ->SetErrorCache(fTEventUp->GetTarget());
   fTEventLow->SetErrorCache(fTEventUp->GetTarget());
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::SVWorkingSet::~SVWorkingSet()
{
   if (fKMatrix   != 0) {delete fKMatrix; fKMatrix = 0;}
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TMVA::SVWorkingSet::ExamineExample( TMVA::SVEvent* jevt )
{
   SVEvent* ievt=0;
   Float_t fErrorC_J = 0.;
   if( jevt->GetIdx()==0) fErrorC_J = jevt->GetErrorCache();
   else{
      Float_t *fKVals = jevt->GetLine();
      fErrorC_J = 0.;
      std::vector<TMVA::SVEvent*>::iterator idIter;

      UInt_t k=0;
      for(idIter = fInputData->begin(); idIter != fInputData->end(); ++idIter){
         if((*idIter)->GetAlpha()>0)
            fErrorC_J += (*idIter)->GetAlpha()*(*idIter)->GetTypeFlag()*fKVals[k];
         k++;
      }


      fErrorC_J -= jevt->GetTypeFlag();
      jevt->SetErrorCache(fErrorC_J);

      if((jevt->GetIdx() == 1) && (fErrorC_J < fB_up )){
         fB_up = fErrorC_J;
         fTEventUp = jevt;
      }
      else if ((jevt->GetIdx() == -1)&&(fErrorC_J > fB_low)) {
         fB_low = fErrorC_J;
         fTEventLow = jevt;
      }
   }
   Bool_t converged = kTRUE;

   if((jevt->GetIdx()>=0) && (fB_low - fErrorC_J > 2*fTolerance)) {
      converged = kFALSE;
      ievt = fTEventLow;
   }

   if((jevt->GetIdx()<=0) && (fErrorC_J - fB_up > 2*fTolerance)) {
      converged = kFALSE;
      ievt = fTEventUp;
   }

   if (converged) return kFALSE;

   if(jevt->GetIdx()==0){
      if(fB_low - fErrorC_J > fErrorC_J - fB_up) ievt = fTEventLow;
      else                                       ievt = fTEventUp;
   }

   if (TakeStep(ievt, jevt)) return kTRUE;
   else                      return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////

Bool_t TMVA::SVWorkingSet::TakeStep(TMVA::SVEvent* ievt,TMVA::SVEvent* jevt )
{
   if (ievt == jevt) return kFALSE;
   std::vector<TMVA::SVEvent*>::iterator idIter;
   const Float_t epsilon = 1e-8; //make it 1-e6 or 1-e5 to make it faster

   Float_t type_I,  type_J;
   Float_t errorC_I,  errorC_J;
   Float_t alpha_I, alpha_J;

   Float_t newAlpha_I, newAlpha_J;
   Int_t   s;

   Float_t l, h, lobj = 0, hobj = 0;
   Float_t eta;

   type_I   = ievt->GetTypeFlag();
   alpha_I  = ievt->GetAlpha();
   errorC_I = ievt->GetErrorCache();

   type_J   = jevt->GetTypeFlag();
   alpha_J  = jevt->GetAlpha();
   errorC_J = jevt->GetErrorCache();

   s = Int_t( type_I * type_J );

   Float_t c_i = ievt->GetCweight();

   Float_t c_j =  jevt->GetCweight();

   // compute l, h objective function

   if (type_I == type_J) {
      Float_t gamma = alpha_I + alpha_J;

      if ( c_i > c_j ) {
         if ( gamma < c_j ) {
            l = 0;
            h = gamma;
         }
         else{
            h = c_j;
            if ( gamma < c_i )
               l = 0;
            else
               l = gamma - c_i;
         }
      }
      else {
         if ( gamma < c_i ){
            l = 0;
            h = gamma;
         }
         else {
            l = gamma - c_i;
            if ( gamma < c_j )
               h = gamma;
            else
               h = c_j;
         }
      }
   }
   else {
      Float_t gamma = alpha_I - alpha_J;
      if (gamma > 0) {
         l = 0;
         if ( gamma >= (c_i - c_j) )
            h = c_i - gamma;
         else
            h = c_j;
      }
      else {
         l = -gamma;
         if ( (c_i - c_j) >= gamma)
            h = c_j;
         else
            h = c_i - gamma;
      }
   }

   if (l == h)  return kFALSE;
   Float_t kernel_II, kernel_IJ, kernel_JJ;

   kernel_II = fKMatrix->GetElement(ievt->GetNs(),ievt->GetNs());
   kernel_IJ = fKMatrix->GetElement(ievt->GetNs(), jevt->GetNs());
   kernel_JJ = fKMatrix->GetElement(jevt->GetNs(),jevt->GetNs());

   eta = 2*kernel_IJ - kernel_II - kernel_JJ;
   if (eta < 0) {
      newAlpha_J = alpha_J + (type_J*( errorC_J - errorC_I ))/eta;
      if      (newAlpha_J < l) newAlpha_J = l;
      else if (newAlpha_J > h) newAlpha_J = h;

   }

   else {

      Float_t c_I = eta/2;
      Float_t c_J = type_J*( errorC_I - errorC_J ) - eta * alpha_J;
      lobj = c_I * l * l + c_J * l;
      hobj = c_I * h * h + c_J * h;

      if      (lobj > hobj + epsilon)  newAlpha_J = l;
      else if (lobj < hobj - epsilon)  newAlpha_J = h;
      else                              newAlpha_J = alpha_J;
   }

   if (TMath::Abs( newAlpha_J - alpha_J ) < ( epsilon * ( newAlpha_J + alpha_J+ epsilon ))){
      return kFALSE;
      //it spends here to much time... it is stupido
   }
   newAlpha_I = alpha_I - s*( newAlpha_J - alpha_J );

   if (newAlpha_I < 0) {
      newAlpha_J += s* newAlpha_I;
      newAlpha_I = 0;
   }
   else if (newAlpha_I > c_i) {
      Float_t temp = newAlpha_I - c_i;
      newAlpha_J += s * temp;
      newAlpha_I = c_i;
   }

   Float_t dL_I = type_I * ( newAlpha_I - alpha_I );
   Float_t dL_J = type_J * ( newAlpha_J - alpha_J );

   Int_t k = 0;
   for(idIter = fInputData->begin(); idIter != fInputData->end(); ++idIter){
      k++;
      if((*idIter)->GetIdx()==0){
         Float_t ii = fKMatrix->GetElement(ievt->GetNs(), (*idIter)->GetNs());
         Float_t jj = fKMatrix->GetElement(jevt->GetNs(), (*idIter)->GetNs());

         (*idIter)->UpdateErrorCache(dL_I * ii + dL_J * jj);
      }
   }
   ievt->SetAlpha(newAlpha_I);
   jevt->SetAlpha(newAlpha_J);
   // set new indexes
   SetIndex(ievt);
   SetIndex(jevt);

   // update error cache
   ievt->SetErrorCache(errorC_I + dL_I*kernel_II + dL_J*kernel_IJ);
   jevt->SetErrorCache(errorC_J + dL_I*kernel_IJ + dL_J*kernel_JJ);

   // compute fI_low, fB_low

   fB_low = -1*1e30;
   fB_up = 1e30;

   for(idIter = fInputData->begin(); idIter != fInputData->end(); ++idIter){
      if((*idIter)->GetIdx()==0){
         if((*idIter)->GetErrorCache()> fB_low){
            fB_low = (*idIter)->GetErrorCache();
            fTEventLow = (*idIter);
         }
         if( (*idIter)->GetErrorCache()< fB_up){
            fB_up =(*idIter)->GetErrorCache();
            fTEventUp = (*idIter);
         }
      }
   }

   // for optimized alfa's
   if (fB_low < TMath::Max(ievt->GetErrorCache(), jevt->GetErrorCache())) {
      if (ievt->GetErrorCache() > fB_low) {
         fB_low = ievt->GetErrorCache();
         fTEventLow = ievt;
      }
      else {
         fB_low = jevt->GetErrorCache();
         fTEventLow = jevt;
      }
   }

   if (fB_up > TMath::Max(ievt->GetErrorCache(), jevt->GetErrorCache())) {
      if (ievt->GetErrorCache()< fB_low) {
         fB_up =ievt->GetErrorCache();
         fTEventUp = ievt;
      }
      else {
         fB_up =jevt->GetErrorCache() ;
         fTEventUp = jevt;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t  TMVA::SVWorkingSet::Terminated()
{
   if((fB_up > fB_low - 2*fTolerance)) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// train the SVM

void TMVA::SVWorkingSet::Train(UInt_t nMaxIter)
{

   Int_t numChanged  = 0;
   Int_t examineAll  = 1;

   Float_t numChangedOld = 0;
   Int_t deltaChanges = 0;
   UInt_t numit    = 0;

   std::vector<TMVA::SVEvent*>::iterator idIter;

   while ((numChanged > 0) || (examineAll > 0)) {
     if (fIPyCurrentIter) *fIPyCurrentIter = numit;
     if (fExitFromTraining && *fExitFromTraining) break;
      numChanged = 0;
      if (examineAll) {
         for (idIter = fInputData->begin(); idIter!=fInputData->end(); ++idIter){
            if(!fdoRegression) numChanged += (UInt_t)ExamineExample(*idIter);
            else numChanged += (UInt_t)ExamineExampleReg(*idIter);
         }
      }
      else {
         for (idIter = fInputData->begin(); idIter!=fInputData->end(); ++idIter) {
            if ((*idIter)->IsInI0()) {
               if(!fdoRegression) numChanged += (UInt_t)ExamineExample(*idIter);
               else numChanged += (UInt_t)ExamineExampleReg(*idIter);
               if (Terminated()) {
                  numChanged = 0;
                  break;
               }
            }
         }
      }

      if      (examineAll == 1) examineAll = 0;
      else if (numChanged == 0 || numChanged < 10 || deltaChanges > 3 ) examineAll = 1;

      if (numChanged == numChangedOld) deltaChanges++;
      else                             deltaChanges = 0;
      numChangedOld = numChanged;
      ++numit;

      if (numit >= nMaxIter) {
         *fLogger << kWARNING
                  << "Max number of iterations exceeded. "
                  << "Training may not be completed. Try use less Cost parameter" << Endl;
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::SVWorkingSet::SetIndex( TMVA::SVEvent* event )
{
   if( (0< event->GetAlpha()) && (event->GetAlpha()< event->GetCweight()))
      event->SetIdx(0);

   if( event->GetTypeFlag() == 1){
      if( event->GetAlpha() == 0)
         event->SetIdx(1);
      else if( event->GetAlpha() == event->GetCweight() )
         event->SetIdx(-1);
   }
   if( event->GetTypeFlag() == -1){
      if( event->GetAlpha() == 0)
         event->SetIdx(-1);
      else if( event->GetAlpha() == event->GetCweight() )
         event->SetIdx(1);
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::SVWorkingSet::PrintStat()
{
   std::vector<TMVA::SVEvent*>::iterator idIter;
   UInt_t counter = 0;
   for( idIter = fInputData->begin(); idIter != fInputData->end(); ++idIter)
      if((*idIter)->GetAlpha() !=0) counter++;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<TMVA::SVEvent*>* TMVA::SVWorkingSet::GetSupportVectors()
{
   std::vector<TMVA::SVEvent*>::iterator idIter;
   if( fSupVec != 0) {delete fSupVec; fSupVec = 0; }
   fSupVec = new std::vector<TMVA::SVEvent*>(0);

   for( idIter = fInputData->begin(); idIter != fInputData->end(); ++idIter){
      if((*idIter)->GetDeltaAlpha() !=0){
         fSupVec->push_back((*idIter));
      }
   }
   return fSupVec;
}

//for regression

Bool_t TMVA::SVWorkingSet::TakeStepReg(TMVA::SVEvent* ievt,TMVA::SVEvent* jevt )
{
   if (ievt == jevt) return kFALSE;
   std::vector<TMVA::SVEvent*>::iterator idIter;
   const Float_t epsilon = 0.001*fTolerance;//TODO

   const Float_t kernel_II = fKMatrix->GetElement(ievt->GetNs(),ievt->GetNs());
   const Float_t kernel_IJ = fKMatrix->GetElement(ievt->GetNs(),jevt->GetNs());
   const Float_t kernel_JJ = fKMatrix->GetElement(jevt->GetNs(),jevt->GetNs());

   //compute eta & gamma
   const Float_t eta = -2*kernel_IJ + kernel_II + kernel_JJ;
   const Float_t gamma = ievt->GetDeltaAlpha() + jevt->GetDeltaAlpha();

   //TODO CHECK WHAT IF ETA <0
   //w.r.t Mercer's conditions it should never happen, but what if?

   Bool_t caseA, caseB, caseC, caseD, terminated;
   caseA = caseB = caseC = caseD = terminated = kFALSE;
   Float_t b_alpha_i, b_alpha_j, b_alpha_i_p, b_alpha_j_p; //temporary Lagrange multipliers
   const Float_t b_cost_i = ievt->GetCweight();
   const Float_t b_cost_j = jevt->GetCweight();

   b_alpha_i   = ievt->GetAlpha();
   b_alpha_j   = jevt->GetAlpha();
   b_alpha_i_p = ievt->GetAlpha_p();
   b_alpha_j_p = jevt->GetAlpha_p();

   //calculate deltafi
   Float_t deltafi = ievt->GetErrorCache()-jevt->GetErrorCache();

   // main loop
   while(!terminated) {
      const Float_t null = 0.; //!!! dummy float null declaration because of problems with TMath::Max/Min(Float_t, Float_t) function
      Float_t low, high;
      Float_t tmp_alpha_i, tmp_alpha_j;
      tmp_alpha_i = tmp_alpha_j = 0.;

      //TODO check this conditions, are they proper
      if((caseA == kFALSE) && (b_alpha_i > 0 || (b_alpha_i_p == 0 && deltafi > 0)) && (b_alpha_j > 0 || (b_alpha_j_p == 0 && deltafi < 0)))
         {
            //compute low, high w.r.t a_i, a_j
            low  = TMath::Max( null, gamma - b_cost_j );
            high = TMath::Min( b_cost_i , gamma);

            if(low<high){
               tmp_alpha_j = b_alpha_j - (deltafi/eta);
               tmp_alpha_j = TMath::Min(tmp_alpha_j,high      );
               tmp_alpha_j = TMath::Max(low        ,tmp_alpha_j);
               tmp_alpha_i = b_alpha_i - (tmp_alpha_j - b_alpha_j);

               //update Li & Lj if change is significant (??)
               if( IsDiffSignificant(b_alpha_j,tmp_alpha_j, epsilon) ||  IsDiffSignificant(b_alpha_i,tmp_alpha_i, epsilon)){
                  b_alpha_j = tmp_alpha_j;
                  b_alpha_i = tmp_alpha_i;
               }

            }
            else
               terminated = kTRUE;

            caseA = kTRUE;
         }
      else if((caseB==kFALSE) && (b_alpha_i>0 || (b_alpha_i_p==0 && deltafi >2*epsilon )) && (b_alpha_j_p>0 || (b_alpha_j==0 && deltafi>2*epsilon)))
         {
            //compute LH w.r.t. a_i, a_j*
            low  = TMath::Max( null, gamma );  //TODO
            high = TMath::Min( b_cost_i , b_cost_j + gamma);


            if(low<high){
               tmp_alpha_j = b_alpha_j_p - ((deltafi-2*epsilon)/eta);
               tmp_alpha_j = TMath::Min(tmp_alpha_j,high);
               tmp_alpha_j = TMath::Max(low,tmp_alpha_j);
               tmp_alpha_i = b_alpha_i - (tmp_alpha_j - b_alpha_j_p);

               //update alphai alphaj_p
               if( IsDiffSignificant(b_alpha_j_p,tmp_alpha_j, epsilon) ||  IsDiffSignificant(b_alpha_i,tmp_alpha_i, epsilon)){
                  b_alpha_j_p = tmp_alpha_j;
                  b_alpha_i   = tmp_alpha_i;
               }
            }
            else
               terminated = kTRUE;

            caseB = kTRUE;
         }
      else if((caseC==kFALSE) && (b_alpha_i_p>0 || (b_alpha_i==0 && deltafi < -2*epsilon )) && (b_alpha_j>0 || (b_alpha_j_p==0 && deltafi< -2*epsilon)))
         {
            //compute LH w.r.t. alphai_p alphaj
            low  = TMath::Max(null, -gamma  );
            high = TMath::Min(b_cost_i, -gamma+b_cost_j);

            if(low<high){
               tmp_alpha_j = b_alpha_j - ((deltafi+2*epsilon)/eta);
               tmp_alpha_j = TMath::Min(tmp_alpha_j,high      );
               tmp_alpha_j = TMath::Max(low        ,tmp_alpha_j);
               tmp_alpha_i = b_alpha_i_p - (tmp_alpha_j - b_alpha_j);

               //update alphai_p alphaj
               if( IsDiffSignificant(b_alpha_j,tmp_alpha_j, epsilon) ||  IsDiffSignificant(b_alpha_i_p,tmp_alpha_i, epsilon)){
                  b_alpha_j     = tmp_alpha_j;
                  b_alpha_i_p   = tmp_alpha_i;
               }
            }
            else
               terminated = kTRUE;

            caseC = kTRUE;
         }
      else if((caseD == kFALSE) &&
              (b_alpha_i_p>0 || (b_alpha_i==0 && deltafi <0 )) &&
              (b_alpha_j_p>0 || (b_alpha_j==0 && deltafi >0 )))
         {
            //compute LH w.r.t. alphai_p alphaj_p
            low  = TMath::Max(null,-gamma - b_cost_j);
            high = TMath::Min(b_cost_i, -gamma);

            if(low<high){
               tmp_alpha_j = b_alpha_j_p + (deltafi/eta);
               tmp_alpha_j = TMath::Min(tmp_alpha_j,high      );
               tmp_alpha_j = TMath::Max(low        ,tmp_alpha_j);
               tmp_alpha_i = b_alpha_i_p - (tmp_alpha_j - b_alpha_j_p);

               if( IsDiffSignificant(b_alpha_j_p,tmp_alpha_j, epsilon) ||  IsDiffSignificant(b_alpha_i_p,tmp_alpha_i, epsilon)){
                  b_alpha_j_p   = tmp_alpha_j;
                  b_alpha_i_p   = tmp_alpha_i;
               }
            }
            else
               terminated = kTRUE;

            caseD = kTRUE;
         }
      else
         terminated = kTRUE;
   }
   // TODO ad commment how it was calculated
   deltafi += ievt->GetDeltaAlpha()*(kernel_II - kernel_IJ) + jevt->GetDeltaAlpha()*(kernel_IJ - kernel_JJ);

   if( IsDiffSignificant(b_alpha_i, ievt->GetAlpha(), epsilon) ||
       IsDiffSignificant(b_alpha_j, jevt->GetAlpha(), epsilon) ||
       IsDiffSignificant(b_alpha_i_p, ievt->GetAlpha_p(), epsilon) ||
       IsDiffSignificant(b_alpha_j_p, jevt->GetAlpha_p(), epsilon) ){

      //TODO check if these conditions might be easier
      //TODO write documentation for this
      const Float_t diff_alpha_i = ievt->GetDeltaAlpha()+b_alpha_i_p - ievt->GetAlpha();
      const Float_t diff_alpha_j = jevt->GetDeltaAlpha()+b_alpha_j_p - jevt->GetAlpha();

      //update error cache
      Int_t k = 0;
      for(idIter = fInputData->begin(); idIter != fInputData->end(); ++idIter){
         k++;
         //there will be some changes in Idx notation
         if((*idIter)->GetIdx()==0){
            Float_t k_ii = fKMatrix->GetElement(ievt->GetNs(), (*idIter)->GetNs());
            Float_t k_jj = fKMatrix->GetElement(jevt->GetNs(), (*idIter)->GetNs());

            (*idIter)->UpdateErrorCache(diff_alpha_i * k_ii + diff_alpha_j * k_jj);
         }
      }

      //store new alphas in SVevents
      ievt->SetAlpha(b_alpha_i);
      jevt->SetAlpha(b_alpha_j);
      ievt->SetAlpha_p(b_alpha_i_p);
      jevt->SetAlpha_p(b_alpha_j_p);

      //TODO update Idexes

      // compute fI_low, fB_low

      fB_low = -1*1e30;
      fB_up =1e30;

      for(idIter = fInputData->begin(); idIter != fInputData->end(); ++idIter){
         if((!(*idIter)->IsInI3()) && ((*idIter)->GetErrorCache()> fB_low)){
            fB_low = (*idIter)->GetErrorCache();
            fTEventLow = (*idIter);

         }
         if((!(*idIter)->IsInI2()) && ((*idIter)->GetErrorCache()< fB_up)){
            fB_up =(*idIter)->GetErrorCache();
            fTEventUp = (*idIter);
         }
      }
      return kTRUE;
   } else return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////

Bool_t TMVA::SVWorkingSet::ExamineExampleReg(TMVA::SVEvent* jevt)
{
   Float_t feps = 1e-7;// TODO check which value is the best
   SVEvent* ievt=0;
   Float_t fErrorC_J = 0.;
   if( jevt->IsInI0()) {
      fErrorC_J = jevt->GetErrorCache();
   }
   else{
      Float_t *fKVals = jevt->GetLine();
      fErrorC_J = 0.;
      std::vector<TMVA::SVEvent*>::iterator idIter;

      UInt_t k=0;
      for(idIter = fInputData->begin(); idIter != fInputData->end(); ++idIter){
         fErrorC_J -= (*idIter)->GetDeltaAlpha()*fKVals[k];
         k++;
      }

      fErrorC_J += jevt->GetTarget();
      jevt->SetErrorCache(fErrorC_J);

      if(jevt->IsInI1()){
         if(fErrorC_J + feps < fB_up ){
            fB_up = fErrorC_J + feps;
            fTEventUp = jevt;
         }
         else if(fErrorC_J -feps > fB_low) {
            fB_low = fErrorC_J - feps;
            fTEventLow = jevt;
         }
      }else if((jevt->IsInI2()) && (fErrorC_J + feps > fB_low)){
         fB_low = fErrorC_J + feps;
         fTEventLow = jevt;
      }else if((jevt->IsInI3()) && (fErrorC_J - feps < fB_up)){
         fB_up = fErrorC_J - feps;
         fTEventUp = jevt;
      }
   }

   Bool_t converged = kTRUE;
   //case 1
   if(jevt->IsInI0a()){
      if( fB_low -fErrorC_J + feps > 2*fTolerance){
         converged = kFALSE;
         ievt = fTEventLow;
         if(fErrorC_J-feps-fB_up > fB_low-fErrorC_J+feps){
            ievt = fTEventUp;
         }
      }else if(fErrorC_J -feps - fB_up > 2*fTolerance){
         converged = kFALSE;
         ievt = fTEventUp;
         if(fB_low - fErrorC_J+feps > fErrorC_J-feps -fB_up){
            ievt = fTEventLow;
         }
      }
   }

   //case 2
   if(jevt->IsInI0b()){
      if( fB_low -fErrorC_J - feps > 2*fTolerance){
         converged = kFALSE;
         ievt = fTEventLow;
         if(fErrorC_J+feps-fB_up > fB_low-fErrorC_J-feps){
            ievt = fTEventUp;
         }
      }else if(fErrorC_J + feps - fB_up > 2*fTolerance){
         converged = kFALSE;
         ievt = fTEventUp;
         if(fB_low - fErrorC_J-feps > fErrorC_J+feps -fB_up){
            ievt = fTEventLow;
         }
      }
   }

   //case 3
   if(jevt->IsInI1()){
      if( fB_low -fErrorC_J - feps > 2*fTolerance){
         converged = kFALSE;
         ievt = fTEventLow;
         if(fErrorC_J+feps-fB_up > fB_low-fErrorC_J-feps){
            ievt = fTEventUp;
         }
      }else if(fErrorC_J - feps - fB_up > 2*fTolerance){
         converged = kFALSE;
         ievt = fTEventUp;
         if(fB_low - fErrorC_J+feps > fErrorC_J-feps -fB_up){
            ievt = fTEventLow;
         }
      }
   }

   //case 4
   if(jevt->IsInI2()){
      if( fErrorC_J + feps -fB_up > 2*fTolerance){
         converged = kFALSE;
         ievt = fTEventUp;
      }
   }

   //case 5
   if(jevt->IsInI3()){
      if(fB_low -fErrorC_J +feps > 2*fTolerance){
         converged = kFALSE;
         ievt = fTEventLow;
      }
   }

   if(converged) return kFALSE;
   if (TakeStepReg(ievt, jevt)) return kTRUE;
   else return kFALSE;
}

Bool_t TMVA::SVWorkingSet::IsDiffSignificant(Float_t a_i, Float_t a_j, Float_t eps)
{
   if( TMath::Abs(a_i - a_j) > eps*(a_i + a_j + eps)) return kTRUE;
   else return kFALSE;
}

