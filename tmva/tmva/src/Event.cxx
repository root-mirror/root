// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss, Jan Therhaag

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Event                                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch> - CERN, Switzerland           *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Jan Therhaag       <Jan.Therhaag@cern.ch>     - U of Bonn, Germany        *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://mva.sourceforge.net/license.txt)                                       *
 **********************************************************************************/

#include "TMVA/Event.h"
#include "TMVA/Tools.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include "TCut.h"

/*! \class TMVA::Event
\ingroup TMVA
*/

Bool_t TMVA::Event::fgIsTraining = kFALSE;
Bool_t TMVA::Event::fgIgnoreNegWeightsInTraining = kFALSE;

////////////////////////////////////////////////////////////////////////////////
/// copy constructor

TMVA::Event::Event()
   : fValues(),
     fValuesDynamic(0),
     fTargets(),
     fSpectators(),
     fVariableArrangement(0),
     fClass(0),
     fWeight(1.0),
     fBoostWeight(1.0),
     fDynamic(kFALSE),
     fDoNotBoost(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::Event::Event( const std::vector<Float_t>& ev,
                    const std::vector<Float_t>& tg,
                    UInt_t cls,
                    Double_t weight,
                    Double_t boostweight )
   : fValues(ev),
     fValuesDynamic(0),
     fTargets(tg),
     fSpectators(0),
     fVariableArrangement(0),
     fClass(cls),
     fWeight(weight),
     fBoostWeight(boostweight),
     fDynamic(kFALSE),
     fDoNotBoost(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::Event::Event( const std::vector<Float_t>& ev,
                    const std::vector<Float_t>& tg,
                    const std::vector<Float_t>& vi,
                    UInt_t cls,
                    Double_t weight,
                    Double_t boostweight )
   : fValues(ev),
     fValuesDynamic(0),
     fTargets(tg),
     fSpectators(vi),
     fVariableArrangement(0),
     fClass(cls),
     fWeight(weight),
     fBoostWeight(boostweight),
     fDynamic(kFALSE),
     fDoNotBoost(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::Event::Event( const std::vector<Float_t>& ev,
                    UInt_t cls,
                    Double_t weight,
                    Double_t boostweight )
   : fValues(ev),
     fValuesDynamic(0),
     fTargets(0),
     fSpectators(0),
     fVariableArrangement(0),
     fClass(cls),
     fWeight(weight),
     fBoostWeight(boostweight),
     fDynamic(kFALSE),
     fDoNotBoost(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor for single events

TMVA::Event::Event( const std::vector<Float_t*>*& evdyn, UInt_t nvar )
   : fValues(nvar),
     fValuesDynamic(0),
     fTargets(0),
     fSpectators(evdyn->size()-nvar),
     fVariableArrangement(0),
     fClass(0),
     fWeight(0),
     fBoostWeight(0),
     fDynamic(true),
     fDoNotBoost(kFALSE)
{
   fValuesDynamic = (std::vector<Float_t*>*) evdyn;
}

////////////////////////////////////////////////////////////////////////////////
/// copy constructor

TMVA::Event::Event( const Event& event )
   : TObject(event),
     fValues(event.fValues),
     fValuesDynamic(event.fValuesDynamic),
     fTargets(event.fTargets),
     fSpectators(event.fSpectators),
     fVariableArrangement(event.fVariableArrangement),
     fClass(event.fClass),
     fWeight(event.fWeight),
     fBoostWeight(event.fBoostWeight),
     fDynamic(event.fDynamic),
     fDoNotBoost(kFALSE)
{
   if (event.fDynamic){
      fValues.clear();
      UInt_t nvar = event.GetNVariables();
      UInt_t idx=0;
      std::vector<Float_t*>::iterator itDyn=event.fValuesDynamic->begin(), itDynEnd=event.fValuesDynamic->end();
      for (; itDyn!=itDynEnd && idx<nvar; ++itDyn){
         Float_t value=*(*itDyn);
         fValues.push_back( value );
         ++idx;
      }
      fSpectators.clear();
      for (; itDyn!=itDynEnd; ++itDyn){
         Float_t value=*(*itDyn);
         fSpectators.push_back( value );
         ++idx;
      }

      fDynamic=kFALSE;
      fValuesDynamic=NULL;
}
}

////////////////////////////////////////////////////////////////////////////////
/// Event destructor

TMVA::Event::~Event()
{
//    delete fValuesDynamic;
}
////////////////////////////////////////////////////////////////////////////////
/// set the variable arrangement

void TMVA::Event::SetVariableArrangement( std::vector<UInt_t>* const m ) const {
   // mapping from global variable index (the position in the vector)
   // to the new index in the subset of variables used by the
   // composite classifier
    if(!m)fVariableArrangement.clear();
    else fVariableArrangement = *m;
}


////////////////////////////////////////////////////////////////////////////////
/// copies only the variable values

void TMVA::Event::CopyVarValues( const Event& other )
{
   fValues      = other.fValues;
   fTargets     = other.fTargets;
   fSpectators  = other.fSpectators;
   if (other.fDynamic){
      UInt_t nvar = other.GetNVariables();
      fValues.clear();
      UInt_t idx=0;
      std::vector<Float_t*>::iterator itDyn=other.fValuesDynamic->begin(), itDynEnd=other.fValuesDynamic->end();
      for (; itDyn!=itDynEnd && idx<nvar; ++itDyn){
         Float_t value=*(*itDyn);
         fValues.push_back( value );
         ++idx;
      }
      fSpectators.clear();
      for (; itDyn!=itDynEnd; ++itDyn){
         Float_t value=*(*itDyn);
         fSpectators.push_back( value );
         ++idx;
      }
   }
   fDynamic     = kFALSE;
   fValuesDynamic = NULL;

   fClass       = other.fClass;
   fWeight      = other.fWeight;
   fBoostWeight = other.fBoostWeight;
}

////////////////////////////////////////////////////////////////////////////////
/// return value of i'th variable

Float_t TMVA::Event::GetValue( UInt_t ivar ) const
{
   Float_t retval;
   if (fVariableArrangement.size()==0) {
      retval = fDynamic ? ( *((fValuesDynamic)->at(ivar)) ) : fValues.at(ivar);
   }
   else {
      UInt_t mapIdx = fVariableArrangement[ivar];
      //   std::cout<< fDynamic ;
      if (fDynamic){
         //     std::cout<< " " << (*fValuesDynamic).size() << " " << fValues.size() << std::endl;
         retval = *((fValuesDynamic)->at(mapIdx));
      }
      else{
         //retval = fValues.at(ivar);
         retval = ( mapIdx<fValues.size() ) ? fValues[mapIdx] : fSpectators[mapIdx-fValues.size()];
      }
   }

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// return spectator content

Float_t TMVA::Event::GetSpectator( UInt_t ivar) const
{
   if (fDynamic) return *(fValuesDynamic->at(GetNVariables()+ivar));
   else          return fSpectators.at(ivar);
}

////////////////////////////////////////////////////////////////////////////////
/// return value vector

const std::vector<Float_t>& TMVA::Event::GetValues() const
{
   if (fVariableArrangement.size()==0) {

      if (fDynamic) {
         fValues.clear();
         for (std::vector<Float_t*>::const_iterator it = fValuesDynamic->begin(), itEnd=fValuesDynamic->end()-GetNSpectators();
              it != itEnd; ++it) {
            Float_t val = *(*it);
            fValues.push_back( val );
         }
      }
   }else{
      UInt_t mapIdx;
      if (fDynamic) {
         fValues.clear();
         for (UInt_t i=0; i< fVariableArrangement.size(); i++){
            mapIdx = fVariableArrangement[i];
            fValues.push_back(*((fValuesDynamic)->at(mapIdx)));
         }
      } else {
         // hmm now you have a problem, as you do not want to mess with the original event variables
         // (change them permanently) ... guess the only way is to add a 'fValuesRearranged' array,
         // and living with the fact that it 'doubles' the Event size :(
         fValuesRearranged.clear();
         for (UInt_t i=0; i< fVariableArrangement.size(); i++){
            mapIdx = fVariableArrangement[i];
            fValuesRearranged.push_back(fValues.at(mapIdx));
         }
         return fValuesRearranged;
      }
   }
   return fValues;
}

////////////////////////////////////////////////////////////////////////////////
/// accessor to the number of variables

UInt_t TMVA::Event::GetNVariables() const
{
   // if variables have to arranged (as it is the case for the
   // composite classifier) the number of the variables changes
   if (fVariableArrangement.size()==0) return fValues.size();
   else                         return fVariableArrangement.size();
}

////////////////////////////////////////////////////////////////////////////////
/// accessor to the number of targets

UInt_t TMVA::Event::GetNTargets() const
{
   return fTargets.size();
}

////////////////////////////////////////////////////////////////////////////////
/// accessor to the number of spectators

UInt_t TMVA::Event::GetNSpectators() const
{
   // if variables have to arranged (as it is the case for the
   // composite classifier) the number of the variables changes

   if (fVariableArrangement.size()==0) return fSpectators.size();
   else                         return fValues.size()-fVariableArrangement.size();
}


////////////////////////////////////////////////////////////////////////////////
/// set variable ivar to val

void TMVA::Event::SetVal( UInt_t ivar, Float_t val )
{
    if ((fDynamic ?( (*fValuesDynamic).size() ) : fValues.size())<=ivar)
        (fDynamic ?( (*fValuesDynamic).resize(ivar+1) ) : fValues.resize(ivar+1));

    (fDynamic ?( *(*fValuesDynamic)[ivar] ) : fValues[ivar])=val;
}

////////////////////////////////////////////////////////////////////////////////
/// print method

void TMVA::Event::Print( std::ostream& o ) const
{
   o << *this << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// set the target value (dimension itgt) to value

void TMVA::Event::SetTarget( UInt_t itgt, Float_t value )
{
   if (fTargets.size() <= itgt) fTargets.resize( itgt+1 );
   fTargets.at(itgt) = value;
}

////////////////////////////////////////////////////////////////////////////////
/// set spectator value (dimension ivar) to value

void TMVA::Event::SetSpectator( UInt_t ivar, Float_t value )
{
   if (fSpectators.size() <= ivar) fSpectators.resize( ivar+1 );
   fSpectators.at(ivar) = value;
}

////////////////////////////////////////////////////////////////////////////////
/// return the event weight - depending on whether the flag
/// *IgnoreNegWeightsInTraining* is or not. If it is set AND it is
/// used for training, then negative event weights are set to zero !
/// NOTE! For events used in Testing, the ORIGINAL possibly negative
/// event weight is used  no matter what

Double_t TMVA::Event::GetWeight() const
{
   return (fgIgnoreNegWeightsInTraining && fgIsTraining && fWeight < 0) ? 0. : fWeight*fBoostWeight;
}

////////////////////////////////////////////////////////////////////////////////
/// when this static function is called, it sets the flag whether
/// events with negative event weight should be ignored in the
/// training, or not.

void TMVA::Event::SetIsTraining(Bool_t b)
{
   fgIsTraining=b;
}
////////////////////////////////////////////////////////////////////////////////
/// when this static function is called, it sets the flag whether
/// events with negative event weight should be ignored in the
/// training, or not.

void TMVA::Event::SetIgnoreNegWeightsInTraining(Bool_t b)
{
   fgIgnoreNegWeightsInTraining=b;
}

////////////////////////////////////////////////////////////////////////////////
/// Outputs the data of an event

std::ostream& TMVA::operator << ( std::ostream& os, const TMVA::Event& event )
{
   os << "Variables [" << event.fValues.size() << "]:";
   for (UInt_t ivar=0; ivar<event.fValues.size(); ++ivar)
      os << " " << std::setw(10) << event.GetValue(ivar);
   os << ", targets [" << event.fTargets.size() << "]:";
   for (UInt_t ivar=0; ivar<event.fTargets.size(); ++ivar)
      os << " " << std::setw(10) << event.GetTarget(ivar);
   os << ", spectators ["<< event.fSpectators.size() << "]:";
   for (UInt_t ivar=0; ivar<event.fSpectators.size(); ++ivar)
      os << " " << std::setw(10) << event.GetSpectator(ivar);
   os << ", weight: " << event.GetWeight();
   os << ", class: " << event.GetClass();
   return os;
}
