// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Peter Speckmayer, Eckhard von Toerne

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VariableNormalizeTransform                                            *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker  <Andreas.Hocker@cern.ch>   - CERN, Switzerland           *
 *      Joerg Stelzer    <Joerg.Stelzer@cern.ch>    - CERN, Switzerland           *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch> - CERN, Switzerland           *
 *      Helge Voss       <Helge.Voss@cern.ch>       - MPI-K Heidelberg, Germany   *
 *      Eckhard v. Toerne     <evt@uni-bonn.de>        - U of Bonn, Germany       *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::VariableNormalizeTransform
\ingroup TMVA
Linear interpolation class
*/

#include "TMVA/VariableNormalizeTransform.h"

#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Event.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/VariableInfo.h"

#include "TMatrixD.h"
#include "TMatrixDBase.h"
#include "TVectorF.h"

#include <iostream>
#include <iomanip>
#include <cfloat>

ClassImp(TMVA::VariableNormalizeTransform);

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::VariableNormalizeTransform::VariableNormalizeTransform( DataSetInfo& dsi, TString strcor )
: VariableTransformBase( dsi, Types::kNormalized, "Norm" ),
   fNoOffset(kFALSE)
{
   if (strcor=="Scale") {fNoOffset = kTRUE;
      SetName("Scale");
   }
}

////////////////////////////////////////////////////////////////////////////////

TMVA::VariableNormalizeTransform::~VariableNormalizeTransform() {
}

////////////////////////////////////////////////////////////////////////////////
/// initialization of the normalization transformation

void TMVA::VariableNormalizeTransform::Initialize()
{
   UInt_t inputSize = fGet.size();
   Int_t numC = GetNClasses()+1;
   if (GetNClasses() <= 1 ) numC = 1;

   fMin.resize( numC );
   fMax.resize( numC );
   for (Int_t i=0; i<numC; i++) {
      fMin.at(i).resize(inputSize);
      fMax.at(i).resize(inputSize);
      fMin.at(i).assign(inputSize, 0);
      fMax.at(i).assign(inputSize, 0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// prepare transformation

Bool_t TMVA::VariableNormalizeTransform::PrepareTransformation (const std::vector<Event*>& events)
{
   if (!IsEnabled() || IsCreated()) return kTRUE;

   Log() << kDEBUG << "\tPreparing the transformation." << Endl;

   Initialize();

   CalcNormalizationParams( events );

   SetCreated( kTRUE );

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// apply the normalization transformation

const TMVA::Event* TMVA::VariableNormalizeTransform::Transform( const TMVA::Event* const ev, Int_t cls ) const
{
   if (!IsCreated()) Log() << kFATAL << "Transformation not yet created" << Endl;

   // if cls (the class chosen by the user) not existing,
   // assume that he wants to have the matrix for all classes together.
   // if (cls < 0 || cls > GetNClasses()) {
   //       if (GetNClasses() > 1 ) cls = GetNClasses();
   //       else cls = (fMin.size()==1?0:2);
   //    }
   // EVT this is a workaround to address the reader problem with transforma and EvaluateMVA(std::vector<float/double> ,...)
   if (cls < 0 || cls >= (int) fMin.size()) cls = fMin.size()-1;
   // EVT workaround end

   FloatVector input; // will be filled with the selected variables, targets, (spectators)
   FloatVector output; // will be filled with the selected variables, targets, (spectators)
   std::vector<Char_t> mask; // entries with kTRUE must not be transformed
   GetInput( ev, input, mask );

   if (fTransformedEvent==0) fTransformedEvent = new Event();

   Float_t min,max;
   const FloatVector& minVector = fMin.at(cls);
   const FloatVector& maxVector = fMax.at(cls);

   UInt_t iidx = 0;
   std::vector<Char_t>::iterator itMask = mask.begin();
   for ( std::vector<Float_t>::iterator itInp = input.begin(), itInpEnd = input.end(); itInp != itInpEnd; ++itInp) { // loop over input variables
      if( (*itMask) ){
         ++iidx;
         ++itMask;
         // don't put any value into output if the value is masked
         continue;
      }

      Float_t val = (*itInp);

      min = minVector.at(iidx);
      max = maxVector.at(iidx);

      Float_t valnorm;
      if (!fNoOffset) {
         Float_t offset = min;
         Float_t scale  = 1.0/(max-min);
         valnorm        = (val-offset)*scale * 2 - 1;
      } else {
         fabs(max)>fabs(min) ? valnorm=val/fabs(max) : valnorm=val/fabs(min);
      }

      output.push_back( valnorm );

      ++iidx;
      ++itMask;
   }

   SetOutput( fTransformedEvent, output, mask, ev );
   return fTransformedEvent;
}

////////////////////////////////////////////////////////////////////////////////
/// apply the inverse transformation

const TMVA::Event* TMVA::VariableNormalizeTransform::InverseTransform(const TMVA::Event* const ev, Int_t cls ) const
{
   if (!IsCreated()) Log() << kFATAL << "Transformation not yet created" << Endl;

   // if cls (the class chosen by the user) not existing,
   // assume that user wants to have the transformation for all classes together.
   if (cls < 0 || cls > GetNClasses()) {
      if (GetNClasses() > 1 ) cls = GetNClasses();
      else cls = 0;
   }

   FloatVector input;  // will be filled with the selected variables, targets, (spectators)
   FloatVector output; // will be filled with the output
   std::vector<Char_t> mask;
   GetInput( ev, input, mask, kTRUE );

   if (fBackTransformedEvent==0) fBackTransformedEvent = new Event( *ev );

   Float_t min,max;
   const FloatVector& minVector = fMin.at(cls);
   const FloatVector& maxVector = fMax.at(cls);

   UInt_t iidx = 0;
   for ( std::vector<Float_t>::iterator itInp = input.begin(), itInpEnd = input.end(); itInp != itInpEnd; ++itInp) { // loop over input variables
      Float_t val = (*itInp);

      min = minVector.at(iidx);
      max = maxVector.at(iidx);

      Float_t valnorm;
      if (!fNoOffset) {
         Float_t offset = min;
         Float_t scale  = 1.0/(max-min);
         valnorm = offset+((val+1)/(scale * 2));
      } else {
         fabs(max)>fabs(min) ? valnorm=val*fabs(max) : valnorm=val*fabs(min);
      }

      output.push_back( valnorm );

      ++iidx;
   }

   SetOutput( fBackTransformedEvent, output, mask, ev, kTRUE );

   return fBackTransformedEvent;
}

////////////////////////////////////////////////////////////////////////////////
/// compute offset and scale from min and max

void TMVA::VariableNormalizeTransform::CalcNormalizationParams( const std::vector< Event*>& events )
{
   if (events.size() <= 1)
      Log() << kFATAL << "Not enough events (found " << events.size() << ") to calculate the normalization" << Endl;

   FloatVector input; // will be filled with the selected variables, targets, (spectators)
   std::vector<Char_t> mask;

   UInt_t inputSize = fGet.size(); // number of input variables

   const UInt_t nCls = GetNClasses();
   Int_t numC = nCls+1;   // prepare the min and max values for each of the classes and additionally for all classes (if more than one)
   Int_t all = nCls; // at idx the min and max values for "all" classes are stored
   if (nCls <= 1 ) {
      numC = 1;
      all = 0;
   }

   for (UInt_t iinp=0; iinp<inputSize; ++iinp) {
      for (Int_t ic = 0; ic < numC; ic++) {
         fMin.at(ic).at(iinp) = FLT_MAX;
         fMax.at(ic).at(iinp) = -FLT_MAX;
      }
   }

   std::vector<Event*>::const_iterator evIt = events.begin();
   for (;evIt!=events.end();++evIt) { // loop over all events
      const TMVA::Event* event = (*evIt);   // get the event

      UInt_t cls = (*evIt)->GetClass(); // get the class of this event

      FloatVector& minVector = fMin.at(cls);
      FloatVector& maxVector = fMax.at(cls);

      FloatVector& minVectorAll = fMin.at(all);
      FloatVector& maxVectorAll = fMax.at(all);

      GetInput(event,input,mask);    // select the input variables for the transformation and get them from the event
      UInt_t iidx = 0;
      for ( std::vector<Float_t>::iterator itInp = input.begin(), itInpEnd = input.end(); itInp != itInpEnd; ++itInp) { // loop over input variables
         Float_t val = (*itInp);

         if( minVector.at(iidx) > val ) minVector.at(iidx) = val;
         if( maxVector.at(iidx) < val ) maxVector.at(iidx) = val;

         if (nCls != 1) { // in case more than one class exists, compute min and max as well for all classes together
            if (minVectorAll.at(iidx) > val) minVectorAll.at(iidx) = val;
            if (maxVectorAll.at(iidx) < val) maxVectorAll.at(iidx) = val;
         }

         ++iidx;
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// creates string with variable transformations applied

std::vector<TString>* TMVA::VariableNormalizeTransform::GetTransformationStrings( Int_t cls ) const
{
   // if cls (the class chosen by the user) not existing, assume that user wants to
   // have the matrix for all classes together.
   if (cls < 0 || cls > GetNClasses()) cls = GetNClasses();

   Float_t min, max;
   const UInt_t size = fGet.size();
   std::vector<TString>* strVec = new std::vector<TString>(size);

   UInt_t iinp = 0;
   for( ItVarTypeIdxConst itGet = fGet.begin(), itGetEnd = fGet.end(); itGet != itGetEnd; ++itGet ) {
      min = fMin.at(cls).at(iinp);
      max = fMax.at(cls).at(iinp);

      Char_t type = (*itGet).first;
      UInt_t idx  = (*itGet).second;
      Float_t offset;
      Float_t scale;
      if (!fNoOffset) {
         offset = min;
         scale  = 1.0/(max-min);
      } else {
         offset = 0.;
         fabs(max)>fabs(min) ? scale=.5/fabs(max) : scale=.5/fabs(min);
      }
      TString str("");
      VariableInfo& varInfo = (type=='v'?fDsi.GetVariableInfo(idx):(type=='t'?fDsi.GetTargetInfo(idx):fDsi.GetSpectatorInfo(idx)));

      if (offset < 0) str = Form( "2*%g*([%s] + %g) - 1", scale, varInfo.GetLabel().Data(), -offset );
      else            str = Form( "2*%g*([%s] - %g) - 1", scale, varInfo.GetLabel().Data(),  offset );
      (*strVec)[iinp] = str;

      ++iinp;
   }

   return strVec;
}

////////////////////////////////////////////////////////////////////////////////
/// write the transformation to the stream

void TMVA::VariableNormalizeTransform::WriteTransformationToStream( std::ostream& o ) const
{
   o << "# min max for all variables for all classes one after the other and as a last entry for all classes together" << std::endl;

   Int_t numC = GetNClasses()+1;
   if (GetNClasses() <= 1 ) numC = 1;

   UInt_t nvars = GetNVariables();
   UInt_t ntgts = GetNTargets();

   for (Int_t icls = 0; icls < numC; icls++ ) {
      o << icls << std::endl;
      for (UInt_t ivar=0; ivar<nvars; ivar++)
         o << std::setprecision(12) << std::setw(20) << fMin.at(icls).at(ivar) << " "
           << std::setprecision(12) << std::setw(20) << fMax.at(icls).at(ivar) << std::endl;
      for (UInt_t itgt=0; itgt<ntgts; itgt++)
         o << std::setprecision(12) << std::setw(20) << fMin.at(icls).at(nvars+itgt) << " "
           << std::setprecision(12) << std::setw(20) << fMax.at(icls).at(nvars+itgt) << std::endl;
   }
   o << "##" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// create XML description of Normalize transformation

void TMVA::VariableNormalizeTransform::AttachXMLTo(void* parent)
{
   void* trfxml = gTools().AddChild(parent, "Transform");
   gTools().AddAttr(trfxml, "Name", "Normalize");
   gTools().AddAttr(trfxml, "UseOffsetOrNot", (fNoOffset?"NoOffset":"UseOffset") );
   VariableTransformBase::AttachXMLTo( trfxml );

   Int_t numC = (GetNClasses()<= 1)?1:GetNClasses()+1;

   for( Int_t icls=0; icls<numC; icls++ ) {
      void* clsxml = gTools().AddChild(trfxml, "Class");
      gTools().AddAttr(clsxml, "ClassIndex", icls);
      void* inpxml = gTools().AddChild(clsxml, "Ranges");
      UInt_t iinp = 0;
      for( ItVarTypeIdx itGet = fGet.begin(), itGetEnd = fGet.end(); itGet != itGetEnd; ++itGet ) {
         void* mmxml = gTools().AddChild(inpxml, "Range");
         gTools().AddAttr(mmxml, "Index", iinp);
         gTools().AddAttr(mmxml, "Min", fMin.at(icls).at(iinp) );
         gTools().AddAttr(mmxml, "Max", fMax.at(icls).at(iinp) );
         ++iinp;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read the transformation matrices from the xml node

void TMVA::VariableNormalizeTransform::ReadFromXML( void* trfnode )
{
   TString UseOffsetOrNot;

   gTools().ReadAttr(trfnode, "UseOffsetOrNot", UseOffsetOrNot );

   if (UseOffsetOrNot == "NoOffset") fNoOffset = kTRUE;
   else                              fNoOffset = kFALSE;

   Bool_t newFormat = kFALSE;

   void* inpnode = NULL;

   inpnode = gTools().GetChild(trfnode, "Selection"); // new xml format
   if( inpnode != NULL )
      newFormat = kTRUE;

   if( newFormat ){
      // ------------- new format --------------------
      // read input
      VariableTransformBase::ReadFromXML( inpnode );

      // read transformation information

      UInt_t size = fGet.size();
      UInt_t classindex, idx;

      void* ch = gTools().GetChild( trfnode, "Class" );
      while(ch) {
         Int_t ci = 0;
         gTools().ReadAttr(ch, "ClassIndex", ci);
         classindex = UInt_t(ci);

         fMin.resize(classindex+1);
         fMax.resize(classindex+1);

         fMin[classindex].resize(size,Float_t(0));
         fMax[classindex].resize(size,Float_t(0));

         void* clch = gTools().GetChild( ch );
         while(clch) {
            TString nodeName(gTools().GetName(clch));
            if(nodeName=="Ranges") {
               void* varch = gTools().GetChild( clch );
               while(varch) {
                  gTools().ReadAttr(varch, "Index", idx);
                  gTools().ReadAttr(varch, "Min",      fMin[classindex][idx]);
                  gTools().ReadAttr(varch, "Max",      fMax[classindex][idx]);
                  varch = gTools().GetNextChild( varch );
               }
            }
            clch = gTools().GetNextChild( clch );
         }
         ch = gTools().GetNextChild( ch );
      }
      SetCreated();
      return;
   }

   // ------------- old format --------------------
   UInt_t classindex, varindex, tgtindex, nvars, ntgts;
   // coverity[tainted_data_argument]
   gTools().ReadAttr(trfnode, "NVariables", nvars);
   // coverity[tainted_data_argument]
   gTools().ReadAttr(trfnode, "NTargets",   ntgts);
   // coverity[tainted_data_argument]

   for( UInt_t ivar = 0; ivar < nvars; ++ivar ){
      fGet.push_back(std::pair<Char_t,UInt_t>('v',ivar));
   }
   for( UInt_t itgt = 0; itgt < ntgts; ++itgt ){
      fGet.push_back(std::pair<Char_t,UInt_t>('t',itgt));
   }
   void* ch = gTools().GetChild( trfnode );
   while(ch) {
      gTools().ReadAttr(ch, "ClassIndex", classindex);

      fMin.resize(classindex+1);
      fMax.resize(classindex+1);
      fMin[classindex].resize(nvars+ntgts,Float_t(0));
      fMax[classindex].resize(nvars+ntgts,Float_t(0));

      void* clch = gTools().GetChild( ch );
      while(clch) {
         TString nodeName(gTools().GetName(clch));
         if(nodeName=="Variables") {
            void* varch = gTools().GetChild( clch );
            while(varch) {
               gTools().ReadAttr(varch, "VarIndex", varindex);
               gTools().ReadAttr(varch, "Min",      fMin[classindex][varindex]);
               gTools().ReadAttr(varch, "Max",      fMax[classindex][varindex]);
               varch = gTools().GetNextChild( varch );
            }
         } else if (nodeName=="Targets") {
            void* tgtch = gTools().GetChild( clch );
            while(tgtch) {
               gTools().ReadAttr(tgtch, "TargetIndex", tgtindex);
               gTools().ReadAttr(tgtch, "Min",      fMin[classindex][nvars+tgtindex]);
               gTools().ReadAttr(tgtch, "Max",      fMax[classindex][nvars+tgtindex]);
               tgtch = gTools().GetNextChild( tgtch );
            }
         }
         clch = gTools().GetNextChild( clch );
      }
      ch = gTools().GetNextChild( ch );
   }
   SetCreated();
}

////////////////////////////////////////////////////////////////////////////////
/// this method is only used when building a normalization transformation
/// from old text files
/// in this case regression didn't exist and there were no targets

void TMVA::VariableNormalizeTransform::BuildTransformationFromVarInfo( const std::vector<TMVA::VariableInfo>& var )
{
   UInt_t nvars = GetNVariables();

   if(var.size() != nvars)
      Log() << kFATAL << "<BuildTransformationFromVarInfo> can't build transformation,"
            << " since the number of variables disagree" << Endl;

   UInt_t numC = (GetNClasses()<=1)?1:GetNClasses()+1;
   fMin.clear();fMin.resize( numC );
   fMax.clear();fMax.resize( numC );


   for(UInt_t cls=0; cls<numC; ++cls) {
      fMin[cls].resize(nvars+GetNTargets(),0);
      fMax[cls].resize(nvars+GetNTargets(),0);
      UInt_t vidx(0);
      for(std::vector<TMVA::VariableInfo>::const_iterator v = var.begin(); v!=var.end(); ++v, ++vidx) {
         fMin[cls][vidx] = v->GetMin();
         fMax[cls][vidx] = v->GetMax();
         fGet.push_back(std::pair<Char_t,UInt_t>('v',vidx));
      }
   }
   SetCreated();
}

////////////////////////////////////////////////////////////////////////////////
/// Read the variable ranges from an input stream

void TMVA::VariableNormalizeTransform::ReadTransformationFromStream( std::istream& istr, const TString& )
{
   UInt_t nvars = GetNVariables();
   UInt_t ntgts = GetNTargets();
   for( UInt_t ivar = 0; ivar < nvars; ++ivar ){
      fGet.push_back(std::pair<Char_t,UInt_t>('v',ivar));
   }
   for( UInt_t itgt = 0; itgt < ntgts; ++itgt ){
      fGet.push_back(std::pair<Char_t,UInt_t>('t',itgt));
   }
   char buf[512];
   char buf2[512];
   istr.getline(buf,512);
   TString strvar, dummy;
   Int_t icls;
   TString test;
   while (!(buf[0]=='#'&& buf[1]=='#')) { // if line starts with ## return
      char* p = buf;
      while (*p==' ' || *p=='\t') p++; // 'remove' leading whitespace
      if (*p=='#' || *p=='\0') {
         istr.getline(buf,512);
         continue; // if comment or empty line, read the next line
      }
      std::stringstream sstr(buf);
      sstr >> icls;
      for (UInt_t ivar=0;ivar<nvars;ivar++) {
         istr.getline(buf2,512); // reading the next line
         std::stringstream sstr2(buf2);
         sstr2 >> fMin[icls][ivar] >> fMax[icls][ivar];
      }
      for (UInt_t itgt=0;itgt<ntgts;itgt++) {
         istr.getline(buf2,512); // reading the next line
         std::stringstream sstr2(buf2);
         sstr2 >> fMin[icls][nvars+itgt] >> fMax[icls][nvars+itgt];
      }
      istr.getline(buf,512); // reading the next line
   }
   SetCreated();
}

////////////////////////////////////////////////////////////////////////////////
/// prints the transformation ranges

void TMVA::VariableNormalizeTransform::PrintTransformation( std::ostream& /* o */ )
{
   Int_t nCls = GetNClasses();
   Int_t numC = nCls+1;
   if (nCls <= 1 ) numC = 1;
   for (Int_t icls = 0; icls < numC; icls++ ) {
      if( icls == nCls )
         Log() << kINFO << "Transformation for all classes based on these ranges:" << Endl;
      else
         Log() << kINFO << "Transformation for class " << icls << " based on these ranges:" << Endl;
      UInt_t iinp = 0;
      for( ItVarTypeIdxConst itGet = fGet.begin(), itGetEnd = fGet.end(); itGet != itGetEnd; ++itGet ){
         Char_t type = (*itGet).first;
         UInt_t idx  = (*itGet).second;

         TString typeString = (type=='v'?"Variable: ": (type=='t'?"Target : ":"Spectator : ") );
         Log() << typeString.Data() << std::setw(20) << fMin[icls][idx] << std::setw(20) << fMax[icls][idx] << Endl;

         ++iinp;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// creates a normalizing function
/// TODO include target-transformation into makefunction

void TMVA::VariableNormalizeTransform::MakeFunction( std::ostream& fout, const TString& fcncName,
                                                     Int_t part, UInt_t trCounter, Int_t )
{
   UInt_t nVar = fGet.size();
   UInt_t numC = fMin.size();
   if (part == 1) {
      fout << std::endl;
      fout << "   double fOff_" << trCounter << "[" << numC << "][" << nVar << "];" << std::endl;
      fout << "   double fScal_" << trCounter << "[" << numC << "][" << nVar << "];" << std::endl;
   }

   if (part == 2) {
      fout << std::endl;
      fout << "//_______________________________________________________________________" << std::endl;
      fout << "inline void " << fcncName << "::InitTransform_" << trCounter << "()" << std::endl;
      fout << "{" << std::endl;
      fout << "   double fMin_" << trCounter << "[" << numC << "][" << nVar << "];" << std::endl;
      fout << "   double fMax_" << trCounter << "[" << numC << "][" << nVar << "];" << std::endl;
      fout << "   // Normalization transformation, initialisation" << std::endl;
      for (UInt_t ivar = 0; ivar < nVar; ivar++) {
         for (UInt_t icls = 0; icls < numC; icls++) {
            Double_t min = TMath::Min(FLT_MAX, fMin.at(icls).at(ivar));
            Double_t max = TMath::Max(-FLT_MAX, fMax.at(icls).at(ivar));
            fout << "   fMin_" << trCounter << "[" << icls << "][" << ivar << "] = " << std::setprecision(12) << min
                 << ";" << std::endl;
            fout << "   fMax_" << trCounter << "[" << icls << "][" << ivar << "] = " << std::setprecision(12) << max
                 << ";" << std::endl;
            fout << "   fScal_" << trCounter << "[" << icls << "][" << ivar << "] = 2.0/(fMax_" << trCounter << "["
                 << icls << "][" << ivar << "]-fMin_" << trCounter << "[" << icls << "][" << ivar << "]);" << std::endl;
            fout << "   fOff_" << trCounter << "[" << icls << "][" << ivar << "] = fMin_" << trCounter << "[" << icls
                 << "][" << ivar << "]*fScal_" << trCounter << "[" << icls << "][" << ivar << "]+1.;" << std::endl;
         }
      }
      fout << "}" << std::endl;
      fout << std::endl;
      fout << "//_______________________________________________________________________" << std::endl;
      fout << "inline void " << fcncName << "::Transform_" << trCounter << "( std::vector<double>& iv, int cls) const"
           << std::endl;
      fout << "{" << std::endl;
      fout << "   // Normalization transformation" << std::endl;
      fout << "   if (cls < 0 || cls > " << GetNClasses() << ") {" << std::endl;
      fout << "   if (" << GetNClasses() << " > 1 ) cls = " << GetNClasses() << ";" << std::endl;
      fout << "      else cls = " << (fMin.size() == 1 ? 0 : 2) << ";" << std::endl;
      fout << "   }" << std::endl;
      fout << "   const int nVar = " << nVar << ";" << std::endl << std::endl;
      fout << "   // get indices of used variables" << std::endl;
      VariableTransformBase::MakeFunction(fout, fcncName, 0, trCounter, 0);
      fout << "   static std::vector<double> dv;"
           << std::endl; // simply made it static so it doesn't need to be re-booked every time
      fout << "   dv.resize(nVar);" << std::endl;
      fout << "   for (int ivar=0; ivar<nVar; ivar++) dv[ivar] = iv[indicesGet.at(ivar)];" << std::endl;

      fout << "   for (int ivar=0;ivar<" << nVar << ";ivar++) {" << std::endl;
      fout << "      double offset = fOff_" << trCounter << "[cls][ivar];" << std::endl;
      fout << "      double scale  = fScal_" << trCounter << "[cls][ivar];" << std::endl;
      fout << "      iv[indicesPut.at(ivar)] = scale*dv[ivar]-offset;" << std::endl;
      fout << "   }" << std::endl;
      fout << "}" << std::endl;
   }
}
