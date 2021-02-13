// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Configurable                                                          *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
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
 *                                                                                *
 **********************************************************************************/

/*! \class TMVA::
\ingroup TMVA

   Base Class for all classes that need option parsing

*/

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

#include "TSystem.h"
#include "TString.h"
#include "TObjString.h"
#include "TQObject.h"
#include "TSpline.h"
#include "TMatrix.h"

#include "TMVA/Configurable.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

// don't change this flag without a good reason ! The FitterBase code won't work anymore !!!
// #define TMVA_Configurable_SanctionUnknownOption kTRUE

ClassImp(TMVA::Configurable);

#ifdef _WIN32
/*Disable warning C4355: 'this' : used in base member initializer list*/
#pragma warning ( disable : 4355 )
#endif

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::Configurable::Configurable( const TString& theOption)
: TNamed("Configurable","Configurable"),
   fOptions                    ( theOption ),
   fLooseOptionCheckingEnabled ( kTRUE ),
   fLastDeclaredOption         ( 0 ),
   fConfigDescription          ( "No description" ),
   fReferenceFile              ( "None" ),
   fLogger                     ( new MsgLogger(this) )
{
   fListOfOptions.SetOwner();

   // check if verbosity "V" set in option
   if (gTools().CheckForVerboseOption( theOption )) Log().SetMinType( kVERBOSE );
}

////////////////////////////////////////////////////////////////////////////////
/// default destructor

TMVA::Configurable::~Configurable()
{
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// splits the option string at ':' and fills the list 'loo' with the primitive strings

void TMVA::Configurable::SplitOptions(const TString& theOpt, TList& loo) const
{
   TString splitOpt(theOpt);
   loo.SetOwner();
   while (splitOpt.Length()>0) {
      if (!splitOpt.Contains(':')) {
         loo.Add(new TObjString(splitOpt));
         splitOpt = "";
      }
      else {
         TString toSave = splitOpt(0,splitOpt.First(':'));
         loo.Add(new TObjString(toSave.Data()));
         splitOpt = splitOpt(splitOpt.First(':')+1,splitOpt.Length());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// resets the IsSet flag for all declare options
/// to be called before options are read from stream

void TMVA::Configurable::ResetSetFlag()
{
   TListIter decOptIt(&fListOfOptions); // declared options
   while (OptionBase* decOpt = (OptionBase*) decOptIt()) { // loop over declared options
      decOpt->fIsSet = kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// options parser

void TMVA::Configurable::ParseOptions()
{
   Log() << kVERBOSE << "Parsing option string: " << Endl;
   TString optionsWithoutTilde(fOptions);
   optionsWithoutTilde.ReplaceAll(TString("~"),TString(""));
   Log() << kVERBOSE << "... \"" << optionsWithoutTilde << "\"" << Endl;

   TList loo; // the List Of Options in the parsed string

   fOptions = fOptions.Strip(TString::kLeading, ':');

   // separate the options by the ':' marker
   SplitOptions(fOptions, loo);
   fOptions = "";

   // loop over the declared options and check for their availability
   std::map<TString, std::vector<std::pair<Int_t, TString> > > arrayTypeOptions;

   TListIter decOptIt(&fListOfOptions); // declared options
   TListIter setOptIt(&loo);   // parsed options
   while (TObjString * os = (TObjString*) setOptIt()) { // loop over parsed options

      TString s = os->GetString();

      // the tilde in the beginning is an indication that the option
      // has been accepted during previous parsing
      //
      // while parsing this option string eventual appearances of the
      // tilde will be preserved, for correctly parsed options a new
      // one will be added (in the end it will be checked if all
      // options were parsed
      Bool_t preserveTilde = s.BeginsWith('~');
      s = s.Strip(TString::kLeading, '~');

      Bool_t paramParsed = kFALSE;
      if (s.Contains('=')) { // desired way of setting an option: "...:optname=optvalue:..."
         TString optname = s(0,s.First('=')); optname.ToLower();
         TString optval = s(s.First('=')+1,s.Length());
         Int_t idx = -1;

         // First check if the optname exists in the list of the
         // objects. This does not depend on the existence of a [] in
         // the optname. Sometimes the [] is part of the optname and
         // does not describe an array
         OptionBase* decOpt = (OptionBase *)fListOfOptions.FindObject(optname);
         if (decOpt==0 && optname.Contains('[')) {
            // now we see if there is an [] and if the optname exists
            // after removing the [idx]
            TString st = optname(optname.First('[')+1,100);
            st.Remove(st.First(']'));
            std::stringstream str(st.Data());
            str >> idx;                              // save the array index
            optname.Remove(optname.First('['));      // and remove [idx] from the option name
            decOpt = (OptionBase *)fListOfOptions.FindObject(optname);
         }

         TListIter optIt(&fListOfOptions);
         if (decOpt!=0) {
            if (decOpt->IsSet())
               Log() << kWARNING << "Value for option " << decOpt->GetName()
                     << " was previously set to " << decOpt->GetValue() << Endl;

            if (!decOpt->HasPreDefinedVal() || (decOpt->HasPreDefinedVal() && decOpt->IsPreDefinedVal(optval)) ) {
               if (decOpt->IsArrayOpt()) { // arrays
                  // if no index was found then we assume the value is to be set for the entire array
                  if (idx==-1) {
                     decOpt->SetValue(optval);
                  }
                  else {
                     // since we don't know what else is coming we just put everthing into a map
                     if (!decOpt->SetValue(optval, idx))
                        Log() << kFATAL << "Index " << idx << " too large for option " << decOpt->TheName()
                              << ", allowed range is [0," << decOpt->GetArraySize()-1 << "]" << Endl;
                  }
               }
               else { // no arrays
                  if (idx!=-1)
                     Log() << kFATAL << "Option " << decOpt->TheName()
                           << " is not an array, but you specified an index" << Endl;
                  decOpt->SetValue(optval);
               }
               paramParsed = kTRUE;
            }
            else Log() << kFATAL << "Option " << decOpt->TheName()
                       << " does not have predefined value: \"" << optval << "\"" << Endl;
         }
      }

      // boolean variables can be specified by just their name (!name),
      // which will set the to true (false):  ...:V:...:!S:..
      Bool_t preserveNotSign = kFALSE;
      if (!paramParsed) {
         Bool_t hasNotSign = kFALSE;
         if (s.BeginsWith("!")) { s.Remove(0,1); preserveNotSign = hasNotSign = kTRUE; }
         TString optname(s); optname.ToLower();
         OptionBase* decOpt = 0;
         Bool_t optionExists = kFALSE;
         TListIter optIt(&fListOfOptions);
         while ((decOpt = (OptionBase*)optIt()) !=0) {
            TString predOptName(decOpt->GetName());
            predOptName.ToLower();
            if (predOptName == optname) optionExists = kTRUE;
            if (dynamic_cast<Option<bool>*>(decOpt)==0) continue; // not a boolean option
            if (predOptName == optname) break;
         }


         if (decOpt != 0) {
            decOpt->SetValue( hasNotSign ? "0" : "1" );
            paramParsed = kTRUE;
         }
         else {
            if (optionExists && hasNotSign) {
               Log() << kFATAL << "Negating a non-boolean variable " << optname
                     << ", please check the options for method: " << GetName() << Endl;
            }
         }
      }


      if (!paramParsed && LooseOptionCheckingEnabled()) {
         // loose options specification, loops through the possible string
         // values any parameter can have not applicable for boolean or floats
         decOptIt.Reset();
         while (OptionBase* decOpt = (OptionBase*) decOptIt()) {
            if (decOpt->HasPreDefinedVal() && decOpt->IsPreDefinedVal(s) ) {
               paramParsed = decOpt->SetValue(s);
               break;
            }
         }
      }

      if (fOptions!="") fOptions += ":";
      if (paramParsed || preserveTilde) fOptions += '~';
      if (preserveNotSign) fOptions += '!';
      fOptions += s;
   }

   // print options summary
   PrintOptions();
   if (gConfig().WriteOptionsReference()) WriteOptionsReferenceToFile();
}

////////////////////////////////////////////////////////////////////////////////
/// checks for unused options in option string

void TMVA::Configurable::CheckForUnusedOptions() const
{
   TString theOpt(fOptions);
   theOpt = theOpt.Strip(TString::kLeading, ':');

   // separate the options by the ':' marker
   TList loo; // the List of Options in the parsed string
   SplitOptions(theOpt, loo);

   TListIter setOptIt(&loo);   // options in a list
   TString unusedOptions("");
   while (TObjString * os = (TObjString*) setOptIt()) { // loop over parsed options

      TString s = os->GetString();
      if (!s.BeginsWith('~')) {
         if (unusedOptions != "") unusedOptions += ':';
         unusedOptions += s;
      }
   }
   if (unusedOptions != "") {
      Log() << kFATAL
            << "The following options were specified, but could not be interpreted: \'"
            << unusedOptions << "\', please check!" << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// prints out the options set in the options string and the defaults

void TMVA::Configurable::PrintOptions() const
{
   Log() << kVERBOSE << "The following options are set:" << Endl;

   TListIter optIt( &fListOfOptions );
   Log() << kVERBOSE << "- By User:" << Endl;
   Bool_t found = kFALSE;
   while (OptionBase* opt = (OptionBase *) optIt()) {
      if (opt->IsSet()) {
          Log() << kVERBOSE << "    ";
          std::ostringstream oss;
          opt->Print(oss);
          Log() << oss.str();
          Log() << Endl; found = kTRUE; }
   }
   if (!found) Log() << kVERBOSE << "    <none>" << Endl;

   optIt.Reset();
   Log() << kVERBOSE << "- Default:" << Endl;
   found = kFALSE;
   while (OptionBase* opt = (OptionBase *) optIt()) {
      if (!opt->IsSet()) {
          Log() << kVERBOSE << "    ";
          std::ostringstream oss;
          opt->Print(oss);
          Log() << oss.str();
          Log() << Endl; found = kTRUE; }
   }
   if (!found) Log() << kVERBOSE << "    <none>" << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// write options to output stream (e.g. in writing the MVA weight files

void TMVA::Configurable::WriteOptionsToStream( std::ostream& o, const TString& prefix ) const
{
   TListIter optIt( &fListOfOptions );
   o << prefix << "# Set by User:" << std::endl;
   while (OptionBase * opt = (OptionBase *) optIt())
      if (opt->IsSet()) { o << prefix; opt->Print(o); o << std::endl; }
   optIt.Reset();
   o << prefix << "# Default:" << std::endl;
   while (OptionBase * opt = (OptionBase *) optIt())
      if (!opt->IsSet()) { o << prefix; opt->Print(o); o << std::endl; }
   o << prefix << "##" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// write options to XML file

void TMVA::Configurable::AddOptionsXMLTo( void* parent ) const
{
   if (!parent) return;
   void* opts = gTools().AddChild(parent, "Options");
   TListIter optIt( &fListOfOptions );
   while (OptionBase * opt = (OptionBase *) optIt()) {
      void* optnode = 0;
      if (opt->IsArrayOpt()) {
         std::stringstream s("");
         s.precision( 16 );
         for(Int_t i=0; i<opt->GetArraySize(); i++) {
            if(i>0) s << " ";
            s << std::scientific << opt->GetValue(i);
         }
         optnode = gTools().AddChild(opts,"Option",s.str().c_str());
      }
      else {
         optnode = gTools().AddChild(opts,"Option", opt->GetValue());
      }
      gTools().AddAttr(optnode, "name", opt->TheName());
      if (opt->IsArrayOpt()) {
         gTools().AddAttr(optnode, "size", opt->GetArraySize());
      }
      gTools().AddAttr(optnode, "modified", (opt->IsSet()?"Yes":"No") );
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::Configurable::ReadOptionsFromXML( void* node )
{
   void* opt = gTools().GetChild(node);
   TString optName, optValue;
   fOptions="";
   while (opt != 0) {
      if (fOptions.Length()!=0) fOptions += ":";
      gTools().ReadAttr(opt, "name", optName);
      optValue = TString( gTools().GetContent(opt) );
      std::stringstream s("");
      s.precision( 16 );
      if (gTools().HasAttr(opt, "size")) {
         UInt_t size;
         gTools().ReadAttr(opt, "size", size);
         std::vector<TString> values = gTools().SplitString(optValue, ' ');
         for(UInt_t i=0; i<size; i++) {
            if(i!=0) s << ":";
            s << std::scientific << optName << "[" << i << "]=" << values[i];
         }
      }
      else {
         s << std::scientific << optName << "=" << optValue;
      }
      fOptions += s.str().c_str();
      opt = gTools().GetNextChild(opt);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// write complete options to output stream

void TMVA::Configurable::WriteOptionsReferenceToFile()
{
   TString dir = gConfig().GetIONames().fOptionsReferenceFileDir;
   gSystem->MakeDirectory( dir );
   fReferenceFile = dir + "/" + GetConfigName() + "_optionsRef.txt";
   std::ofstream o( fReferenceFile );
   if (!o.good()) { // file could not be opened --> Error
      Log() << kFATAL << "<WriteOptionsToInfoFile> Unable to open output file: " << fReferenceFile << Endl;
   }

   TListIter optIt( &fListOfOptions );
   o << "# List of options:" << std::endl;
   o << "# Configurable: " << GetConfigName() << std::endl;
   o << "# Description: " << GetConfigDescription() << std::endl;
   while (OptionBase * opt = (OptionBase *) optIt()) {
      opt->Print( o, 1 );
      o << std::endl << "# ------------------------------------------------" << std::endl;
   }

   o.close();
   Log() << kVERBOSE << "Wrote options reference file: \"" << fReferenceFile << "\"" << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// read option back from the weight file

void TMVA::Configurable::ReadOptionsFromStream(std::istream& istr)
{
   // first set the IsSet flag of all declared options to false
   // that is only necessary in our factory, when we test right
   // after the training
   ResetSetFlag();
   fOptions = "";
   char buf[512];
   istr.getline(buf,512);
   TString stropt, strval;
   while (istr.good() && !istr.eof() && !(buf[0]=='#' && buf[1]=='#')) { // if line starts with ## return
      char *p = buf;
      while (*p==' ' || *p=='\t') p++; // 'remove' leading whitespace
      if (*p=='#' || *p=='\0') {
         istr.getline(buf,512); // reading the next line
         continue; // if comment or empty line, read the next line
      }
      std::stringstream sstr(buf);
      sstr >> stropt >> strval;
      stropt.ReplaceAll(':','=');
      strval.ReplaceAll("\"","");
      if (fOptions.Length()!=0) fOptions += ":";
      fOptions += stropt;
      fOptions += strval;
      istr.getline(buf,512); // reading the next line
   }
}

