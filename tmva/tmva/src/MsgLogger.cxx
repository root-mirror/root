// @(#)root/tmva $Id$
// Author: Attila Krasznahorkay, Andreas Hoecker, Joerg Stelzer, Eckhard von Toerne

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
 *      Attila Krasznahorkay <Attila.Krasznahorkay@cern.ch> - CERN, Switzerland   *
 *      Andreas Hoecker       <Andreas.Hocker@cern.ch> - CERN, Switzerland        *
 *      Joerg Stelzer         <stelzer@cern.ch>        - DESY, Germany            *
 *      Eckhard v. Toerne     <evt@uni-bonn.de>        - U of Bonn, Germany       *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MsgLogger
\ingroup TMVA
ostringstream derivative to redirect and format output
*/

// Local include(s):
#include "TMVA/MsgLogger.h"

#include "TMVA/Config.h"
#include "TMVA/Types.h"

// ROOT include(s):
#include "Rtypes.h"
#include "TObject.h"

// STL include(s):
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <memory>


ClassImp(TMVA::MsgLogger);

// declaration of global variables
// this is the hard-coded maximum length of the source names
const UInt_t                           TMVA::MsgLogger::fgMaxSourceSize = 25;

const std::string                      TMVA::MsgLogger::fgPrefix = "";
const std::string                      TMVA::MsgLogger::fgSuffix = ": ";
#if __cplusplus > 199711L
std::atomic<Bool_t>                                       TMVA::MsgLogger::fgInhibitOutput{kFALSE};
std::atomic<const std::map<TMVA::EMsgType, std::string>*> TMVA::MsgLogger::fgTypeMap{0};
std::atomic<const std::map<TMVA::EMsgType, std::string>*> TMVA::MsgLogger::fgColorMap{0};
#else
Bool_t                                       TMVA::MsgLogger::fgInhibitOutput = kFALSE;
const std::map<TMVA::EMsgType, std::string>* TMVA::MsgLogger::fgTypeMap  = 0;
const std::map<TMVA::EMsgType, std::string>* TMVA::MsgLogger::fgColorMap = 0;
#endif
static std::unique_ptr<const std::map<TMVA::EMsgType, std::string> > gOwnTypeMap;
static std::unique_ptr<const std::map<TMVA::EMsgType, std::string> > gOwnColorMap;


void   TMVA::MsgLogger::InhibitOutput() { fgInhibitOutput = kTRUE;  }
void   TMVA::MsgLogger::EnableOutput()  { fgInhibitOutput = kFALSE; }
////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::MsgLogger::MsgLogger( const TObject* source, EMsgType minType )
   : fObjSource ( source ),
     fStrSource ( "" ),
     fActiveType( kINFO ),
     fMinType   ( minType )
{
   InitMaps();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::MsgLogger::MsgLogger( const std::string& source, EMsgType minType )
   : fObjSource ( 0 ),
     fStrSource ( source ),
     fActiveType( kINFO ),
     fMinType   ( minType )
{
   InitMaps();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::MsgLogger::MsgLogger( EMsgType minType )
   : fObjSource ( 0 ),
     fStrSource ( "Unknown" ),
     fActiveType( kINFO ),
     fMinType   ( minType )
{
   InitMaps();
}

////////////////////////////////////////////////////////////////////////////////
/// copy constructor

TMVA::MsgLogger::MsgLogger( const MsgLogger& parent )
   : std::basic_ios<MsgLogger::char_type, MsgLogger::traits_type>(),
     std::ostringstream(),
     TObject(),
     fObjSource(0)
{
   InitMaps();
   *this = parent;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MsgLogger::~MsgLogger()
{
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

TMVA::MsgLogger& TMVA::MsgLogger::operator= ( const MsgLogger& parent )
{
   if (&parent != this) {
      fObjSource  = parent.fObjSource;
      fStrSource  = parent.fStrSource;
      fActiveType = parent.fActiveType;
      fMinType    = parent.fMinType;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// make sure the source name is no longer than fgMaxSourceSize:

std::string TMVA::MsgLogger::GetFormattedSource() const
{
   std::string source_name;
   if (fActiveType == kHEADER)
   {
       source_name = fStrSource;
   }
   if (fActiveType == kWARNING)
     {
       source_name ="<WARNING>";
     }
   if (source_name.size() > fgMaxSourceSize) {
      source_name = source_name.substr( 0, fgMaxSourceSize - 3 );
      source_name += "...";
   }

   return source_name;
}

////////////////////////////////////////////////////////////////////////////////
/// returns the maximum source size

UInt_t TMVA::MsgLogger::GetMaxSourceSize()
{
   return static_cast<UInt_t>(fgMaxSourceSize);
}

////////////////////////////////////////////////////////////////////////////////
/// the full logger prefix

std::string TMVA::MsgLogger::GetPrintedSource() const
{
   std::string source_name = GetFormattedSource();
   if (source_name.size() < fgMaxSourceSize)
      for (std::string::size_type i=source_name.size(); i<fgMaxSourceSize; i++) source_name.push_back( ' ' );

   return fgPrefix + source_name + fgSuffix;
}

////////////////////////////////////////////////////////////////////////////////
/// activates the logger writer

void TMVA::MsgLogger::Send()
{
   // make sure the source name is no longer than fgMaxSourceSize:
   std::string source_name = GetFormattedSource();

   std::string message = this->str();
   std::string::size_type previous_pos = 0, current_pos = 0;

   // slice the message into lines:
   while (kTRUE) {
      current_pos = message.find( '\n', previous_pos );
      std::string line = message.substr( previous_pos, current_pos - previous_pos );

      std::ostringstream message_to_send;
      // must call the modifiers like this, otherwise g++ get's confused with the operators...
      message_to_send.setf( std::ios::adjustfield, std::ios::left );
      message_to_send.width( fgMaxSourceSize );
      message_to_send << source_name << fgSuffix << line;
      std::string msg = message_to_send.str();
      this->WriteMsg( fActiveType, msg );

      if (current_pos == message.npos) break;
      previous_pos = current_pos + 1;
   }

   // reset the stream buffer:
   this->str( "" );
   fActiveType = kINFO; // To always print messages that have no level specified...
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// putting the output string, the message type, and the color
/// switcher together into a single string

void TMVA::MsgLogger::WriteMsg( EMsgType type, const std::string& line ) const
{
  if ( (type < fMinType || fgInhibitOutput) && type!=kFATAL ) return; // no output

  std::map<EMsgType, std::string>::const_iterator stype;

  if ((stype = fgTypeMap.load()->find( type )) != fgTypeMap.load()->end()) {
    if (!gConfig().IsSilent() || type==kFATAL) {
      if (gConfig().UseColor()) {
   // no text for INFO or VERBOSE
   if (type == kHEADER || type ==kWARNING)
     std::cout << fgPrefix << line << std::endl;
   else if (type == kINFO || type == kVERBOSE)
     //std::cout << fgPrefix << line << std::endl; // no color for info
     std::cout << line << std::endl;
   else{
     //std::cout<<"prefix='"<<fgPrefix<<"'"<<std::endl;
     std::cout << fgColorMap.load()->find( type )->second << "<" << stype->second << ">" << line << "\033[0m" << std::endl;
}
      }

      else {
   if (type == kINFO) std::cout << fgPrefix << line << std::endl;
   else               std::cout << fgPrefix << "<" << stype->second << "> " << line << std::endl;
      }
    }
  }

   // take decision to stop if fatal error
   if (type == kFATAL) {
      std::cout << "***> abort program execution" << std::endl;
      throw std::runtime_error("FATAL error");

      //std::exit(1);
      //assert(false);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// end line

TMVA::MsgLogger& TMVA::MsgLogger::Endmsg( MsgLogger& logger )
{
   logger.Send();
   return logger;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the message type and color maps

void TMVA::MsgLogger::InitMaps()
{
   if(!fgTypeMap) {
      std::map<TMVA::EMsgType, std::string>*tmp  = new std::map<TMVA::EMsgType, std::string>();

      (*tmp)[kVERBOSE]  = std::string("VERBOSE");
      (*tmp)[kDEBUG]    = std::string("DEBUG");
      (*tmp)[kINFO]     = std::string("INFO");
      (*tmp)[kWARNING]  = std::string("WARNING");
      (*tmp)[kERROR]    = std::string("ERROR");
      (*tmp)[kFATAL]    = std::string("FATAL");
      (*tmp)[kSILENT]   = std::string("SILENT");
      (*tmp)[kHEADER]   = std::string("HEADER");
      const std::map<TMVA::EMsgType, std::string>* expected=0;
      if(fgTypeMap.compare_exchange_strong(expected,tmp)) {
         //Have the global own this
         gOwnTypeMap.reset(tmp);
      } else {
         //Another thread beat us in creating the instance
         delete tmp;
      }
   }

   if(!fgColorMap) {
      std::map<TMVA::EMsgType, std::string>*tmp  = new std::map<TMVA::EMsgType, std::string>();

      (*tmp)[kVERBOSE] = std::string("");
      (*tmp)[kDEBUG]   = std::string("\033[34m");
      (*tmp)[kINFO]    = std::string("");
      (*tmp)[kWARNING] = std::string("\033[1;31m");
      (*tmp)[kERROR]   = std::string("\033[31m");
      (*tmp)[kFATAL]   = std::string("\033[37;41;1m");
      (*tmp)[kSILENT]  = std::string("\033[30m");

      const std::map<TMVA::EMsgType, std::string>* expected=0;
      if(fgColorMap.compare_exchange_strong(expected,tmp)) {
         //Have the global own this
         gOwnColorMap.reset(tmp);
      } else {
         //Another thread beat us in creating the instance
         delete tmp;
      }
   }
}
