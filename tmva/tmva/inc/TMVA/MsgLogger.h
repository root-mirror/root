// @(#)root/tmva $Id$
// Author: Attila Krasznahorkay, Andreas Hoecker, Joerg Stelzer, Eckhard von Toerne

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MsgLogger                                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      TMVA output logger class producing nicely formatted log messages          *
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

#ifndef ROOT_TMVA_MsgLogger
#define ROOT_TMVA_MsgLogger

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MsgLogger                                                            //
//                                                                      //
// ostringstream derivative to redirect and format output               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

// STL include(s):
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#if __cplusplus > 199711L
#include <atomic>
#endif

// ROOT include(s)
#include "TObject.h"

#include "TMVA/Types.h"

// Local include(s):

namespace TMVA {

   class MsgLogger : public std::ostringstream, public TObject {

   public:

      MsgLogger( const TObject* source, EMsgType minType = kINFO );
      MsgLogger( const std::string& source, EMsgType minType = kINFO );
      MsgLogger( EMsgType minType = kINFO );
      MsgLogger( const MsgLogger& parent );
      ~MsgLogger();

      // Accessors
      void        SetSource ( const std::string& source ) { fStrSource = source; }
      EMsgType    GetMinType()                      const { return fMinType; }
      void        SetMinType( EMsgType minType )          { fMinType = minType; }
      std::string GetSource()          const              { return fStrSource; }
      std::string GetPrintedSource()   const;
      std::string GetFormattedSource() const;

      static UInt_t GetMaxSourceSize();

      // Needed for copying
      MsgLogger& operator= ( const MsgLogger& parent );

      // Stream modifier(s)
      static MsgLogger& Endmsg( MsgLogger& logger );

      // Accept stream modifiers
      MsgLogger& operator<< ( MsgLogger& ( *_f )( MsgLogger& ) );
      MsgLogger& operator<< ( std::ostream& ( *_f )( std::ostream& ) );
      MsgLogger& operator<< ( std::ios& ( *_f )( std::ios& ) );

      // Accept message type specification
      MsgLogger& operator<< ( EMsgType type );

      // For all the "conventional" inputs
      template <class T> MsgLogger& operator<< ( T arg ) {
         *(std::ostringstream*)this << arg;
         return *this;
      }

      // Temporally disables all the loggers (Caution! Use with care !)
      static void  InhibitOutput();
      static void  EnableOutput();

   private:

      // private utility routines
      void Send();
      void InitMaps();
      void WriteMsg( EMsgType type, const std::string& line ) const;

      const TObject*           fObjSource;        // the source TObject (used for name)
      std::string              fStrSource;        // alternative string source
      static const std::string fgPrefix;          // the prefix of the source name
      static const std::string fgSuffix;          // suffix following source name
      EMsgType                 fActiveType;       // active type
      static const UInt_t      fgMaxSourceSize;   // maximum length of source name
#if __cplusplus > 199711L
      static std::atomic<Bool_t> fgOutputSupressed; // disable the output globally (used by generic booster)
      static std::atomic<Bool_t> fgInhibitOutput;   // flag to suppress all output

      static std::atomic<const std::map<EMsgType, std::string>*> fgTypeMap;   // matches output types with strings
      static std::atomic<const std::map<EMsgType, std::string>*> fgColorMap;  // matches output types with terminal colors
#else
      static Bool_t            fgOutputSupressed; // disable the output globally (used by generic booster)
      static Bool_t            fgInhibitOutput;   // flag to suppress all output

      static const std::map<EMsgType, std::string>* fgTypeMap;   // matches output types with strings
      static const std::map<EMsgType, std::string>* fgColorMap;  // matches output types with terminal colors
#endif
      EMsgType                                fMinType;    // minimum type for output

      ClassDef(MsgLogger,0) // Ostringstream derivative to redirect and format logging output
   }; // class MsgLogger

   inline MsgLogger& MsgLogger::operator<< ( MsgLogger& (*_f)( MsgLogger& ) )
   {
      return (_f)(*this);
   }

   inline MsgLogger& MsgLogger::operator<< ( std::ostream& (*_f)( std::ostream& ) )
   {
      (_f)(*this);
      return *this;
   }

   inline MsgLogger& MsgLogger::operator<< ( std::ios& ( *_f )( std::ios& ) )
   {
      (_f)(*this);
      return *this;
   }

   inline MsgLogger& MsgLogger::operator<< ( EMsgType type )
   {
      fActiveType = type;
      return *this;
   }

   // Shortcut
   inline MsgLogger& Endl(MsgLogger& ml) { return MsgLogger::Endmsg(ml); }

}

#endif // TMVA_MsgLogger
