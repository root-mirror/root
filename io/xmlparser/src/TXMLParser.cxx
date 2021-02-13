// @(#)root/xmlparser:$Id$
// Author: Jose Lo   12/1/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
\class TXMLParser
\ingroup IO

TXMLParser is an abstract class which interfaces with Libxml2.
Libxml2 is the XML C parser and toolkit developed for the Gnome
project.
The libxml library provides two interfaces to the parser, a DOM
style tree interface and a SAX style event based interface.
TXMLParser is parent class of TSAXParser and TDOMParser, which are
a SAX interface and DOM interface of libxml.
*/

/*************************************************************************
  This source is based on libxml++, a C++ wrapper for the libxml XML
  parser library. Copyright (C) 2000 by Ari Johnson.

  libxml++ are copyright (C) 2000 by Ari Johnson, and are covered by the
  GNU Lesser General Public License, which should be included with
  libxml++ as the file COPYING.
 *************************************************************************/

#include "TXMLParser.h"

#include <libxml/parser.h>


namespace {
   // See https://lists.fedoraproject.org/pipermail/devel/2010-January/129117.html :
   // "That function might delete TLS fields that belong to other libraries
   // [...] if called twice."
   // The same (though with less dramatic consequences) holds for xmlInitParser().
   struct InitAndCleanupTheXMLParserOnlyOnceCommaEver {
      InitAndCleanupTheXMLParserOnlyOnceCommaEver() {
         xmlInitParser();
      }
      ~InitAndCleanupTheXMLParserOnlyOnceCommaEver() {
         xmlCleanupParser();
      }
   } gInitAndCleanupTheXMLParserOnlyOnceCommaEver;
}

ClassImp(TXMLParser);

////////////////////////////////////////////////////////////////////////////////
/// Initializes parser variables.

TXMLParser::TXMLParser()
   : fContext(0), fValidate(kTRUE), fReplaceEntities(kFALSE), fStopError(kFALSE), fParseCode(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup.

TXMLParser::~TXMLParser()
{
   ReleaseUnderlying();
   fParseCode = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// The parser will validate the xml file if val = true.

void TXMLParser::SetValidate(Bool_t val)
{
   fValidate = val;
}

////////////////////////////////////////////////////////////////////////////////
/// The parser will replace/expand entities.

void TXMLParser::SetReplaceEntities(Bool_t val)
{
   fReplaceEntities = val;
}

////////////////////////////////////////////////////////////////////////////////
/// To release any existing document.

void TXMLParser::ReleaseUnderlying()
{
   if (fContext) {
      fContext->_private = 0;
      xmlFreeParserCtxt(fContext);
      fContext = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This function is called when an error from the parser has occured.
/// Message is the parse error.

void TXMLParser::OnValidateError(const TString& message)
{
   fValidateError += message;
}

////////////////////////////////////////////////////////////////////////////////
/// This function is called when a warning from the parser has occured.
/// Message is the parse error.

void TXMLParser::OnValidateWarning(const TString& message)
{
   fValidateWarning += message;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the parse code message.

const char *TXMLParser::GetParseCodeMessage(Int_t parseCode) const
{
   switch (parseCode) {
      case -1:
         return "Attempt to parse a second file while a parse is in progress";
         break;
      case -2:
         return "Parse context is not created";
         break;
      case -3:
         return "An error occured while parsing file";
         break;
      case -4:
         return "A fatal error occured while parsing file";
         break;
      case -5:
         return "Document is not well-formed";
         break;
      case -6:
         return "Document is not valid";
         break;
      default:
         return "Parse code does not exist";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize parser parameters, such as, disactivate non-standards libxml1
/// features, on/off validation, clear error and warning messages.

void TXMLParser::InitializeContext()
{
   fContext->linenumbers = 1; // TRUE - This is the default anyway.
   fContext->validate = fValidate ? 1 : 0;
   fContext->replaceEntities = fReplaceEntities ? 1 : 0;
   fContext->_private = this;

   fValidateError = "";
   fValidateWarning = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Stops parsing.

void TXMLParser::StopParser()
{
   if (fContext)
      xmlStopParser(fContext);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the parse code:
///   - \b 0: Parse successful
///   - \b -1: Attempt to parse a second file while a parse is in progress
///   - \b -2: Parse context is not created
///   - \b -3: An error occured while parsing file
///   - \b -4: A fatal error occured while parsing file
///   - \b -5: Document is not well-formed

void TXMLParser::SetParseCode(Int_t errorcode)
{
   fParseCode = errorcode;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parser stops in case of error:
///   - \b stop = true, stops on error
///   - \b stop = false, continue parsing on error...

void TXMLParser::SetStopOnError(Bool_t stop)
{
   fStopError = stop;
}
