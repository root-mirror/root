// @(#)root/xmlparser:$Id$
// Author: Jose Lo   12/4/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TXMLParser
#define ROOT_TXMLParser

#include "TQObject.h"

#include "TObject.h"

#include "TString.h"

struct   _xmlParserCtxt;


class TXMLParser : public TObject, public TQObject {

private:
   TXMLParser(const TXMLParser&) = delete;
   TXMLParser& operator=(const TXMLParser&) = delete;

protected:
   _xmlParserCtxt     *fContext;          ///< Parse the xml file
   Bool_t              fValidate;         ///< To validate the parse context
   Bool_t              fReplaceEntities;  ///< Replace entities
   Bool_t              fStopError;        ///< Stop when parse error occurs
   TString             fValidateError;    ///< Parse error
   TString             fValidateWarning;  ///< Parse warning
   Int_t               fParseCode;        ///< To keep track of the errorcodes

   virtual void        InitializeContext();
   virtual void        ReleaseUnderlying();
   virtual void        OnValidateError(const TString& message);
   virtual void        OnValidateWarning(const TString& message);
   virtual void        SetParseCode(Int_t code);

public:
   TXMLParser();
   virtual ~TXMLParser();

   void                SetValidate(Bool_t val = kTRUE);
   Bool_t              GetValidate() const { return fValidate; }

   void                SetReplaceEntities(Bool_t val = kTRUE);
   Bool_t              GetReplaceEntities() const { return fReplaceEntities; }

   virtual Int_t       ParseFile(const char *filename) = 0;
   virtual Int_t       ParseBuffer(const char *contents, Int_t len) = 0;
   virtual void        StopParser();

   Int_t               GetParseCode() const { return fParseCode; }

   const char         *GetParseCodeMessage(Int_t parseCode) const;

   void                SetStopOnError(Bool_t stop = kTRUE);
   Bool_t              GetStopOnError() const { return fStopError; }

   const char         *GetValidateError() const { return fValidateError; }
   const char         *GetValidateWarning() const { return fValidateWarning; }

   ClassDef(TXMLParser,0);  // XML SAX parser
};

#endif
