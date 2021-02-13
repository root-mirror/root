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
\class TSAXParser
\ingroup IO

TSAXParser is a subclass of TXMLParser, it is a wraper class to
libxml library.
SAX (Simple API for XML) is an event based interface, which doesn't
maintain the DOM tree in memory, in other words, it's much more
efficient for large document.
TSAXParserCallback contains a number of callback routines to the
parser in a xmlSAXHandler structure. The parser will then parse the
document and call the appropriate callback when certain conditions
occur.
*/

/*************************************************************************
  This source is based on libxml++, a C++ wrapper for the libxml XML
  parser library.Copyright (C) 2000 by Ari Johnson

  libxml++ are copyright (C) 2000 by Ari Johnson, and are covered by the
  GNU Lesser General Public License, which should be included with
  libxml++ as the file COPYING.
 *************************************************************************/

#include "TSAXParser.h"
#include "TXMLAttr.h"
#include "Varargs.h"
#include "strlcpy.h"
#include "TList.h"
#include "TClass.h"

#include <libxml/parser.h>
#include <libxml/parserInternals.h>


class TSAXParserCallback {
public:
   static void StartDocument(void *fParser);
   static void EndDocument(void *fParser);
   static void StartElement(void *fParser, const xmlChar *name, const xmlChar **p);
   static void EndElement(void *fParser, const xmlChar *name);
   static void Characters(void *fParser, const xmlChar *ch, Int_t len);
   static void Comment(void *fParser, const xmlChar *value);
   static void CdataBlock(void *fParser, const xmlChar *value, Int_t len);
   static void Warning(void *fParser, const char *fmt, ...);
   static void Error(void *fParser, const char *fmt, ...);
   static void FatalError(void *fParser, const char *fmt, ...);
};


ClassImp(TSAXParser);

////////////////////////////////////////////////////////////////////////////////
/// Create SAX parser.

TSAXParser::TSAXParser()
{
   fSAXHandler = new xmlSAXHandler;
   memset(fSAXHandler, 0, sizeof(xmlSAXHandler));

   fSAXHandler->startDocument =
                   (startDocumentSAXFunc)TSAXParserCallback::StartDocument;
   fSAXHandler->endDocument   =
                   (endDocumentSAXFunc)TSAXParserCallback::EndDocument;
   fSAXHandler->startElement  =
                   (startElementSAXFunc)TSAXParserCallback::StartElement;
   fSAXHandler->endElement    =
                   (endElementSAXFunc)TSAXParserCallback::EndElement;
   fSAXHandler->characters    =
                   (charactersSAXFunc)TSAXParserCallback::Characters;
   fSAXHandler->comment       =
                   (commentSAXFunc)TSAXParserCallback::Comment;
   fSAXHandler->cdataBlock    =
                   (cdataBlockSAXFunc)TSAXParserCallback::CdataBlock;
   fSAXHandler->warning       =
                   (warningSAXFunc)TSAXParserCallback::Warning;
   fSAXHandler->error         =
                   (errorSAXFunc)TSAXParserCallback::Error;
   fSAXHandler->fatalError    =
                   (fatalErrorSAXFunc)TSAXParserCallback::FatalError;
}

////////////////////////////////////////////////////////////////////////////////
/// TSAXParser desctructor

TSAXParser::~TSAXParser()
{
   ReleaseUnderlying();

   delete fSAXHandler;
}

////////////////////////////////////////////////////////////////////////////////
/// Emit a signal for OnStartDocument.

void TSAXParser::OnStartDocument()
{
   Emit("OnStartDocument()");
}

////////////////////////////////////////////////////////////////////////////////
/// Emit a signal for OnEndDocument.

void TSAXParser::OnEndDocument()
{
   Emit("OnEndDocument()");
}

////////////////////////////////////////////////////////////////////////////////
/// Emit a signal for OnStarElement, where name is the Element's name and
/// attribute is a TList of (TObjString*, TObjString *) TPair's.
/// The TPair's key is the attribute's name and value is the attribute's
/// value.

void TSAXParser::OnStartElement(const char *name, const TList *attributes)
{
   Long_t args[2];
   args[0] = (Long_t)name;
   args[1] = (Long_t)attributes;

   Emit("OnStartElement(const char *, const TList *)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit a signal for OnEndElement, where name is the Element's name.

void TSAXParser::OnEndElement(const char *name)
{
   Emit("OnEndElement(const char *)", name);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit a signal for OnCharacters, where characters are the characters
/// outside of tags.

void TSAXParser::OnCharacters(const char *characters)
{
   Emit("OnCharacters(const char *)", characters);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit a signal for OnComment, where text is the comment.

void TSAXParser::OnComment(const char *text)
{
   Emit("OnComment(const char *)", text);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit a signal for OnWarning, where text is the warning.

void TSAXParser::OnWarning(const char *text)
{
   Emit("OnWarning(const char *)", text);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit a signal for OnError, where text is the error and it returns the
/// Parse Error Code, see TXMLParser.

Int_t TSAXParser::OnError(const char *text)
{
   TString errmsg;
   errmsg.Form("Source line: %d  %s", fContext->input->line, text);

   Emit("OnError(const char *)", errmsg.Data());
   return -3;
}

////////////////////////////////////////////////////////////////////////////////
/// Emit a signal for OnFactalError, where text is the error and it
/// returns the Parse Error Code, see TXMLParser.

Int_t TSAXParser::OnFatalError(const char *text)
{
   TString errmsg;
   errmsg.Form("Source line: %d  %s", fContext->input->line, text);

   Emit("OnFatalError(const char *)", errmsg.Data());
   return -4;
}

////////////////////////////////////////////////////////////////////////////////
/// Emit a signal for OnCdataBlock.

void TSAXParser::OnCdataBlock(const char *text, Int_t len)
{
   Long_t args[2];
   args[0] = (Long_t)text;
   args[1] = len;

   Emit("OnCdataBlock(const char *, Int_t)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// This function parses the xml file, by initializing the parser and checks
/// whether the parse context is created or not, it will check as well
/// whether the document is well formated.
/// It returns the parse error code, see TXMLParser.

Int_t TSAXParser::Parse()
{
   if (!fContext) {
      return -2;
   }

   xmlSAXHandlerPtr oldSAX = fContext->sax;
   fContext->sax = fSAXHandler;
   fContext->userData = this;

   InitializeContext();

   xmlParseDocument(fContext);

   fContext->sax = oldSAX;

   if (!fContext->wellFormed && fParseCode == 0) {
      fParseCode = -5;
   }

   ReleaseUnderlying();

   return fParseCode;
}

////////////////////////////////////////////////////////////////////////////////
/// It creates the parse context of the xml file, where the xml file name is
/// filename. If context is created sucessfully, it will call Parse()
/// It returns parse error code, see TXMLParser.

Int_t TSAXParser::ParseFile(const char *filename)
{
   // Attempt to parse a second file while a parse is in progress.
   if (fContext) {
      return -1;
   }

   fContext = xmlCreateFileParserCtxt(filename);
   return Parse();
}

////////////////////////////////////////////////////////////////////////////////
/// It parse the contents, instead of a file.
/// It will return error if is attempted to parse a second file while
/// a parse is in progres.
/// It returns parse code error, see TXMLParser.

Int_t TSAXParser::ParseBuffer(const char *contents, Int_t len)
{
   // Attempt to parse a second file while a parse is in progress.
   if (fContext) {
      return -1;
   }

   fContext = xmlCreateMemoryParserCtxt(contents, len);
   return Parse();
}


//--- TSAXParserCallback -------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// StartDocument Callback function.

void TSAXParserCallback::StartDocument(void *fParser)
{
   TSAXParser *parser = (TSAXParser*)fParser;
   parser->OnStartDocument();
}

////////////////////////////////////////////////////////////////////////////////
/// EndDocument callback function.

void TSAXParserCallback::EndDocument(void *fParser)
{
   TSAXParser *parser = (TSAXParser*)fParser;
   parser->OnEndDocument();
}

////////////////////////////////////////////////////////////////////////////////
/// StartElement callback function, where name is the name of the element
/// and p contains the attributes for the start tag.

void TSAXParserCallback::StartElement(void *fParser, const xmlChar *name,
                                      const xmlChar **p)
{
   TSAXParser *parser = (TSAXParser*)fParser;
   TList *attributes = new TList;

   if (p) {
      for (const xmlChar **cur = p; cur && *cur; cur += 2) {
         attributes->Add(new TXMLAttr((const char*)*cur,
                                      (const char*)*(cur + 1)));
      }
   }

   parser->OnStartElement((const char*) name, attributes);

   attributes->Delete();
   delete attributes;
}

////////////////////////////////////////////////////////////////////////////////
/// EndElement callback function, where name is the name of the element.

void TSAXParserCallback::EndElement(void *fParser, const xmlChar *name)
{
   TSAXParser *parser = (TSAXParser*)fParser;
   parser->OnEndElement((const char*) name);
}

////////////////////////////////////////////////////////////////////////////////
/// Character callback function. It is called when there are characters that
/// are outside of tags get parsed and the context will be stored in ch,
/// len is the length of ch.

void TSAXParserCallback::Characters(void *fParser, const xmlChar *ch,
                                    Int_t len)
{
   TSAXParser *parser = (TSAXParser*)fParser;

   char *str = new char[len+1];
   strlcpy(str, (const char*) ch, len+1);
   str[len] = '\0';

   parser->OnCharacters(str);

   delete [] str;
}

////////////////////////////////////////////////////////////////////////////////
/// Comment callback function.
/// Comment of the xml file will be parsed to value.

void TSAXParserCallback::Comment(void *fParser, const xmlChar *value)
{
   TSAXParser *parser = (TSAXParser*)fParser;
   parser->OnComment((const char*) value);
}

////////////////////////////////////////////////////////////////////////////////
/// Warning callback function. Warnings while parsing a xml file will
/// be stored at fmt.

void TSAXParserCallback::Warning(void * fParser, const char *va_(fmt), ...)
{
   TSAXParser *parser = (TSAXParser*)fParser;

   va_list arg;
   char buffer[2048];

   va_start(arg, va_(fmt));
   vsnprintf(buffer, 2048, va_(fmt), arg);
   va_end(arg);

   TString buff(buffer);

   parser->OnWarning(buff.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Error callback function. Errors while parsing a xml file will be stored
/// at fmt.

void TSAXParserCallback::Error(void *fParser, const char *va_(fmt), ...)
{
   Int_t errorcode;
   TSAXParser *parser = (TSAXParser*)fParser;

   va_list arg;
   char buffer[2048];

   va_start(arg, va_(fmt));
   vsnprintf(buffer, 2048, va_(fmt), arg);
   va_end(arg);

   TString buff(buffer);

   errorcode = parser->OnError(buff.Data());
   if (errorcode < 0) { //When error occurs, write fErrorCode
      parser->SetParseCode(errorcode);
   }

   if (errorcode < 0 && parser->GetStopOnError()) {
      //When GetStopOnError is enabled, stop the parse when an error occurs
      parser->StopParser();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// FactalError callback function. Factal errors while parsing a xml file
/// will be stored at fmt.

void TSAXParserCallback::FatalError(void *fParser, const char *va_(fmt), ...)
{
   Int_t errorcode;
   TSAXParser *parser = (TSAXParser*)fParser;

   va_list arg;
   char buffer[2048];

   va_start(arg, va_(fmt));
   vsnprintf(buffer, 2048, va_(fmt), arg);
   va_end(arg);

   TString buff(buffer);

   errorcode = parser->OnFatalError(buff);
   if (errorcode < 0) {
      parser->SetParseCode(errorcode);
      parser->StopParser();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// CdataBlock Callback function.

void TSAXParserCallback::CdataBlock(void *fParser, const xmlChar *value,
                                    Int_t len)
{
   TSAXParser *parser = (TSAXParser*)fParser;
   parser->OnCdataBlock((const char*)value, len);
}

////////////////////////////////////////////////////////////////////////////////
/// A default TSAXParser to a user-defined Handler connection function.
/// This function makes connection between various function from TSAXParser
/// with the user-define SAX Handler, whose functions has to be exactly the
/// same as in TSAXParser.
///
/// \param[in] handler Name User-defined SAX Handler class name
/// \param[in] handler Pointer to the user-defined SAX Handler
///
/// See SAXHandler.C tutorial.

void TSAXParser::ConnectToHandler(const char *handlerName, void *handler)
{
   const TString kFunctionsName [] = {
      "OnStartDocument()",
      "OnEndDocument()",
      "OnStartElement(const char *, const TList *)",
      "OnEndElement(const char *)",
      "OnCharacters(const char *)",
      "OnComment(const char *)",
      "OnWarning(const char *)",
      "OnError(const char *)",
      "OnFatalError(const char *)",
      "OnCdataBlock(const char *, Int_t)"
   };

   TClass *cl = TClass::GetClass(handlerName);

   for (Int_t i = 0; i < 10; i++) {
      if (CheckConnectArgs(this, this->IsA(), kFunctionsName[i],
                           cl, kFunctionsName[i]) != -1)
         Connect(kFunctionsName[i], handlerName, handler, kFunctionsName[i]);
   }
}
