// @(#)root/xmlparser:$Id$
// Author: Jose Lo   12/4/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TXMLDocument
#define ROOT_TXMLDocument

#include "TObject.h"

struct _xmlDoc;
class TXMLNode;


class TXMLDocument : public TObject {

private:
   TXMLDocument(const TXMLDocument&) = delete;
   TXMLDocument& operator=(const TXMLDocument&) = delete;

   _xmlDoc  *fXMLDoc;           // libxml xml doc
   TXMLNode *fRootNode;         // the root node

public:
   TXMLDocument(_xmlDoc *doc);
   virtual ~TXMLDocument();

   TXMLNode   *GetRootNode() const;

   const char *Version() const;
   const char *Encoding() const;
   const char *URL() const;

   ClassDef(TXMLDocument,0)  // XML document created by the DOM parser
};

#endif
