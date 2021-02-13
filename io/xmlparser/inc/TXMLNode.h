// @(#)root/xmlparser:$Id$
// Author: Jose Lo   12/4/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TXMLNode
#define ROOT_TXMLNode

#include "TObject.h"

class TList;
struct _xmlNode;

class TXMLNode : public TObject {

private:
   TXMLNode(const TXMLNode&) = delete;
   TXMLNode& operator=(const TXMLNode&) = delete;

   _xmlNode *fXMLNode;        ///< libxml node

   TXMLNode *fParent;         ///< Parent node
   TXMLNode *fChildren;       ///< Children node
   TXMLNode *fNextNode;       ///< Next sibling node
   TXMLNode *fPreviousNode;   ///< Previous sibling node
   TList    *fAttrList;       ///< List of Attributes

public:
   /// This enum is based on libxml tree Enum xmlElementType
   enum EXMLElementType {
      kXMLElementNode = 1,
      kXMLAttributeNode = 2,
      kXMLTextNode = 3,
      kXMLCommentNode = 8
   };

   TXMLNode(_xmlNode *node, TXMLNode* parent=0, TXMLNode* previous=0);

   virtual ~TXMLNode();

   EXMLElementType GetNodeType() const;
   const char *GetNodeName() const;
   TXMLNode   *GetChildren();
   TXMLNode   *GetParent() const;
   TXMLNode   *GetNextNode();
   TXMLNode   *GetPreviousNode() const;
   const char *GetContent() const;
   const char *GetText() const;
   TList      *GetAttributes();

   Bool_t      HasChildren() const;
   Bool_t      HasNextNode() const;
   Bool_t      HasParent() const;
   Bool_t      HasPreviousNode() const;
   Bool_t      HasAttributes() const;

   const char *GetNamespaceHref() const;
   const char *GetNamespacePrefix() const;

   ClassDef(TXMLNode,0);  // XML node under DOM tree
};

#endif
