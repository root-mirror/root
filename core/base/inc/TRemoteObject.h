// @(#)root/base:$Id$
// Author: Bertrand Bellenot   19/06/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TRemoteObject
#define ROOT_TRemoteObject

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRemoteObject                                                        //
//                                                                      //
// The TRemoteObject class provides protocol for browsing ROOT objects  //
// from a remote ROOT session.                                          //
// It contains information on the real remote object as:               //
//  - Object Properties (i.e. file stat if the object is a TSystemFile) //
//  - Object Name                                                       //
//  - Class Name                                                        //
//  - TKey Object Name (if the remote object is a TKey)                 //
//  - TKey Class Name (if the remote object is a TKey)                  //
//  - Remote object address                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TSystem.h"

class TList;

class TRemoteObject : public TNamed {

protected:
   FileStat_t  fFileStat;        // file status
   Bool_t      fIsFolder;        // is folder flag
   Long64_t    fRemoteAddress;   // remote address
   TString     fClassName;       // real object class name
   TString     fKeyObjectName;   // key object name
   TString     fKeyClassName;    // key object class name

public:
   TRemoteObject();
   TRemoteObject(const char *name, const char *title, const char *classname);

   virtual ~TRemoteObject();

   virtual void            Browse(TBrowser *b);
   Bool_t                  IsFolder() const { return fIsFolder; }
   TList                  *Browse();
   Bool_t                  GetFileStat(FileStat_t *sbuf);
   const char             *GetClassName() const { return fClassName.Data(); }
   const char             *GetKeyObjectName() const { return fKeyObjectName.Data(); }
   const char             *GetKeyClassName() const { return fKeyClassName.Data(); }
   void                    SetFolder(Bool_t isFolder) { fIsFolder = isFolder; }
   void                    SetKeyObjectName(const char *name) { fKeyObjectName = name; }
   void                    SetKeyClassName(const char *name) { fKeyClassName = name; }
   void                    SetRemoteAddress(Long_t addr) { fRemoteAddress = addr; }

   ClassDef(TRemoteObject,0)  //A remote object
};

#endif

