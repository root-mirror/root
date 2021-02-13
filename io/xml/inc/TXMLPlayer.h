// @(#)root/xml:$Id$
// Author: Sergey Linev  10.05.2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TXMLPlayer
#define ROOT_TXMLPlayer

#include "TObject.h"

#include "TXMLSetup.h"

class TStreamerInfo;
class TStreamerElement;
class TStreamerSTL;
class TDataMember;
class TList;

class TXMLPlayer : public TObject {
public:
   TXMLPlayer();
   virtual ~TXMLPlayer();

   Bool_t ProduceCode(TList *cllist, const char *filename);

protected:
   TString GetStreamerName(TClass *cl);

   const char *ElementGetter(TClass *cl, const char *membername, int specials = 0);
   const char *ElementSetter(TClass *cl, const char *membername, char *endch);

   TString GetMemberTypeName(TDataMember *member);
   TString GetBasicTypeName(TStreamerElement *el);
   TString GetBasicTypeReaderMethodName(Int_t type, const char *realname);
   void ProduceStreamerSource(std::ostream &fs, TClass *cl, TList *cllist);

   void ReadSTLarg(std::ostream &fs, TString &argname, int argtyp, Bool_t isargptr, TClass *argcl, TString &tname,
                   TString &ifcond);
   void WriteSTLarg(std::ostream &fs, const char *accname, int argtyp, Bool_t isargptr, TClass *argcl);
   Bool_t ProduceSTLstreamer(std::ostream &fs, TClass *cl, TStreamerSTL *el, Bool_t isWriting);

   TString fGetterName; //!  buffer for name of getter method
   TString fSetterName; //!  buffer for name of setter method
   TXMLSetup fXmlSetup; //!  buffer for xml names conversion

   ClassDef(TXMLPlayer, 1) // Generation of external xml streamers
};

#endif
