// @(#)root/base:$Id$
// Author: Rene Brun   16/08/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMacro
#define ROOT_TMacro

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMacro                                                               //
//                                                                      //
// Class supporting a collection of lines with C++ code.                //
// A TMacro can be executed, saved to a ROOT file, edited, etc.         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TList;
class TObjString;
class TMD5;
class TCollection;


class TMacro : public TNamed {

protected:
   TList         *fLines;      //collection of lines
   TString        fParams;     //default string of macro parameters

   void           SaveSource(FILE *fp);

public:
   TMacro();
   TMacro(const TMacro&);
   TMacro(const char *name, const char *title="");
   virtual ~TMacro();
   TMacro& operator=(const TMacro&);
   virtual TObjString  *AddLine(const char *text);
   virtual void         Browse(TBrowser *b);
   virtual TMD5        *Checksum();
   virtual TObjString  *GetLineWith(const char *text) const;
   virtual Bool_t       Load() const; //*MENU*
   virtual Long_t       Exec(const char *params = 0, Int_t* error = 0); //*MENU*
   TList               *GetListOfLines() const {return fLines;}
   virtual const char  *GetParams() const { return fParams; }
   virtual Int_t        Merge(TCollection *list);
   virtual void         Paint(Option_t *option="");
   virtual void         Print(Option_t *option="") const;  //*MENU*
   virtual Int_t        ReadFile(const char *filename);
   virtual void         SaveSource(const char *filename);  //*MENU*
   virtual void         SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void         SetParams(const char *params=0); //*MENU*

   ClassDef(TMacro,1)  // Class supporting a collection of lines with C++ code.
};

#endif
