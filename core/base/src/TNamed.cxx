// @(#)root/base:$Id$
// Author: Rene Brun   26/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TNamed
\ingroup Base

The TNamed class is the base class for all named ROOT classes.

A TNamed contains the essential elements (name, title)
to identify a derived object in containers, directories and files.
Most member functions defined in this base class are in general
overridden by the derived classes.
*/

#include <iostream>
#include "Strlen.h"
#include "TNamed.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "TClass.h"

ClassImp(TNamed);

////////////////////////////////////////////////////////////////////////////////
/// TNamed copy ctor.

TNamed::TNamed(const TNamed &named) : TObject(named),fName(named.fName),fTitle(named.fTitle)
{
}

////////////////////////////////////////////////////////////////////////////////
/// TNamed destructor.

TNamed::~TNamed()
{
   // Required since we overload TObject::Hash.
   ROOT::CallRecursiveRemoveIfNeeded(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// TNamed assignment operator.

TNamed& TNamed::operator=(const TNamed& rhs)
{
   if (this != &rhs) {
      TObject::operator=(rhs);
      fName  = rhs.fName;
      fTitle = rhs.fTitle;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Set name and title to empty strings ("").

void TNamed::Clear(Option_t *)
{
   fName  = "";
   fTitle = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of an object using the Streamer facility.
/// If newname is specified, this will be the name of the new object.

TObject *TNamed::Clone(const char *newname) const
{
   TNamed *named = (TNamed*)TObject::Clone(newname);
   if (newname && strlen(newname)) named->SetName(newname);
   return named;
}

////////////////////////////////////////////////////////////////////////////////
/// Compare two TNamed objects. Returns 0 when equal, -1 when this is
/// smaller and +1 when bigger (like strcmp).

Int_t TNamed::Compare(const TObject *obj) const
{
   if (this == obj) return 0;
   return fName.CompareTo(obj->GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this to obj.

void TNamed::Copy(TObject &obj) const
{
   TObject::Copy(obj);
   ((TNamed&)obj).fName  = fName;
   ((TNamed&)obj).fTitle = fTitle;
}

////////////////////////////////////////////////////////////////////////////////
/// Encode TNamed into output buffer.

void TNamed::FillBuffer(char *&buffer)
{
   fName.FillBuffer(buffer);
   fTitle.FillBuffer(buffer);
}

////////////////////////////////////////////////////////////////////////////////
/// List TNamed name and title.

void TNamed::ls(Option_t *opt) const
{
   TROOT::IndentLevel();
   if (opt && strstr(opt,"noaddr")) {
      std::cout <<"OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << " : "
                << Int_t(TestBit(kCanDelete)) << std::endl;
   } else {
      std::cout <<"OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << " : "
                 << Int_t(TestBit(kCanDelete)) << " at: "<<this<< std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print TNamed name and title.

void TNamed::Print(Option_t *) const
{
   std::cout <<"OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the name of the TNamed.
///
/// WARNING: if the object is a member of a THashTable or THashList container
/// the container must be Rehash()'ed after SetName(). For example the list
/// of objects in the current directory is a THashList.

void TNamed::SetName(const char *name)
{
   fName = name;
   if (gPad && TestBit(kMustCleanup)) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Set all the TNamed parameters (name and title).
//
/// WARNING: if the name is changed and the object is a member of a
/// THashTable or THashList container the container must be Rehash()'ed
/// after SetName(). For example the list of objects in the current
/// directory is a THashList.

void TNamed::SetNameTitle(const char *name, const char *title)
{
   fName  = name;
   fTitle = title;
   if (gPad && TestBit(kMustCleanup)) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the title of the TNamed.

void TNamed::SetTitle(const char *title)
{
   fTitle = title;
   if (gPad && TestBit(kMustCleanup)) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Return size of the TNamed part of the TObject.

Int_t TNamed::Sizeof() const
{
   Int_t nbytes = fName.Sizeof() + fTitle.Sizeof();
   return nbytes;
}
