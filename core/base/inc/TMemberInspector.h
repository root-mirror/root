// @(#)root/base:$Id$
// Author: Fons Rademakers   15/07/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMemberInspector
#define ROOT_TMemberInspector

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMemberInspector                                                     //
//                                                                      //
// Abstract base class for accessing the datamembers of a class.        //
// Classes derived from this class can be given as argument to the      //
// ShowMembers() methods of ROOT classes. This feature facilitates      //
// the writing of class browsers and inspectors.                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

class TObject;
class TClass;

class TMemberInspector {
public:
   enum EObjectPointerState {
      kUnset, // No Inspect() call has been seen yet.
      kNoObjectGiven, // No object was given to the initial Inspect() call.
      kValidObjectGiven // The address points to an actual object.
   };
private:
   class TParentBuf;
   TParentBuf* fParent; // current inspection "path"
   EObjectPointerState fObjectPointerState; // whether the address is valid or only an offset

   TMemberInspector(const TMemberInspector &) = delete;
   TMemberInspector &operator=(const TMemberInspector &) = delete;

public:
   TMemberInspector();
   virtual ~TMemberInspector();

   EObjectPointerState GetObjectValidity() const { return fObjectPointerState; }
   void SetObjectValidity(EObjectPointerState val) { fObjectPointerState = val; }
   virtual void Inspect(TClass *cl, const char *parent, const char *name, const void *addr);
   virtual void Inspect(TClass *cl, const char *parent, const char *name, const void *addr, Bool_t /* isTransient */) { Inspect(cl,parent,name,addr); }

   const char* GetParent() const;
   Ssiz_t GetParentLen() const;
   void AddToParent(const char* name);
   void RemoveFromParent(Ssiz_t startingAt);

   virtual Bool_t IsTreatingNonAccessibleTypes() {return kTRUE;}

   template <class T>
   void InspectMember(const T& obj, const char* name, Bool_t isTransient) {
      Ssiz_t len = GetParentLen();
      AddToParent(name);
      obj.IsA()->CallShowMembers(&obj, *this, isTransient);
      RemoveFromParent(len);
   }

   void InspectMember(const TObject& obj, const char* name, Bool_t isTransient);
   void InspectMember(const char* topclassname, const void* pobj, const char* name,
                      Bool_t transient);
   void InspectMember(TClass* cl, const void* pobj, const char* name,
                      Bool_t isTransient);

   void GenericShowMembers(const char *topClassName, const void *obj,
                           Bool_t transientMember);

   ClassDef(TMemberInspector,0)  //ABC for inspecting class data members
};

#endif
