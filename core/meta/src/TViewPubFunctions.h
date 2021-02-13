// @(#)root/meta:$Id$
// Author: Philippe Canal October 2013

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TViewPubFunctions
#define ROOT_TViewPubFunctions

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TViewPubFunctions                                                    //
//                                                                      //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TList.h"

class TClass;


class TViewPubFunctions : public TList {

protected:
   TList   fClasses; // list of the all the (base) classes for which we list methods.

private:
   TViewPubFunctions(const TViewPubFunctions&) = delete;
   TViewPubFunctions& operator=(const TViewPubFunctions&) = delete;

public:
   TViewPubFunctions(TClass *cl = nullptr);
   virtual    ~TViewPubFunctions();

   TObject   *FindObject(const char *name) const;
   TObject   *FindObject(const TObject *obj) const;

   virtual TObject  *At(Int_t idx) const;
   virtual const TList *GetListOfClasses() const { return &fClasses; }
   virtual Int_t        GetSize() const;
   virtual TIterator   *MakeIterator(Bool_t dir = kIterForward) const;

   void       Load();

   // All the following routines are explicitly disallow/unsupported for
   // a view
protected:
   void       Clear(Option_t *option="");
   void       Delete(Option_t *option="");

   void       AddFirst(TObject *obj);
   void       AddFirst(TObject *obj, Option_t *opt);
   void       AddLast(TObject *obj);
   void       AddLast(TObject *obj, Option_t *opt);
   void       AddAt(TObject *obj, Int_t idx);
   void       AddAfter(const TObject *after, TObject *obj);
   void       AddAfter(TObjLink *after, TObject *obj);
   void       AddBefore(const TObject *before, TObject *obj);
   void       AddBefore(TObjLink *before, TObject *obj);

   virtual TObject  *After(const TObject *obj) const;
   virtual TObject  *Before(const TObject *obj) const;
   virtual TObject  *First() const;
   virtual TObjLink *FirstLink() const;
   virtual TObject **GetObjectRef(const TObject *obj) const;
   virtual TObject  *Last() const;
   virtual TObjLink *LastLink() const;

   void       RecursiveRemove(TObject *obj);
   TObject   *Remove(TObject *obj);
   TObject   *Remove(TObjLink *lnk);

public:
   ClassDefInline(TViewPubFunctions, 0) // Doubly linked list with hashtable for lookup
};

// Preventing warnings with -Weffc++ in GCC since it is a false positive for the TListIter destructor.
#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) >= 40600
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TViewPubFunctionsIter                                                //
//                                                                      //
// Iterator of view of linked list.      `1234                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TViewPubFunctionsIter : public TIterator,
                              public std::iterator<std::bidirectional_iterator_tag,
                                                   TObject*, std::ptrdiff_t,
                                                    const TObject**, const TObject*&>
{
protected:
   const TList *fView;   //View we are iterating over.
   TIter        fClassIter;    //iterator over the classes
   TIter        fFuncIter;     //iterator over the method of the current class
   Bool_t       fStarted;      //iteration started
   Bool_t       fDirection;    //iteration direction

   TViewPubFunctionsIter() : fView(0), fClassIter((TCollection *)0), fFuncIter((TCollection *)0), fStarted(kFALSE), fDirection(kIterForward) { }

public:
   TViewPubFunctionsIter(const TViewPubFunctions *l, Bool_t dir = kIterForward);
   TViewPubFunctionsIter(const TViewPubFunctionsIter &iter);
   ~TViewPubFunctionsIter() { }
   TIterator &operator=(const TIterator &rhs);
   TViewPubFunctionsIter &operator=(const TViewPubFunctionsIter &rhs);

   const TCollection *GetCollection() const { return fView; }
   TObject           *Next();
   void               Reset();
   Bool_t             operator!=(const TIterator &aIter) const;
   Bool_t             operator!=(const TViewPubFunctionsIter &aIter) const;
   TObject           *operator*() const { return *fFuncIter; }

   // ClassDefInline does not yet support non default constructible classes
   // ClassDefInline(TViewPubFunctionsIter,0)  //Linked list iterator
};

#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) >= 40600
#pragma GCC diagnostic pop
#endif

#endif // ROOT_TViewPubFunctions
