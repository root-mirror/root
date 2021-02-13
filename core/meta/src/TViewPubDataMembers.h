// @(#)root/meta:$Id$
// Author: Philippe Canal October 2013

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TViewPubDataMembers
#define ROOT_TViewPubDataMembers

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TViewPubDataMembers                                                  //
//                                                                      //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TList.h"

class TClass;


class TViewPubDataMembers : public TList {

protected:
   TList   fClasses; // list of the all the (base) classes for which we list methods.

private:
   TViewPubDataMembers(const TViewPubDataMembers&) = delete;
   TViewPubDataMembers& operator=(const TViewPubDataMembers&) = delete;

public:
   TViewPubDataMembers(TClass *cl = nullptr);
   virtual    ~TViewPubDataMembers();

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
   ClassDefInline(TViewPubDataMembers, 0)
};

// Preventing warnings with -Weffc++ in GCC since it is a false positive for the TListIter destructor.
#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) >= 40600
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TViewPubDataMembersIter                                              //
//                                                                      //
// Iterator of view of linked list.      `                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TViewPubDataMembersIter : public TIterator,
                                public std::iterator<std::bidirectional_iterator_tag, TObject *, std::ptrdiff_t,
                                                     const TObject **, const TObject *&> {
protected:
   const TList *fView;   //View we are iterating over.
   TIter        fClassIter;    //iterator over the classes
   TIter        fIter;         //iterator over the members of the current class
   Bool_t       fStarted;      //iteration started
   Bool_t       fDirection;    //iteration direction

   TViewPubDataMembersIter() : fView(0), fClassIter((TCollection *)0), fIter((TCollection *)0), fStarted(kFALSE), fDirection(kIterForward) { }

public:
   TViewPubDataMembersIter(const TViewPubDataMembers *l, Bool_t dir = kIterForward);
   TViewPubDataMembersIter(const TViewPubDataMembersIter &iter);
   ~TViewPubDataMembersIter() { }
   TIterator &operator=(const TIterator &rhs);
   TViewPubDataMembersIter &operator=(const TViewPubDataMembersIter &rhs);

   const TCollection *GetCollection() const { return fView; }
   TObject           *Next();
   void               Reset();
   Bool_t             operator!=(const TIterator &aIter) const;
   Bool_t             operator!=(const TViewPubDataMembersIter &aIter) const;
   TObject           *operator*() const { return *fIter; }

   // ClassDefInline does not yet support non default constructible classes
   //    ClassDefInline(TViewPubDataMembersIter,0)
};

#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) >= 40600
#pragma GCC diagnostic pop
#endif

#endif // ROOT_TViewPubDataMembers
