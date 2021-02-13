// @(#)root/cont:$Id$
// Author: Fons Rademakers   10/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TList
#define ROOT_TList


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TList                                                                //
//                                                                      //
// A doubly linked list. All classes inheriting from TObject can be     //
// inserted in a TList.                                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TSeqCollection.h"
#include "TString.h"

#include <iterator>
#include <memory>

#if (__GNUC__ >= 3) && !defined(__INTEL_COMPILER)
// Prevent -Weffc++ from complaining about the inheritance
// TListIter from std::iterator.
#pragma GCC system_header
#endif

const Bool_t kSortAscending  = kTRUE;
const Bool_t kSortDescending = !kSortAscending;

class TObjLink;
class TListIter;


class TList : public TSeqCollection {

friend  class TListIter;

protected:
   using TObjLinkPtr_t = std::shared_ptr<TObjLink>;
   using TObjLinkWeakPtr_t = std::weak_ptr<TObjLink>;

   TObjLinkPtr_t     fFirst;     //! pointer to first entry in linked list
   TObjLinkPtr_t     fLast;      //! pointer to last entry in linked list
   TObjLinkWeakPtr_t fCache;     //! cache to speedup sequential calling of Before() and After() functions
   Bool_t     fAscending; //! sorting order (when calling Sort() or for TSortedList)

   TObjLink          *LinkAt(Int_t idx) const;
   TObjLink          *FindLink(const TObject *obj, Int_t &idx) const;

   TObjLinkPtr_t *DoSort(TObjLinkPtr_t *head, Int_t n);

   Bool_t         LnkCompare(const TObjLinkPtr_t &l1, const TObjLinkPtr_t &l2);
   TObjLinkPtr_t  NewLink(TObject *obj, const TObjLinkPtr_t &prev = nullptr);
   TObjLinkPtr_t  NewOptLink(TObject *obj, Option_t *opt,  const TObjLinkPtr_t &prev = nullptr);
   TObjLinkPtr_t  NewLink(TObject *obj, TObjLink *prev);
   TObjLinkPtr_t  NewOptLink(TObject *obj, Option_t *opt,  TObjLink *prev);
   // virtual void       DeleteLink(TObjLink *lnk);

   void InsertAfter(const TObjLinkPtr_t &newlink, const TObjLinkPtr_t &prev);

private:
   TList(const TList&) = delete;
   TList& operator=(const TList&) = delete;

public:
   typedef TListIter Iterator_t;

   TList() : fAscending(kTRUE) { }
   TList(TObject *) : fAscending(kTRUE) { } // for backward compatibility, don't use
   virtual           ~TList();
   virtual void      Clear(Option_t *option="");
   virtual void      Delete(Option_t *option="");
   virtual TObject  *FindObject(const char *name) const;
   virtual TObject  *FindObject(const TObject *obj) const;
   virtual TIterator *MakeIterator(Bool_t dir = kIterForward) const;

   virtual void      Add(TObject *obj) { AddLast(obj); }
   virtual void      Add(TObject *obj, Option_t *opt) { AddLast(obj, opt); }
   virtual void      AddFirst(TObject *obj);
   virtual void      AddFirst(TObject *obj, Option_t *opt);
   virtual void      AddLast(TObject *obj);
   virtual void      AddLast(TObject *obj, Option_t *opt);
   virtual void      AddAt(TObject *obj, Int_t idx);
   virtual void      AddAfter(const TObject *after, TObject *obj);
   virtual void      AddAfter(TObjLink *after, TObject *obj);
   virtual void      AddBefore(const TObject *before, TObject *obj);
   virtual void      AddBefore(TObjLink *before, TObject *obj);
   virtual TObject  *Remove(TObject *obj);
   virtual TObject  *Remove(TObjLink *lnk);
           TObject  *Remove(const TObjLinkPtr_t &lnk) { return Remove(lnk.get()); }
   virtual void      RemoveLast();
   virtual void      RecursiveRemove(TObject *obj);

   virtual TObject  *At(Int_t idx) const;
   virtual TObject  *After(const TObject *obj) const;
   virtual TObject  *Before(const TObject *obj) const;
   virtual TObject  *First() const;
   virtual TObjLink *FirstLink() const { return fFirst.get(); }
   virtual TObject **GetObjectRef(const TObject *obj) const;
   virtual TObject  *Last() const;
   virtual TObjLink *LastLink() const { return fLast.get(); }

   virtual void      Sort(Bool_t order = kSortAscending);
   Bool_t            IsAscending() { return fAscending; }

   ClassDef(TList,5)  //Doubly linked list
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TObjLink                                                             //
//                                                                      //
// Wrapper around a TObject so it can be stored in a TList.             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TObjLink : public std::enable_shared_from_this<TObjLink> {

friend class TList;

private:
   using TObjLinkPtr_t = std::shared_ptr<TObjLink>;
   using TObjLinkWeakPtr_t = std::weak_ptr<TObjLink>;

   TObjLinkPtr_t     fNext;
   TObjLinkWeakPtr_t fPrev;

   TObject    *fObject; // should be atomic ...

   TObjLink(const TObjLink&) = delete;
   TObjLink& operator=(const TObjLink&) = delete;
   TObjLink() = delete;

public:

   TObjLink(TObject *obj) : fObject(obj) { }
   virtual ~TObjLink() { }

   TObject                *GetObject() const { return fObject; }
   TObject               **GetObjectRef() { return &fObject; }
   void                    SetObject(TObject *obj) { fObject = obj; }
   virtual Option_t       *GetAddOption() const { return ""; }
   virtual Option_t       *GetOption() const { return fObject->GetOption(); }
   virtual void            SetOption(Option_t *) { }
   TObjLink               *Next() { return fNext.get(); }
   TObjLink               *Prev() { return fPrev.lock().get(); }
   TObjLinkPtr_t           NextSP() { return fNext; }
   TObjLinkPtr_t           PrevSP() { return fPrev.lock(); }
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TObjOptLink                                                          //
//                                                                      //
// Wrapper around a TObject so it can be stored in a TList including    //
// an option string.                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TObjOptLink : public TObjLink {

private:
   TString   fOption;

public:
   TObjOptLink(TObject *obj, Option_t *opt) : TObjLink(obj), fOption(opt) { }
   ~TObjOptLink() { }
   Option_t        *GetAddOption() const { return fOption.Data(); }
   Option_t        *GetOption() const { return fOption.Data(); }
   void             SetOption(Option_t *option) { fOption = option; }
};


// Preventing warnings with -Weffc++ in GCC since it is a false positive for the TListIter destructor.
#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) >= 40600
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TListIter                                                            //
//                                                                      //
// Iterator of linked list.                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TListIter : public TIterator,
                  public std::iterator<std::bidirectional_iterator_tag,
                                       TObject*, std::ptrdiff_t,
                                       const TObject**, const TObject*&> {

protected:
   using TObjLinkPtr_t = std::shared_ptr<TObjLink>;

   const TList   *fList;         //list being iterated
   TObjLinkPtr_t  fCurCursor;    //current position in list
   TObjLinkPtr_t  fCursor;       //next position in list
   Bool_t         fDirection;    //iteration direction
   Bool_t         fStarted;      //iteration started

   TListIter() : fList(nullptr), fCurCursor(), fCursor(), fDirection(kIterForward),
                 fStarted(kFALSE) { }

public:
   TListIter(const TList *l, Bool_t dir = kIterForward);
   TListIter(const TListIter &iter);
   ~TListIter() { }
   TIterator &operator=(const TIterator &rhs);
   TListIter &operator=(const TListIter &rhs);

   const TCollection *GetCollection() const { return fList; }
   Option_t          *GetOption() const;
   void               SetOption(Option_t *option);
   TObject           *Next();
   void               Reset();
   Bool_t             operator!=(const TIterator &aIter) const;
   Bool_t             operator!=(const TListIter &aIter) const;
   TObject           *operator*() const { return (fCurCursor ? fCurCursor->GetObject() : nullptr); }

   ClassDef(TListIter,0)  //Linked list iterator
};

inline bool operator==(TObjOptLink *l, const std::shared_ptr<TObjLink> &r) {
   return l == r.get();
}

inline bool operator==(const std::shared_ptr<TObjLink> &l, TObjOptLink *r) {
    return r == l;
}

inline TList::TObjLinkPtr_t  TList::NewLink(TObject *obj, TObjLink *prev) {
   return NewLink(obj, prev ? prev->shared_from_this() : nullptr);
}
inline TList::TObjLinkPtr_t  TList::NewOptLink(TObject *obj, Option_t *opt,  TObjLink *prev) {
   return NewOptLink(obj, opt, prev ? prev->shared_from_this() : nullptr);
}


#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) >= 40600
#pragma GCC diagnostic pop
#endif

#endif
