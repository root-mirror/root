// @(#)root/tree:$Id$
// Author: Axel Naumann   14/10/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBranchBrowsable
#define ROOT_TBranchBrowsable

#include "TNamed.h"

#include "TList.h"

#include <list>

class TMethod;
class TBowser;
class TClass;
class TBranch;
class TBranchElement;
class TString;
class TStreamerElement;

class TVirtualBranchBrowsable: public TNamed {
public:

   // these methods are registered in RegisterGenerator, and
   // called to create the list of browsables. See e.g.
   // TMethodBrowsable::Register
   typedef Int_t (*MethodCreateListOfBrowsables_t)
      (TList &, const TBranch *branch, const TVirtualBranchBrowsable *parent);

   ~TVirtualBranchBrowsable();

   void Browse(TBrowser *b) override;

   /** return icon shown when browsing a TVirtualBranchBrowsable */
   const char* GetIconName() const override
   {
      if (IsFolder())
         return "TBranchElement-folder";
      else
         return "TBranchElement-leaf";
   }

   void GetScope(TString &scope) const;

   /** check whether we have sub-elements */
   Bool_t IsFolder() const override { return (GetLeaves() && GetLeaves()->GetSize()); }

   static Int_t FillListOfBrowsables(TList &list, const TBranch *branch,
         const TVirtualBranchBrowsable *parent = nullptr);

   /** return the parent branch (might be many levels up) */
   const TBranch* GetBranch() const { return fBranch; }

   /** return the parent TVirtualBranchBrowsable */
   const TVirtualBranchBrowsable* GetParent() const { return fParent; }

   /** return the type of this browsable object */
   TClass *GetClassType() const { return fClass; }

   /** return whether the type of this browsable object is a pointer */
   Bool_t TypeIsPointer() const { return fTypeIsPointer; }

   TList *GetLeaves() const;

   // static void Register()   has to be implemented for all derived classes!
   // static void Unregister() has to be implemented for all derived classes!

protected:
   TVirtualBranchBrowsable(const TBranch *b, TClass *type, Bool_t typeIsPointer,
         const TVirtualBranchBrowsable *parent = nullptr);
   static TClass* GetCollectionContainedType(const TBranch *b,
         const TVirtualBranchBrowsable *parent, TClass *&contained);
   static std::list<MethodCreateListOfBrowsables_t>& GetRegisteredGenerators();
   static void RegisterGenerator(MethodCreateListOfBrowsables_t generator);
   static void UnregisterGenerator(MethodCreateListOfBrowsables_t generator);
   /** sets the type of this browsable object */
   void SetType(TClass* type) { fClass = type; }

   /** sets whether the type of this browsable object is a pointer */
   void SetTypeIsPointer(Bool_t set=kTRUE) { fTypeIsPointer = set; }

private:
   static void RegisterDefaultGenerators();
   const TBranch    *fBranch{nullptr}; ///< pointer to the branch element representing the top object
   const TVirtualBranchBrowsable *fParent{nullptr}; ///< parent method if this method is member of a returned class
   TList            *fLeaves{nullptr}; ///< pointer to leaves
   TClass           *fClass{nullptr};  ///< pointer to TClass representing our type (i.e. return type for methods), 0 if basic type
   Bool_t            fTypeIsPointer{kFALSE}; ///< return type is pointer to class
   static std::list<MethodCreateListOfBrowsables_t> fgGenerators; ///< list of MethodCreateListOfBrowsables_t called by CreateListOfBrowsables
   static Bool_t     fgGeneratorsSet; ///< have we set the generators yet? empty is not good enough - user might have removed them
   ClassDefOverride(TVirtualBranchBrowsable, 0); ///< Base class for helper objects used for browsing
};


class TMethodBrowsable: public TVirtualBranchBrowsable {
public:
   ~TMethodBrowsable() {};

   static Int_t GetBrowsables(TList &list, const TBranch *branch,
         const TVirtualBranchBrowsable *parent = nullptr);

   /** return our special icons */
   const char* GetIconName() const override
   {
      if (IsFolder())
         return "TMethodBrowsable-branch";
      return "TMethodBrowsable-leaf";
   }
   static Bool_t IsMethodBrowsable(const TMethod *m);
   static void Register();
   static void Unregister();

protected:
   static void GetBrowsableMethodsForClass(TClass *cl, TList &list);
   TMethodBrowsable(const TBranch *branch, TMethod *m,
         const TVirtualBranchBrowsable *parent = nullptr);

private:
   TMethod         *fMethod{nullptr}; // pointer to a method
   ClassDefOverride(TMethodBrowsable,0); // Helper object to browse methods
};


class TNonSplitBrowsable: public TVirtualBranchBrowsable {
public:
   ~TNonSplitBrowsable() {}

   static Int_t GetBrowsables(TList &list, const TBranch *branch,
         const TVirtualBranchBrowsable *parent = nullptr);
   static void Register();
   static void Unregister();

protected:
   TNonSplitBrowsable(const TStreamerElement *element, const TBranch *branch,
         const TVirtualBranchBrowsable *parent = nullptr);

private:
   ClassDef(TNonSplitBrowsable, 0); // Helper object to browse unsplit objects
};


class TCollectionPropertyBrowsable: public TVirtualBranchBrowsable {
public:
   ~TCollectionPropertyBrowsable() {}

   void Browse(TBrowser *b) override;
   static Int_t GetBrowsables(TList &list, const TBranch *branch,
         const TVirtualBranchBrowsable *parent = nullptr);
   /** return the string passed to TTree::Draw */
   const char *GetDraw() const { return fDraw.Data(); }
   static void Register();
   static void Unregister();

protected:
   /** constructor, which sets the name and title according to the parameters
     * (and thus differently than our base class TVirtualBranchBrowsable) */
   TCollectionPropertyBrowsable(const char *name, const char *title,
         const char *draw, const TBranch *branch,
         const TVirtualBranchBrowsable *parent = nullptr) :
         TVirtualBranchBrowsable(branch, 0, kFALSE, parent), fDraw(draw)
   {
      SetNameTitle(name, title);
   }

private:
   TString fDraw; // string to send to TTree::Draw(), NOT by GetScope()!
   ClassDefOverride(TCollectionPropertyBrowsable, 0); // Helper object to add browsable collection properties
};

class TCollectionMethodBrowsable: public TMethodBrowsable {
public:
   ~TCollectionMethodBrowsable() {};

   static Int_t GetBrowsables(TList &list, const TBranch *branch,
         const TVirtualBranchBrowsable *parent = nullptr);
   static void Register();
   static void Unregister();

protected:
   TCollectionMethodBrowsable(const TBranch* branch, TMethod* m,
      const TVirtualBranchBrowsable* parent = nullptr);

   ClassDef(TCollectionMethodBrowsable,0); // Helper object to browse a collection's methods
};

#endif // defined ROOT_TBranchBrowsable
