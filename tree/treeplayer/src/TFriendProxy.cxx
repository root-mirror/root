// @(#)root/base:$Id$
// Author: Philippe Canal  13/05/2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TFriendProxy
Concrete implementation of the proxy around a Friend Tree.
*/

#include "TFriendProxy.h"
#include "TTree.h"
#include "TList.h"
#include "TFriendElement.h"

ClassImp(ROOT::Internal::TFriendProxy);

namespace ROOT {
namespace Internal {

/////////////////////////////////////////////////////////////////////////////

TFriendProxy::TFriendProxy() : fDirector(0,-1), fIndex(-1)
{
}

   /////////////////////////////////////////////////////////////////////////////
   /// Constructor.

   TFriendProxy::TFriendProxy(TBranchProxyDirector *director, TTree *main, Int_t index) :
      fDirector(0,-1), fIndex(index)
   {
      // The list of friends needs to be accessed via GetTree()->GetListOfFriends()
      // (and not directly GetListOfFriends()), otherwise when `main` is a TChain we
      // might not recover the list correctly (see #6993 for the TTreeReader issue
      // and #6741 for a more complete discussion/explanation).
      if (main && main->GetTree()->GetListOfFriends()) {
         TObject *obj = main->GetTree()->GetListOfFriends()->At(fIndex);
         TFriendElement *element = dynamic_cast<TFriendElement*>( obj );
         if (element) fDirector.SetTree(element->GetTree());
      }
      director->Attach(this);
   }

   /////////////////////////////////////////////////////////////////////////////
   /// Return the entry number currently being looked at.

   Long64_t TFriendProxy::GetReadEntry() const
   {
      return fDirector.GetReadEntry();
   }

   /////////////////////////////////////////////////////////////////////////////
   /// Refresh the cached read entry number from the original tree.

   void TFriendProxy::ResetReadEntry()
   {
      // The 2nd call to GetTree is to insure we get the 'local' tree's entry in the case of a
      // chain.
      if (fDirector.GetTree()) fDirector.SetReadEntry(fDirector.GetTree()->GetTree()->GetReadEntry());
   }

  //////////////////////////////////////////////////////////////////////////////
  /// Update the address of the underlying tree.

   void TFriendProxy::Update(TTree *newmain)
   {
      if (newmain && newmain->GetTree()->GetListOfFriends()) {
         TObject *obj = newmain->GetTree()->GetListOfFriends()->At(fIndex);
         TFriendElement *element = dynamic_cast<TFriendElement*>( obj );
         if (element) fDirector.SetTree(element->GetTree());
         else fDirector.SetTree(0);
      } else {
         fDirector.SetTree(0);
      }
   }

} // namespace Internal
} // namespace ROOT
