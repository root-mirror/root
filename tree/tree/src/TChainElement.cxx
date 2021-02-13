// @(#)root/tree:$Id$
// Author: Rene Brun   11/02/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TChainElement
\ingroup tree

A TChainElement describes a component of a TChain.
*/

#include "TChainElement.h"
#include "TTree.h"
#include "TROOT.h"
#include <iostream>

ClassImp(TChainElement);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for a chain element.

TChainElement::TChainElement() : TNamed(),fBaddress(0),fBaddressType(0),
   fBaddressIsPtr(kFALSE), fDecomposedObj(kFALSE), fCheckedType(kFALSE), fBranchPtr(0), fLoadResult(0)
{
   fNPackets   = 0;
   fPackets    = 0;
   fEntries    = 0;
   fPacketSize = 100;
   fStatus     = -1;
   ResetBit(kHasBeenLookedUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a chain element.

TChainElement::TChainElement(const char *name, const char *title)
   :TNamed(name,title),fBaddress(0),fBaddressType(0),
    fBaddressIsPtr(kFALSE), fDecomposedObj(kFALSE), fCheckedType(kFALSE), fBranchPtr(0), fLoadResult(0)
{
   fNPackets   = 0;
   fPackets    = 0;
   fEntries    = 0;
   fPacketSize = 100;
   fStatus     = -1;
   ResetBit(kHasBeenLookedUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for a chain element.

TChainElement::~TChainElement()
{
   delete [] fPackets;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the packet descriptor string.

void TChainElement::CreatePackets()
{
   fNPackets = 1 + Int_t(fEntries/fPacketSize);
   delete [] fPackets;
   fPackets = new char[fNPackets+1];
   for (Int_t i=0;i<fNPackets;i++) fPackets[i] = ' ';
   fPackets[fNPackets] = 0;

}

////////////////////////////////////////////////////////////////////////////////
/// List files in the chain.

void TChainElement::ls(Option_t *) const
{
   TROOT::IndentLevel();
   std::cout << GetTitle() << "tree:" << GetName() << " entries=";
   if (fEntries == TTree::kMaxEntries)
      std::cout << "<not calculated>";
   else
      std::cout << fEntries;
   std::cout << '\n';
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of entries per packet for parallel root.

void TChainElement::SetPacketSize(Int_t size)
{
   fPacketSize = size;
}

////////////////////////////////////////////////////////////////////////////////
/// Set/Reset the looked-up bit

void TChainElement::SetLookedUp(Bool_t y)
{
   if (y)
      SetBit(kHasBeenLookedUp);
   else
      ResetBit(kHasBeenLookedUp);
}

