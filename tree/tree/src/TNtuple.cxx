// @(#)root/tree:$Id$
// Author: Rene Brun   06/04/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TNtuple.h"
#include "TBuffer.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TBrowser.h"
#include "TreeUtils.h"
#include "strlcpy.h"

#include <string>

ClassImp(TNtuple);

/** \class TNtuple
\ingroup tree

A simple TTree restricted to a list of float variables only.

Each variable goes to a separate branch.

A Ntuple is created via
~~~ {.cpp}
    TNtuple(name,title,varlist,bufsize)
~~~
It is filled via:
~~~ {.cpp}
    TNtuple::Fill(*x)  or
    TNtuple::Fill(v1,v2,v3.....)
~~~
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for Ntuple.

TNtuple::TNtuple(): TTree()
{
   fNvar = 0;
   fArgs = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create an Ntuple.
///
/// The parameter varlist describes the list of the ntuple variables
/// separated by a colon:
///
/// Example:  `x:y:z:energy`
///
/// For each variable in the list a separate branch is created.
///
/// NOTE:
///  - Use TTree to create branches with variables of different data types.
///  - Use TTree when the number of branches is large (> 100).

TNtuple::TNtuple(const char *name, const char *title, const char *varlist, Int_t bufsize)
       :TTree(name,title)
{
   Int_t i;
   fNvar = 0;
   fArgs = 0;

//   Count number of variables (separated by :)
   Int_t nch = strlen(varlist);
   if (nch == 0) return;
   char *vars = new char[nch+1];
   strlcpy(vars,varlist,nch+1);
   Int_t *pvars = new Int_t[nch+1];
   fNvar = 1;
   pvars[0] = 0;
   for (i=1;i<nch;i++) {
      if (vars[i] == ':') {
         pvars[fNvar] = i+1;
         vars[i] = 0;
         fNvar++;
      }
   }
   fArgs = new Float_t[fNvar];

//  Create one branch for each variable
   for (i=0;i<fNvar;i++) {
      Int_t pv = pvars[i];
      TTree::Branch(&vars[pv],&fArgs[i],&vars[pv],bufsize);
   }

   delete [] vars;
   delete [] pvars;
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for an Ntuple.

TNtuple::~TNtuple()
{
   delete [] fArgs;
   fArgs = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a clone of this tree and copy nentries.
///
/// By default copy all entries.
/// Note that only active branches are copied.
/// The compression level of the cloned tree is set to the destination file's
/// compression level.
///
/// See TTree::CloneTree for more details.

TTree* TNtuple::CloneTree(Long64_t nentries /* = -1 */, Option_t* option /* = "" */)
{
   TNtuple *newtuple = dynamic_cast<TNtuple*> (TTree::CloneTree(nentries,option) );
   if (newtuple) {
      // To deal with the cases of some of the branches where dropped.
      newtuple->fNvar = newtuple->fBranches.GetEntries();
   }
   return newtuple;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the branch addresses to the internal fArgs array. Use this
/// method when the addresses were changed via calls to SetBranchAddress().

void TNtuple::ResetBranchAddress(TBranch *branch)
{
   if (branch) {
      Int_t index = fBranches.IndexOf(branch);
      if (index>=0) {
         branch->SetAddress(&fArgs[index]);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the branch addresses to the internal fArgs array. Use this
/// method when the addresses were changed via calls to SetBranchAddress().

void TNtuple::ResetBranchAddresses()
{
   for (Int_t i = 0; i < fNvar; i++) {
      TBranch *branch = (TBranch*)fBranches.UncheckedAt(i);
      if (branch) branch->SetAddress(&fArgs[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Browse content of the ntuple

void TNtuple::Browse(TBrowser *b)
{
   fLeaves.Browse( b );
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Ntuple with current values in fArgs.
///
/// Note that this function is protected.
/// Currently called only by TChain::Merge

Int_t TNtuple::Fill()
{
   return TTree::Fill();
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Ntuple with an array of floats

Int_t TNtuple::Fill(const Float_t *x)
{

   // Store array x into buffer
   for (Int_t i=0;i<fNvar;i++)  {
      fArgs[i] = x[i];
   }

   return TTree::Fill();
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a Ntuple: Each Ntuple item is an argument

Int_t TNtuple::Fill(Float_t x0,Float_t x1,Float_t x2,Float_t x3,Float_t x4
              ,Float_t x5,Float_t x6,Float_t x7,Float_t x8,Float_t x9
              ,Float_t x10,Float_t x11,Float_t x12,Float_t x13,Float_t x14)
{
   if (fNvar >  0) fArgs[0]  = x0;
   if (fNvar >  1) fArgs[1]  = x1;
   if (fNvar >  2) fArgs[2]  = x2;
   if (fNvar >  3) fArgs[3]  = x3;
   if (fNvar >  4) fArgs[4]  = x4;
   if (fNvar >  5) fArgs[5]  = x5;
   if (fNvar >  6) fArgs[6]  = x6;
   if (fNvar >  7) fArgs[7]  = x7;
   if (fNvar >  8) fArgs[8]  = x8;
   if (fNvar >  9) fArgs[9]  = x9;
   if (fNvar > 10) fArgs[10] = x10;
   if (fNvar > 11) fArgs[11] = x11;
   if (fNvar > 12) fArgs[12] = x12;
   if (fNvar > 13) fArgs[13] = x13;
   if (fNvar > 14) fArgs[14] = x14;

   return TTree::Fill();
}

////////////////////////////////////////////////////////////////////////////////
/// Read from filename as many columns as variables in the ntuple
/// the function returns the number of rows found in the file
/// The second argument "branchDescriptor" is currently not used.
/// Lines in the input file starting with "#" are ignored.

Long64_t TNtuple::ReadStream(std::istream &inputStream, const char * /*branchDescriptor*/, char delimiter)
{
   /*
   Long64_t nlines = 0;
   char newline = GetNewlineValue(inputStream);
   while (1) {
      if ( inputStream.peek() != '#' ) {
         for (Int_t i=0;i<fNvar;i++) {
            inputStream >> fArgs[i];
            if (inputStream.peek() == delimiter) {
               inputStream.get(); // skip delimiter.
            }
         }
         if (!inputStream.good()) break;
         TTree::Fill();
         ++nlines;
      }
      inputStream.ignore(8192,newline);
   }
   return nlines;
   */

   //The last argument - true == strict mode.
   return ROOT::TreeUtils::FillNtupleFromStream<Float_t, TNtuple>(inputStream, *this, delimiter, true);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TNtuple::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         b.ReadClassBuffer(TNtuple::Class(), this, R__v, R__s, R__c);
      } else {
         //====process old versions before automatic schema evolution
         TTree::Streamer(b);
         b >> fNvar;
         b.CheckByteCount(R__s, R__c, TNtuple::IsA());
         //====end of old versions
      }
      if (fNvar <= 0) return;
      fArgs = new Float_t[fNvar];
      for (Int_t i=0;i<fNvar;i++) {
         TBranch *branch = (TBranch*)fBranches.UncheckedAt(i);
         if (branch) branch->SetAddress(&fArgs[i]);
      }
   } else {
      b.WriteClassBuffer(TNtuple::Class(),this);
   }
}
