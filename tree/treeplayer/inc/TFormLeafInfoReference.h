// @(#)root/treeplayer:$Id$
// Author: Markus Frank 01/02/2006

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers and al.        *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFormLeafInfoReference
#define ROOT_TFormLeafInfoReference

#include "TFormLeafInfo.h"

// Forward declarations
class TVirtualRefProxy;

class TFormLeafInfoReference : public TFormLeafInfo {
   typedef TVirtualRefProxy Proxy;
public:
   Proxy*      fProxy;         //! Cached pointer to reference proxy
   TBranch*    fBranch;        //! Cached pointer to branch object
public:
   // Initializing constructor
   TFormLeafInfoReference(TClass* classptr, TStreamerElement* element, int off);
   // Copy constructor
   TFormLeafInfoReference(const TFormLeafInfoReference& orig);
   // Default destructor
   virtual ~TFormLeafInfoReference();
   // Exception safe swap.
   void Swap(TFormLeafInfoReference &other);
   // Exception safe assignment operator.
   TFormLeafInfoReference &operator=(const TFormLeafInfoReference &orig);
   // Virtual copy constructor
   virtual TFormLeafInfo* DeepCopy()  const;

   // Access to the info's proxy
   Proxy*           GetProxy()        const      {  return fProxy;        }
   // Access to the info's connected branch
   TBranch*         GetBranch()       const      {  return fBranch;       }
   // Access to the info's connected branch
   void             SetBranch(TBranch* branch)
   {  fBranch = branch; if ( fNext ) fNext->SetBranch(branch);            }
   // Access to the offset of the data
   virtual Int_t    GetOffset()       const     {  return fOffset;       }
   // Return true only if the underlying data is an integral value
   virtual Bool_t   IsInteger()       const     {  return kFALSE;         }
   // Return true only if the underlying data is a string
   virtual Bool_t   IsString()        const     {  return kFALSE;         }
   // Return true only if the underlying data is a reference
   virtual Bool_t   IsReference()     const     {  return kTRUE;          }
   // Access to target class pointer (if available)
   virtual TClass*  GetClass()        const;
   // Access to the value class of the reference proxy
   virtual TClass*  GetValueClass(TLeaf* from);
   // Access to the value class from the object pointer
   virtual TClass*  GetValueClass(void* from);
   // Return the address of the local value
   virtual void    *GetLocalValuePointer( TLeaf *from, Int_t instance = 0);
   // Return the address of the local value
   virtual void    *GetLocalValuePointer(char *from, Int_t instance = 0);
   // Return true if any of underlying data has a array size counter
   virtual Bool_t HasCounter() const;
   // Return the size of the underlying array for the current entry in the TTree.
   virtual Int_t ReadCounterValue(char *where);
   // Return the current size of the array container
   virtual Int_t GetCounterValue(TLeaf* leaf);

   // Access value of referenced object (macro from TFormLeafInfo.g)
   DECLARE_GETVAL;
   // Read value of referenced object
   DECLARE_READVAL;
   // TFormLeafInfo overload: Update (and propagate) cached information
   virtual Bool_t   Update();
};

#endif /* ROOT_TFormLeafInfoReference */
