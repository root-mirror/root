// @(#)root/meta:$Id$
// Author: Markus Frank 20/05/2005

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TRefProxy
#define ROOT_TRefProxy

#include "TVirtualRefProxy.h"

#include "TClassRef.h"

// Forward declarations
class TFormLeafInfoReference;

//______________________________________________________________________________
//
// TRefProxy is a reference proxy, which allows to access ROOT references (TRef)
// stored contained in other objects from TTree::Draw
//______________________________________________________________________________

class TRefProxy : public TVirtualRefProxy  {
protected:
   TClassRef fClass;    //! Pointer to the reference class (TRef::Class())

public:
   /// Default constructor
   TRefProxy() : fClass("TRef") {}
   /// Copy constructor
   TRefProxy(const TRefProxy& c) : TVirtualRefProxy(), fClass(c.fClass) {}
   /// Assignement operator
   TRefProxy &operator=(const TRefProxy& c) { fClass =c.fClass; return *this; }

   /// TVirtualRefProxy overload: Release the reference proxy (virtual destructor)
   virtual void Release()                         { delete this;                }
   /// TVirtualRefProxy overload: Clone the reference proxy (virtual constructor)
   virtual TVirtualRefProxy* Clone() const        { return new TRefProxy(*this);}
   /// TVirtualRefProxy overload: Setter of reference class (executed when the proxy is adopted)
   virtual void SetClass(TClass *cl)              { fClass = cl;                }
   /// TVirtualRefProxy overload: Getter of reference class (executed when the proxy is adopted)
   virtual TClass * GetClass() const              { return fClass;              }
   /// TVirtualRefProxy overload: Access to value class
   virtual TClass* GetValueClass(void* data) const;
   /// TVirtualRefProxy overload: Prepare reused reference object (e.g. ZERO data pointers)
   virtual void* GetPreparedReference(void* data) {  return data;               }
   /// TVirtualRefProxy overload: Update (and propagate) cached information
   virtual Bool_t Update();
   /// TVirtualRefProxy overload: Flag to indicate if this is a container reference
   virtual Bool_t HasCounter()  const             { return kFALSE;              }
   /// TVirtualRefProxy overload: Access to container size (if container reference (ie TRefArray) etc)
   virtual Int_t  GetCounterValue(TFormLeafInfoReference* /* info */, void* /* data */)
   {  return 0;                                                                 }
   /// TVirtualRefProxy overload: Access referenced object(-data)
   virtual void* GetObject(TFormLeafInfoReference* info, void* data, int instance);
};
#endif // ROOT_TRefProxy
