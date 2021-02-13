// @(#)root/meta: $Id$
// Author: Markus Frank 20/05/2005

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TVirtualRefProxy
#define ROOT_TVirtualRefProxy

// Framework include files
#include "RtypesCore.h"

// Forward declarations
class TClass;
class TFormLeafInfoReference;

//______________________________________________________________________________
//
//   Abstract proxy definition to follow reference objects.
//
//
//   Generic Mechanism for Object References
//   =======================================
//
//   References are a well known mechanism to support persistency
//   of entities, which in C++ typically are represented as
//   pointers. The generic mechanism allows clients to supply
//   hooks to the ROOT framework in interactive mode in order to
//   dereference these objects and access the objects pointed to by
//   the reference objects.
//
//   Implementations are supplied for ROOT own reference mechanism
//   based on instances of the TRef and the TRefArray classes.
//
//   To support generality this mechanism was implemented using a
//   proxy mechanism, which shields the concrete implementation of the
//   reference classes from ROOT. Hence, this mechanism also works for
//   references as they are supported by the POOL persistency framework
//   and by frameworks like Gaudi.
//
//   To enable reference support a concrete sub-class instance of
//   the TVirtualRefProxy base class must be attached to the TClass
//   instance representing the reference itself. Please see the
//   header- and implementation file TRefProxy.h/cxx for details.
//   For ROOT's own references this is done simply by a call like:
//
//        #include "TROOT.h"
//        #include "TClass.h"
//        #include "TRefProxy.h"
//
//        ...
//        gROOT->GetClass("TRef")->AdoptReferenceProxy(new TRefProxy());
//
//      - GetObject() must return the pointer to the referenced
//        object. TTreeFormula then figures out how to access the
//        value to be plotted.
//        Hence, the actual work is done inside a call to:
//
//        void* TRefProxy::GetObject(TFormLeafInfoReference* info, void* data, int)
//        {
//          if ( data )  {
//            TRef*      ref    = (TRef*)((char*)data + info->GetOffset());
//            // Dereference TRef and return pointer to object
//            void* obj = ref->GetObject();
//            if ( obj )  {         return obj;      }
//
//            ... else handle error or implement failover ....
//
//
//   The type of the referenced object must either be known at compilation
//      time or it must be possible to guess it reading the first TTree entry.
//      In this case the following conditions must be met:
//      - GetValueClass() must return the TClass to the referenced
//        objects (or a base class)
//
//______________________________________________________________________________
class TVirtualRefProxy  {
public:
   // Virtual Destructor
   virtual ~TVirtualRefProxy() {};

   // Release the reference proxy (virtual destructor)
   virtual void Release() = 0;

   // Clone the reference proxy (virtual constructor)
   virtual TVirtualRefProxy* Clone() const = 0;

   // Setter of reference class (executed when the proxy is adopted)
   // Setup the reference when it is adopted by the TClass structure
   //
   // classptr [IN]    Pointer to the reference class.
   virtual void SetClass(TClass *classptr) = 0;

   // Getter of reference class.
   // The function returns the class description of the reference class
   // ie. in the case of TRef TRef::Class
   virtual TClass * GetClass() const = 0;

   // Access to the target class.
   // In the event the value class cannot be specified from the reference
   // itself, because the object behind the reference requires a cast,
   // the return value must be NULL.
   //
   // data   [IN]   Resolved pointer to the referenced object
   virtual TClass* GetValueClass(void* data) const = 0;

   // Update (and propagate) cached information
   virtual Bool_t Update() = 0;

   // Flag to indicate if this is a container reference
   virtual Bool_t HasCounter()  const = 0;

   // Access to container size (if container reference (ie TRefArray) etc)
   //
   // info    [IN]   Pointer to the structure called by TTree::Draw
   //                to extract the required object information.
   // data    [IN]   Pointer to the reference object
   //
   //  return value: The prepared pointer to the reference.
   virtual Int_t  GetCounterValue(TFormLeafInfoReference* info, void *data) = 0;

   // Prepare reused reference object (e.g. ZERO data pointers)
   // Because TTree::Draw reuses objects some reference implementations
   // require setup. For example the pointer to the object the reference points to
   // needs to be ZEROed.
   //
   // data   [IN]   Pointer to the reference object
   //
   // return value: The prepared pointer to the reference.
   virtual void* GetPreparedReference(void* data) = 0;

   // Access referenced object(-data)
   //
   // info     [IN]   Pointer to the structure called by TTree::Draw
   //                 to extract the required object information.
   // data     [IN]   Pointer to the referenced object
   // instance [IN]   Item number if ref collection
   //
   // return value: Pointer to the requested information
   virtual void* GetObject(TFormLeafInfoReference* info, void* data, int instance) = 0;
};
#endif // ROOT_TVirtualRefProxy
