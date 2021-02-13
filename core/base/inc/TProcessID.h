// @(#)root/cont:$Id$
// Author: Rene Brun   28/09/2001

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TProcessID
#define ROOT_TProcessID


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TProcessID                                                           //
//                                                                      //
// Process Identifier object                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TNamed.h"

#ifdef R__LESS_INCLUDES
class TObjArray;
#else
#include "TObjArray.h"
#endif

#include <atomic>
#include <type_traits>

class TExMap;

namespace ROOT {
   namespace Internal {
     /**
      * \class ROOT::Internal::TAtomicPointer
      * \brief Helper class to manage atomic pointers.
      * \tparam T Pointer type to be made atomic
      *
      * Helper class to manage atomic pointers. The class enforces that the templated type
      * is a pointer.
      */
      template <typename T> class TAtomicPointer {
         private:
            std::atomic<T> fAtomic;

         public:
            TAtomicPointer() : fAtomic(nullptr)
            {
               static_assert(std::is_pointer<T>::value, "Only pointer types supported");
            }

            ~TAtomicPointer() { delete fAtomic.load(); }

            T operator->() const { return fAtomic; }

            operator T() const { return fAtomic; }

            T operator=(const T& t)
            {
               fAtomic = t;
               return t;
            }
      };
   } // End of namespace Internal
} // End of namespace ROOT


class TProcessID : public TNamed {

private:
   TProcessID(const TProcessID &ref);            // TProcessID are not copiable.
   TProcessID& operator=(const TProcessID &ref); // TProcessID are not copiable.

protected:
   std::atomic_int    fCount;                           //!Reference count to this object (from TFile)
   ROOT::Internal::TAtomicPointer<TObjArray*> fObjects; //!Array pointing to the referenced objects
   std::atomic_flag   fLock;                            //!Spin lock for initialization of fObjects

   static TProcessID *fgPID;      //Pointer to current session ProcessID
   static TObjArray  *fgPIDs;     //Table of ProcessIDs
   static TExMap     *fgObjPIDs;  //Table pointer to pids

   static std::atomic_uint      fgNumber;   //Referenced objects count

public:
   TProcessID();
   virtual ~TProcessID();
   void             CheckInit();
   virtual void     Clear(Option_t *option="");
   Int_t            DecrementCount();
   Int_t            IncrementCount();
   Int_t            GetCount() const {return fCount;}
   TObjArray       *GetObjects() const {return fObjects;}
   TObject         *GetObjectWithID(UInt_t uid);
   void             PutObjectWithID(TObject *obj, UInt_t uid=0);
   virtual void     RecursiveRemove(TObject *obj);

   static TProcessID  *AddProcessID();
   static UInt_t       AssignID(TObject *obj);
   static void         Cleanup();
   static UInt_t       GetNProcessIDs();
   static TProcessID  *GetPID();
   static TObjArray   *GetPIDs();
   static TProcessID  *GetProcessID(UShort_t pid);
   static TProcessID  *GetProcessWithUID(const TObject *obj);
   static TProcessID  *GetProcessWithUID(UInt_t uid,const void *obj);
   static TProcessID  *GetSessionProcessID();
   static  UInt_t      GetObjectCount();
   static  Bool_t      IsValid(TProcessID *pid);
   static  void        SetObjectCount(UInt_t number);

   ClassDef(TProcessID,1)  //Process Unique Identifier in time and space
};

#endif
