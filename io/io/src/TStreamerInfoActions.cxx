// @(#)root/io:$Id$
// Author: Philippe Canal 05/2010

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TStreamerInfo.h"
#include "TStreamerInfoActions.h"
#include "TROOT.h"
#include "TStreamerElement.h"
#include "TVirtualMutex.h"
#include "TInterpreter.h"
#include "TError.h"
#include "TVirtualArray.h"
#include "TBufferFile.h"
#include "TBufferText.h"
#include "TMemberStreamer.h"
#include "TClassEdit.h"
#include "TVirtualCollectionIterators.h"
#include "TProcessID.h"
#include "TFile.h"

static const Int_t kRegrouped = TStreamerInfo::kOffsetL;

// More possible optimizations:
// Avoid call the virtual version of TBuffer::ReadInt and co.
// Merge the Reading of the version and the looking up or the StreamerInfo
// Avoid if (bytecnt) inside the CheckByteCount routines and avoid multiple (mostly useless nested calls)
// Try to avoid if statement on onfile class being set (TBufferFile::ReadClassBuffer).

using namespace TStreamerInfoActions;

#ifdef _AIX
# define INLINE_TEMPLATE_ARGS
#else
# define INLINE_TEMPLATE_ARGS inline
#endif


namespace TStreamerInfoActions
{
   bool IsDefaultVector(TVirtualCollectionProxy &proxy)
   {
      const auto props = proxy.GetProperties();
      const bool isVector = proxy.GetCollectionType() == ROOT::kSTLvector;
      const bool hasDefaultAlloc = !(props & TVirtualCollectionProxy::kCustomAlloc);
      const bool isEmulated = props & TVirtualCollectionProxy::kIsEmulated;

      return isEmulated || (isVector && hasDefaultAlloc);
   }

   template <typename From>
   struct WithFactorMarker {
      typedef From Value_t;
   };

   template <typename From>
   struct NoFactorMarker {
      typedef From Value_t;
   };

   struct BitsMarker {
      typedef UInt_t Value_t;
   };

   void TConfiguration::AddToOffset(Int_t delta)
   {
      // Add the (potentially negative) delta to all the configuration's offset.  This is used by
      // TBranchElement in the case of split sub-object.

      if (fOffset != TVirtualStreamerInfo::kMissing)
         fOffset += delta;
   }

  void TConfiguration::SetMissing()
   {
      // Add the (potentially negative) delta to all the configuration's offset.  This is used by
      // TBranchElement in the case of split sub-object.

      fOffset = TVirtualStreamerInfo::kMissing;
   }

   void TConfiguredAction::PrintDebug(TBuffer &buf, void *addr) const
   {
      // Inform the user what we are about to stream.

      // Idea, we should find a way to print the name of the function
      if (fConfiguration) fConfiguration->PrintDebug(buf,addr);
   }

   void TConfiguration::Print() const
   {
      // Inform the user what we are about to stream.

      TStreamerInfo *info = (TStreamerInfo*)fInfo;
      TStreamerElement *aElement = fCompInfo->fElem;
      TString sequenceType;
      aElement->GetSequenceType(sequenceType);

      printf("StreamerInfoAction, class:%s, name=%s, fType[%d]=%d,"
             " %s, offset=%d (%s)\n",
             info->GetClass()->GetName(), aElement->GetName(), fElemId, fCompInfo->fType,
             aElement->ClassName(), fOffset, sequenceType.Data());
   }

   void TConfiguration::PrintDebug(TBuffer &buf, void *addr) const
   {
      // Inform the user what we are about to stream.

      if (gDebug > 1) {
         // Idea: We should print the name of the action function.
         TStreamerInfo *info = (TStreamerInfo*)fInfo;
         TStreamerElement *aElement = fCompInfo->fElem;
         TString sequenceType;
         aElement->GetSequenceType(sequenceType);

         printf("StreamerInfoAction, class:%s, name=%s, fType[%d]=%d,"
                " %s, bufpos=%d, arr=%p, offset=%d (%s)\n",
                info->GetClass()->GetName(), aElement->GetName(), fElemId, fCompInfo->fType,
                aElement->ClassName(), buf.Length(), addr, fOffset, sequenceType.Data());
      }
   }

   void TLoopConfiguration::Print() const
   {
      // Inform the user what we are about to stream.

      printf("TLoopConfiguration: unconfigured\n");
   }


   struct TGenericConfiguration : TConfiguration {
      // Configuration of action using the legacy code.
      // Mostly to cancel out the PrintDebug.
   public:
      TGenericConfiguration(TVirtualStreamerInfo *info, UInt_t id, TCompInfo_t *compinfo, Int_t offset = 0) : TConfiguration(info,id,compinfo,offset) {};
      void PrintDebug(TBuffer &, void *) const {
         // Since we call the old code, it will print the debug statement.
      }

      virtual TConfiguration *Copy() { return new TGenericConfiguration(*this); }
  };

   struct TBitsConfiguration : TConfiguration {
      // Configuration of action handling kBits.
      // In this case we need to know both the location
      // of the member (fBits) and the start of the object
      // (its TObject part to be exact).

      Int_t  fObjectOffset;  // Offset of the TObject part within the object

      TBitsConfiguration(TVirtualStreamerInfo *info, UInt_t id, TCompInfo_t *compinfo, Int_t offset = 0) : TConfiguration(info,id,compinfo,offset),fObjectOffset(0) {};
      void PrintDebug(TBuffer &, void *) const {
         TStreamerInfo *info = (TStreamerInfo*)fInfo;
         TStreamerElement *aElement = fCompInfo->fElem;
         TString sequenceType;
         aElement->GetSequenceType(sequenceType);

         printf("StreamerInfoAction, class:%s, name=%s, fType[%d]=%d,"
                " %s, offset=%d (%s)\n",
                info->GetClass()->GetName(), aElement->GetName(), fElemId, fCompInfo->fType,
                aElement->ClassName(), fOffset, sequenceType.Data());
      }

      void AddToOffset(Int_t delta)
      {
         // Add the (potentially negative) delta to all the configuration's offset.  This is used by
         // TBranchElement in the case of split sub-object.

         if (fOffset != TVirtualStreamerInfo::kMissing)
            fOffset += delta;
         fObjectOffset = 0;
      }

      void SetMissing()
      {
         fOffset = TVirtualStreamerInfo::kMissing;
         fObjectOffset = 0;
      }

      virtual TConfiguration *Copy() { return new TBitsConfiguration(*this); }

   };

   Int_t GenericReadAction(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      char *obj = (char*)addr;
      TGenericConfiguration *conf = (TGenericConfiguration*)config;
      return ((TStreamerInfo*)conf->fInfo)->ReadBuffer(buf, &obj, &(conf->fCompInfo), /*first*/ 0, /*last*/ 1, /*narr*/ 1, config->fOffset, 2);
   }

   Int_t GenericWriteAction(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      char *obj = (char*)addr;
      TGenericConfiguration *conf = (TGenericConfiguration*)config;
      return ((TStreamerInfo*)conf->fInfo)->WriteBufferAux(buf, &obj, &(conf->fCompInfo), /*first*/ 0, /*last*/ 1, /*narr*/ 1, config->fOffset, 2);
   }

   template <typename T>
   INLINE_TEMPLATE_ARGS Int_t ReadBasicType(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      T *x = (T*)( ((char*)addr) + config->fOffset );
      // Idea: Implement buf.ReadBasic/Primitive to avoid the return value
      buf >> *x;
      return 0;
   }

   void HandleReferencedTObject(TBuffer &buf, void *addr, const TConfiguration *config) {
      TBitsConfiguration *conf = (TBitsConfiguration*)config;
      UShort_t pidf;
      buf >> pidf;
      pidf += buf.GetPidOffset();
      TProcessID *pid = buf.ReadProcessID(pidf);
      if (pid!=0) {
         TObject *obj = (TObject*)( ((char*)addr) + conf->fObjectOffset);
         UInt_t gpid = pid->GetUniqueID();
         UInt_t uid;
         if (gpid>=0xff) {
            uid = obj->GetUniqueID() | 0xff000000;
         } else {
            uid = ( obj->GetUniqueID() & 0xffffff) + (gpid<<24);
         }
         obj->SetUniqueID(uid);
         pid->PutObjectWithID(obj);
      }
   }

   template <>
   INLINE_TEMPLATE_ARGS Int_t ReadBasicType<BitsMarker>(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      UInt_t *x = (UInt_t*)( ((char*)addr) + config->fOffset );
      // Idea: Implement buf.ReadBasic/Primitive to avoid the return value
      // Idea: This code really belongs inside TBuffer[File]
      buf >> *x;

      if ((*x & kIsReferenced) != 0) {
         HandleReferencedTObject(buf,addr,config);
      }
      return 0;
   }

   template <typename T>
   INLINE_TEMPLATE_ARGS Int_t WriteBasicType(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      T *x = (T *)(((char *)addr) + config->fOffset);
      // Idea: Implement buf.ReadBasic/Primitive to avoid the return value
      buf << *x;
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t WriteTextTNamed(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      void *x = (void *)(((char *)addr) + config->fOffset);
      // Idea: Implement buf.ReadBasic/Primitive to avoid the return value
      buf.StreamObject(x, TNamed::Class(), TNamed::Class());
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t WriteTextTObject(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      void *x = (void *)(((char *)addr) + config->fOffset);
      // Idea: Implement buf.ReadBasic/Primitive to avoid the return value
      buf.StreamObject(x, TObject::Class(), TObject::Class());
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t WriteTextBaseClass(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      void *x = (void *)(((char *)addr) + config->fOffset);
      // Idea: Implement buf.ReadBasic/Primitive to avoid the return value
      ((TBufferText *)&buf)->WriteBaseClass(x, (TStreamerBase *)config->fCompInfo->fElem);
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t WriteTextStreamer(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      void *x = (void *)(((char *)addr) + config->fOffset);
      TMemberStreamer *pstreamer = config->fCompInfo->fStreamer;
      UInt_t pos = buf.WriteVersion(config->fInfo->IsA(), kTRUE);
      (*pstreamer)(buf, x, config->fCompInfo->fLength);
      buf.SetByteCount(pos, kTRUE);
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t ReadTextObject(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      void *x = (void *)(((char *)addr) + config->fOffset);
      buf.ReadFastArray(x, config->fCompInfo->fClass, config->fCompInfo->fLength, config->fCompInfo->fStreamer);
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t ReadTextTObject(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      void *x = (void *)(((char *)addr) + config->fOffset);
      buf.StreamObject(x, TObject::Class(), TObject::Class());
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t ReadTextBaseClass(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      void *x = (void *)(((char *)addr) + config->fOffset);
      // Idea: Implement buf.ReadBasic/Primitive to avoid the return value
      ((TBufferText *)&buf)->ReadBaseClass(x, (TStreamerBase *)config->fCompInfo->fElem);
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t ReadTextStreamer(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      void *x = (void *)(((char *)addr) + config->fOffset);
      TMemberStreamer *pstreamer = config->fCompInfo->fStreamer;

      UInt_t start, count;
      /* Version_t v = */ buf.ReadVersion(&start, &count, config->fCompInfo->fClass);
      (*pstreamer)(buf, x, config->fCompInfo->fLength);
      buf.CheckByteCount(start, count, config->fCompInfo->fElem->GetFullName());
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t ReadTextTObjectBase(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      // action required to call custom code for TObject as base class
      void *x = (void *)(((char *)addr) + config->fOffset);
      buf.ReadClassBuffer(TObject::Class(), x, TObject::Class());
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t ReadTextTNamed(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      void *x = (void *)(((char *)addr) + config->fOffset);
      buf.StreamObject(x, TNamed::Class(), TNamed::Class());
      return 0;
   }

   /** Direct copy of code from TStreamerInfo::WriteBufferAux,
    * potentially can be used later for non-text streaming */
   template<bool kIsTextT>
   INLINE_TEMPLATE_ARGS Int_t WriteSTLp(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      TClass *cl                 = config->fCompInfo->fClass;
      TMemberStreamer *pstreamer = config->fCompInfo->fStreamer;
      TVirtualCollectionProxy *proxy = cl->GetCollectionProxy();
      TClass* vClass = proxy ? proxy->GetValueClass() : 0;
      UInt_t eoffset = 0; // extra parameter of TStreamerInfo::WriteBufferAux, 0 for all kind of objects writing
      UInt_t ioffset = eoffset + config->fOffset;

      if (!buf.TestBit(TBuffer::kCannotHandleMemberWiseStreaming)
          && proxy && vClass
          && config->fInfo->GetStreamMemberWise()
          && cl->CanSplit()
          && !(strspn(config->fCompInfo->fElem->GetTitle(),"||") == 2)
          && !(vClass->TestBit(TClass::kHasCustomStreamerMember)) ) {
         // Let's save the collection member-wise.

         UInt_t pos = buf.WriteVersionMemberWise(config->fInfo->IsA(),kTRUE);
         buf.WriteVersion( vClass, kFALSE );

         // TODO: subinfo used for WriteBufferSTL call, which is private for the moment
         //TStreamerInfo *subinfo = (TStreamerInfo*)vClass->GetStreamerInfo();

         //for (int k = 0; k < narr; ++k) {
            char **contp = (char **)((char *)addr + ioffset);
            for(int j=0;j<config->fCompInfo->fLength;++j) {
               char *cont = contp[j];
               TVirtualCollectionProxy::TPushPop helper( proxy, cont );
               Int_t nobjects = cont ? proxy->Size() : 0;
               buf << nobjects;

               // TODO: method is private, should be made accesible from here
               // subinfo->WriteBufferSTL(buf,proxy,nobjects);
            }
         //}
         buf.SetByteCount(pos,kTRUE);
         return 0;
      }
      UInt_t pos = buf.WriteVersion(config->fInfo->IsA(), kTRUE);
      if (kIsTextT) {
         // use same method which is used in kSTL
         buf.WriteFastArray((void **)((char *)addr + ioffset), cl, config->fCompInfo->fLength, kFALSE, pstreamer);
      } else if (pstreamer == nullptr) {
         // for (int k = 0; k < narr; ++k) {
         char **contp = (char **)((char *)addr + ioffset);
         for (int j = 0; j < config->fCompInfo->fLength; ++j) {
            char *cont = contp[j];
            cl->Streamer(cont, buf);
         }
         // }
      } else {
         // for (int k = 0; k < narr; ++k) {
         (*pstreamer)(buf, (char *)addr + ioffset, config->fCompInfo->fLength);
         //}
      }
      buf.SetByteCount(pos, kTRUE);
      return 0;
   }


   /** Direct copy of code from TStreamerInfo::WriteBufferAux,
    * potentially can be used later for non-text streaming */
   template<bool kIsTextT>
   INLINE_TEMPLATE_ARGS Int_t ReadSTLp(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      TClass *cle                = config->fCompInfo->fClass;
      TStreamerElement * aElement  = (TStreamerElement*) config->fCompInfo->fElem;
      TMemberStreamer *pstreamer = config->fCompInfo->fStreamer;
      //TVirtualCollectionProxy *proxy = cl->GetCollectionProxy();
      //TClass* vClass = proxy ? proxy->GetValueClass() : 0;

      UInt_t eoffset = 0; // extra parameter of TStreamerInfo::WriteBufferAux, 0 for all kind of objects writing
      UInt_t ioffset = eoffset + config->fOffset;


      UInt_t start,count;
      Version_t vers = buf.ReadVersion(&start, &count, cle);

      if (!kIsTextT && (vers & TBufferFile::kStreamedMemberWise) ) {
         // Collection was saved member-wise

         vers &= ~( TBufferFile::kStreamedMemberWise );

         TClass *newClass = aElement->GetNewClass();
         TClass *oldClass = aElement->GetClassPointer();
         if( vers < 9 && newClass && newClass!=oldClass ) {
            Error( "ReadBuffer", "Unfortunately, version %d of TStreamerInfo (used in %s) did not record enough information to convert a %s into a %s.",
                  vers, buf.GetParent() ? buf.GetParent()->GetName() : "memory/socket", oldClass->GetName(), newClass->GetName() );
            return 0;
         }

         Version_t vClVersion = 0; // For vers less than 9, we have to use the current version.
         if( vers >= 9 ) {
            vClVersion = buf.ReadVersionForMemberWise( cle->GetCollectionProxy()->GetValueClass() );
         }

         TVirtualCollectionProxy *newProxy = (newClass ? newClass->GetCollectionProxy() : nullptr);
         TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
         TStreamerInfo *subinfo = nullptr;

         if( newProxy ) {
            // coverity[dereference] oldProxy->GetValueClass() can not be null since this was streamed memberwise.
            subinfo = (TStreamerInfo*)newProxy->GetValueClass()->GetConversionStreamerInfo( oldProxy->GetValueClass(), vClVersion );
         } else {
            subinfo = (TStreamerInfo*)oldProxy->GetValueClass()->GetStreamerInfo( vClVersion );
            newProxy = oldProxy;
         }
         if (subinfo) {
            // DOLOOP {
               void* env;
               void **contp = (void**)((char *) addr + ioffset);
               for(int j=0;j<config->fCompInfo->fLength;j++) {
                  void *cont = contp[j];
                  if (cont==nullptr) {
                     contp[j] = cle->New();
                     cont = contp[j];
                  }
                  TVirtualCollectionProxy::TPushPop helper( newProxy, cont );
                  Int_t nobjects;
                  buf >> nobjects;
                  env = newProxy->Allocate(nobjects,true);
                  subinfo->ReadBufferSTL(buf,newProxy,nobjects,/* offset */ 0, vers>=7 );
                  newProxy->Commit(env);
               }
           // } // DOLOOP
         }
         buf.CheckByteCount(start,count,aElement->GetFullName());
         return 0;
      }

      if (kIsTextT) {
         // use same method which is used in kSTL
         buf.ReadFastArray((void **)((char *)addr + ioffset), cle, config->fCompInfo->fLength, kFALSE, pstreamer);
      } else if (pstreamer == nullptr) {
         // DOLOOP {
         void **contp = (void **)((char *)addr + ioffset);
         for (int j = 0; j < config->fCompInfo->fLength; j++) {
            void *cont = contp[j];
            if (cont == nullptr) {
               // int R__n;
               // b >> R__n;
               // b.SetOffset(b.GetOffset()-4); // rewind to the start of the int
               // if (R__n) continue;
               contp[j] = cle->New();
               cont = contp[j];
            }
            cle->Streamer(cont, buf);
         }
         // }
      } else {
         (*pstreamer)(buf, (char *)addr + ioffset, config->fCompInfo->fLength);
      }
      buf.CheckByteCount(start, count, aElement->GetFullName());

      return 0;
   }

   /** Direct copy of code from TStreamerInfo::WriteBufferAux,
    * potentially can be used later for non-text streaming */
   template<bool kIsTextT>
   INLINE_TEMPLATE_ARGS Int_t WriteStreamerLoop(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      UInt_t eoffset = 0; // extra parameter of TStreamerInfo::WriteBufferAux, 0 for all kind of objects writing
      UInt_t ioffset = eoffset + config->fOffset;

      if (!kIsTextT && config->fCompInfo->fStreamer) {
         // Get any private streamer which was set for the data member.
         TMemberStreamer* pstreamer = config->fCompInfo->fStreamer;
         // -- We have a private streamer.
         UInt_t pos = buf.WriteVersion(config->fInfo->IsA(), kTRUE);
         // Loop over the entries in the clones array or the STL container.
         //for (int k = 0; k < narr; ++k) {
            // Get a pointer to the counter for the varying length array.
            Int_t* counter = (Int_t*) ((char *) addr /*entry pointer*/ + eoffset /*entry offset*/ + config->fCompInfo->fMethod /*counter offset*/);

            // And call the private streamer, passing it the buffer, the object, and the counter.
            (*pstreamer)(buf, (char *) addr /*entry pointer*/ + ioffset /*object offset*/, *counter);
         //} for k
         buf.SetByteCount(pos, kTRUE);
         // We are done, next streamer element.
         return 0;
      }

      // Get the class of the data member.
      TClass* cl = config->fCompInfo->fClass;
      // Which are we, an array of objects or an array of pointers to objects?
      Bool_t isPtrPtr = (strstr(config->fCompInfo->fElem->GetTypeName(), "**") != 0);

      // By default assume the file version is the newest.
      Int_t fileVersion = kMaxInt;

      if (!kIsTextT) {
         // At this point we do *not* have a private streamer.
         // Get the version of the file we are writing to.
         TFile* file = (TFile*) buf.GetParent();
         if (file) {
            fileVersion = file->GetVersion();
         }
      }
      // Write the class version to the buffer.
      UInt_t pos = buf.WriteVersion(config->fInfo->IsA(), kTRUE);
      if (fileVersion > 51508) {
         // -- Newer versions allow polymorphic pointers to objects.
         // Loop over the entries in the clones array or the STL container.
         //for (int k = 0; k < narr; ++k) {
            // Get the counter for the varying length array.
            Int_t vlen = *((Int_t*) ((char *) addr /*entry pointer*/ + eoffset /*entry offset*/ + config->fCompInfo->fMethod /*counter offset*/));

            //b << vlen;
            if (vlen) {
               // Get a pointer to the array of pointers.
               char** pp = (char**) ((char *) addr /*entry pointer*/ + ioffset /*object offset*/);
               // Loop over each element of the array of pointers to varying-length arrays.
               for (Int_t ndx = 0; ndx < config->fCompInfo->fLength; ++ndx) {
                  if (!pp[ndx]) {
                     // -- We do not have a pointer to a varying-length array.
                     // Error("WriteBufferAux", "The pointer to element %s::%s type %d (%s) is null\n", GetName(), aElement->GetFullName(), compinfo[i]->fType, aElement->GetTypeName());
                     // ::ErrorHandler(kError, "::WriteStreamerLoop", Form("The pointer to element %s::%s type %d (%s) is null\n", config->fInfo->GetName(), config->fCompInfo->fElem->GetFullName(), config->fCompInfo->fType, config->fCompInfo->fElem->GetTypeName()));
                     printf("WriteStreamerLoop - The pointer to element %s::%s type %d (%s) is null\n", config->fInfo->GetName(), config->fCompInfo->fElem->GetFullName(), config->fCompInfo->fType, config->fCompInfo->fElem->GetTypeName());
                     continue;
                  }
                  if (!isPtrPtr) {
                     // -- We are a varying-length array of objects.
                     // Write the entire array of objects to the buffer.
                     // Note: Polymorphism is not allowed here.
                     buf.WriteFastArray(pp[ndx], cl, vlen, nullptr);
                  } else {
                     // -- We are a varying-length array of pointers to objects.
                     // Write the entire array of object pointers to the buffer.
                     // Note: The object pointers are allowed to be polymorphic.
                     buf.WriteFastArray((void **)pp[ndx], cl, vlen, kFALSE, nullptr);
                  } // isPtrPtr
               } // ndx
            } else // vlen
            if (kIsTextT) {
               // special handling for the text-based streamers
               for (Int_t ndx = 0; ndx < config->fCompInfo->fLength; ++ndx)
                  buf.WriteFastArray((void *)nullptr, cl, -1, nullptr);
            }
            //} // k
      }
      else {
         // -- Older versions do *not* allow polymorphic pointers to objects.
         // Loop over the entries in the clones array or the STL container.
         //for (int k = 0; k < narr; ++k) {
            // Get the counter for the varying length array.
            Int_t vlen = *((Int_t*) ((char *) addr /*entry pointer*/ + eoffset /*entry offset*/ + config->fCompInfo->fMethod /*counter offset*/));
            //b << vlen;
            if (vlen) {
               // Get a pointer to the array of pointers.
               char** pp = (char**) ((char *) addr /*entry pointer*/ + ioffset /*object offset*/);
               // -- Older versions do *not* allow polymorphic pointers to objects.
               // Loop over each element of the array of pointers to varying-length arrays.
               for (Int_t ndx = 0; ndx < config->fCompInfo->fLength; ++ndx) {
                  if (!pp[ndx]) {
                     // -- We do not have a pointer to a varying-length array.
                     //Error("WriteBufferAux", "The pointer to element %s::%s type %d (%s) is null\n", GetName(), aElement->GetFullName(), compinfo[i]->fType, aElement->GetTypeName());
                     // ::ErrorHandler(kError, "::WriteTextStreamerLoop", Form("The pointer to element %s::%s type %d (%s) is null\n", config->fInfo->GetName(), config->fCompInfo->fElem->GetFullName(), config->fCompInfo->fType, config->fCompInfo->fElem->GetTypeName()));
                     printf("WriteStreamerLoop - The pointer to element %s::%s type %d (%s) is null\n", config->fInfo->GetName(), config->fCompInfo->fElem->GetFullName(), config->fCompInfo->fType, config->fCompInfo->fElem->GetTypeName());
                     continue;
                  }
                  if (!isPtrPtr) {
                     // -- We are a varying-length array of objects.
                     // Loop over the elements of the varying length array.
                     for (Int_t v = 0; v < vlen; ++v) {
                        // Write the object to the buffer.
                        cl->Streamer(pp[ndx] + (v * cl->Size()), buf);
                     } // v
                  }
                  else {
                     // -- We are a varying-length array of pointers to objects.
                     // Loop over the elements of the varying length array.
                     for (Int_t v = 0; v < vlen; ++v) {
                        // Get a pointer to the object pointer.
                        char** r = (char**) pp[ndx];
                        // Write the object to the buffer.
                        cl->Streamer(r[v], buf);
                     } // v
                  } // isPtrPtr
               } // ndx
            } // vlen
         //} // k
      } // fileVersion
      // Backpatch the byte count into the buffer.
      buf.SetByteCount(pos, kTRUE);

      return 0;
   }


   /** Direct copy of code from TStreamerInfo::WriteBufferAux,
    * potentially can be used later for non-text streaming */
   template<bool kIsTextT>
   INLINE_TEMPLATE_ARGS Int_t ReadStreamerLoop(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      UInt_t eoffset = 0; // extra parameter of TStreamerInfo::WriteBufferAux, 0 for all kind of objects writing
      UInt_t ioffset = eoffset + config->fOffset;

      // Get the class of the data member.
      TClass* cl = config->fCompInfo->fClass;

      // Check for a private streamer.
      if (!kIsTextT && config->fCompInfo->fStreamer) {
         // Get any private streamer which was set for the data member.
         TMemberStreamer* pstreamer = config->fCompInfo->fStreamer;
         // -- We have a private streamer.
         // Read the class version and byte count from the buffer.
         UInt_t start = 0;
         UInt_t count = 0;
         buf.ReadVersion(&start, &count, cl);
         // Loop over the entries in the clones array or the STL container.
         //for (Int_t k = 0; k < narr; ++k) {

            Int_t* counter = (Int_t*) ((char *) addr /*entry pointer*/ + eoffset /*entry offset*/ + config->fCompInfo->fMethod /*counter offset*/);
            // And call the private streamer, passing it the buffer, the object, and the counter.
            (*pstreamer)(buf, (char *) addr /*entry pointer*/ + ioffset /*object offset*/, *counter);

         // } // for k
         buf.CheckByteCount(start, count, config->fCompInfo->fElem->GetFullName());
         // We are done, next streamer element.
         return 0;
      }

      // Which are we, an array of objects or an array of pointers to objects?
      Bool_t isPtrPtr = (strstr(config->fCompInfo->fElem->GetTypeName(), "**") != 0);

      // By default assume the file version is the newest.
      Int_t fileVersion = kMaxInt;
      if (!kIsTextT) {
         // At this point we do *not* have a private streamer.
         // Get the version of the file we are reading from.
         TFile* file = (TFile*) buf.GetParent();
         if (file) {
            fileVersion = file->GetVersion();
         }
      }
      // Read the class version and byte count from the buffer.
      UInt_t start = 0;
      UInt_t count = 0;
      buf.ReadVersion(&start, &count, cl);
      if (fileVersion > 51508) {
         // -- Newer versions allow polymorphic pointers.
         // Loop over the entries in the clones array or the STL container.
         // for (Int_t k = 0; k < narr; ++k) {
         // Get the counter for the varying length array.
         Int_t vlen = *((Int_t *)((char *)addr /*entry pointer*/ + eoffset /*entry offset*/ +
                                  config->fCompInfo->fMethod /*counter offset*/));
         // Int_t realLen;
         // b >> realLen;
         // if (realLen != vlen) {
         //   fprintf(stderr, "read vlen: %d  realLen: %s\n", vlen, realLen);
         //}
         // Get a pointer to the array of pointers.
         char **pp = (char **)((char *)addr /*entry pointer*/ + ioffset /*object offset*/);
         // Loop over each element of the array of pointers to varying-length arrays.
         // if (!pp) {
         //   continue;
         // }

         if (pp) // SL: place it here instead of continue, which is related to for(k) loop
            for (Int_t ndx = 0; ndx < config->fCompInfo->fLength; ++ndx) {
               // if (!pp[ndx]) {
               // -- We do not have a pointer to a varying-length array.
               // Error("ReadBuffer", "The pointer to element %s::%s type %d (%s) is null\n", thisVar->GetName(),
               // aElement->GetFullName(), compinfo[i]->fType, aElement->GetTypeName());
               // continue;
               //}
               // Delete any memory at pp[ndx].
               if (!isPtrPtr) {
                  cl->DeleteArray(pp[ndx]);
                  pp[ndx] = 0;
               } else {
                  // Using vlen is wrong here because it has already
                  // been overwritten with the value needed to read
                  // the current record.  Fixing this will require
                  // doing a pass over the object at the beginning
                  // of the I/O and releasing all the buffer memory
                  // for varying length arrays before we overwrite
                  // the counter values.
                  //
                  // For now we will just leak memory, just as we
                  // have always done in the past.  Fix this.
                  //
                  // char** r = (char**) pp[ndx];
                  // if (r) {
                  //   for (Int_t v = 0; v < vlen; ++v) {
                  //      cl->Destructor(r[v]);
                  //      r[v] = 0;
                  //   }
                  //}
                  delete[] pp[ndx];
                  pp[ndx] = 0;
               }
               if (!vlen) {
                  if (kIsTextT) {
                     // special handling for the text-based streamers - keep calling to shift array index
                     buf.ReadFastArray((void *)nullptr, cl, -1, nullptr);
                  }
                  continue;
               }
               // Note: We now have pp[ndx] is null.
               // Allocate memory to read into.
               if (!isPtrPtr) {
                  // -- We are a varying-length array of objects.
                  // Note: Polymorphism is not allowed here.
                  // Allocate a new array of objects to read into.
                  pp[ndx] = (char *)cl->NewArray(vlen);
                  if (!pp[ndx]) {
                     Error("ReadBuffer", "Memory allocation failed!\n");
                     continue;
                  }
               } else {
                  // -- We are a varying-length array of pointers to objects.
                  // Note: The object pointers are allowed to be polymorphic.
                  // Allocate a new array of pointers to objects to read into.
                  pp[ndx] = (char *)new char *[vlen];
                  if (!pp[ndx]) {
                     Error("ReadBuffer", "Memory allocation failed!\n");
                     continue;
                  }
                  // And set each pointer to null.
                  memset(pp[ndx], 0, vlen * sizeof(char *)); // This is the right size we really have a char**: pp[ndx]
                                                             // = (char*) new char*[vlen];
               }
               if (!isPtrPtr) {
                  // -- We are a varying-length array of objects.
                  buf.ReadFastArray(pp[ndx], cl, vlen, nullptr);
               } else {
                  // -- We are a varying-length array of object pointers.
                  buf.ReadFastArray((void **)pp[ndx], cl, vlen, kFALSE, nullptr);
               } // isPtrPtr
            }    // ndx
                 // } // k
      } else {
         // -- Older versions do *not* allow polymorphic pointers.
         // Loop over the entries in the clones array or the STL container.
         // for (Int_t k = 0; k < narr; ++k) {
         // Get the counter for the varying length array.
         Int_t vlen = *((Int_t *)((char *)addr /*entry pointer*/ + eoffset /*entry offset*/ +
                                  config->fCompInfo->fMethod /*counter offset*/));
         // Int_t realLen;
         // b >> realLen;
         // if (realLen != vlen) {
         //   fprintf(stderr, "read vlen: %d  realLen: %s\n", vlen, realLen);
         //}
         // Get a pointer to the array of pointers.
         char **pp = (char **)((char *)addr /*entry pointer*/ + ioffset /*object offset*/);
         // if (!pp) {
         //   continue;
         //}

         if (pp) // SL: place it here instead of continue, which is related to for(k) loop

            // Loop over each element of the array of pointers to varying-length arrays.
            for (Int_t ndx = 0; ndx < config->fCompInfo->fLength; ++ndx) {
               // if (!pp[ndx]) {
               // -- We do not have a pointer to a varying-length array.
               // Error("ReadBuffer", "The pointer to element %s::%s type %d (%s) is null\n", thisVar->GetName(),
               // aElement->GetFullName(), compinfo[i]->fType, aElement->GetTypeName());
               // continue;
               //}
               // Delete any memory at pp[ndx].
               if (!isPtrPtr) {
                  cl->DeleteArray(pp[ndx]);
                  pp[ndx] = 0;
               } else {
                  // Using vlen is wrong here because it has already
                  // been overwritten with the value needed to read
                  // the current record.  Fixing this will require
                  // doing a pass over the object at the beginning
                  // of the I/O and releasing all the buffer memory
                  // for varying length arrays before we overwrite
                  // the counter values.
                  //
                  // For now we will just leak memory, just as we
                  // have always done in the past.  Fix this.
                  //
                  // char** r = (char**) pp[ndx];
                  // if (r) {
                  //   for (Int_t v = 0; v < vlen; ++v) {
                  //      cl->Destructor(r[v]);
                  //      r[v] = 0;
                  //   }
                  //}
                  delete[] pp[ndx];
                  pp[ndx] = 0;
               }
               if (!vlen) {
                  continue;
               }
               // Note: We now have pp[ndx] is null.
               // Allocate memory to read into.
               if (!isPtrPtr) {
                  // -- We are a varying-length array of objects.
                  // Note: Polymorphism is not allowed here.
                  // Allocate a new array of objects to read into.
                  pp[ndx] = (char *)cl->NewArray(vlen);
                  if (!pp[ndx]) {
                     Error("ReadBuffer", "Memory allocation failed!\n");
                     continue;
                  }
               } else {
                  // -- We are a varying-length array of pointers to objects.
                  // Note: The object pointers are allowed to be polymorphic.
                  // Allocate a new array of pointers to objects to read into.
                  pp[ndx] = (char *)new char *[vlen];
                  if (!pp[ndx]) {
                     Error("ReadBuffer", "Memory allocation failed!\n");
                     continue;
                  }
                  // And set each pointer to null.
                  memset(pp[ndx], 0, vlen * sizeof(char *)); // This is the right size we really have a char**: pp[ndx]
                                                             // = (char*) new char*[vlen];
               }
               if (!isPtrPtr) {
                  // -- We are a varying-length array of objects.
                  // Loop over the elements of the varying length array.
                  for (Int_t v = 0; v < vlen; ++v) {
                     // Read the object from the buffer.
                     cl->Streamer(pp[ndx] + (v * cl->Size()), buf);
                  } // v
               } else {
                  // -- We are a varying-length array of object pointers.
                  // Get a pointer to the object pointer array.
                  char **r = (char **)pp[ndx];
                  // Loop over the elements of the varying length array.
                  for (Int_t v = 0; v < vlen; ++v) {
                     // Allocate an object to read into.
                     r[v] = (char *)cl->New();
                     if (!r[v]) {
                        // Do not print a second error message here.
                        // Error("ReadBuffer", "Memory allocation failed!\n");
                        continue;
                     }
                     // Read the object from the buffer.
                     cl->Streamer(r[v], buf);
                  } // v
               }    // isPtrPtr
            }       // ndx
         // } // k
      } // fileVersion
      buf.CheckByteCount(start, count, config->fCompInfo->fElem->GetFullName());
      return 0;
   }

   class TConfWithFactor : public TConfiguration {
      // Configuration object for the Float16/Double32 where a factor has been specified.
   public:
      Double_t fFactor;
      Double_t fXmin;
      TConfWithFactor(TVirtualStreamerInfo *info, UInt_t id, TCompInfo_t *compinfo, Int_t offset, Double_t factor, Double_t xmin) : TConfiguration(info,id,compinfo,offset),fFactor(factor),fXmin(xmin) {};
      virtual TConfiguration *Copy() { return new TConfWithFactor(*this); }
   };

   template <typename T>
   INLINE_TEMPLATE_ARGS Int_t ReadBasicType_WithFactor(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      // Stream a Float16 or Double32 where a factor has been specified.
      //a range was specified. We read an integer and convert it back to a double.

      TConfWithFactor *conf = (TConfWithFactor *)config;
      buf.ReadWithFactor((T*)( ((char*)addr) + config->fOffset ), conf->fFactor, conf->fXmin);
      return 0;
   }

   class TConfNoFactor : public TConfiguration {
      // Configuration object for the Float16/Double32 where a factor has been specified.
   public:
      Int_t fNbits;
      TConfNoFactor(TVirtualStreamerInfo *info, UInt_t id, TCompInfo_t *compinfo, Int_t offset, Int_t nbits) : TConfiguration(info,id,compinfo,offset),fNbits(nbits) {};
      virtual TConfiguration *Copy() { return new TConfNoFactor(*this); }
   };

   template <typename T>
   INLINE_TEMPLATE_ARGS Int_t ReadBasicType_NoFactor(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      // Stream a Float16 or Double32 where a factor has not been specified.

      TConfNoFactor *conf = (TConfNoFactor *)config;
      Int_t nbits = conf->fNbits;

      buf.ReadWithNbits( (T*)( ((char*)addr) + config->fOffset ), nbits );
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t ReadTString(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      // Read in a TString object.

      // Idea: We could separate the TString Streamer in its two parts and
      // avoid the if (buf.IsReading()) and try having it inlined.
      ((TString*)(((char*)addr)+config->fOffset))->TString::Streamer(buf);
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t ReadTObject(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      // Read in a TObject object part.

      // Idea: We could separate the TObject Streamer in its two parts and
      // avoid the if (buf.IsReading()).
      ((TObject*)(((char*)addr)+config->fOffset))->TObject::Streamer(buf);
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t ReadTNamed(TBuffer &buf, void *addr, const TConfiguration *config)
   {
      // Read in a TNamed object part.
      // Since the TNamed streamer is solely delegating back to the StreamerInfo we
      // can skip the streamer.

      // Idea: We could extract the code from ReadClassBuffer and avoid one function
      // code.
      static const TClass *TNamed_cl = TNamed::Class();
      return buf.ReadClassBuffer(TNamed_cl,(((char*)addr)+config->fOffset));
   }

   class TConfigSTL : public TConfiguration {
      // Configuration object for the kSTL case
   private:
      void Init() {
         TVirtualCollectionProxy *proxy = fNewClass->GetCollectionProxy();
         if (proxy) {
            fCreateIterators = proxy->GetFunctionCreateIterators();
            fCreateWriteIterators = proxy->GetFunctionCreateIterators(kFALSE);
            fCopyIterator = proxy->GetFunctionCopyIterator();
            fDeleteIterator = proxy->GetFunctionDeleteIterator();
            fDeleteTwoIterators = proxy->GetFunctionDeleteTwoIterators();
         }
      }

   public:
      TClass          *fOldClass;   // Class of the content on file
      TClass          *fNewClass;   // Class of the content in memory.
      TMemberStreamer *fStreamer;
      const char      *fTypeName;   // Type name of the member as typed by ther user.
      Bool_t          fIsSTLBase;  // aElement->IsBase() && aElement->IsA()!=TStreamerBase::Class()

      TVirtualCollectionProxy::CreateIterators_t    fCreateIterators;
      TVirtualCollectionProxy::CreateIterators_t    fCreateWriteIterators;
      TVirtualCollectionProxy::CopyIterator_t       fCopyIterator;
      TVirtualCollectionProxy::DeleteIterator_t     fDeleteIterator;
      TVirtualCollectionProxy::DeleteTwoIterators_t fDeleteTwoIterators;

      TConfigSTL(TVirtualStreamerInfo *info, UInt_t id, TCompInfo_t *compinfo, Int_t offset, UInt_t length, TClass *oldClass, const char *type_name, Bool_t isbase) :
         TConfiguration(info,id,compinfo,offset,length), fOldClass(oldClass), fNewClass(oldClass), fStreamer(0), fTypeName(type_name), fIsSTLBase(isbase),
         fCreateIterators(0), fCopyIterator(0), fDeleteIterator(0), fDeleteTwoIterators(0) { Init(); }

      TConfigSTL(TVirtualStreamerInfo *info, UInt_t id, TCompInfo_t *compinfo, Int_t offset, UInt_t length, TClass *oldClass, TClass *newClass, const char *type_name, Bool_t isbase) :
         TConfiguration(info,id,compinfo,offset,length), fOldClass(oldClass), fNewClass(newClass), fStreamer(0), fTypeName(type_name), fIsSTLBase(isbase),
         fCreateIterators(0), fCopyIterator(0), fDeleteIterator(0), fDeleteTwoIterators(0) { Init(); }

      TConfigSTL(TVirtualStreamerInfo *info, UInt_t id, TCompInfo_t *compinfo, Int_t offset, UInt_t length, TClass *oldClass, TMemberStreamer* streamer, const char *type_name, Bool_t isbase) :
         TConfiguration(info,id,compinfo,offset,length), fOldClass(oldClass), fNewClass(oldClass), fStreamer(streamer), fTypeName(type_name), fIsSTLBase(isbase),
         fCreateIterators(0), fCopyIterator(0), fDeleteIterator(0), fDeleteTwoIterators(0) { Init(); }

      TConfigSTL(TVirtualStreamerInfo *info, UInt_t id, TCompInfo_t *compinfo, Int_t offset, UInt_t length, TClass *oldClass, TClass *newClass, TMemberStreamer* streamer, const char *type_name, Bool_t isbase) :
         TConfiguration(info,id,compinfo,offset,length), fOldClass(oldClass), fNewClass(newClass), fStreamer(streamer), fTypeName(type_name), fIsSTLBase(isbase),
         fCreateIterators(0), fCopyIterator(0), fDeleteIterator(0), fDeleteTwoIterators(0) { Init(); }

      virtual TConfiguration *Copy() { return new TConfigSTL(*this); }
   };

   class TConfSTLWithFactor : public TConfigSTL {
      // Configuration object for the Float16/Double32 where a factor has been specified.
   public:
      Double_t fFactor;
      Double_t fXmin;
      TConfSTLWithFactor(TConfigSTL *orig, Double_t factor, Double_t xmin) : TConfigSTL(*orig),fFactor(factor),fXmin(xmin) {};
      virtual TConfiguration *Copy() { return new TConfSTLWithFactor(*this); }
   };

   class TConfSTLNoFactor : public TConfigSTL {
      // Configuration object for the Float16/Double32 where a factor has been specified.
   public:
      Int_t fNbits;
      TConfSTLNoFactor(TConfigSTL *orig, Int_t nbits) : TConfigSTL(*orig),fNbits(nbits) {};
      virtual TConfiguration *Copy() { return new TConfSTLNoFactor(*this); }
   };

   class TVectorLoopConfig : public TLoopConfiguration {
      // Base class of the Configurations used in member wise streaming.
   protected:
   public:
      Long_t fIncrement; // Either a value to increase the cursor by and
   public:
      TVectorLoopConfig(TVirtualCollectionProxy *proxy, Long_t increment, Bool_t /* read */) : TLoopConfiguration(proxy), fIncrement(increment) {};
      //virtual void PrintDebug(TBuffer &buffer, void *);
      virtual ~TVectorLoopConfig() {};
      void Print() const
      {
         printf("TVectorLoopConfig: increment=%ld\n",fIncrement);
      }

      void* GetFirstAddress(void *start, const void * /* end */) const
      {
         // Return the address of the first element of the collection.

         return start;
      }

      virtual TLoopConfiguration* Copy() const { return new TVectorLoopConfig(*this); }
   };

   class TAssocLoopConfig : public TLoopConfiguration {
      // Base class of the Configurations used in member wise streaming.
   public:
      TAssocLoopConfig(TVirtualCollectionProxy *proxy, Bool_t /* read */) : TLoopConfiguration(proxy) {};
      //virtual void PrintDebug(TBuffer &buffer, void *);
      virtual ~TAssocLoopConfig() {};
      void Print() const
      {
         printf("TAssocLoopConfig: proxy=%s\n",fProxy->GetCollectionClass()->GetName());
      }
      virtual TLoopConfiguration* Copy() const { return new TAssocLoopConfig(*this); }

      void* GetFirstAddress(void *start, const void * /* end */) const
      {
         // Return the address of the first element of the collection.

         R__ASSERT(0);
//         char iterator[TVirtualCollectionProxy::fgIteratorArenaSize];
//         void *iter = genloopconfig->fCopyIterator(&iterator,start_collection);
//         arr0 = genloopconfig->fNext(iter,end_collection);
//         if (iter != &iterator[0]) {
//            genloopconfig->fDeleteIterator(iter);
//         }
         return start;
      }
   };

   class TGenericLoopConfig : public TLoopConfiguration {
      // Configuration object for the generic case of member wise streaming looping.
   private:
      void Init(Bool_t read) {
         if (fProxy) {
            if (fProxy->HasPointers()) {
               fNext = TVirtualCollectionPtrIterators::Next;
               fCopyIterator = TVirtualCollectionPtrIterators::CopyIterator;
               fDeleteIterator = TVirtualCollectionPtrIterators::DeleteIterator;
            } else {
               fNext = fProxy->GetFunctionNext(read);
               fCopyIterator = fProxy->GetFunctionCopyIterator(read);
               fDeleteIterator = fProxy->GetFunctionDeleteIterator(read);
            }
         }
      }
   public:
      TVirtualCollectionProxy::Next_t               fNext;
      TVirtualCollectionProxy::CopyIterator_t       fCopyIterator;
      TVirtualCollectionProxy::DeleteIterator_t     fDeleteIterator;

      TGenericLoopConfig(TVirtualCollectionProxy *proxy, Bool_t read) : TLoopConfiguration(proxy), fNext(0), fCopyIterator(0), fDeleteIterator(0)
      {
         Init(read);
      }
      virtual ~TGenericLoopConfig() {};
      void Print() const
      {
         printf("TGenericLoopConfig: proxy=%s\n",fProxy->GetCollectionClass()->GetName());
      }
      virtual TLoopConfiguration* Copy() const { return new TGenericLoopConfig(*this); }

      void* GetFirstAddress(void *start_collection, const void *end_collection) const
      {
         // Return the address of the first element of the collection.

         char iterator[TVirtualCollectionProxy::fgIteratorArenaSize];
         void *iter = fCopyIterator(&iterator,start_collection);
         void *arr0 = fNext(iter,end_collection);
         if (iter != &iterator[0]) {
            fDeleteIterator(iter);
         }
         return arr0;
      }
   };

   INLINE_TEMPLATE_ARGS void ReadSTLMemberWiseSameClass(TBuffer &buf, void *addr, const TConfiguration *conf, Version_t vers)
   {
      // Collection was saved member-wise

      TConfigSTL *config = (TConfigSTL*)conf;
      vers &= ~( TBufferFile::kStreamedMemberWise );

      if( vers >= 8 ) {

         TClass *oldClass = config->fOldClass;

         TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
         if (!oldProxy) {
            // Missing information, broken file ... give up
            return;
         }
         TClass *valueClass = oldProxy->GetValueClass();
         Version_t vClVersion = buf.ReadVersionForMemberWise( valueClass );

         TVirtualCollectionProxy::TPushPop helper( oldProxy, (char*)addr );
         Int_t nobjects;
         buf.ReadInt(nobjects);
         void* alternative = oldProxy->Allocate(nobjects,true);
         if (nobjects) {
            TActionSequence *actions = oldProxy->GetReadMemberWiseActions( vClVersion );

            char startbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
            char endbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
            void *begin = &(startbuf[0]);
            void *end = &(endbuf[0]);
            config->fCreateIterators(alternative, &begin, &end, oldProxy);
            // We can not get here with a split vector of pointer, so we can indeed assume
            // that actions->fConfiguration != null.
            buf.ApplySequence(*actions, begin, end);
            if (begin != &(startbuf[0])) {
               // assert(end != endbuf);
               config->fDeleteTwoIterators(begin,end);
            }
         }
         oldProxy->Commit(alternative);

      } else {

         TClass *oldClass = config->fOldClass;

         TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
         if (!oldProxy) {
            // Missing information, broken file ... give up
            return;
         }

         TVirtualCollectionProxy::TPushPop helper( oldProxy, (char*)addr );
         Int_t nobjects;
         buf.ReadInt(nobjects);
         void* env = oldProxy->Allocate(nobjects,true);

         if (nobjects || vers < 7 ) {
            // coverity[dereference] since this is a member streaming action by definition the collection contains objects.
            TStreamerInfo *subinfo = (TStreamerInfo*)oldProxy->GetValueClass()->GetStreamerInfo( 0 );

            subinfo->ReadBufferSTL(buf, oldProxy, nobjects, /* offset */ 0, /* v7 */ kFALSE);
         }
         oldProxy->Commit(env);
      }
   }

   INLINE_TEMPLATE_ARGS void ReadArraySTLMemberWiseSameClass(TBuffer &buf, void *addr, const TConfiguration *conf, Version_t vers)
   {
      // Collection was saved member-wise

      TConfigSTL *config = (TConfigSTL*)conf;
      vers &= ~( TBufferFile::kStreamedMemberWise );

      if( vers >= 8 ) {

         TClass *oldClass = config->fOldClass;

         TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
         if (!oldProxy) {
            // Missing information, broken file ... give up
            return;
         }
         TClass *valueClass = oldProxy->GetValueClass();
         Version_t vClVersion = buf.ReadVersionForMemberWise( valueClass );

         TActionSequence *actions = oldProxy->GetReadMemberWiseActions( vClVersion );

         int objectSize = oldClass->Size();
         char *obj = (char*)addr;
         char *endobj = obj + conf->fLength*objectSize;

         for(; obj<endobj; obj+=objectSize) {
            Int_t nobjects;
            buf.ReadInt(nobjects);
            TVirtualCollectionProxy::TPushPop helper( oldProxy, (char*)obj );
            void* alternative = oldProxy->Allocate(nobjects,true);
            if (nobjects) {
               char startbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
               char endbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
               void *begin = &(startbuf[0]);
               void *end = &(endbuf[0]);
               config->fCreateIterators(alternative, &begin, &end, oldProxy);
               // We can not get here with a split vector of pointer, so we can indeed assume
               // that actions->fConfiguration != null.
               buf.ApplySequence(*actions, begin, end);
               if (begin != &(startbuf[0])) {
                  // assert(end != endbuf);
                  config->fDeleteTwoIterators(begin,end);
               }
            }
            oldProxy->Commit(alternative);
         }

      } else {

         TClass *oldClass = config->fOldClass;

         TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
         if (!oldProxy) {
            // Missing information, broken file ... give up
            return;
         }

         int objectSize = oldClass->Size();
         char *obj = (char*)addr;
         char *endobj = obj + conf->fLength*objectSize;

         for(; obj<endobj; obj+=objectSize) {
            TVirtualCollectionProxy::TPushPop helper( oldProxy, (char*)obj );
            Int_t nobjects;
            buf.ReadInt(nobjects);
            void* env = oldProxy->Allocate(nobjects,true);

            if (nobjects || vers < 7 ) {
               // coverity[dereference] since this is a member streaming action by definition the collection contains objects.
               TStreamerInfo *subinfo = (TStreamerInfo*)oldProxy->GetValueClass()->GetStreamerInfo( 0 );

               subinfo->ReadBufferSTL(buf, oldProxy, nobjects, /* offset */ 0, /* v7 */ kFALSE);
            }
            oldProxy->Commit(env);
         }
      }
   }

   INLINE_TEMPLATE_ARGS void ReadSTLMemberWiseChangedClass(TBuffer &buf, void *addr, const TConfiguration *conf, Version_t vers)
   {
      // Collection was saved member-wise

      TConfigSTL *config = (TConfigSTL*)conf;

      vers &= ~( TBufferFile::kStreamedMemberWise );

      TClass *newClass = config->fNewClass;
      TClass *oldClass = config->fOldClass;

      if( vers < 8 ) {
         Error( "ReadSTLMemberWiseChangedClass", "Unfortunately, version %d of TStreamerInfo (used in %s) did not record enough information to convert a %s into a %s.",
               vers, buf.GetParent() ? buf.GetParent()->GetName() : "memory/socket", oldClass->GetName(), newClass->GetName() );
      } else {

         Version_t vClVersion = buf.ReadVersionForMemberWise( oldClass->GetCollectionProxy()->GetValueClass() );

         TVirtualCollectionProxy *newProxy = newClass->GetCollectionProxy();
         TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();

         TVirtualCollectionProxy::TPushPop helper( newProxy, (char*)addr );
         Int_t nobjects;
         buf.ReadInt(nobjects);
         void* alternative = newProxy->Allocate(nobjects,true);
         if (nobjects) {
            TActionSequence *actions = newProxy->GetConversionReadMemberWiseActions( oldProxy->GetValueClass(), vClVersion );
            char startbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
            char endbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
            void *begin = &(startbuf[0]);
            void *end = &(endbuf[0]);
            config->fCreateIterators( alternative, &begin, &end, newProxy);
            // We can not get here with a split vector of pointer, so we can indeed assume
            // that actions->fConfiguration != null.
            buf.ApplySequence(*actions, begin, end);
            if (begin != &(startbuf[0])) {
               // assert(end != endbuf);
               config->fDeleteTwoIterators(begin,end);
            }
         }
         newProxy->Commit(alternative);
      }
   }

   INLINE_TEMPLATE_ARGS void ReadArraySTLMemberWiseChangedClass(TBuffer &buf, void *addr, const TConfiguration *conf, Version_t vers)
   {
      // Collection was saved member-wise

      TConfigSTL *config = (TConfigSTL*)conf;

      vers &= ~( TBufferFile::kStreamedMemberWise );

      TClass *newClass = config->fNewClass;
      TClass *oldClass = config->fOldClass;

      if( vers < 8 ) {
         Error( "ReadSTLMemberWiseChangedClass", "Unfortunately, version %d of TStreamerInfo (used in %s) did not record enough information to convert a %s into a %s.",
               vers, buf.GetParent() ? buf.GetParent()->GetName() : "memory/socket", oldClass->GetName(), newClass->GetName() );
      } else {

         Version_t vClVersion = buf.ReadVersionForMemberWise( oldClass->GetCollectionProxy()->GetValueClass() );

         TVirtualCollectionProxy *newProxy = newClass->GetCollectionProxy();
         TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();

         int objectSize = newClass->Size();
         char *obj = (char*)addr;
         char *endobj = obj + conf->fLength*objectSize;

         for(; obj<endobj; obj+=objectSize) {
            TVirtualCollectionProxy::TPushPop helper( newProxy, (char*)obj );
            Int_t nobjects;
            buf.ReadInt(nobjects);
            void* alternative = newProxy->Allocate(nobjects,true);
            if (nobjects) {
               TActionSequence *actions = newProxy->GetConversionReadMemberWiseActions( oldProxy->GetValueClass(), vClVersion );
               char startbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
               char endbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
               void *begin = &(startbuf[0]);
               void *end = &(endbuf[0]);
               config->fCreateIterators( alternative, &begin, &end, newProxy);
               // We can not get here with a split vector of pointer, so we can indeed assume
               // that actions->fConfiguration != null.
               buf.ApplySequence(*actions, begin, end);
               if (begin != &(startbuf[0])) {
                  // assert(end != endbuf);
                  config->fDeleteTwoIterators(begin,end);
               }
            }
            newProxy->Commit(alternative);
         }
      }
   }


   INLINE_TEMPLATE_ARGS void ReadSTLObjectWiseFastArray(TBuffer &buf, void *addr, const TConfiguration *conf, Version_t /* vers */, UInt_t /* start */)
   {
      TConfigSTL *config = (TConfigSTL*)conf;
      // Idea: This needs to be unrolled, it currently calls the TGenCollectionStreamer ....
      buf.ReadFastArray(addr,config->fNewClass,conf->fLength,(TMemberStreamer*)0,config->fOldClass);
   }
   INLINE_TEMPLATE_ARGS void ReadSTLObjectWiseStreamer(TBuffer &buf, void *addr, const TConfiguration *conf, Version_t /* vers */, UInt_t /* start */)
   {
      TConfigSTL *config = (TConfigSTL*)conf;
      (*config->fStreamer)(buf,addr,conf->fLength);
   }
   INLINE_TEMPLATE_ARGS void ReadSTLObjectWiseFastArrayV2(TBuffer &buf, void *addr, const TConfiguration *conf, Version_t vers, UInt_t start)
   {
      // case of old TStreamerInfo

      TConfigSTL *config = (TConfigSTL*)conf;
      //  Backward compatibility. Some TStreamerElement's where without
      //  Streamer but were not removed from element list
      if (config->fIsSTLBase || vers == 0) {
         buf.SetBufferOffset(start);  //there is no byte count
      }
      // Idea: This needs to be unrolled, it currently calls the TGenCollectionStreamer ....
      buf.ReadFastArray(addr,config->fNewClass,conf->fLength,(TMemberStreamer*)0,config->fOldClass);
   }
   INLINE_TEMPLATE_ARGS void ReadSTLObjectWiseStreamerV2(TBuffer &buf, void *addr, const TConfiguration *conf, Version_t vers, UInt_t start)
   {
      // case of old TStreamerInfo

      TConfigSTL *config = (TConfigSTL*)conf;
      //  Backward compatibility. Some TStreamerElement's where without
      //  Streamer but were not removed from element list
      if (config->fIsSTLBase || vers == 0) {
         buf.SetBufferOffset(start);  //there is no byte count
      }
      (*config->fStreamer)(buf,addr,conf->fLength);
   }

   template <void (*memberwise)(TBuffer&,void *,const TConfiguration*, Version_t),
             void (*objectwise)(TBuffer&,void *,const TConfiguration*, Version_t, UInt_t)>
   INLINE_TEMPLATE_ARGS Int_t ReadSTL(TBuffer &buf, void *addr, const TConfiguration *conf)
   {
      TConfigSTL *config = (TConfigSTL*)conf;
      UInt_t start, count;
      Version_t vers = buf.ReadVersion(&start, &count, config->fOldClass);
      if ( vers & TBufferFile::kStreamedMemberWise ) {
         memberwise(buf,((char*)addr)+config->fOffset,config, vers);
      } else {
         objectwise(buf,((char*)addr)+config->fOffset,config, vers, start);
      }
      buf.CheckByteCount(start,count,config->fTypeName);
      return 0;
   }

   template <typename From, typename To>
   struct ConvertBasicType {
      static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *addr, const TConfiguration *config)
      {
         // Simple conversion from a 'From' on disk to a 'To' in memory.
         From temp;
         buf >> temp;
         *(To*)( ((char*)addr) + config->fOffset ) = (To)temp;
         return 0;
      }
   };

   template <typename To>
   struct ConvertBasicType<BitsMarker,To> {
      static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *addr, const TConfiguration *config)
      {
         // Simple conversion from a 'From' on disk to a 'To' in memory
         UInt_t temp;
         buf >> temp;

         if ((temp & kIsReferenced) != 0) {
            HandleReferencedTObject(buf,addr,config);
         }

         *(To*)( ((char*)addr) + config->fOffset ) = (To)temp;
         return 0;
      }
   };

   template <typename From, typename To>
   struct ConvertBasicType<WithFactorMarker<From>,To> {
      static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *addr, const TConfiguration *config)
      {
         // Simple conversion from a 'From' on disk to a 'To' in memory.
         TConfWithFactor *conf = (TConfWithFactor *)config;
         From temp;
         buf.ReadWithFactor(&temp, conf->fFactor, conf->fXmin);
         *(To*)( ((char*)addr) + config->fOffset ) = (To)temp;
         return 0;
      }
   };

   template <typename From, typename To>
   struct ConvertBasicType<NoFactorMarker<From>,To> {
      static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *addr, const TConfiguration *config)
      {
         // Simple conversion from a 'From' on disk to a 'To' in memory.
         TConfNoFactor *conf = (TConfNoFactor *)config;
         From temp;
         buf.ReadWithNbits(&temp, conf->fNbits);
         *(To*)( ((char*)addr) + config->fOffset ) = (To)temp;
         return 0;
      }
   };

   class TConfigurationPushDataCache : public TConfiguration {
      // Configuration object for the PushDataCache case.
   public:
      TVirtualArray *fOnfileObject;

      TConfigurationPushDataCache(TVirtualStreamerInfo *info, TVirtualArray *onfileObject, Int_t offset) :
         TConfiguration(info, -1, nullptr, offset), fOnfileObject(onfileObject)
      {}

      virtual void Print() const {
         TStreamerInfo *info = (TStreamerInfo*)fInfo;
         if (fOnfileObject)
            printf("StreamerInfoAction, class:%s, PushDataCache offset=%d\n",
                   info->GetClass()->GetName(), fOffset);
         else
            printf("StreamerInfoAction, class:%s, PopDataCache offset=%d\n",
                   info->GetClass()->GetName(), fOffset);
      }
      virtual void PrintDebug(TBuffer &buffer, void *object) const {
         if (gDebug > 1) {
            TStreamerInfo *info = (TStreamerInfo*)fInfo;
            printf("StreamerInfoAction, class:%s, %sDataCache, bufpos=%d, arr=%p, offset=%d, onfileObject=%p\n",
                  info->GetClass()->GetName(), fOnfileObject ? "Push" : "Pop", buffer.Length(), object, fOffset, fOnfileObject);

         }
      }
   };

   Int_t PushDataCache(TBuffer &b, void *, const TConfiguration *conf)
   {
      TConfigurationPushDataCache *config = (TConfigurationPushDataCache*)conf;
      auto onfileObject = config->fOnfileObject;

      // onfileObject->SetSize(1);
      b.PushDataCache( onfileObject );

      return 0;
   }

   Int_t PushDataCacheVectorPtr(TBuffer &b, void *, const void *, const TConfiguration *conf)
   {
      TConfigurationPushDataCache *config = (TConfigurationPushDataCache*)conf;
      auto onfileObject = config->fOnfileObject;

      // onfileObject->SetSize(n);
      b.PushDataCache( onfileObject );

      return 0;
   }

   Int_t PushDataCacheGenericCollection(TBuffer &b, void *, const void *, const TLoopConfiguration *loopconfig, const TConfiguration *conf)
   {
      TConfigurationPushDataCache *config = (TConfigurationPushDataCache*)conf;
      auto onfileObject = config->fOnfileObject;

      TVirtualCollectionProxy *proxy = ((TGenericLoopConfig*)loopconfig)->fProxy;
      UInt_t n = proxy->Size();

      onfileObject->SetSize(n);
      b.PushDataCache( onfileObject );

      return 0;
   }

   Int_t PopDataCache(TBuffer &b, void *, const TConfiguration *)
   {
      b.PopDataCache();
      return 0;
   }

   Int_t PopDataCacheVectorPtr(TBuffer &b, void *, const void *, const TConfiguration *)
   {
      b.PopDataCache();
      return 0;
   }

   Int_t PopDataCacheGenericCollection(TBuffer &b, void *, const void *, const TLoopConfiguration *, const TConfiguration *)
   {
      b.PopDataCache();
      return 0;
   }

   class TConfigurationUseCache : public TConfiguration {
      // Configuration object for the UseCache case.
   public:
      TConfiguredAction fAction;
      Bool_t            fNeedRepeat;

      TConfigurationUseCache(TVirtualStreamerInfo *info, TConfiguredAction &action, Bool_t repeat) :
              TConfiguration(info,action.fConfiguration->fElemId,action.fConfiguration->fCompInfo,action.fConfiguration->fOffset),fAction(action),fNeedRepeat(repeat) {};
      virtual void PrintDebug(TBuffer &b, void *addr) const
      {
         if (gDebug > 1) {
            // Idea: We should print the name of the action function.
            TStreamerInfo *info = (TStreamerInfo*)fInfo;
            TStreamerElement *aElement = fCompInfo->fElem;
            fprintf(stdout,"StreamerInfoAction, class:%s, name=%s, fType[%d]=%d,"
                   " %s, bufpos=%d, arr=%p, eoffset=%d, Redirect=%p\n",
                   info->GetClass()->GetName(),aElement->GetName(),fElemId,fCompInfo->fType,
                   aElement->ClassName(),b.Length(),addr, 0,b.PeekDataCache() ? b.PeekDataCache()->GetObjectAt(0) : 0);
         }

      }
      virtual ~TConfigurationUseCache() {};
      virtual TConfiguration *Copy() {
         TConfigurationUseCache *copy = new TConfigurationUseCache(*this);
         fAction.fConfiguration = copy->fAction.fConfiguration->Copy(); // since the previous allocation did a 'move' of fAction we need to fix it.
         return copy;
      }
   };

   INLINE_TEMPLATE_ARGS Int_t UseCache(TBuffer &b, void *addr, const TConfiguration *conf)
   {
      TConfigurationUseCache *config = (TConfigurationUseCache*)conf;

      Int_t bufpos = b.Length();
      TVirtualArray *cached = b.PeekDataCache();
      if (cached==0) {
         TStreamerElement *aElement = conf->fCompInfo->fElem;
         TStreamerInfo *info = (TStreamerInfo*)conf->fInfo;
         Warning("ReadBuffer","Skipping %s::%s because the cache is missing.",info->GetName(),aElement->GetName());
         char *ptr = (char*)addr;
         info->ReadBufferSkip(b,&ptr,config->fCompInfo,config->fCompInfo->fType+TStreamerInfo::kSkip,aElement,1,0);
      } else {
         config->fAction(b, (*cached)[0]);
      }
      // Idea: Factor out this 'if' to a UseCacheRepeat function
      if (config->fNeedRepeat) {
         b.SetBufferOffset(bufpos);
      }
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t UseCacheVectorPtrLoop(TBuffer &b, void *start, const void *end, const TConfiguration *conf)
   {
      TConfigurationUseCache *config = (TConfigurationUseCache*)conf;
      Int_t bufpos = b.Length();

      TVirtualArray *cached = b.PeekDataCache();
      if (cached==0) {
         TStreamerElement *aElement = config->fCompInfo->fElem;
         TStreamerInfo *info = (TStreamerInfo*)config->fInfo;
         Warning("ReadBuffer","Skipping %s::%s because the cache is missing.",info->GetName(),aElement->GetName());
         char *ptr = (char*)start;
         UInt_t n = (((void**)end)-((void**)start));
         info->ReadBufferSkip(b,&ptr,config->fCompInfo,conf->fCompInfo->fType+TStreamerInfo::kSkip,aElement,n,0);
      } else {
         TVectorLoopConfig cached_config( nullptr, cached->fClass->Size(), /* read */ kTRUE );
         void *cached_start = (*cached)[0];
         void *cached_end = ((char*)cached_start) + cached->fSize * cached_config.fIncrement;
         config->fAction(b,cached_start,cached_end,&cached_config);
      }
      // Idea: Factor out this 'if' to a UseCacheRepeat function
      if (config->fNeedRepeat) {
         b.SetBufferOffset(bufpos);
      }
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t UseCacheVectorLoop(TBuffer &b, void *start, const void *end, const TLoopConfiguration *loopconf, const TConfiguration *conf)
   {
      TConfigurationUseCache *config = (TConfigurationUseCache*)conf;

      Int_t bufpos = b.Length();
      TVirtualArray *cached = b.PeekDataCache();
      if (cached==0) {
         TStreamerElement *aElement = config->fCompInfo->fElem;
         TStreamerInfo *info = (TStreamerInfo*)config->fInfo;
         Warning("ReadBuffer","Skipping %s::%s because the cache is missing.",info->GetName(),aElement->GetName());
         char *ptr = (char*)start;
         UInt_t n = (((char*)end)-((char*)start))/((TVectorLoopConfig*)loopconf)->fIncrement;
         info->ReadBufferSkip(b,&ptr,config->fCompInfo,config->fCompInfo->fType+TStreamerInfo::kSkip,aElement,n,0);
      } else {
         TVectorLoopConfig cached_config( nullptr, cached->fClass->Size(), /* read */ kTRUE );
         void *cached_start = (*cached)[0];
         void *cached_end = ((char*)cached_start) + cached->fSize * cached_config.fIncrement;
         config->fAction(b,cached_start,cached_end,&cached_config);
      }
      // Idea: Factor out this 'if' to a UseCacheRepeat function
      if (config->fNeedRepeat) {
         b.SetBufferOffset(bufpos);
      }
      return 0;
   }

   INLINE_TEMPLATE_ARGS Int_t UseCacheGenericCollection(TBuffer &b, void *, const void *, const TLoopConfiguration *loopconfig, const TConfiguration *conf)
   {
      TConfigurationUseCache *config = (TConfigurationUseCache*)conf;

      Int_t bufpos = b.Length();
      TVirtualArray *cached = b.PeekDataCache();
      if (cached==0) {
         TStreamerElement *aElement = config->fCompInfo->fElem;
         TStreamerInfo *info = (TStreamerInfo*)config->fInfo;

         TVirtualCollectionProxy *proxy = ((TGenericLoopConfig*)loopconfig)->fProxy;
         Warning("ReadBuffer","Skipping %s::%s because the cache is missing.",info->GetName(),aElement->GetName());
         UInt_t n = proxy->Size();
         info->ReadBufferSkip(b, *proxy,config->fCompInfo,config->fCompInfo->fType+TStreamerInfo::kSkip,aElement,n,0);
      } else {
         TVectorLoopConfig cached_config( nullptr, cached->fClass->Size(), /* read */ kTRUE );
         void *cached_start = (*cached)[0];
         void *cached_end = ((char*)cached_start) + cached->fSize * cached_config.fIncrement;
         config->fAction(b,cached_start,cached_end,&cached_config);
      }
      // Idea: Factor out this 'if' to a UseCacheRepeat function
      if (config->fNeedRepeat) {
         b.SetBufferOffset(bufpos);
      }
      return 0;
   }

   // Support for collections.

   Int_t ReadLoopInvalid(TBuffer &, void *, const void *, const TConfiguration *config)
   {
      Fatal("ApplySequence","The sequence of actions to read %s:%d member-wise was not initialized.",config->fInfo->GetName(),config->fInfo->GetClassVersion());
      return 0;
   }

   Int_t WriteLoopInvalid(TBuffer &, void *, const void *, const TConfiguration *config)
   {
      Fatal("ApplySequence","The sequence of actions to write %s:%d member-wise was not initialized.",config->fInfo->GetName(),config->fInfo->GetClassVersion());
      return 0;
   }

   enum ESelectLooper { kVectorLooper, kVectorPtrLooper, kAssociativeLooper, kGenericLooper };

   ESelectLooper SelectLooper(TVirtualCollectionProxy &proxy)
   {
      if ( (proxy.GetCollectionType() == ROOT::kSTLvector) || (proxy.GetProperties() & TVirtualCollectionProxy::kIsEmulated) ) {
         if (proxy.GetProperties() & TVirtualCollectionProxy::kCustomAlloc)
            return kGenericLooper;
         else
            return kVectorLooper;
      } else if (proxy.GetCollectionType() == ROOT::kSTLset || proxy.GetCollectionType() == ROOT::kSTLunorderedset
                 || proxy.GetCollectionType() == ROOT::kSTLmultiset || proxy.GetCollectionType() == ROOT::kSTLunorderedmultiset
                 || proxy.GetCollectionType() == ROOT::kSTLmap || proxy.GetCollectionType() == ROOT::kSTLmultimap
                 || proxy.GetCollectionType() == ROOT::kSTLunorderedmap || proxy.GetCollectionType() == ROOT::kSTLunorderedmultimap
                 || proxy.GetCollectionType() == ROOT::kSTLbitset) {
         return kAssociativeLooper;
      } else {
         return kGenericLooper;
      }
   }

   struct VectorLooper {

      template <typename T>
      static INLINE_TEMPLATE_ARGS Int_t ReadBasicType(TBuffer &buf, void *iter, const void *end, const TLoopConfiguration *loopconfig, const TConfiguration *config)
      {
         const Int_t incr = ((TVectorLoopConfig*)loopconfig)->fIncrement;
         iter = (char*)iter + config->fOffset;
         end = (char*)end + config->fOffset;
         for(; iter != end; iter = (char*)iter + incr ) {
            T *x = (T*) ((char*) iter);
            buf >> *x;
         }
         return 0;
      }

      template <typename From, typename To>
      struct ConvertBasicType {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *iter, const void *end, const TLoopConfiguration *loopconfig, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.
            From temp;
            const Int_t incr = ((TVectorLoopConfig*)loopconfig)->fIncrement;
            iter = (char*)iter + config->fOffset;
            end = (char*)end + config->fOffset;
            for(; iter != end; iter = (char*)iter + incr ) {
               buf >> temp;
               *(To*)( ((char*)iter) ) = (To)temp;
            }
            return 0;
         }
      };

      template <typename To>
      struct ConvertBasicType<BitsMarker,To> {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *iter, const void *end, const TLoopConfiguration *loopconfig, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.
            UInt_t temp;
            const Int_t incr = ((TVectorLoopConfig*)loopconfig)->fIncrement;
            iter = (char*)iter + config->fOffset;
            end = (char*)end + config->fOffset;
            for(; iter != end; iter = (char*)iter + incr ) {
               buf >> temp;

               if ((temp & kIsReferenced) != 0) {
                  HandleReferencedTObject(buf, (char*)iter - config->fOffset, config);
               }

               *(To*)( ((char*)iter) ) = (To)temp;
            }
            return 0;
         }
      };

      template <typename From, typename To>
      struct ConvertBasicType<WithFactorMarker<From>,To> {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *iter, const void *end, const TLoopConfiguration *loopconfig, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.
            TConfWithFactor *conf = (TConfWithFactor *)config;
            From temp;
            const Int_t incr = ((TVectorLoopConfig*)loopconfig)->fIncrement;
            iter = (char*)iter + config->fOffset;
            end = (char*)end + config->fOffset;
            for(; iter != end; iter = (char*)iter + incr ) {
               buf.ReadWithFactor(&temp, conf->fFactor, conf->fXmin);
               *(To*)( ((char*)iter) ) = (To)temp;
            }
            return 0;
         }
      };

      template <typename From, typename To>
      struct ConvertBasicType<NoFactorMarker<From>,To> {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *iter, const void *end, const TLoopConfiguration *loopconfig, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.
            TConfNoFactor *conf = (TConfNoFactor *)config;
            From temp;
            const Int_t incr = ((TVectorLoopConfig*)loopconfig)->fIncrement;
            iter = (char*)iter + config->fOffset;
            end = (char*)end + config->fOffset;
            for(; iter != end; iter = (char*)iter + incr ) {
               buf.ReadWithNbits(&temp, conf->fNbits);
               *(To*)( ((char*)iter) ) = (To)temp;
            }
            return 0;
         }
      };

      template <typename T>
      static INLINE_TEMPLATE_ARGS Int_t WriteBasicType(TBuffer &buf, void *iter, const void *end, const TLoopConfiguration *loopconfig, const TConfiguration *config)
      {
         const Int_t incr = ((TVectorLoopConfig*)loopconfig)->fIncrement;
         iter = (char*)iter + config->fOffset;
         end = (char*)end + config->fOffset;
         for(; iter != end; iter = (char*)iter + incr ) {
            T *x = (T*) ((char*) iter);
            buf << *x;
         }
         return 0;
      }

      template <Int_t (*iter_action)(TBuffer&,void *,const TConfiguration*)>
      static INLINE_TEMPLATE_ARGS Int_t ReadAction(TBuffer &buf, void *start, const void *end, const TLoopConfiguration *loopconfig, const TConfiguration *config)
      {
         const Int_t incr = ((TVectorLoopConfig*)loopconfig)->fIncrement;
         //Idea: can we factor out the addition of fOffset
         //  iter = (char*)iter + config->fOffset;
         for(void *iter = start; iter != end; iter = (char*)iter + incr ) {
            iter_action(buf, iter, config);
         }
         return 0;
      }

      static INLINE_TEMPLATE_ARGS Int_t ReadBase(TBuffer &buf, void *start, const void *end, const TLoopConfiguration * loopconfig, const TConfiguration *config)
      {
         // Well the implementation is non trivial since we do not have a proxy for the container of _only_ the base class.  For now
         // punt.

         UInt_t incr = ((TVectorLoopConfig*)loopconfig)->fIncrement;
         UInt_t n = (((char*)end)-((char*)start))/incr;
         char **arrptr = new char*[n];
         UInt_t i = 0;
         for(void *iter = start; iter != end; iter = (char*)iter + incr, ++i ) {
            arrptr[i] = (char*)iter;
         }
         ((TStreamerInfo*)config->fInfo)->ReadBuffer(buf, arrptr, &(config->fCompInfo), /*first*/ 0, /*last*/ 1, /*narr*/ n, config->fOffset, 1|2 );
         delete [] arrptr;

         //      // Idea: need to cache this result!
         //      TStreamerInfo *info = (TStreamerInfo*)config->fInfo;
         //      TStreamerElement *aElement = (TStreamerElement*)info->GetElem(config->fElemId);
         //
         //      *Int_t clversion = ((TStreamerBase*)aElement)->Get BaseVersion();
         //      *TClass *cle = aElement->GetNewBaseClass();
         //      *(TSequence *actions = CreateReadMemberWiseActions( cle->GetStreamerInfo(clversion), ???? );
         //
         //      TSequence *actions = CreateReadMemberWiseActions( ((TStreamerBase*)aElement)->GetBaseStreamerInfo(), ???? );
         //
         //      actions->ReadBuffer(b,start,end);
         //      delete actions;

         //      const Int_t incr = ((TVectorLoopConfig*)loopconfig)->fIncrement;
         //      for(void *iter = start; iter != end; iter = (char*)iter + incr )
         //      {
         //         ((TStreamerInfo*)(((TStreamerBase*)aElement)->GetBaseStreamerInfo())->ReadBuffer(b,arr,-1,narr,ioffset,arrayMode);
         //
         //         ((TStreamerInfo*)config->fInfo)->ReadBuffer(buf, (char**)&iter, config->fElemId, 1, config->fOffset, 1|2 );
         //      }
         return 0;
      }

      static INLINE_TEMPLATE_ARGS Int_t GenericRead(TBuffer &buf, void *start, const void *end, const TLoopConfiguration * loopconfig, const TConfiguration *config)
      {
         // Well the implementation is non trivial. For now punt.

         UInt_t incr = ((TVectorLoopConfig*)loopconfig)->fIncrement;
         UInt_t n = (((char*)end)-((char*)start))/incr;
         char **arrptr = new char*[n];
         UInt_t i = 0;
         for(void *iter = start; iter != end; iter = (char*)iter + incr, ++i ) {
            arrptr[i] = (char*)iter;
         }
         ((TStreamerInfo*)config->fInfo)->ReadBuffer(buf, arrptr, &(config->fCompInfo), /*first*/ 0, /*last*/ 1, /*narr*/ n, config->fOffset, 1|2 );
         delete [] arrptr;
         return 0;
      }

      static INLINE_TEMPLATE_ARGS Int_t GenericWrite(TBuffer &buf, void *start, const void *end, const TLoopConfiguration * loopconfig, const TConfiguration *config)
      {
         // Well the implementation is non trivial. For now punt.

         UInt_t incr = ((TVectorLoopConfig*)loopconfig)->fIncrement;
         UInt_t n = (((char*)end)-((char*)start))/incr;
         char **arrptr = new char*[n];
         UInt_t i = 0;
         for(void *iter = start; iter != end; iter = (char*)iter + incr, ++i ) {
            arrptr[i] = (char*)iter;
         }
         ((TStreamerInfo*)config->fInfo)->WriteBufferAux(buf, arrptr, &(config->fCompInfo), /*first*/ 0, /*last*/ 1, n, config->fOffset, 1|2 );
         delete [] arrptr;
         return 0;
      }

      template <typename T>
      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionBasicType(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         // Collection of numbers.  Memberwise or not, it is all the same.

         TConfigSTL *config = (TConfigSTL*)conf;
         UInt_t start, count;
         /* Version_t vers = */ buf.ReadVersion(&start, &count, config->fOldClass);

         std::vector<T> *const vec = (std::vector<T>*)(((char*)addr)+config->fOffset);
         Int_t nvalues;
         buf.ReadInt(nvalues);
         vec->resize(nvalues);

#ifdef R__VISUAL_CPLUSPLUS
         if (nvalues <= 0) {
            buf.CheckByteCount(start,count,config->fTypeName);
            return 0;
         }
#endif
         T *begin = &(*vec->begin());
         buf.ReadFastArray(begin, nvalues);

         buf.CheckByteCount(start,count,config->fTypeName);
         return 0;
      }

      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionBool(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         // Collection of numbers.  Memberwise or not, it is all the same.

         TConfigSTL *config = (TConfigSTL*)conf;
         UInt_t start, count;
         /* Version_t vers = */ buf.ReadVersion(&start, &count, config->fOldClass);

         std::vector<bool> *const vec = (std::vector<bool>*)(((char*)addr)+config->fOffset);
         Int_t nvalues;
         buf.ReadInt(nvalues);
         vec->resize(nvalues);

         bool *items = new bool[nvalues];
         buf.ReadFastArray(items, nvalues);
         for(Int_t i = 0 ; i < nvalues; ++i) {
            (*vec)[i] = items[i];
         }
         delete [] items;

         // We could avoid the call to ReadFastArray, and we could
         // the following, however this breaks TBufferXML ...
         // for(Int_t i = 0 ; i < nvalues; ++i) {
         //    bool tmp; buf >> tmp;
         //    (*vec)[i] = tmp;
         // }

         buf.CheckByteCount(start,count,config->fTypeName);
         return 0;
      }

      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionFloat16(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         // Collection of numbers.  Memberwise or not, it is all the same.

         TConfigSTL *config = (TConfigSTL*)conf;
         UInt_t start, count;
         /* Version_t vers = */ buf.ReadVersion(&start, &count, config->fOldClass);

         std::vector<float> *const vec = (std::vector<float>*)(((char*)addr)+config->fOffset);
         Int_t nvalues;
         buf.ReadInt(nvalues);
         vec->resize(nvalues);

#ifdef R__VISUAL_CPLUSPLUS
         if (nvalues <= 0) {
            buf.CheckByteCount(start,count,config->fTypeName);
            return 0;
         }
#endif
         float *begin = &(*vec->begin());
         buf.ReadFastArrayFloat16(begin, nvalues);

         buf.CheckByteCount(start,count,config->fTypeName);
         return 0;
      }

      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionDouble32(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         // Collection of numbers.  Memberwise or not, it is all the same.

         TConfigSTL *config = (TConfigSTL*)conf;
         UInt_t start, count;
         /* Version_t vers = */ buf.ReadVersion(&start, &count, config->fOldClass);

         std::vector<double> *const vec = (std::vector<double>*)(((char*)addr)+config->fOffset);
         Int_t nvalues;
         buf.ReadInt(nvalues);
         vec->resize(nvalues);

#ifdef R__VISUAL_CPLUSPLUS
         if (nvalues <= 0) {
            buf.CheckByteCount(start,count,config->fTypeName);
            return 0;
         }
#endif
         double *begin = &(*vec->begin());
         buf.ReadFastArrayDouble32(begin, nvalues);

         buf.CheckByteCount(start,count,config->fTypeName);
         return 0;
      }

      template <typename From, typename To>
      struct ConvertCollectionBasicType {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *addr, const TConfiguration *conf)
         {
            // Collection of numbers.  Memberwise or not, it is all the same.

            TConfigSTL *config = (TConfigSTL*)conf;
            UInt_t start, count;
            /* Version_t vers = */ buf.ReadVersion(&start, &count, config->fOldClass);

            std::vector<To> *const vec = (std::vector<To>*)(((char*)addr)+config->fOffset);
            Int_t nvalues;
            buf.ReadInt(nvalues);
            vec->resize(nvalues);

            From *temp = new From[nvalues];
            buf.ReadFastArray(temp, nvalues);
            for(Int_t ind = 0; ind < nvalues; ++ind) {
               (*vec)[ind] = (To)temp[ind];
            }
            delete [] temp;

            buf.CheckByteCount(start,count,config->fTypeName);
            return 0;
         }
      };

      template <typename From, typename To>
      struct ConvertCollectionBasicType<NoFactorMarker<From>,To> {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *addr, const TConfiguration *conf)
         {
            // Collection of numbers.  Memberwise or not, it is all the same.

            TConfigSTL *config = (TConfigSTL*)conf;
            UInt_t start, count;
            /* Version_t vers = */ buf.ReadVersion(&start, &count, config->fOldClass);

            std::vector<To> *const vec = (std::vector<To>*)(((char*)addr)+config->fOffset);
            Int_t nvalues;
            buf.ReadInt(nvalues);
            vec->resize(nvalues);

            From *temp = new From[nvalues];
            buf.ReadFastArrayWithNbits(temp, nvalues, 0);
            for(Int_t ind = 0; ind < nvalues; ++ind) {
               (*vec)[ind] = (To)temp[ind];
            }
            delete [] temp;

            buf.CheckByteCount(start,count,config->fTypeName);
            return 0;
         }
      };

      template <typename To>
      static INLINE_TEMPLATE_ARGS Int_t ConvertCollectionDouble32(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         // Collection of numbers.  Memberwise or not, it is all the same.

         TConfigSTL *config = (TConfigSTL*)conf;
         UInt_t start, count;
         /* Version_t vers = */ buf.ReadVersion(&start, &count, config->fOldClass);

         std::vector<To> *const vec = (std::vector<To>*)(((char*)addr)+config->fOffset);
         Int_t nvalues;
         buf.ReadInt(nvalues);
         vec->resize(nvalues);

         Double32_t *temp = new Double32_t[nvalues];
         buf.ReadFastArrayDouble32(temp, nvalues);
         for(Int_t ind = 0; ind < nvalues; ++ind) {
            (*vec)[ind] = (To)temp[ind];
         }
         delete [] temp;

         buf.CheckByteCount(start,count,config->fTypeName);
         return 0;
      }

   };

   struct VectorPtrLooper {

      template <typename T>
      static INLINE_TEMPLATE_ARGS Int_t ReadBasicType(TBuffer &buf, void *iter, const void *end, const TConfiguration *config)
      {
         const Int_t offset = config->fOffset;

         for(; iter != end; iter = (char*)iter + sizeof(void*) ) {
            T *x = (T*)( ((char*) (*(void**)iter) ) + offset );
            buf >> *x;
         }
         return 0;
      }

      template <typename From, typename To>
      struct ConvertBasicType {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *iter, const void *end, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.
            From temp;
            const Int_t offset = config->fOffset;
            for(; iter != end; iter = (char*)iter + sizeof(void*) ) {
               buf >> temp;
               To *x = (To*)( ((char*) (*(void**)iter) ) + offset );
               *x = (To)temp;
            }
            return 0;
         }
      };

      template <typename To>
      struct ConvertBasicType<BitsMarker,To> {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *iter, const void *end, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.
            UInt_t temp;
            const Int_t offset = config->fOffset;
            for(; iter != end; iter = (char*)iter + sizeof(void*) ) {
               buf >> temp;

               if ((temp & kIsReferenced) != 0) {
                  HandleReferencedTObject(buf,*(void**)iter,config);
               }

               To *x = (To*)( ((char*) (*(void**)iter) ) + offset );
               *x = (To)temp;
            }
            return 0;
         }
      };

      template <typename From, typename To>
      struct ConvertBasicType<WithFactorMarker<From>,To> {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *iter, const void *end, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.
            TConfWithFactor *conf = (TConfWithFactor *)config;
            From temp;
            const Int_t offset = config->fOffset;
            for(; iter != end; iter = (char*)iter + sizeof(void*) ) {
               buf.ReadWithFactor(&temp, conf->fFactor, conf->fXmin);
               To *x = (To*)( ((char*) (*(void**)iter) ) + offset );
               *x = (To)temp;
            }
            return 0;
         }
      };

      template <typename From, typename To>
      struct ConvertBasicType<NoFactorMarker<From>,To> {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *iter, const void *end, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.
            TConfNoFactor *conf = (TConfNoFactor *)config;
            From temp;
            const Int_t offset = config->fOffset;
            for(; iter != end; iter = (char*)iter + sizeof(void*) ) {
               buf.ReadWithNbits(&temp, conf->fNbits);
               To *x = (To*)( ((char*) (*(void**)iter) ) + offset );
               *x = (To)temp;
            }
            return 0;
         }
      };

      template <typename T>
      static INLINE_TEMPLATE_ARGS Int_t WriteBasicType(TBuffer &buf, void *iter, const void *end, const TConfiguration *config)
      {
         const Int_t offset = config->fOffset;

         for(; iter != end; iter = (char*)iter + sizeof(void*) ) {
            T *x = (T*)( ((char*) (*(void**)iter) ) + offset );
            buf << *x;
         }
         return 0;
      }

      template <Int_t (*action)(TBuffer&,void *,const TConfiguration*)>
      static INLINE_TEMPLATE_ARGS Int_t ReadAction(TBuffer &buf, void *start, const void *end, const TConfiguration *config)
      {
         for(void *iter = start; iter != end; iter = (char*)iter + sizeof(void*) ) {
            action(buf, *(void**)iter, config);
         }
         return 0;
      }

      static INLINE_TEMPLATE_ARGS Int_t ReadBase(TBuffer &buf, void *start, const void *end, const TConfiguration *config)
      {
         // Well the implementation is non trivial since we do not have a proxy for the container of _only_ the base class.  For now
         // punt.

         return GenericRead(buf,start,end,config);
      }

      static INLINE_TEMPLATE_ARGS Int_t GenericRead(TBuffer &buf, void *iter, const void *end, const TConfiguration *config)
      {
         Int_t n = ( ((void**)end) - ((void**)iter) );
         char **arr = (char**)iter;
         return ((TStreamerInfo*)config->fInfo)->ReadBuffer(buf, arr, &(config->fCompInfo), /*first*/ 0, /*last*/ 1, /*narr*/ n, config->fOffset, 1|2 );
      }

      static INLINE_TEMPLATE_ARGS Int_t GenericWrite(TBuffer &buf, void *iter, const void *end, const TConfiguration *config)
      {
         Int_t n = ( ((void**)end) - ((void**)iter) );
         char **arr = (char**)iter;
         return ((TStreamerInfo*)config->fInfo)->WriteBufferAux(buf, arr, &(config->fCompInfo), /*first*/ 0, /*last*/ 1, n, config->fOffset, 1|2 );
      }

   };

   struct AssociativeLooper {

      template <typename T>
      static INLINE_TEMPLATE_ARGS void SimpleRead(TBuffer &buf, void *addr, Int_t nvalues)
      {
         buf.ReadFastArray((T*)addr, nvalues);
      }

      static INLINE_TEMPLATE_ARGS void SimpleReadFloat16(TBuffer &buf, void *addr, Int_t nvalues)
      {
         buf.ReadFastArrayFloat16((float*)addr, nvalues);
      }

      static INLINE_TEMPLATE_ARGS void SimpleReadDouble32(TBuffer &buf, void *addr, Int_t nvalues)
      {
         buf.ReadFastArrayDouble32((double*)addr, nvalues);
      }

      template <typename T,void (*action)(TBuffer&,void *,Int_t)>
      static INLINE_TEMPLATE_ARGS Int_t ReadNumericalCollection(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         // Collection of numbers.  Memberwise or not, it is all the same.

         TConfigSTL *config = (TConfigSTL*)conf;
         UInt_t start, count;
         /* Version_t vers = */ buf.ReadVersion(&start, &count, config->fOldClass);

         TClass *newClass = config->fNewClass;
         TVirtualCollectionProxy *newProxy = newClass->GetCollectionProxy();
         TVirtualCollectionProxy::TPushPop helper( newProxy, ((char*)addr)+config->fOffset );

         Int_t nvalues;
         buf.ReadInt(nvalues);
         void* alternative = newProxy->Allocate(nvalues,true);
         if (nvalues) {
            char startbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
            char endbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
            void *begin = &(startbuf[0]);
            void *end = &(endbuf[0]);
            config->fCreateIterators(alternative, &begin, &end, newProxy);
            // We can not get here with a split vector of pointer, so we can indeed assume
            // that actions->fConfiguration != null.

            action(buf,begin,nvalues);

            if (begin != &(startbuf[0])) {
               // assert(end != endbuf);
               config->fDeleteTwoIterators(begin,end);
            }
         }
         newProxy->Commit(alternative);

         buf.CheckByteCount(start,count,config->fTypeName);
         return 0;
      }

      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionBool(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         return ReadNumericalCollection<bool,SimpleRead<bool> >(buf,addr,conf);
      }

      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionFloat16(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         return ReadNumericalCollection<Float_t,SimpleReadFloat16 >(buf,addr,conf);
      }

      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionDouble32(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         return ReadNumericalCollection<Double_t,SimpleReadDouble32 >(buf,addr,conf);
      }

      template <typename T>
      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionBasicType(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         return ReadNumericalCollection<T,SimpleRead<T> >(buf,addr,conf);
      }

      template <typename From, typename To>
      struct ConvertRead {
         static INLINE_TEMPLATE_ARGS void Action(TBuffer &buf, void *addr, Int_t nvalues)
         {
            From *temp = new From[nvalues];
            buf.ReadFastArray(temp, nvalues);
            To *vec = (To*)addr;
            for(Int_t ind = 0; ind < nvalues; ++ind) {
               vec[ind] = (To)temp[ind];
            }
            delete [] temp;
         }
      };

      template <typename From, typename To>
      struct ConvertRead<NoFactorMarker<From>,To> {
         static INLINE_TEMPLATE_ARGS void Action(TBuffer &buf, void *addr, Int_t nvalues)
         {
            From *temp = new From[nvalues];
            buf.ReadFastArrayWithNbits(temp, nvalues,0);
            To *vec = (To*)addr;
            for(Int_t ind = 0; ind < nvalues; ++ind) {
               vec[ind] = (To)temp[ind];
            }
            delete [] temp;
         }
      };

      template <typename From, typename To>
      struct ConvertRead<WithFactorMarker<From>,To> {
         static INLINE_TEMPLATE_ARGS void Action(TBuffer &buf, void *addr, Int_t nvalues)
         {
            From *temp = new From[nvalues];
            double factor,min; // needs to be initialized.
            buf.ReadFastArrayWithFactor(temp, nvalues, factor, min);
            To *vec = (To*)addr;
            for(Int_t ind = 0; ind < nvalues; ++ind) {
               vec[ind] = (To)temp[ind];
            }
            delete [] temp;
         }
      };

      template <typename From, typename To>
      struct ConvertCollectionBasicType {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *addr, const TConfiguration *conf)
         {
            return ReadNumericalCollection<To,ConvertRead<From,To>::Action >(buf,addr,conf);
         }
      };

    };

   struct GenericLooper {

      template <typename T>
      static INLINE_TEMPLATE_ARGS Int_t ReadBasicType(TBuffer &buf, void *start, const void *end, const TLoopConfiguration *loopconf, const TConfiguration *config)
      {
         TGenericLoopConfig *loopconfig = (TGenericLoopConfig*)loopconf;

         Next_t next = loopconfig->fNext;
         const Int_t offset = config->fOffset;

         char iterator[TVirtualCollectionProxy::fgIteratorArenaSize];
         void *iter = loopconfig->fCopyIterator(iterator,start);
         void *addr;
         while( (addr = next(iter,end)) ) {
            T *x =  (T*)( ((char*)addr) + offset );
            buf >> *x;
         }
         if (iter != &iterator[0]) {
            loopconfig->fDeleteIterator(iter);
         }
         return 0;
      }

      template <typename T>
      static INLINE_TEMPLATE_ARGS Int_t WriteBasicType(TBuffer &buf, void *start, const void *end, const TLoopConfiguration *loopconf, const TConfiguration *config)
      {
         TGenericLoopConfig *loopconfig = (TGenericLoopConfig*)loopconf;

         Next_t next = loopconfig->fNext;
         const Int_t offset = config->fOffset;

         char iterator[TVirtualCollectionProxy::fgIteratorArenaSize];
         void *iter = loopconfig->fCopyIterator(iterator,start);
         void *addr;
         while( (addr = next(iter,end)) ) {
            T *x =  (T*)( ((char*)addr) + offset );
            buf << *x;
         }
         if (iter != &iterator[0]) {
            loopconfig->fDeleteIterator(iter);
         }
         return 0;
      }

      template <Int_t (*iter_action)(TBuffer&,void *,const TConfiguration*)>
      static INLINE_TEMPLATE_ARGS Int_t ReadAction(TBuffer &buf, void *start, const void *end, const TLoopConfiguration *loopconf, const TConfiguration *config)
      {
         TGenericLoopConfig *loopconfig = (TGenericLoopConfig*)loopconf;

         // const Int_t offset = config->fOffset;
         Next_t next = loopconfig->fNext;

         char iterator[TVirtualCollectionProxy::fgIteratorArenaSize];
         void *iter = loopconfig->fCopyIterator(&iterator,start);
         void *addr;
         while( (addr = next(iter,end)) ) {
            iter_action(buf, addr, config);
         }
         if (iter != &iterator[0]) {
            loopconfig->fDeleteIterator(iter);
         }
         return 0;
      }

      template <typename From, typename To>
      struct Generic {
         static void ConvertAction(From *items, void *start, const void *end, const TLoopConfiguration *loopconf, const TConfiguration *config)
         {
            TGenericLoopConfig *loopconfig = (TGenericLoopConfig*)loopconf;

            const Int_t offset = config->fOffset;
            Next_t next = loopconfig->fNext;

            char iterator[TVirtualCollectionProxy::fgIteratorArenaSize];
            void *iter = loopconfig->fCopyIterator(&iterator,start);
            void *addr;
            while( (addr = next(iter,end)) ) {
               To *x =  (To*)( ((char*)addr) + offset );
               *x = (To)(*items);
               ++items;
            }
            if (iter != &iterator[0]) {
               loopconfig->fDeleteIterator(iter);
            }
         }
      };

      template <typename From, typename To>
      struct Numeric {
         static void ConvertAction(From *items, void *start, const void *end, const TLoopConfiguration *loopconf, const TConfiguration * /* config */)
         {
            // The difference with ConvertAction is that we can modify the start
            // iterator and skip the copy.  We also never have an offset.

            TGenericLoopConfig *loopconfig = (TGenericLoopConfig*)loopconf;
            Next_t next = loopconfig->fNext;

            void *iter = start;
            void *addr;
            while( (addr = next(iter,end)) ) {
               To *x =  (To*)(addr);
               *x = (To)(*items);
               ++items;
            }
         }
      };

      template <typename From, typename To, template <typename F, typename T> class Converter = Generic >
      struct ConvertBasicType {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *start, const void *end, const TLoopConfiguration *loopconf, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.

            TGenericLoopConfig *loopconfig = (TGenericLoopConfig*)loopconf;
            TVirtualCollectionProxy *proxy = loopconfig->fProxy;
            Int_t nvalues = proxy->Size();

            From *items = new From[nvalues];
            buf.ReadFastArray(items, nvalues);
            Converter<From,To>::ConvertAction(items,start,end,loopconfig,config);
            delete [] items;
            return 0;
         }
      };

      template <typename To>
      struct ConvertBasicType<BitsMarker, To, Generic> {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *start, const void *end, const TLoopConfiguration *loopconf, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.

            TGenericLoopConfig *loopconfig = (TGenericLoopConfig*)loopconf;
            TVirtualCollectionProxy *proxy = loopconfig->fProxy;
            Int_t nvalues = proxy->Size();

            UInt_t *items_storage = new UInt_t[nvalues];
            UInt_t *items = items_storage;

            const Int_t offset = config->fOffset;
            Next_t next = loopconfig->fNext;

            char iterator[TVirtualCollectionProxy::fgIteratorArenaSize];
            void *iter = loopconfig->fCopyIterator(&iterator,start);
            void *addr;
            while( (addr = next(iter,end)) ) {
               buf >> (*items);
               if (((*items) & kIsReferenced) != 0) {
                  HandleReferencedTObject(buf, addr, config);
               }
               To *x =  (To*)( ((char*)addr) + offset );
               *x = (To)(*items);
               ++items;
            }
            if (iter != &iterator[0]) {
               loopconfig->fDeleteIterator(iter);
            }

            delete [] items_storage;
            return 0;
         }
      };

      template <typename From, typename To, template <typename F, typename T> class Converter >
      struct ConvertBasicType<WithFactorMarker<From>,To,Converter > {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *start, const void *end, const TLoopConfiguration *loopconf, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.

            TGenericLoopConfig *loopconfig = (TGenericLoopConfig*)loopconf;
            TVirtualCollectionProxy *proxy = loopconfig->fProxy;
            Int_t nvalues = proxy->Size();

            TConfSTLWithFactor *conf = (TConfSTLWithFactor *)config;

            From *items = new From[nvalues];
            buf.ReadFastArrayWithFactor(items, nvalues, conf->fFactor, conf->fXmin);
            Converter<From,To>::ConvertAction(items,start,end,loopconfig,config);
            delete [] items;
            return 0;
         }
      };

      template <typename From, typename To, template <typename F, typename T> class Converter >
      struct ConvertBasicType<NoFactorMarker<From>,To,Converter > {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *start, const void *end, const TLoopConfiguration *loopconf, const TConfiguration *config)
         {
            // Simple conversion from a 'From' on disk to a 'To' in memory.

            TGenericLoopConfig *loopconfig = (TGenericLoopConfig*)loopconf;
            TVirtualCollectionProxy *proxy = loopconfig->fProxy;
            Int_t nvalues = proxy->Size();

            TConfSTLNoFactor *conf = (TConfSTLNoFactor *)config;

            From *items = new From[nvalues];
            buf.ReadFastArrayWithNbits(items, nvalues, conf->fNbits);
            Converter<From,To>::ConvertAction(items,start,end,loopconfig,config);
            delete [] items;
            return 0;
         }
      };

      static INLINE_TEMPLATE_ARGS Int_t ReadBase(TBuffer &buf, void *start, const void *end, const TLoopConfiguration * loopconfig, const TConfiguration *config)
      {
         // Well the implementation is non trivial since we do not have a proxy for the container of _only_ the base class.  For now
         // punt.

         return GenericRead(buf,start,end,loopconfig, config);
      }

      static INLINE_TEMPLATE_ARGS Int_t GenericRead(TBuffer &buf, void *, const void *, const TLoopConfiguration * loopconf, const TConfiguration *config)
      {
         TGenericLoopConfig *loopconfig = (TGenericLoopConfig*)loopconf;
         TVirtualCollectionProxy *proxy = loopconfig->fProxy;
         return ((TStreamerInfo*)config->fInfo)->ReadBuffer(buf, *proxy, &(config->fCompInfo), /*first*/ 0, /*last*/ 1, /*narr*/ proxy->Size(), config->fOffset, 1|2 );
      }

      static INLINE_TEMPLATE_ARGS Int_t GenericWrite(TBuffer &buf, void *, const void *, const TLoopConfiguration * loopconf, const TConfiguration *config)
      {
         TGenericLoopConfig *loopconfig = (TGenericLoopConfig*)loopconf;
         TVirtualCollectionProxy *proxy = loopconfig->fProxy;
         return ((TStreamerInfo*)config->fInfo)->WriteBufferAux(buf, *proxy, &(config->fCompInfo), /*first*/ 0, /*last*/ 1, proxy->Size(), config->fOffset, 1|2 );
      }

      template <typename T>
      static INLINE_TEMPLATE_ARGS void SimpleRead(TBuffer &buf, void *addr)
      {
         buf >> *(T*)addr;
      }

      static INLINE_TEMPLATE_ARGS void SimpleReadFloat16(TBuffer &buf, void *addr)
      {
         buf.ReadWithNbits((float*)addr,12);
      }

      static INLINE_TEMPLATE_ARGS void SimpleReadDouble32(TBuffer &buf, void *addr)
      {
         //we read a float and convert it to double
         Float_t afloat;
         buf >> afloat;
         *(double*)addr = (Double_t)afloat;
      }

      template <typename ActionHolder>
      static INLINE_TEMPLATE_ARGS Int_t ReadNumericalCollection(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         // Collection of numbers.  Memberwise or not, it is all the same.

         TConfigSTL *config = (TConfigSTL*)conf;
         UInt_t start, count;
         /* Version_t vers = */ buf.ReadVersion(&start, &count, config->fOldClass);

         TClass *newClass = config->fNewClass;
         TVirtualCollectionProxy *newProxy = newClass->GetCollectionProxy();
         TVirtualCollectionProxy::TPushPop helper( newProxy, ((char*)addr)+config->fOffset );

         Int_t nvalues;
         buf.ReadInt(nvalues);
         void* alternative = newProxy->Allocate(nvalues,true);
         if (nvalues) {
            char startbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
            char endbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
            void *begin = &(startbuf[0]);
            void *end = &(endbuf[0]);
            config->fCreateIterators(alternative, &begin, &end, newProxy);
            // We can not get here with a split vector of pointer, so we can indeed assume
            // that actions->fConfiguration != null.

            TGenericLoopConfig loopconf(newProxy, /* read */ kTRUE);
            ActionHolder::Action(buf,begin,end,&loopconf,config);

            if (begin != &(startbuf[0])) {
               // assert(end != endbuf);
               config->fDeleteTwoIterators(begin,end);
            }
         }
         newProxy->Commit(alternative);

         buf.CheckByteCount(start,count,config->fTypeName);
         return 0;
      }

      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionBool(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         return ReadNumericalCollection<ConvertBasicType<bool,bool,Numeric > >(buf,addr,conf);
      }

      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionFloat16(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         return ReadNumericalCollection<ConvertBasicType<NoFactorMarker<float>,float,Numeric > >(buf,addr,conf);
      }

      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionDouble32(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         return ReadNumericalCollection<ConvertBasicType<float,double,Numeric > >(buf,addr,conf);
         // Could also use:
         // return ReadNumericalCollection<ConvertBasicType<NoFactorMarker<double>,double,Numeric > >(buf,addr,conf);
      }

      template <typename T>
      static INLINE_TEMPLATE_ARGS Int_t ReadCollectionBasicType(TBuffer &buf, void *addr, const TConfiguration *conf)
      {
         return ReadNumericalCollection<ConvertBasicType<T,T,Numeric > >(buf,addr,conf);
      }

      template <typename From, typename To>
      struct ConvertCollectionBasicType {
         static INLINE_TEMPLATE_ARGS Int_t Action(TBuffer &buf, void *addr, const TConfiguration *conf)
         {
            // return ReadNumericalCollection<To,ConvertRead<From,To>::Action >(buf,addr,conf);
            return ReadNumericalCollection<ConvertBasicType<From,To,Numeric > >(buf,addr,conf);
         }
      };

   };
}

template <typename Looper, typename From>
static TConfiguredAction GetCollectionReadConvertAction(Int_t newtype, TConfiguration *conf)
{
   switch (newtype) {
      case TStreamerInfo::kBool:    return TConfiguredAction( Looper::template ConvertBasicType<From,bool>::Action, conf ); break;
      case TStreamerInfo::kChar:    return TConfiguredAction( Looper::template ConvertBasicType<From,char>::Action, conf ); break;
      case TStreamerInfo::kShort:   return TConfiguredAction( Looper::template ConvertBasicType<From,short>::Action, conf );  break;
      case TStreamerInfo::kInt:     return TConfiguredAction( Looper::template ConvertBasicType<From,Int_t>::Action, conf ); break;
      case TStreamerInfo::kLong:    return TConfiguredAction( Looper::template ConvertBasicType<From,Long_t>::Action, conf ); break;
      case TStreamerInfo::kLong64:  return TConfiguredAction( Looper::template ConvertBasicType<From,Long64_t>::Action, conf ); break;
      case TStreamerInfo::kFloat:   return TConfiguredAction( Looper::template ConvertBasicType<From,float>::Action, conf ); break;
      case TStreamerInfo::kFloat16: return TConfiguredAction( Looper::template ConvertBasicType<From,float>::Action, conf ); break;
      case TStreamerInfo::kDouble:  return TConfiguredAction( Looper::template ConvertBasicType<From,double>::Action, conf ); break;
      case TStreamerInfo::kDouble32:return TConfiguredAction( Looper::template ConvertBasicType<From,double>::Action, conf ); break;
      case TStreamerInfo::kUChar:   return TConfiguredAction( Looper::template ConvertBasicType<From,UChar_t>::Action, conf ); break;
      case TStreamerInfo::kUShort:  return TConfiguredAction( Looper::template ConvertBasicType<From,UShort_t>::Action, conf ); break;
      case TStreamerInfo::kUInt:    return TConfiguredAction( Looper::template ConvertBasicType<From,UInt_t>::Action, conf ); break;
      case TStreamerInfo::kULong:   return TConfiguredAction( Looper::template ConvertBasicType<From,ULong_t>::Action, conf ); break;
      case TStreamerInfo::kULong64: return TConfiguredAction( Looper::template ConvertBasicType<From,ULong64_t>::Action, conf );  break;
      case TStreamerInfo::kBits:    return TConfiguredAction( Looper::template ConvertBasicType<From,UInt_t>::Action, conf ); break;
      default:
         return TConfiguredAction( Looper::GenericRead, conf );
         break;
   }
   R__ASSERT(0); // We should never be here
   return TConfiguredAction();
}

template <class Looper>
static TConfiguredAction GetNumericCollectionReadAction(Int_t type, TConfigSTL *conf)
{
   // If we ever support std::vector<Double32_t> fValues; //[...] we would get the info from the StreamerElement for fValues.

   switch (type) {
      // Read basic types.
      case TStreamerInfo::kBool:    return TConfiguredAction( Looper::ReadCollectionBool, conf );    break;
      case TStreamerInfo::kChar:    return TConfiguredAction( Looper::template ReadCollectionBasicType<Char_t>, conf );    break;
      case TStreamerInfo::kShort:   return TConfiguredAction( Looper::template ReadCollectionBasicType<Short_t>,conf );   break;
      case TStreamerInfo::kInt:     return TConfiguredAction( Looper::template ReadCollectionBasicType<Int_t>,  conf );     break;
      case TStreamerInfo::kLong:    return TConfiguredAction( Looper::template ReadCollectionBasicType<Long_t>, conf );    break;
      case TStreamerInfo::kLong64:  return TConfiguredAction( Looper::template ReadCollectionBasicType<Long64_t>, conf );  break;
      case TStreamerInfo::kFloat:   return TConfiguredAction( Looper::template ReadCollectionBasicType<Float_t>,  conf );   break;
      case TStreamerInfo::kDouble:  return TConfiguredAction( Looper::template ReadCollectionBasicType<Double_t>, conf );  break;
      case TStreamerInfo::kUChar:   return TConfiguredAction( Looper::template ReadCollectionBasicType<UChar_t>,  conf );   break;
      case TStreamerInfo::kUShort:  return TConfiguredAction( Looper::template ReadCollectionBasicType<UShort_t>, conf );  break;
      case TStreamerInfo::kUInt:    return TConfiguredAction( Looper::template ReadCollectionBasicType<UInt_t>,   conf );    break;
      case TStreamerInfo::kULong:   return TConfiguredAction( Looper::template ReadCollectionBasicType<ULong_t>,  conf );   break;
      case TStreamerInfo::kULong64: return TConfiguredAction( Looper::template ReadCollectionBasicType<ULong64_t>, conf ); break;
      case TStreamerInfo::kBits:    Error("GetNumericCollectionReadAction","There is no support for kBits outside of a TObject."); break;
      case TStreamerInfo::kFloat16: {
         TConfigSTL *alternate = new TConfSTLNoFactor(conf,12);
         delete conf;
         return TConfiguredAction( Looper::ReadCollectionFloat16, alternate );
         // if (element->GetFactor() != 0) {
         //    return TConfiguredAction( Looper::template ReadAction<ReadBasicType_WithFactor<float> >, new TConfWithFactor(info,i,compinfo,offset,element->GetFactor(),element->GetXmin()) );
         // } else {
         //    Int_t nbits = (Int_t)element->GetXmin();
         //    if (!nbits) nbits = 12;
         //    return TConfiguredAction( Looper::template ReadAction<ReadBasicType_NoFactor<float> >, new TConfNoFactor(info,i,compinfo,offset,nbits) );
         // }
         break;
      }
      case TStreamerInfo::kDouble32: {
         TConfigSTL *alternate = new TConfSTLNoFactor(conf,0);
         delete conf;
         return TConfiguredAction( Looper::ReadCollectionDouble32, alternate );
         // if (element->GetFactor() != 0) {
         //    return TConfiguredAction( Looper::template ReadAction<ReadBasicType_WithFactor<double> >, new TConfWithFactor(info,i,compinfo,offset,element->GetFactor(),element->GetXmin()) );
         // } else {
         //    Int_t nbits = (Int_t)element->GetXmin();
         //    if (!nbits) {
         //       return TConfiguredAction( Looper::template ReadAction<ConvertBasicType<float,double> >, new TConfiguration(info,i,compinfo,offset) );
         //    } else {
         //       return TConfiguredAction( Looper::template ReadAction<ReadBasicType_NoFactor<double> >, new TConfNoFactor(info,i,compinfo,offset,nbits) );
         //    }
         // }
         break;
      }
   }
   Fatal("GetNumericCollectionReadAction","Is confused about %d",type);
   R__ASSERT(0); // We should never be here
   return TConfiguredAction();
}

template <typename Looper, typename From>
static TConfiguredAction GetConvertCollectionReadActionFrom(Int_t newtype, TConfiguration *conf)
{
   switch (newtype) {
      case TStreamerInfo::kBool:    return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,bool>::Action, conf ); break;
      case TStreamerInfo::kChar:    return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,char>::Action, conf ); break;
      case TStreamerInfo::kShort:   return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,short>::Action, conf );  break;
      case TStreamerInfo::kInt:     return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,Int_t>::Action, conf ); break;
      case TStreamerInfo::kLong:    return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,Long_t>::Action, conf ); break;
      case TStreamerInfo::kLong64:  return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,Long64_t>::Action, conf ); break;
      case TStreamerInfo::kFloat:   return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,float>::Action, conf ); break;
      case TStreamerInfo::kFloat16: return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,float>::Action, conf ); break;
      case TStreamerInfo::kDouble:  return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,double>::Action, conf ); break;
      case TStreamerInfo::kDouble32:return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,double>::Action, conf ); break;
      case TStreamerInfo::kUChar:   return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,UChar_t>::Action, conf ); break;
      case TStreamerInfo::kUShort:  return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,UShort_t>::Action, conf ); break;
      case TStreamerInfo::kUInt:    return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,UInt_t>::Action, conf ); break;
      case TStreamerInfo::kULong:   return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,ULong_t>::Action, conf ); break;
      case TStreamerInfo::kULong64: return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,ULong64_t>::Action, conf );  break;
      case TStreamerInfo::kBits:    return TConfiguredAction( Looper::template ConvertCollectionBasicType<From,UInt_t>::Action, conf );  break;
      default:
         break;
   }
   R__ASSERT(0); // We should never be here
   return TConfiguredAction();
}

template <typename Looper>
static TConfiguredAction GetConvertCollectionReadAction(Int_t oldtype, Int_t newtype, TConfiguration *conf)
{
   switch (oldtype) {
      case TStreamerInfo::kBool:
         return GetConvertCollectionReadActionFrom<Looper,Bool_t>(newtype, conf );
         break;
      case TStreamerInfo::kChar:
         return GetConvertCollectionReadActionFrom<Looper,Char_t>(newtype, conf );
         break;
      case TStreamerInfo::kShort:
         return GetConvertCollectionReadActionFrom<Looper,Short_t>(newtype, conf );
         break;
      case TStreamerInfo::kInt:
         return GetConvertCollectionReadActionFrom<Looper,Int_t>(newtype, conf );
         break;
      case TStreamerInfo::kLong:
         return GetConvertCollectionReadActionFrom<Looper,Long_t>(newtype, conf );
         break;
      case TStreamerInfo::kLong64:
         return GetConvertCollectionReadActionFrom<Looper,Long64_t>(newtype, conf );
         break;
      case TStreamerInfo::kFloat:
         return GetConvertCollectionReadActionFrom<Looper,Float_t>( newtype, conf );
         break;
      case TStreamerInfo::kDouble:
         return GetConvertCollectionReadActionFrom<Looper,Double_t>(newtype, conf );
         break;
      case TStreamerInfo::kUChar:
         return GetConvertCollectionReadActionFrom<Looper,UChar_t>(newtype, conf );
         break;
      case TStreamerInfo::kUShort:
         return GetConvertCollectionReadActionFrom<Looper,UShort_t>(newtype, conf );
         break;
      case TStreamerInfo::kUInt:
         return GetConvertCollectionReadActionFrom<Looper,UInt_t>(newtype, conf );
         break;
      case TStreamerInfo::kULong:
         return GetConvertCollectionReadActionFrom<Looper,ULong_t>(newtype, conf );
         break;
      case TStreamerInfo::kULong64:
         return GetConvertCollectionReadActionFrom<Looper,ULong64_t>(newtype, conf );
         break;
      case TStreamerInfo::kFloat16:
         return GetConvertCollectionReadActionFrom<Looper,NoFactorMarker<Float16_t> >( newtype, conf );
         break;
      case TStreamerInfo::kDouble32:
         return GetConvertCollectionReadActionFrom<Looper,NoFactorMarker<Double32_t> >( newtype, conf );
         break;
      case TStreamerInfo::kBits:
         Error("GetConvertCollectionReadAction","There is no support for kBits outside of a TObject.");
         break;
      default:
         break;
   }
   R__ASSERT(0); // We should never be here
   return TConfiguredAction();
}

template <class Looper>
static TConfiguredAction GetCollectionReadAction(TVirtualStreamerInfo *info, TStreamerElement *element, Int_t type, UInt_t i, TStreamerInfo::TCompInfo_t *compinfo, Int_t offset)
{
   switch (type) {
      // Read basic types.
      case TStreamerInfo::kBool:    return TConfiguredAction( Looper::template ReadBasicType<Bool_t>, new TConfiguration(info,i,compinfo,offset) );    break;
      case TStreamerInfo::kChar:    return TConfiguredAction( Looper::template ReadBasicType<Char_t>, new TConfiguration(info,i,compinfo,offset) );    break;
      case TStreamerInfo::kShort:   return TConfiguredAction( Looper::template ReadBasicType<Short_t>,new TConfiguration(info,i,compinfo,offset) );   break;
      case TStreamerInfo::kInt:     return TConfiguredAction( Looper::template ReadBasicType<Int_t>,  new TConfiguration(info,i,compinfo,offset) );     break;
      case TStreamerInfo::kLong:    return TConfiguredAction( Looper::template ReadBasicType<Long_t>, new TConfiguration(info,i,compinfo,offset) );    break;
      case TStreamerInfo::kLong64:  return TConfiguredAction( Looper::template ReadBasicType<Long64_t>, new TConfiguration(info,i,compinfo,offset) );  break;
      case TStreamerInfo::kFloat:   return TConfiguredAction( Looper::template ReadBasicType<Float_t>,  new TConfiguration(info,i,compinfo,offset) );   break;
      case TStreamerInfo::kDouble:  return TConfiguredAction( Looper::template ReadBasicType<Double_t>, new TConfiguration(info,i,compinfo,offset) );  break;
      case TStreamerInfo::kUChar:   return TConfiguredAction( Looper::template ReadBasicType<UChar_t>,  new TConfiguration(info,i,compinfo,offset) );   break;
      case TStreamerInfo::kUShort:  return TConfiguredAction( Looper::template ReadBasicType<UShort_t>, new TConfiguration(info,i,compinfo,offset) );  break;
      case TStreamerInfo::kUInt:    return TConfiguredAction( Looper::template ReadBasicType<UInt_t>,   new TConfiguration(info,i,compinfo,offset) );    break;
      case TStreamerInfo::kULong:   return TConfiguredAction( Looper::template ReadBasicType<ULong_t>,  new TConfiguration(info,i,compinfo,offset) );   break;
      case TStreamerInfo::kULong64: return TConfiguredAction( Looper::template ReadBasicType<ULong64_t>, new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kBits: return TConfiguredAction( Looper::template ReadAction<TStreamerInfoActions::ReadBasicType<BitsMarker> > , new TBitsConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kFloat16: {
         if (element->GetFactor() != 0) {
            return TConfiguredAction( Looper::template ReadAction<ReadBasicType_WithFactor<float> >, new TConfWithFactor(info,i,compinfo,offset,element->GetFactor(),element->GetXmin()) );
         } else {
            Int_t nbits = (Int_t)element->GetXmin();
            if (!nbits) nbits = 12;
            return TConfiguredAction( Looper::template ReadAction<ReadBasicType_NoFactor<float> >, new TConfNoFactor(info,i,compinfo,offset,nbits) );
         }
         break;
      }
      case TStreamerInfo::kDouble32: {
         if (element->GetFactor() != 0) {
            return TConfiguredAction( Looper::template ReadAction<ReadBasicType_WithFactor<double> >, new TConfWithFactor(info,i,compinfo,offset,element->GetFactor(),element->GetXmin()) );
         } else {
            Int_t nbits = (Int_t)element->GetXmin();
            if (!nbits) {
               return TConfiguredAction( Looper::template ReadAction<ConvertBasicType<float,double>::Action >, new TConfiguration(info,i,compinfo,offset) );
            } else {
               return TConfiguredAction( Looper::template ReadAction<ReadBasicType_NoFactor<double> >, new TConfNoFactor(info,i,compinfo,offset,nbits) );
            }
         }
         break;
      }
      case TStreamerInfo::kTNamed:  return TConfiguredAction( Looper::template ReadAction<ReadTNamed >, new TConfiguration(info,i,compinfo,offset) );    break;
         // Idea: We should calculate the CanIgnoreTObjectStreamer here and avoid calling the
         // Streamer alltogether.
      case TStreamerInfo::kTObject: return TConfiguredAction( Looper::template ReadAction<ReadTObject >, new TConfiguration(info,i,compinfo,offset) );    break;
      case TStreamerInfo::kTString: return TConfiguredAction( Looper::template ReadAction<ReadTString >, new TConfiguration(info,i,compinfo,offset) );    break;
      case TStreamerInfo::kArtificial:
      case TStreamerInfo::kCacheNew:
      case TStreamerInfo::kCacheDelete:
      case TStreamerInfo::kSTL:  return TConfiguredAction( Looper::GenericRead, new TGenericConfiguration(info,i,compinfo) ); break;
      case TStreamerInfo::kBase: return TConfiguredAction( Looper::ReadBase, new TGenericConfiguration(info,i,compinfo) ); break;

      // Conversions.
      case TStreamerInfo::kConv + TStreamerInfo::kBool:
         return GetCollectionReadConvertAction<Looper,Bool_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kChar:
         return GetCollectionReadConvertAction<Looper,Char_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kShort:
         return GetCollectionReadConvertAction<Looper,Short_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kInt:
         return GetCollectionReadConvertAction<Looper,Int_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kLong:
         return GetCollectionReadConvertAction<Looper,Long_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kLong64:
         return GetCollectionReadConvertAction<Looper,Long64_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kFloat:
         return GetCollectionReadConvertAction<Looper,Float_t>( element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kDouble:
         return GetCollectionReadConvertAction<Looper,Double_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kUChar:
         return GetCollectionReadConvertAction<Looper,UChar_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kUShort:
         return GetCollectionReadConvertAction<Looper,UShort_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kUInt:
         return GetCollectionReadConvertAction<Looper,UInt_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kULong:
         return GetCollectionReadConvertAction<Looper,ULong_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kULong64:
         return GetCollectionReadConvertAction<Looper,ULong64_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kBits:
         return GetCollectionReadConvertAction<Looper,BitsMarker>(element->GetNewType(), new TBitsConfiguration(info,i,compinfo,offset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kFloat16: {
         if (element->GetFactor() != 0) {
            return GetCollectionReadConvertAction<Looper,WithFactorMarker<float> >(element->GetNewType(), new TConfWithFactor(info,i,compinfo,offset,element->GetFactor(),element->GetXmin()) );
         } else {
            Int_t nbits = (Int_t)element->GetXmin();
            if (!nbits) nbits = 12;
            return GetCollectionReadConvertAction<Looper,NoFactorMarker<float> >(element->GetNewType(), new TConfNoFactor(info,i,compinfo,offset,nbits) );
         }
         break;
      }
      case TStreamerInfo::kConv + TStreamerInfo::kDouble32: {
         if (element->GetFactor() != 0) {
            return GetCollectionReadConvertAction<Looper,WithFactorMarker<double> >(element->GetNewType(), new TConfWithFactor(info,i,compinfo,offset,element->GetFactor(),element->GetXmin()) );
         } else {
            Int_t nbits = (Int_t)element->GetXmin();
            if (!nbits) {
               return GetCollectionReadConvertAction<Looper,Float_t>(element->GetNewType(), new TConfiguration(info,i,compinfo,offset) );
            } else {
               return GetCollectionReadConvertAction<Looper,NoFactorMarker<double> >(element->GetNewType(), new TConfNoFactor(info,i,compinfo,offset,nbits) );
            }
         }
         break;
      }
      default:
         return TConfiguredAction( Looper::GenericRead, new TGenericConfiguration(info,i,compinfo) );
         break;
   }
   R__ASSERT(0); // We should never be here
   return TConfiguredAction();
}

template <class Looper>
static TConfiguredAction GetCollectionWriteAction(TVirtualStreamerInfo *info, TStreamerElement * /*element*/, Int_t type, UInt_t i, TStreamerInfo::TCompInfo_t *compinfo, Int_t offset) {
   switch (type) {
      // read basic types
      case TStreamerInfo::kBool:    return TConfiguredAction( Looper::template WriteBasicType<Bool_t>,   new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kChar:    return TConfiguredAction( Looper::template WriteBasicType<Char_t>,   new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kShort:   return TConfiguredAction( Looper::template WriteBasicType<Short_t>,  new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kInt:     return TConfiguredAction( Looper::template WriteBasicType<Int_t>,    new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kLong:    return TConfiguredAction( Looper::template WriteBasicType<Long_t>,   new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kLong64:  return TConfiguredAction( Looper::template WriteBasicType<Long64_t>, new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kFloat:   return TConfiguredAction( Looper::template WriteBasicType<Float_t>,  new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kDouble:  return TConfiguredAction( Looper::template WriteBasicType<Double_t>, new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kUChar:   return TConfiguredAction( Looper::template WriteBasicType<UChar_t>,  new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kUShort:  return TConfiguredAction( Looper::template WriteBasicType<UShort_t>, new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kUInt:    return TConfiguredAction( Looper::template WriteBasicType<UInt_t>,   new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kULong:   return TConfiguredAction( Looper::template WriteBasicType<ULong_t>,  new TConfiguration(info,i,compinfo,offset) ); break;
      case TStreamerInfo::kULong64: return TConfiguredAction( Looper::template WriteBasicType<ULong64_t>,new TConfiguration(info,i,compinfo,offset) ); break;
      // the simple type missing are kBits and kCounter.
      default:
         return TConfiguredAction( Looper::GenericWrite, new TConfiguration(info,i,compinfo,0 /* 0 because we call the legacy code */) );
   }
   R__ASSERT(0); // We should never be here
   return TConfiguredAction();
}


////////////////////////////////////////////////////////////////////////////////
/// loop on the TStreamerElement list
/// regroup members with same type
/// Store predigested information into local arrays. This saves a huge amount
/// of time compared to an explicit iteration on all elements.

void TStreamerInfo::Compile()
{
   if (IsCompiled()) {
      //Error("Compile","can only be called once; this first call generates both the optimized and memberwise actions.");
      return;
   }
   R__LOCKGUARD(gInterpreterMutex);

   // fprintf(stderr,"Running Compile for %s %d %d req=%d,%d\n",GetName(),fClassVersion,fOptimized,CanOptimize(),TestBit(kCannotOptimize));

   // if (IsCompiled() && (!fOptimized || (CanOptimize() && !TestBit(kCannotOptimize)))) return;
   fOptimized = kFALSE;
   fNdata = 0;
   fNfulldata = 0;

   TObjArray* infos = (TObjArray*) gROOT->GetListOfStreamerInfo();
   if (fNumber < 0) {
      ++fgCount;
      fNumber = fgCount;
   }
   if (fNumber >= infos->GetSize()) {
      infos->AddAtAndExpand(this, fNumber);
   } else {
      if (!infos->At(fNumber)) {
         infos->AddAt(this, fNumber);
      }
   }

   assert(fComp == 0 && fCompFull == 0 && fCompOpt == 0);


   Int_t ndata = fElements->GetEntriesFast();


   if (fReadObjectWise) fReadObjectWise->fActions.clear();
   else fReadObjectWise = new TStreamerInfoActions::TActionSequence(this,ndata);

   if (fWriteObjectWise) fWriteObjectWise->fActions.clear();
   else fWriteObjectWise = new TStreamerInfoActions::TActionSequence(this,ndata);

   if (fReadMemberWise) fReadMemberWise->fActions.clear();
   else fReadMemberWise = new TStreamerInfoActions::TActionSequence(this,ndata);

   if (fReadText) fReadText->fActions.clear();
   else fReadText = new TStreamerInfoActions::TActionSequence(this,ndata);

   if (fWriteMemberWise) fWriteMemberWise->fActions.clear();
   else fWriteMemberWise = new TStreamerInfoActions::TActionSequence(this,ndata);

   if (fReadMemberWiseVecPtr) fReadMemberWiseVecPtr->fActions.clear();
   else fReadMemberWiseVecPtr = new TStreamerInfoActions::TActionSequence(this, ndata, kTRUE);

   if (fWriteMemberWiseVecPtr) fWriteMemberWiseVecPtr->fActions.clear();
   else fWriteMemberWiseVecPtr = new TStreamerInfoActions::TActionSequence(this, ndata, kTRUE);

   if (fWriteText) fWriteText->fActions.clear();
   else fWriteText = new TStreamerInfoActions::TActionSequence(this,ndata);

   if (!ndata) {
      // This may be the case for empty classes (e.g., TAtt3D).
      // We still need to properly set the size of emulated classes (i.e. add the virtual table)
      if (fClass->TestBit(TClass::kIsEmulation) && fNVirtualInfoLoc!=0) {
         fSize = sizeof(TStreamerInfo*);
      }
      fComp = new TCompInfo[1];
      fCompFull = new TCompInfo*[1];
      fCompOpt  = new TCompInfo*[1];
      fCompOpt[0] = fCompFull[0] = &(fComp[0]);
      SetIsCompiled();
      return;
   }

   // At most half of the elements can be used to hold optimized versions.
   // We use the bottom to hold the optimized-into elements and the non-optimized elements
   // and the top to hold the original copy of the optimized out elements.
   fNslots = ndata + ndata/2 + 1;
   Int_t optiOut = 0;

   fComp = new TCompInfo[fNslots];
   fCompFull = new TCompInfo*[ndata];
   fCompOpt  = new TCompInfo*[ndata];

   TStreamerElement* element;
   TStreamerElement* previous = 0;
   Int_t keep = -1;
   Int_t i;

   if (!CanOptimize()) {
      SetBit(kCannotOptimize);
   }

   Bool_t isOptimized = kFALSE;
   Bool_t previousOptimized = kFALSE;

   for (i = 0; i < ndata; ++i) {
      element = (TStreamerElement*) fElements->At(i);
      if (!element) {
         break;
      }

      Int_t asize = element->GetSize();
      if (element->GetArrayLength()) {
         asize /= element->GetArrayLength();
      }
      fComp[fNdata].fType = element->GetType();
      fComp[fNdata].fNewType = element->GetNewType();
      fComp[fNdata].fOffset = element->GetOffset();
      fComp[fNdata].fLength = element->GetArrayLength();
      fComp[fNdata].fElem = element;
      fComp[fNdata].fMethod = element->GetMethod();
      fComp[fNdata].fClass = element->GetClassPointer();
      fComp[fNdata].fNewClass = element->GetNewClass();
      fComp[fNdata].fClassName = TString(element->GetTypeName()).Strip(TString::kTrailing, '*');
      fComp[fNdata].fStreamer = element->GetStreamer();

      // try to group consecutive members of the same type
      if (!TestBit(kCannotOptimize)
          && (keep >= 0)
          && (element->GetType() > 0)
          && (element->GetType() < 10)
          && (fComp[fNdata].fType == fComp[fNdata].fNewType)
          && (fComp[keep].fMethod == 0)
          && (element->GetArrayDim() == 0)
          && (fComp[keep].fType < kObject)
          && (fComp[keep].fType != kCharStar) /* do not optimize char* */
          && (element->GetType() == (fComp[keep].fType%kRegrouped))
          && ((element->GetOffset()-fComp[keep].fOffset) == (fComp[keep].fLength)*asize)
          && ((fOldVersion<6) || !previous || /* In version of TStreamerInfo less than 6, the Double32_t were merged even if their annotation (aka factor) were different */
              ((element->GetFactor() == previous->GetFactor())
               && (element->GetXmin() == previous->GetXmin())
               && (element->GetXmax() == previous->GetXmax())
               )
              )
          && (element->TestBit(TStreamerElement::kCache) == previous->TestBit(TStreamerElement::kCache))
          && (element->TestBit(TStreamerElement::kWrite) == previous->TestBit(TStreamerElement::kWrite))
          // kWholeObject and kDoNotDelete do not apply to numerical elements.
          )
      {
         if (!previousOptimized) {
            // The element was not yet optimized we first need to copy it into
            // the set of original copies.
            fComp[fNslots - (++optiOut) ] = fComp[keep];   // Copy the optimized out elements.
            fCompFull[fNfulldata-1] = &(fComp[fNslots - optiOut]); // Reset the pointer in the full list.
         }
         fComp[fNslots - (++optiOut) ] = fComp[fNdata]; // Copy the optimized out elements.
         fCompFull[fNfulldata] = &(fComp[fNslots - optiOut]);

         R__ASSERT( keep < (fNslots - optiOut) );

         if (fComp[keep].fLength == 0) {
            fComp[keep].fLength++;
         }
         fComp[keep].fLength++;
         fComp[keep].fType = element->GetType() + kRegrouped;
         isOptimized = kTRUE;
         previousOptimized = kTRUE;
      } else if (element->GetType() < 0) {

         // -- Deal with an ignored TObject base class.
         // Note: The only allowed negative value here is -1,
         // and signifies that Build() has found a TObject
         // base class and TClass::IgnoreTObjectStreamer() was
         // called.  In this case the compiled version of the
         // elements omits the TObject base class element,
         // which has to be compensated for by TTree::Bronch()
         // when it is making branches for a split object.
         fComp[fNslots - (++optiOut) ] = fComp[fNdata]; // Copy the 'ignored' element.
         fCompFull[fNfulldata] = &(fComp[fNslots - optiOut]);
         keep = -1;
         previousOptimized = kFALSE;

      } else {
         if (fComp[fNdata].fNewType != fComp[fNdata].fType) {
            if (fComp[fNdata].fNewType > 0) {
               if ( (fComp[fNdata].fNewType == kObjectp || fComp[fNdata].fNewType == kAnyp
                     || fComp[fNdata].fNewType == kObject || fComp[fNdata].fNewType == kAny
                     || fComp[fNdata].fNewType == kTObject || fComp[fNdata].fNewType == kTNamed || fComp[fNdata].fNewType == kTString)
                   && (fComp[fNdata].fType == kObjectp || fComp[fNdata].fType == kAnyp
                       || fComp[fNdata].fType == kObject || fComp[fNdata].fType == kAny
                       || fComp[fNdata].fType == kTObject || fComp[fNdata].fType == kTNamed || fComp[fNdata].fType == kTString )
                   ) {
                  fComp[fNdata].fType = fComp[fNdata].fNewType;
               } else if (fComp[fNdata].fType != kCounter) {
                  fComp[fNdata].fType += kConv;
               }
            } else {
               if (fComp[fNdata].fType == kCounter) {
                  Warning("Compile", "Counter %s should not be skipped from class %s", element->GetName(), GetName());
               }
               fComp[fNdata].fType += kSkip;
            }
         }
         fCompOpt[fNdata] = &(fComp[fNdata]);
         fCompFull[fNfulldata] = &(fComp[fNdata]);

         R__ASSERT( fNdata < (fNslots - optiOut) );

         keep = fNdata;
         if (fComp[keep].fLength == 0) {
            fComp[keep].fLength = 1;
         }
         fNdata++;
         previousOptimized = kFALSE;
      }
      // The test 'fMethod[keep] == 0' fails to detect a variable size array
      // if the counter happens to have an offset of zero, so let's explicitly
      // prevent for here.
      if (element->HasCounter()) keep = -1;
      ++fNfulldata;
      previous = element;
   }

   for (i = 0; i < fNdata; ++i) {
      if (!fCompOpt[i]->fElem || fCompOpt[i]->fElem->GetType()< 0) {
         continue;
      }
      AddReadAction(fReadObjectWise, i, fCompOpt[i]);
      AddWriteAction(fWriteObjectWise, i, fCompOpt[i]);
   }
   for (i = 0; i < fNfulldata; ++i) {
      if (!fCompFull[i]->fElem || fCompFull[i]->fElem->GetType()< 0) {
         continue;
      }
      AddReadAction(fReadMemberWise, i, fCompFull[i]);
      AddWriteAction(fWriteMemberWise, i, fCompFull[i]);
      AddReadMemberWiseVecPtrAction(fReadMemberWiseVecPtr, i, fCompFull[i]);
      AddWriteMemberWiseVecPtrAction(fWriteMemberWiseVecPtr, i, fCompFull[i]);

      AddReadTextAction(fReadText, i, fCompFull[i]);
      AddWriteTextAction(fWriteText, i, fCompFull[i]);
   }
   ComputeSize();

   fOptimized = isOptimized;
   SetIsCompiled();

   if (gDebug > 0) {
      ls();
   }
}

template <typename From>
static void AddReadConvertAction(TStreamerInfoActions::TActionSequence *sequence, Int_t newtype, TConfiguration *conf)
{
   switch (newtype) {
      case TStreamerInfo::kBool:    sequence->AddAction( ConvertBasicType<From,bool>::Action,  conf ); break;
      case TStreamerInfo::kChar:    sequence->AddAction( ConvertBasicType<From,char>::Action,  conf ); break;
      case TStreamerInfo::kShort:   sequence->AddAction( ConvertBasicType<From,short>::Action, conf );  break;
      case TStreamerInfo::kInt:     sequence->AddAction( ConvertBasicType<From,Int_t>::Action, conf ); break;
      case TStreamerInfo::kLong:    sequence->AddAction( ConvertBasicType<From,Long_t>::Action,conf ); break;
      case TStreamerInfo::kLong64:  sequence->AddAction( ConvertBasicType<From,Long64_t>::Action, conf ); break;
      case TStreamerInfo::kFloat:   sequence->AddAction( ConvertBasicType<From,float>::Action,    conf ); break;
      case TStreamerInfo::kFloat16: sequence->AddAction( ConvertBasicType<From,float>::Action,    conf ); break;
      case TStreamerInfo::kDouble:  sequence->AddAction( ConvertBasicType<From,double>::Action,   conf ); break;
      case TStreamerInfo::kDouble32:sequence->AddAction( ConvertBasicType<From,double>::Action,   conf ); break;
      case TStreamerInfo::kUChar:   sequence->AddAction( ConvertBasicType<From,UChar_t>::Action,  conf ); break;
      case TStreamerInfo::kUShort:  sequence->AddAction( ConvertBasicType<From,UShort_t>::Action, conf ); break;
      case TStreamerInfo::kUInt:    sequence->AddAction( ConvertBasicType<From,UInt_t>::Action,   conf ); break;
      case TStreamerInfo::kULong:   sequence->AddAction( ConvertBasicType<From,ULong_t>::Action,  conf ); break;
      case TStreamerInfo::kULong64: sequence->AddAction( ConvertBasicType<From,ULong64_t>::Action,conf );  break;
      case TStreamerInfo::kBits:    sequence->AddAction( ConvertBasicType<From,UInt_t>::Action,   conf ); break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add a read action for the given element.

void TStreamerInfo::AddReadAction(TStreamerInfoActions::TActionSequence *readSequence, Int_t i, TStreamerInfo::TCompInfo *compinfo)
{
   TStreamerElement *element = compinfo->fElem;

   if (element->TestBit(TStreamerElement::kWrite)) return;

   switch (compinfo->fType) {
      // read basic types
      case TStreamerInfo::kBool:    readSequence->AddAction( ReadBasicType<Bool_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
      case TStreamerInfo::kChar:    readSequence->AddAction( ReadBasicType<Char_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
      case TStreamerInfo::kShort:   readSequence->AddAction( ReadBasicType<Short_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );   break;
      case TStreamerInfo::kInt:     readSequence->AddAction( ReadBasicType<Int_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );     break;
      case TStreamerInfo::kLong:    readSequence->AddAction( ReadBasicType<Long_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
      case TStreamerInfo::kLong64:  readSequence->AddAction( ReadBasicType<Long64_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );  break;
      case TStreamerInfo::kFloat:   readSequence->AddAction( ReadBasicType<Float_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );   break;
      case TStreamerInfo::kDouble:  readSequence->AddAction( ReadBasicType<Double_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );  break;
      case TStreamerInfo::kUChar:   readSequence->AddAction( ReadBasicType<UChar_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );   break;
      case TStreamerInfo::kUShort:  readSequence->AddAction( ReadBasicType<UShort_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );  break;
      case TStreamerInfo::kUInt:    readSequence->AddAction( ReadBasicType<UInt_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
      case TStreamerInfo::kULong:   readSequence->AddAction( ReadBasicType<ULong_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );   break;
      case TStreamerInfo::kULong64: readSequence->AddAction( ReadBasicType<ULong64_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) ); break;
      case TStreamerInfo::kBits:    readSequence->AddAction( ReadBasicType<BitsMarker>, new TBitsConfiguration(this,i,compinfo,compinfo->fOffset) );     break;
      case TStreamerInfo::kFloat16: {
         if (element->GetFactor() != 0) {
            readSequence->AddAction( ReadBasicType_WithFactor<float>, new TConfWithFactor(this,i,compinfo,compinfo->fOffset,element->GetFactor(),element->GetXmin()) );
         } else {
            Int_t nbits = (Int_t)element->GetXmin();
            if (!nbits) nbits = 12;
            readSequence->AddAction( ReadBasicType_NoFactor<float>, new TConfNoFactor(this,i,compinfo,compinfo->fOffset,nbits) );
         }
         break;
      }
      case TStreamerInfo::kDouble32: {
         if (element->GetFactor() != 0) {
            readSequence->AddAction( ReadBasicType_WithFactor<double>, new TConfWithFactor(this,i,compinfo,compinfo->fOffset,element->GetFactor(),element->GetXmin()) );
         } else {
            Int_t nbits = (Int_t)element->GetXmin();
            if (!nbits) {
               readSequence->AddAction( ConvertBasicType<float,double>::Action, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
            } else {
               readSequence->AddAction( ReadBasicType_NoFactor<double>, new TConfNoFactor(this,i,compinfo,compinfo->fOffset,nbits) );
            }
         }
         break;
      }
      case TStreamerInfo::kTNamed:  readSequence->AddAction( ReadTNamed, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
         // Idea: We should calculate the CanIgnoreTObjectStreamer here and avoid calling the
         // Streamer alltogether.
      case TStreamerInfo::kTObject: readSequence->AddAction( ReadTObject, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
      case TStreamerInfo::kTString: readSequence->AddAction( ReadTString, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
      case TStreamerInfo::kSTL: {
         TClass *newClass = element->GetNewClass();
         TClass *oldClass = element->GetClassPointer();
         Bool_t isSTLbase = element->IsBase() && element->IsA()!=TStreamerBase::Class();

         if (element->GetArrayLength() <= 1) {
            if (fOldVersion<3){   // case of old TStreamerInfo
               if (newClass && newClass != oldClass) {
                  if (element->GetStreamer()) {
                     readSequence->AddAction(ReadSTL<ReadSTLMemberWiseChangedClass,ReadSTLObjectWiseStreamerV2>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,newClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
                  } else {
                     readSequence->AddAction(ReadSTL<ReadSTLMemberWiseChangedClass,ReadSTLObjectWiseFastArrayV2>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,newClass,element->GetTypeName(),isSTLbase));
                  }
               } else {
                  if (element->GetStreamer()) {
                     readSequence->AddAction(ReadSTL<ReadSTLMemberWiseSameClass,ReadSTLObjectWiseStreamerV2>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
                  } else {
                     readSequence->AddAction(ReadSTL<ReadSTLMemberWiseSameClass,ReadSTLObjectWiseFastArrayV2>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,element->GetTypeName(),isSTLbase));
                  }
               }
            } else {
               if (newClass && newClass != oldClass) {
                  if (element->GetStreamer()) {
                     readSequence->AddAction(ReadSTL<ReadSTLMemberWiseChangedClass,ReadSTLObjectWiseStreamer>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,newClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
                  } else {
                     if (oldClass->GetCollectionProxy() == 0 || oldClass->GetCollectionProxy()->GetValueClass() || oldClass->GetCollectionProxy()->HasPointers() ) {
                        readSequence->AddAction(ReadSTL<ReadSTLMemberWiseChangedClass,ReadSTLObjectWiseFastArray>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,newClass,element->GetTypeName(),isSTLbase));
                     } else {
                        switch (SelectLooper(*newClass->GetCollectionProxy())) {
                        case kVectorLooper:
                           readSequence->AddAction(GetConvertCollectionReadAction<VectorLooper>(oldClass->GetCollectionProxy()->GetType(), newClass->GetCollectionProxy()->GetType(), new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,newClass,element->GetTypeName(),isSTLbase)));
                           break;
                        case kAssociativeLooper:
                           readSequence->AddAction(GetConvertCollectionReadAction<AssociativeLooper>(oldClass->GetCollectionProxy()->GetType(), newClass->GetCollectionProxy()->GetType(), new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,newClass,element->GetTypeName(),isSTLbase)));
                           break;
                        case kVectorPtrLooper:
                        case kGenericLooper:
                        default:
                           // For now TBufferXML would force use to allocate the data buffer each time and copy into the real thing.
                           readSequence->AddAction(GetConvertCollectionReadAction<GenericLooper>(oldClass->GetCollectionProxy()->GetType(), newClass->GetCollectionProxy()->GetType(), new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,newClass,element->GetTypeName(),isSTLbase)));
                           break;
                        }
                     }
                  }
               } else {
                  if (element->GetStreamer()) {
                     readSequence->AddAction(ReadSTL<ReadSTLMemberWiseSameClass,ReadSTLObjectWiseStreamer>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
                  } else {
                     if (oldClass->GetCollectionProxy() == 0 || oldClass->GetCollectionProxy()->GetValueClass() || oldClass->GetCollectionProxy()->HasPointers() ) {
                        readSequence->AddAction(ReadSTL<ReadSTLMemberWiseSameClass,ReadSTLObjectWiseFastArray>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,element->GetTypeName(),isSTLbase));
                     } else {
                        switch (SelectLooper(*oldClass->GetCollectionProxy())) {
                        case kVectorLooper:
                           readSequence->AddAction(GetNumericCollectionReadAction<VectorLooper>(oldClass->GetCollectionProxy()->GetType(), new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,element->GetTypeName(),isSTLbase)));
                           break;
                        case kAssociativeLooper:
                           readSequence->AddAction(GetNumericCollectionReadAction<AssociativeLooper>(oldClass->GetCollectionProxy()->GetType(), new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,element->GetTypeName(),isSTLbase)));
                           break;
                        case kVectorPtrLooper:
                        case kGenericLooper:
                        default:
                           // For now TBufferXML would force use to allocate the data buffer each time and copy into the real thing.
                           readSequence->AddAction(GetNumericCollectionReadAction<GenericLooper>(oldClass->GetCollectionProxy()->GetType(), new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,element->GetTypeName(),isSTLbase)));
                           break;
                        }
                     }
                  }
               }
            }
         } else {
            if (fOldVersion<3){   // case of old TStreamerInfo
               if (newClass && newClass != oldClass) {
                  if (element->GetStreamer()) {
                     readSequence->AddAction(ReadSTL<ReadArraySTLMemberWiseChangedClass,ReadSTLObjectWiseStreamerV2>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,newClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
                  } else {
                     readSequence->AddAction(ReadSTL<ReadArraySTLMemberWiseChangedClass,ReadSTLObjectWiseFastArrayV2>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,newClass,element->GetTypeName(),isSTLbase));
                  }
               } else {
                  if (element->GetStreamer()) {
                     readSequence->AddAction(ReadSTL<ReadArraySTLMemberWiseSameClass,ReadSTLObjectWiseStreamerV2>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
                  } else {
                     readSequence->AddAction(ReadSTL<ReadArraySTLMemberWiseSameClass,ReadSTLObjectWiseFastArrayV2>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,element->GetTypeName(),isSTLbase));
                  }
               }
            } else {
               if (newClass && newClass != oldClass) {
                  if (element->GetStreamer()) {
                     readSequence->AddAction(ReadSTL<ReadArraySTLMemberWiseChangedClass,ReadSTLObjectWiseStreamer>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,newClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
                  } else {
                     readSequence->AddAction(ReadSTL<ReadArraySTLMemberWiseChangedClass,ReadSTLObjectWiseFastArray>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,newClass,element->GetTypeName(),isSTLbase));
                  }
               } else {
                  if (element->GetStreamer()) {
                     readSequence->AddAction(ReadSTL<ReadArraySTLMemberWiseSameClass,ReadSTLObjectWiseStreamer>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
                  } else {
                     readSequence->AddAction(ReadSTL<ReadArraySTLMemberWiseSameClass,ReadSTLObjectWiseFastArray>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,element->GetTypeName(),isSTLbase));
                  }
               }
            }
         }
         break;
      }

      case TStreamerInfo::kConv + TStreamerInfo::kBool:
         AddReadConvertAction<Bool_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kChar:
         AddReadConvertAction<Char_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kShort:
         AddReadConvertAction<Short_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kInt:
         AddReadConvertAction<Int_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kLong:
         AddReadConvertAction<Long_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kLong64:
         AddReadConvertAction<Long64_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kFloat:
         AddReadConvertAction<Float_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kDouble:
         AddReadConvertAction<Double_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kUChar:
         AddReadConvertAction<UChar_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kUShort:
         AddReadConvertAction<UShort_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kUInt:
         AddReadConvertAction<UInt_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kULong:
         AddReadConvertAction<ULong_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kULong64:
         AddReadConvertAction<ULong64_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kBits:
         AddReadConvertAction<BitsMarker>(readSequence, compinfo->fNewType, new TBitsConfiguration(this,i,compinfo,compinfo->fOffset) );
         break;
      case TStreamerInfo::kConv + TStreamerInfo::kFloat16: {
         if (element->GetFactor() != 0) {
            AddReadConvertAction<WithFactorMarker<float> >(readSequence, compinfo->fNewType, new TConfWithFactor(this,i,compinfo,compinfo->fOffset,element->GetFactor(),element->GetXmin()) );
         } else {
            Int_t nbits = (Int_t)element->GetXmin();
            if (!nbits) nbits = 12;
            AddReadConvertAction<NoFactorMarker<float> >(readSequence, compinfo->fNewType, new TConfNoFactor(this,i,compinfo,compinfo->fOffset,nbits) );
         }
         break;
      }
      case TStreamerInfo::kConv + TStreamerInfo::kDouble32: {
         if (element->GetFactor() != 0) {
            AddReadConvertAction<WithFactorMarker<double> >(readSequence, compinfo->fNewType, new TConfWithFactor(this,i,compinfo,compinfo->fOffset,element->GetFactor(),element->GetXmin()) );
         } else {
            Int_t nbits = (Int_t)element->GetXmin();
            if (!nbits) {
               AddReadConvertAction<Float_t>(readSequence, compinfo->fNewType, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
            } else {
               AddReadConvertAction<NoFactorMarker<double> >(readSequence, compinfo->fNewType, new TConfNoFactor(this,i,compinfo,compinfo->fOffset,nbits) );
            }
         }
         break;
      }
      default:
         readSequence->AddAction( GenericReadAction, new TGenericConfiguration(this,i,compinfo) );
         break;
   }
   if (element->TestBit(TStreamerElement::kCache)) {
      TConfiguredAction action( readSequence->fActions.back() );  // Action is moved, we must pop it next.
      readSequence->fActions.pop_back();
      readSequence->AddAction( UseCache, new TConfigurationUseCache(this,action,element->TestBit(TStreamerElement::kRepeat)) );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add a read text action for the given element.

void TStreamerInfo::AddReadTextAction(TStreamerInfoActions::TActionSequence *readSequence, Int_t i, TStreamerInfo::TCompInfo *compinfo)
{
   TStreamerElement *element = compinfo->fElem;

   if (element->TestBit(TStreamerElement::kWrite))
      return;

   Bool_t generic = kFALSE, isBase = kFALSE;

   switch (compinfo->fType) {
   case TStreamerInfo::kTObject:
      if (element->IsBase())
         isBase = kTRUE;
      // readSequence->AddAction( ReadTextTObjectBase, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
      else
         readSequence->AddAction(ReadTextTObject, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;

   case TStreamerInfo::kTNamed:
      if (element->IsBase())
         isBase = kTRUE;
      // generic = kTRUE; // for the base class one cannot call TClass::Streamer() as performed for the normal object
      else
         readSequence->AddAction(ReadTextTNamed, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;

   case TStreamerInfo::kObject: // Class      derived from TObject
   case TStreamerInfo::kAny:    // Class  NOT derived from TObject
   case TStreamerInfo::kOffsetL + TStreamerInfo::kObject:
   case TStreamerInfo::kAny + TStreamerInfo::kOffsetL:
      readSequence->AddAction(ReadTextObject, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;

   case TStreamerInfo::kSTLp: // Pointer to container with no virtual table (stl) and no comment
   case TStreamerInfo::kSTLp +
      TStreamerInfo::kOffsetL: // array of pointers to container with no virtual table (stl) and no comment
      readSequence->AddAction(ReadSTLp<true>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;

   case TStreamerInfo::kStreamLoop:
   case TStreamerInfo::kOffsetL + TStreamerInfo::kStreamLoop:
      readSequence->AddAction(ReadStreamerLoop<true>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;

   case TStreamerInfo::kBase: isBase = kTRUE; break;

   case TStreamerInfo::kStreamer:
      readSequence->AddAction(ReadTextStreamer, new TGenericConfiguration(this, i, compinfo));
      break;

   default: generic = kTRUE; break;
   }

   if (isBase) {
      if (compinfo->fStreamer) {
         readSequence->AddAction(ReadTextStreamer, new TGenericConfiguration(this, i, compinfo));
      } else {
         readSequence->AddAction(ReadTextBaseClass, new TGenericConfiguration(this, i, compinfo));
      }
   } else if (generic)
      readSequence->AddAction(GenericReadAction, new TGenericConfiguration(this, i, compinfo));
}

////////////////////////////////////////////////////////////////////////////////
/// Add a read action for the given element.
/// This is for streaming via a TClonesArray (or a vector of pointers of this type).

void TStreamerInfo::AddReadMemberWiseVecPtrAction(TStreamerInfoActions::TActionSequence *readSequence, Int_t i, TStreamerInfo::TCompInfo *compinfo)
{
   TStreamerElement *element = compinfo->fElem;

   if (element->TestBit(TStreamerElement::kWrite)) return;

   if (element->TestBit(TStreamerElement::kCache)) {
      TConfiguredAction action( GetCollectionReadAction<VectorLooper>(this,element,compinfo->fType,i,compinfo,compinfo->fOffset) );
      readSequence->AddAction( UseCacheVectorPtrLoop, new TConfigurationUseCache(this,action,element->TestBit(TStreamerElement::kRepeat)) );
   } else {
      readSequence->AddAction( GetCollectionReadAction<VectorPtrLooper>(this,element,compinfo->fType,i,compinfo,compinfo->fOffset) );
   }
}

////////////////////////////////////////////////////////////////////////////////

void TStreamerInfo::AddWriteAction(TStreamerInfoActions::TActionSequence *writeSequence, Int_t i, TStreamerInfo::TCompInfo *compinfo)
{
   TStreamerElement *element = compinfo->fElem;
   if (element->TestBit(TStreamerElement::kCache) && !element->TestBit(TStreamerElement::kWrite)) {
      // Skip element cached for reading purposes.
      return;
   }
   if (element->GetType() >= kArtificial &&  !element->TestBit(TStreamerElement::kWrite)) {
      // Skip artificial element used for reading purposes.
      return;
   }
   switch (compinfo->fType) {
      // write basic types
      case TStreamerInfo::kBool:    writeSequence->AddAction( WriteBasicType<Bool_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
      case TStreamerInfo::kChar:    writeSequence->AddAction( WriteBasicType<Char_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
      case TStreamerInfo::kShort:   writeSequence->AddAction( WriteBasicType<Short_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );   break;
      case TStreamerInfo::kInt:     writeSequence->AddAction( WriteBasicType<Int_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );     break;
      case TStreamerInfo::kLong:    writeSequence->AddAction( WriteBasicType<Long_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
      case TStreamerInfo::kLong64:  writeSequence->AddAction( WriteBasicType<Long64_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );  break;
      case TStreamerInfo::kFloat:   writeSequence->AddAction( WriteBasicType<Float_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );   break;
      case TStreamerInfo::kDouble:  writeSequence->AddAction( WriteBasicType<Double_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );  break;
      case TStreamerInfo::kUChar:   writeSequence->AddAction( WriteBasicType<UChar_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );   break;
      case TStreamerInfo::kUShort:  writeSequence->AddAction( WriteBasicType<UShort_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );  break;
      case TStreamerInfo::kUInt:    writeSequence->AddAction( WriteBasicType<UInt_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
      case TStreamerInfo::kULong:   writeSequence->AddAction( WriteBasicType<ULong_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );   break;
      case TStreamerInfo::kULong64: writeSequence->AddAction( WriteBasicType<ULong64_t>, new TConfiguration(this,i,compinfo,compinfo->fOffset) ); break;
       // case TStreamerInfo::kBits:    writeSequence->AddAction( WriteBasicType<BitsMarker>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
     /*case TStreamerInfo::kFloat16: {
         if (element->GetFactor() != 0) {
            writeSequence->AddAction( WriteBasicType_WithFactor<float>, new TConfWithFactor(this,i,compinfo,compinfo->fOffset,element->GetFactor(),element->GetXmin()) );
         } else {
            Int_t nbits = (Int_t)element->GetXmin();
            if (!nbits) nbits = 12;
            writeSequence->AddAction( WriteBasicType_NoFactor<float>, new TConfNoFactor(this,i,compinfo,compinfo->fOffset,nbits) );
         }
         break;
      } */
     /*case TStreamerInfo::kDouble32: {
        if (element->GetFactor() != 0) {
           writeSequence->AddAction( WriteBasicType_WithFactor<double>, new TConfWithFactor(this,i,compinfo,compinfo->fOffset,element->GetFactor(),element->GetXmin()) );
        } else {
           Int_t nbits = (Int_t)element->GetXmin();
           if (!nbits) {
              writeSequence->AddAction( ConvertBasicType<float,double>, new TConfiguration(this,i,compinfo,compinfo->fOffset) );
           } else {
              writeSequence->AddAction( WriteBasicType_NoFactor<double>, new TConfNoFactor(this,i,compinfo,compinfo->fOffset,nbits) );
           }
        }
        break;
     } */
     //case TStreamerInfo::kTNamed:  writeSequence->AddAction( WriteTNamed, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
        // Idea: We should calculate the CanIgnoreTObjectStreamer here and avoid calling the
        // Streamer alltogether.
     //case TStreamerInfo::kTObject: writeSequence->AddAction( WriteTObject, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
     //case TStreamerInfo::kTString: writeSequence->AddAction( WriteTString, new TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
     /*case TStreamerInfo::kSTL: {
        TClass *newClass = element->GetNewClass();
        TClass *oldClass = element->GetClassPointer();
        Bool_t isSTLbase = element->IsBase() && element->IsA()!=TStreamerBase::Class();

        if (element->GetArrayLength() <= 1) {
           if (newClass && newClass != oldClass) {
              if (element->GetStreamer()) {
                 writeSequence->AddAction(WriteSTL<WriteSTLMemberWiseChangedClass,WriteSTLObjectWiseStreamer>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,newClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
              } else {
                 writeSequence->AddAction(WriteSTL<WriteSTLMemberWiseChangedClass,WriteSTLObjectWiseFastArray>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,newClass,element->GetTypeName(),isSTLbase));
              }
           } else {
              if (element->GetStreamer()) {
                 writeSequence->AddAction(WriteSTL<WriteSTLMemberWiseSameClass,WriteSTLObjectWiseStreamer>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
              } else {
                 writeSequence->AddAction(WriteSTL<WriteSTLMemberWiseSameClass,WriteSTLObjectWiseFastArray>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,element->GetTypeName(),isSTLbase));
              }
           }
        } else {
           if (newClass && newClass != oldClass) {
              if (element->GetStreamer()) {
                 writeSequence->AddAction(WriteSTL<WriteArraySTLMemberWiseChangedClass,WriteSTLObjectWiseStreamer>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,newClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
              } else {
                 writeSequence->AddAction(WriteSTL<WriteArraySTLMemberWiseChangedClass,WriteSTLObjectWiseFastArray>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,newClass,element->GetTypeName(),isSTLbase));
              }
           } else {
              if (element->GetStreamer()) {
                 writeSequence->AddAction(WriteSTL<WriteArraySTLMemberWiseSameClass,WriteSTLObjectWiseStreamer>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
              } else {
                 writeSequence->AddAction(WriteSTL<WriteArraySTLMemberWiseSameClass,WriteSTLObjectWiseFastArray>, new TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,element->GetTypeName(),isSTLbase));
              }
           }
        }
        break;
     } */
      default:
         writeSequence->AddAction( GenericWriteAction, new TGenericConfiguration(this,i,compinfo) );
         break;
   }
#if defined(CDJ_NO_COMPILE)
   if (element->TestBit(TStreamerElement::kCache)) {
      TConfiguredAction action( writeSequence->fActions.back() );  // Action is moved, we must pop it next.
      writeSequence->fActions.pop_back();
      writeSequence->AddAction( UseCache, new TConfigurationUseCache(this,action,element->TestBit(TStreamerElement::kRepeat)) );
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////

void TStreamerInfo::AddWriteTextAction(TStreamerInfoActions::TActionSequence *writeSequence, Int_t i, TStreamerInfo::TCompInfo *compinfo)
{
   TStreamerElement *element = compinfo->fElem;
   if (element->TestBit(TStreamerElement::kCache) && !element->TestBit(TStreamerElement::kWrite)) {
      // Skip element cached for reading purposes.
      return;
   }
   if (element->GetType() >= kArtificial && !element->TestBit(TStreamerElement::kWrite)) {
      // Skip artificial element used for reading purposes.
      return;
   }

   Bool_t generic = kFALSE, isBase = kFALSE;

   switch (compinfo->fType) {
   // write basic types
   case TStreamerInfo::kBool:
      writeSequence->AddAction(WriteBasicType<Bool_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kChar:
      writeSequence->AddAction(WriteBasicType<Char_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kShort:
      writeSequence->AddAction(WriteBasicType<Short_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kInt:
      writeSequence->AddAction(WriteBasicType<Int_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kLong:
      writeSequence->AddAction(WriteBasicType<Long_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kLong64:
      writeSequence->AddAction(WriteBasicType<Long64_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kFloat:
      writeSequence->AddAction(WriteBasicType<Float_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kDouble:
      writeSequence->AddAction(WriteBasicType<Double_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kUChar:
      writeSequence->AddAction(WriteBasicType<UChar_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kUShort:
      writeSequence->AddAction(WriteBasicType<UShort_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kUInt:
      writeSequence->AddAction(WriteBasicType<UInt_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kULong:
      writeSequence->AddAction(WriteBasicType<ULong_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;
   case TStreamerInfo::kULong64:
      writeSequence->AddAction(WriteBasicType<ULong64_t>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;

   case TStreamerInfo::kTObject:
      if (element->IsBase())
         isBase = kTRUE;
      else
         writeSequence->AddAction(WriteTextTObject, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;

   case TStreamerInfo::kTNamed:
      if (element->IsBase())
         isBase = kTRUE;
      else
         writeSequence->AddAction(WriteTextTNamed, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;

   case TStreamerInfo::kSTLp: // Pointer to container with no virtual table (stl) and no comment
   case TStreamerInfo::kSTLp +
      TStreamerInfo::kOffsetL: // array of pointers to container with no virtual table (stl) and no comment
      writeSequence->AddAction(WriteSTLp<true>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;

   case TStreamerInfo::kStreamLoop:
   case TStreamerInfo::kOffsetL + TStreamerInfo::kStreamLoop:
      writeSequence->AddAction(WriteStreamerLoop<true>, new TConfiguration(this, i, compinfo, compinfo->fOffset));
      break;

   case TStreamerInfo::kBase: isBase = kTRUE; break;

   case TStreamerInfo::kStreamer:
      writeSequence->AddAction(WriteTextStreamer, new TGenericConfiguration(this, i, compinfo));
      break;

   // case TStreamerInfo::kBits:    writeSequence->AddAction( WriteBasicType<BitsMarker>, new
   // TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
   /*case TStreamerInfo::kFloat16: {
       if (element->GetFactor() != 0) {
          writeSequence->AddAction( WriteBasicType_WithFactor<float>, new
    TConfWithFactor(this,i,compinfo,compinfo->fOffset,element->GetFactor(),element->GetXmin()) );
       } else {
          Int_t nbits = (Int_t)element->GetXmin();
          if (!nbits) nbits = 12;
          writeSequence->AddAction( WriteBasicType_NoFactor<float>, new
    TConfNoFactor(this,i,compinfo,compinfo->fOffset,nbits) );
       }
       break;
    } */
   /*case TStreamerInfo::kDouble32: {
      if (element->GetFactor() != 0) {
         writeSequence->AddAction( WriteBasicType_WithFactor<double>, new
   TConfWithFactor(this,i,compinfo,compinfo->fOffset,element->GetFactor(),element->GetXmin()) );
      } else {
         Int_t nbits = (Int_t)element->GetXmin();
         if (!nbits) {
            writeSequence->AddAction( ConvertBasicType<float,double>, new
   TConfiguration(this,i,compinfo,compinfo->fOffset) );
         } else {
            writeSequence->AddAction( WriteBasicType_NoFactor<double>, new
   TConfNoFactor(this,i,compinfo,compinfo->fOffset,nbits) );
         }
      }
      break;
   } */
   // case TStreamerInfo::kTNamed:  writeSequence->AddAction( WriteTNamed, new
   // TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
   // Idea: We should calculate the CanIgnoreTObjectStreamer here and avoid calling the
   // Streamer alltogether.
   // case TStreamerInfo::kTObject: writeSequence->AddAction( WriteTObject, new
   // TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
   // case TStreamerInfo::kTString: writeSequence->AddAction( WriteTString, new
   // TConfiguration(this,i,compinfo,compinfo->fOffset) );    break;
   /*case TStreamerInfo::kSTL: {
      TClass *newClass = element->GetNewClass();
      TClass *oldClass = element->GetClassPointer();
      Bool_t isSTLbase = element->IsBase() && element->IsA()!=TStreamerBase::Class();

      if (element->GetArrayLength() <= 1) {
         if (newClass && newClass != oldClass) {
            if (element->GetStreamer()) {
               writeSequence->AddAction(WriteSTL<WriteSTLMemberWiseChangedClass,WriteSTLObjectWiseStreamer>, new
   TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,newClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
            } else {
               writeSequence->AddAction(WriteSTL<WriteSTLMemberWiseChangedClass,WriteSTLObjectWiseFastArray>, new
   TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,newClass,element->GetTypeName(),isSTLbase));
            }
         } else {
            if (element->GetStreamer()) {
               writeSequence->AddAction(WriteSTL<WriteSTLMemberWiseSameClass,WriteSTLObjectWiseStreamer>, new
   TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
            } else {
               writeSequence->AddAction(WriteSTL<WriteSTLMemberWiseSameClass,WriteSTLObjectWiseFastArray>, new
   TConfigSTL(this,i,compinfo,compinfo->fOffset,1,oldClass,element->GetTypeName(),isSTLbase));
            }
         }
      } else {
         if (newClass && newClass != oldClass) {
            if (element->GetStreamer()) {
               writeSequence->AddAction(WriteSTL<WriteArraySTLMemberWiseChangedClass,WriteSTLObjectWiseStreamer>, new
   TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,newClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
            } else {
               writeSequence->AddAction(WriteSTL<WriteArraySTLMemberWiseChangedClass,WriteSTLObjectWiseFastArray>, new
   TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,newClass,element->GetTypeName(),isSTLbase));
            }
         } else {
            if (element->GetStreamer()) {
               writeSequence->AddAction(WriteSTL<WriteArraySTLMemberWiseSameClass,WriteSTLObjectWiseStreamer>, new
   TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,element->GetStreamer(),element->GetTypeName(),isSTLbase));
            } else {
               writeSequence->AddAction(WriteSTL<WriteArraySTLMemberWiseSameClass,WriteSTLObjectWiseFastArray>, new
   TConfigSTL(this,i,compinfo,compinfo->fOffset,element->GetArrayLength(),oldClass,element->GetTypeName(),isSTLbase));
            }
         }
      }
      break;
   } */
   default: generic = kTRUE; break;
   }

   if (isBase) {
      if (compinfo->fStreamer) {
         writeSequence->AddAction(WriteTextStreamer, new TGenericConfiguration(this, i, compinfo));
      } else {
         writeSequence->AddAction(WriteTextBaseClass, new TGenericConfiguration(this, i, compinfo));
      }

   } else

      // use generic write action when special handling is not provided
      if (generic)
      writeSequence->AddAction(GenericWriteAction, new TGenericConfiguration(this, i, compinfo));

#if defined(CDJ_NO_COMPILE)
   if (element->TestBit(TStreamerElement::kCache)) {
      TConfiguredAction action(writeSequence->fActions.back()); // Action is moved, we must pop it next.
      writeSequence->fActions.pop_back();
      writeSequence->AddAction(UseCache,
                               new TConfigurationUseCache(this, action, element->TestBit(TStreamerElement::kRepeat)));
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// This is for streaming via a TClonesArray (or a vector of pointers of this type).

void TStreamerInfo::AddWriteMemberWiseVecPtrAction(TStreamerInfoActions::TActionSequence *writeSequence, Int_t i, TStreamerInfo::TCompInfo *compinfo)
{
   TStreamerElement *element = compinfo->fElem;
   if (element->TestBit(TStreamerElement::kCache) && !element->TestBit(TStreamerElement::kWrite)) {
      // Skip element cached for reading purposes.
      return;
   }
   if (element->GetType() >= kArtificial &&  !element->TestBit(TStreamerElement::kWrite)) {
      // Skip artificial element used for reading purposes.
      return;
   }

#if defined(CDJ_NO_COMPILE)
   if (element->TestBit(TStreamerElement::kCache)) {
      TConfiguredAction action( GetCollectionWriteAction<VectorLooper>(this,element,compinfo->fType,i,compinfo,compinfo->fOffset) );
      writeSequence->AddAction( UseCacheVectorPtrLoop, new TConfigurationUseCache(this,action,element->TestBit(TStreamerElement::kRepeat)) );
   } else {
      writeSequence->Addaction( GetCollectionWriteAction<VectorPtrLooper>(this,element,compinfo->fType,i,compinfo,compinfo->fOffset) );
   }
#else
   writeSequence->AddAction( VectorPtrLooper::GenericWrite, new TGenericConfiguration(this,i,compinfo) );
#endif

}

////////////////////////////////////////////////////////////////////////////////
/// Create the bundle of the actions necessary for the streaming memberwise of the content described by 'info' into the collection described by 'proxy'

TStreamerInfoActions::TActionSequence *TStreamerInfoActions::TActionSequence::CreateReadMemberWiseActions(TVirtualStreamerInfo *info, TVirtualCollectionProxy &proxy)
{
   if (info == 0) {
      return new TStreamerInfoActions::TActionSequence(0,0);
   }

   TStreamerInfo *sinfo = static_cast<TStreamerInfo*>(info);

   UInt_t ndata = info->GetElements()->GetEntriesFast();
   TStreamerInfoActions::TActionSequence *sequence = new TStreamerInfoActions::TActionSequence(info,ndata);
   if (IsDefaultVector(proxy))
   {
      if (proxy.HasPointers()) {
         // Instead of the creating a new one let's copy the one from the StreamerInfo.
         delete sequence;

         sequence = sinfo->GetReadMemberWiseActions(kTRUE)->CreateCopy();

         return sequence;
      }

      // We can speed up the iteration in case of vector.  We also know that all emulated collection are stored internally as a vector.
      Long_t increment = proxy.GetIncrement();
      sequence->fLoopConfig = new TVectorLoopConfig(&proxy, increment, /* read */ kTRUE);
   } else if (proxy.GetCollectionType() == ROOT::kSTLset || proxy.GetCollectionType() == ROOT::kSTLunorderedset
              || proxy.GetCollectionType() == ROOT::kSTLmultiset || proxy.GetCollectionType() == ROOT::kSTLunorderedmultiset
              || proxy.GetCollectionType() == ROOT::kSTLmap || proxy.GetCollectionType() == ROOT::kSTLmultimap
              || proxy.GetCollectionType() == ROOT::kSTLunorderedmap || proxy.GetCollectionType() == ROOT::kSTLunorderedmultimap)
   {
      Long_t increment = proxy.GetIncrement();
      sequence->fLoopConfig = new TVectorLoopConfig(&proxy, increment, /* read */ kTRUE);
      // sequence->fLoopConfig = new TAssocLoopConfig(proxy);
   } else {
      sequence->fLoopConfig = new TGenericLoopConfig(&proxy, /* read */ kTRUE);
   }
   for (UInt_t i = 0; i < ndata; ++i) {
      TStreamerElement *element = (TStreamerElement*) info->GetElements()->At(i);
      if (!element) {
         break;
      }
      if (element->GetType() < 0) {
         // -- Skip an ignored TObject base class.
         // Note: The only allowed negative value here is -1, and signifies that Build() has found a TObject
         // base class and TClass::IgnoreTObjectStreamer() was called.  In this case the compiled version of the
         // elements omits the TObject base class element, which has to be compensated for by TTree::Bronch()
         // when it is making branches for a split object.
         continue;
      }
      if (element->TestBit(TStreamerElement::kWrite)) {
         // Skip element that only for writing.
         continue;
      }
      TStreamerBase *baseEl = dynamic_cast<TStreamerBase*>(element);
      if (baseEl) {
         if (!baseEl->TestBit(TStreamerElement::kWarned) && baseEl->GetErrorMessage()[0]) {
            // There was a problem with the checksum, the user likely did not
            // increment the version number of the derived class when the
            // base class changed.  Since we will be member wise streaming
            // this class, let's warn the user that something is wrong.
            ::Warning("CreateReadMemberWiseActions","%s",
                      baseEl->GetErrorMessage());
            baseEl->SetBit(TStreamerElement::kWarned);
         }
      }

      TStreamerInfo::TCompInfo_t *compinfo = sinfo->fCompFull[i];

      Int_t asize = element->GetSize();
      if (element->GetArrayLength()) {
         asize /= element->GetArrayLength();
      }
      Int_t oldType = element->GetType();
      Int_t newType = element->GetNewType();

      Int_t offset = element->GetOffset();
      if (newType != oldType) {
         if (newType > 0) {
            if (oldType != TVirtualStreamerInfo::kCounter) {
               oldType += TVirtualStreamerInfo::kConv;
            }
         } else {
            oldType += TVirtualStreamerInfo::kSkip;
         }
      }
      switch (SelectLooper(proxy)) {
      case kAssociativeLooper:
//         } else if (proxy.GetCollectionType() == ROOT::kSTLset || proxy.GetCollectionType() == ROOT::kSTLmultiset
//                    || proxy.GetCollectionType() == ROOT::kSTLmap || proxy.GetCollectionType() == ROOT::kSTLmultimap) {
//            sequence->AddAction( GenericAssocCollectionAction, new TConfigSTL(info,i,compinfo,offset,0,proxy.GetCollectionClass(),0,0) );
      case kVectorLooper:
      case kVectorPtrLooper:
         // We can speed up the iteration in case of vector.  We also know that all emulated collection are stored internally as a vector.
         if (element->TestBit(TStreamerElement::kCache)) {
            TConfiguredAction action( GetCollectionReadAction<VectorLooper>(info,element,oldType,i,compinfo,offset) );
            sequence->AddAction( UseCacheVectorLoop,  new TConfigurationUseCache(info,action,element->TestBit(TStreamerElement::kRepeat)) );
         } else {
            sequence->AddAction( GetCollectionReadAction<VectorLooper>(info,element,oldType,i,compinfo,offset));
         }
         break;
      case kGenericLooper:
      default:
         // The usual collection case.
         if (element->TestBit(TStreamerElement::kCache)) {
            TConfiguredAction action( GetCollectionReadAction<VectorLooper>(info,element,oldType,i,compinfo,offset) );
            sequence->AddAction( UseCacheGenericCollection, new TConfigurationUseCache(info,action,element->TestBit(TStreamerElement::kRepeat)) );
         } else {
            sequence->AddAction( GetCollectionReadAction<GenericLooper>(info,element,oldType,i,compinfo,offset) );
         }
         break;
      }
   }
   return sequence;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the bundle of the actions necessary for the streaming memberwise of the content described by 'info' into the collection described by 'proxy'

TStreamerInfoActions::TActionSequence *TStreamerInfoActions::TActionSequence::CreateWriteMemberWiseActions(TVirtualStreamerInfo *info, TVirtualCollectionProxy &proxy)
{
      if (info == 0) {
         return new TStreamerInfoActions::TActionSequence(0,0);
      }

      UInt_t ndata = info->GetElements()->GetEntriesFast();
      TStreamerInfo *sinfo = static_cast<TStreamerInfo*>(info);
      TStreamerInfoActions::TActionSequence *sequence = new TStreamerInfoActions::TActionSequence(info,ndata);

      if (IsDefaultVector(proxy))
      {
         if (proxy.HasPointers()) {
            // Instead of the creating a new one let's copy the one from the StreamerInfo.
            delete sequence;

            sequence = sinfo->GetWriteMemberWiseActions(kTRUE)->CreateCopy();

            return sequence;
         }

         // We can speed up the iteration in case of vector.  We also know that all emulated collection are stored internally as a vector.
         Long_t increment = proxy.GetIncrement();
         sequence->fLoopConfig = new TVectorLoopConfig(&proxy, increment, /* read */ kFALSE);
      /*} else if (proxy.GetCollectionType() == ROOT::kSTLset || proxy.GetCollectionType() == ROOT::kSTLmultiset
                 || proxy.GetCollectionType() == ROOT::kSTLmap || proxy.GetCollectionType() == ROOT::kSTLmultimap)
      {
         Long_t increment = proxy.GetIncrement();
         sequence->fLoopConfig = new TVectorLoopConfig(increment);
         // sequence->fLoopConfig = new TAssocLoopConfig(proxy); */
      } else {
         sequence->fLoopConfig = new TGenericLoopConfig(&proxy, /* read */ kFALSE);
      }
      for (UInt_t i = 0; i < ndata; ++i) {
         TStreamerElement *element = (TStreamerElement*) info->GetElements()->At(i);
         if (!element) {
            break;
         }
         if (element->GetType() < 0) {
            // -- Skip an ignored TObject base class.
            // Note: The only allowed negative value here is -1, and signifies that Build() has found a TObject
            // base class and TClass::IgnoreTObjectStreamer() was called.  In this case the compiled version of the
            // elements omits the TObject base class element, which has to be compensated for by TTree::Bronch()
            // when it is making branches for a split object.
            continue;
         }
         if (element->TestBit(TStreamerElement::kCache) && !element->TestBit(TStreamerElement::kWrite)) {
            // Skip element cached for reading purposes.
            continue;
         }
         if (element->GetType() >= TVirtualStreamerInfo::kArtificial &&  !element->TestBit(TStreamerElement::kWrite)) {
            // Skip artificial element used for reading purposes.
            continue;
         }
         TStreamerInfo::TCompInfo *compinfo = sinfo->fCompFull[i];
         Int_t asize = element->GetSize();
         if (element->GetArrayLength()) {
            asize /= element->GetArrayLength();
         }
         Int_t oldType = element->GetType();
         Int_t offset = element->GetOffset();
#if defined(CDJ_NO_COMPILE)
         Int_t newType = element->GetNewType();

         if (newType != oldType) {
            if (newType > 0) {
               if (oldType != TVirtualStreamerInfo::kCounter) {
                  oldType += TVirtualStreamerInfo::kConv;
               }
            } else {
               oldType += TVirtualStreamerInfo::kSkip;
            }
         }
         if ( IsDefaultVector(proxy)
               /*|| (proxy.GetCollectionType() == ROOT::kSTLset || proxy.GetCollectionType() == ROOT::kSTLmultiset
               || proxy.GetCollectionType() == ROOT::kSTLmap || proxy.GetCollectionType() == ROOT::kSTLmultimap) */ )
         {

            // We can speed up the iteration in case of vector.  We also know that all emulated collection are stored internally as a vector.
            if (element->TestBit(TStreamerElement::kCache)) {
               TConfiguredAction action( GetCollectionWriteAction<VectorLooper>(info,element,oldType,i,compinfo,offset) );
               sequence->AddAction( UseCacheVectorLoop,  new TConfigurationUseCache(info,action,element->TestBit(TStreamerElement::kRepeat)) );
            } else {
               sequence->AddAction(GetCollectionWriteAction<VectorLooper>(info,element,oldType,i,compinfo,offset));
            }

   //         } else if (proxy.GetCollectionType() == ROOT::kSTLset || proxy.GetCollectionType() == ROOT::kSTLmultiset
   //                    || proxy.GetCollectionType() == ROOT::kSTLmap || proxy.GetCollectionType() == ROOT::kSTLmultimap) {
   //            sequence->AddAction( GenericAssocCollectionAction, new TConfigSTL(info,i,compinfo,offset,0,proxy.GetCollectionClass(),0,0) );
         } else {
            // The usual collection case.
            if (element->TestBit(TStreamerElement::kCache)) {
               TConfiguredAction action( GetWriteAction<VectorLooper>(info,element,oldType,i,compinfo,offset) );
               sequence->AddAction( UseCacheGenericCollection, new TConfigurationUseCache(info,action,element->TestBit(TStreamerElement::kRepeat)) );
            } else {
               switch (oldType) {
                     // read basic types
                  case TVirtualStreamerInfo::kBool:    sequence->AddAction( WriteBasicTypeGenericLoop<Bool_t>, new TConfiguration(info,i,compinfo,offset) );    break;
                  case TVirtualStreamerInfo::kChar:    sequence->AddAction( WriteBasicTypeGenericLoop<Char_t>, new TConfiguration(info,i,compinfo,offset) );    break;
                  case TVirtualStreamerInfo::kShort:   sequence->AddAction( WriteBasicTypeGenericLoop<Short_t>, new TConfiguration(info,i,compinfo,offset) );   break;
                  case TVirtualStreamerInfo::kInt:     sequence->AddAction( WriteBasicTypeGenericLoop<Int_t>, new TConfiguration(info,i,compinfo,offset) );     break;
                  case TVirtualStreamerInfo::kLong:    sequence->AddAction( WriteBasicTypeGenericLoop<Long_t>, new TConfiguration(info,i,compinfo,offset) );    break;
                  case TVirtualStreamerInfo::kLong64:  sequence->AddAction( WriteBasicTypeGenericLoop<Long64_t>, new TConfiguration(info,i,compinfo,offset) );  break;
                  case TVirtualStreamerInfo::kFloat:   sequence->AddAction( WriteBasicTypeGenericLoop<Float_t>, new TConfiguration(info,i,compinfo,offset) );   break;
                  case TVirtualStreamerInfo::kDouble:  sequence->AddAction( WriteBasicTypeGenericLoop<Double_t>, new TConfiguration(info,i,compinfo,offset) );  break;
                  case TVirtualStreamerInfo::kUChar:   sequence->AddAction( WriteBasicTypeGenericLoop<UChar_t>, new TConfiguration(info,i,compinfo,offset) );   break;
                  case TVirtualStreamerInfo::kUShort:  sequence->AddAction( WriteBasicTypeGenericLoop<UShort_t>, new TConfiguration(info,i,compinfo,offset) );  break;
                  case TVirtualStreamerInfo::kUInt:    sequence->AddAction( WriteBasicTypeGenericLoop<UInt_t>, new TConfiguration(info,i,compinfo,offset) );    break;
                  case TVirtualStreamerInfo::kULong:   sequence->AddAction( WriteBasicTypeGenericLoop<ULong_t>, new TConfiguration(info,i,compinfo,offset) );   break;
                  case TVirtualStreamerInfo::kULong64: sequence->AddAction( WriteBasicTypeGenericLoop<ULong64_t>, new TConfiguration(info,i,compinfo,offset) ); break;
                  // case TVirtualStreamerInfo::kBits:    sequence->AddAction( WriteBasicTypeGenericLoop<BitsMarker>, new TConfiguration(info,i,compinfo,offset) );    break;
                  case TVirtualStreamerInfo::kFloat16: {
                     if (element->GetFactor() != 0) {
                        sequence->AddAction( GenericLooper<WriteBasicType_WithFactor<float> >, new TConfWithFactor(info,i,compinfo,offset,element->GetFactor(),element->GetXmin()) );
                     } else {
                        Int_t nbits = (Int_t)element->GetXmin();
                        if (!nbits) nbits = 12;
                        sequence->AddAction( GenericLooper<WriteBasicType_NoFactor<float> >, new TConfNoFactor(info,i,compinfo,offset,nbits) );
                     }
                     break;
                  }
                  case TVirtualStreamerInfo::kDouble32: {
                     if (element->GetFactor() != 0) {
                        sequence->AddAction( GenericLooper<WriteBasicType_WithFactor<double> >, new TConfWithFactor(info,i,compinfo,offset,element->GetFactor(),element->GetXmin()) );
                     } else {
                        Int_t nbits = (Int_t)element->GetXmin();
                        if (!nbits) {
                           sequence->AddAction( GenericLooper<ConvertBasicType<float,double> >, new TConfiguration(info,i,compinfo,offset) );
                        } else {
                           sequence->AddAction( GenericLooper<WriteBasicType_NoFactor<double> >, new TConfNoFactor(info,i,compinfo,offset,nbits) );
                        }
                     }
                     break;
                  }
                  case TVirtualStreamerInfo::kTNamed:  sequence->AddAction( GenericLooper<WriteTNamed >, new TConfiguration(info,i,compinfo,offset) );    break;
                     // Idea: We should calculate the CanIgnoreTObjectStreamer here and avoid calling the
                     // Streamer alltogether.
                  case TVirtualStreamerInfo::kTObject: sequence->AddAction( GenericLooper<WriteTObject >, new TConfiguration(info,i,compinfo,offset) );    break;
                  case TVirtualStreamerInfo::kTString: sequence->AddAction( GenericLooper<WriteTString >, new TConfiguration(info,i,compinfo,offset) );    break;
                  default:
                     sequence->AddAction( GenericCollectionWriteAction, new TConfigSTL(info,i,0 /* the offset will be used from TStreamerInfo */,0,proxy.GetCollectionClass(),0,0) );
                     break;
               }
            }
         }
#else
         if ( IsDefaultVector(proxy)
               /*|| (proxy.GetCollectionType() == ROOT::kSTLset || proxy.GetCollectionType() == ROOT::kSTLmultiset
                || proxy.GetCollectionType() == ROOT::kSTLmap || proxy.GetCollectionType() == ROOT::kSTLmultimap)*/ )
         {
            sequence->AddAction( GetCollectionWriteAction<VectorLooper>(info,element,oldType,i,compinfo,offset) );
         } else {
            // NOTE: TBranch::FillLeavesCollection[Member] is not yet ready to handle the sequence
            // as it does not create/use a TStaging as expected ... but then again it might
            // not be the right things to expect ...
            // sequence->AddAction( GetCollectionWriteAction<GenericLooper>(info,element,oldType,i,compinfo,offset) );
            sequence->AddAction( GenericLooper::GenericWrite, new TConfigSTL(info,i,compinfo,0 /* the offset will be used from TStreamerInfo */,0,proxy.GetCollectionClass(),0,0) );
         }
#endif
      }
      return sequence;
}

void TStreamerInfoActions::TActionSequence::AddToOffset(Int_t delta)
{
   // Add the (potentially negative) delta to all the configuration's offset.  This is used by
   // TBranchElement in the case of split sub-object.

   TStreamerInfoActions::ActionContainer_t::iterator end = fActions.end();
   for(TStreamerInfoActions::ActionContainer_t::iterator iter = fActions.begin();
       iter != end;
       ++iter)
   {
      if (!iter->fConfiguration->fInfo->GetElements()->At(iter->fConfiguration->fElemId)->TestBit(TStreamerElement::kCache))
         iter->fConfiguration->AddToOffset(delta);
   }
}

void TStreamerInfoActions::TActionSequence::SetMissing()
{
   // Add the (potentially negative) delta to all the configuration's offset.  This is used by
   // TBranchElement in the case of split sub-object.

   TStreamerInfoActions::ActionContainer_t::iterator end = fActions.end();
   for(TStreamerInfoActions::ActionContainer_t::iterator iter = fActions.begin();
       iter != end;
       ++iter)
   {
      if (!iter->fConfiguration->fInfo->GetElements()->At(iter->fConfiguration->fElemId)->TestBit(TStreamerElement::kCache))
         iter->fConfiguration->SetMissing();
   }
}

TStreamerInfoActions::TActionSequence *TStreamerInfoActions::TActionSequence::CreateCopy()
{
   // Create a copy of this sequence.

   TStreamerInfoActions::TActionSequence *sequence = new TStreamerInfoActions::TActionSequence(fStreamerInfo, fActions.size(), IsForVectorPtrLooper());

   sequence->fLoopConfig = fLoopConfig ? fLoopConfig->Copy() : 0;

   TStreamerInfoActions::ActionContainer_t::iterator end = fActions.end();
   for(TStreamerInfoActions::ActionContainer_t::iterator iter = fActions.begin();
       iter != end;
       ++iter)
   {
      TConfiguration *conf = iter->fConfiguration->Copy();
      sequence->AddAction( iter->fAction, conf );
   }
   return sequence;
}

void TStreamerInfoActions::TActionSequence::AddToSubSequence(TStreamerInfoActions::TActionSequence *sequence,
      const TStreamerInfoActions::TIDs &element_ids,
      Int_t offset,
      TStreamerInfoActions::TActionSequence::SequenceGetter_t create)
{
   for(UInt_t id = 0; id < element_ids.size(); ++id) {
      if ( element_ids[id].fElemID < 0 ) {
         if (element_ids[id].fNestedIDs) {
            auto original = create(element_ids[id].fNestedIDs->fInfo,
                                   sequence->fLoopConfig ? sequence->fLoopConfig->GetCollectionProxy() : nullptr,
                                   nullptr);
            if (element_ids[id].fNestedIDs->fOnfileObject) {
               auto conf = new TConfigurationPushDataCache(element_ids[id].fNestedIDs->fInfo, element_ids[id].fNestedIDs->fOnfileObject, offset);
               if ( sequence->fLoopConfig )
                  sequence->AddAction( PushDataCacheGenericCollection, conf );
               else if ( sequence->IsForVectorPtrLooper() )
                  sequence->AddAction( PushDataCacheVectorPtr, conf );
               else
                  sequence->AddAction( PushDataCache, conf );
            }

            original->AddToSubSequence(sequence, element_ids[id].fNestedIDs->fIDs, element_ids[id].fNestedIDs->fOffset, create);

            if (element_ids[id].fNestedIDs->fOnfileObject) {
               auto conf =
                  new TConfigurationPushDataCache(element_ids[id].fNestedIDs->fInfo, nullptr, element_ids[id].fNestedIDs->fOffset);
               if ( sequence->fLoopConfig )
                  sequence->AddAction( PopDataCacheGenericCollection, conf );
               else if ( sequence->IsForVectorPtrLooper() )
                  sequence->AddAction( PopDataCacheVectorPtr, conf );
               else
                  sequence->AddAction( PopDataCache, conf );
            }
         } else {
            TStreamerInfoActions::ActionContainer_t::iterator end = fActions.end();
            for(TStreamerInfoActions::ActionContainer_t::iterator iter = fActions.begin();
               iter != end;
               ++iter)
            {
               TConfiguration *conf = iter->fConfiguration->Copy();
               if (!iter->fConfiguration->fInfo->GetElements()->At(iter->fConfiguration->fElemId)->TestBit(TStreamerElement::kCache))
                  conf->AddToOffset(offset);
               sequence->AddAction( iter->fAction, conf );
            }
         }
      } else {
         int localIndex = 0;
         TStreamerInfoActions::ActionContainer_t::iterator end = fActions.end();
         for(TStreamerInfoActions::ActionContainer_t::iterator iter = fActions.begin();
             iter != end;
             ++iter) {
            // fprintf(stderr, "With element_ids[%d] For %s comparing act[%d/%zu] %d to %d  for %p vs %p %s\n",
            //         id,
            //         iter->fConfiguration->fInfo->GetName(),
            //         localIndex, fActions.size(),
            //         iter->fConfiguration->fElemId,
            //         (UInt_t)element_ids[id].fElemID, iter->fConfiguration->fInfo,
            //         element_ids[id].fInfo,
            //         element_ids[id].fInfo ? element_ids[id].fInfo->GetName() : "nullptr" );
            ++localIndex;
            if ( iter->fConfiguration->fElemId == (UInt_t)element_ids[id].fElemID ) {
               TConfiguration *conf = iter->fConfiguration->Copy();
               if (!iter->fConfiguration->fInfo->GetElements()->At(iter->fConfiguration->fElemId)->TestBit(TStreamerElement::kCache))
                  conf->AddToOffset(offset);
               sequence->AddAction( iter->fAction, conf );
            }
         }
      }
   }
}

TStreamerInfoActions::TActionSequence *TStreamerInfoActions::TActionSequence::CreateSubSequence(const TIDs &element_ids, size_t offset,
      TStreamerInfoActions::TActionSequence::SequenceGetter_t create)
{
   // Create a sequence containing the subset of the action corresponding to the SteamerElement whose ids is contained in the vector.
   // 'offset' is the location of this 'class' within the object (address) that will be passed to ReadBuffer when using this sequence.

   TStreamerInfoActions::TActionSequence *sequence = new TStreamerInfoActions::TActionSequence(fStreamerInfo, element_ids.size(), IsForVectorPtrLooper());

   sequence->fLoopConfig = fLoopConfig ? fLoopConfig->Copy() : 0;

   AddToSubSequence(sequence, element_ids, offset, create);

   return sequence;
}

TStreamerInfoActions::TActionSequence *TStreamerInfoActions::TActionSequence::CreateSubSequence(const std::vector<Int_t> &element_ids, size_t offset)
{
   // Create a sequence containing the subset of the action corresponding to the SteamerElement whose ids is contained in the vector.
   // 'offset' is the location of this 'class' within the object (address) that will be passed to ReadBuffer when using this sequence.

   TStreamerInfoActions::TActionSequence *sequence = new TStreamerInfoActions::TActionSequence(fStreamerInfo, element_ids.size(), IsForVectorPtrLooper());

   sequence->fLoopConfig = fLoopConfig ? fLoopConfig->Copy() : 0;

   for(UInt_t id = 0; id < element_ids.size(); ++id) {
      if ( element_ids[id] < 0 ) {
         TStreamerInfoActions::ActionContainer_t::iterator end = fActions.end();
         for(TStreamerInfoActions::ActionContainer_t::iterator iter = fActions.begin();
             iter != end;
             ++iter)
         {
            TConfiguration *conf = iter->fConfiguration->Copy();
            if (!iter->fConfiguration->fInfo->GetElements()->At(iter->fConfiguration->fElemId)->TestBit(TStreamerElement::kCache))
               conf->AddToOffset(offset);
            sequence->AddAction( iter->fAction, conf );
         }
      } else {
         TStreamerInfoActions::ActionContainer_t::iterator end = fActions.end();
         for(TStreamerInfoActions::ActionContainer_t::iterator iter = fActions.begin();
             iter != end;
             ++iter) {
            if ( iter->fConfiguration->fElemId == (UInt_t)element_ids[id] ) {
               TConfiguration *conf = iter->fConfiguration->Copy();
               if (!iter->fConfiguration->fInfo->GetElements()->At(iter->fConfiguration->fElemId)->TestBit(TStreamerElement::kCache))
                  conf->AddToOffset(offset);
               sequence->AddAction( iter->fAction, conf );
            }
         }
      }
   }
   return sequence;
}

#if !defined(R__WIN32) && !defined(_AIX)

#include <dlfcn.h>

#endif

typedef void (*voidfunc)();
static const char *R__GetSymbolName(voidfunc func)
{
#if defined(R__WIN32) || defined(__CYGWIN__) || defined(_AIX)
   return "not available on this platform";
#if 0
   MEMORY_BASIC_INFORMATION mbi;
   if (!VirtualQuery (func, &mbi, sizeof (mbi)))
   {
      return 0;
   }

   HMODULE hMod = (HMODULE) mbi.AllocationBase;
   static char moduleName[MAX_PATH];

   if (!GetModuleFileNameA (hMod, moduleName, sizeof (moduleName)))
   {
      return 0;
   }
   return moduleName;
#endif
#else
   Dl_info info;
   if (dladdr((void*)func,&info)==0) {
      // Not in a known share library, let's give up
      return "name not found";
   } else {
      //fprintf(stdout,"Found address in %s\n",info.dli_fname);
      return info.dli_sname;
   }
#endif
}

void TStreamerInfoActions::TActionSequence::Print(Option_t *opt) const
{
   // Add the (potentially negative) delta to all the configuration's offset.  This is used by
   // TTBranchElement in the case of split sub-object.
   // If opt contains 'func', also print the (mangled) name of the function that will be executed.

   if (fLoopConfig) {
      fLoopConfig->Print();
   }
   TStreamerInfoActions::ActionContainer_t::const_iterator end = fActions.end();
   for(TStreamerInfoActions::ActionContainer_t::const_iterator iter = fActions.begin();
       iter != end;
       ++iter)
   {
      iter->fConfiguration->Print();
      if (strstr(opt,"func")) {
         printf("StreamerInfoAction func: %s\n",R__GetSymbolName((voidfunc)iter->fAction));
      }
   }
}


