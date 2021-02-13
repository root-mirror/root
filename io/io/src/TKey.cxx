// @(#)root/io:$Id$
// Author: Rene Brun   28/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
\class TKey
\ingroup IO

 Book space in a file, create I/O buffers, to fill them, (un)compress them.

 The TKey class includes functions to book space in a file, to create I/O
 buffers, to fill these buffers, to compress/uncompress data buffers.
 Before saving (making persistent) an object in a file, a key must
 be created. The key structure contains all the information to
 uniquely identify a persistent object in a file.
 | Data Member | Explanation |
 |-------------|-------------|
 |  fNbytes    | Number of bytes for the compressed object and key. |
 |  fObjlen    | Length of uncompressed object. |
 |  fDatime    | Date/Time when the object was written. |
 |  fKeylen    | Number of bytes for the key structure. |
 |  fCycle     | Cycle number of the object. |
 |  fSeekKey   | Address of the object on file (points to fNbytes). This is a redundant information used to cross-check the data base integrity. |
 |  fSeekPdir  | Pointer to the directory supporting this object.|
 |  fClassName | Object class name. |
 |  fName      | Name of the object. |
 |  fTitle     | Title of the object. |

 In the 16 highest bits of fSeekPdir is encoded a pid offset.  This
 offset is to be added to the pid index stored in the TRef object
 and the referenced TObject.

 The TKey class is used by ROOT to:
   - Write an object in the current directory
   - Write a new ntuple buffer

 The structure of a file is shown in TFile::TFile.
 The structure of a directory is shown in TDirectoryFile::TDirectoryFile.
 The TKey class is used by the TBasket class.
 See also TTree.
*/

#include <atomic>
#include <iostream>

#include "TROOT.h"
#include "TClass.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TKey.h"
#include "TBufferFile.h"
#include "TFree.h"
#include "TBrowser.h"
#include "Bytes.h"
#include "TInterpreter.h"
#include "TError.h"
#include "TVirtualStreamerInfo.h"
#include "TSchemaRuleSet.h"
#include "ThreadLocalStorage.h"

#include "RZip.h"

const Int_t kTitleMax = 32000;
#if 0
const Int_t kMAXFILEBUFFER = 262144;
#endif

#if !defined(_MSC_VER) || (_MSC_VER>1300)
const ULong64_t kPidOffsetMask = 0xffffffffffffULL;
#else
const ULong64_t kPidOffsetMask = 0xffffffffffffUL;
#endif
const UChar_t kPidOffsetShift = 48;

const static TString gTDirectoryString("TDirectory");
std::atomic<UInt_t> keyAbsNumber{0};

ClassImp(TKey);

////////////////////////////////////////////////////////////////////////////////
/// TKey default constructor.

TKey::TKey() : TNamed(), fDatime((UInt_t)0)
{
   Build(0, "", 0);

   fKeylen     = Sizeof();

   keyAbsNumber++; SetUniqueID(keyAbsNumber);
}

////////////////////////////////////////////////////////////////////////////////
/// TKey default constructor.

TKey::TKey(TDirectory* motherDir) : TNamed(), fDatime((UInt_t)0)
{
   Build(motherDir, "", 0);

   fKeylen     = Sizeof();

   keyAbsNumber++; SetUniqueID(keyAbsNumber);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a TKey from its original directory to the new 'motherDir'

TKey::TKey(TDirectory* motherDir, const TKey &orig, UShort_t pidOffset) : TNamed(), fDatime((UInt_t)0)
{
   fMotherDir  = motherDir;

   fPidOffset  = orig.fPidOffset + pidOffset;
   fNbytes     = orig.fNbytes;
   fObjlen     = orig.fObjlen;
   fClassName  = orig.fClassName;
   fName       = orig.fName;
   fTitle      = orig.fTitle;

   fCycle      = fMotherDir->AppendKey(this);
   fSeekPdir   = 0;
   fSeekKey    = 0;
   fLeft       = 0;

   fVersion    = TKey::Class_Version();
   Long64_t filepos = GetFile()->GetEND();
   if (filepos > TFile::kStartBigFile || fPidOffset) fVersion += 1000;

   fKeylen     = Sizeof();  // fVersion must be set.

   UInt_t bufferDecOffset = 0;
   UInt_t bufferIncOffset = 0;
   UInt_t alloc = fNbytes + sizeof(Int_t);  // The extra Int_t is for any free space information.
   if (fKeylen < orig.fKeylen) {
      bufferDecOffset = orig.fKeylen - fKeylen;
      fNbytes -= bufferDecOffset;
   } else if (fKeylen > orig.fKeylen) {
      bufferIncOffset = fKeylen - orig.fKeylen;
      alloc += bufferIncOffset;
      fNbytes += bufferIncOffset;
   }

   fBufferRef  = new TBufferFile(TBuffer::kWrite, alloc);
   fBuffer     = fBufferRef->Buffer();

   // Steal the data from the old key.

   TFile* f = orig.GetFile();
   if (f) {
      Int_t nsize = orig.fNbytes;
      f->Seek(orig.fSeekKey);
      if( f->ReadBuffer(fBuffer+bufferIncOffset,nsize) )
      {
         Error("ReadFile", "Failed to read data.");
         return;
      }
      if (gDebug) {
         std::cout << "TKey Reading "<<nsize<< " bytes at address "<<fSeekKey<<std::endl;
      }
   }
   fBuffer += bufferDecOffset; // Reset the buffer to be appropriate for this key.
   Int_t nout = fNbytes - fKeylen;
   Create(nout);
   fBufferRef->SetBufferOffset(bufferDecOffset);
   Streamer(*fBufferRef);         //write key itself again
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TKey object to read keys.
/// Constructor called by TDirectoryFile::ReadKeys and by TFile::TFile.
/// A TKey object is created to read the keys structure itself.

TKey::TKey(Long64_t pointer, Int_t nbytes, TDirectory* motherDir) : TNamed()
{
   Build(motherDir, "", pointer);

   fSeekKey    = pointer;
   fNbytes     = nbytes;
   fBuffer     = new char[nbytes];
   keyAbsNumber++; SetUniqueID(keyAbsNumber);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TKey object with the specified name, title for the given class.
///
///  WARNING: in name avoid special characters like '^','$','.' that are used
///  by the regular expression parser (see TRegexp).

TKey::TKey(const char *name, const char *title, const TClass *cl, Int_t nbytes, TDirectory* motherDir)
      : TNamed(name,title)
{
   Build(motherDir, cl->GetName(), -1);

   fKeylen     = Sizeof();
   fObjlen     = nbytes;
   Create(nbytes);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TKey object with the specified name, title for the given class.
///
///  WARNING: in name avoid special characters like '^','$','.' that are used
///  by the regular expression parser (see TRegexp).

TKey::TKey(const TString &name, const TString &title, const TClass *cl, Int_t nbytes, TDirectory* motherDir)
      : TNamed(name,title)
{
   Build(motherDir, cl->GetName(), -1);

   fKeylen     = Sizeof();
   fObjlen     = nbytes;
   Create(nbytes);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TKey object for a TObject* and fill output buffer
///
///  WARNING: in name avoid special characters like '^','$','.' that are used
///  by the regular expression parser (see TRegexp).

TKey::TKey(const TObject *obj, const char *name, Int_t bufsize, TDirectory* motherDir)
     : TNamed(name, obj->GetTitle())
{
   R__ASSERT(obj);

   if (!obj->IsA()->HasDefaultConstructor()) {
      Warning("TKey", "since %s has no public constructor\n"
              "\twhich can be called without argument, objects of this class\n"
              "\tcan not be read with the current library. You will need to\n"
              "\tadd a default constructor before attempting to read it.",
              obj->ClassName());
   }

   Build(motherDir, obj->ClassName(), -1);

   Int_t lbuf, nout, noutot, bufmax, nzip;
   fBufferRef = new TBufferFile(TBuffer::kWrite, bufsize);
   fBufferRef->SetParent(GetFile());
   fCycle     = fMotherDir->AppendKey(this);

   Streamer(*fBufferRef);         //write key itself
   fKeylen    = fBufferRef->Length();
   fBufferRef->MapObject(obj);    //register obj in map in case of self reference
   ((TObject*)obj)->Streamer(*fBufferRef);    //write object
   lbuf       = fBufferRef->Length();
   fObjlen    = lbuf - fKeylen;

   Int_t cxlevel = GetFile() ? GetFile()->GetCompressionLevel() : 0;
   ROOT::RCompressionSetting::EAlgorithm::EValues cxAlgorithm = static_cast<ROOT::RCompressionSetting::EAlgorithm::EValues>(GetFile() ? GetFile()->GetCompressionAlgorithm() : 0);
   if (cxlevel > 0 && fObjlen > 256) {
      Int_t nbuffers = 1 + (fObjlen - 1)/kMAXZIPBUF;
      Int_t buflen = TMath::Max(512,fKeylen + fObjlen + 9*nbuffers + 28); //add 28 bytes in case object is placed in a deleted gap
      fBuffer = new char[buflen];
      char *objbuf = fBufferRef->Buffer() + fKeylen;
      char *bufcur = &fBuffer[fKeylen];
      noutot = 0;
      nzip   = 0;
      for (Int_t i = 0; i < nbuffers; ++i) {
         if (i == nbuffers - 1) bufmax = fObjlen - nzip;
         else               bufmax = kMAXZIPBUF;
         R__zipMultipleAlgorithm(cxlevel, &bufmax, objbuf, &bufmax, bufcur, &nout, cxAlgorithm);
         if (nout == 0 || nout >= fObjlen) { //this happens when the buffer cannot be compressed
            fBuffer = fBufferRef->Buffer();
            Create(fObjlen);
            fBufferRef->SetBufferOffset(0);
            Streamer(*fBufferRef);         //write key itself again
            return;
         }
         bufcur += nout;
         noutot += nout;
         objbuf += kMAXZIPBUF;
         nzip   += kMAXZIPBUF;
      }
      Create(noutot);
      fBufferRef->SetBufferOffset(0);
      Streamer(*fBufferRef);         //write key itself again
      memcpy(fBuffer,fBufferRef->Buffer(),fKeylen);
      delete fBufferRef; fBufferRef = 0;
   } else {
      fBuffer = fBufferRef->Buffer();
      Create(fObjlen);
      fBufferRef->SetBufferOffset(0);
      Streamer(*fBufferRef);         //write key itself again
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TKey object for any object obj of class cl d and fill
/// output buffer.
///
///  WARNING: in name avoid special characters like '^','$','.' that are used
///  by the regular expression parser (see TRegexp).

TKey::TKey(const void *obj, const TClass *cl, const char *name, Int_t bufsize, TDirectory* motherDir)
     : TNamed(name, "object title")
{
   R__ASSERT(obj && cl);

   if (!cl->HasDefaultConstructor()) {
      Warning("TKey", "since %s has no public constructor\n"
              "\twhich can be called without argument, objects of this class\n"
              "\tcan not be read with the current library. You will need to\n"
              "\tadd a default constructor before attempting to read it.",
              cl->GetName());
   }

   TClass *clActual = cl->GetActualClass(obj);
   const void* actualStart;
   if (clActual) {
      const char *temp = (const char*) obj;
      // clActual->GetStreamerInfo();
      Int_t offset = (cl != clActual) ?
                     clActual->GetBaseClassOffset(cl) : 0;
      temp -= offset;
      actualStart = temp;
   } else {
      // We could not determine the real type of this object,
      // let's assume it is the one given by the caller.
      clActual = const_cast<TClass*>(cl);
      actualStart = obj;
   }

   Build(motherDir, clActual->GetName(), -1);

   fBufferRef = new TBufferFile(TBuffer::kWrite, bufsize);
   fBufferRef->SetParent(GetFile());
   fCycle     = fMotherDir->AppendKey(this);

   Streamer(*fBufferRef);         //write key itself
   fKeylen    = fBufferRef->Length();

   Int_t lbuf, nout, noutot, bufmax, nzip;

   fBufferRef->MapObject(actualStart,clActual);         //register obj in map in case of self reference
   clActual->Streamer((void*)actualStart, *fBufferRef); //write object
   lbuf       = fBufferRef->Length();
   fObjlen    = lbuf - fKeylen;

   Int_t cxlevel = GetFile() ? GetFile()->GetCompressionLevel() : 0;
   ROOT::RCompressionSetting::EAlgorithm::EValues cxAlgorithm = static_cast<ROOT::RCompressionSetting::EAlgorithm::EValues>(GetFile() ? GetFile()->GetCompressionAlgorithm() : 0);
   if (cxlevel > 0 && fObjlen > 256) {
      Int_t nbuffers = 1 + (fObjlen - 1)/kMAXZIPBUF;
      Int_t buflen = TMath::Max(512,fKeylen + fObjlen + 9*nbuffers + 28); //add 28 bytes in case object is placed in a deleted gap
      fBuffer = new char[buflen];
      char *objbuf = fBufferRef->Buffer() + fKeylen;
      char *bufcur = &fBuffer[fKeylen];
      noutot = 0;
      nzip   = 0;
      for (Int_t i = 0; i < nbuffers; ++i) {
         if (i == nbuffers - 1) bufmax = fObjlen - nzip;
         else               bufmax = kMAXZIPBUF;
         R__zipMultipleAlgorithm(cxlevel, &bufmax, objbuf, &bufmax, bufcur, &nout, cxAlgorithm);
         if (nout == 0 || nout >= fObjlen) { //this happens when the buffer cannot be compressed
            fBuffer = fBufferRef->Buffer();
            Create(fObjlen);
            fBufferRef->SetBufferOffset(0);
            Streamer(*fBufferRef);         //write key itself again
            return;
         }
         bufcur += nout;
         noutot += nout;
         objbuf += kMAXZIPBUF;
         nzip   += kMAXZIPBUF;
      }
      Create(noutot);
      fBufferRef->SetBufferOffset(0);
      Streamer(*fBufferRef);         //write key itself again
      memcpy(fBuffer,fBufferRef->Buffer(),fKeylen);
      delete fBufferRef; fBufferRef = 0;
   } else {
      fBuffer = fBufferRef->Buffer();
      Create(fObjlen);
      fBufferRef->SetBufferOffset(0);
      Streamer(*fBufferRef);         //write key itself again
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Method used in all TKey constructor to initialize basic data fields.
///
/// The member filepos is used to calculate correct version number of key
/// if filepos==-1, end of file position is used.

void TKey::Build(TDirectory* motherDir, const char* classname, Long64_t filepos)
{
   fMotherDir = motherDir;

   fPidOffset  = 0;
   fNbytes     = 0;
   fBuffer     = 0;
   fKeylen     = 0;
   fObjlen     = 0;
   fBufferRef  = 0;
   fCycle      = 0;
   fSeekPdir   = 0;
   fSeekKey    = 0;
   fLeft       = 0;

   fClassName = classname;
   //the following test required for forward and backward compatibility
   if (fClassName == "TDirectoryFile") SetBit(kIsDirectoryFile);

   fVersion = TKey::Class_Version();

   if ((filepos==-1) && GetFile()) filepos = GetFile()->GetEND();
   if (filepos > TFile::kStartBigFile) fVersion += 1000;

   if (fTitle.Length() > kTitleMax) fTitle.Resize(kTitleMax);

   if (GetFile() && GetFile()->TestBit(TFile::kReproducible))
      SetBit(TKey::kReproducible);
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from disk and call its Browse() method.
///
/// If object with same name already exist in memory delete it (like
/// TDirectoryFile::Get() is doing), except when the key references a
/// folder in which case we don't want to re-read the folder object
/// since it might contain new objects not yet saved.

void TKey::Browse(TBrowser *b)
{
   if (fMotherDir==0) return;

   TClass *objcl = TClass::GetClass(GetClassName());

   void* obj = fMotherDir->GetList()->FindObject(GetName());
   if (obj && objcl->IsTObject()) {
      TObject *tobj = (TObject*) objcl->DynamicCast(TObject::Class(), obj);
      if (!tobj->IsFolder()) {
         if (tobj->InheritsFrom(TCollection::Class()))
            tobj->Delete();   // delete also collection elements
         delete tobj;
         obj = 0;
      }
   }

   if (!obj)
      obj = ReadObj();

   if (b && obj) {
      objcl->Browse(obj,b);
      b->SetRefreshFlag(kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TKey object of specified size.
///
/// If externFile!=0, key will be allocated in specified file, otherwise file
/// of mother directory will be used.

void TKey::Create(Int_t nbytes, TFile* externFile)
{
   keyAbsNumber++; SetUniqueID(keyAbsNumber);

   TFile *f = externFile;
   if (!f) f = GetFile();
   if (!f) {
      Error("Create","Cannot create key without file");
      return;
   }

   Int_t nsize      = nbytes + fKeylen;
   TList *lfree     = f->GetListOfFree();
   TFree *f1        = (TFree*)lfree->First();
//*-*-------------------find free segment
//*-*                    =================
   TFree *bestfree  = f1->GetBestFree(lfree,nsize);
   if (bestfree == 0) {
      Error("Create","Cannot allocate %d bytes for ID = %s Title = %s",
            nsize,GetName(),GetTitle());
      return;
   }

   if (f->TestBit(TFile::kReproducible))
      SetBit(TKey::kReproducible);

   fDatime.Set();
   fSeekKey  = bestfree->GetFirst();
//*-*----------------- Case Add at the end of the file
   if (fSeekKey >= f->GetEND()) {
      f->SetEND(fSeekKey+nsize);
      bestfree->SetFirst(fSeekKey+nsize);
      if (f->GetEND() > bestfree->GetLast()) {
         bestfree->SetLast(bestfree->GetLast() + 1000000000);
      }
      fLeft   = -1;
      if (!fBuffer) fBuffer = new char[nsize];
   } else {
      fLeft = Int_t(bestfree->GetLast() - fSeekKey - nsize + 1);
   }
//*-*----------------- Case where new object fills exactly a deleted gap
   fNbytes = nsize;
   if (fLeft == 0) {
      if (!fBuffer) {
         fBuffer = new char[nsize];
      }
      lfree->Remove(bestfree);
      delete bestfree;
   }
//*-*----------------- Case where new object is placed in a deleted gap larger than itself
   if (fLeft > 0) {    // found a bigger segment
      if (!fBuffer) {
         fBuffer = new char[nsize+sizeof(Int_t)];
      }
      char *buffer  = fBuffer+nsize;
      Int_t nbytesleft = -fLeft;  // set header of remaining record
      tobuf(buffer, nbytesleft);
      bestfree->SetFirst(fSeekKey+nsize);
   }

   fSeekPdir = externFile ? externFile->GetSeekDir() : fMotherDir->GetSeekDir();
}

////////////////////////////////////////////////////////////////////////////////
/// TKey default destructor.

TKey::~TKey()
{
   if (fMotherDir && fMotherDir->GetListOfKeys())
      fMotherDir->GetListOfKeys()->Remove(this);
   TKey::DeleteBuffer();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete an object from the file.
///
/// Note: the key is not deleted. You still have to call "delete key".
/// This is different from the behaviour of TObject::Delete()!

void TKey::Delete(Option_t *option)
{
   if (TestBit(kIsDirectoryFile)) {
      // TDirectoryFile assumes that its location on file never change (for example updates are partial)
      // and never checks if the space might have been released and thus over-write any data that might
      // have been written there.
      if (option && option[0] == 'v')
         printf("Rejected attempt to delete TDirectoryFile key: %s at address %lld, nbytes = %d\n",
                GetName(), fSeekKey, fNbytes);
      return;
   }
   if (option && option[0] == 'v') printf("Deleting key: %s at address %lld, nbytes = %d\n",GetName(),fSeekKey,fNbytes);
   Long64_t first = fSeekKey;
   Long64_t last  = fSeekKey + fNbytes -1;
   if (GetFile()) GetFile()->MakeFree(first, last);  // release space used by this key
   fMotherDir->GetListOfKeys()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete key buffer(s).

void TKey::DeleteBuffer()
{
   if (fBufferRef) {
      delete fBufferRef;
      fBufferRef = 0;
   } else {
      // We only need to delete fBuffer if fBufferRef is zero because
      // if fBufferRef exists, we delegate ownership of fBuffer to fBufferRef.
      if (fBuffer) {
         delete [] fBuffer;
      }
   }
   fBuffer = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return cycle number associated to this key.

Short_t TKey::GetCycle() const
{
   return ((fCycle >0) ? fCycle : -fCycle);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns file to which key belong.

TFile *TKey::GetFile() const
{
   return fMotherDir!=0 ? fMotherDir->GetFile() : gFile;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the "KEEP" status.

Short_t TKey::GetKeep() const
{
   return ((fCycle >0) ? 0 : 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Encode key header into output buffer.

void TKey::FillBuffer(char *&buffer)
{
   tobuf(buffer, fNbytes);
   Version_t version = fVersion;
   tobuf(buffer, version);

   tobuf(buffer, fObjlen);
   if (TestBit(TKey::kReproducible))
      TDatime((UInt_t) 1).FillBuffer(buffer);
   else
      fDatime.FillBuffer(buffer);
   tobuf(buffer, fKeylen);
   tobuf(buffer, fCycle);
   if (fVersion > 1000) {
      tobuf(buffer, fSeekKey);

      // We currently store in the 16 highest bit of fSeekPdir the value of
      // fPidOffset.  This offset is used when a key (or basket) is transfered from one
      // file to the other.  In this case the TRef and TObject might have stored a
      // pid index (to retrieve TProcessIDs) which refered to their order on the original
      // file, the fPidOffset is to be added to those values to correctly find the
      // TProcessID.  This fPidOffset needs to be increment if the key/basket is copied
      // and need to be zero for new key/basket.
      Long64_t pdir = (((Long64_t)fPidOffset)<<kPidOffsetShift) | (kPidOffsetMask & fSeekPdir);
      tobuf(buffer, pdir);
   } else {
      tobuf(buffer, (Int_t)fSeekKey);
      tobuf(buffer, (Int_t)fSeekPdir);
   }
   if (TestBit(kIsDirectoryFile)) {
      // We want to record "TDirectory" instead of TDirectoryFile so that the file can be read by ancient version of ROOT.
      gTDirectoryString.FillBuffer(buffer);
   } else {
      fClassName.FillBuffer(buffer);
   }
   fName.FillBuffer(buffer);
   fTitle.FillBuffer(buffer);
}

////////////////////////////////////////////////////////////////////////////////
/// Increment fPidOffset by 'offset'.
///
/// This offset is used when a key (or basket) is transfered from one file to
/// the other.  In this case the TRef and TObject might have stored a pid
/// index (to retrieve TProcessIDs) which refered to their order on the
/// original file, the fPidOffset is to be added to those values to correctly
/// find the TProcessID.  This fPidOffset needs to be increment if the
/// key/basket is copied and need to be zero for new key/basket.

void TKey::IncrementPidOffset(UShort_t offset)
{
   fPidOffset += offset;
   if (fPidOffset) {
      // We currently store fPidOffset in the 16 highest bit of fSeekPdir, which
      // need to be store as a 64 bit integer.  So we require this key to be
      // a 'large file' key.
      if (fVersion<1000) fVersion += 1000;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if object referenced by the key is a folder.

Bool_t TKey::IsFolder() const
{
   Bool_t ret = kFALSE;

   TClass *classPtr = TClass::GetClass((const char *) fClassName);
   if (classPtr && classPtr->GetState() > TClass::kEmulated && classPtr->IsTObject()) {
      TObject *obj = (TObject *) classPtr->DynamicCast(TObject::Class(), classPtr->New(TClass::kDummyNew));
      if (obj) {
         ret = obj->IsFolder();
         delete obj;
      }
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the "KEEP" status.
///
/// When the KEEP flag is set to 1 the object cannot be purged.

void TKey::Keep()
{
   if (fCycle >0)  fCycle = -fCycle;
}

////////////////////////////////////////////////////////////////////////////////
/// List Key contents.
/// Add indicator of whether it is the current item or a backup copy.

void TKey::ls(Bool_t current) const
{
   TROOT::IndentLevel();
   std::cout <<"KEY: "<<fClassName<<"\t"<<GetName()<<";"<<GetCycle()<<"\t"<<GetTitle();
   std::cout << (current ? " [current cycle]" : " [backup cycle]");
   std::cout <<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// List Key contents.

void TKey::ls(Option_t *) const
{
   TROOT::IndentLevel();
   std::cout <<"KEY: "<<fClassName<<"\t"<<GetName()<<";"<<GetCycle()<<"\t"<<GetTitle()<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Print key contents.

void TKey::Print(Option_t *) const
{
   printf("TKey Name = %s, Title = %s, Cycle = %d\n",GetName(),GetTitle(),GetCycle());
}

////////////////////////////////////////////////////////////////////////////////
/// To read a TObject* from the file.
///
/// The object associated to this key is read from the file into memory
/// Once the key structure is read (via Streamer) the class identifier
/// of the object is known.
/// Using the class identifier we find the TClass object for this class.
/// A TClass object contains a full description (i.e. dictionary) of the
/// associated class. In particular the TClass object can create a new
/// object of the class type it describes. This new object now calls its
/// Streamer function to rebuilt itself.
///
/// Use TKey::ReadObjectAny to read any object non-derived from TObject
///
/// ### Note
/// A C style cast can only be used in the case where the final class
/// of this object derives from TObject as a first inheritance, otherwise
/// one must use a dynamic_cast.
///
/// #### Example1: simplified case
///     class MyClass : public TObject, public AnotherClass
/// then on return, one get away with using:
///     MyClass *obj = (MyClass*)key->ReadObj();
///
/// #### Example2: Usual case (recommended unless performance is critical)
///     MyClass *obj = dynamic_cast<MyClass*>(key->ReadObj());
/// which support also the more complex inheritance like:
///     class MyClass : public AnotherClass, public TObject
///
/// Of course, dynamic_cast<> can also be used in the example 1.

TObject *TKey::ReadObj()
{
   TClass *cl = TClass::GetClass(fClassName.Data());
   if (!cl) {
      Error("ReadObj", "Unknown class %s", fClassName.Data());
      return 0;
   }
   if (!cl->IsTObject()) {
      // in principle user should call TKey::ReadObjectAny!
      return (TObject*)ReadObjectAny(0);
   }

   TBufferFile bufferRef(TBuffer::kRead, fObjlen+fKeylen);
   if (!bufferRef.Buffer()) {
      Error("ReadObj", "Cannot allocate buffer: fObjlen = %d", fObjlen);
      return 0;
   }
   if (GetFile()==0) return 0;
   bufferRef.SetParent(GetFile());
   bufferRef.SetPidOffset(fPidOffset);

   std::unique_ptr<char []> compressedBuffer;
   auto storeBuffer = fBuffer;
   if (fObjlen > fNbytes-fKeylen) {
      compressedBuffer.reset(new char[fNbytes]);
      fBuffer = compressedBuffer.get();
      if( !ReadFile() )                    //Read object structure from file
      {
        fBuffer = 0;
        return 0;
      }
      memcpy(bufferRef.Buffer(),fBuffer,fKeylen);
   } else {
      fBuffer = bufferRef.Buffer();
      if( !ReadFile() ) {                   //Read object structure from file

         fBuffer = 0;
         return 0;
      }
   }
   fBuffer = storeBuffer;

   // get version of key
   bufferRef.SetBufferOffset(sizeof(fNbytes));
   Version_t kvers = bufferRef.ReadVersion();

   bufferRef.SetBufferOffset(fKeylen);
   TObject *tobj = 0;
   // Create an instance of this class

   char *pobj = (char*)cl->New();
   if (!pobj) {
      Error("ReadObj", "Cannot create new object of class %s", fClassName.Data());
      return 0;
   }
   Int_t baseOffset = cl->GetBaseClassOffset(TObject::Class());
   if (baseOffset==-1) {
      // cl does not inherit from TObject.
      // Since this is not possible yet, the only reason we could reach this code
      // is because something is screw up in the ROOT code.
      Fatal("ReadObj","Incorrect detection of the inheritance from TObject for class %s.\n",
            fClassName.Data());
   }
   tobj = (TObject*)(pobj+baseOffset);
   if (kvers > 1)
      bufferRef.MapObject(pobj,cl);  //register obj in map to handle self reference

   if (fObjlen > fNbytes-fKeylen) {
      char *objbuf = bufferRef.Buffer() + fKeylen;
      UChar_t *bufcur = (UChar_t *)&compressedBuffer[fKeylen];
      Int_t nin, nout = 0, nbuf;
      Int_t noutot = 0;
      while (1) {
         Int_t hc = R__unzip_header(&nin, bufcur, &nbuf);
         if (hc!=0) break;
         R__unzip(&nin, bufcur, &nbuf, (unsigned char*) objbuf, &nout);
         if (!nout) break;
         noutot += nout;
         if (noutot >= fObjlen) break;
         bufcur += nin;
         objbuf += nout;
      }
      compressedBuffer.reset(nullptr);
      if (nout) {
         tobj->Streamer(bufferRef); //does not work with example 2 above
      } else {
         // Even-though we have a TObject, if the class is emulated the virtual
         // table may not be 'right', so let's go via the TClass.
         cl->Destructor(pobj);
         return nullptr;
      }
   } else {
      tobj->Streamer(bufferRef);
   }

   if (gROOT->GetForceStyle()) tobj->UseCurrentStyle();

   if (cl->InheritsFrom(TDirectoryFile::Class())) {
      TDirectory *dir = static_cast<TDirectoryFile*>(tobj);
      dir->SetName(GetName());
      dir->SetTitle(GetTitle());
      dir->SetMother(fMotherDir);
      fMotherDir->Append(dir);
   }

   // Append the object to the directory if requested:
   {
      ROOT::DirAutoAdd_t addfunc = cl->GetDirectoryAutoAdd();
      if (addfunc) {
         addfunc(pobj, fMotherDir);
      }
   }

   return tobj;
}

////////////////////////////////////////////////////////////////////////////////
/// To read a TObject* from bufferRead.
///
/// This function is identical to TKey::ReadObj, but it reads directly from
/// bufferRead instead of reading from a file.
/// The object associated to this key is read from the buffer into memory
/// Using the class identifier we find the TClass object for this class.
/// A TClass object contains a full description (i.e. dictionary) of the
/// associated class. In particular the TClass object can create a new
/// object of the class type it describes. This new object now calls its
/// Streamer function to rebuilt itself.
///
/// ### Note
/// This function is called only internally by ROOT classes.
/// Although being public it is not supposed to be used outside ROOT.
/// If used, you must make sure that the bufferRead is large enough to
/// accomodate the object being read.

TObject *TKey::ReadObjWithBuffer(char *bufferRead)
{

   TClass *cl = TClass::GetClass(fClassName.Data());
   if (!cl) {
      Error("ReadObjWithBuffer", "Unknown class %s", fClassName.Data());
      return 0;
   }
   if (!cl->IsTObject()) {
      // in principle user should call TKey::ReadObjectAny!
      return (TObject*)ReadObjectAny(0);
   }

   TBufferFile bufferRef(TBuffer::kRead, fObjlen+fKeylen);
   if (!bufferRef.Buffer()) {
      Error("ReadObjWithBuffer", "Cannot allocate buffer: fObjlen = %d", fObjlen);
      return 0;
   }
   if (GetFile()==0) return 0;
   bufferRef.SetParent(GetFile());
   bufferRef.SetPidOffset(fPidOffset);

   auto storeBuffer = fBuffer;
   if (fObjlen > fNbytes-fKeylen) {
      fBuffer = bufferRead;
      memcpy(bufferRef.Buffer(),fBuffer,fKeylen);
   } else {
      fBuffer = bufferRef.Buffer();
      ReadFile();                    //Read object structure from file
   }
   fBuffer = storeBuffer;

   // get version of key
   bufferRef.SetBufferOffset(sizeof(fNbytes));
   Version_t kvers = bufferRef.ReadVersion();

   bufferRef.SetBufferOffset(fKeylen);
   TObject *tobj = 0;
   // Create an instance of this class

   char *pobj = (char*)cl->New();
   if (!pobj) {
      Error("ReadObjWithBuffer", "Cannot create new object of class %s", fClassName.Data());
      return 0;
   }
   Int_t baseOffset = cl->GetBaseClassOffset(TObject::Class());
   if (baseOffset==-1) {
      // cl does not inherit from TObject.
      // Since this is not possible yet, the only reason we could reach this code
      // is because something is screw up in the ROOT code.
      Fatal("ReadObjWithBuffer","Incorrect detection of the inheritance from TObject for class %s.\n",
            fClassName.Data());
   }
   tobj = (TObject*)(pobj+baseOffset);

   if (kvers > 1)
      bufferRef.MapObject(pobj,cl);  //register obj in map to handle self reference

   if (fObjlen > fNbytes-fKeylen) {
      char *objbuf = bufferRef.Buffer() + fKeylen;
      UChar_t *bufcur = (UChar_t *)&bufferRead[fKeylen];
      Int_t nin, nout = 0, nbuf;
      Int_t noutot = 0;
      while (1) {
         Int_t hc = R__unzip_header(&nin, bufcur, &nbuf);
         if (hc!=0) break;
         R__unzip(&nin, bufcur, &nbuf, (unsigned char*) objbuf, &nout);
         if (!nout) break;
         noutot += nout;
         if (noutot >= fObjlen) break;
         bufcur += nin;
         objbuf += nout;
      }
      if (nout) {
         tobj->Streamer(bufferRef); //does not work with example 2 above
      } else {
         // Even-though we have a TObject, if the class is emulated the virtual
         // table may not be 'right', so let's go via the TClass.
         cl->Destructor(pobj);
         return nullptr;
      }
   } else {
      tobj->Streamer(bufferRef);
   }

   if (gROOT->GetForceStyle()) tobj->UseCurrentStyle();

   if (cl->InheritsFrom(TDirectoryFile::Class())) {
      TDirectory *dir = static_cast<TDirectoryFile*>(tobj);
      dir->SetName(GetName());
      dir->SetTitle(GetTitle());
      dir->SetMother(fMotherDir);
      fMotherDir->Append(dir);
   }

   // Append the object to the directory if requested:
   {
      ROOT::DirAutoAdd_t addfunc = cl->GetDirectoryAutoAdd();
      if (addfunc) {
         addfunc(pobj, fMotherDir);
      }
   }

   return tobj;
}

////////////////////////////////////////////////////////////////////////////////
/// To read an object (non deriving from TObject) from the file.
///
/// If expectedClass is not null, we checked that that actual class of the
/// object stored is suitable to be stored in a pointer pointing to an object
/// of class 'expectedClass'.  We also adjust the value of the returned address
/// so that it is suitable to be cast (C-Style)
/// a pointer pointing to an object of class 'expectedClass'.
///
/// So for example if the class Bottom inherits from Top and the object
/// stored is of type Bottom you can safely do:
/// ~~~{.cpp}
/// auto TopClass = TClass::GetClass("Top");
/// auto ptr = (Top*) key->ReadObjectAny( TopClass );
/// if (ptr==0) printError("the object stored in the key is not of the expected type\n");
/// ~~~
/// The object associated to this key is read from the file into memory.
/// Once the key structure is read (via Streamer) the class identifier
/// of the object is known.
/// Using the class identifier we find the TClass object for this class.
/// A TClass object contains a full description (i.e. dictionary) of the
/// associated class. In particular the TClass object can create a new
/// object of the class type it describes. This new object now calls its
/// Streamer function to rebuilt itself.

void *TKey::ReadObjectAny(const TClass* expectedClass)
{
   TBufferFile bufferRef(TBuffer::kRead, fObjlen+fKeylen);
   if (!bufferRef.Buffer()) {
      Error("ReadObj", "Cannot allocate buffer: fObjlen = %d", fObjlen);
      return 0;
   }
   if (GetFile()==0) return 0;
   bufferRef.SetParent(GetFile());
   bufferRef.SetPidOffset(fPidOffset);

   std::unique_ptr<char []> compressedBuffer;
   auto storeBuffer = fBuffer;
   if (fObjlen > fNbytes-fKeylen) {
      compressedBuffer.reset(new char[fNbytes]);
      fBuffer = compressedBuffer.get();
      ReadFile();                    //Read object structure from file
      memcpy(bufferRef.Buffer(),fBuffer,fKeylen);
   } else {
      fBuffer = bufferRef.Buffer();
      ReadFile();                    //Read object structure from file
   }
   fBuffer = storeBuffer;

   // get version of key
   bufferRef.SetBufferOffset(sizeof(fNbytes));
   Version_t kvers = bufferRef.ReadVersion();

   bufferRef.SetBufferOffset(fKeylen);
   TClass *cl = TClass::GetClass(fClassName.Data());
   TClass *clOnfile = 0;
   if (!cl) {
      Error("ReadObjectAny", "Unknown class %s", fClassName.Data());
      return 0;
   }
   Int_t baseOffset = 0;
   if (expectedClass) {
       // baseOffset will be -1 if cl does not inherit from expectedClass
      baseOffset = cl->GetBaseClassOffset(expectedClass);
      if (baseOffset == -1) {
         // The 2 classes are unrelated, maybe there is a converter between the 2.

         if (!expectedClass->GetSchemaRules() ||
             !expectedClass->GetSchemaRules()->HasRuleWithSourceClass(cl->GetName()))
         {
            // There is no converter
            return 0;
         }
         baseOffset = 0; // For now we do not support requesting from a class that is the base of one of the class for which there is transformation to ....
         clOnfile = cl;
         cl = const_cast<TClass*>(expectedClass);
         Info("ReadObjectAny","Using Converter StreamerInfo from %s to %s",clOnfile->GetName(),expectedClass->GetName());
      }
      if (cl->GetState() > TClass::kEmulated && expectedClass->GetState() <= TClass::kEmulated) {
         //we cannot mix a compiled class with an emulated class in the inheritance
         Warning("ReadObjectAny",
                 "Trying to read an emulated class (%s) to store in a compiled pointer (%s)",
                 cl->GetName(),expectedClass->GetName());
      }
   }
   // Create an instance of this class

   void *pobj = cl->New();
   if (!pobj) {
      Error("ReadObjectAny", "Cannot create new object of class %s", fClassName.Data());
      return 0;
   }

   if (kvers > 1)
      bufferRef.MapObject(pobj,cl);  //register obj in map to handle self reference

   if (fObjlen > fNbytes-fKeylen) {
      char *objbuf = bufferRef.Buffer() + fKeylen;
      UChar_t *bufcur = (UChar_t *)&compressedBuffer[fKeylen];
      Int_t nin, nout = 0, nbuf;
      Int_t noutot = 0;
      while (1) {
         Int_t hc = R__unzip_header(&nin, bufcur, &nbuf);
         if (hc!=0) break;
         R__unzip(&nin, bufcur, &nbuf, (unsigned char*) objbuf, &nout);
         if (!nout) break;
         noutot += nout;
         if (noutot >= fObjlen) break;
         bufcur += nin;
         objbuf += nout;
      }
      if (nout) {
         cl->Streamer((void*)pobj, bufferRef, clOnfile);    //read object
      } else {
         cl->Destructor(pobj);
         return nullptr;
      }
   } else {
      cl->Streamer((void*)pobj, bufferRef, clOnfile);    //read object
   }

   if (cl->IsTObject()) {
      auto tobjBaseOffset = cl->GetBaseClassOffset(TObject::Class());
      if (tobjBaseOffset == -1) {
         Fatal("ReadObj","Incorrect detection of the inheritance from TObject for class %s.\n",
               fClassName.Data());
      }
      TObject *tobj = (TObject*)( ((char*)pobj) + tobjBaseOffset);

      // See similar adjustments in ReadObj
      if (gROOT->GetForceStyle()) tobj->UseCurrentStyle();

      if (cl->InheritsFrom(TDirectoryFile::Class())) {
         TDirectory *dir = static_cast<TDirectoryFile*>(tobj);
         dir->SetName(GetName());
         dir->SetTitle(GetTitle());
         dir->SetMother(fMotherDir);
         fMotherDir->Append(dir);
      }
   }

   {
      // Append the object to the directory if requested:
      ROOT::DirAutoAdd_t addfunc = cl->GetDirectoryAutoAdd();
      if (addfunc) {
         addfunc(pobj, fMotherDir);
      }
   }

   return ( ((char*)pobj) + baseOffset );
}

////////////////////////////////////////////////////////////////////////////////
/// To read an object from the file.
///
/// The object associated to this key is read from the file into memory.
/// Before invoking this function, obj has been created via the
/// default constructor.

Int_t TKey::Read(TObject *obj)
{
   if (!obj || (GetFile()==0)) return 0;

   TBufferFile bufferRef(TBuffer::kRead, fObjlen+fKeylen);
   bufferRef.SetParent(GetFile());
   bufferRef.SetPidOffset(fPidOffset);

   if (fVersion > 1)
      bufferRef.MapObject(obj);  //register obj in map to handle self reference

   std::unique_ptr<char []> compressedBuffer;
   auto storeBuffer = fBuffer;
   if (fObjlen > fNbytes-fKeylen) {
      compressedBuffer.reset(new char[fNbytes]);
      fBuffer = compressedBuffer.get();
      ReadFile();                    //Read object structure from file
      memcpy(bufferRef.Buffer(),fBuffer,fKeylen);
   } else {
      fBuffer = bufferRef.Buffer();
      ReadFile();                    //Read object structure from file
   }
   fBuffer = storeBuffer;

   bufferRef.SetBufferOffset(fKeylen);
   if (fObjlen > fNbytes-fKeylen) {
      char *objbuf = bufferRef.Buffer() + fKeylen;
      UChar_t *bufcur = (UChar_t *)&compressedBuffer[fKeylen];
      Int_t nin, nout = 0, nbuf;
      Int_t noutot = 0;
      while (1) {
         Int_t hc = R__unzip_header(&nin, bufcur, &nbuf);
         if (hc!=0) break;
         R__unzip(&nin, bufcur, &nbuf, (unsigned char*) objbuf, &nout);
         if (!nout) break;
         noutot += nout;
         if (noutot >= fObjlen) break;
         bufcur += nin;
         objbuf += nout;
      }
      if (nout) obj->Streamer(bufferRef);
   } else {
      obj->Streamer(bufferRef);
   }

   // Append the object to the directory if requested:
   {
      ROOT::DirAutoAdd_t addfunc = obj->IsA()->GetDirectoryAutoAdd();
      if (addfunc) {
         addfunc(obj, fMotherDir);
      }
   }

   return fNbytes;
}

////////////////////////////////////////////////////////////////////////////////
/// Decode input buffer.
///
/// In some situation will add key to gDirectory.

void TKey::ReadBuffer(char *&buffer)
{
   ReadKeyBuffer(buffer);

   if (!gROOT->ReadingObject() && gDirectory) {
      if (fSeekPdir != gDirectory->GetSeekDir()) gDirectory->AppendKey(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Decode input buffer.

void TKey::ReadKeyBuffer(char *&buffer)
{
   frombuf(buffer, &fNbytes);
   Version_t version;
   frombuf(buffer,&version);
   fVersion = (Int_t)version;
   frombuf(buffer, &fObjlen);
   fDatime.ReadBuffer(buffer);
   frombuf(buffer, &fKeylen);
   frombuf(buffer, &fCycle);
   if (fVersion > 1000) {
      frombuf(buffer, &fSeekKey);

      // We currently store in the 16 highest bit of fSeekPdir the value of
      // fPidOffset.  This offset is used when a key (or basket) is transfered from one
      // file to the other.  In this case the TRef and TObject might have stored a
      // pid index (to retrieve TProcessIDs) which refered to their order on the original
      // file, the fPidOffset is to be added to those values to correctly find the
      // TProcessID.  This fPidOffset needs to be increment if the key/basket is copied
      // and need to be zero for new key/basket.
      Long64_t pdir;
      frombuf(buffer, &pdir);
      fPidOffset = pdir >> kPidOffsetShift;
      fSeekPdir = pdir & kPidOffsetMask;
   } else {
      UInt_t seekkey,seekdir;
      frombuf(buffer, &seekkey); fSeekKey = (Long64_t)seekkey;
      frombuf(buffer, &seekdir); fSeekPdir= (Long64_t)seekdir;
   }
   fClassName.ReadBuffer(buffer);
   //the following test required for forward and backward compatibility
   if (fClassName == "TDirectory") {
      fClassName = "TDirectoryFile";
      SetBit(kIsDirectoryFile);
   }

   fName.ReadBuffer(buffer);
   fTitle.ReadBuffer(buffer);
}

////////////////////////////////////////////////////////////////////////////////
/// Read the key structure from the file

Bool_t TKey::ReadFile()
{
   TFile* f = GetFile();
   if (f==0) return kFALSE;

   Int_t nsize = fNbytes;
   f->Seek(fSeekKey);
#if 0
   for (Int_t i = 0; i < nsize; i += kMAXFILEBUFFER) {
      int nb = kMAXFILEBUFFER;
      if (i+nb > nsize) nb = nsize - i;
      f->ReadBuffer(fBuffer+i,nb);
   }
#else
   if( f->ReadBuffer(fBuffer,nsize) )
   {
      Error("ReadFile", "Failed to read data.");
      return kFALSE;
   }
#endif
   if (gDebug) {
      std::cout << "TKey Reading "<<nsize<< " bytes at address "<<fSeekKey<<std::endl;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parent in key buffer.

void TKey::SetParent(const TObject *parent)
{
   if (fBufferRef) fBufferRef->SetParent((TObject*)parent);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the key as it had not been 'filled' yet.

void TKey::Reset()
{
   fPidOffset  = 0;
   fNbytes     = 0;
   fBuffer     = 0;
   fObjlen     = 0;
   fCycle      = 0;
   fSeekPdir   = 0;
   fSeekKey    = 0;
   fLeft       = 0;
   fDatime     = (UInt_t)0;

   // fBufferRef and fKeylen intentionally not reset/changed

   keyAbsNumber++; SetUniqueID(keyAbsNumber);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the size in bytes of the key header structure.
///
/// An explaination about the nbytes (Int_t nbytes) variable used in the
/// function. The size of fSeekKey and fSeekPdir is 8 instead of 4 if version is
/// greater than 1000.
/// | Component         | Sizeof |
/// |-------------------|--------|
/// | fNbytes           | 4      |
/// | sizeof(Version_t) | 2      |
/// | fObjlen           | 4      |
/// | fKeylen           | 2      |
/// | fCycle            | 2      |
/// | fSeekKey          | 4 or 8 |
/// | fSeekPdir         | 4 or 8 |
/// | **TOTAL**         |   22   |

Int_t TKey::Sizeof() const
{
   Int_t nbytes = 22; if (fVersion > 1000) nbytes += 8;
   nbytes      += fDatime.Sizeof();
   if (TestBit(kIsDirectoryFile)) {
      nbytes   += 11; // strlen("TDirectory")+1
   } else {
      nbytes   += fClassName.Sizeof();
   }
   nbytes      += fName.Sizeof();
   nbytes      += fTitle.Sizeof();
   return nbytes;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TKey::Streamer(TBuffer &b)
{
   Version_t version;
   if (b.IsReading()) {
      b >> fNbytes;
      b >> version; fVersion = (Int_t)version;
      b >> fObjlen;
      fDatime.Streamer(b);
      b >> fKeylen;
      b >> fCycle;
      if (fVersion > 1000) {
         b >> fSeekKey;

         // We currently store in the 16 highest bit of fSeekPdir the value of
         // fPidOffset.  This offset is used when a key (or basket) is transfered from one
         // file to the other.  In this case the TRef and TObject might have stored a
         // pid index (to retrieve TProcessIDs) which refered to their order on the original
         // file, the fPidOffset is to be added to those values to correctly find the
         // TProcessID.  This fPidOffset needs to be increment if the key/basket is copied
         // and need to be zero for new key/basket.
         Long64_t pdir;
         b >> pdir;
         fPidOffset = pdir >> kPidOffsetShift;
         fSeekPdir = pdir & kPidOffsetMask;
      } else {
         UInt_t seekkey, seekdir;
         b >> seekkey; fSeekKey = (Long64_t)seekkey;
         b >> seekdir; fSeekPdir= (Long64_t)seekdir;
      }
      fClassName.Streamer(b);
      //the following test required for forward and backward compatibility
      if (fClassName == "TDirectory") {
         fClassName = "TDirectoryFile";
         SetBit(kIsDirectoryFile);
      }
      fName.Streamer(b);
      fTitle.Streamer(b);
      if (fKeylen < 0) {
         Error("Streamer","The value of fKeylen is incorrect (%d) ; trying to recover by setting it to zero",fKeylen);
         MakeZombie();
         fKeylen = 0;
      }
      if (fObjlen < 0) {
         Error("Streamer","The value of fObjlen is incorrect (%d) ; trying to recover by setting it to zero",fObjlen);
         MakeZombie();
         fObjlen = 0;
      }
      if (fNbytes < 0) {
         Error("Streamer","The value of fNbytes is incorrect (%d) ; trying to recover by setting it to zero",fNbytes);
         MakeZombie();
         fNbytes = 0;
      }

   } else {
      b << fNbytes;
      version = (Version_t)fVersion;
      b << version;
      b << fObjlen;
      if (fDatime.Get() == 0) fDatime.Set();
      if (TestBit(TKey::kReproducible))
         TDatime((UInt_t) 1).Streamer(b);
      else
         fDatime.Streamer(b);
      b << fKeylen;
      b << fCycle;
      if (fVersion > 1000) {
         b << fSeekKey;

         // We currently store in the 16 highest bit of fSeekPdir the value of
         // fPidOffset.  This offset is used when a key (or basket) is transfered from one
         // file to the other.  In this case the TRef and TObject might have stored a
         // pid index (to retrieve TProcessIDs) which refered to their order on the original
         // file, the fPidOffset is to be added to those values to correctly find the
         // TProcessID.  This fPidOffset needs to be increment if the key/basket is copied
         // and need to be zero for new key/basket.
         Long64_t pdir = (((Long64_t)fPidOffset)<<kPidOffsetShift) | (kPidOffsetMask & fSeekPdir);
         b << pdir;
      } else {
         b << (Int_t)fSeekKey;
         b << (Int_t)fSeekPdir;
      }
      if (TestBit(kIsDirectoryFile)) {
         // We want to record "TDirectory" instead of TDirectoryFile so that the file can be read by ancient version of ROOT.
         b.WriteTString(gTDirectoryString);
      } else {
         fClassName.Streamer(b);
      }
      fName.Streamer(b);
      fTitle.Streamer(b);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write the encoded object supported by this key.
/// The function returns the number of bytes committed to the file.
/// If a write error occurs, the number of bytes returned is -1.

Int_t TKey::WriteFile(Int_t cycle, TFile* f)
{
   if (!f) f = GetFile();
   if (!f) return -1;

   Int_t nsize  = fNbytes;
   char *buffer = fBuffer;
   if (cycle) {
      fCycle = cycle;
      FillBuffer(buffer);
      buffer = fBuffer;
   }

   if (fLeft > 0) nsize += sizeof(Int_t);
   f->Seek(fSeekKey);
#if 0
   for (Int_t i=0;i<nsize;i+=kMAXFILEBUFFER) {
      Int_t nb = kMAXFILEBUFFER;
      if (i+nb > nsize) nb = nsize - i;
      f->WriteBuffer(buffer,nb);
      buffer += nb;
   }
#else
   Bool_t result = f->WriteBuffer(buffer,nsize);
#endif
   //f->Flush(); Flushing takes too much time.
   //            Let user flush the file when they want.
   if (gDebug) {
      std::cout <<"   TKey Writing "<<nsize<< " bytes at address "<<fSeekKey
           <<" for ID= " <<GetName()<<" Title= "<<GetTitle()<<std::endl;
   }

   DeleteBuffer();
   return result==kTRUE ? -1 : nsize;
}

////////////////////////////////////////////////////////////////////////////////
/// Write the encoded object supported by this key.
/// The function returns the number of bytes committed to the file.
/// If a write error occurs, the number of bytes returned is -1.

Int_t TKey::WriteFileKeepBuffer(TFile *f)
{
   if (!f) f = GetFile();
   if (!f) return -1;

   Int_t nsize  = fNbytes;
   char *buffer = fBuffer;

   if (fLeft > 0) nsize += sizeof(Int_t);
   f->Seek(fSeekKey);
#if 0
   for (Int_t i=0;i<nsize;i+=kMAXFILEBUFFER) {
      Int_t nb = kMAXFILEBUFFER;
      if (i+nb > nsize) nb = nsize - i;
      f->WriteBuffer(buffer,nb);
      buffer += nb;
   }
#else
   Bool_t result = f->WriteBuffer(buffer,nsize);
#endif
   //f->Flush(); Flushing takes too much time.
   //            Let user flush the file when they want.
   if (gDebug) {
      std::cout <<"   TKey Writing "<<nsize<< " bytes at address "<<fSeekKey
      <<" for ID= " <<GetName()<<" Title= "<<GetTitle()<<std::endl;
   }

   return result==kTRUE ? -1 : nsize;
}

////////////////////////////////////////////////////////////////////////////////
/// Title can keep 32x32 xpm thumbnail/icon of the parent object.

const char *TKey::GetIconName() const
{
   return (!fTitle.IsNull() && fTitle.BeginsWith("/* ") ?  fTitle.Data() : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns title (title can contain 32x32 xpm thumbnail/icon).

const char *TKey::GetTitle() const
{
   if (!fTitle.IsNull() && fTitle.BeginsWith("/* ")) { // title contains xpm thumbnail
      static TString ret;
      int start = fTitle.Index("/*") + 3;
      int stop = fTitle.Index("*/") - 1;
      ret = fTitle(start, stop - start);
      return ret.Data();
   }
   return fTitle.Data();
}
