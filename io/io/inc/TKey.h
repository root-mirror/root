// @(#)root/io
// Author: Rene Brun   28/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TKey
#define ROOT_TKey

#include "TNamed.h"
#include "TDatime.h"
#include "TClass.h"
#ifdef R__LESS_INCLUDES
class TBuffer;
#else
#include "TBuffer.h"
#endif

class TBrowser;
class TDirectory;
class TFile;

class TKey : public TNamed {

private:
   enum EStatusBits {
      kIsDirectoryFile = BIT(14),
      kReproducible = BIT(15)
   };
   TKey(const TKey&) = delete;            // TKey objects are not copiable.
   TKey& operator=(const TKey&) = delete; // TKey objects are not copiable.

protected:
   Int_t       fVersion;     ///< Key version identifier
   Int_t       fNbytes;      ///< Number of bytes for the object on file
   Int_t       fObjlen;      ///< Length of uncompressed object in bytes
   TDatime     fDatime;      ///< Date/Time of insertion in file
   Short_t     fKeylen;      ///< Number of bytes for the key itself
   Short_t     fCycle;       ///< Cycle number
   Long64_t    fSeekKey;     ///< Location of object on file
   Long64_t    fSeekPdir;    ///< Location of parent directory on file
   TString     fClassName;   ///< Object Class name
   Int_t       fLeft;        ///< Number of bytes left in current segment
   char       *fBuffer;      ///< Object buffer
   TBuffer    *fBufferRef;   ///< Pointer to the TBuffer object
   UShort_t    fPidOffset;   ///<!Offset to be added to the pid index in this key/buffer.  This is actually saved in the high bits of fSeekPdir
   TDirectory *fMotherDir;   ///<!pointer to mother directory

   virtual Int_t    Read(const char *name) { return TObject::Read(name); }
   virtual void     Create(Int_t nbytes, TFile* f = 0);
           void     Build(TDirectory* motherDir, const char* classname, Long64_t filepos);
           void     Reset(); // Currently only for the use of TBasket.
   virtual Int_t    WriteFileKeepBuffer(TFile *f = 0);


 public:
   TKey();
   TKey(TDirectory* motherDir);
   TKey(TDirectory* motherDir, const TKey &orig, UShort_t pidOffset);
   TKey(const char *name, const char *title, const TClass *cl, Int_t nbytes, TDirectory* motherDir);
   TKey(const TString &name, const TString &title, const TClass *cl, Int_t nbytes, TDirectory* motherDir);
   TKey(const TObject *obj, const char *name, Int_t bufsize, TDirectory* motherDir);
   TKey(const void *obj, const TClass *cl, const char *name, Int_t bufsize, TDirectory* motherDir);
   TKey(Long64_t pointer, Int_t nbytes, TDirectory* motherDir = 0);
   virtual ~TKey();

   virtual void        Browse(TBrowser *b);
   virtual void        Delete(Option_t *option="");
   virtual void        DeleteBuffer();
   virtual void        FillBuffer(char *&buffer);
   virtual const char *GetClassName() const {return fClassName.Data();}
   virtual const char *GetIconName() const;
   virtual const char *GetTitle() const;
   virtual char       *GetBuffer() const {return fBuffer+fKeylen;}
           TBuffer    *GetBufferRef() const {return fBufferRef;}
           Short_t     GetCycle() const;
   const   TDatime    &GetDatime() const   {return fDatime;}
           TFile      *GetFile() const;
           Short_t     GetKeep() const;
           Int_t       GetKeylen() const   {return fKeylen;}
           TDirectory* GetMotherDir() const { return fMotherDir; }
           Int_t       GetNbytes() const   {return fNbytes;}
           Int_t       GetObjlen() const   {return fObjlen;}
           Int_t       GetVersion() const  {return fVersion;}
   virtual Long64_t    GetSeekKey() const  {return fSeekKey;}
   virtual Long64_t    GetSeekPdir() const {return fSeekPdir;}
   virtual void        IncrementPidOffset(UShort_t offset);
           Bool_t      IsFolder() const;
   virtual void        Keep();
   virtual void        ls(Bool_t current) const;
   virtual void        ls(Option_t *option="") const;
   virtual void        Print(Option_t *option="") const;
   virtual Int_t       Read(TObject *obj);
   virtual TObject    *ReadObj();
   virtual TObject    *ReadObjWithBuffer(char *bufferRead);
   /// To read an object (non deriving from TObject) from the file.
   /// This is more user friendly version of TKey::ReadObjectAny.
   /// See TKey::ReadObjectAny for more details.
   template <typename T> T *ReadObject() {
      return reinterpret_cast<T*>(ReadObjectAny(TClass::GetClass<T>()));
   }
   virtual void       *ReadObjectAny(const TClass *expectedClass);
   virtual void        ReadBuffer(char *&buffer);
           void        ReadKeyBuffer(char *&buffer);
   virtual Bool_t      ReadFile();
   virtual void        SetBuffer() { DeleteBuffer(); fBuffer = new char[fNbytes];}
   virtual void        SetParent(const TObject *parent);
           void        SetMotherDir(TDirectory* dir) { fMotherDir = dir; }
   virtual Int_t       Sizeof() const;
   virtual Int_t       WriteFile(Int_t cycle=1, TFile* f = 0);

   ClassDef(TKey,4); //Header description of a logical record on file.
};

#endif
