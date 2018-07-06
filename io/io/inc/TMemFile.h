// @(#)root/io:$Id$
// Author: Philippe Canal, May 2011

/*************************************************************************
 * Copyright (C) 1995-2009, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMemFile
#define ROOT_TMemFile

#include "TFile.h"
#include <vector>
#include <memory>

class TMemFile : public TFile {
public:
   using ExternalDataPtr_t = std::shared_ptr<const std::vector<char>>;

private:
   struct TMemBlock {
   private:
      TMemBlock(const TMemBlock&);            // Not implemented
      TMemBlock &operator=(const TMemBlock&); // Not implemented.
   public:
      TMemBlock();
      TMemBlock(Long64_t size, TMemBlock *previous = 0);
      TMemBlock(UChar_t* externalBuffer, Long64_t size);
      ~TMemBlock();

      void CreateNext(Long64_t size);

      TMemBlock *fPrevious;
      TMemBlock *fNext;
      UChar_t   *fBuffer;
      Long64_t   fSize;
   };
   TMemBlock    fBlockList;               ///< Collection of memory blocks of size fgDefaultBlockSize
   const ExternalDataPtr_t fExternalData; ///< shared file data / content
   Long64_t     fSize;                    ///< Total file size (sum of the size of the chunks)
   Long64_t     fSysOffset;               ///< Seek offset in file
   TMemBlock   *fBlockSeek;               ///< Pointer to the block we seeked to.
   Long64_t     fBlockOffset;             ///< Seek offset within the block

   static Long64_t fgDefaultBlockSize;

   Long64_t MemRead(Int_t fd, void *buf, Long64_t len) const;

   // Overload TFile interfaces.
   Int_t    SysOpen(const char *pathname, Int_t flags, UInt_t mode);
   Int_t    SysClose(Int_t fd);
   Int_t    SysReadImpl(Int_t fd, void *buf, Long64_t len);
   Int_t    SysWriteImpl(Int_t fd, const void *buf, Long64_t len);
   Int_t    SysRead(Int_t fd, void *buf, Int_t len);
   Int_t    SysWrite(Int_t fd, const void *buf, Int_t len);
   Long64_t SysSeek(Int_t fd, Long64_t offset, Int_t whence);
   Int_t    SysStat(Int_t fd, Long_t *id, Long64_t *size, Long_t *flags, Long_t *modtime);
   Int_t    SysSync(Int_t fd);

   void ResetObjects(TDirectoryFile *, TFileMergeInfo *) const;

   enum EMode : int {
      kUnmanagedMask     = 0x10000,
      kCreate            = 0x00001,
      kRecreate          = 0x00010,
      kUpdate            = 0x00100,
      kRead              = 0x01000,
      kCreateUnmanaged   = 0x10001,
      kRecreateUnmanaged = 0x10010,
      kUpdateUnmanaged   = 0x10100,
      kReadUnmanaged     = 0x11000,
   };

   bool NeedsToWrite(EMode mode) const { return !(mode & EMode::kRead); }
   bool NeedsExistingFile(EMode mode) const { return (mode & EMode::kUpdate) || (mode & EMode::kRead);}

   EMode ParseOption(Option_t *option);

   TMemFile &operator=(const TMemFile&); // Not implemented.

public:
   TMemFile(const char *name, Option_t *option="", const char *ftitle="", Int_t compress=4);
   TMemFile(const char *name, char *buffer, Long64_t size, Option_t *option="", const char *ftitle="", Int_t compress=4);
   TMemFile(const char *name, ExternalDataPtr_t data);
   TMemFile(const TMemFile &orig);
   virtual ~TMemFile();

   virtual Long64_t CopyTo(void *to, Long64_t maxsize) const;
   virtual void     CopyTo(TBuffer &tobuf) const;
   virtual Long64_t GetSize() const;

   void ResetAfterMerge(TFileMergeInfo *);
   void ResetErrno() const;

   virtual void        Print(Option_t *option="") const;

   ClassDef(TMemFile, 0) // A ROOT file that reads/writes on a chunk of memory
};

#endif
