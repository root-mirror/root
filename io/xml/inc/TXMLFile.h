// @(#)root/xml:$Id$
// Author: Sergey Linev  10.05.2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TXMLFile
#define ROOT_TXMLFile

#include "TFile.h"
#include "TXMLSetup.h"
#include "TXMLEngine.h"
#include "Compression.h"
#include <memory>

class TKeyXML;
class TStreamerElement;
class TStreamerInfo;

class TXMLFile final : public TFile, public TXMLSetup {

protected:
   void InitXmlFile(Bool_t create);
   // Interface to basic system I/O routines
   Int_t SysOpen(const char *, Int_t, UInt_t) final { return 0; }
   Int_t SysClose(Int_t) final { return 0; }
   Int_t SysRead(Int_t, void *, Int_t) final { return 0; }
   Int_t SysWrite(Int_t, const void *, Int_t) final { return 0; }
   Long64_t SysSeek(Int_t, Long64_t, Int_t) final { return 0; }
   Int_t SysStat(Int_t, Long_t *, Long64_t *, Long_t *, Long_t *) final { return 0; }
   Int_t SysSync(Int_t) final { return 0; }

   // Overwrite methods for directory I/O
   Long64_t DirCreateEntry(TDirectory *) final;
   Int_t DirReadKeys(TDirectory *) final;
   void DirWriteKeys(TDirectory *) final;
   void DirWriteHeader(TDirectory *) final;

   InfoListRet GetStreamerInfoListImpl(bool lookupSICache) final;

private:
   TXMLFile(const TXMLFile &) = delete;            // TXMLFile cannot be copied, not implemented
   void operator=(const TXMLFile &) = delete;      // TXMLFile cannot be copied, not implemented

public:
   TXMLFile() {} // NOLINT: not allowed to use = default because of TObject::kIsOnHeap detection, see ROOT-10300
   TXMLFile(const char *filename, Option_t *option = "read", const char *title = "title", Int_t compression = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault);
   virtual ~TXMLFile();

   void Close(Option_t *option = "") final; // *MENU*
   TKey *CreateKey(TDirectory *mother, const TObject *obj, const char *name, Int_t bufsize) final;
   TKey *CreateKey(TDirectory *mother, const void *obj, const TClass *cl, const char *name, Int_t bufsize) final;
   void DrawMap(const char * = "*", Option_t * = "") final {}
   void FillBuffer(char *&) final {}
   void Flush() final {}

   Long64_t GetEND() const final { return 0; }
   Int_t GetErrno() const final { return 0; }
   void ResetErrno() const final {}

   Int_t GetNfree() const final { return 0; }
   Int_t GetNbytesInfo() const final { return 0; }
   Int_t GetNbytesFree() const final { return 0; }
   Long64_t GetSeekFree() const final { return 0; }
   Long64_t GetSeekInfo() const final { return 0; }
   Long64_t GetSize() const final { return 0; }

   Int_t GetIOVersion() const { return fIOVersion; }

   Bool_t IsOpen() const final;

   void MakeFree(Long64_t, Long64_t) final {}
   void MakeProject(const char *, const char * = "*", Option_t * = "new") final {} // *MENU*
   void Map(Option_t *) final {}                                                   //
   void Map() final {}                                                             //
   void Paint(Option_t * = "") final {}
   void Print(Option_t * = "") const final {}
   Bool_t ReadBuffer(char *, Int_t) final { return kFALSE; }
   Bool_t ReadBuffer(char *, Long64_t, Int_t) final { return kFALSE; }
   void ReadFree() final {}
   Int_t Recover() final { return 0; }
   Int_t ReOpen(Option_t *mode) final;
   void Seek(Long64_t, ERelativeTo = kBeg) final {}

   void SetEND(Long64_t) final {}
   Int_t Sizeof() const final { return 0; }

   Bool_t WriteBuffer(const char *, Int_t) final { return kFALSE; }
   Int_t Write(const char * = nullptr, Int_t = 0, Int_t = 0) final { return 0; }
   Int_t Write(const char * = nullptr, Int_t = 0, Int_t = 0) const final { return 0; }
   void WriteFree() final {}
   void WriteHeader() final {}
   void WriteStreamerInfo() final;

   // XML specific functions

   void SetXmlLayout(EXMLLayout layout) final;
   void SetStoreStreamerInfos(Bool_t iConvert = kTRUE) final;
   void SetUsedDtd(Bool_t use = kTRUE) final;
   void SetUseNamespaces(Bool_t iUseNamespaces = kTRUE) final;

   Bool_t AddXmlComment(const char *comment);
   Bool_t AddXmlStyleSheet(const char *href, const char *type = "text/css", const char *title = nullptr, int alternate = -1,
                           const char *media = nullptr, const char *charset = nullptr);
   Bool_t AddXmlLine(const char *line);

   TXMLEngine *XML() { return fXML.get(); }

protected:
   // functions to store streamer infos

   void StoreStreamerElement(XMLNodePointer_t node, TStreamerElement *elem);
   void ReadStreamerElement(XMLNodePointer_t node, TStreamerInfo *info);

   Bool_t ReadFromFile();
   Int_t ReadKeysList(TDirectory *dir, XMLNodePointer_t topnode);
   TKeyXML *FindDirKey(TDirectory *dir);
   TDirectory *FindKeyDir(TDirectory *mother, Long64_t keyid);
   void CombineNodesTree(TDirectory *dir, XMLNodePointer_t topnode, Bool_t dolink);

   void SaveToFile();

   static void ProduceFileNames(const char *filename, TString &fname, TString &dtdname);

   XMLDocPointer_t fDoc{nullptr}; //!

   XMLNodePointer_t fStreamerInfoNode{nullptr}; //!  pointer of node with streamer info data

   std::unique_ptr<TXMLEngine> fXML; //! object for interface with xml library

   Int_t fIOVersion{0}; //! indicates format of ROOT xml file

   Long64_t fKeyCounter{0}; //! counter of created keys, used for keys id

   ClassDefOverride(TXMLFile, 3) // ROOT file in XML format
};

#endif
