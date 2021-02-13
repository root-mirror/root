// $Id$
// Author: Sergey Linev   22/12/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TRootSniffer
#define ROOT_TRootSniffer

#include "TNamed.h"
#include "TList.h"
#include <memory>
#include <string>

class TFolder;
class TKey;
class TBufferFile;
class TDataMember;
class THttpCallArg;
class TRootSnifferStore;
class TRootSniffer;

class TRootSnifferScanRec {

   friend class TRootSniffer;

protected:
   // different bits used to scan hierarchy
   enum {
      kScan = 0x0001,        ///< normal scan of hierarchy
      kExpand = 0x0002,      ///< expand of specified item - allowed to scan object members
      kSearch = 0x0004,      ///< search for specified item (only objects and collections)
      kCheckChilds = 0x0008, ///< check if there childs, very similar to search
      kOnlyFields = 0x0010,  ///< if set, only fields for specified item will be set (but all fields)
      kActions = 0x001F      ///< mask for actions, only actions copied to child rec
   };

   TRootSnifferScanRec *fParent{nullptr}; ///<! pointer on parent record
   UInt_t fMask{0};                       ///<! defines operation kind
   const char *fSearchPath{nullptr};      ///<! current path searched
   Int_t fLevel{0};                       ///<! current level of hierarchy
   TString fItemName;                     ///<! name of current item
   TList fItemsNames;                     ///<! list of created items names, need to avoid duplication
   Int_t fRestriction{0};                 ///<! restriction 0 - default, 1 - read-only, 2 - full access

   TRootSnifferStore *fStore{nullptr}; ///<! object to store results
   Bool_t fHasMore{kFALSE};            ///<! indicates that potentially there are more items can be found
   Bool_t fNodeStarted{kFALSE};        ///<! indicate if node was started
   Int_t fNumFields{0};                ///<! number of fields
   Int_t fNumChilds{0};                ///<! number of childs

public:
   TRootSnifferScanRec();
   virtual ~TRootSnifferScanRec();

   void CloseNode();

   /** return true when fields could be set to the hierarchy item */
   Bool_t CanSetFields() const { return (fMask & kScan) && (fStore != nullptr); }

   /** return true when only fields are scanned by the sniffer */
   Bool_t ScanOnlyFields() const { return (fMask & kOnlyFields) && (fMask & kScan); }

   /** Starts new node, must be closed at the end */
   void CreateNode(const char *_node_name);

   void BeforeNextChild();

   /** Set item field only when creating is specified */
   void SetField(const char *name, const char *value, Bool_t with_quotes = kTRUE);

   /** Mark item with ROOT class and correspondent streamer info */
   void SetRootClass(TClass *cl);

   /** Returns true when item can be expanded */
   Bool_t CanExpandItem();

   /** Checks if result will be accepted. Used to verify if sniffer should read object from the file */
   Bool_t IsReadyForResult() const;

   /** Obsolete, use SetFoundResult instead */
   Bool_t SetResult(void *obj, TClass *cl, TDataMember *member = nullptr);

   /** Set found element with class and datamember (optional) */
   Bool_t SetFoundResult(void *obj, TClass *cl, TDataMember *member = nullptr);

   /** Returns depth of hierarchy */
   Int_t Depth() const;

   /** Method indicates that scanning can be interrupted while result is set */
   Bool_t Done() const;

   /** Construct item name, using object name as basis */
   void MakeItemName(const char *objname, TString &itemname);

   /** Produces full name for the current item */
   void BuildFullName(TString &buf, TRootSnifferScanRec *prnt = nullptr);

   /** Returns read-only flag for current item */
   Bool_t IsReadOnly(Bool_t dflt = kTRUE);

   Bool_t GoInside(TRootSnifferScanRec &super, TObject *obj, const char *obj_name = nullptr,
                   TRootSniffer *sniffer = nullptr);

   ClassDef(TRootSnifferScanRec, 0) // Scan record for objects sniffer
};

//_______________________________________________________________________

class TRootSniffer : public TNamed {
   enum {
      kItemField = BIT(21) // item property stored as TNamed
   };

protected:
   TString fObjectsPath;    ///<! default path for registered objects
   Bool_t fReadOnly{kTRUE}; ///<! indicate if sniffer allowed to change ROOT structures - like read objects from file
   Bool_t fScanGlobalDir{kTRUE};       ///<! when enabled (default), scan gROOT for histograms, canvases, open files
   std::unique_ptr<TFolder> fTopFolder; ///<! own top TFolder object, used for registering objects
   THttpCallArg *fCurrentArg{nullptr}; ///<! current http arguments (if any)
   Int_t fCurrentRestrict{0};          ///<! current restriction for last-found object
   TString fCurrentAllowedMethods;     ///<! list of allowed methods, extracted when analyzed object restrictions
   TList fRestrictions;                ///<! list of restrictions for different locations
   TString fAutoLoad;                  ///<! scripts names, which are add as _autoload parameter to h.json request

   void ScanObjectMembers(TRootSnifferScanRec &rec, TClass *cl, char *ptr);

   virtual void ScanObjectProperties(TRootSnifferScanRec &rec, TObject *obj);

   virtual void ScanKeyProperties(TRootSnifferScanRec &rec, TKey *key, TObject *&obj, TClass *&obj_class);

   virtual void ScanObjectChilds(TRootSnifferScanRec &rec, TObject *obj);

   void
   ScanCollection(TRootSnifferScanRec &rec, TCollection *lst, const char *foldername = nullptr, TCollection *keys_lst = nullptr);

   virtual void ScanRoot(TRootSnifferScanRec &rec);

   TString DecodeUrlOptionValue(const char *value, Bool_t remove_quotes = kTRUE);

   TObject *GetItem(const char *fullname, TFolder *&parent, Bool_t force = kFALSE, Bool_t within_objects = kTRUE);

   TFolder *GetSubFolder(const char *foldername, Bool_t force = kFALSE);

   const char *GetItemField(TFolder *parent, TObject *item, const char *name);

   Bool_t IsItemField(TObject *obj) const;

   Bool_t AccessField(TFolder *parent, TObject *item, const char *name, const char *value, TNamed **only_get = nullptr);

   Int_t WithCurrentUserName(const char *option);

   virtual Bool_t CanDrawClass(TClass *) { return kFALSE; }

   virtual Bool_t HasStreamerInfo() const { return kFALSE; }

   virtual Bool_t ProduceJson(const std::string &path, const std::string &options, std::string &res);

   virtual Bool_t ProduceXml(const std::string &path, const std::string &options, std::string &res);

   virtual Bool_t ProduceBinary(const std::string &path, const std::string &options, std::string &res);

   virtual Bool_t ProduceImage(Int_t kind, const std::string &path, const std::string &options, std::string &res);

   virtual Bool_t ProduceExe(const std::string &path, const std::string &options, Int_t reskind, std::string &res);

   virtual Bool_t ExecuteCmd(const std::string &path, const std::string &options, std::string &res);

   virtual Bool_t
   ProduceItem(const std::string &path, const std::string &options, std::string &res, Bool_t asjson = kTRUE);

   virtual Bool_t
   ProduceMulti(const std::string &path, const std::string &options, std::string &res, Bool_t asjson = kTRUE);

public:
   TRootSniffer(const char *name, const char *objpath = "Objects");
   virtual ~TRootSniffer();

   /** When readonly on (default), sniffer is not allowed to change ROOT structures
     * For instance, it is not allowed to read new objects from files */
   void SetReadOnly(Bool_t on = kTRUE) { fReadOnly = on; }

   /** Returns readonly mode */
   Bool_t IsReadOnly() const { return fReadOnly; }

   void Restrict(const char *path, const char *options);

   Bool_t HasRestriction(const char *item_name);

   Int_t CheckRestriction(const char *item_name);

   void CreateOwnTopFolder();

   TFolder *GetTopFolder(Bool_t force = kFALSE);

   /** When enabled (default), sniffer scans gROOT for files, canvases, histograms */
   void SetScanGlobalDir(Bool_t on = kTRUE) { fScanGlobalDir = on; }

   void SetAutoLoad(const char *scripts = "");

   const char *GetAutoLoad() const;

   /** Returns true when sniffer allowed to scan global directories */
   Bool_t IsScanGlobalDir() const { return fScanGlobalDir; }

   Bool_t RegisterObject(const char *subfolder, TObject *obj);

   Bool_t UnregisterObject(TObject *obj);

   Bool_t RegisterCommand(const char *cmdname, const char *method, const char *icon);

   Bool_t CreateItem(const char *fullname, const char *title);

   Bool_t SetItemField(const char *fullname, const char *name, const char *value);

   const char *GetItemField(const char *fullname, const char *name);

   void SetCurrentCallArg(THttpCallArg *arg);

   /** Method scans normal objects, registered in ROOT */
   void ScanHierarchy(const char *topname, const char *path, TRootSnifferStore *store, Bool_t only_fields = kFALSE);

   TObject *FindTObjectInHierarchy(const char *path);

   virtual void *
   FindInHierarchy(const char *path, TClass **cl = nullptr, TDataMember **member = nullptr, Int_t *chld = nullptr);

   Bool_t CanDrawItem(const char *path);

   Bool_t CanExploreItem(const char *path);

   virtual Bool_t IsStreamerInfoItem(const char *) { return kFALSE; }

   virtual ULong_t GetStreamerInfoHash() { return 0; }

   virtual ULong_t GetItemHash(const char *itemname);

   Bool_t Produce(const std::string &path, const std::string &file, const std::string &options, std::string &res);

   ClassDefOverride(TRootSniffer, 0) // Sniffer of ROOT objects (basic version)
};

#endif
