// @(#)root/base:$Id: 65b4f3646f4e5b2fa77218ba786b7fe4e16e27be $
// Author: Rene Brun   28/11/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#include <cstdlib>

#include "Strlen.h"
#include "strlcpy.h"
#include "TDirectory.h"
#include "TBuffer.h"
#include "TClassTable.h"
#include "TInterpreter.h"
#include "THashList.h"
#include "TBrowser.h"
#include "TROOT.h"
#include "TError.h"
#include "TClass.h"
#include "TRegexp.h"
#include "TSystem.h"
#include "TVirtualMutex.h"
#include "TThreadSlots.h"
#include "TMethod.h"

#include "TSpinLockGuard.h"

Bool_t TDirectory::fgAddDirectory = kTRUE;

const Int_t  kMaxLen = 2048;

/** \class TDirectory
\ingroup Base

Describe directory structure in memory.
*/

ClassImp(TDirectory);

////////////////////////////////////////////////////////////////////////////////
/// Directory default constructor.

TDirectory::TDirectory() : TNamed()
{
   // MSVC doesn't support fSpinLock=ATOMIC_FLAG_INIT; in the class definition
   std::atomic_flag_clear( &fSpinLock );
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new Directory.
///
///  A new directory with name,title is created in the current directory
///  The directory header information is immediately saved in the file
///  A new key is added in the parent directory
///
///  When this constructor is called from a class directly derived
///  from TDirectory, the third argument classname MUST be specified.
///  In this case, classname must be the name of the derived class.
///
///  Note that the directory name cannot contain slashes.

TDirectory::TDirectory(const char *name, const char *title, Option_t * /*classname*/, TDirectory* initMotherDir)
   : TNamed(name, title)
{
   // MSVC doesn't support fSpinLock=ATOMIC_FLAG_INIT; in the class definition
   std::atomic_flag_clear( &fSpinLock );

   if (!initMotherDir) initMotherDir = gDirectory;

   if (strchr(name,'/')) {
      ::Error("TDirectory::TDirectory","directory name (%s) cannot contain a slash", name);
      gDirectory = nullptr;
      return;
   }
   if (strlen(GetName()) == 0) {
      ::Error("TDirectory::TDirectory","directory name cannot be \"\"");
      gDirectory = nullptr;
      return;
   }

   BuildDirectory(initMotherDir ? initMotherDir->GetFile() : nullptr, initMotherDir);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TDirectory::~TDirectory()
{
   if (!gROOT) {
      delete fList;
      return; //when called by TROOT destructor
   }

   if (fList) {
      if (!fList->IsUsingRWLock())
         Fatal("~TDirectory","In %s:%p the fList (%p) is not using the RWLock\n",
               GetName(),this,fList);
      fList->Delete("slow");
      SafeDelete(fList);
   }

   TDirectory::CleanTargets();

   TDirectory* mom = GetMotherDir();

   if (mom) {
      mom->Remove(this);
   }

   if (gDebug) {
      Info("~TDirectory", "dtor called for %s", GetName());
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Set the current directory to null.
/// This is called from the TContext destructor.  Since the destructor is
/// inline, we do not want to have it directly use a global variable.

void TDirectory::TContext::CdNull()
{
   gDirectory = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.
///
/// Reset the current directory to its previous state.

TDirectory::TContext::~TContext()
{
   fActiveDestructor = true;
   if (fDirectory) {
      // UnregisterContext must not be virtual to allow
      // this to work even with fDirectory set to nullptr.
      (*fDirectory).UnregisterContext(this);
      // While we were waiting for the lock, the TDirectory
      // may have been deleted by another thread, so
      // we need to recheck the value of fDirectory.
      if (fDirectory)
         (*fDirectory).cd();
      else
         CdNull();
   } else {
      CdNull();
   }
   fActiveDestructor = false;
   while(fDirectoryWait);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the flag controlling the automatic add objects like histograms, TGraph2D, etc
/// in memory
///
/// By default (fAddDirectory = kTRUE), these objects are automatically added
/// to the list of objects in memory.
/// Note that in the classes like TH1, TGraph2D supporting this facility,
/// one object can be removed from its support directory
/// by calling object->SetDirectory(nullptr) or object->SetDirectory(dir) to add it
/// to the list of objects in the directory dir.
///
///  NOTE that this is a static function. To call it, use:
/// ~~~ {.cpp}
///     TDirectory::AddDirectory
/// ~~~

void TDirectory::AddDirectory(Bool_t add)
{
   fgAddDirectory = add;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function: see TDirectory::AddDirectory for more comments.

Bool_t TDirectory::AddDirectoryStatus()
{
   return fgAddDirectory;
}

////////////////////////////////////////////////////////////////////////////////
/// Append object to this directory.
///
/// If `replace` is true:
///   remove any existing objects with the same name (if the name is not "")

void TDirectory::Append(TObject *obj, Bool_t replace /* = kFALSE */)
{
   if (!obj || !fList) return;

   if (replace && obj->GetName() && obj->GetName()[0]) {
      TObject *old;
      while (nullptr != (old = GetList()->FindObject(obj->GetName()))) {
         Warning("Append","Replacing existing %s: %s (Potential memory leak).",
                 obj->IsA()->GetName(),obj->GetName());
         ROOT::DirAutoAdd_t func = old->IsA()->GetDirectoryAutoAdd();
         if (func) {
            func(old,nullptr);
         } else {
            Remove(old);
         }
      }
   }

   fList->Add(obj);
   obj->SetBit(kMustCleanup);
}

////////////////////////////////////////////////////////////////////////////////
/// Browse the content of the directory.

void TDirectory::Browse(TBrowser *b)
{
   if (b) {
      TObject *obj = nullptr;
      TIter nextin(fList);

      cd();

      //Add objects that are only in memory
      while ((obj = nextin())) {
         b->Add(obj, obj->GetName());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialise directory to defaults.
///
/// If directory is created via default ctor (when dir is read from file)
/// don't add it here to the directory since its name is not yet known.
/// It will be added to the directory in TKey::ReadObj().

void TDirectory::BuildDirectory(TFile* /*motherFile*/, TDirectory* motherDir)
{
   fList       = new THashList(100,50);
   fList->UseRWLock();
   fMother     = motherDir;
   SetBit(kCanDelete);

   // Build is done and is the last part of the constructor (and is not
   // being called from the derived classes) so we can publish.
   if (motherDir && strlen(GetName()) != 0) motherDir->Append(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Clean the pointers to this object (gDirectory, TContext, etc.).

void TDirectory::CleanTargets()
{
   std::vector<TContext*> extraWait;

   {
      ROOT::Internal::TSpinLockGuard slg(fSpinLock);

      while (fContext) {
         const auto next = fContext->fNext;
         const auto ctxt = fContext;
         ctxt->fDirectoryWait = true;

         ctxt->fDirectory = nullptr;

         if (ctxt->fActiveDestructor) {
            extraWait.push_back(fContext);
         } else {
            ctxt->fDirectoryWait = false;
         }
         fContext = next;
      }

      for(auto ptr : fGDirectories) {
         if (ptr->load() == this) {
            (*ptr) = nullptr;
         }
      }
   }
   for(auto &&context : extraWait) {
      // Wait until the TContext is done spinning
      // over the lock.
      while(context->fActiveDestructor);
      context->fDirectoryWait = false;
   }

   if (gDirectory == this) {
      TDirectory *cursav = GetMotherDir();
      if (cursav && cursav != this) {
         cursav->cd();
      } else {
         if (this == gROOT) {
            gDirectory = nullptr;
         } else {
            gROOT->cd();
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fast execution of 'new TBufferFile(TBuffer::kWrite,10000), without having
/// a compile time circular dependency ... alternatively we could (should?)
/// introduce yet another abstract interface.

static TBuffer* R__CreateBuffer()
{
   typedef void (*tcling_callfunc_Wrapper_t)(void*, int, void**, void*);
   static tcling_callfunc_Wrapper_t creator = nullptr;
   if (!creator) {
      R__LOCKGUARD(gROOTMutex);
      TClass *c = TClass::GetClass("TBufferFile");
      TMethod *m = c->GetMethodWithPrototype("TBufferFile","TBuffer::EMode,Int_t",kFALSE,ROOT::kExactMatch);
      creator = (tcling_callfunc_Wrapper_t)( m->InterfaceMethod() );
   }
   TBuffer::EMode mode = TBuffer::kWrite;
   Int_t size = 10000;
   void *args[] = { &mode, &size };
   TBuffer *result;
   creator(0,2,args,&result);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Clone an object.
/// This function is called when the directory is not a TDirectoryFile.
/// This version has to load the I/O package, hence via Cling.
///
/// If autoadd is true and if the object class has a
/// DirectoryAutoAdd function, it will be called at the end of the
/// function with the parameter gDirector.  This usually means that
/// the object will be appended to the current ROOT directory.

TObject *TDirectory::CloneObject(const TObject *obj, Bool_t autoadd /* = kTRUE */)
{
   // if no default ctor return immediately (error issued by New())
   char *pobj = (char*)obj->IsA()->New();
   if (!pobj) {
     Fatal("CloneObject","Failed to create new object");
     return nullptr;
   }

   Int_t baseOffset = obj->IsA()->GetBaseClassOffset(TObject::Class());
   if (baseOffset==-1) {
      // cl does not inherit from TObject.
      // Since this is not supported in this function, the only reason we could reach this code
      // is because something is screwed up in the ROOT code.
      Fatal("CloneObject","Incorrect detection of the inheritance from TObject for class %s.\n",
            obj->IsA()->GetName());
   }
   TObject *newobj = (TObject*)(pobj+baseOffset);

   //create a buffer where the object will be streamed
   //We are forced to go via the I/O package (ie TBufferFile).
   //Invoking TBufferFile via CINT will automatically load the I/O library
   TBuffer *buffer = R__CreateBuffer();
   if (!buffer) {
      Fatal("CloneObject","Not able to create a TBuffer!");
      return nullptr;
   }
   buffer->MapObject(obj);  //register obj in map to handle self reference
   const_cast<TObject*>(obj)->Streamer(*buffer);

   // read new object from buffer
   buffer->SetReadMode();
   buffer->ResetMap();
   buffer->SetBufferOffset(0);
   buffer->MapObject(newobj);  //register obj in map to handle self reference
   newobj->Streamer(*buffer);
   newobj->ResetBit(kIsReferenced);
   newobj->ResetBit(kCanDelete);

   delete buffer;
   if (autoadd) {
      ROOT::DirAutoAdd_t func = obj->IsA()->GetDirectoryAutoAdd();
      if (func) {
         func(newobj,this);
      }
   }
   return newobj;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current directory for the current thread.

std::atomic<TDirectory*> &TDirectory::CurrentDirectory()
{
   static std::atomic<TDirectory*> currentDirectory{nullptr};
   if (!gThreadTsd)
      return currentDirectory;
   else
      return *(std::atomic<TDirectory*>*)(*gThreadTsd)(&currentDirectory,ROOT::kDirectoryThreadSlot);
}

////////////////////////////////////////////////////////////////////////////////
/// Find a directory using apath.
/// It apath is null or empty, returns "this" directory.
/// Otherwise use apath to find a directory.
/// The absolute path syntax is: `file.root:/dir1/dir2`
///
/// where file.root is the file and /dir1/dir2 the desired subdirectory
/// in the file. Relative syntax is relative to "this" directory. E.g: `../aa`.
/// Returns 0 in case path does not exist.
/// If printError is true, use Error with 'funcname' to issue an error message.

TDirectory *TDirectory::GetDirectory(const char *apath,
                                     Bool_t printError, const char *funcname)
{
   Int_t nch = 0;
   if (apath) nch = strlen(apath);
   if (!nch) {
      return this;
   }

   if (funcname==0 || strlen(funcname)==0) funcname = "GetDirectory";

   TDirectory *result = this;

   char *path = new char[nch+1]; path[0] = 0;
   if (nch) strlcpy(path,apath,nch+1);
   char *s = (char*)strrchr(path, ':');
   if (s) {
      *s = '\0';
      R__LOCKGUARD(gROOTMutex);
      TDirectory *f = (TDirectory *)gROOT->GetListOfFiles()->FindObject(path);
      if (!f && !strcmp(gROOT->GetName(), path)) f = gROOT;
      if (s) *s = ':';
      if (f) {
         result = f;
         if (s && *(s+1)) result = f->GetDirectory(s+1,printError,funcname);
         delete [] path; return result;
      } else {
         if (printError) Error(funcname, "No such file %s", path);
         delete [] path; return nullptr;
      }
   }

   // path starts with a slash (assumes current file)
   if (path[0] == '/') {
      TDirectory *td = gROOT;
      result = td->GetDirectory(path+1,printError,funcname);
      delete [] path; return result;
   }

   TObject *obj;
   char *slash = (char*)strchr(path,'/');
   if (!slash) {                     // we are at the lowest level
      if (!strcmp(path, "..")) {
         result = GetMotherDir();
         delete [] path; return result;
      }
      obj = Get(path);
      if (!obj) {
         if (printError) Error(funcname,"Unknown directory %s", path);
         delete [] path; return nullptr;
      }

      //Check return object is a directory
      if (!obj->InheritsFrom(TDirectory::Class())) {
         if (printError) Error(funcname,"Object %s is not a directory", path);
         delete [] path; return nullptr;
      }
      delete [] path; return (TDirectory*)obj;
   }

   TString subdir(path);
   slash = (char*)strchr(subdir.Data(),'/');
   *slash = 0;
   //Get object with path from current directory/file
   if (!strcmp(subdir, "..")) {
      TDirectory* mom = GetMotherDir();
      if (mom)
         result = mom->GetDirectory(slash+1,printError,funcname);
      delete [] path; return result;
   }
   obj = Get(subdir);
   if (!obj) {
      if (printError) Error(funcname,"Unknown directory %s", subdir.Data());
      delete [] path; return nullptr;
   }

   //Check return object is a directory
   if (!obj->InheritsFrom(TDirectory::Class())) {
      if (printError) Error(funcname,"Object %s is not a directory", subdir.Data());
      delete [] path; return nullptr;
   }
   result = ((TDirectory*)obj)->GetDirectory(slash+1,printError,funcname);
   delete [] path; return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Change current directory to "this" directory.
///
/// Using path one can change the current directory to "path". The absolute path
/// syntax is: `file.root:/dir1/dir2`
/// where `file.root` is the file and `/dir1/dir2` the desired subdirectory
/// in the file.
///
/// Relative syntax is relative to "this" directory. E.g: `../aa`.
///
/// Returns kTRUE in case of success.

Bool_t TDirectory::cd(const char *path)
{
   return cd1(path);
}

////////////////////////////////////////////////////////////////////////////////
/// Change current directory to "this" directory.
///
/// Using path one can
/// change the current directory to "path". The absolute path syntax is:
/// `file.root:/dir1/dir2`
/// where `file.root` is the file and `/dir1/dir2` the desired subdirectory
/// in the file.
///
/// Relative syntax is relative to "this" directory. E.g: `../aa`.
///
/// Returns kFALSE in case path does not exist.

Bool_t TDirectory::cd1(const char *apath)
{
   Int_t nch = 0;
   if (apath) nch = strlen(apath);
   if (!nch) {
      auto &global = CurrentDirectory();
      RegisterGDirectory(&global);
      global = this;
      return kTRUE;
   }

   TDirectory *where = GetDirectory(apath,kTRUE,"cd");
   if (where) {
      where->cd();
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Change current directory to "path". The absolute path syntax is:
/// `file.root:/dir1/dir2`
/// where file.root is the file and `/dir1/dir2 the desired subdirectory
/// in the file.
/// Relative syntax is relative to the current directory `gDirectory`, e.g.: `../aa`.
///
/// Returns kTRUE in case of success.

Bool_t TDirectory::Cd(const char *path)
{
   return Cd1(path);
}

////////////////////////////////////////////////////////////////////////////////
/// Change current directory to "path". The path syntax is:
/// `file.root:/dir1/dir2`
/// where file.root is the file and `/dir1/dir2` the desired subdirectory
/// in the file.
///
/// Returns kFALSE in case path does not exist.

Bool_t TDirectory::Cd1(const char *apath)
{
   // null path is always true (i.e. stay in the current directory)
   Int_t nch = 0;
   if (apath) nch = strlen(apath);
   if (!nch) return kTRUE;

   TDirectory *where = gDirectory->GetDirectory(apath,kTRUE,"Cd");
   if (where) {
      where->cd();
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete all objects from a Directory list.

void TDirectory::Clear(Option_t *)
{
   if (fList) fList->Clear();

}

////////////////////////////////////////////////////////////////////////////////
/// Delete all objects from memory and directory structure itself.
/// if option is "slow", iterate through the containers in a way to can handle
///    'external' modification (induced by recursions)
/// if option is "nodelete", write the TDirectory but do not delete the contained
///    objects.
void TDirectory::Close(Option_t *option)
{
   if (!fList) {
      return;
   }

   // Save the directory key list and header
   Save();

   Bool_t nodelete = option ? (!strcmp(option, "nodelete") ? kTRUE : kFALSE) : kFALSE;

   if (!nodelete) {
      Bool_t slow = option ? (!strcmp(option, "slow") ? kTRUE : kFALSE) : kFALSE;
      if (!slow) {
         // Check if it is wise to use the fast deletion path.
         TObjLink *lnk = fList->FirstLink();
         while (lnk) {
            if (lnk->GetObject()->IsA() == TDirectory::Class()) {
               slow = kTRUE;
               break;
            }
            lnk = lnk->Next();
         }
      }

      // Delete objects from directory list, this in turn, recursively closes all
      // sub-directories (that were allocated on the heap)
      // if this dir contains subdirs, we must use the slow option for Delete!
      // we must avoid "slow" as much as possible, in particular Delete("slow")
      // with a large number of objects (eg >10^5) would take for ever.
      if (slow) fList->Delete("slow");
      else      fList->Delete();
   }

   TDirectory::CleanTargets();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete all objects from memory.

void TDirectory::DeleteAll(Option_t *)
{
   fList->Delete("slow");
}

////////////////////////////////////////////////////////////////////////////////
/// Delete Objects or/and keys in a directory.
///
///  - namecycle has the format name;cycle
///  - namecycle = "" same as namecycle ="T*"
///  - name  = * means all
///  - cycle = * means all cycles (memory and keys)
///  - cycle = "" or cycle = 9999 ==> apply to a memory object
///    When name=* use T* to delete subdirectories also
///
/// To delete one directory, you must specify the directory cycle,
/// eg.  `file.Delete("dir1;1");`
///
/// examples:
///  - foo   : delete object named foo in memory
///  - foo*  : delete all objects with a name starting with foo
///  - foo;1 : delete cycle 1 of foo on file
///  - foo;* : delete all cycles of foo on file and also from memory
///  - *;2   : delete all objects on file having the cycle 2
///  - *;*   : delete all objects from memory and file
///  - T*;*  : delete all objects from memory and file and all subdirectories

void TDirectory::Delete(const char *namecycle)
{
   if (gDebug)
     Info("Delete","Call for this = %s namecycle = %s",
               GetName(), (namecycle ? namecycle : "null"));

   TDirectory::TContext ctxt(this);
   Short_t  cycle;
   char     name[kMaxLen];
   DecodeNameCycle(namecycle, name, cycle, kMaxLen);

   Int_t deleteall    = 0;
   Int_t deletetree   = 0;
   if(strcmp(name,"*") == 0)   deleteall = 1;
   if(strcmp(name,"*T") == 0){ deleteall = 1; deletetree = 1;}
   if(strcmp(name,"T*") == 0){ deleteall = 1; deletetree = 1;}
   if(namecycle==0 || !namecycle[0]){ deleteall = 1; deletetree = 1;}
   TRegexp re(name,kTRUE);
   TString s;
   Int_t deleteOK = 0;

//*-*---------------------Case of Object in memory---------------------
//                        ========================
   if (cycle >= 9999 ) {
      TNamed *idcur;
      TIter   next(fList);
      while ((idcur = (TNamed *) next())) {
         deleteOK = 0;
         s = idcur->GetName();
         if (deleteall || s.Index(re) != kNPOS) {
            deleteOK = 1;
            if (idcur->IsA() == TDirectory::Class()) {
               deleteOK = 2;
               if (!deletetree && deleteall) deleteOK = 0;
            }
         }
         if (deleteOK != 0) {
            fList->Remove(idcur);
            if (deleteOK==2) {
               // read subdirectories to correctly delete them
               if (deletetree)
                  ((TDirectory*) idcur)->ReadAll("dirs");
               idcur->Delete(deletetree ? "T*;*" : "*");
               delete idcur;
            } else
               idcur->Delete(name);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill Graphics Structure and Paint.
///
/// Loop on all objects (memory or file) and all subdirectories

void TDirectory::Draw(Option_t *option)
{
   fList->R__FOR_EACH(TObject,Draw)(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Find object in the list of memory objects.

TObject *TDirectory::FindObject(const TObject *obj) const
{
   return fList->FindObject(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Find object by name in the list of memory objects.

TObject *TDirectory::FindObject(const char *name) const
{
   return fList->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Find object by name in the list of memory objects of the current
/// directory or its sub-directories.
/// After this call the current directory is not changed.
/// To automatically set the current directory where the object is found,
/// use FindKeyAny(aname)->ReadObj().

TObject *TDirectory::FindObjectAny(const char *aname) const
{
   //object may be already in the list of objects in memory
   TObject *obj = fList->FindObject(aname);
   if (obj) return obj;

   //try with subdirectories
   TIter next(fList);
   while( (obj = next()) ) {
      if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
         TDirectory* subdir = static_cast<TDirectory*>(obj);
         TObject *subobj = subdir->TDirectory::FindObjectAny(aname); // Explicitly recurse into _this_ exact function.
         if (subobj) {
            return subobj;
         }
      }
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to object identified by namecycle.
///
///   namecycle has the format name;cycle
///  - name  = * is illegal, cycle = * is illegal
///  - cycle = "" or cycle = 9999 ==> apply to a memory object
///
/// examples:
///  - foo   : get object named foo in memory
///                if object is not in memory, try with highest cycle from file
///  - foo;1 : get cycle 1 of foo on file
///
///  The retrieved object should in principle derive from TObject.
///  If not, the function TDirectory::GetObject should be called.
///  However, this function will still work for a non-TObject, providing that
///  the calling application cast the return type to the correct type (which
///  is the actual type of the object).
///
///  NOTE:
///
///  The method GetObject offer better protection and avoid the need
///  for any cast:
/// ~~~ {.cpp}
///      MyClass *obj;
///      directory->GetObject("some object",obj);
///      if (obj) { ... the object exist and inherits from MyClass ... }
/// ~~~
///
///  VERY IMPORTANT NOTE:
///
///  In case the class of this object derives from TObject but not
///  as a first inheritance, one must use dynamic_cast<>().
/// #### Example 1: Normal case:
/// ~~~ {.cpp}
///      class MyClass : public TObject, public AnotherClass
/// ~~~
///   then on return, one can do:
/// ~~~ {.cpp}
///      MyClass *obj = (MyClass*)directory->Get("some object of MyClass");
/// ~~~
/// #### Example 2: Special case:
/// ~~~ {.cpp}
///      class MyClass : public AnotherClass, public TObject
/// ~~~
///  then on return, one must do:
/// ~~~ {.cpp}
///      MyClass *obj = dynamic_cast<MyClass*>(directory->Get("some object of MyClass"));
/// ~~~
///  Of course, dynamic_cast<> can also be used in the example 1.

TObject *TDirectory::Get(const char *namecycle)
{
   Short_t  cycle;
   char     name[kMaxLen];

   DecodeNameCycle(namecycle, name, cycle, kMaxLen);
   char *namobj = name;
   Int_t nch = strlen(name);
   for (Int_t i = nch-1; i > 0; i--) {
      if (name[i] == '/') {
         name[i] = 0;
         TDirectory* dirToSearch=GetDirectory(name);
         namobj = name + i + 1;
         name[i] = '/';
         return dirToSearch ? dirToSearch->Get(namobj) : nullptr;
      }
   }

//*-*---------------------Case of Object in memory---------------------
//                        ========================
   TObject *idcur = fList->FindObject(namobj);
   if (idcur) {
      if (idcur==this && strlen(namobj)!=0) {
         // The object has the same name has the directory and
         // that's what we picked-up!  We just need to ignore
         // it ...
         idcur = nullptr;
      } else if (cycle == 9999) {
         return idcur;
      } else {
         if (idcur->InheritsFrom(TCollection::Class()))
            idcur->Delete();  // delete also list elements
         delete idcur;
         idcur = nullptr;
      }
   }
   return idcur;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to object identified by namecycle.
/// The returned object may or may not derive from TObject.
///
///  - namecycle has the format name;cycle
///  - name  = * is illegal, cycle = * is illegal
///  - cycle = "" or cycle = 9999 ==> apply to a memory object
///
///  VERY IMPORTANT NOTE:
///
///  The calling application must cast the returned object to
///  the final type, e.g.
/// ~~~ {.cpp}
///      MyClass *obj = (MyClass*)directory->GetObject("some object of MyClass");
/// ~~~

void *TDirectory::GetObjectUnchecked(const char *namecycle)
{
   return GetObjectChecked(namecycle,(TClass *)nullptr);
}

////////////////////////////////////////////////////////////////////////////////
/// See documentation of TDirectory::GetObjectCheck(const char *namecycle, const TClass *cl)

void *TDirectory::GetObjectChecked(const char *namecycle, const char* classname)
{
   return GetObjectChecked(namecycle, TClass::GetClass(classname));
}


////////////////////////////////////////////////////////////////////////////////
/// Return pointer to object identified by namecycle if and only if the actual
/// object is a type suitable to be stored as a pointer to a "expectedClass"
/// If expectedClass is null, no check is performed.
///
/// namecycle has the format `name;cycle`
///  - name  = * is illegal, cycle = * is illegal
///  - cycle = "" or cycle = 9999 ==> apply to a memory object
///
///  VERY IMPORTANT NOTE:
///
///  The calling application must cast the returned pointer to
///  the type described by the 2 arguments (i.e. cl):
/// ~~~ {.cpp}
///      MyClass *obj = (MyClass*)directory->GetObjectChecked("some object of MyClass","MyClass"));
/// ~~~
///  Note: We recommend using the method TDirectory::GetObject:
/// ~~~ {.cpp}
///      MyClass *obj = nullptr;
///      directory->GetObject("some object inheriting from MyClass",obj);
///      if (obj) { ... we found what we are looking for ... }
/// ~~~

void *TDirectory::GetObjectChecked(const char *namecycle, const TClass* expectedClass)
{
   Short_t  cycle;
   char     name[kMaxLen];

   DecodeNameCycle(namecycle, name, cycle, kMaxLen);
   char *namobj = name;
   Int_t nch = strlen(name);
   for (Int_t i = nch-1; i > 0; i--) {
      if (name[i] == '/') {
         name[i] = 0;
         TDirectory* dirToSearch=GetDirectory(name);
         namobj = name + i + 1;
         name[i] = '/';
         if (dirToSearch) {
            return dirToSearch->GetObjectChecked(namobj, expectedClass);
         } else {
            return nullptr;
         }
      }
   }

//*-*---------------------Case of Object in memory---------------------
//                        ========================
   if (!expectedClass || expectedClass->IsTObject()) {
      TObject *objcur = fList->FindObject(namobj);
      if (objcur) {
         if (objcur==this && strlen(namobj)!=0) {
            // The object has the same name has the directory and
            // that's what we picked-up!  We just need to ignore
            // it ...
            objcur = nullptr;
         } else if (cycle == 9999) {
            // Check type
            if (expectedClass && objcur->IsA()->GetBaseClassOffset(expectedClass) == -1) return nullptr;
            else return objcur;
         } else {
            if (objcur->InheritsFrom(TCollection::Class()))
               objcur->Delete();  // delete also list elements
            delete objcur;
            objcur = nullptr;
         }
      }
   }

   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the full path of the directory. E.g. `file:/dir1/dir2`.
/// The returned path will be re-used by the next call to GetPath().

const char *TDirectory::GetPathStatic() const
{
   static char *path = nullptr;
   const int kMAXDEPTH = 128;
   const TDirectory *d[kMAXDEPTH];
   const TDirectory *cur = this;
   int depth = 0, len = 0;

   d[depth++] = cur;
   len = strlen(cur->GetName()) + 1;  // +1 for the /

   while (cur->fMother && depth < kMAXDEPTH) {
      cur = (TDirectory *)cur->fMother;
      d[depth++] = cur;
      len += strlen(cur->GetName()) + 1;
   }

   if (path) delete [] path;
   path = new char[len+2];

   for (int i = depth-1; i >= 0; i--) {
      if (i == depth-1) {    // file or TROOT name
         strlcpy(path, d[i]->GetName(),len+2);
         strlcat(path, ":",len+2);
         if (i == 0) strlcat(path, "/",len+2);
      } else {
         strlcat(path, "/",len+2);
         strlcat(path, d[i]->GetName(),len+2);
      }
   }

   return path;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the full path of the directory. E.g. `file:/dir1/dir2`.
/// The returned path will be re-used by the next call to GetPath().

const char *TDirectory::GetPath() const
{
   FillFullPath(fPathBuffer);

   if (!GetMotherDir()) // case of file
      fPathBuffer.Append("/");

   return fPathBuffer.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Recursive method to fill full path for directory.

void TDirectory::FillFullPath(TString& buf) const
{
   TDirectory* mom = GetMotherDir();
   if (mom) {
      mom->FillFullPath(buf);
      buf += "/";
      buf += GetName();
   } else {
      buf = GetName();
      buf += ":";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a sub-directory "a" or a hierarchy of sub-directories "a/b/c/...".
///
/// Returns 0 in case of error or if a sub-directory (hierarchy) with the requested
/// name already exists.
/// returnExistingDirectory returns a pointer to an already existing sub-directory with the same name.
/// Returns a pointer to the created sub-directory or to the top sub-directory of
/// the hierarchy (in the above example, the returned TDirectory * always points
/// to "a").
/// In particular, the steps to create first a/b/c and then a/b/d without receiving
/// errors are:
/// ~~~ {.cpp}
///    TFile * file = new TFile("afile","RECREATE");
///    file->mkdir("a");
///    file->cd("a");
///    gDirectory->mkdir("b/c");
///    gDirectory->cd("b");
///    gDirectory->mkdir("d");
/// ~~~
/// or
/// ~~~ {.cpp}
///    TFile * file = new TFile("afile","RECREATE");
///    file->mkdir("a");
///    file->cd("a");
///    gDirectory->mkdir("b/c");
///    gDirectory->mkdir("b/d", "", true);
/// ~~~

TDirectory *TDirectory::mkdir(const char *name, const char *title, Bool_t returnExistingDirectory)
{
   if (returnExistingDirectory) {
      auto existingdir = GetDirectory(name);
      if (existingdir)
        return existingdir;
   }
   if (!name || !title || !name[0]) return nullptr;
   if (!title[0]) title = name;
   TDirectory *newdir = nullptr;
   if (const char *slash = strchr(name,'/')) {
      Long_t size = Long_t(slash-name);
      char *workname = new char[size+1];
      strncpy(workname, name, size);
      workname[size] = 0;
      TDirectory *tmpdir;
      GetObject(workname,tmpdir);
      if (!tmpdir)
         tmpdir = mkdir(workname,title);
      delete[] workname;
      if (!tmpdir) return nullptr;
      if (!newdir) newdir = tmpdir;
      return tmpdir->mkdir(slash+1);
   }

   TDirectory::TContext ctxt(this);

   newdir = new TDirectory(name, title, "", this);

   return newdir;
}

////////////////////////////////////////////////////////////////////////////////
/// List Directory contents.
///
/// Indentation is used to identify the directory tree
/// Subdirectories are listed first, then objects in memory.
///
/// The option can has the following format:
///
///      [<regexp>]
///
/// The `<regexp>` will be used to match the name of the objects.
/// By default memory and disk objects are listed.

void TDirectory::ls(Option_t *option) const
{
   TROOT::IndentLevel();
   TROOT::IncreaseDirLevel();

   TString opta = option;
   TString opt  = opta.Strip(TString::kBoth);
   Bool_t memobj  = kTRUE;
   TString reg = "*";
   if (opt.BeginsWith("-m")) {
      if (opt.Length() > 2)
         reg = opt(2,opt.Length());
   } else if (opt.BeginsWith("-d")) {
      memobj  = kFALSE;
      if (opt.Length() > 2)
         reg = opt(2,opt.Length());
   } else if (!opt.IsNull())
      reg = opt;

   TRegexp re(reg, kTRUE);

   if (memobj) {
      TObject *obj;
      TIter nextobj(fList);
      while ((obj = (TObject *) nextobj())) {
         TString s = obj->GetName();
         if (s.Index(re) == kNPOS) continue;
         obj->ls(option);            //*-* Loop on all the objects in memory
      }
   }
   TROOT::DecreaseDirLevel();
}

////////////////////////////////////////////////////////////////////////////////
/// Paint all objects in the directory.

void TDirectory::Paint(Option_t *option)
{
   fList->R__FOR_EACH(TObject,Paint)(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Print all objects in the directory.

void TDirectory::Print(Option_t *option) const
{
   fList->R__FOR_EACH(TObject,Print)(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Print the path of the directory.

void TDirectory::pwd() const
{
   Printf("%s", GetPath());
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively remove object from a Directory.

void TDirectory::RecursiveRemove(TObject *obj)
{
   fList->RecursiveRemove(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove an object from the in-memory list.

TObject *TDirectory::Remove(TObject* obj)
{
   TObject *p = 0;
   if (fList) {
      p = fList->Remove(obj);
   }
   return p;
}

////////////////////////////////////////////////////////////////////////////////
/// Removes subdirectory from the directory
/// When directory is deleted, all keys in all subdirectories will be
/// read first and deleted from file (if exists)
/// Equivalent call is Delete("name;*");

void TDirectory::rmdir(const char *name)
{
   if ((name==0) || (*name==0)) return;

   TString mask(name);
   mask+=";*";
   Delete(mask);
}

////////////////////////////////////////////////////////////////////////////////
/// Save object in filename,
/// if filename is 0 or "", a file with "objectname.root" is created.
/// The name of the key is the object name.
/// If the operation is successful, it returns the number of bytes written to the file
/// otherwise it returns 0.
/// By default a message is printed. Use option "q" to not print the message.
/// If filename contains ".json" extension, JSON representation of the object
/// will be created and saved in the text file. Such file can be used in
/// JavaScript ROOT (https://root.cern.ch/js/) to display object in web browser
/// When creating JSON file, option string may contain compression level from 0 to 3 (default 0)

Int_t TDirectory::SaveObjectAs(const TObject *obj, const char *filename, Option_t *option) const
{
   if (!obj) return 0;
   Int_t nbytes = 0;
   TString fname = filename;
   if (!filename || !filename[0]) {
      fname.Form("%s.root",obj->GetName());
   }
   TString cmd;
   if (fname.Index(".json") > 0) {
      cmd.Form("TBufferJSON::ExportToFile(\"%s\",(TObject*) %s, \"%s\");", fname.Data(), TString::LLtoa((Long_t)obj, 10).Data(), (option ? option : ""));
      nbytes = gROOT->ProcessLine(cmd);
   } else {
      cmd.Form("TFile::Open(\"%s\",\"recreate\");",fname.Data());
      TContext ctxt; // The TFile::Open will change the current directory.
      TDirectory *local = (TDirectory*)gROOT->ProcessLine(cmd);
      if (!local) return 0;
      nbytes = obj->Write();
      delete local;
   }
   TString opt(option);
   opt.ToLower();
   if (!opt.Contains("q")) {
      if (!gSystem->AccessPathName(fname.Data())) obj->Info("SaveAs", "ROOT file %s has been created", fname.Data());
   }
   return nbytes;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the name for directory
/// If the directory name is changed after the directory was written once,
/// ROOT currently would NOT change the name of correspondent key in the
/// mother directory.
/// DO NOT use this method to 'rename a directory'.
/// Renaming a directory is currently NOT supported.

void TDirectory::SetName(const char* newname)
{
   TNamed::SetName(newname);
}

////////////////////////////////////////////////////////////////////////////////
/// Encode the name and cycle into buffer like: "aap;2".

void TDirectory::EncodeNameCycle(char *buffer, const char *name, Short_t cycle)
{
   if (cycle == 9999)
      strcpy(buffer, name);
   else
      sprintf(buffer, "%s;%d", name, cycle);
}

////////////////////////////////////////////////////////////////////////////////
/// Decode a namecycle "aap;2" into name "aap" and cycle "2". Destination
/// buffer size for name (including string terminator) should be specified in
/// namesize.

void TDirectory::DecodeNameCycle(const char *buffer, char *name, Short_t &cycle,
                                 const size_t namesize)
{
   size_t len = 0;
   const char *ni = strchr(buffer, ';');

   if (ni) {
      // Found ';'
      len = ni - buffer;
      ++ni;
   } else {
      // No ';' found
      len = strlen(buffer);
      ni = &buffer[len];
   }

   if (namesize) {
      if (len > namesize-1ul) len = namesize-1;  // accommodate string terminator
   } else {
      ::Warning("TDirectory::DecodeNameCycle",
         "Using unsafe version: invoke this metod by specifying the buffer size");
   }

   strncpy(name, buffer, len);
   name[len] = '\0';

   if (*ni == '*')
      cycle = 10000;
   else if (isdigit(*ni)) {
      long parsed = strtol(ni,nullptr,10);
      if (parsed >= (long) std::numeric_limits<Short_t>::max())
         cycle = 0;
      else
         cycle = (Short_t)parsed;
   } else
      cycle = 9999;
}

///////////////////////////////////////////////////////////////////////////////
/// Register a TContext pointing to this TDirectory object

void TDirectory::RegisterContext(TContext *ctxt) {
   ROOT::Internal::TSpinLockGuard slg(fSpinLock);

   if (fContext) {
      TContext *current = fContext;
      while(current->fNext) {
         current = current->fNext;
      }
      current->fNext = ctxt;
      ctxt->fPrevious = current;
   } else {
      fContext = ctxt;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Register a std::atomic<TDirectory*> pointing to this TDirectory object

void TDirectory::RegisterGDirectory(std::atomic<TDirectory*> *globalptr)
{
   ROOT::Internal::TSpinLockGuard slg(fSpinLock);

   auto oldvalue = globalptr->load();

   auto iter = std::find(fGDirectories.begin(), fGDirectories.end(), globalptr);
   if (iter != fGDirectories.begin())
      fGDirectories.push_back(globalptr);

   if (oldvalue && oldvalue != this) {
      iter = std::find(oldvalue->fGDirectories.begin(), oldvalue->fGDirectories.end(), globalptr);
      if (iter != oldvalue->fGDirectories.begin())
         oldvalue->fGDirectories.erase(iter);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// \copydoc TDirectory::WriteTObject().

Int_t TDirectory::WriteTObject(const TObject *obj, const char *name, Option_t * /*option*/, Int_t /*bufsize*/)
{
   const char *objname = "no name specified";
   if (name) objname = name;
   else if (obj) objname = obj->GetName();
   Error("WriteTObject","The current directory (%s) is not associated with a file. The object (%s) has not been written.",GetName(),objname);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// UnRegister a TContext pointing to this TDirectory object

void TDirectory::UnregisterContext(TContext *ctxt) {

   ROOT::Internal::TSpinLockGuard slg(fSpinLock);

   // Another thread already unregistered the TContext.
   if (ctxt->fDirectory == nullptr)
      return;

   if (ctxt==fContext) {
      fContext = ctxt->fNext;
      if (fContext) fContext->fPrevious = nullptr;
      ctxt->fPrevious = ctxt->fNext = nullptr;
   } else {
      TContext *next = ctxt->fNext;
      ctxt->fPrevious->fNext = next;
      if (next) next->fPrevious = ctxt->fPrevious;
      ctxt->fPrevious = ctxt->fNext = nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TDirectory Streamer.
void TDirectory::Streamer(TBuffer &R__b)
{
   // Stream an object of class TDirectory.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> fMother;
      R__b >> fList;
      fList->UseRWLock();
      fUUID.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TDirectory::IsA());
   } else {
      R__c = R__b.WriteVersion(TDirectory::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << fMother;
      R__b << fList;
      fUUID.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}
