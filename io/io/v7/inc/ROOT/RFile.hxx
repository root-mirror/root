/// \file ROOT/RFile.hxx
/// \ingroup Base ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-07-31
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RFile
#define ROOT7_RFile

#include "ROOT/RDirectory.hxx"
#include "ROOT/RMakeUnique.hxx"
#include "ROOT/RStringView.hxx"

#include "TClass.h"
#include <memory>

namespace ROOT {
namespace Experimental {

class RFilePtr;

namespace Internal {
class RFileStorageInterface;
class RFileSharedPtrCtor;
} // namespace Internal

/** \class ROOT::Experimental::RFile
 A ROOT file.

 A ROOT file is an object store: it can serialize any
 object for which ROOT I/O is available (generally: an object which has a
 dictionary), and it stores the object's data under a key name.

 */
class RFile: public RDirectory {
private:
   std::unique_ptr<Internal::RFileStorageInterface> fStorage; ///< Storage backend.

   RFile(std::unique_ptr<Internal::RFileStorageInterface> &&storage);

   /// Serialize the object at address, using the object's TClass.
   // FIXME: what about `cl` "pointing" to a base class?
   void WriteMemoryWithType(std::string_view name, const void *address, TClass *cl);

   friend Internal::RFileSharedPtrCtor;

public:
   /// Options for RFile construction.
   struct Options_t {
      /// Default constructor needed for member inits.
      Options_t() {}

      /// Whether the file should be opened asynchronously, if available.
      bool fAsynchronousOpen = false;

      /// Timeout for asynchronous opening.
      int fAsyncTimeout = 0;

      /// Whether the file should be cached before reading. Only available for
      /// "remote" file protocols. If the download fails, the file will be opened
      /// remotely.
      bool fCachedRead = false;

      /// Where to cache the file. If empty, defaults to RFilePtr::GetCacheDir().
      std::string fCacheDir;
   };

   ///\name Generator functions
   ///\{

   /// Open a file with `name` for reading.
   ///
   /// \note: Synchronizes multi-threaded accesses through locks.
   static RFilePtr Open(std::string_view name, const Options_t &opts = Options_t());

   /// Open an existing file with `name` for reading and writing. If a file with
   /// that name does not exist, an invalid RFilePtr will be returned.
   ///
   /// \note: Synchronizes multi-threaded accesses through locks.
   static RFilePtr OpenForUpdate(std::string_view name, const Options_t &opts = Options_t());

   /// Open a file with `name` for reading and writing. Fail (return an invalid
   /// `RFilePtr`) if a file with this name already exists.
   ///
   /// \note: Synchronizes multi-threaded accesses through locks.
   static RFilePtr Create(std::string_view name, const Options_t &opts = Options_t());

   /// Open a file with `name` for reading and writing. If a file with this name
   /// already exists, delete it and create a new one. Else simply create a new file.
   ///
   /// \note: Synchronizes multi-threaded accesses through locks.
   static RFilePtr Recreate(std::string_view name, const Options_t &opts = Options_t());

   ///\}

   /// Set the new directory used for cached reads, returns the old directory.
   ///
   /// \note: Synchronizes multi-threaded accesses through locks.
   static std::string SetCacheDir(std::string_view path);

   /// Get the directory used for cached reads.
   static std::string GetCacheDir();

   /// Must not call Write() of all attached objects:
   /// some might not be needed to be written or writing might be aborted due to
   /// an exception; require explicit Write().
   ~RFile();

   /// Save all objects associated with this directory (including file header) to
   /// the storage medium.
   void Flush();

   /// Flush() and make the file non-writable: close it.
   void Close();

   /// Read the object for a key. `T` must be the object's type.
   /// This will re-read the object for each call, returning a new copy; whether
   /// the `RDirectory` is managing an object attached to this key or not.
   /// \returns a `unique_ptr` to the object.
   /// \throws RDirectoryUnknownKey if no object is stored under this name.
   /// \throws RDirectoryTypeMismatch if the object stored under this name is of
   ///   a type different from `T`.
   template <class T>
   std::unique_ptr<T> Read(std::string_view name)
   {
      // FIXME: need separate collections for a RDirectory's key/value and registered objects. Here, we want to emit a
      // read and must look through the key/values without attaching an object to the RDirectory.
      // FIXME: do not register read object in RDirectory
      // FIXME: implement actual read
      // FIXME: for now, copy out of whatever the RDirectory manages.
      return std::make_unique<T>(*Get<T>(name));
   }

   /// Write an object that is not lifetime managed by this RFileImplBase.
   template <class T>
   void Write(std::string_view name, const T &obj)
   {
      WriteMemoryWithType(name, &obj, TClass::GetClass<T>());
   }

   /// Write an object that is not lifetime managed by this RFileImplBase.
   template <class T>
   void Write(std::string_view name, const T *obj)
   {
      WriteMemoryWithType(name, obj, TClass::GetClass<T>());
   }

   /// Write an object that is already lifetime managed by this RFileImplBase.
   void Write(std::string_view name)
   {
      auto dep = Find(name);
      WriteMemoryWithType(name, dep.GetPointer().get(), dep.GetType());
   }

   /// Hand over lifetime management of an object to this RFileImplBase, and
   /// write it.
   template <class T>
   void Write(std::string_view name, std::shared_ptr<T> &&obj)
   {
      Add(name, obj);
      // FIXME: use an iterator from the insertion to write instead of a second name lookup.
      Write(name);
   }
};

/**
 \class RFilePtr
 \brief Points to an object that stores or reads objects in ROOT's binary
 format.

 FIXME: implement async open; likely using std::future, possibly removing the
 Option_t element.

 */

class RFilePtr {
private:
   std::shared_ptr<RFile> fFile;

   /// Constructed by Open etc.
   RFilePtr(std::shared_ptr<RFile> &&);

   friend class RFile;

public:
   /// Dereference the file pointer, giving access to the RFileImplBase object.
   RFile *operator->() { return fFile.get(); }

   /// Dereference the file pointer, giving access to the RFileImplBase object.
   /// const overload.
   const RFile *operator->() const { return fFile.get(); }

   /// Check the validity of the file pointer.
   operator bool() const { return fFile.get(); }
};

} // namespace Experimental
} // namespace ROOT
#endif
