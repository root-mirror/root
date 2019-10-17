/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RBrowsable
#define ROOT7_RBrowsable


#include <ROOT/RBrowserItem.hxx>

#include <memory>
#include <string>
#include <map>
#include <vector>


class TClass;
class TObject;

namespace ROOT {
namespace Experimental {

namespace Browsable {

class RLevelIter;

/** \class RElement
\ingroup rbrowser
\brief Basic element of RBrowsable hierarchy. Provides access to data, creates iterator if any
\author Sergey Linev <S.Linev@gsi.de>
\date 2019-10-14
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RElement {
public:
   virtual ~RElement() = default;

   /** Name of RBrowsable, must be provided in derived classes */
   virtual std::string GetName() const = 0;

   /** Title of RBrowsable (optional) */
   virtual std::string GetTitle() const { return ""; }

   /** Create iterator for childs elements if any */
   virtual std::unique_ptr<RLevelIter> GetChildsIter() { return nullptr; }

   virtual bool HasTextContent() const { return false; }

   virtual std::string GetTextContent() { return ""; }

   /** Temporary solution, later better interface should be provided */
   virtual TObject *GetObjectToDraw() { return nullptr; }
};

/** \class RLevelIter
\ingroup rbrowser
\brief Iterator over single level hierarchy like TList
\author Sergey Linev <S.Linev@gsi.de>
\date 2019-10-14
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/


class RLevelIter {
public:
   virtual ~RLevelIter() = default;

   /** Shift to next element */
   virtual bool Next() { return false; }

   /** Reset iterator to the first element, returns false if not supported */
   virtual bool Reset() { return false; }

   /** Is there current element  */
   virtual bool HasItem() const { return false; }

   /** Returns current element name  */
   virtual std::string GetName() const { return ""; }

   virtual bool Find(const std::string &name);

   /** If element may have childs: 0 - no, >0 - yes, -1 - maybe */
   virtual int CanHaveChilds() const { return 0; }

   virtual std::unique_ptr<RBrowserItem> CreateBrowserItem()
   {
      return std::make_unique<RBrowserItem>(GetName(), CanHaveChilds());
   }

   /** Returns full information for current element */
   virtual std::shared_ptr<RElement> GetElement() { return nullptr; }
};


/** \class RProvider
\ingroup rbrowser
\brief Provider of different browsing methods for supported classes
\author Sergey Linev <S.Linev@gsi.de>
\date 2019-10-14
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/


class RProvider {

   using Map_t = std::map<const TClass*, std::shared_ptr<RProvider>>;

   static Map_t &GetMap();

public:
   virtual ~RProvider() = default;

   /** Returns supported class */
   virtual const TClass *GetSupportedClass() const = 0;

   /** Returns true if derived classes supported as well */
   virtual bool SupportDerivedClasses() const { return false; }


   static void Register(std::shared_ptr<RProvider> provider);
   static std::shared_ptr<RProvider> GetProvider(const TClass *cl, bool check_base = true);

};

} // namespace Browsable


/** \class RBrowsable
\ingroup rbrowser
\brief Way to browse (hopefully) everything in ROOT
\author Sergey Linev <S.Linev@gsi.de>
\date 2019-10-14
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/


class RBrowsable {

   struct RLevel {
      std::string name;
      std::unique_ptr<Browsable::RLevelIter> iter;
      std::shared_ptr<Browsable::RElement> item;
      RLevel(const std::string &_name) : name(_name) {}
   };

   std::shared_ptr<Browsable::RElement> fItem; ///<! top-level item to browse
   std::vector<RLevel> fLevels;           ///<! navigated levels

   bool Navigate(const std::vector<std::string> &path);

   bool DecomposePath(const std::string &path, std::vector<std::string> &arr);

public:
   RBrowsable() = default;

   RBrowsable(std::shared_ptr<Browsable::RElement> item)
   {
      fItem = item;
   }

   virtual ~RBrowsable() = default;


   void SetTopItem(std::shared_ptr<Browsable::RElement> item)
   {
      fLevels.clear();
      fItem = item;
   }

   bool ProcessRequest(const RBrowserRequest &request, RBrowserReplyNew &reply);

   std::shared_ptr<Browsable::RElement> GetElement(const std::string &path);
};


} // namespace Experimental
} // namespace ROOT

#endif
