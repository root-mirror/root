// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2019-10-14
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RBrowserData
#define ROOT7_RBrowserData

#include <ROOT/Browsable/RElement.hxx>

#include <ROOT/RBrowserRequest.hxx>
#include <ROOT/RBrowserReply.hxx>

#include <memory>
#include <string>
#include <vector>
#include <map>

namespace ROOT {
namespace Experimental {

class RLogChannel;
/// Log channel for Browser diagnostics.
RLogChannel &BrowserLog();


class RBrowserData {

   std::shared_ptr<Browsable::RElement> fTopElement;    ///<! top element

   Browsable::RElementPath_t  fWorkingPath;             ///<! path showed in Breadcrumb

   std::map<Browsable::RElementPath_t, std::shared_ptr<Browsable::RElement>> fCache; ///<! already requested elements

   Browsable::RElementPath_t fLastPath;                  ///<! path to last used element
   std::shared_ptr<Browsable::RElement> fLastElement;    ///<! last element used in request
   std::vector<std::unique_ptr<Browsable::RItem>> fLastItems; ///<! created browser items - used in requests
   bool fLastAllChilds{false};                           ///<! if all chlds were extracted
   std::vector<const Browsable::RItem *> fLastSortedItems;   ///<! sorted child items, used in requests
   std::string fLastSortMethod;                          ///<! last sort method
   bool fLastSortReverse{false};                         ///<! last request reverse order

   void ResetLastRequestData(bool with_element);

   bool ProcessBrowserRequest(const RBrowserRequest &request, RBrowserReply &reply);

public:
   RBrowserData() = default;

   RBrowserData(std::shared_ptr<Browsable::RElement> elem) { SetTopElement(elem); }

   virtual ~RBrowserData() = default;

   void SetTopElement(std::shared_ptr<Browsable::RElement> elem);

   void SetWorkingPath(const Browsable::RElementPath_t &path);

   const Browsable::RElementPath_t &GetWorkingPath() const { return fWorkingPath; }

   std::string ProcessRequest(const RBrowserRequest &request);

   std::shared_ptr<Browsable::RElement> GetElement(const std::string &str);
   std::shared_ptr<Browsable::RElement> GetElementFromTop(const Browsable::RElementPath_t &path);

   Browsable::RElementPath_t DecomposePath(const std::string &path, bool relative_to_work_element);
   std::shared_ptr<Browsable::RElement> GetSubElement(const Browsable::RElementPath_t &path);


};


} // namespace Experimental
} // namespace ROOT

#endif
