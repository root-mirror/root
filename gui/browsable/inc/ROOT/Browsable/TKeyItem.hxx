/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_Browsable_TKeyItem
#define ROOT7_Browsable_TKeyItem

#include <ROOT/Browsable/RItem.hxx>

class TDirectory;

namespace ROOT {
namespace Experimental {
namespace Browsable {

/** \class TKeyItem
\ingroup rbrowser
Representation of single item in the file browser for object from TKey
*/

class TKeyItem : public RItem {
   std::string className; ///< class name

public:

   TKeyItem() = default;

   TKeyItem(const std::string &_name, int _nchilds) : RItem(_name, _nchilds) {}

   // should be here, one needs virtual table for correct streaming of RRootBrowserReply
   virtual ~TKeyItem() = default;

   void SetClassName(const std::string &_className) { className = _className; }
};

} // namespace Browsable
} // namespace Experimental
} // namespace ROOT


#endif
