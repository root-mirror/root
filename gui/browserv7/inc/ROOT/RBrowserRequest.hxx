// Authors: Bertrand Bellenot <bertrand.bellenot@cern.ch> Sergey Linev <S.Linev@gsi.de>
// Date: 2019-02-28
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RBrowserRequest
#define ROOT7_RBrowserRequest

#include <string>
#include <vector>

namespace ROOT {
namespace Experimental {

/** \class ROOT::Experimental::RBrowserRequest
\ingroup rbrowser
\brief Request send from client to get content of path element
*/

class RBrowserRequest {
public:
   std::vector<std::string> path; ///< requested path
   int first{0};          ///< first child to request
   int number{0};         ///< number of childs to request, 0 - all childs
   std::string sort;      ///< kind of sorting
   bool reverse{false};   ///< reverse item order
   bool hidden{false};    ///< show hidden files
   bool reload{false};    ///< force items reload
   std::string regex;     ///< applied regex
};

} // namespace Experimental
} // namespace ROOT

#endif


