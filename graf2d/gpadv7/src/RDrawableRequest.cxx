/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RDrawableRequest.hxx"

#include "ROOT/RDrawable.hxx"
#include "ROOT/RCanvas.hxx"

using namespace ROOT::Experimental;


/////////////////////////////////////////////////////////////////////////////
/// destructor, pin vtable

RDrawableReply::~RDrawableReply() = default;

/////////////////////////////////////////////////////////////////////////////
/// destructor, pin vtable

RDrawableRequest::~RDrawableRequest() = default;


/////////////////////////////////////////////////////////////////////////////
/// Execute method of the drawable

std::unique_ptr<RDrawableReply> RDrawableExecRequest::Process()
{
   if (!exec.empty() && GetContext().GetDrawable()) {
      std::string buf = exec;

      // many operations can be separated by ";;" string
      // TODO: exclude potential mistake if such symbols appears inside quotes

      while (!buf.empty()) {
         std::string sub = buf;
         auto pos = buf.find(";;");
         if (pos == std::string::npos) {
            sub = buf;
            buf.clear();
         } else {
            sub = buf.substr(0,pos);
            buf = buf.substr(pos+2);
         }
         if (!sub.empty())
            GetContext().GetDrawable()->Execute(sub);
      }
   }

   if (GetContext().GetCanvas())
      GetContext().GetCanvas()->Modified();

   return nullptr;
}
