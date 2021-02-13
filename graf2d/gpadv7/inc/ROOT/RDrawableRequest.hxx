/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RDrawableRequest
#define ROOT7_RDrawableRequest

#include <string>
#include <cstdint>
#include <memory>

#include <ROOT/RDrawable.hxx>

namespace ROOT {
namespace Experimental {


/** \class RDrawableReply
\ingroup GpadROOT7
\brief Base class for replies on RDrawableRequest
\author Sergey Linev <s.linev@gsi.de>
\date 2020-04-14
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RDrawableReply {
   uint64_t reqid{0}; ///< request id

public:

   void SetRequestId(uint64_t _reqid) { reqid = _reqid; }
   uint64_t GetRequestId() const { return reqid; }

   virtual ~RDrawableReply();
};


/** \class RDrawableRequest
\ingroup GpadROOT7
\brief Base class for requests which can be submitted from the clients
\author Sergey Linev <s.linev@gsi.de>
\date 2020-04-14
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RDrawableRequest {
   std::string id; ///< drawable id
   uint64_t reqid{0}; ///< request id

   RDrawable::RDisplayContext fContext; ///<! display context

   //const RCanvas *fCanvas{nullptr}; ///<! pointer on canvas, can be used in Process
   //const RPadBase *fPad{nullptr};   ///<! pointer on pad with drawable, can be used in Process
   //RDrawable *fDrawable{nullptr};   ///<! pointer on drawable, can be used in Process

public:
   const std::string &GetId() const { return id; }
   uint64_t GetRequestId() const { return reqid; }

   const RDrawable::RDisplayContext &GetContext() const { return fContext; }
   RDrawable::RDisplayContext &GetContext() { return fContext; }

   //void SetCanvas(const RCanvas *canv) { fCanvas = canv; }
   //void SetPad(const RPadBase *pad) { fPad = pad; }
   //void SetDrawable(RDrawable *dr) { fDrawable = dr; }

   virtual ~RDrawableRequest();

   bool ShouldBeReplyed() const { return GetRequestId() > 0; }

   virtual std::unique_ptr<RDrawableReply> Process() { return nullptr; }

   virtual bool NeedCanvasUpdate() const { return false; }
};


/** \class RDrawableExecRequest
\ingroup GpadROOT7
\brief Request execution of method of referenced drawable, no reply
\author Sergey Linev <s.linev@gsi.de>
\date 2020-04-14
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RDrawableExecRequest : public RDrawableRequest {
   std::string exec; ///< that to execute
public:
   std::unique_ptr<RDrawableReply> Process() override;
};

} // namespace Experimental
} // namespace ROOT


#endif
