/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RDrawable
#define ROOT7_RDrawable

#include <memory>
#include <string>
#include <vector>

#include <ROOT/RAttrMap.hxx>
#include <ROOT/RStyle.hxx>

namespace ROOT {
namespace Experimental {

class RMenuItems;
class RPadBase;
class RAttrBase;
class RDisplayItem;
class RDrawableDisplayItem;
class RIndirectDisplayItem;
class RLegend;
class RCanvas;
class RChangeAttrRequest;
class RDrawableMenuRequest;
class RDrawableExecRequest;

namespace Internal {

/** \class RIOSharedBase
\ingroup GpadROOT7
\author Sergey Linev <s.linev@gsi.de>
\date 2019-09-24
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RIOSharedBase {
public:
   virtual const void *GetIOPtr() const = 0;
   virtual bool HasShared() const = 0;
   virtual void *MakeShared() = 0;
   virtual void SetShared(void *shared) = 0;
   virtual ~RIOSharedBase() = default;
};

using RIOSharedVector_t = std::vector<RIOSharedBase *>;

template<class T>
class RIOShared final : public RIOSharedBase {
   std::shared_ptr<T>  fShared;  ///<!   holder of object
   T* fIO{nullptr};              ///<    plain pointer for IO
public:
   const void *GetIOPtr() const final { return fIO; }
   bool HasShared() const final { return fShared.get() != nullptr; }
   void *MakeShared() final { fShared.reset(fIO); return &fShared; }
   void SetShared(void *shared) final { fShared = *((std::shared_ptr<T> *) shared); }

   RIOShared() = default;

   RIOShared(const std::shared_ptr<T> &ptr) : RIOSharedBase()
   {
      fShared = ptr;
      fIO = ptr.get();
   }

   RIOShared &operator=(const std::shared_ptr<T> &ptr)
   {
      fShared = ptr;
      fIO = ptr.get();
      return *this;
   }

   operator bool() const { return !!fShared || !!fIO; }

   const T *get() const { return fShared ? fShared.get() : fIO; }
   T *get() { return fShared ? fShared.get() : fIO; }

   const T *operator->() const { return get(); }
   T *operator->() { return get(); }

   std::shared_ptr<T> get_shared() const { return fShared; }

   void reset() { fShared.reset(); fIO = nullptr; }

   // reset IO pointer, object will not be stored in normal output operation
   void reset_io() { fIO = nullptr; }

   // restore IO pointer, object will be stored in normal output operation
   void restore_io() { fIO = fShared.get(); }
};

} // namespace Internal

/** \class RDrawable
\ingroup GpadROOT7
\brief Base class for drawable entities: objects that can be painted on a `RPad`.
\author Axel Naumann <axel@cern.ch>
\author Sergey Linev <s.linev@gsi.de>
\date 2015-08-07
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RDrawable {

friend class RPadBase; // to access Display method and IsFrameRequired
friend class RAttrBase;
friend class RStyle;
friend class RLegend; // to access CollectShared method
friend class RDrawableDisplayItem;  // to call OnDisplayItemDestroyed
friend class RIndirectDisplayItem;  // to access attributes and other members
friend class RChangeAttrRequest; // access SetDrawableVersion and AttrMap
friend class RDrawableMenuRequest; // access PopulateMenu method
friend class RDrawableExecRequest; // access Execute() method

public:

   using Version_t = uint64_t;

   class RDisplayContext {
      RCanvas *fCanvas{nullptr};     ///<! canvas where drawable is displayed
      RPadBase *fPad{nullptr};       ///<! subpad where drawable is displayed
      RDrawable *fDrawable{nullptr}; ///<! reference on the drawable
      Version_t fLastVersion{0};     ///<! last displayed version
      unsigned fIndex{0};            ///<! index in list of primitives
      unsigned fConnId{0};           ///<! connection id
      bool fMainConn{false};         ///<! is main connection

   public:

      RDisplayContext() = default;

      RDisplayContext(RCanvas *canv, RPadBase *pad, Version_t vers = 0) :
         fCanvas(canv), fPad(pad), fLastVersion(vers)
      {
      }

      /** Set canvas */
      void SetCanvas(RCanvas *canv) { fCanvas = canv; }
      /** Set pad */
      void SetPad(RPadBase *pad) { fPad = pad; }
      /** Set drawable and its index in list of primitives */
      void SetDrawable(RDrawable *dr, unsigned indx)
      {
         fDrawable = dr;
         fIndex = indx;
      }
      /** Set connection id and ismain flag for connection */
      void SetConnection(unsigned connid, bool ismain)
      {
         fConnId = connid;
         fMainConn = ismain;
      }

      RCanvas *GetCanvas() const { return fCanvas; }
      RPadBase *GetPad() const { return fPad; }
      RDrawable *GetDrawable() const { return fDrawable; }
      unsigned GetIndex() const { return fIndex; }

      Version_t GetLastVersion() const { return fLastVersion; }

      unsigned GetConnId() const { return fConnId; }
      bool IsMainConn() const { return fMainConn; }
   };

private:
   RAttrMap fAttr;               ///< attributes values
   std::weak_ptr<RStyle> fStyle; ///<! style applied for RDrawable, not stored when canvas is saved
   std::string fCssType;         ///<! drawable type, not stored in the root file, must be initialized in constructor
   std::string fCssClass;        ///< user defined drawable class, can later go inside map
   std::string fId;              ///< optional object identifier, may be used in CSS as well
   Version_t fVersion{1};        ///<! drawable version, changed from the canvas

protected:

   virtual void CollectShared(Internal::RIOSharedVector_t &) {}

   virtual bool IsFrameRequired() const { return false; }

   RAttrMap &GetAttrMap() { return fAttr; }
   const RAttrMap &GetAttrMap() const { return fAttr; }

   bool MatchSelector(const std::string &selector) const;

   virtual std::unique_ptr<RDisplayItem> Display(const RDisplayContext &);

   virtual void OnDisplayItemDestroyed(RDisplayItem *) const {}

   virtual void SetDrawableVersion(Version_t vers) { fVersion = vers; }
   Version_t GetVersion() const { return fVersion; }

   virtual void PopulateMenu(RMenuItems &);

   virtual void Execute(const std::string &);

   RDrawable(const RDrawable &) = delete;
   RDrawable &operator=(const RDrawable &) = delete;

public:

   explicit RDrawable(const std::string &type) : fCssType(type) {}

   virtual ~RDrawable();

   virtual void UseStyle(const std::shared_ptr<RStyle> &style) { fStyle = style; }
   void ClearStyle() { UseStyle(nullptr); }

   void SetCssClass(const std::string &cl) { fCssClass = cl; }
   const std::string &GetCssClass() const { return fCssClass; }

   const std::string &GetCssType() const { return fCssType; }

   const std::string &GetId() const { return fId; }
   void SetId(const std::string &id) { fId = id; }
};

/// Central method to insert drawable in list of pad primitives
/// By default drawable placed as is.

template <class DRAWABLE, std::enable_if_t<std::is_base_of<RDrawable, DRAWABLE>{}>* = nullptr>
inline auto GetDrawable(const std::shared_ptr<DRAWABLE> &drawable)
{
   return drawable;
}

} // namespace Experimental
} // namespace ROOT

#endif
