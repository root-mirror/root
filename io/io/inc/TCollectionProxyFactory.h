// @(#)root/io:$Id$
// Author: Markus Frank  28/10/04

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TCollectionProxyFactory
#define ROOT_TCollectionProxyFactory

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  Small helper to save proxy environment in the event of
//  recursive calls.
//
//////////////////////////////////////////////////////////////////////////

#include "TCollectionProxyInfo.h"

#include "TClassStreamer.h"

#include "TMemberStreamer.h"

#include "TGenCollectionProxy.h"

// Forward declarations
class TBuffer;
class TGenCollectionProxy;
class TGenCollectionStreamer;
class TVirtualCollectionProxy;
class TEmulatedCollectionProxy;

#if defined(_WIN32)
   #if _MSC_VER<1300
      #define TYPENAME
      #define R__VCXX6
   #else
      #define TYPENAME typename
   #endif
#else
   #define TYPENAME typename
#endif


/** \class TCollectionProxyFactory TCollectionProxyFactory.h
 TCollectionProxyFactory
 Interface to collection proxy and streamer generator.
 Proxy around an arbitrary container, which implements basic
 functionality and iteration. The purpose of this implementation
 is to shield any generated dictionary implementation from the
 underlying streamer/proxy implementation and only expose
 the creation functions.

 In particular this is used to implement splitting and abstract
 element access of any container. Access to compiled code is necessary
 to implement the abstract iteration sequence and functionality like
 size(), clear(), resize(). resize() may be a void operation.

 \author  M.Frank
 \version 1.0
*/
class TCollectionProxyFactory  {
public:

   typedef TVirtualCollectionProxy Proxy_t;
#ifdef R__HPUX
   typedef const std::type_info&      Info_t;
#else
   typedef const std::type_info& Info_t;
#endif

   /// Generate emulated collection proxy for a given class
   static TVirtualCollectionProxy* GenEmulatedProxy(const char* class_name, Bool_t silent);

   /// Generate emulated class streamer for a given collection class
   static TClassStreamer* GenEmulatedClassStreamer(const char* class_name, Bool_t silent);

   /// Generate emulated member streamer for a given collection class
   static TMemberStreamer* GenEmulatedMemberStreamer(const char* class_name, Bool_t silent);


   /// Generate proxy from static functions
   static Proxy_t* GenExplicitProxy( const ::ROOT::TCollectionProxyInfo &info, TClass *cl );

   /// Generate proxy from template
   template <class T> static Proxy_t* GenProxy(const T &arg, TClass *cl)  {
      return GenExplicitProxy( ::ROOT::TCollectionProxyInfo::Get(arg), cl );
   }

   /// Generate streamer from static functions
   static TGenCollectionStreamer*
      GenExplicitStreamer( const ::ROOT::TCollectionProxyInfo &info, TClass *cl );

   /// Generate class streamer from static functions
   static TClassStreamer*
      GenExplicitClassStreamer( const ::ROOT::TCollectionProxyInfo &info, TClass *cl );

   /// Generate class streamer from template
   template <class T> static TClassStreamer* GenClassStreamer(const T &arg, TClass *cl)  {
      return GenExplicitClassStreamer(::ROOT::TCollectionProxyInfo::Get(arg), cl);
   }

   /// Generate member streamer from static functions
   static TMemberStreamer*
      GenExplicitMemberStreamer(const ::ROOT::TCollectionProxyInfo &info, TClass *cl);

   /// Generate member streamer from template
   template <class T> static TMemberStreamer* GenMemberStreamer(const T &arg, TClass *cl)  {
      return GenExplicitMemberStreamer(::ROOT::TCollectionProxyInfo::Get(arg), cl);
   }
};

/**
 \class TCollectionStreamer TCollectionProxyFactory.h
 \ingroup IO

 Class streamer object to implement TClassStreamer functionality for I/O emulation.

 @author  M.Frank
 @version 1.0
*/
class TCollectionStreamer   {
private:
   TCollectionStreamer& operator=(const TCollectionStreamer&);   // not implemented

protected:
   TGenCollectionProxy* fStreamer;   ///< Pointer to worker streamer

   /// Issue Error about invalid proxy
   void InvalidProxyError();

public:
   /// Initializing constructor
   TCollectionStreamer();
   /// Copy constructor
   TCollectionStreamer(const TCollectionStreamer& c);
   /// Standard destructor
   virtual ~TCollectionStreamer();
   /// Attach worker proxy
   void AdoptStreamer(TGenCollectionProxy* streamer);
   /// Streamer for I/O handling
   void Streamer(TBuffer &refBuffer, void *obj, int siz, TClass *onFileClass );
};

/**
 \class TCollectionClassStreamer TCollectionProxyFactory.h
 \ingroup IO

 Class streamer object to implement TClassStreamer functionality
 for I/O emulation.
 \author  M.Frank
 \version 1.0
*/
class TCollectionClassStreamer : public TClassStreamer, public TCollectionStreamer {
 protected:
   TCollectionClassStreamer &operator=(const TCollectionClassStreamer &rhs); // Not implemented.
   /// Copy constructor
   TCollectionClassStreamer(const TCollectionClassStreamer& c)
      : TClassStreamer(c), TCollectionStreamer(c)      {                        }

public:
   /// Initializing constructor
   TCollectionClassStreamer() : TClassStreamer(0)     {                        }
   /// Standard destructor
   virtual ~TCollectionClassStreamer()                {                        }
   /// Streamer for I/O handling
   virtual void operator()(TBuffer &buff, void *obj ) { Streamer(buff,obj,0,fOnFileClass); }

   virtual void Stream(TBuffer &b, void *obj, const TClass *onfileClass)
   {
      if (b.IsReading()) {
         TGenCollectionProxy *proxy = TCollectionStreamer::fStreamer;
         if (onfileClass==0 || onfileClass == proxy->GetCollectionClass()) {
            proxy->ReadBuffer(b,obj);
         } else {
            proxy->ReadBuffer(b,obj,onfileClass);
         }
      } else {
         // fStreamer->WriteBuffer(b,objp,onfileClass);
         Streamer(b,obj,0,(TClass*)onfileClass);
      }
   }

   /// Virtual copy constructor.
   virtual TClassStreamer *Generate() const {
      return new TCollectionClassStreamer(*this);
   }

   TGenCollectionProxy *GetXYZ() { return TCollectionStreamer::fStreamer; }

};

/**
 \class TCollectionMemberStreamer TCollectionProxyFactory.h
 \ingroup IO

 Class streamer object to implement TMemberStreamer functionality
 for I/O emulation.
 \author  M.Frank
 \version 1.0
  */
class TCollectionMemberStreamer : public TMemberStreamer, public TCollectionStreamer {
private:
   TCollectionMemberStreamer &operator=(const TCollectionMemberStreamer &rhs); // Not implemented.
public:
   /// Initializing constructor
   TCollectionMemberStreamer() : TMemberStreamer(0) { }
   /// Copy constructor
   TCollectionMemberStreamer(const TCollectionMemberStreamer& c)
      : TMemberStreamer(c), TCollectionStreamer(c)   { }
   /// Standard destructor
   virtual ~TCollectionMemberStreamer()             { }
   /// Streamer for I/O handling
   virtual void operator()(TBuffer &buff,void *obj,Int_t siz=0)
   { Streamer(buff, obj, siz, 0); /* FIXME */ }
};

#endif
