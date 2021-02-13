/* @(#)root/multiproc:$Id$ */
// Author: Enrico Guiraud July 2015

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_MPSendRecv
#define ROOT_MPSendRecv

#include "TBufferFile.h"
#include "TClass.h"
#include "TError.h"
#include "TSocket.h"
#include <memory> //unique_ptr
#include <type_traits> //enable_if
#include <typeinfo> //typeid
#include <utility> //pair
#include <string>

//////////////////////////////////////////////////////////////////////////
/// An std::pair that wraps the code and optional object contained in a message.
/// \param first message code
/// \param second a smart pointer to a TBufferFile that contains the message object\n
/// The smart pointer is null if the message does not contain an object
/// but only consists of a code. See MPRecv() description on how to
/// retrieve the object from the TBufferFile.
using MPCodeBufPair = std::pair<unsigned, std::unique_ptr<TBufferFile>>;


/************ FUNCTIONS' DECLARATIONS *************/

// There are several versions of this function: this is one sends a
// message with a code and no object. The templated versions are used
// to send a code and an object of any non-pointer type.
int MPSend(TSocket *s, unsigned code);

template<class T, typename std::enable_if<std::is_class<T>::value>::type * = nullptr>
int MPSend(TSocket *s, unsigned code, T obj);

template < class T, typename std::enable_if < !std::is_class<T>::value  &&!std::is_pointer<T>::value >::type * = nullptr >
int MPSend(TSocket *s, unsigned code, T obj);

template<class T, typename std::enable_if<std::is_same<const char *, T>::value>::type * = nullptr>
int MPSend(TSocket *s, unsigned code, T obj);

template < class T, typename std::enable_if < std::is_pointer<T>::value  &&std::is_constructible<TObject *, T>::value >::type * = nullptr >
int MPSend(TSocket *s, unsigned code, T obj);

MPCodeBufPair MPRecv(TSocket *s);


//this version reads classes from the message
template<class T, typename std::enable_if<std::is_class<T>::value>::type * = nullptr>
T ReadBuffer(TBufferFile *buf);

//this version reads built-in types from the message
template < class T, typename std::enable_if < !std::is_class<T>::value  &&!std::is_pointer<T>::value >::type * = nullptr >
T ReadBuffer(TBufferFile *buf);

//this version reads std::string and c-strings from the message
template<class T, typename std::enable_if<std::is_same<const char *, T>::value>::type * = nullptr>
T ReadBuffer(TBufferFile *buf);

//this version reads a TObject* from the message
template < class T, typename std::enable_if < std::is_pointer<T>::value  &&std::is_constructible<TObject *, T>::value >::type * = nullptr >
T ReadBuffer(TBufferFile *buf);


/************ TEMPLATE FUNCTIONS' IMPLEMENTATIONS *******************/

//////////////////////////////////////////////////////////////////////////
/// Send a message with a code and an object to socket s.
/// The number of bytes sent is returned, as per TSocket::SendRaw.
/// This standalone function can be used to send a code and possibly
/// an object on a given socket. This function does not check whether the
/// socket connection is in a valid state. MPRecv() must be used to
/// retrieve the contents of the message.\n
/// **Note:** only objects the headers of which have been parsed by
/// cling can be sent using MPSend(). User-defined types can be made available to
/// cling via a call like `gSystem->ProcessLine("#include \"header.h\"")`.
/// Pointer types cannot be sent via MPSend() (with the exception of const char*).
/// \param s a pointer to a valid TSocket. No validity checks are performed\n
/// \param code the code to be sent
/// \param obj the object to be sent
/// \return the number of bytes sent, as per TSocket::SendRaw
template<class T, typename std::enable_if<std::is_class<T>::value>::type *>
int MPSend(TSocket *s, unsigned code, T obj)
{
   TClass *c = TClass::GetClass<T>();
   if (!c) {
      Error("MPSend", "[E] Could not find cling definition for class %s\n", typeid(T).name());
      return -1;
   }
   TBufferFile objBuf(TBuffer::kWrite);
   objBuf.WriteObjectAny(&obj, c);
   TBufferFile wBuf(TBuffer::kWrite);
   wBuf.WriteUInt(code);
   wBuf.WriteULong(objBuf.Length());
   wBuf.WriteBuf(objBuf.Buffer(), objBuf.Length());
   return s->SendRaw(wBuf.Buffer(), wBuf.Length());
}

/// \cond
// send a built-in type that is not a pointer (under the hypothesis that
// TBuffer's operator<< works with any built-in type that is not a pointer)
template < class T, typename std::enable_if < !std::is_class<T>::value  &&!std::is_pointer<T>::value >::type * >
int MPSend(TSocket *s, unsigned code, T obj)
{
   TBufferFile wBuf(TBuffer::kWrite);
   ULong_t size = sizeof(T);
   wBuf << code << size << obj;
   return s->SendRaw(wBuf.Buffer(), wBuf.Length());
}

// send an null-terminated c-string or an std::string (which is converted to a c-string)
//TODO can this become a partial specialization instead?
template<class T, typename std::enable_if<std::is_same<const char *, T>::value>::type *>
int MPSend(TSocket *s, unsigned code, T obj)
{
   TBufferFile wBuf(TBuffer::kWrite);
   wBuf.WriteUInt(code);
   wBuf.WriteULong(strlen(obj) + 1); //strlen does not count the trailing \0
   wBuf.WriteString(obj);
   return s->SendRaw(wBuf.Buffer(), wBuf.Length());
}

// send a TObject*. Allows polymorphic behaviour and pters to derived classes
template < class T, typename std::enable_if < std::is_pointer<T>::value && std::is_constructible<TObject *, T>::value >::type * >
int MPSend(TSocket *s, unsigned code, T obj)
{
   //find out the size of the object
   TBufferFile objBuf(TBuffer::kWrite);
   if(obj != nullptr)
      objBuf.WriteObjectAny(obj, obj->IsA());

   //write everything together in a buffer
   TBufferFile wBuf(TBuffer::kWrite);
   wBuf.WriteUInt(code);
   wBuf.WriteULong(objBuf.Length());
   if(objBuf.Length())
      wBuf.WriteBuf(objBuf.Buffer(), objBuf.Length());
   return s->SendRaw(wBuf.Buffer(), wBuf.Length());
}

/// \endcond

//////////////////////////////////////////////////////////////////////////
/// One of the template functions used to read objects from messages.
/// Different implementations are provided for different types of objects:
/// classes, non-pointer built-ins and const char*. Reading pointers is
/// not implemented (at the time of writing, sending pointers is not either).
template<class T, typename std::enable_if<std::is_class<T>::value>::type *>
T ReadBuffer(TBufferFile *buf)
{
   TClass *c = TClass::GetClass(typeid(T));
   T *objp = (T *)buf->ReadObjectAny(c);
   T obj = *objp; //this is slow, but couldn't find a better way of returning a T without leaking memory
   delete objp;
   return obj;
}

/// \cond
template < class T, typename std::enable_if < !std::is_class<T>::value  &&!std::is_pointer<T>::value >::type * >
T ReadBuffer(TBufferFile *buf)
{
   //read built-in type
   T obj;
   *(buf) >> obj;
   return obj;
}

template<class T, typename std::enable_if<std::is_same<const char *, T>::value>::type *>
T ReadBuffer(TBufferFile *buf)
{
   //read c-string
   char *c = new char[buf->BufferSize()];
   buf->ReadString(c, buf->BufferSize());
   return c;
}

template < class T, typename std::enable_if < std::is_pointer<T>::value  &&std::is_constructible<TObject *, T>::value >::type * >
T ReadBuffer(TBufferFile *buf)
{
   //read TObject*
   using objType = typename std::remove_pointer<T>::type;
   return (T)buf->ReadObjectAny(objType::Class());
}
/// \endcond

#endif
