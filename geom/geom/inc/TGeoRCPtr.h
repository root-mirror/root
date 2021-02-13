// @(#)root/geom:$Id$
// Author: Andrei Gheata   29/05/2013

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoRCPtr
#define ROOT_TGeoRCPtr

/** \class TGeoRCPtr
\ingroup Geometry_classes

A reference counting-managed pointer for classes derived from TGeoExtension
which can be used as C pointer. Based on CodeProject implementation example

### Example:

~~~ {.cpp}
class MyExtension : public TGeoExtension {
public:
   MyExtension() : TGeoExtension(), fRC(0) {printf("Created MyExtension\n");}
   virtual ~MyExtension() {printf("Deleted MyExtension\n");}

   virtual TGeoExtension *Grab() const {fRC++; return (TGeoExtension*)this;}
   virtual void Release() const {assert(fRC > 0); fRC--; if (fRC ==0) delete this;}
   void print() const {printf("MyExtension object %p\n", this);}
private:
   mutable Int_t        fRC;           // Reference counter
   ClassDef(MyExtension,1)
};
~~~

### Usage:

~~~ {.cpp}
 // Module 1 creates an object
 TGeoRCPtr<MyExtension> a2 = new MyExtension(); //fRC=1

 // Module 2 grabs object
 TGeoRCPtr<MyExtension> ptr2 = a2; //fRC=2

 // Module 2 invokes a method
 ptr2->Print();
 (*ptr2).Print();

 // Module 1 no longer needs object
  a2 = 0;      //RC=1

 // Module 2 no longer needs object
  ptr2 = 0;    //object will be destroyed here
~~~

### Note:

Event if one forgets to call ptr2 = 0, the object gets delete when the method
using ptr2 gets out of scope.
*/

template<class T>
class TGeoRCPtr
{
public:
   //Construct using a C pointer, e.g. TGeoRCPtr<T> x = new T();
   TGeoRCPtr(T* ptr = 0)
       : fPtr(ptr)
   {
      if(ptr != 0) ptr->Grab();
   }

   //Copy constructor
   TGeoRCPtr(const TGeoRCPtr &ptr)
       : fPtr(ptr.fPtr)
   {
      if(fPtr != 0) fPtr->Grab();
   }

   ~TGeoRCPtr()
   {
      if(fPtr != 0) fPtr->Release();
   }

   //Assign a pointer, e.g. x = new T();
   TGeoRCPtr &operator=(T* ptr)
   {
      if(ptr != 0) ptr->Grab();
      if(fPtr != 0) fPtr->Release();
      fPtr = ptr;
      return (*this);
   }

   //Assign another TGeoRCPtr
   TGeoRCPtr &operator=(const TGeoRCPtr &ptr)
   {
      return (*this) = ptr.fPtr;
   }

   //Retrieve actual pointer
   T* Get() const
   {
      return fPtr;
   }

   //Some overloaded operators to facilitate dealing with an TGeoRCPtr as a conventional C pointer.
   //Without these operators, one can still use the less transparent Get() method to access the pointer.
   T* operator->() const {return fPtr;}   //x->member
   T &operator*() const {return *fPtr;}   //*x, (*x).member
   operator T*() const {return fPtr;}      //T* y = x;
   operator bool() const {return fPtr != 0;}   //if(x) {/*x is not NULL*/}
   bool operator==(const TGeoRCPtr &ptr) {return fPtr == ptr.fPtr;}
   bool operator==(const T *ptr) {return fPtr == ptr;}

private:
   T *fPtr;   //Actual pointer
};

#endif
