/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooRealProxy.h,v 1.23 2007/07/12 20:30:28 wouter Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_TEMPLATE_PROXY
#define ROO_TEMPLATE_PROXY

#include "RooAbsReal.h"
#include "RooArgProxy.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategory.h"
#include "RooMsgService.h"
#include <string>

/**
\class RooTemplateProxy
\ingroup Roofitcore

## Introduction
A RooTemplateProxy is used to hold references to other RooFit objects in an expression tree.
A `RooGaussian(..., x, mean, sigma)` can e.g. store references to `x, mean, sigma` as
```
RooTemplateProxy<RooAbsReal> _x;
RooTemplateProxy<RooAbsReal> _mean;
RooTemplateProxy<RooAbsReal> _sigma;
```
Now, the values of these three can be accessed, and the template argument ensures that only objects that evaluate
to real numbers (RooAbsReal) can be stored in such a proxy. These can e.g. be variables, PDFs and functions.
To store an object that's a `RooCategory`, one would, for example, use
```
RooTemplateProxy<RooCategory> _category;
```

Since %ROOT 6.22, the proxy can be used like a pointer to an instance of the template argument.
For this, it provides `operator*` and `operator->`, e.g.
```
double oldValue = _x->getVal(normalisationSet);
*_x = 17.;
```

RooTemplateProxy's base class RooArgProxy registers the proxied objects as "servers" of the object
that holds the proxy. When the value of the proxied object is changed, the owner is
notified, and can recalculate its own value. Renaming or exchanging objects that
serve values to the owner of the proxy is handled automatically.

## Modernisation of proxies in %ROOT 6.22
In ROOT 6.22, the classes RooRealProxy and RooCategoryProxy were replaced by RooTemplateProxy<class T>.

Two typedefs have been defined for backward compatibility:
- `RooRealProxy = RooTemplateProxy<RooAbsReal>`. Any generic object that converts to a real value.
- `RooCategoryProxy = RooTemplateProxy<RooAbsCategory>`. Any category object.

To modernise a class, one can change the template argument of the proxy to the most appropriate type,
and increment the class version of the owner.

<table>
<tr><th> %RooFit before %ROOT 6.22 <th> %RooFit starting with %ROOT 6.22
<tr><td>
~~~{.cpp}
// In .h: Declare member
RooRealProxy pdfProxy;

ClassDef(MyPdf, 1)
};

// In .cxx: Initialise proxy in constructor
// The proxy will accept any RooAbsArg, so the type of
// "thePdf" has to be checked manually.
MyPdf::MyPdf(name, title, ...) :
  pdfProxy("pdfProxy", "Proxy holding a PDF", this, thePdf) {
  [ Extra checking here ... ]
}


// In .cxx: Accessing the proxy
RooAbsArg* absArg = pdfProxy.absArg();
RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(absArg);
assert(pdf); // Manual type checking ...
pdf->fitTo(...);
~~~
<td>
~~~{.cpp}
// In .h: Declare member
RooTemplateProxy<RooAbsPdf> pdfProxy;

ClassDef(MyPdf, 2)
};

// In .cxx: Initialise proxy in constructor
// The program will not compile if "thePdf" is not a
// type deriving from RooAbsPdf
MyPdf::MyPdf(name, title, ...) :
  pdfProxy("pdfProxy", "Proxy holding a PDF", this, thePdf) {

}


// In .cxx: Accessing the proxy



pdfProxy->fitTo(...);
~~~
</table>


### How to modernise old code

1. Choose the proper template argument for the proxy.
 - If a PDF is stored: `RooTemplateProxy<RooAbsPdf>`.
 - If a real-valued object is stored: `RooTemplateProxy<RooAbsReal>`.
 - If a category is stored: `RooTemplateProxy<RooCategory>`.
 - If a variable is stored (i.e. one wants to be able to assign values to it): `RooTemplateProxy<RooRealVar>`
 Other template arguments are possible, as long as they derive from RooAbsArg.
2. Increment the class version of the owning class.
3. Make sure that the right type is passed in the constructor of the proxy.
4. Always use `proxy->` and `*proxy` to work with the stored object. No need to cast.
5. **Only if necessary** If errors about missing symbols connected to RooTemplateProxy appear at link time,
   a specific template instantiation for RooTemplateProxy is not yet in ROOT's dictionaries.
   These two lines should be added to the LinkDef.h of the project:
   ~~~{.cpp}
       #pragma link C++ class RooTemplateProxy<RooMultiCategory>+;
       #pragma read sourceClass="RooCategoryProxy" targetClass="RooTemplateProxy<RooMultiCategory>"
   ~~~
   Replace `RooMultiCategory` by the proper type. If the proxy was holding a real-valued object, use `sourceClass="RooRealProxy"`.

   The first line adds the proxy class to the dictionary, the second line enables reading a legacy
   `RooCategoryProxy` from a file, and converting it to the new type-safe proxy. If no old proxies
   have to be read from files, this line can be omitted.

   If the template instantiation that triggered the missing symbols seems to be a very common instantiation,
   request for it to be added to RooFit by creating a pull request for ROOT. If it is rather uncommon,
   it is sufficient to add it to the LinkDef.h of the local project only.

**/

template<class T>
class RooTemplateProxy : public RooArgProxy {
public:

  RooTemplateProxy() {} ;

  ////////////////////////////////////////////////////////////////////////////////
  /// Constructor with owner.
  /// \param[in] theName Name of this proxy (for printing).
  /// \param[in] desc Description what this proxy should act as.
  /// \param[in] owner The object that owns the proxy. This is important for tracking
  ///            of client-server dependencies.
  /// \param[in] valueServer Notify the owner if value changes.
  /// \param[in] shapeServer Notify the owner if shape (e.g. binning) changes.
  /// \param[in] proxyOwnsArg Proxy will delete the payload if owning.
  RooTemplateProxy(const char* theName, const char* desc, RooAbsArg* owner,
      Bool_t valueServer=true, Bool_t shapeServer=false, Bool_t proxyOwnsArg=false)
  : RooArgProxy(theName, desc, owner, valueServer, shapeServer, proxyOwnsArg) { }

  ////////////////////////////////////////////////////////////////////////////////
  /// Constructor with owner and proxied object.
  /// \param[in] theName Name of this proxy (for printing).
  /// \param[in] desc Description what this proxy should act as.
  /// \param[in] owner The object that owns the proxy. This is important for tracking
  ///            of client-server dependencies.
  /// \param[in] ref Reference to the object that the proxy should hold.
  /// \param[in] valueServer Notify the owner if value changes.
  /// \param[in] shapeServer Notify the owner if shape (e.g. binning) changes.
  /// \param[in] proxyOwnsArg Proxy will delete the payload if owning.
  RooTemplateProxy(const char* theName, const char* desc, RooAbsArg* owner, T& ref,
      Bool_t valueServer=true, Bool_t shapeServer=false, Bool_t proxyOwnsArg=false) :
        RooArgProxy(theName, desc, owner, ref, valueServer, shapeServer, proxyOwnsArg) { }


  ////////////////////////////////////////////////////////////////////////////////
  /// Copy from an existing proxy.
  /// It will accept any RooTemplateProxy instance, and attempt a dynamic_cast on its payload.
  /// \param[in] theName Name of this proxy.
  /// \param[in] owner Pointer to the owner this proxy should be registered to.
  /// \param[in] other Instance of a differen proxy whose payload should be copied.
  /// \param[in] allowWrongTypes Instead of throwing a std::invalid_argument, only issue an
  /// error message when payload with wrong type is found. This is unsafe, but may be necessary
  /// when reading back legacy types. Defaults to false.
  /// \throw std::invalid_argument if the types of the payloads are incompatible.
  template<typename U>
  RooTemplateProxy(const char* theName, RooAbsArg* owner, const RooTemplateProxy<U>& other, bool allowWrongTypes = false) :
    RooArgProxy(theName, owner, other) {
    if (_arg && !dynamic_cast<const T*>(_arg)) {
      if (allowWrongTypes) {
        coutE(InputArguments) << "Error trying to copy an argument from a proxy with an incompatible payload." << std::endl;
      } else {
        throw std::invalid_argument("Tried to construct a RooTemplateProxy with incompatible payload.");
      }
    }
  }

  virtual TObject* Clone(const char* newName=0) const { return new RooTemplateProxy<T>(newName,_owner,*this); }
 

  /// Return reference to the proxied object.
  T& operator*() const {
    return static_cast<T&>(*_arg);
  }

  /// Member access operator to proxied object.
  T* operator->() const {
    return static_cast<T*>(_arg);
  }


  /// Convert the proxy into a number.
  /// \return A category proxy will return the index state, real proxies the result of RooAbsReal::getVal(normSet).
  operator typename T::value_type() const {
    return retrieveValue(arg());
  }


  ////////////////////////////////////////////////////////////////////////////////
  /// Change object held in proxy into newRef
  bool setArg(T& newRef) {
    if (_arg) {
      if (std::string(arg().GetName()) != newRef.GetName()) {
        newRef.setAttribute(Form("ORIGNAME:%s", arg().GetName())) ;
      }
      return changePointer(RooArgSet(newRef), true);
    } else {
      return changePointer(RooArgSet(newRef), false, true);
    }
  }


  /// \name Legacy interface
  /// In ROOT versions before 6.22, RooFit didn't have this typed proxy. Therefore, a number of functions
  /// for forwarding calls to the proxied objects were necessary. The functions in this group can all be
  /// replaced by directly accessing the proxied objects using e.g. the member access operator like
  /// `proxy->function()` or by dereferencing like `*proxy = value`.
  /// For this to work, choose the template argument appropriately. That is, if the
  /// proxy stores a PDF, use `RooTemplateProxy<RooAbsPdf>`, *etc.*.
  /// @{

  /// Get the label of the current category state. This function only makes sense for category proxies.
  const char* label() const {
    return arg().getCurrentLabel();
  }

  /// Check if the stored object has a range with the given name.
  bool hasRange(const char* rangeName) const {
    return arg().hasRange(rangeName);
  }

  /// Return reference to object held in proxy.
  const T& arg() const { return static_cast<const T&>(*_arg); }

  /// Assign a new value to the object pointed to by the proxy.
  /// This requires the payload to be assignable (RooAbsRealLValue or derived, RooAbsCategoryLValue).
  RooTemplateProxy<T>& operator=(typename T::value_type value) {
    lvptr(static_cast<T*>(nullptr))->operator=(value);
    return *this;
  }
  /// Set a category state using its state name. This function can only work for category-type proxies.
  RooTemplateProxy<T>& operator=(const std::string& newState) {
    static_assert(std::is_base_of<RooAbsCategory, T>::value, "Strings can only be assigned to category proxies.");
    lvptr(static_cast<RooAbsCategoryLValue*>(nullptr))->operator=(newState.c_str());
    return *this;
  }

  /// Query lower limit of range. This requires the payload to be RooAbsRealLValue or derived.
  double min(const char* rname=0) const  { return lvptr(static_cast<const T*>(nullptr))->getMin(rname) ; }
  /// Query upper limit of range. This requires the payload to be RooAbsRealLValue or derived.
  double max(const char* rname=0) const  { return lvptr(static_cast<const T*>(nullptr))->getMax(rname) ; }
  /// Check if the range has a lower bound. This requires the payload to be RooAbsRealLValue or derived.
  bool hasMin(const char* rname=0) const { return lvptr(static_cast<const T*>(nullptr))->hasMin(rname) ; }
  /// Check if the range has a upper bound. This requires the payload to be RooAbsRealLValue or derived.
  bool hasMax(const char* rname=0) const { return lvptr(static_cast<const T*>(nullptr))->hasMax(rname) ; }

  /// @}


private:
  /// Are we a real-valued proxy or a category proxy?
  using LValue_t = typename std::conditional<std::is_base_of<RooAbsReal, T>::value,
      RooAbsRealLValue, RooAbsCategoryLValue>::type;

  ////////////////////////////////////////////////////////////////////////////////
  /// Return l-value pointer to contents. If the contents derive from RooAbsLValue or RooAbsCategoryLValue,
  /// the conversion is safe, and the function directly returns the pointer using a static_cast.
  /// If the template parameter of this proxy is not an LValue type, then
  /// - in a debug build, a dynamic_cast with an assertion is used.
  /// - in a release build, a static_cast is forced, irrespective of what the type of the object actually is. This
  /// is dangerous, but equivalent to the behaviour before refactoring the RooFit proxies.
  /// \deprecated This function is unneccessary if the template parameter is RooAbsRealLValue (+ derived types) or
  /// RooAbsCategoryLValue (+derived types), as arg() will always return the correct type.
  const LValue_t* lvptr(const LValue_t*) const {
    return static_cast<const LValue_t*>(_arg);
  }
  /// \copydoc lvptr(const LValue_t*) const
  LValue_t* lvptr(LValue_t*) {
    return static_cast<LValue_t*>(_arg);
  }
  /// \copydoc lvptr(const LValue_t*) const
  const LValue_t* lvptr(const RooAbsArg*) const
  R__SUGGEST_ALTERNATIVE("The template argument of RooTemplateProxy needs to derive from RooAbsRealLValue or RooAbsCategoryLValue to safely call this function.") {
#ifdef NDEBUG
    return static_cast<const LValue_t*>(_arg);
#else
    auto theArg = dynamic_cast<const LValue_t*>(_arg);
    assert(theArg);
    return theArg;
#endif
  }
  /// \copydoc lvptr(const LValue_t*) const
  LValue_t* lvptr(RooAbsArg*)
  R__SUGGEST_ALTERNATIVE("The template argument of RooTemplateProxy needs to derive from RooAbsRealLValue or RooAbsCategoryLValue to safely call this function.") {
#ifdef NDEBUG
    return static_cast<LValue_t*>(_arg);
#else
    auto theArg = dynamic_cast<LValue_t*>(_arg);
    assert(theArg);
    return theArg;
#endif
  }


  /// Retrieve index state from a category.
  typename T::value_type retrieveValue(const RooAbsCategory& cat) const {
    return cat.getCurrentIndex();
  }

  /// Retrieve value from a real-valued object.
  typename T::value_type retrieveValue(const RooAbsReal& real) const {
    return real.getVal(_nset);
  }

  ClassDef(RooTemplateProxy,1) // Proxy for a RooAbsReal object
};

#endif
