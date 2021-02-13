/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
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

/// \class RooMappedCategory
/// RooMappedCategory provides a category-to-category mapping defined
/// by pattern matching on their state labels.
///
/// The mapping function consists of a series of wild card regular expressions.
/// Each expression is matched to the input categories' state labels, and an associated
/// output state label.

#include "RooMappedCategory.h"

#include "RooFit.h"
#include "RooStreamParser.h"
#include "RooMsgService.h"
#include "Riostream.h"
#include "RooAbsCache.h"

#include "TBuffer.h"
#include "TString.h"
#include "TRegexp.h"


class RooMappedCategoryCache : public RooAbsCache {
  public:
    RooMappedCategoryCache(RooAbsArg* owner = 0) : RooAbsCache(owner)
  { initialise(); }
    RooMappedCategoryCache(const RooAbsCache& other, RooAbsArg* owner = 0) :
      RooAbsCache(other, owner)
    { initialise(); }

    // look up our parent's output based on our parent's input category index
    RooAbsCategory::value_type lookup(Int_t idx) const
    { return _map[idx]; }

    virtual void wireCache()
    { _map.clear(); initialise(); }

    virtual Bool_t redirectServersHook(const RooAbsCollection& /*newServerList*/, Bool_t /*mustReplaceAll*/, Bool_t /*nameChange*/, Bool_t /*isRecursive*/)
    { _map.clear(); initialise(); return kFALSE; }

  private:
    mutable std::map<Int_t, RooAbsCategory::value_type> _map;

    // pre-map categories of input category on something easily searchable
    // like the index (not the name!)
    void initialise()
    {
      const RooMappedCategory& parent = *static_cast<const RooMappedCategory*>(_owner);

      for (const auto& inCat : *parent._inputCat) {
        const std::string& inKey = inCat.first;
        // Scan array of regexps
        bool found = false;
        for (const auto& strAndEntry : parent._mapArray) {
          if (strAndEntry.second.match(inKey.c_str())) {
            found = true;
            _map[inCat.second] = strAndEntry.second.outCat();
            break;
          }
        }

        if (!found)
          _map[inCat.second] = parent._defCat;
      }
    }
};

RooMappedCategory::RooMappedCategory(const char *name, const char *title, RooAbsCategory& inputCat, const char* defOut, Int_t defOutIdx) :
  RooAbsCategory(name, title), _inputCat("input","Input category",this,inputCat),
  _mapcache(0)
{
  // Constructor with input category and name of default output state, which is assigned
  // to all input category states that do not follow any mapping rule.
  if (defOutIdx==NoCatIdx) {
    _defCat = defineState(defOut).second;
  } else {
    _defCat = defineState(defOut,defOutIdx).second;
  }
}


RooMappedCategory::RooMappedCategory(const RooMappedCategory& other, const char *name) :
  RooAbsCategory(other,name), _inputCat("input",this,other._inputCat), _mapArray(other._mapArray),
  _mapcache(0)
{
  _defCat = lookupIndex(other.lookupName(other._defCat));
}



RooMappedCategory::~RooMappedCategory()
{
    // Destructor
    delete _mapcache;
}



Bool_t RooMappedCategory::map(const char* inKeyRegExp, const char* outKey, Int_t outIdx)
{
  // Add mapping rule: any input category state label matching the 'inKeyRegExp'
  // wildcard expression will be mapped to an output state with name 'outKey'
  //
  // Rules are evaluated in the order they were added. In case an input state
  // matches more than one rule, the first rules output state will be assigned

  if (!inKeyRegExp || !outKey) return kTRUE ;

  // Check if pattern is already registered
  if (_mapArray.find(inKeyRegExp)!=_mapArray.end()) {
    coutE(InputArguments) << "RooMappedCategory::map(" << GetName() << "): ERROR expression "
                          << inKeyRegExp << " already mapped" << std::endl ;
    return kTRUE ;
  }

  // Check if output type exists, if not register
  value_type catIdx = lookupIndex(outKey);
  if (catIdx == invalidCategory().second) {
    if (outIdx==NoCatIdx) {
      catIdx = defineState(outKey).second;
    } else {
      catIdx = defineState(outKey,outIdx).second;
    }
  }

  if (catIdx == invalidCategory().second) {
    coutE(InputArguments) << "RooMappedCategory::map(" << GetName()
                          << "): ERROR, unable to define category for output type " << outKey << std::endl ;
    return true;
  }

  // Create new map entry ;
  Entry e(inKeyRegExp, catIdx);
  if (!e.ok()) {
    coutE(InputArguments) << "RooMappedCategory::map(" << GetName()
                          << "): ERROR, expression " << inKeyRegExp << " didn't compile" << std::endl ;
    return kTRUE ;
  }

  _mapArray[inKeyRegExp] = e ;
  return kFALSE ;
}



RooAbsCategory::value_type RooMappedCategory::evaluate() const
{
    const RooMappedCategoryCache* cache = getOrCreateCache();
    return cache->lookup(_inputCat->getCurrentIndex());
}

const RooMappedCategoryCache* RooMappedCategory::getOrCreateCache() const
{
    if (!_mapcache) _mapcache = new RooMappedCategoryCache(
            const_cast<RooMappedCategory*>(this));
    return _mapcache;
}

void RooMappedCategory::printMultiline(std::ostream& os, Int_t content, Bool_t verbose, TString indent) const
{
  // Print info about this mapped category to the specified stream. In addition to the info
  // from RooAbsCategory::printStream() we add:
  //
  //  Standard : input category
  //     Shape : default value
  //   Verbose : list of mapping rules

  RooAbsCategory::printMultiline(os,content,verbose,indent);

  if (verbose) {
    os << indent << "--- RooMappedCategory ---" << std::endl
       << indent << "  Maps from " ;
    _inputCat.arg().printStream(os,0,kStandard);

    os << indent << "  Default value is " << lookupName(_defCat) << " = " << _defCat << '\n';

    os << indent << "  Mapping rules:" << std::endl;
    for (const auto& strAndEntry : _mapArray) {
      os << indent << "  " << strAndEntry.first << " -> " << strAndEntry.second.outCat() << std::endl ;
    }
  }
}


Bool_t RooMappedCategory::readFromStream(std::istream& is, Bool_t compact, Bool_t /*verbose*/)
{
  // Read object contents from given stream
   if (compact) {
     coutE(InputArguments) << "RooMappedCategory::readFromSteam(" << GetName() << "): can't read in compact mode" << std::endl ;
     return kTRUE ;
   } else {

     //Clear existing definitions, but preserve default output
     TString defCatName(lookupName(_defCat));
     _mapArray.clear() ;
     delete _mapcache;
     _mapcache = 0;
     clearTypes() ;
     _defCat = defineState(defCatName.Data()).second;

     TString token,errorPrefix("RooMappedCategory::readFromStream(") ;
     errorPrefix.Append(GetName()) ;
     errorPrefix.Append(")") ;
     RooStreamParser parser(is,errorPrefix) ;
     parser.setPunctuation(":,") ;

     TString destKey,srcKey ;
     Bool_t readToken(kTRUE) ;

    // Loop over definition sequences
     while(1) {
       if (readToken) token=parser.readToken() ;
       if (token.IsNull()) break ;
       readToken=kTRUE ;

       destKey = token ;
       if (parser.expectToken(":",kTRUE)) return kTRUE ;

       // Loop over list of sources for this destination
       while(1) {
         srcKey = parser.readToken() ;
         token = parser.readToken() ;

         // Map a value
         if (map(srcKey,destKey)) return kTRUE ;

         // Unless next token is ',' current token
         // is destination part of next sequence
         if (token.CompareTo(",")) {
           readToken = kFALSE ;
           break ;
         }
       }
     }
     return kFALSE ;
   }
   //return kFALSE ; // statement unreachable (OSF)
}


////////////////////////////////////////////////////////////////////////////////
/// Customized printing of arguments of a RooMappedCategory to more intuitively reflect the contents of the
/// product operator construction

void RooMappedCategory::printMetaArgs(std::ostream& os) const
{
  // Scan array of regexps
  RooAbsCategory::value_type prevOutCat = invalidCategory().second;
  Bool_t first(kTRUE) ;
  os << "map=(" ;
  for (const auto& iter : _mapArray) {
    if (iter.second.outCat() != prevOutCat) {
      if (!first) { os << " " ; }
      first=kFALSE ;

      os << iter.second.outCat() << ":" << iter.first ;
      prevOutCat = iter.second.outCat();
    } else {
      os << "," << iter.first ;
    }
  }

  if (!first) { os << " " ; }
  os << lookupName(_defCat) << ":*" ;

  os << ") " ;
}




void RooMappedCategory::writeToStream(std::ostream& os, Bool_t compact) const
{
  // Write object contents to given stream
  if (compact) {
    // Write value only
    os << getCurrentLabel() ;
  } else {
    // Write mapping expression

    // Scan array of regexps
    RooAbsCategory::value_type prevOutCat = invalidCategory().second;
    Bool_t first(kTRUE) ;
    for (const auto& iter : _mapArray) {
      if (iter.second.outCat() != prevOutCat) {
        if (!first) { os << " " ; }
        first=kFALSE ;

        os << iter.second.outCat() << "<-" << iter.first ;
        prevOutCat = iter.second.outCat();
      } else {
        os << "," << iter.first ;
      }
    }

    if (!first) { os << " " ; }
    os << lookupName(_defCat) << ":*" ;
  }
}


/// When the input category changes states, the cached state mappings are invalidated
void RooMappedCategory::recomputeShape() {
  // There is no need to recompute _stateNames and _insertionOrder, as only defining new
  // mappings has an effect on these. When the input category changes it shape, it is sufficient
  // to clear the cached state mappings.
  if (_mapcache) {
    _mapcache->wireCache();
  }
}


RooMappedCategory::Entry::Entry(const char* exp, RooAbsCategory::value_type cat) :
    _expr(exp), _regexp(nullptr), _catIdx(cat) {}
RooMappedCategory::Entry::Entry(const Entry& other) :
    _expr(other._expr), _regexp(nullptr), _catIdx(other._catIdx) {}
RooMappedCategory::Entry::~Entry() { delete _regexp; }
bool RooMappedCategory::Entry::ok() { return (const_cast<TRegexp*>(regexp())->Status()==TRegexp::kOK) ; }

////////////////////////////////////////////////////////////////////////////////

RooMappedCategory::Entry& RooMappedCategory::Entry::operator=(const RooMappedCategory::Entry& other)
{
  if (&other==this) return *this ;

  _expr = other._expr ;
  _catIdx = other._catIdx;

  if (_regexp) {
    delete _regexp ;
  }

  return *this;
}

bool RooMappedCategory::Entry::match(const char* testPattern) const {
  return (TString(testPattern).Index(*regexp())>=0);
}

////////////////////////////////////////////////////////////////////////////////
/// Mangle name : escape regexp character '+'

TString RooMappedCategory::Entry::mangle(const char* exp) const
{
  TString t ;
  const char *c = exp ;
  while(*c) {
    if (*c=='+') t.Append('\\') ;
    t.Append(*c) ;
    c++ ;
  }
  return t ;
}

const TRegexp* RooMappedCategory::Entry::regexp() const {
  if (!_regexp) {
    _regexp = new TRegexp(mangle(_expr), true);
  }

  return _regexp;
}
