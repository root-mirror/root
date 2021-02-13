/******************************************************
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

//////////////////////////////////////////////////////////////////////////////
/**  \class RooAbsArg
     \ingroup Roofitcore

RooAbsArg is the common abstract base class for objects that
represent a value and a "shape" in RooFit. Values or shapes usually depend on values
or shapes of other RooAbsArg instances. Connecting several RooAbsArg in
a computation graph models an expression tree that can be evaluated.

### Building a computation graph of RooFit objects
Therefore, RooAbsArg provides functionality to connect objects of type RooAbsArg into
a computation graph to pass values between those objects.
A value can e.g. be a real-valued number, (instances of RooAbsReal), or an integer, that is,
catgory index (instances of RooAbsCategory). The third subclass of RooAbsArg is RooStringVar,
but it is rarely used.

The "shapes" that a RooAbsArg can possess can e.g. be the definition
range of an observable, or how many states a category object has. In computations,
values are expected to change often, while shapes remain mostly constant
(unless e.g. a new range is set for an observable).

Nodes of a computation graph are connected using instances of RooAbsProxy.
If Node B declares a member `RooTemplateProxy<TypeOfNodeA>`, Node A will be
registered as a server of values to Node B, and Node B will know that it is
a client of node A. Using functions like dependsOn(), or getObservables()
/ getParameters(), the relation of `A --> B` can be queried. Using graphVizTree(),
one can create a visualisation of the expression tree.


An instance of RooAbsArg can have named attributes. It also has flags
to indicate that either its value or its shape were changed (= it is dirty).
RooAbsArg provides functionality to manage client/server relations in
a computation graph (\ref clientServerInterface), and helps propagating
value/shape changes through the graph. RooAbsArg implements interfaces
for inspecting client/server relationships (\ref clientServerInterface) and
setting/clearing/querying named attributes.

### Caching of values
The values of nodes in the computation graph are cached in RooFit. If
a value is used in two nodes of a graph, it doesn't need to be recomputed. If
a node acquires a new value, it notifies its consumers ("clients") that
their cached values are dirty. See the functions in \ref optimisationInterface
for details.
A node uses its isValueDirty() and isShapeDirty() functions to decide if a
computation is necessary. Caching can be vetoed globally by setting a
bit using setDirtyInhibit(). This will make computations slower, but all the
nodes of the computation graph will be evaluated irrespective of whether their
state is clean or dirty. Using setOperMode(), caching can also be enabled/disabled
for single nodes.

*/

#include "RooFit.h"

#include "TBuffer.h"
#include "TClass.h"
#include "TVirtualStreamerInfo.h"
#include "strlcpy.h"

#include "RooSecondMoment.h"
#include "RooNameSet.h"
#include "RooWorkspace.h"

#include "RooMsgService.h"
#include "RooAbsArg.h"
#include "RooArgSet.h"
#include "RooArgProxy.h"
#include "RooSetProxy.h"
#include "RooListProxy.h"
#include "RooAbsData.h"
#include "RooAbsCategoryLValue.h"
#include "RooAbsRealLValue.h"
#include "RooTrace.h"
#include "RooRealIntegral.h"
#include "RooConstVar.h"
#include "RooExpensiveObjectCache.h"
#include "RooAbsDataStore.h"
#include "RooResolutionModel.h"
#include "RooVectorDataStore.h"
#include "RooTreeDataStore.h"
#include "ROOT/RMakeUnique.hxx"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>

using namespace std;

#if (__GNUC__==3&&__GNUC_MINOR__==2&&__GNUC_PATCHLEVEL__==3)
char* operator+( streampos&, char* );
#endif

ClassImp(RooAbsArg);
;

Bool_t RooAbsArg::_verboseDirty(kFALSE) ;
Bool_t RooAbsArg::_inhibitDirty(kFALSE) ;
Bool_t RooAbsArg::inhibitDirty() const { return _inhibitDirty && !_localNoInhibitDirty; }

std::map<RooAbsArg*,std::unique_ptr<TRefArray>> RooAbsArg::_ioEvoList;
std::stack<RooAbsArg*> RooAbsArg::_ioReadStack ;


////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooAbsArg::RooAbsArg()
   : TNamed(), _deleteWatch(kFALSE), _valueDirty(kTRUE), _shapeDirty(kTRUE), _operMode(Auto), _fast(kFALSE), _ownedComponents(nullptr),
     _prohibitServerRedirect(kFALSE), _namePtr(0), _isConstant(kFALSE), _localNoInhibitDirty(kFALSE),
     _myws(0)
{
  _namePtr = (TNamed*) RooNameReg::instance().constPtr(GetName()) ;

}

////////////////////////////////////////////////////////////////////////////////
/// Create an object with the specified name and descriptive title.
/// The newly created object has no clients or servers and has its
/// dirty flags set.

RooAbsArg::RooAbsArg(const char *name, const char *title)
   : TNamed(name, title), _deleteWatch(kFALSE), _valueDirty(kTRUE), _shapeDirty(kTRUE), _operMode(Auto), _fast(kFALSE),
     _ownedComponents(0), _prohibitServerRedirect(kFALSE), _namePtr(0), _isConstant(kFALSE),
     _localNoInhibitDirty(kFALSE), _myws(0)
{
  if (name == nullptr || strlen(name) == 0) {
    throw std::logic_error("Each RooFit object needs a name. "
        "Objects representing the same entity (e.g. an observable 'x') are identified using their name.");
  }
  _namePtr = (TNamed*) RooNameReg::instance().constPtr(GetName()) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor transfers all boolean and string properties of the original
/// object. Transient properties and client-server links are not copied

RooAbsArg::RooAbsArg(const RooAbsArg &other, const char *name)
   : TNamed(other.GetName(), other.GetTitle()), RooPrintable(other), _boolAttrib(other._boolAttrib),
     _stringAttrib(other._stringAttrib), _deleteWatch(other._deleteWatch), _operMode(Auto), _fast(kFALSE),
     _ownedComponents(0), _prohibitServerRedirect(kFALSE), _namePtr(other._namePtr),
     _isConstant(other._isConstant), _localNoInhibitDirty(other._localNoInhibitDirty), _myws(0)
{
  // Use name in argument, if supplied
  if (name) {
    TNamed::SetName(name) ;
    _namePtr = (TNamed*) RooNameReg::instance().constPtr(name) ;
  } else {
    // Same name, don't recalculate name pointer (expensive)
    TNamed::SetName(other.GetName()) ;
    _namePtr = other._namePtr ;
  }

  // Copy server list by hand
  Bool_t valueProp, shapeProp ;
  for (const auto server : other._serverList) {
    valueProp = server->_clientListValue.containsByNamePtr(&other);
    shapeProp = server->_clientListShape.containsByNamePtr(&other);
    addServer(*server,valueProp,shapeProp) ;
  }

  setValueDirty() ;
  setShapeDirty() ;

  //setAttribute(Form("CloneOf(%08x)",&other)) ;
  //cout << "RooAbsArg::cctor(" << this << ") #bools = " << _boolAttrib.size() << " #strings = " << _stringAttrib.size() << endl ;

}

////////////////////////////////////////////////////////////////////////////////
/// Assign all boolean and string properties of the original
/// object. Transient properties and client-server links are not assigned.
RooAbsArg& RooAbsArg::operator=(const RooAbsArg& other) {
  TNamed::operator=(other);
  RooPrintable::operator=(other);
  _boolAttrib = other._boolAttrib;
  _stringAttrib = other._stringAttrib;
  _deleteWatch = other._deleteWatch;
  _operMode = other._operMode;
  _fast = other._fast;
  _ownedComponents = nullptr;
  _prohibitServerRedirect = other._prohibitServerRedirect;
  _namePtr = other._namePtr;
  _isConstant = other._isConstant;
  _localNoInhibitDirty = other._localNoInhibitDirty;
  _myws = nullptr;

  bool valueProp, shapeProp;
  for (const auto server : other._serverList) {
    valueProp = server->_clientListValue.containsByNamePtr(&other);
    shapeProp = server->_clientListShape.containsByNamePtr(&other);
    addServer(*server,valueProp,shapeProp) ;
  }

  setValueDirty();
  setShapeDirty();

  return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

RooAbsArg::~RooAbsArg()
{
  // Notify all servers that they no longer need to serve us
  while (!_serverList.empty()) {
    removeServer(*_serverList.containedObjects().back(), kTRUE);
  }

  // Notify all clients that they are in limbo
  std::vector<RooAbsArg*> clientListTmp(_clientList.begin(), _clientList.end()); // have to copy, as we invalidate iterators
  Bool_t first(kTRUE) ;
  for (auto client : clientListTmp) {
    client->setAttribute("ServerDied") ;
    TString attr("ServerDied:");
    attr.Append(GetName());
    attr.Append(Form("(%lx)",(ULong_t)this)) ;
    client->setAttribute(attr.Data());
    client->removeServer(*this,kTRUE);

    if (_verboseDirty) {

      if (first) {
	cxcoutD(Tracing) << "RooAbsArg::dtor(" << GetName() << "," << this << ") DeleteWatch: object is being destroyed" << endl ;
	first = kFALSE ;
      }

      cxcoutD(Tracing)  << fName << "::" << ClassName() << ":~RooAbsArg: dependent \""
		       << client->GetName() << "\" should have been deleted first" << endl ;
    }
  }

  if (_ownedComponents) {
    delete _ownedComponents ;
    _ownedComponents = 0 ;
  }

}


////////////////////////////////////////////////////////////////////////////////
/// Control global dirty inhibit mode. When set to true no value or shape dirty
/// flags are propagated and cache is always considered to be dirty.

void RooAbsArg::setDirtyInhibit(Bool_t flag)
{
  _inhibitDirty = flag ;
}


////////////////////////////////////////////////////////////////////////////////
/// Activate verbose messaging related to dirty flag propagation

void RooAbsArg::verboseDirty(Bool_t flag)
{
  _verboseDirty = flag ;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if this object was created as a clone of 'other'

Bool_t RooAbsArg::isCloneOf(const RooAbsArg& other) const
{
  return (getAttribute(Form("CloneOf(%lx)",(ULong_t)&other)) ||
	  other.getAttribute(Form("CloneOf(%lx)",(ULong_t)this))) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Set (default) or clear a named boolean attribute of this object.

void RooAbsArg::setAttribute(const Text_t* name, Bool_t value)
{
  // Preserve backward compatibility - any strong
  if(string("Constant")==name) {
    _isConstant = value ;
  }

  if (value) {
    _boolAttrib.insert(name) ;
  } else {
    set<string>::iterator iter = _boolAttrib.find(name) ;
    if (iter != _boolAttrib.end()) {
      _boolAttrib.erase(iter) ;
    }

  }

}


////////////////////////////////////////////////////////////////////////////////
/// Check if a named attribute is set. By default, all attributes are unset.

Bool_t RooAbsArg::getAttribute(const Text_t* name) const
{
  return (_boolAttrib.find(name) != _boolAttrib.end()) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Associate string 'value' to this object under key 'key'

void RooAbsArg::setStringAttribute(const Text_t* key, const Text_t* value)
{
  if (value) {
    _stringAttrib[key] = value ;
  } else {
    _stringAttrib.erase(key) ;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Get string attribute mapped under key 'key'. Returns null pointer
/// if no attribute exists under that key

const Text_t* RooAbsArg::getStringAttribute(const Text_t* key) const
{
  map<string,string>::const_iterator iter = _stringAttrib.find(key) ;
  if (iter!=_stringAttrib.end()) {
    return iter->second.c_str() ;
  } else {
    return 0 ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Set (default) or clear a named boolean attribute of this object.

void RooAbsArg::setTransientAttribute(const Text_t* name, Bool_t value)
{
  if (value) {

    _boolAttribTransient.insert(name) ;

  } else {

    set<string>::iterator iter = _boolAttribTransient.find(name) ;
    if (iter != _boolAttribTransient.end()) {
      _boolAttribTransient.erase(iter) ;
    }

  }

}


////////////////////////////////////////////////////////////////////////////////
/// Check if a named attribute is set. By default, all attributes
/// are unset.

Bool_t RooAbsArg::getTransientAttribute(const Text_t* name) const
{
  return (_boolAttribTransient.find(name) != _boolAttribTransient.end()) ;
}




////////////////////////////////////////////////////////////////////////////////
/// Register another RooAbsArg as a server to us, ie, declare that
/// we depend on it.
/// \param server The server to be registered.
/// \param valueProp In addition to the basic client-server relationship, declare dependence on the server's value.
/// \param valueProp In addition to the basic client-server relationship, declare dependence on the server's shape.
/// \param refCount Optionally add with higher reference count (if multiple components depend on it)

void RooAbsArg::addServer(RooAbsArg& server, Bool_t valueProp, Bool_t shapeProp, std::size_t refCount)
{
  if (_prohibitServerRedirect) {
    cxcoutF(LinkStateMgmt) << "RooAbsArg::addServer(" << this << "," << GetName()
			   << "): PROHIBITED SERVER ADDITION REQUESTED: adding server " << server.GetName()
			   << "(" << &server << ") for " << (valueProp?"value ":"") << (shapeProp?"shape":"") << endl ;
    throw std::logic_error("PROHIBITED SERVER ADDITION REQUESTED in RooAbsArg::addServer");
  }

  cxcoutD(LinkStateMgmt) << "RooAbsArg::addServer(" << this << "," << GetName() << "): adding server " << server.GetName()
			 << "(" << &server << ") for " << (valueProp?"value ":"") << (shapeProp?"shape":"") << endl ;

  if (server.operMode()==ADirty && operMode()!=ADirty && valueProp) {
    setOperMode(ADirty) ;
  }


  // LM: use hash tables for larger lists
//  if (_serverList.GetSize() > 999 && _serverList.getHashTableSize() == 0) _serverList.setHashTableSize(1000);
//  if (server._clientList.GetSize() > 999 && server._clientList.getHashTableSize() == 0) server._clientList.setHashTableSize(1000);
//  if (server._clientListValue.GetSize() >  999 && server._clientListValue.getHashTableSize() == 0) server._clientListValue.setHashTableSize(1000);

  // Add server link to given server
  _serverList.Add(&server, refCount) ;

  server._clientList.Add(this, refCount);
  if (valueProp) server._clientListValue.Add(this, refCount);
  if (shapeProp) server._clientListShape.Add(this, refCount);
}



////////////////////////////////////////////////////////////////////////////////
/// Register a list of RooAbsArg as servers to us by calling
/// addServer() for each arg in the list

void RooAbsArg::addServerList(RooAbsCollection& serverList, Bool_t valueProp, Bool_t shapeProp)
{
  _serverList.reserve(_serverList.size() + serverList.size());

  for (const auto arg : serverList) {
    addServer(*arg,valueProp,shapeProp) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Unregister another RooAbsArg as a server to us, ie, declare that
/// we no longer depend on its value and shape.

void RooAbsArg::removeServer(RooAbsArg& server, Bool_t force)
{
  if (_prohibitServerRedirect) {
    cxcoutF(LinkStateMgmt) << "RooAbsArg::addServer(" << this << "," << GetName() << "): PROHIBITED SERVER REMOVAL REQUESTED: removing server "
			   << server.GetName() << "(" << &server << ")" << endl ;
    assert(0) ;
  }

  if (_verboseDirty) {
    cxcoutD(LinkStateMgmt) << "RooAbsArg::removeServer(" << GetName() << "): removing server "
			   << server.GetName() << "(" << &server << ")" << endl ;
  }

  // Remove server link to given server
  _serverList.Remove(&server, force) ;

  server._clientList.Remove(this, force) ;
  server._clientListValue.Remove(this, force) ;
  server._clientListShape.Remove(this, force) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Replace 'oldServer' with 'newServer'

void RooAbsArg::replaceServer(RooAbsArg& oldServer, RooAbsArg& newServer, Bool_t propValue, Bool_t propShape)
{
  Int_t count = _serverList.refCount(&oldServer);
  removeServer(oldServer, kTRUE);
  addServer(newServer, propValue, propShape, count);
}


////////////////////////////////////////////////////////////////////////////////
/// Change dirty flag propagation mask for specified server

void RooAbsArg::changeServer(RooAbsArg& server, Bool_t valueProp, Bool_t shapeProp)
{
  if (!_serverList.containsByNamePtr(&server)) {
    coutE(LinkStateMgmt) << "RooAbsArg::changeServer(" << GetName() << "): Server "
	 << server.GetName() << " not registered" << endl ;
    return ;
  }

  // This condition should not happen, but check anyway
  if (!server._clientList.containsByNamePtr(this)) {
    coutE(LinkStateMgmt) << "RooAbsArg::changeServer(" << GetName() << "): Server "
			 << server.GetName() << " doesn't have us registered as client" << endl ;
    return ;
  }

  // Remove all propagation links, then reinstall requested ones ;
  Int_t vcount = server._clientListValue.refCount(this) ;
  Int_t scount = server._clientListShape.refCount(this) ;
  server._clientListValue.RemoveAll(this) ;
  server._clientListShape.RemoveAll(this) ;
  if (valueProp) {
    server._clientListValue.Add(this, vcount) ;
  }
  if (shapeProp) {
    server._clientListShape.Add(this, scount) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Fill supplied list with all leaf nodes of the arg tree, starting with
/// ourself as top node. A leaf node is node that has no servers declared.

void RooAbsArg::leafNodeServerList(RooAbsCollection* list, const RooAbsArg* arg, Bool_t recurseNonDerived) const
{
  treeNodeServerList(list,arg,kFALSE,kTRUE,kFALSE,recurseNonDerived) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Fill supplied list with all branch nodes of the arg tree starting with
/// ourself as top node. A branch node is node that has one or more servers declared.

void RooAbsArg::branchNodeServerList(RooAbsCollection* list, const RooAbsArg* arg, Bool_t recurseNonDerived) const
{
  treeNodeServerList(list,arg,kTRUE,kFALSE,kFALSE,recurseNonDerived) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Fill supplied list with nodes of the arg tree, following all server links,
/// starting with ourself as top node.
/// \param[in] list Output list
/// \param[in] arg Start searching at this element of the tree.
/// \param[in] doBranch Add branch nodes to the list.
/// \param[in] doLeaf Add leaf nodes to the list.
/// \param[in] valueOnly Only check if an element is a value server (no shape server).
/// \param[in] recurseFundamental

void RooAbsArg::treeNodeServerList(RooAbsCollection* list, const RooAbsArg* arg, Bool_t doBranch, Bool_t doLeaf, Bool_t valueOnly, Bool_t recurseFundamental) const
{
//   if (arg==0) {
//     cout << "treeNodeServerList(" << GetName() << ") doBranch=" << (doBranch?"T":"F") << " doLeaf = " << (doLeaf?"T":"F") << " valueOnly=" << (valueOnly?"T":"F") << endl ;
//   }

  if (!arg) {
    list->reserve(10);
    arg=this ;
  }

  // Decide if to add current node
  if ((doBranch&&doLeaf) ||
      (doBranch&&arg->isDerived()) ||
      (doLeaf&&arg->isFundamental()&&(!(recurseFundamental&&arg->isDerived()))) ||
      (doLeaf && !arg->isFundamental() && !arg->isDerived())) {

    list->add(*arg,kTRUE) ;
  }

  // Recurse if current node is derived
  if (arg->isDerived() && (!arg->isFundamental() || recurseFundamental)) {
    for (const auto server : arg->_serverList) {

      // Skip non-value server nodes if requested
      Bool_t isValueSrv = server->_clientListValue.containsByNamePtr(arg);
      if (valueOnly && !isValueSrv) {
        continue ;
      }
      treeNodeServerList(list,server,doBranch,doLeaf,valueOnly,recurseFundamental) ;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Create a list of leaf nodes in the arg tree starting with
/// ourself as top node that don't match any of the names of the variable list
/// of the supplied data set (the dependents). The caller of this
/// function is responsible for deleting the returned argset.
/// The complement of this function is getObservables()

RooArgSet* RooAbsArg::getParameters(const RooAbsData* set, Bool_t stripDisconnected) const
{
  return getParameters(set?set->get():0,stripDisconnected) ;
}


////////////////////////////////////////////////////////////////////////////////
/// INTERNAL helper function for getParameters()

void RooAbsArg::addParameters(RooArgSet& params, const RooArgSet* nset,Bool_t stripDisconnected) const
{
  // INTERNAL helper function for getParameters()

  RooArgSet nodeParamServers ;
  RooArgSet nodeBranchServers ;
  for (const auto server : _serverList) {
    if (server->isValueServer(*this)) {
      if (server->isFundamental()) {
        if (!nset || !server->dependsOn(*nset)) {
          nodeParamServers.add(*server) ;
        }
      } else {
        nodeBranchServers.add(*server) ;
      }
    }
  }

  // Allow pdf to strip parameters from list before adding it
  getParametersHook(nset,&nodeParamServers,stripDisconnected) ;

  // Add parameters of this node to the combined list
  params.add(nodeParamServers,kTRUE) ;

  // Now recurse into branch servers
  for (const auto server : nodeBranchServers) {
    server->addParameters(params,nset) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Create a list of leaf nodes in the arg tree starting with
/// ourself as top node that don't match any of the names the args in the
/// supplied argset. The caller of this function is responsible
/// for deleting the returned argset. The complement of this function
/// is getObservables()

RooArgSet* RooAbsArg::getParameters(const RooArgSet* nset, Bool_t stripDisconnected) const
{

   // Check for cached parameter set
   if (_myws) {
      RooNameSet nsetObs(nset ? *nset : RooArgSet());
      const RooArgSet *paramSet = _myws->set(Form("CACHE_PARAMS_OF_PDF_%s_FOR_OBS_%s", GetName(), nsetObs.content()));
      if (paramSet) {
         // cout << " restoring parameter cache from workspace for pdf " << IsA()->GetName() << "::" << GetName() <<
         // endl ;
         return new RooArgSet(*paramSet);
      }
   }

   RooArgSet *parList = new RooArgSet("parameters");

   addParameters(*parList, nset, stripDisconnected);

   parList->sort();

   // Cache parameter set
   if (_myws && parList->getSize() > 10) {
      RooNameSet nsetObs(nset ? *nset : RooArgSet());
      _myws->defineSetInternal(Form("CACHE_PARAMS_OF_PDF_%s_FOR_OBS_%s", GetName(), nsetObs.content()), *parList);
      // cout << " caching parameters in workspace for pdf " << IsA()->GetName() << "::" << GetName() << endl ;
   }

   return parList;
}



////////////////////////////////////////////////////////////////////////////////
/// Create a list of leaf nodes in the arg tree starting with
/// ourself as top node that match any of the names of the variable list
/// of the supplied data set (the dependents). The caller of this
/// function is responsible for deleting the returned argset.
/// The complement of this function is getParameters().

RooArgSet* RooAbsArg::getObservables(const RooAbsData* set) const
{
  if (!set) return new RooArgSet ;

  return getObservables(set->get()) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a list of leaf nodes in the arg tree starting with
/// ourself as top node that match any of the names the args in the
/// supplied argset. The caller of this function is responsible
/// for deleting the returned argset. The complement of this function
/// is getParameters().

RooArgSet* RooAbsArg::getObservables(const RooArgSet* dataList, Bool_t valueOnly) const
{
  //cout << "RooAbsArg::getObservables(" << GetName() << ")" << endl ;

  RooArgSet* depList = new RooArgSet("dependents") ;
  if (!dataList) return depList ;

  // Make iterator over tree leaf node list
  RooArgSet leafList("leafNodeServerList") ;
  treeNodeServerList(&leafList,0,kFALSE,kTRUE,valueOnly) ;

  if (valueOnly) {
    for (const auto arg : leafList) {
      if (arg->dependsOnValue(*dataList) && arg->isLValue()) {
        depList->add(*arg) ;
      }
    }
  } else {
    for (const auto arg : leafList) {
      if (arg->dependsOn(*dataList) && arg->isLValue()) {
        depList->add(*arg) ;
      }
    }
  }

  return depList ;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a RooArgSet with all components (branch nodes) of the
/// expression tree headed by this object.
RooArgSet* RooAbsArg::getComponents() const
{
  TString name(GetName()) ;
  name.Append("_components") ;

  RooArgSet* set = new RooArgSet(name) ;
  branchNodeServerList(set) ;

  return set ;
}



////////////////////////////////////////////////////////////////////////////////
/// Overloadable function in which derived classes can implement
/// consistency checks of the variables. If this function returns
/// true, indicating an error, the fitter or generator will abort.

Bool_t RooAbsArg::checkObservables(const RooArgSet*) const
{
  return kFALSE ;
}


////////////////////////////////////////////////////////////////////////////////
/// Recursively call checkObservables on all nodes in the expression tree

Bool_t RooAbsArg::recursiveCheckObservables(const RooArgSet* nset) const
{
  RooArgSet nodeList ;
  treeNodeServerList(&nodeList) ;
  RooFIter iter = nodeList.fwdIterator() ;

  RooAbsArg* arg ;
  Bool_t ret(kFALSE) ;
  while((arg=iter.next())) {
    if (arg->getAttribute("ServerDied")) {
      coutE(LinkStateMgmt) << "RooAbsArg::recursiveCheckObservables(" << GetName() << "): ERROR: one or more servers of node "
			   << arg->GetName() << " no longer exists!" << endl ;
      arg->Print("v") ;
      ret = kTRUE ;
    }
    ret |= arg->checkObservables(nset) ;
  }

  return ret ;
}


////////////////////////////////////////////////////////////////////////////////
/// Test whether we depend on (ie, are served by) any object in the
/// specified collection. Uses the dependsOn(RooAbsArg&) member function.

Bool_t RooAbsArg::dependsOn(const RooAbsCollection& serverList, const RooAbsArg* ignoreArg, Bool_t valueOnly) const
{
  // Test whether we depend on (ie, are served by) any object in the
  // specified collection. Uses the dependsOn(RooAbsArg&) member function.

  for (auto server : serverList) {
    if (dependsOn(*server,ignoreArg,valueOnly)) {
      return kTRUE;
    }
  }
  return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Test whether we depend on (ie, are served by) the specified object.
/// Note that RooAbsArg objects are considered equivalent if they have
/// the same name.

Bool_t RooAbsArg::dependsOn(const RooAbsArg& testArg, const RooAbsArg* ignoreArg, Bool_t valueOnly) const
{
  if (this==ignoreArg) return kFALSE ;

  // First check if testArg is self
  //if (!TString(testArg.GetName()).CompareTo(GetName())) return kTRUE ;
  if (testArg.namePtr()==namePtr()) return kTRUE ;


  // Next test direct dependence
  RooAbsArg* foundServer = findServer(testArg) ;
  if (foundServer) {

    // Return true if valueOnly is FALSE or if server is value server, otherwise keep looking
    if ( !valueOnly || foundServer->isValueServer(*this)) {
      return kTRUE ;
    }
  }

  // If not, recurse
  for (const auto server : _serverList) {
    if ( !valueOnly || server->isValueServer(*this)) {
      if (server->dependsOn(testArg,ignoreArg,valueOnly)) {
        return kTRUE ;
      }
    }
  }

  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Test if any of the nodes of tree are shared with that of the given tree

Bool_t RooAbsArg::overlaps(const RooAbsArg& testArg, Bool_t valueOnly) const
{
  RooArgSet list("treeNodeList") ;
  treeNodeServerList(&list) ;

  return valueOnly ? testArg.dependsOnValue(list) : testArg.dependsOn(list) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Test if any of the dependents of the arg tree (as determined by getObservables)
/// overlaps with those of the testArg.

Bool_t RooAbsArg::observableOverlaps(const RooAbsData* dset, const RooAbsArg& testArg) const
{
  return observableOverlaps(dset->get(),testArg) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Test if any of the dependents of the arg tree (as determined by getObservables)
/// overlaps with those of the testArg.

Bool_t RooAbsArg::observableOverlaps(const RooArgSet* nset, const RooAbsArg& testArg) const
{
  RooArgSet* depList = getObservables(nset) ;
  Bool_t ret = testArg.dependsOn(*depList) ;
  delete depList ;
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Mark this object as having changed its value, and propagate this status
/// change to all of our clients. If the object is not in automatic dirty
/// state propagation mode, this call has no effect.

void RooAbsArg::setValueDirty(const RooAbsArg* source)
{
  _allBatchesDirty = true;

  if (_operMode!=Auto || _inhibitDirty) return ;

  // Handle no-propagation scenarios first
  if (_clientListValue.size() == 0) {
    _valueDirty = kTRUE ;
    return ;
  }

  // Cyclical dependency interception
  if (source==0) {
    source=this ;
  } else if (source==this) {
    // Cyclical dependency, abort
    coutE(LinkStateMgmt) << "RooAbsArg::setValueDirty(" << GetName()
			 << "): cyclical dependency detected, source = " << source->GetName() << endl ;
    //assert(0) ;
    return ;
  }

  // Propagate dirty flag to all clients if this is a down->up transition
  if (_verboseDirty) {
    cxcoutD(LinkStateMgmt) << "RooAbsArg::setValueDirty(" << (source?source->GetName():"self") << "->" << GetName() << "," << this
			   << "): dirty flag " << (_valueDirty?"already ":"") << "raised" << endl ;
  }

  _valueDirty = kTRUE ;


  for (auto client : _clientListValue) {
    client->setValueDirty(source) ;
  }


}


////////////////////////////////////////////////////////////////////////////////
/// Mark this object as having changed its shape, and propagate this status
/// change to all of our clients.

void RooAbsArg::setShapeDirty(const RooAbsArg* source)
{
  if (_verboseDirty) {
    cxcoutD(LinkStateMgmt) << "RooAbsArg::setShapeDirty(" << GetName()
			   << "): dirty flag " << (_shapeDirty?"already ":"") << "raised" << endl ;
  }

  if (_clientListShape.empty()) {
    _shapeDirty = kTRUE ;
    return ;
  }

  // Set 'dirty' shape state for this object and propagate flag to all its clients
  if (source==0) {
    source=this ;
  } else if (source==this) {
    // Cyclical dependency, abort
    coutE(LinkStateMgmt) << "RooAbsArg::setShapeDirty(" << GetName()
	 << "): cyclical dependency detected" << endl ;
    return ;
  }

  // Propagate dirty flag to all clients if this is a down->up transition
  _shapeDirty=kTRUE ;

  for (auto client : _clientListShape) {
    client->setShapeDirty(source) ;
    client->setValueDirty(source) ;
  }

}



////////////////////////////////////////////////////////////////////////////////
/// Replace all direct servers of this object with the new servers in `newServerList`.
/// This substitutes objects that we receive values from with new objects that have the same name.
/// \see recursiveRedirectServers() Use recursive version if servers that are only indirectly serving this object should be replaced as well.
/// \see redirectServers() If only the direct servers of an object need to be replaced.
///
/// Note that changing the types of objects is generally allowed, but can be wrong if the interface of an object changes.
/// For example, one can reparametrise a model by substituting a variable with a function:
/// \f[
///   f(x\, |\, a) = a \cdot x \rightarrow f(x\, |\, b) = (2.1 \cdot b) \cdot x
/// \f]
/// If an object, however, expects a PDF, and this is substituted with a function that isn't normalised, wrong results might be obtained
/// or it might even crash the program. The types of the objects being substituted are not checked.
///
/// \param[in] newSetOrig Set of new servers that should be used instead of the current servers.
/// \param[in] mustReplaceAll A warning is printed and error status is returned if not all servers could be
/// substituted successfully.
/// \param[in] nameChange If false, an object named "x" is replaced with an object named "x" in `newSetOrig`.
/// If the object in `newSet` is called differently, set `nameChange` to true and use setStringAttribute on the x object:
/// ```
/// objectToReplaceX.setStringAttribute("ORIGNAME", "x")
/// ```
/// \param[in] isRecursionStep Internal switch used when called from recursiveRedirectServers().
Bool_t RooAbsArg::redirectServers(const RooAbsCollection& newSetOrig, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursionStep)
{
  // Trivial case, no servers
  if (_serverList.empty()) return kFALSE ;
  if (newSetOrig.getSize()==0) return kFALSE ;

  // Strip any non-matching removal nodes from newSetOrig
  RooAbsCollection* newSet ;

  if (nameChange) {
    newSet = new RooArgSet ;
    for (auto arg : newSetOrig) {

      if (string("REMOVAL_DUMMY")==arg->GetName()) {

        if (arg->getAttribute("REMOVE_ALL")) {
          newSet->add(*arg) ;
        } else if (arg->getAttribute(Form("REMOVE_FROM_%s",getStringAttribute("ORIGNAME")))) {
          newSet->add(*arg) ;
        }
      } else {
        newSet->add(*arg) ;
      }
    }
  } else {
    newSet = (RooAbsCollection*) &newSetOrig ;
  }

  // Replace current servers with new servers with the same name from the given list
  Bool_t ret(kFALSE) ;

  //Copy original server list to not confuse the iterator while deleting
  std::vector<RooAbsArg*> origServerList, origServerValue, origServerShape;
  auto origSize = _serverList.size();
  origServerList.reserve(origSize);
  origServerValue.reserve(origSize);

  for (const auto oldServer : _serverList) {
    origServerList.push_back(oldServer) ;

    // Retrieve server side link state information
    if (oldServer->_clientListValue.containsByNamePtr(this)) {
      origServerValue.push_back(oldServer) ;
    }
    if (oldServer->_clientListShape.containsByNamePtr(this)) {
      origServerShape.push_back(oldServer) ;
    }
  }

  // Delete all previously registered servers
  for (auto oldServer : origServerList) {

    RooAbsArg * newServer= oldServer->findNewServer(*newSet, nameChange);

    if (newServer && _verboseDirty) {
      cxcoutD(LinkStateMgmt) << "RooAbsArg::redirectServers(" << (void*)this << "," << GetName() << "): server " << oldServer->GetName()
			         << " redirected from " << oldServer << " to " << newServer << endl ;
    }

    if (!newServer) {
      if (mustReplaceAll) {
        coutE(LinkStateMgmt) << "RooAbsArg::redirectServers(" << (void*)this << "," << GetName() << "): server " << oldServer->GetName()
			           << " (" << (void*)oldServer << ") not redirected" << (nameChange?"[nameChange]":"") << endl ;
        ret = kTRUE ;
      }
      continue ;
    }

    auto findByNamePtr = [&oldServer](const RooAbsArg * item) {
      return oldServer->namePtr() == item->namePtr();
    };
    bool propValue = std::any_of(origServerValue.begin(), origServerValue.end(), findByNamePtr);
    bool propShape = std::any_of(origServerShape.begin(), origServerShape.end(), findByNamePtr);

    if (newServer != this) {
      replaceServer(*oldServer,*newServer,propValue,propShape) ;
    }
  }


  setValueDirty() ;
  setShapeDirty() ;

  // Process the proxies
  for (int i=0 ; i<numProxies() ; i++) {
    RooAbsProxy* p = getProxy(i) ;
    if (!p) continue ;
    Bool_t ret2 = p->changePointer(*newSet,nameChange,kFALSE) ;

    if (mustReplaceAll && !ret2) {
      auto ap = dynamic_cast<const RooArgProxy*>(p);
      coutE(LinkStateMgmt) << "RooAbsArg::redirectServers(" << GetName()
              << "): ERROR, proxy '" << p->name()
              << "' with arg '" << (ap ? ap->absArg()->GetName() : "<could not cast>") << "' could not be adjusted" << endl;
      ret = kTRUE ;
    }
  }


  // Optional subclass post-processing
  for (Int_t i=0 ;i<numCaches() ; i++) {
    ret |= getCache(i)->redirectServersHook(*newSet,mustReplaceAll,nameChange,isRecursionStep) ;
  }
  ret |= redirectServersHook(*newSet,mustReplaceAll,nameChange,isRecursionStep) ;

  if (nameChange) {
    delete newSet ;
  }

  return ret ;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the new server in the specified set that matches the old server.
/// Allow a name change if nameChange is kTRUE, in which case the new
/// server is selected by searching for a new server with an attribute
/// of "ORIGNAME:<oldName>". Return zero if there is not a unique match.

RooAbsArg *RooAbsArg::findNewServer(const RooAbsCollection &newSet, Bool_t nameChange) const
{
  RooAbsArg *newServer = 0;
  if (!nameChange) {
    newServer = newSet.find(*this) ;
  }
  else {
    // Name changing server redirect:
    // use 'ORIGNAME:<oldName>' attribute instead of name of new server
    TString nameAttrib("ORIGNAME:") ;
    nameAttrib.Append(GetName()) ;

    RooArgSet* tmp = (RooArgSet*) newSet.selectByAttrib(nameAttrib,kTRUE) ;
    if(0 != tmp) {

      // Check if any match was found
      if (tmp->getSize()==0) {
        delete tmp ;
        return 0 ;
      }

      // Check if match is unique
      if(tmp->getSize()>1) {
        coutF(LinkStateMgmt) << "RooAbsArg::redirectServers(" << GetName() << "): FATAL Error, " << tmp->getSize() << " servers with "
            << nameAttrib << " attribute" << endl ;
        tmp->Print("v") ;
        assert(0) ;
      }

      // use the unique element in the set
      newServer= tmp->first();
      delete tmp ;
    }
  }
  return newServer;
}


////////////////////////////////////////////////////////////////////////////////
/// Recursively replace all servers with the new servers in `newSet`.
/// This substitutes objects that we receive values from (also indirectly through other objects) with new objects that have the same name.
///
/// *Copied from redirectServers:*
///
/// \copydetails RooAbsArg::redirectServers
Bool_t RooAbsArg::recursiveRedirectServers(const RooAbsCollection& newSet, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t recurseInNewSet)
{
  // Cyclic recursion protection
  static std::set<const RooAbsArg*> callStack;
  {
    std::set<const RooAbsArg*>::iterator it = callStack.lower_bound(this);
    if (it != callStack.end() && this == *it) {
      return kFALSE;
    } else {
      callStack.insert(it, this);
    }
  }

  // Do not recurse into newset if not so specified
//   if (!recurseInNewSet && newSet.contains(*this)) {
//     return kFALSE ;
//   }


  // Apply the redirectServers function recursively on all branch nodes in this argument tree.
  Bool_t ret(kFALSE) ;

  cxcoutD(LinkStateMgmt) << "RooAbsArg::recursiveRedirectServers(" << this << "," << GetName() << ") newSet = " << newSet << " mustReplaceAll = "
			 << (mustReplaceAll?"T":"F") << " nameChange = " << (nameChange?"T":"F") << " recurseInNewSet = " << (recurseInNewSet?"T":"F") << endl ;

  // Do redirect on self (identify operation as recursion step)
  ret |= redirectServers(newSet,mustReplaceAll,nameChange,kTRUE) ;

  // Do redirect on servers
  for (const auto server : _serverList){
    ret |= server->recursiveRedirectServers(newSet,mustReplaceAll,nameChange,recurseInNewSet) ;
  }

  callStack.erase(this);
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Register an RooArgProxy in the proxy list. This function is called by owned
/// proxies upon creation. After registration, this arg wil forward pointer
/// changes from serverRedirects and updates in cached normalization sets
/// to the proxies immediately after they occur. The proxied argument is
/// also added as value and/or shape server

void RooAbsArg::registerProxy(RooArgProxy& proxy)
{
  // Every proxy can be registered only once
  if (_proxyList.FindObject(&proxy)) {
    coutE(LinkStateMgmt) << "RooAbsArg::registerProxy(" << GetName() << "): proxy named "
			 << proxy.GetName() << " for arg " << proxy.absArg()->GetName()
			 << " already registered" << endl ;
    return ;
  }

//   cout << (void*)this << " " << GetName() << ": registering proxy "
//        << (void*)&proxy << " with name " << proxy.name() << " in mode "
//        << (proxy.isValueServer()?"V":"-") << (proxy.isShapeServer()?"S":"-") << endl ;

  // Register proxied object as server
  if (proxy.absArg()) {
    addServer(*proxy.absArg(),proxy.isValueServer(),proxy.isShapeServer()) ;
  }

  // Register proxy itself
  _proxyList.Add(&proxy) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Remove proxy from proxy list. This functions is called by owned proxies
/// upon their destruction.

void RooAbsArg::unRegisterProxy(RooArgProxy& proxy)
{
  _proxyList.Remove(&proxy) ;
  _proxyList.Compress() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Register an RooSetProxy in the proxy list. This function is called by owned
/// proxies upon creation. After registration, this arg wil forward pointer
/// changes from serverRedirects and updates in cached normalization sets
/// to the proxies immediately after they occur.

void RooAbsArg::registerProxy(RooSetProxy& proxy)
{
  // Every proxy can be registered only once
  if (_proxyList.FindObject(&proxy)) {
    coutE(LinkStateMgmt) << "RooAbsArg::registerProxy(" << GetName() << "): proxy named "
			 << proxy.GetName() << " already registered" << endl ;
    return ;
  }

  // Register proxy itself
  _proxyList.Add(&proxy) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Remove proxy from proxy list. This functions is called by owned proxies
/// upon their destruction.

void RooAbsArg::unRegisterProxy(RooSetProxy& proxy)
{
  _proxyList.Remove(&proxy) ;
  _proxyList.Compress() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Register an RooListProxy in the proxy list. This function is called by owned
/// proxies upon creation. After registration, this arg wil forward pointer
/// changes from serverRedirects and updates in cached normalization sets
/// to the proxies immediately after they occur.

void RooAbsArg::registerProxy(RooListProxy& proxy)
{
  // Every proxy can be registered only once
  if (_proxyList.FindObject(&proxy)) {
    coutE(LinkStateMgmt) << "RooAbsArg::registerProxy(" << GetName() << "): proxy named "
			 << proxy.GetName() << " already registered" << endl ;
    return ;
  }

  // Register proxy itself
  Int_t nProxyOld = _proxyList.GetEntries() ;
  _proxyList.Add(&proxy) ;
  if (_proxyList.GetEntries()!=nProxyOld+1) {
    cout << "RooAbsArg::registerProxy(" << GetName() << ") proxy registration failure! nold=" << nProxyOld << " nnew=" << _proxyList.GetEntries() << endl ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Remove proxy from proxy list. This functions is called by owned proxies
/// upon their destruction.

void RooAbsArg::unRegisterProxy(RooListProxy& proxy)
{
  _proxyList.Remove(&proxy) ;
  _proxyList.Compress() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the nth proxy from the proxy list.

RooAbsProxy* RooAbsArg::getProxy(Int_t index) const
{
  // Cross cast: proxy list returns TObject base pointer, we need
  // a RooAbsProxy base pointer. C++ standard requires
  // a dynamic_cast for this.
  return dynamic_cast<RooAbsProxy*> (_proxyList.At(index)) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the number of registered proxies.

Int_t RooAbsArg::numProxies() const
{
   return _proxyList.GetEntriesFast();
}



////////////////////////////////////////////////////////////////////////////////
/// Forward a change in the cached normalization argset
/// to all the registered proxies.

void RooAbsArg::setProxyNormSet(const RooArgSet* nset)
{
  for (int i=0 ; i<numProxies() ; i++) {
    RooAbsProxy* p = getProxy(i) ;
    if (!p) continue ;
    getProxy(i)->changeNormSet(nset) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Overloadable function for derived classes to implement
/// attachment as branch to a TTree

void RooAbsArg::attachToTree(TTree& ,Int_t)
{
  coutE(Contents) << "RooAbsArg::attachToTree(" << GetName()
		  << "): Cannot be attached to a TTree" << endl ;
}



////////////////////////////////////////////////////////////////////////////////
/// WVE (08/21/01) Probably obsolete now

Bool_t RooAbsArg::isValid() const
{
  return kTRUE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Print object name

void RooAbsArg::printName(ostream& os) const
{
  os << GetName() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print object title

void RooAbsArg::printTitle(ostream& os) const
{
  os << GetTitle() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print object class name

void RooAbsArg::printClassName(ostream& os) const
{
  os << IsA()->GetName() ;
}


void RooAbsArg::printAddress(ostream& os) const
{
  // Print addrss of this RooAbsArg
  os << this ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print object arguments, ie its proxies

void RooAbsArg::printArgs(ostream& os) const
{
  // Print nothing if there are no dependencies
  if (numProxies()==0) return ;

  os << "[ " ;
  for (Int_t i=0 ; i<numProxies() ; i++) {
    RooAbsProxy* p = getProxy(i) ;
    if (p==0) continue ;
    if (!TString(p->name()).BeginsWith("!")) {
      p->print(os) ;
      os << " " ;
    }
  }
  printMetaArgs(os) ;
  os << "]" ;
}



////////////////////////////////////////////////////////////////////////////////
/// Define default contents to print

Int_t RooAbsArg::defaultPrintContents(Option_t* /*opt*/) const
{
  return kName|kClassName|kValue|kArgs ;
}



////////////////////////////////////////////////////////////////////////////////
/// Implement multi-line detailed printing

void RooAbsArg::printMultiline(ostream& os, Int_t /*contents*/, Bool_t /*verbose*/, TString indent) const
{
  os << indent << "--- RooAbsArg ---" << endl;
  // dirty state flags
  os << indent << "  Value State: " ;
  switch(_operMode) {
  case ADirty: os << "FORCED DIRTY" ; break ;
  case AClean: os << "FORCED clean" ; break ;
  case Auto: os << (isValueDirty() ? "DIRTY":"clean") ; break ;
  }
  os << endl
     << indent << "  Shape State: " << (isShapeDirty() ? "DIRTY":"clean") << endl;
  // attribute list
  os << indent << "  Attributes: " ;
  printAttribList(os) ;
  os << endl ;
  // our memory address (for x-referencing with client addresses of other args)
  os << indent << "  Address: " << (void*)this << endl;
  // client list
  os << indent << "  Clients: " << endl;
  for (const auto client : _clientList) {
    os << indent << "    (" << (void*)client  << ","
       << (_clientListValue.containsByNamePtr(client)?"V":"-")
       << (_clientListShape.containsByNamePtr(client)?"S":"-")
       << ") " ;
    client->printStream(os,kClassName|kTitle|kName,kSingleLine);
  }

  // server list
  os << indent << "  Servers: " << endl;
  for (const auto server : _serverList) {
    os << indent << "    (" << (void*)server << ","
       << (server->_clientListValue.containsByNamePtr(this)?"V":"-")
       << (server->_clientListShape.containsByNamePtr(this)?"S":"-")
       << ") " ;
    server->printStream(os,kClassName|kName|kTitle,kSingleLine);
  }

  // proxy list
  os << indent << "  Proxies: " << endl ;
  for (int i=0 ; i<numProxies() ; i++) {
    RooAbsProxy* proxy=getProxy(i) ;
    if (!proxy) continue ;
    if (proxy->IsA()->InheritsFrom(RooArgProxy::Class())) {
      os << indent << "    " << proxy->name() << " -> " ;
      RooAbsArg* parg = ((RooArgProxy*)proxy)->absArg() ;
      if (parg) {
	parg->printStream(os,kName,kSingleLine) ;
      } else {
	os << " (empty)" << endl ; ;
      }
    } else {
      os << indent << "    " << proxy->name() << " -> " ;
      os << endl ;
      TString moreIndent(indent) ;
      moreIndent.Append("    ") ;
      ((RooSetProxy*)proxy)->printStream(os,kName,kStandard,moreIndent.Data()) ;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Print object tree structure

void RooAbsArg::printTree(ostream& os, TString /*indent*/) const
{
  const_cast<RooAbsArg*>(this)->printCompactTree(os) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Ostream operator

ostream& operator<<(ostream& os, RooAbsArg &arg)
{
  arg.writeToStream(os,kTRUE) ;
  return os ;
}

////////////////////////////////////////////////////////////////////////////////
/// Istream operator

istream& operator>>(istream& is, RooAbsArg &arg)
{
  arg.readFromStream(is,kTRUE,kFALSE) ;
  return is ;
}

////////////////////////////////////////////////////////////////////////////////
/// Print the attribute list

void RooAbsArg::printAttribList(ostream& os) const
{
  set<string>::const_iterator iter = _boolAttrib.begin() ;
  Bool_t first(kTRUE) ;
  while (iter != _boolAttrib.end()) {
    os << (first?" [":",") << *iter ;
    first=kFALSE ;
    ++iter ;
  }
  if (!first) os << "] " ;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace server nodes with names matching the dataset variable names
/// with those data set variables, making this PDF directly dependent on the dataset.

void RooAbsArg::attachDataSet(const RooAbsData &data)
{
  const RooArgSet* set = data.get() ;
  RooArgSet branches ;
  branchNodeServerList(&branches,0,kTRUE) ;

  RooFIter iter = branches.fwdIterator() ;
  RooAbsArg* branch ;
  while((branch=iter.next())) {
    branch->redirectServers(*set,kFALSE,kFALSE) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Replace server nodes with names matching the dataset variable names
/// with those data set variables, making this PDF directly dependent on the dataset

void RooAbsArg::attachDataStore(const RooAbsDataStore &dstore)
{
  const RooArgSet* set = dstore.get() ;
  RooArgSet branches ;
  branchNodeServerList(&branches,0,kTRUE) ;

  RooFIter iter = branches.fwdIterator() ;
  RooAbsArg* branch ;
  while((branch=iter.next())) {
    branch->redirectServers(*set,kFALSE,kFALSE) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Utility function used by TCollection::Sort to compare contained TObjects
/// We implement comparison by name, resulting in alphabetical sorting by object name.

Int_t RooAbsArg::Compare(const TObject* other) const
{
  return strcmp(GetName(),other->GetName()) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Print information about current value dirty state information.
/// If depth flag is true, information is recursively printed for
/// all nodes in this arg tree.

void RooAbsArg::printDirty(Bool_t depth) const
{
  if (depth) {

    RooArgSet branchList ;
    branchNodeServerList(&branchList) ;
    RooFIter bIter = branchList.fwdIterator() ;
    RooAbsArg* branch ;
    while((branch=bIter.next())) {
      branch->printDirty(kFALSE) ;
    }

  } else {
    cout << GetName() << " : " ;
    switch (_operMode) {
    case AClean: cout << "FORCED clean" ; break ;
    case ADirty: cout << "FORCED DIRTY" ; break ;
    case Auto:   cout << "Auto  " << (isValueDirty()?"DIRTY":"clean") ;
    }
    cout << endl ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Activate cache mode optimization with given definition of observables.
/// The cache operation mode of all objects in the expression tree will
/// modified such that all nodes that depend directly or indirectly on
/// any of the listed observables will be set to ADirty, as they are
/// expected to change every time. This save change tracking overhead for
/// nodes that are a priori known to change every time

void RooAbsArg::optimizeCacheMode(const RooArgSet& observables)
{
  RooLinkedList proc;
  RooArgSet opt ;
  optimizeCacheMode(observables,opt,proc) ;

  coutI(Optimization) << "RooAbsArg::optimizeCacheMode(" << GetName() << ") nodes " << opt << " depend on observables, "
			<< "changing cache operation mode from change tracking to unconditional evaluation" << endl ;
}


////////////////////////////////////////////////////////////////////////////////
/// Activate cache mode optimization with given definition of observables.
/// The cache operation mode of all objects in the expression tree will
/// modified such that all nodes that depend directly or indirectly on
/// any of the listed observables will be set to ADirty, as they are
/// expected to change every time. This save change tracking overhead for
/// nodes that are a priori known to change every time

void RooAbsArg::optimizeCacheMode(const RooArgSet& observables, RooArgSet& optimizedNodes, RooLinkedList& processedNodes)
{
  // Optimization applies only to branch nodes, not to leaf nodes
  if (!isDerived()) {
    return ;
  }


  // Terminate call if this node was already processed (tree structure may be cyclical)
  if (processedNodes.findArg(this)) {
    return ;
  } else {
    processedNodes.Add(this) ;
  }

  // Set cache mode operator to 'AlwaysDirty' if we depend on any of the given observables
  if (dependsOnValue(observables)) {

    if (dynamic_cast<RooRealIntegral*>(this)) {
      cxcoutI(Integration) << "RooAbsArg::optimizeCacheMode(" << GetName() << ") integral depends on value of one or more observables and will be evaluated for every event" << endl ;
    }
    optimizedNodes.add(*this,kTRUE) ;
    if (operMode()==AClean) {
    } else {
      setOperMode(ADirty,kTRUE) ; // WVE propagate flag recursively to top of tree
    }
  } else {
  }
  // Process any RooAbsArgs contained in any of the caches of this object
  for (Int_t i=0 ;i<numCaches() ; i++) {
    getCache(i)->optimizeCacheMode(observables,optimizedNodes,processedNodes) ;
  }

  // Forward calls to all servers
  for (const auto server : _serverList) {
    server->optimizeCacheMode(observables,optimizedNodes,processedNodes) ;
  }

}

////////////////////////////////////////////////////////////////////////////////
/// Find branch nodes with all-constant parameters, and add them to the list of
/// nodes that can be cached with a dataset in a test statistic calculation

Bool_t RooAbsArg::findConstantNodes(const RooArgSet& observables, RooArgSet& cacheList)
{
  RooLinkedList proc ;
  Bool_t ret = findConstantNodes(observables,cacheList,proc) ;

  // If node can be optimized and hasn't been identified yet, add it to the list
  coutI(Optimization) << "RooAbsArg::findConstantNodes(" << GetName() << "): components "
			<< cacheList << " depend exclusively on constant parameters and will be precalculated and cached" << endl ;

  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Find branch nodes with all-constant parameters, and add them to the list of
/// nodes that can be cached with a dataset in a test statistic calculation

Bool_t RooAbsArg::findConstantNodes(const RooArgSet& observables, RooArgSet& cacheList, RooLinkedList& processedNodes)
{
  // Caching only applies to branch nodes
  if (!isDerived()) {
    return kFALSE;
  }

  // Terminate call if this node was already processed (tree structure may be cyclical)
  if (processedNodes.findArg(this)) {
    return kFALSE ;
  } else {
    processedNodes.Add(this) ;
  }

  // Check if node depends on any non-constant parameter
  Bool_t canOpt(kTRUE) ;
  RooArgSet* paramSet = getParameters(observables) ;
  RooFIter iter = paramSet->fwdIterator() ;
  RooAbsArg* param ;
  while((param = iter.next())) {
    if (!param->isConstant()) {
      canOpt=kFALSE ;
      break ;
    }
  }
  delete paramSet ;


  if (getAttribute("NeverConstant")) {
    canOpt = kFALSE ;
  }

  if (canOpt) {
    setAttribute("ConstantExpression") ;
  }

  // If yes, list node eligible for caching, if not test nodes one level down
  if (canOpt||getAttribute("CacheAndTrack")) {

    if (!cacheList.find(*this) && dependsOnValue(observables) && !observables.find(*this) ) {

      // Add to cache list
      cxcoutD(Optimization) << "RooAbsArg::findConstantNodes(" << GetName() << ") adding self to list of constant nodes" << endl ;

      if (canOpt) setAttribute("ConstantExpressionCached") ;
      cacheList.add(*this,kFALSE) ;
    }
  }

  if (!canOpt) {

    // If not, see if next level down can be cached
    for (const auto server : _serverList) {
      if (server->isDerived()) {
        server->findConstantNodes(observables,cacheList,processedNodes) ;
      }
    }
  }

  // Forward call to all cached contained in current object
  for (Int_t i=0 ;i<numCaches() ; i++) {
    getCache(i)->findConstantNodes(observables,cacheList,processedNodes) ;
  }

  return kFALSE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Interface function signaling a request to perform constant term
/// optimization. This default implementation takes no action other than to
/// forward the calls to all servers

void RooAbsArg::constOptimizeTestStatistic(ConstOpCode opcode, Bool_t doAlsoTrackingOpt)
{
  for (const auto server : _serverList) {
    server->constOptimizeTestStatistic(opcode,doAlsoTrackingOpt) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Change cache operation mode to given mode. If recurseAdirty
/// is true, then a mode change to AlwaysDirty will automatically
/// be propagated recursively to all client nodes

void RooAbsArg::setOperMode(OperMode mode, Bool_t recurseADirty)
{
  // Prevent recursion loops
  if (mode==_operMode) return ;

  _operMode = mode ;
  _fast = ((mode==AClean) || dynamic_cast<RooRealVar*>(this)!=0 || dynamic_cast<RooConstVar*>(this)!=0 ) ;
  for (Int_t i=0 ;i<numCaches() ; i++) {
    getCache(i)->operModeHook() ;
  }
  operModeHook() ;

  // Propagate to all clients
  if (mode==ADirty && recurseADirty) {
    for (auto clientV : _clientListValue) {
      clientV->setOperMode(mode) ;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Print tree structure of expression tree on stdout, or to file if filename is specified.
/// If namePat is not "*", only nodes with names matching the pattern will be printed.
/// The client argument is used in recursive calls to properly display the value or shape nature
/// of the client-server links. It should be zero in calls initiated by users.

void RooAbsArg::printCompactTree(const char* indent, const char* filename, const char* namePat, RooAbsArg* client)
{
  if (filename) {
    ofstream ofs(filename) ;
    printCompactTree(ofs,indent,namePat,client) ;
  } else {
    printCompactTree(cout,indent,namePat,client) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Print tree structure of expression tree on given ostream.
/// If namePat is not "*", only nodes with names matching the pattern will be printed.
/// The client argument is used in recursive calls to properly display the value or shape nature
/// of the client-server links. It should be zero in calls initiated by users.

void RooAbsArg::printCompactTree(ostream& os, const char* indent, const char* namePat, RooAbsArg* client)
{
  if ( !namePat || TString(GetName()).Contains(namePat)) {
    os << indent << this ;
    if (client) {
      os << "/" ;
      if (isValueServer(*client)) os << "V" ; else os << "-" ;
      if (isShapeServer(*client)) os << "S" ; else os << "-" ;
    }
    os << " " ;

    os << IsA()->GetName() << "::" << GetName() <<  " = " ;
    printValue(os) ;

    if (!_serverList.empty()) {
      switch(operMode()) {
      case Auto:   os << " [Auto," << (isValueDirty()?"Dirty":"Clean") << "] "  ; break ;
      case AClean: os << " [ACLEAN] " ; break ;
      case ADirty: os << " [ADIRTY] " ; break ;
      }
    }
    os << endl ;

    for (Int_t i=0 ;i<numCaches() ; i++) {
      getCache(i)->printCompactTreeHook(os,indent) ;
    }
    printCompactTreeHook(os,indent) ;
  }

  TString indent2(indent) ;
  indent2 += "  " ;
  for (const auto arg : _serverList) {
    arg->printCompactTree(os,indent2,namePat,this) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Print tree structure of expression tree on given ostream, only branch nodes are printed.
/// Lead nodes (variables) will not be shown
///
/// If namePat is not "*", only nodes with names matching the pattern will be printed.

void RooAbsArg::printComponentTree(const char* indent, const char* namePat, Int_t nLevel)
{
  if (nLevel==0) return ;
  if (isFundamental()) return ;
  RooResolutionModel* rmodel = dynamic_cast<RooResolutionModel*>(this) ;
  if (rmodel && rmodel->isConvolved()) return ;
  if (InheritsFrom("RooConstVar")) return ;

  if ( !namePat || TString(GetName()).Contains(namePat)) {
    cout << indent ;
    Print() ;
  }

  TString indent2(indent) ;
  indent2 += "  " ;
  for (const auto arg : _serverList) {
    arg->printComponentTree(indent2.Data(),namePat,nLevel-1) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Construct a mangled name from the actual name that
/// is free of any math symbols that might be interpreted by TTree

TString RooAbsArg::cleanBranchName() const
{
  // Check for optional alternate name of branch for this argument
  TString rawBranchName = GetName() ;
  if (getStringAttribute("BranchName")) {
    rawBranchName = getStringAttribute("BranchName") ;
  }

  TString cleanName(rawBranchName) ;
  cleanName.ReplaceAll("/","D") ;
  cleanName.ReplaceAll("-","M") ;
  cleanName.ReplaceAll("+","P") ;
  cleanName.ReplaceAll("*","X") ;
  cleanName.ReplaceAll("[","L") ;
  cleanName.ReplaceAll("]","R") ;
  cleanName.ReplaceAll("(","L") ;
  cleanName.ReplaceAll(")","R") ;
  cleanName.ReplaceAll("{","L") ;
  cleanName.ReplaceAll("}","R") ;

  return cleanName;
}


////////////////////////////////////////////////////////////////////////////////
/// Hook function interface for object to insert additional information
/// when printed in the context of a tree structure. This default
/// implementation prints nothing

void RooAbsArg::printCompactTreeHook(ostream&, const char *)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Register RooAbsCache with this object. This function is called
/// by RooAbsCache constructors for objects that are a datamember
/// of this RooAbsArg. By registering itself the RooAbsArg is aware
/// of all its cache data members and will forward server change
/// and cache mode change calls to the cache objects, which in turn
/// can forward them their contents

void RooAbsArg::registerCache(RooAbsCache& cache)
{
  _cacheList.push_back(&cache) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Unregister a RooAbsCache. Called from the RooAbsCache destructor

void RooAbsArg::unRegisterCache(RooAbsCache& cache)
{
  _cacheList.erase(std::remove(_cacheList.begin(), _cacheList.end(), &cache),
	  _cacheList.end());
}


////////////////////////////////////////////////////////////////////////////////
/// Return number of registered caches

Int_t RooAbsArg::numCaches() const
{
  return _cacheList.size() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return registered cache object by index

RooAbsCache* RooAbsArg::getCache(Int_t index) const
{
  return _cacheList[index] ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return RooArgSet with all variables (tree leaf nodes of expresssion tree)

RooArgSet* RooAbsArg::getVariables(Bool_t stripDisconnected) const
{
  return getParameters(RooArgSet(),stripDisconnected) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return ancestors in cloning chain of this RooAbsArg. NOTE: Returned pointers
/// are not guaranteed to be 'live', so do not dereference without proper caution

RooLinkedList RooAbsArg::getCloningAncestors() const
{
  RooLinkedList retVal ;

  set<string>::const_iterator iter= _boolAttrib.begin() ;
  while(iter != _boolAttrib.end()) {
    if (TString(*iter).BeginsWith("CloneOf(")) {
      char buf[128] ;
      strlcpy(buf,iter->c_str(),128) ;
      strtok(buf,"(") ;
      char* ptrToken = strtok(0,")") ;
      RooAbsArg* ptr = (RooAbsArg*) strtol(ptrToken,0,16) ;
      retVal.Add(ptr) ;
    }
  }

  return retVal ;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a GraphViz .dot file visualizing the expression tree headed by
/// this RooAbsArg object. Use the GraphViz tool suite to make e.g. a gif
/// or ps file from the .dot file
///
/// Based on concept developed by Kyle Cranmer

void RooAbsArg::graphVizTree(const char* fileName, const char* delimiter, bool useTitle, bool useLatex)
{
  ofstream ofs(fileName) ;
  if (!ofs) {
    coutE(InputArguments) << "RooAbsArg::graphVizTree() ERROR: Cannot open graphViz output file with name " << fileName << endl ;
    return ;
  }
  graphVizTree(ofs, delimiter, useTitle, useLatex) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Write the GraphViz representation of the expression tree headed by
/// this RooAbsArg object to the given ostream.
///
/// Based on concept developed by Kyle Cranmer

void RooAbsArg::graphVizTree(ostream& os, const char* delimiter, bool useTitle, bool useLatex)
{
  if (!os) {
    coutE(InputArguments) << "RooAbsArg::graphVizTree() ERROR: output stream provided as input argument is in invalid state" << endl ;
  }

  // Write header
  os << "digraph \"" << GetName() << "\"{" << endl ;

  // First list all the tree nodes
  RooArgSet nodeSet ;
  treeNodeServerList(&nodeSet) ;
  RooFIter iter = nodeSet.fwdIterator() ;
  RooAbsArg* node ;

  // iterate over nodes
  while((node=iter.next())) {
    string nodeName = node->GetName();
    string nodeTitle = node->GetTitle();
    string nodeLabel = (useTitle && !nodeTitle.empty()) ? nodeTitle : nodeName;

    // if using latex, replace ROOT's # with normal latex backslash
    string::size_type position = nodeLabel.find("#") ;
    while(useLatex && position!=nodeLabel.npos){
      nodeLabel.replace(position, 1, "\\");
    }

    string typeFormat = "\\texttt{";
    string nodeType = (useLatex) ? typeFormat+node->IsA()->GetName()+"}" : node->IsA()->GetName();

    os << "\"" << nodeName << "\" [ color=" << (node->isFundamental()?"blue":"red")
       << ", label=\"" << nodeType << delimiter << nodeLabel << "\"];" << endl ;

  }

  // Get set of all server links
  set<pair<RooAbsArg*,RooAbsArg*> > links ;
  graphVizAddConnections(links) ;

  // And write them out
  for(auto const& link : links) {
    os << "\"" << link.first->GetName() << "\" -> \"" << link.second->GetName() << "\";" << endl ;
  }

  // Write trailer
  os << "}" << endl ;

}

////////////////////////////////////////////////////////////////////////////////
/// Utility function that inserts all point-to-point client-server connections
/// between any two RooAbsArgs in the expression tree headed by this object
/// in the linkSet argument.

void RooAbsArg::graphVizAddConnections(set<pair<RooAbsArg*,RooAbsArg*> >& linkSet)
{
  for (const auto server : _serverList) {
    linkSet.insert(make_pair(this,server)) ;
    server->graphVizAddConnections(linkSet) ;
  }
}



// //_____________________________________________________________________________
// TGraphStruct* RooAbsArg::graph(Bool_t useFactoryTag, Double_t textSize)
// {
//   // Return a TGraphStruct object filled with the tree structure of the pdf object

//   TGraphStruct* theGraph = new TGraphStruct() ;

//   // First list all the tree nodes
//   RooArgSet nodeSet ;
//   treeNodeServerList(&nodeSet) ;
//   TIterator* iter = nodeSet.createIterator() ;
//   RooAbsArg* node ;


//   // iterate over nodes
//   while((node=(RooAbsArg*)iter->Next())) {

//     // Skip node that represent numeric constants
//     if (node->IsA()->InheritsFrom(RooConstVar::Class())) continue ;

//     string nodeName ;
//     if (useFactoryTag && node->getStringAttribute("factory_tag")) {
//       nodeName  = node->getStringAttribute("factory_tag") ;
//     } else {
//       if (node->isFundamental()) {
// 	nodeName = node->GetName();
//       } else {
// 	ostringstream oss ;
// 	node->printStream(oss,(node->defaultPrintContents(0)&(~kValue)),node->defaultPrintStyle(0)) ;
// 	nodeName= oss.str() ;
// // 	nodeName = Form("%s::%s",node->IsA()->GetName(),node->GetName());

//       }
//     }
//     if (strncmp(nodeName.c_str(),"Roo",3)==0) {
//       nodeName = string(nodeName.c_str()+3) ;
//     }
//     node->setStringAttribute("graph_name",nodeName.c_str()) ;

//     TGraphNode* gnode = theGraph->AddNode(nodeName.c_str(),nodeName.c_str()) ;
//     gnode->SetLineWidth(2) ;
//     gnode->SetTextColor(node->isFundamental()?kBlue:kRed) ;
//     gnode->SetTextSize(textSize) ;
//   }
//   delete iter ;

//   // Get set of all server links
//   set<pair<RooAbsArg*,RooAbsArg*> > links ;
//   graphVizAddConnections(links) ;

//   // And insert them into the graph
//   set<pair<RooAbsArg*,RooAbsArg*> >::iterator liter = links.begin() ;
//   for( ; liter != links.end() ; ++liter ) {

//     TGraphNode* n1 = (TGraphNode*)theGraph->GetListOfNodes()->FindObject(liter->first->getStringAttribute("graph_name")) ;
//     TGraphNode* n2 = (TGraphNode*)theGraph->GetListOfNodes()->FindObject(liter->second->getStringAttribute("graph_name")) ;
//     if (n1 && n2) {
//       TGraphEdge* edge = theGraph->AddEdge(n1,n2) ;
//       edge->SetLineWidth(2) ;
//     }
//   }

//   return theGraph ;
// }



// //_____________________________________________________________________________
// Bool_t RooAbsArg::inhibitDirty()
// {
//   // Return current status of the inhibitDirty global flag. If true
//   // no dirty state change tracking occurs and all caches are considered
//   // to be always dirty.
//   return _inhibitDirty ;
// }


////////////////////////////////////////////////////////////////////////////////
/// Take ownership of the contents of 'comps'

Bool_t RooAbsArg::addOwnedComponents(const RooArgSet& comps)
{
  if (!_ownedComponents) {
    _ownedComponents = new RooArgSet("owned components") ;
  }
  return _ownedComponents->addOwned(comps) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Clone tree expression of objects. All tree nodes will be owned by
/// the head node return by cloneTree()

RooAbsArg* RooAbsArg::cloneTree(const char* newname) const
{
  // Clone tree using snapshot
  RooArgSet* clonedNodes = (RooArgSet*) RooArgSet(*this).snapshot(kTRUE) ;

  // Find the head node in the cloneSet
  RooAbsArg* head = clonedNodes->find(*this) ;
  assert(head);

  // Remove the head node from the cloneSet
  // To release it from the set ownership
  clonedNodes->remove(*head) ;

  // Add the set as owned component of the head
  head->addOwnedComponents(*clonedNodes) ;

  // Delete intermediate container
  clonedNodes->releaseOwnership() ;
  delete clonedNodes ;

  // Adjust name of head node if requested
  if (newname) {
    head->TNamed::SetName(newname) ;
    head->_namePtr = (TNamed*) RooNameReg::instance().constPtr(newname) ;
  }

  // Return the head
  return head ;
}



////////////////////////////////////////////////////////////////////////////////

void RooAbsArg::attachToStore(RooAbsDataStore& store)
{
  if (dynamic_cast<RooTreeDataStore*>(&store)) {
    attachToTree(((RooTreeDataStore&)store).tree()) ;
  } else if (dynamic_cast<RooVectorDataStore*>(&store)) {
    attachToVStore((RooVectorDataStore&)store) ;
  }
}



////////////////////////////////////////////////////////////////////////////////

RooExpensiveObjectCache& RooAbsArg::expensiveObjectCache() const
{
  if (_eocache) {
    return *_eocache ;
  } else {
    return RooExpensiveObjectCache::instance() ;
  }
}


////////////////////////////////////////////////////////////////////////////////

const char* RooAbsArg::aggregateCacheUniqueSuffix() const
{
  string suffix ;

  RooArgSet branches ;
  branchNodeServerList(&branches) ;
  RooFIter iter = branches.fwdIterator();
  RooAbsArg* arg ;
  while((arg=iter.next())) {
    const char* tmp = arg->cacheUniqueSuffix() ;
    if (tmp) suffix += tmp ;
  }
  return Form("%s",suffix.c_str()) ;
}


////////////////////////////////////////////////////////////////////////////////

void RooAbsArg::wireAllCaches()
{
  RooArgSet branches ;
  branchNodeServerList(&branches) ;
  RooFIter iter = branches.fwdIterator() ;
  RooAbsArg* arg ;
  while((arg=iter.next())) {
//     cout << "wiring caches on node " << arg->GetName() << endl ;
    for (deque<RooAbsCache*>::iterator iter2 = arg->_cacheList.begin() ; iter2 != arg->_cacheList.end() ; ++iter2) {
      (*iter2)->wireCache() ;
    }
  }
}



////////////////////////////////////////////////////////////////////////////////

void RooAbsArg::SetName(const char* name)
{
  TNamed::SetName(name) ;
  TNamed* newPtr = (TNamed*) RooNameReg::instance().constPtr(GetName()) ;
  if (newPtr != _namePtr) {
    //cout << "Rename '" << _namePtr->GetName() << "' to '" << name << "' (set flag in new name)" << endl;
    _namePtr = newPtr;
    _namePtr->SetBit(RooNameReg::kRenamedArg);
  }
}




////////////////////////////////////////////////////////////////////////////////

void RooAbsArg::SetNameTitle(const char *name, const char *title)
{
  TNamed::SetNameTitle(name,title) ;
  TNamed* newPtr = (TNamed*) RooNameReg::instance().constPtr(GetName()) ;
  if (newPtr != _namePtr) {
    //cout << "Rename '" << _namePtr->GetName() << "' to '" << name << "' (set flag in new name)" << endl;
    _namePtr = newPtr;
    _namePtr->SetBit(RooNameReg::kRenamedArg);
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooAbsArg.

void RooAbsArg::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
     _ioReadStack.push(this) ;
     R__b.ReadClassBuffer(RooAbsArg::Class(),this);
     _ioReadStack.pop() ;
     _namePtr = (TNamed*) RooNameReg::instance().constPtr(GetName()) ;
     _isConstant = getAttribute("Constant") ;
   } else {
     R__b.WriteClassBuffer(RooAbsArg::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Method called by workspace container to finalize schema evolution issues
/// that cannot be handled in a single ioStreamer pass.
///
/// A second pass is typically needed when evolving data member of RooAbsArg-derived
/// classes that are container classes with references to other members, which may
/// not yet be 'live' in the first ioStreamer() evolution pass.
///
/// Classes may overload this function, but must call the base method in the
/// overloaded call to ensure base evolution is handled properly

void RooAbsArg::ioStreamerPass2()
{
  // Handling of v5-v6 migration (TRefArray _proxyList --> RooRefArray _proxyList)
  auto iter = _ioEvoList.find(this);
  if (iter != _ioEvoList.end()) {

    // Transfer contents of saved TRefArray to RooRefArray now
    if (!_proxyList.GetEntriesFast())
       _proxyList.Expand(iter->second->GetEntriesFast());
    for (int i = 0; i < iter->second->GetEntriesFast(); i++) {
       _proxyList.Add(iter->second->At(i));
    }
    // Delete TRefArray and remove from list
    _ioEvoList.erase(iter);
  }
}




////////////////////////////////////////////////////////////////////////////////
/// Method called by workspace container to finalize schema evolution issues
/// that cannot be handled in a single ioStreamer pass. This static finalize method
/// is called after ioStreamerPass2() is called on each directly listed object
/// in the workspace. It's purpose is to complete schema evolution of any
/// objects in the workspace that are not directly listed as content elements
/// (e.g. analytical convolution tokens )

void RooAbsArg::ioStreamerPass2Finalize()
{
  // Handling of v5-v6 migration (TRefArray _proxyList --> RooRefArray _proxyList)
  for (const auto& iter : _ioEvoList) {

    // Transfer contents of saved TRefArray to RooRefArray now
    if (!iter.first->_proxyList.GetEntriesFast())
       iter.first->_proxyList.Expand(iter.second->GetEntriesFast());
    for (int i = 0; i < iter.second->GetEntriesFast(); i++) {
       iter.first->_proxyList.Add(iter.second->At(i));
    }
  }

  _ioEvoList.clear();
}


RooAbsArg::RefCountListLegacyIterator_t *
RooAbsArg::makeLegacyIterator(const RooAbsArg::RefCountList_t& list) const {
  return new RefCountListLegacyIterator_t(list.containedObjects());
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooRefArray.

void RooRefArray::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {

      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }

      // Make temporary refArray and read that from the streamer
      auto refArray = std::make_unique<TRefArray>();
      refArray->Streamer(R__b) ;
      R__b.CheckByteCount(R__s, R__c, refArray->IsA());

      // Schedule deferred processing of TRefArray into proxy list
      RooAbsArg::_ioEvoList[RooAbsArg::_ioReadStack.top()] = std::move(refArray);

   } else {

     R__c = R__b.WriteVersion(RooRefArray::IsA(), kTRUE);

     // Make a temporary refArray and write that to the streamer
     TRefArray refArray(GetEntriesFast());
     TIterator* iter = MakeIterator() ;
     TObject* tmpObj ; while ((tmpObj = iter->Next())) {
       refArray.Add(tmpObj) ;
     }
     delete iter ;

     refArray.Streamer(R__b) ;
     R__b.SetByteCount(R__c, kTRUE) ;

   }
}

/// Print at the prompt
namespace cling {
std::string printValue(RooAbsArg *raa)
{
   std::stringstream s;
   if (0 == *raa->GetName() && 0 == *raa->GetTitle()) {
      s << "An instance of " << raa->ClassName() << ".";
      return s.str();
   }
   raa->printStream(s, raa->defaultPrintContents(""), raa->defaultPrintStyle(""));
   return s.str();
}
} // namespace cling

