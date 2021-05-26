/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooCustomizer.h,v 1.11 2007/05/11 09:11:30 verkerke Exp $
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

#ifndef ROO_PDF_CUSTOMIZER
#define ROO_PDF_CUSTOMIZER

#include "TList.h"
#include "TNamed.h"
#include "TString.h"
#include "RooArgSet.h"
#include "RooPrintable.h"
#include "RooFactoryWSTool.h"

#include <vector>
#include <string>

class RooAbsCategoryLValue ;
class RooAbsCategory ;
class RooAbsArg ;
class RooAbsPdf ;

class RooCustomizer : public TNamed, public RooPrintable {

public:

  // Constructors, assignment etc
  RooCustomizer(const RooAbsArg& pdf, const RooAbsCategoryLValue& masterCat, RooArgSet& splitLeafListOwned, RooArgSet* splitLeafListAll=0) ;
  RooCustomizer(const RooAbsArg& pdf, const char* name) ;
  virtual ~RooCustomizer() ;

  void setOwning(Bool_t flag) {
    // If flag is true, make customizer own all created components
    _owning = flag ;
  }

  void splitArgs(const RooArgSet& argSet, const RooAbsCategory& splitCat) ;
  void splitArg(const RooAbsArg& arg, const RooAbsCategory& splitCat) ;
  void replaceArg(const RooAbsArg& orig, const RooAbsArg& subst) ;
  RooAbsArg* build(const char* masterCatState, Bool_t verbose=kFALSE) ;
  RooAbsArg* build(Bool_t verbose=kFALSE) ;

  const RooArgSet& cloneBranchList() const {
    // Return list of cloned branch nodes
    return *_cloneBranchList ;
  }
  const RooArgSet& cloneLeafList() const {
    // Return list of cloned leaf nodes
    return *_cloneNodeListOwned ;
  }

  // Printing interface
  virtual void printName(std::ostream& os) const ;
  virtual void printTitle(std::ostream& os) const ;
  virtual void printClassName(std::ostream& os) const ;
  virtual void printArgs(std::ostream& os) const ;
  virtual void printMultiline(std::ostream& os, Int_t content, Bool_t verbose=kFALSE, TString indent= "") const;

  inline virtual void Print(Option_t *options= 0) const {
    // Printing interface
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  // Releases ownership of list of cloned branch nodes
  void setCloneBranchSet(RooArgSet& cloneBranchSet) ;

  // Factory interface
  class CustIFace : public RooFactoryWSTool::IFace {
  public:
    virtual ~CustIFace() {} ;
    std::string create(RooFactoryWSTool& ft, const char* typeName, const char* instanceName, std::vector<std::string> args) ;
  } ;

protected:

  RooCustomizer(const RooCustomizer&) ;
  void initialize() ;

  RooAbsArg* doBuild(const char* masterCatState, Bool_t verbose) ;

  Bool_t _sterile ; // If true we do not have as associated master category
  Bool_t _owning ;  // If true we own all created components
  TString _name ;   // Name of this object

  TList _splitArgList ; // List of RooAbsArgs to be split
  TList _splitCatList ; // List of categories to be used for above splits

  TList _replaceArgList ; // List of RooAbsArgs to be replaced
  TList _replaceSubList ; // List of replacement RooAbsArgs

  // Master nodes are not owned
  RooAbsArg* _masterPdf ;             // Pointer to input p.d.f
  RooAbsCategoryLValue* _masterCat ;  // Pointer to input master category

  RooArgSet  _masterBranchList ;      // List of branch nodes
  RooArgSet  _masterLeafList ;        // List of leaf nodes

  RooArgSet  _internalCloneBranchList ; // List of branches of internal clone
  RooArgSet* _cloneBranchList ;         // Pointer to list of cloned branches used

  // Cloned leafs are owned by the user supplied list in the ctor
  RooArgSet* _cloneNodeListAll ;        // List of all cloned nodes
  RooArgSet* _cloneNodeListOwned ;      // List of owned cloned nodes

  ClassDef(RooCustomizer,0) // Editing tool for RooAbsArg composite object expressions
} ;

#endif
