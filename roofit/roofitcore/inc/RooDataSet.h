/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooDataSet.h,v 1.59 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_DATA_SET
#define ROO_DATA_SET

class TDirectory;
class RooAbsRealLValue;
class RooRealVar;
class RooDataHist;

#include "RooAbsData.h"
#include "RooDirItem.h"

#include "ROOT/RStringView.hxx"

#include <list>


#define USEMEMPOOLFORDATASET
template <class RooSet_t, size_t>
class MemPoolForRooSets;

class RooDataSet : public RooAbsData, public RooDirItem {
public:

#ifdef USEMEMPOOLFORDATASET
  void* operator new (size_t bytes);
  void operator delete (void *ptr);
#endif
 

  // Constructors, factory methods etc.
  RooDataSet() ; 

  // Empty constructor 
  RooDataSet(std::string_view name, std::string_view title, const RooArgSet& vars, const char* wgtVarName=0) ;

  // Universal constructor
  RooDataSet(std::string_view name, std::string_view title, const RooArgSet& vars, const RooCmdArg& arg1, const RooCmdArg& arg2=RooCmdArg(), 
	     const RooCmdArg& arg3=RooCmdArg(), const RooCmdArg& arg4=RooCmdArg(),const RooCmdArg& arg5=RooCmdArg(),
	     const RooCmdArg& arg6=RooCmdArg(),const RooCmdArg& arg7=RooCmdArg(),const RooCmdArg& arg8=RooCmdArg()) ; 

    // Constructor for subset of existing dataset
  RooDataSet(std::string_view name, std::string_view title, RooDataSet *data, const RooArgSet& vars, 
             const char *cuts=0, const char* wgtVarName=0);
  RooDataSet(std::string_view name, std::string_view title, RooDataSet *data, const RooArgSet& vars,  
	     const RooFormulaVar& cutVar, const char* wgtVarName=0) ;  


  // Constructor importing data from external ROOT Tree
  RooDataSet(std::string_view name, std::string_view title, TTree *tree, const RooArgSet& vars,
	     const char *cuts=0, const char* wgtVarName=0); 
  RooDataSet(std::string_view name, std::string_view title, TTree *tree, const RooArgSet& vars,
	     const RooFormulaVar& cutVar, const char* wgtVarName=0) ;  

  WRITE_TSTRING_COMPATIBLE_CONSTRUCTOR(RooDataSet)

  RooDataSet(RooDataSet const & other, const char* newname=0) ;  
  virtual TObject* Clone(const char* newname = "") const override {
    return new RooDataSet(*this, newname && newname[0] != '\0' ? newname : GetName());
  }
  virtual ~RooDataSet() ;

  virtual RooAbsData* emptyClone(const char* newName=0, const char* newTitle=0, const RooArgSet* vars=0, const char* wgtVarName=0) const override;

  RooDataHist* binnedClone(const char* newName=0, const char* newTitle=0) const ;

  virtual Double_t sumEntries() const override;
  virtual Double_t sumEntries(const char* cutSpec, const char* cutRange=0) const override;

  virtual RooPlot* plotOnXY(RooPlot* frame, 
			    const RooCmdArg& arg1=RooCmdArg::none(), const RooCmdArg& arg2=RooCmdArg::none(),
			    const RooCmdArg& arg3=RooCmdArg::none(), const RooCmdArg& arg4=RooCmdArg::none(),
			    const RooCmdArg& arg5=RooCmdArg::none(), const RooCmdArg& arg6=RooCmdArg::none(),
			    const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none()) const ;


  // Read data from a text file and create a dataset from it.
  // The possible options are: (D)ebug, (Q)uiet.
  static RooDataSet *read(const char *filename, const RooArgList &variables,
			  const char *opts= "", const char* commonPath="",
			  const char *indexCatName=0) ;
  Bool_t write(const char* filename) const;
  Bool_t write(std::ostream & ofs) const;

/*   void setWeightVar(const char* name=0) ; */
/*   void setWeightVar(const RooAbsArg& arg) {  */
/*     // Interpret given argument as event weight */
/*     setWeightVar(arg.GetName()) ;  */
/*   } */
  virtual Bool_t isWeighted() const override;
  virtual Bool_t isNonPoissonWeighted() const override;

  virtual Double_t weight() const override;
  virtual Double_t weightSquared() const override;
  virtual void weightError(double& lo, double& hi,ErrorType etype=SumW2) const override;
  double weightError(ErrorType etype=SumW2) const override;

  virtual const RooArgSet* get(Int_t index) const override;
  virtual const RooArgSet* get() const override;

  void getBatches(RooBatchCompute::RunContext& evalData,
      std::size_t first = 0, std::size_t len = std::numeric_limits<std::size_t>::max()) const override;
  virtual RooSpan<const double> getWeightBatch(std::size_t first, std::size_t len) const override;

  // Add one ore more rows of data
  virtual void add(const RooArgSet& row, Double_t weight=1.0, Double_t weightError=0) override;
  virtual void add(const RooArgSet& row, Double_t weight, Double_t weightErrorLo, Double_t weightErrorHi);

  virtual void addFast(const RooArgSet& row, Double_t weight=1.0, Double_t weightError=0);

  void append(RooDataSet& data) ;
  Bool_t merge(RooDataSet* data1, RooDataSet* data2=0, RooDataSet* data3=0,  
 	       RooDataSet* data4=0, RooDataSet* data5=0, RooDataSet* data6=0) ; 
  Bool_t merge(std::list<RooDataSet*> dsetList) ;

  virtual RooAbsArg* addColumn(RooAbsArg& var, Bool_t adjustRange=kTRUE) ;
  virtual RooArgSet* addColumns(const RooArgList& varList) ;

  // Plot the distribution of a real valued arg
  using RooAbsData::createHistogram ;
  TH2F* createHistogram(const RooAbsRealLValue& var1, const RooAbsRealLValue& var2, const char* cuts="", 
			const char *name= "hist") const;	 
  TH2F* createHistogram(const RooAbsRealLValue& var1, const RooAbsRealLValue& var2, Int_t nx, Int_t ny,
                        const char* cuts="", const char *name="hist") const;

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent="") const override;
  virtual void printArgs(std::ostream& os) const override;
  virtual void printValue(std::ostream& os) const override;

  void SetName(const char *name) override;
  void SetNameTitle(const char *name, const char* title) override;

  static void cleanup();

  void convertToTreeStore() override;

protected:

  virtual RooAbsData* cacheClone(const RooAbsArg* newCacheOwner, const RooArgSet* newCacheVars, const char* newName=0) override;

  friend class RooProdGenContext ;

  void initialize(const char* wgtVarName) ;
  
  // Cache copy feature is not publicly accessible
  RooAbsData* reduceEng(const RooArgSet& varSubset, const RooFormulaVar* cutVar, const char* cutRange=0, 
	                std::size_t nStart=0, std::size_t nStop = std::numeric_limits<std::size_t>::max(), Bool_t copyCache=kTRUE) override;
  RooDataSet(std::string_view name, std::string_view title, RooDataSet *ntuple, 
	     const RooArgSet& vars, const RooFormulaVar* cutVar, const char* cutRange, std::size_t nStart, std::size_t nStop, Bool_t copyCache, const char* wgtVarName=0);
  
  RooArgSet addWgtVar(const RooArgSet& origVars, const RooAbsArg* wgtVar) ; 
  
  RooArgSet _varsNoWgt ;   // Vars without weight variable 
  RooRealVar* _wgtVar ;    // Pointer to weight variable (if set) 

private:
#ifdef USEMEMPOOLFORDATASET
  typedef MemPoolForRooSets<RooDataSet, 5*150> MemPool; // 150 = about 100kb
  static MemPool * memPool();
#endif
  unsigned short _errorMsgCount{0}; //! Counter to silence error messages when filling dataset.
  bool _doWeightErrorCheck{true}; //! When adding events with weights, check that weights can actually be stored.

  ClassDefOverride(RooDataSet,2) // Unbinned data set
};

#endif
