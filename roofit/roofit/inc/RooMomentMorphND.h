/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef ROOMOMENTMORPHND
#define ROOMOMENTMORPHND

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooSetProxy.h"
#include "RooListProxy.h"
#include "RooArgList.h"
#include "RooBinning.h"

#include "TMatrixD.h"
#include "TMap.h"

#include <vector>
#include <map>

class RooChangeTracker;

class RooMomentMorphND : public RooAbsPdf {

public:
   class Grid {
   public:
      Grid(){};
      Grid(const Grid &other);
      Grid(const RooAbsBinning &binning_x) { _grid.push_back(binning_x.clone()); };
      Grid(const RooAbsBinning &binning_x, const RooAbsBinning &binning_y)
      {
         _grid.push_back(binning_x.clone());
         _grid.push_back(binning_y.clone());
      };
      Grid(const RooAbsBinning &binning_x, const RooAbsBinning &binning_y, const RooAbsBinning &binning_z)
      {
         _grid.push_back(binning_x.clone());
         _grid.push_back(binning_y.clone());
         _grid.push_back(binning_z.clone());
      };
      Grid(const std::vector<RooAbsBinning *> binnings)
      {
         for (unsigned int i = 0; i < binnings.size(); i++) {
            _grid.push_back(binnings[i]->clone());
         }
      };

      ~Grid();

      void addPdf(const RooAbsPdf &pdf, int bin_x);
      void addPdf(const RooAbsPdf &pdf, int bin_x, int bin_y);
      void addPdf(const RooAbsPdf &pdf, int bin_x, int bin_y, int bin_z);
      void addPdf(const RooAbsPdf &pdf, std::vector<int> bins);
      void addBinning(const RooAbsBinning &binning) { _grid.push_back(binning.clone()); };

      mutable std::vector<RooAbsBinning *> _grid;
      mutable RooArgList _pdfList;
      mutable std::map<std::vector<int>, int> _pdfMap;

      mutable std::vector<std::vector<double>> _nref;
      mutable std::vector<int> _nnuis;
   };

protected:
   class CacheElem : public RooAbsCacheElement {
   public:
      CacheElem(RooAbsPdf &sumPdf, RooChangeTracker &tracker, const RooArgList &flist)
         : _sumPdf(&sumPdf), _tracker(&tracker)
      {
         _frac.add(flist);
      };
      virtual ~CacheElem();
      virtual RooArgList containedArgs(Action);
      RooAbsPdf *_sumPdf;
      RooChangeTracker *_tracker;
      RooArgList _frac;

      RooRealVar *frac(int i);
      const RooRealVar *frac(int i) const;
      void calculateFractions(const RooMomentMorphND &self, Bool_t verbose = kTRUE) const;
   };

public:
   enum Setting { Linear, SineLinear, NonLinear, NonLinearPosFractions, NonLinearLinFractions };

   RooMomentMorphND();
   RooMomentMorphND(const char *name, const char *title, RooAbsReal &_m, const RooArgList &varList,
                    const RooArgList &pdfList, const RooArgList &mrefList, Setting setting);
   RooMomentMorphND(const char *name, const char *title, const RooArgList &parList, const RooArgList &obsList,
                    const Grid &referenceGrid, const Setting &setting);
   RooMomentMorphND(const RooMomentMorphND &other, const char *name = 0);
   RooMomentMorphND(const char *name, const char *title, RooAbsReal &_m, const RooArgList &varList,
                    const RooArgList &pdfList, const TVectorD &mrefpoints, Setting setting);
   virtual ~RooMomentMorphND();
   virtual TObject *clone(const char *newname) const { return new RooMomentMorphND(*this, newname); }

   void setMode(const Setting &setting) { _setting = setting; }
   virtual Bool_t selfNormalized() const { return kTRUE; }
   Bool_t setBinIntegrator(RooArgSet &allVars);
   void useHorizontalMorphing(Bool_t val) { _useHorizMorph = val; }

   Double_t evaluate() const;
   virtual Double_t getVal(const RooArgSet *set = 0) const;

protected:
   void initialize();
   void initializeParameters(const RooArgList &parList);
   void initializeObservables(const RooArgList &obsList);

   RooAbsPdf *sumPdf(const RooArgSet *nset);
   CacheElem *getCache(const RooArgSet *nset) const;

   void findShape(const std::vector<double> &x) const;

   friend class CacheElem;
   friend class Grid;

   mutable RooObjCacheManager _cacheMgr;
   mutable RooArgSet *_curNormSet;

   RooListProxy _parList;
   RooSetProxy _obsList;
   // RooListProxy _pdfList ;
   TIterator *_parItr; //! Do not persist
   TIterator *_obsItr; //! Do not persist
   mutable Grid _referenceGrid;
   RooListProxy _pdfList;

   mutable TMatrixD *_M;
   mutable TMatrixD *_MSqr;
   mutable std::vector<std::vector<double>> _squareVec;
   mutable std::vector<int> _squareIdx;

   Setting _setting;
   Bool_t _useHorizMorph;

   inline int sij(const int &i, const int &j) const { return (i * _obsList.getSize() + j); }

   ClassDef(RooMomentMorphND, 1)
};

#endif
