#ifndef RooNLLVarNew_h
#define RooNLLVarNew_h

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooTemplateProxy.h"

#include "RooBatchCompute.h"

class RooNLLVarNew : public RooAbsReal {

public:
  RooNLLVarNew(){};
  RooNLLVarNew(const char *name, const char *title, RooAbsPdf &pdf,
               RooArgSet const& observables,
               RooAbsReal* weight, RooAbsReal* constraints,
               bool isExtended);
  RooNLLVarNew(const RooNLLVarNew &other, const char *name = 0);
  virtual TObject *clone(const char *newname) const override { return new RooNLLVarNew(*this, newname); }

  using RooAbsArg::getParameters;
  bool getParameters(const RooArgSet *depList, RooArgSet& outSet, bool stripDisconnected=true) const override;

  virtual Double_t defaultErrorLevel() const override
  {
    // Return default level for MINUIT error analysis
    return 0.5;
  }

  inline RooAbsPdf* getPdf() const { return &*_pdf; }
  void computeBatch(double* output, size_t nEvents, rbc::DataMap& dataMap) const override;
  double reduce(const double* input, size_t nEvents) const;

protected:
  RooTemplateProxy<RooAbsPdf> _pdf;
  RooArgSet const* _observables = nullptr;
  std::unique_ptr<RooTemplateProxy<RooAbsReal>> _weight = nullptr;
  RooAbsReal const* _constraints = nullptr;
  mutable double _sumWeight = 0.0; //!
  bool _isExtended;

  double getValV(const RooArgSet *normalisationSet = nullptr) const override;

  double evaluate() const override;

  RooSpan<double> evaluateSpan(RooBatchCompute::RunContext &evalData, const RooArgSet *normSet) const override;

  RooSpan<const double>
  getValues(RooBatchCompute::RunContext &evalData, const RooArgSet *normSet = nullptr) const override;

private:
  ClassDefOverride(RooNLLVarNew, 1)
};

#endif
