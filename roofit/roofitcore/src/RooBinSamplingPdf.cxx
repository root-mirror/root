// Authors: Stephan Hageboeck, CERN; Andrea Sciandra, SCIPP-UCSC/Atlas; Nov 2020

/*****************************************************************************
 * RooFit
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2020, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
 * \class RooBinSamplingPdf
 * The RooBinSamplingPdf is supposed to be used as an adapter between a continuous PDF
 * and a binned distribution.
 * When RooFit is used to fit binned data, and the PDF is continuous, it takes the probability density
 * at the bin centre as a proxy for the probability averaged (integrated) over the entire bin. This is
 * correct only if the second derivative of the function vanishes, though.
 *
 * For PDFs that have larger curvatures, the RooBinSamplingPdf can be used. It integrates the PDF in each
 * bin using an adaptive integrator. This usually requires 21 times more function evaluations, but significantly
 * reduces biases due to better sampling of the PDF. The integrator can be accessed from the outside
 * using integrator(). This can be used to change the integration rules, so less/more function evaluations are
 * performed. The target precision of the integrator can be set in the constructor.
 *
 * \note This feature is currently limited to one-dimensional PDFs.
 *
 * ### How to use it
 * There are two ways to use this class:
 * - Manually wrap a PDF:
 * ```
 *   RooBinSamplingPdf binSampler("<name>", "title", <binned observable of PDF>, <original PDF> [, <precision for integrator>]);
 *   binSampler.fitTo(data);
 * ```
 *   When a PDF is wrapped with a RooBinSamplingPDF, just use the bin sampling PDF instead of the original one for fits
 *   or plotting etc. Note that the binning will be taken from the observable.
 * - Instruct test statistics to carry out this wrapping automatically:
 * ```
 *   pdf.fitTo(data, IntegrateBins(<precision>));
 * ```
 *   This method is especially useful when used with a simultaneous PDF, since each component will automatically be wrapped,
 *   depending on the value of `precision`:
 *   - `precision < 0.`: None of the PDFs are touched, bin sampling is off.
 *   - `precision = 0.`: Continuous PDFs that are fit to a RooDataHist are wrapped into a RooBinSamplingPdf. The target precision
 *      forwarded to the integrator is 1.E-4 (the default argument of the constructor).
 *   - `precision > 0.`: All continuous PDFs are automatically wrapped into a RooBinSamplingPdf. The `'precision'` is used for all
 *      integrators.
 *
 * ### Simulating a binned fit using RooDataSet
 *   Some frameworks use unbinned data (RooDataSet) to simulate binned datasets. By adding one entry for each bin centre with the
 *   appropriate weight, one can achieve the same result as fitting with RooDataHist. In this case, however, RooFit cannot
 *   auto-detect that a binned fit is running, and that an integration over the bin is desired (note that there are no bins to
 *   integrate over in this kind of dataset).
 *
 *   In this case, `IntegrateBins(>0.)` needs to be used, and the desired binning needs to be assigned to the observable
 *   of the dataset:
 *   ```
 *     RooRealVar x("x", "x", 0., 5.);
 *     x.setBins(10);
 *
 *     // <create dataset and model>
 *
 *     model.fitTo(data, IntegrateBins(>0.));
 *   ```
 *
 *   \see RooAbsPdf::fitTo()
 *   \see IntegrateBins()
 *
 *
 * \htmlonly <style>div.image img[src="RooBinSamplingPdf_OFF.png"]{width:15cm;}</style> \endhtmlonly
 * \htmlonly <style>div.image img[src="RooBinSamplingPdf_ON.png" ]{width:15cm;}</style> \endhtmlonly
 * <table>
 * <tr><td>
 * \image html RooBinSamplingPdf_OFF.png "Binned fit without RooBinSamplingPdf"
 * <td>
 * \image html RooBinSamplingPdf_ON.png "Binned fit with RooBinSamplingPdf"
 * </table>
 *
 */


#include "RooBinSamplingPdf.h"

#include "RooHelpers.h"
#include "RooRealBinding.h"
#include "RunContext.h"

#include "Math/Integrator.h"

#include <algorithm>

////////////////////////////////////////////////////////////////////////////////
/// Construct a new RooBinSamplingPdf.
/// \param[in] name A name to identify this object.
/// \param[in] title Title (for e.g. plotting)
/// \param[in] observable Observable to integrate over (the one that is binned).
/// \param[in] inputPdf A PDF whose bins should be sampled with higher precision.
/// \param[in] epsilon Relative precision for the integrator, which is used to sample the bins.
/// Note that ROOT's default is to use an adaptive integrator, which in its first iteration usually reaches
/// relative precision of 1.E-4 or better. Therefore, asking for lower precision rarely has an effect.
RooBinSamplingPdf::RooBinSamplingPdf(const char *name, const char *title, RooAbsRealLValue& observable,
    RooAbsPdf& inputPdf, double epsilon) :
      RooAbsPdf(name, title),
      _pdf("inputPdf", "Function to be converted into a PDF", this, inputPdf),
      _observable("observable", "Observable to integrate over", this, observable, true, true),
      _relEpsilon(epsilon) {
  if (!_pdf->dependsOn(*_observable)) {
    throw std::invalid_argument(std::string("RooBinSamplingPDF(") + GetName()
        + "): The PDF " + _pdf->GetName() + " needs to depend on the observable "
        + _observable->GetName());
  }
}


 ////////////////////////////////////////////////////////////////////////////////
 /// Copy a RooBinSamplingPdf.
 /// \param[in] other PDF to copy.
 /// \param[in] name Optionally rename the copy.
 RooBinSamplingPdf::RooBinSamplingPdf(const RooBinSamplingPdf& other, const char* name) :
   RooAbsPdf(other, name),
   _pdf("inputPdf", this, other._pdf),
   _observable("observable", this, other._observable),
   _relEpsilon(other._relEpsilon) { }


////////////////////////////////////////////////////////////////////////////////
/// Integrate the PDF over the current bin of the observable.
double RooBinSamplingPdf::evaluate() const {
  const unsigned int bin = _observable->getBin();
  const double low = _observable->getBinning().binLow(bin);
  const double high = _observable->getBinning().binHigh(bin);

  const double oldX = _observable->getVal();
  double result;
  {
    // Important: When the integrator samples x, caching of sub-tree values needs to be off.
    RooHelpers::DisableCachingRAII disableCaching(inhibitDirty());
    result = integrate(_normSet, low, high) / (high-low);
  }

  _observable->setVal(oldX);

  return result;
}


////////////////////////////////////////////////////////////////////////////////
/// Integrate the PDF over all its bins, and return a batch with those values.
/// \param[in/out] evalData Struct with evaluation data.
/// \param[in] normSet Normalisation set that's used to evaluate the PDF.
RooSpan<double> RooBinSamplingPdf::evaluateSpan(RooBatchCompute::RunContext& evalData, const RooArgSet* normSet) const {
  // Retrieve binning, which we need to compute the probabilities
  auto boundaries = binBoundaries();
  auto xValues = _observable->getValues(evalData, normSet);
  auto results = evalData.makeBatch(this, xValues.size());

  // Important: When the integrator samples x, caching of sub-tree values needs to be off.
  RooHelpers::DisableCachingRAII disableCaching(inhibitDirty());

  // Now integrate PDF in each bin:
  for (unsigned int i=0; i < xValues.size(); ++i) {
    const double x = xValues[i];
    const auto upperIt = std::upper_bound(boundaries.begin(), boundaries.end(), x);
    const unsigned int bin = std::distance(boundaries.begin(), upperIt) - 1;
    assert(bin < boundaries.size());

    results[i] = integrate(normSet, boundaries[bin], boundaries[bin+1]) / (boundaries[bin+1]-boundaries[bin]);
  }

  return results;
}


////////////////////////////////////////////////////////////////////////////////
/// Get the bin boundaries for the observable.
/// These will be recomputed whenever the shape of this object is dirty.
RooSpan<const double> RooBinSamplingPdf::binBoundaries() const {
  if (isShapeDirty() || _binBoundaries.empty()) {
    _binBoundaries.clear();
    const RooAbsBinning& binning = _observable->getBinning(nullptr);
    const double* boundaries = binning.array();

    for (int i=0; i < binning.numBoundaries(); ++i) {
      _binBoundaries.push_back(boundaries[i]);
    }

    assert(std::is_sorted(_binBoundaries.begin(), _binBoundaries.end()));

    clearShapeDirty();
  }

  return {_binBoundaries};
}


////////////////////////////////////////////////////////////////////////////////
/// Return a list of all bin boundaries, so the PDF is plotted correctly.
/// \param[in] obs Observable to generate the boundaries for.
/// \param[in] xlo Beginning of range to create list of boundaries for.
/// \param[in] xhi End of range to create to create list of boundaries for.
/// \return Pointer to a list to be deleted by caller.
std::list<double>* RooBinSamplingPdf::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const {
  if (obs.namePtr() != _observable->namePtr()) {
    coutE(Plotting) << "RooBinSamplingPdf::binBoundaries(" << GetName() << "): observable '" << obs.GetName()
        << "' is not the observable of this PDF ('" << _observable->GetName() << "')." << std::endl;
    return nullptr;
  }

  auto list = new std::list<double>;
  for (double val : binBoundaries()) {
    if (xlo <= val && val < xhi)
      list->push_back(val);
  }

  return list;
}


////////////////////////////////////////////////////////////////////////////////
/// Return a list of all bin edges, so the PDF is plotted as a step function.
/// \param[in] obs Observable to generate the sampling hint for.
/// \param[in] xlo Beginning of range to create sampling hint for.
/// \param[in] xhi End of range to create sampling hint for.
/// \return Pointer to a list to be deleted by caller.
std::list<double>* RooBinSamplingPdf::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const {
  if (obs.namePtr() != _observable->namePtr()) {
    coutE(Plotting) << "RooBinSamplingPdf::plotSamplingHint(" << GetName() << "): observable '" << obs.GetName()
        << "' is not the observable of this PDF ('" << _observable->GetName() << "')." << std::endl;
    return nullptr;
  }

  auto binEdges = new std::list<double>;
  const auto& binning = obs.getBinning();

  for (unsigned int bin=0, numBins = static_cast<unsigned int>(binning.numBins()); bin < numBins; ++bin) {
    const double low  = std::max(binning.binLow(bin), xlo);
    const double high = std::min(binning.binHigh(bin), xhi);
    const double width = high - low;

    // Check if this bin is in plotting range at all
    if (low >= high)
      continue;

    // Move support points slightly inside the bin, so step function is plotted correctly.
    binEdges->push_back(low  + 0.001 * width);
    binEdges->push_back(high - 0.001 * width);
  }

  return binEdges;
}


////////////////////////////////////////////////////////////////////////////////
/// Return reference to a unique_ptr holding the integrator that's used to sample the bins.
/// This can be used to change options such as sampling accuracy or to entirely exchange the integrator.
/// \see ROOT::Math::IntegratorOneDim::SetOptions to alter integration options.
/// \note When integration options are altered, these changes are not saved to files.
std::unique_ptr<ROOT::Math::IntegratorOneDim>& RooBinSamplingPdf::integrator() const {
  if (!_integrator) {
    _integrator.reset(new ROOT::Math::IntegratorOneDim(*this,
        ROOT::Math::IntegrationOneDim::kADAPTIVE, // GSL Integrator. Will really get it only if MathMore enabled.
        -1., _relEpsilon, // Abs epsilon = default, rel epsilon set by us.
        0, // We don't limit the sub-intervals. Steer run time via _relEpsilon.
        2 // This should be ROOT::Math::Integration::kGAUSS21, but this is in MathMore, so cannot include.
        ));
  }

  return _integrator;
}


////////////////////////////////////////////////////////////////////////////////
/// Binding used by the integrator to evaluate the PDF.
double RooBinSamplingPdf::operator()(double x) const {
  _observable->setVal(x);
  return _pdf->getVal(_normSetForIntegrator);
}


////////////////////////////////////////////////////////////////////////////////
/// Integrate the wrapped PDF using our current integrator, with given norm set and limits.
double RooBinSamplingPdf::integrate(const RooArgSet* normSet, double low, double high) const {
  // Need to set this because operator() only takes one argument.
  _normSetForIntegrator = normSet;

  return integrator()->Integral(low, high);
}

