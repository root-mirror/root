/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// This macro demonstrates how to set up a fit in two ranges
/// such that it does not only fit the shapes in each region, but also
/// takes into account the relative normalization of the two.
///
/// ### 1. Shape fits (plain likelihood)
///
/// If you perform a fit in two ranges in RooFit, e.g. `pdf->fitTo(data,Range("Range1,Range2"))`
/// it will construct a simple simultaneous fit of the two regions.
///
/// In case the pdf is not extended, i.e., a shape fit, it will only fit the shapes in the
/// two selected ranges, and not take into account the relative predicted yields in those ranges.
///
/// In certain models (like exponential decays) and configurations (e.g. narrow ranges that are far apart),
/// the relative normalization of the ranges may carry much more information about the function parameters
/// than the shape of the distribution inside those ranges. Therefore, it is important to take that into
/// account.
///
/// This is particularly important for cases where the 2-range fit is meant to be representative of
/// a full-range fit, but with a blinded signal region inside it.
///
///
/// ### 2. Shape+rate fits (extended likelihood)
///
/// Also if your pdf is already extended, i.e. measuring both the distribution in the observable as well
/// as the event count in the fitted region, some intervention is needed to make fits in ranges
/// work in a way that corresponds to intuition.
///
/// If an extended fit is performed in a sub-range, the observed yield is only that of the subrange, hence
/// the expected event count will converge to a number that is smaller than what's visible in a plot.
/// In such cases, it is often preferred to interpret the extended term with respect to the full range
/// that's plotted, i.e., apply a correction to the extended likelihood term in such a way
/// that the interpretation of the expected event count remains that of the full range. This can
/// be done by applying a correcion factor (equal to the fraction of the pdf that is contained in the
/// fitted range) in the Poisson term that represents the extended likelihood term.
///
/// If an extended likelihood fit is performed over *two* sub-ranges, this correction is
/// even more important: without it, each component likelihood would have a different interpretation
/// of the expected event count (each corresponding to the count in its own region), and a joint
/// fit of these regions with different interpretations of the same model parameter results
/// in a number that is not easily interpreted.
///
/// If both regions correct their interpretatin such that N_expected refers to the full range,
/// it is interpreted easily, and consistent in both regions.
///
/// This requires that the likelihood model is extended using RooAddPdf in the
/// form SumPdf = Nsig * sigPdf + Nbkg * bkgPdf.
///
/// \macro_image
/// \macro_code
/// \macro_output
///
/// \authors Stephan Hageboeck, Wouter Verkerke

#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooExtendPdf.h"
#include "TCanvas.h"
using namespace RooFit;

void rf204b_extendedLikelihood_rangedFit()
{

 // PART 1: Background-only fits
 // ----------------------------

 // Build plain exponential model
 RooRealVar x("x", "x", 10, 100);
 RooRealVar alpha("alpha", "alpha", -0.04, -0.1, -0.0);
 RooExponential model("model", "Exponential model", x, alpha);

 // Define side band regions and full range
 x.setRange("LEFT",10,20);
 x.setRange("RIGHT",60,100);

 x.setRange("FULL",10,100);

 RooDataSet* data = model.generate(x, 10000);

 // Construct an extended pdf, which measures the event count N **on the full range**.
 // If the actual domain of x that is fitted is identical to FULL, this has no affect.
 //
 // If the fitted domain is a subset of `FULL`, though, the expected event count is divided by
 // \f[
 //   \mathrm{frac} = \frac{
 //     \int_{\mathrm{Fit range}} \mathrm{model}(x)  \; \mathrm{d}x }{
 //     \int_{\mathrm{Full range}} \mathrm{model}(x) \; \mathrm{d}x }.
 // \f]
 // `N` will therefore return the count extrapolated to the full range instead of the fit range.
 //
 // **Note**: When using a RooAddPdf for extending the likelihood, the same effect can be achieved with
 // [RooAddPdf::fixCoefRange()](https://root.cern.ch/doc/master/classRooAddPdf.html#ab631caf4b59e4c4221f8967aecbf2a65),

 RooRealVar N("N", "Extended term", 0, 20000);
 RooExtendPdf extmodel("extmodel", "Extended model", model, N, "FULL");


 // It can be instructive to fit the above model to either the LEFT or RIGHT range. `N` should approximately converge to the expected number
 // of events in the full range. One may try to leave out `"FULL"` in the constructor, o the the interpretation of `N` changes.
 extmodel.fitTo(*data, Range("LEFT"), PrintLevel(-1));
 N.Print();


 // If we now do a simultaneous fit to the extended model, instead of the original model, the LEFT and RIGHT range will each correct their local
 // `N` such that it refers to the `FULL` range.
 //
 // This joint fit of the extmodel will include (w.r.t. the plain model fit) a product of extended terms
 // \f[
 //     L_\mathrm{ext} = L
 //       \cdot \mathrm{Poisson} \left( N_\mathrm{obs}^\mathrm{LEFT}  | N_\mathrm{exp} / \mathrm{frac LEFT} \right)
 //       \cdot \mathrm{Poisson} \left( N_\mathrm{obs}^\mathrm{RIGHT} | N_\mathrm{exp} / \mathrm{frac RIGHT} \right)
 // \f]
 // that will introduce additional sensitivity of the likelihood to the slope parameter alpha of the exponential model through the `frac_LEFT` and `frac_RIGHT` integrals.
 //
 // In the extreme case of an exponential function and a fit in narrow LEFT and RIGHT ranges, this sensitivity may actually be larger
 // than from the shapes.
 //
 // This is also nicely demonstrated in the example below where the uncertainty on alpha is almost 5x smaller if the extended term is included.


 TCanvas* c = new TCanvas("c", "c", 2100, 700);
 c->Divide(3);
 c->cd(1);

 RooFitResult* r = model.fitTo(*data, Range("LEFT,RIGHT"), Save());

 RooPlot* frame = x.frame();
 data->plotOn(frame);
 model.plotOn(frame, VisualizeError(*r));
 model.plotOn(frame);
 model.paramOn(frame, Label("Bkg fit. Large errors since\nnormalisation ignored"));
 frame->Draw();

 c->cd(2);

 RooFitResult* r2 = extmodel.fitTo(*data, Range("LEFT,RIGHT"), Save());
 RooPlot* frame2 = x.frame();
 data->plotOn(frame2);
 extmodel.plotOn(frame2);
 extmodel.plotOn(frame2, VisualizeError(*r2));
 extmodel.paramOn(frame2, Label("Bkg fit. Normalisation\nincluded"), Layout(0.4,0.95));
 frame2->Draw();

 // PART 2: Extending with RooAddPdf
 // --------------------------------
 //
 // Now we repeat the above exercise, but instead of explicitly adding an extended term to a single shape pdf (RooExponential),
 // we assume that we already have an extended likelihood model in the form of a RooAddPdf constructed in the form `Nsig * sigPdf + Nbkg * bkgPdf`.
 //
 // We add a Gaussian to the previously defined exponential background.
 // The signal shape parameters are chosen constant, since the signal region is entirely in the blinded region, i.e., the fit has no sensitivity.

 c->cd(3);

 RooRealVar Nsig("Nsig", "Number of signal events", 1000, 0, 2000);
 RooRealVar Nbkg("Nbkg", "Number of background events", 10000, 0, 20000);

 RooRealVar mean("mean", "Mean of signal model", 40.);
 RooRealVar width("width", "Width of signal model", 5.);
 RooGaussian sig("sig", "Signal model", x, mean, width);

 RooAddPdf modelsum("modelsum", "NSig*signal + NBkg*background", RooArgSet(sig, model), RooArgSet(Nsig, Nbkg));

 // This model will automatically insert the correction factor for the reinterpretation of Nsig and Nnbkg in the full ranges.
 //
 // When this happens, it reports this with lines like the following:
 // ```
 // [#1] INFO:Fitting -- RooAbsOptTestStatistic::ctor(nll_modelsum_modelsumData_LEFT) fixing interpretation of coefficients of any RooAddPdf to full domain of observables
 // [#1] INFO:Fitting -- RooAbsOptTestStatistic::ctor(nll_modelsum_modelsumData_RIGHT) fixing interpretation of coefficients of any RooAddPdf to full domain of observables
 // ```

 RooFitResult* r3 = modelsum.fitTo(*data, Range("LEFT,RIGHT"), Save());

 RooPlot* frame3 = x.frame();
 data->plotOn(frame3);
 modelsum.plotOn(frame3);
 modelsum.plotOn(frame3, VisualizeError(*r3));
 modelsum.paramOn(frame3, Label("S+B fit with RooAddPdf"), Layout(0.3,0.95));
 frame3->Draw();

 c->Draw();
}
