/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Likelihood and minimization: Recover from regions where the function is not defined.
///
/// We demonstrate improved recovery from disallowed parameters. For this, we use a polynomial PDF of the form
/// \f[
///   \mathrm{Pol2} = \mathcal{N} \left( c + a_1 \cdot x + a_2 \cdot x^2 + 0.01 \cdot x^3 \right),
/// \f]
/// where \f$ \mathcal{N} \f$ is a normalisation factor. Unless the parameters are chosen carefully,
/// this function can be negative, and hence, it cannot be used as a PDF. In this case, RooFit passes
/// an error to the minimiser, which might try to recover.
///
/// Before ROOT 6.24, RooFit always passed the highest function value that was encountered during the minimisation
/// to the minimiser. If a parameter is far in a disallowed region, the minimiser has to blindly test various
/// values of the parameters. It might find the correct values by chance, but often be unable to recover from bad
/// starting values. Here, we use a model with such bad values.
///
/// Starting with ROOT 6.24, the minimiser receives more information. For example, when a PDF is negative,
/// the magnitude of the "undershoot" is passed to the minimiser. The minimiser can use this to compute a
/// gradient, which will eventually lead it out of the disallowed region. The steepness of this gradient
/// can be chosen using RooFit::RecoverFromUndefinedRegions(double). A value of zero is equivalent to RooFit
/// before ROOT 6.24. Positive values activate the recovery. Values between 1. and 10. were found to be a
/// good default. If no argument is passed, RooFit uses 10.
///
/// \macro_output
/// \macro_code
///
/// \date 10/2020
/// \author Stephan Hageboeck

#include <RooRealVar.h>
#include <RooPolynomial.h>
#include <RooPlot.h>
#include <RooDataSet.h>
#include <RooGlobalFunc.h>
#include <RooFitResult.h>
#include <RooMsgService.h>

#include <TCanvas.h>
#include <TLegend.h>

void rf612_recoverFromInvalidParameters() {

  // Create a fit model:
  // The polynomial is notoriously unstable, because it can quickly go negative.
  // Since PDFs need to be positive, one often ends up with an unstable fit model.
  RooRealVar x("x", "x", -15, 15);
  RooRealVar a1("a1", "a1", -0.5, -10., 20.);
  RooRealVar a2("a2", "a2", 0.2, -10., 20.);
  RooRealVar a3("a3", "a3", 0.01);
  RooPolynomial pdf("pol3", "c + a1 * x + a2 * x*x + 0.01 * x*x*x", x, RooArgSet(a1, a2, a3));

  // Create toy data with all-positive coefficients:
  std::unique_ptr<RooDataSet> data(pdf.generate(x, 10000));

  // For plotting.
  // We create pointers to the plotted objects. We want these objects to leak out of the function,
  // so we can still see them after it returns.
  TCanvas* c = new TCanvas();
  RooPlot* frame = x.frame();
  data->plotOn(frame, RooFit::Name("data"));

  // Plotting a PDF with disallowed parameters doesn't work. We would get a lot of error messages.
  // Therefore, we disable plotting messages in RooFit's message streams:
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);


  // RooFit before ROOT 6.24
  // --------------------------------
  // Before 6.24, RooFit wasn't able to recover from invalid parameters. The minimiser just errs around
  // the starting values of the parameters without finding any improvement.

  // Set up the parameters such that the PDF would come out negative. The PDF is now undefined.
  a1.setVal(10.);
  a2.setVal(-1.);

  // Perform a fit:
  RooFitResult* fitWithoutRecovery = pdf.fitTo(*data, RooFit::Save(),
      RooFit::RecoverFromUndefinedRegions(0.), // This is how RooFit behaved prior to ROOT 6.24
      RooFit::PrintEvalErrors(-1), // We are expecting a lot of evaluation errors. -1 switches off printing.
      RooFit::PrintLevel(-1));

  pdf.plotOn(frame, RooFit::LineColor(kRed), RooFit::Name("noRecovery"));



  // RooFit since ROOT 6.24
  // --------------------------------
  // The minimiser gets information about the "badness" of the violation of the function definition. It uses this
  // to find its way out of the disallowed parameter regions.
  std::cout << "\n\n\n-------------- Starting second fit ---------------\n\n" << std::endl;

  // Reset the parameters such that the PDF is again undefined.
  a1.setVal(10.);
  a2.setVal(-1.);

  // Fit again, but pass recovery information to the minimiser:
  RooFitResult* fitWithRecovery = pdf.fitTo(*data, RooFit::Save(),
      RooFit::RecoverFromUndefinedRegions(1.), // The magnitude of the recovery information can be chosen here.
                                                // Higher values mean more aggressive recovery.
      RooFit::PrintEvalErrors(-1), // We are still expecting a few evaluation errors.
      RooFit::PrintLevel(0));

  pdf.plotOn(frame, RooFit::LineColor(kBlue), RooFit::Name("recovery"));



  // Collect results and plot.
  // --------------------------------
  // We print the two fit results, and plot the fitted curves.
  // The curve of the fit without recovery cannot be plotted, because the PDF is undefined if a2 < 0.
  fitWithoutRecovery->Print();
  std::cout << "Without recovery, the fitter encountered " << fitWithoutRecovery->numInvalidNLL()
      << " invalid function values. The parameters are unchanged." << std::endl;

  fitWithRecovery->Print();
  std::cout << "With recovery, the fitter encountered " << fitWithRecovery->numInvalidNLL()
      << " invalid function values, but the parameters are fitted." << std::endl;

  TLegend* legend = new TLegend(0.5, 0.7, 0.9, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(frame->findObject("data"), "Data", "P");
  legend->AddEntry(frame->findObject("noRecovery"), "Without recovery (cannot be plotted)", "L");
  legend->AddEntry(frame->findObject("recovery"), "With recovery", "L");
  frame->Draw();
  legend->Draw();
  c->Draw();
}

