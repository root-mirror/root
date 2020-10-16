// Tests for the RooRealBinding, and numeric integrators that use it.
// Author: Stephan Hageboeck, CERN  05/2020

#include "RooRealBinding.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooIntegrator1D.h"
#include "RooNumIntConfig.h"
#include "RooHelpers.h"
#include <Math/ProbFuncMathCore.h>
#include <Math/SpecFuncMathCore.h>
#include <Math/Integrator.h>
#include <TMath.h>
#include <TStopwatch.h>

#include "../src/RooFitLegacy/OldRooIntegrator1D.h"

#include "gtest/gtest.h"

#include <numeric>
#include <algorithm>

///
TEST(RooRealBinding, BatchEvalFeature) {
  RooRealVar a("a", "a", -100, 100);
  RooRealVar b("b", "b", -100, 100);
  RooFormulaVar formula("formula", "1.3*a + 1.4*b", RooArgList(a, b));

  std::vector<double> as;
  std::vector<double> bs;
  std::generate_n(std::back_inserter(as), 10, [](){ static double val = 0; return val += 0.3;});
  std::generate_n(std::back_inserter(bs), 10, [](){ static double val = 0; return val += 0.4;});

  std::vector<RooSpan<const double>> data;
  data.emplace_back(as);
  data.emplace_back(bs);

  RooRealBinding binding(formula, RooArgSet(a, b));
  auto result = binding.getValues(data);
  for (unsigned int i=0; i < 10; ++i) {
    EXPECT_DOUBLE_EQ(result[i], 1.3 * as[i] + 1.4 * bs[i]) << "result[" << i << "] a=" << as[i] << " b=" << bs[i];
  }
}


TEST(Roo1DIntegrator, RunFormulaVar_Trapezoid) {

  RooRealVar x("x", "x", -100, 100);
  RooRealVar a("a", "a", 0.2, -100, 100);
  RooRealVar b("b", "b", 0.3, -100, 100);
  RooFormulaVar formula("formula", "0.1 + x*(a + b*x)", RooArgList(x, a, b));
  auto solution = [&x](double a, double b){
    auto indefInt = [=](double y){
      return y*(0.1 + y*(1./2.*a + 1./3.*b * y));
    };
    return indefInt(x.getMax()) - indefInt(x.getMin());
  };
  RooRealBinding binding(formula, RooArgSet(x, a, b));

  // The integrators will warn, since we let them run until maxSteps
  RooHelpers::HijackMessageStream hijack(RooFit::WARNING, RooFit::Integration);

  // Test the recursion anchors of the Romberg integration
  {
    RooIntegrator1D oneStep(binding, RooIntegrator1D::Trapezoid, 1, 1.E-15);
    EXPECT_DOUBLE_EQ(oneStep.integral(), 0.5*200.*(2*0.1 + 2.*0.3*10000.));
    x = -100.;
    const double left = formula.getVal();
    x = 100.;
    const double right = formula.getVal();
    // Run integral again, also to make sure that setting x has no effect:
    EXPECT_DOUBLE_EQ(oneStep.integral(), 0.5*200.*(left + right));

    RooIntegrator1D twoStep(binding, RooIntegrator1D::Trapezoid, 2, 1.E-15);
    x = 0.;
    const double middle = formula.getVal();
    twoStep.applySeriesAcceleration(false);
    const double noAccel = twoStep.integral();
    EXPECT_DOUBLE_EQ(noAccel, 0.25*200.*(left + right) + 0.5*200.*middle);

    twoStep.applySeriesAcceleration(true);
    const double accel = twoStep.integral();
    EXPECT_LT(fabs(accel - solution(a.getVal(), b.getVal())), 0.8 * fabs(noAccel - solution(a.getVal(), b.getVal())))
        << "Expect with acceleration to be better than without.";
  }

  // Now run many steps
  {
    constexpr unsigned int nSteps = 25;
    constexpr double relEps = 1.E-50;
    RooIntegrator1D integrator(binding, RooIntegrator1D::Trapezoid, nSteps, relEps);
    double inputs[] = {1., 3.123};
    double target = solution(1., 3.123);
    EXPECT_LT(fabs(integrator.integral(inputs) - target)/target, 1.E-14);

    target = solution(a.getVal(), b.getVal());
    EXPECT_LT(fabs(integrator.integral() - target)/target, 1.E-14);

    inputs[0] = 4.; inputs[1] = 5.;
    target = solution(4., 5.);
    EXPECT_LT(fabs(integrator.integral(inputs) - target)/target, 1.E-14);
  }
}


TEST(Roo1DIntegrator, RunQuarticFormulaVar) {
  constexpr unsigned int nSteps = 25;
  constexpr double relEps = 1.E-50;
  RooRealVar x("x", "x", -50, 50);
  RooRealVar a("a", "a", 0.2, -100, 100);
  RooRealVar b("b", "b", 0.3, -100, 100);
  RooRealVar c("c", "c", 0.4, -100, 100);
  RooRealVar d("d", "d", 0.5, -100, 100);
  RooFormulaVar formula("formula", "0.1 + x*(a + x*(b + x*(c + d * x)))", RooArgList(x, a, b, c, d));
  auto solution = [&x](double a, double b, double c, double d){
    auto indefInt = [=](double y){
      return y*(0.1 + y*(1./2.*a + y*(1./3.*b + y*(1./4.*c + 1./5.*d * y))));
    };
    return indefInt(x.getMax()) - indefInt(x.getMin());
  };
  RooRealBinding binding(formula, RooArgSet(x, a, b, c, d));
  RooIntegrator1D integrator(binding, RooIntegrator1D::Trapezoid, nSteps, relEps);

  double target = solution(0.2, 0.3, 0.4, 0.5);
  EXPECT_LT(fabs(integrator.integral() - target)/target, 1.E-13);
}


/// Run numeric integrations using RooIntegrator1D and ROOT's adaptive integrator. Ensure that
/// they reach the requested precision.
void testConvergenceSettings(const RooFormulaVar& formula, const RooArgSet& observables, const RooArgSet& funcParameters,
    std::function<double(const double* pars, unsigned int nPar)> solution, const std::string& name) {
  constexpr unsigned int nSteps = 25;
  constexpr bool printDetails = false;
  SCOPED_TRACE(std::string("Function to integrate: ") + name + "\t" + formula.GetTitle());

  for (double relEps : {1.E-3, 1.E-6, 1.E-8}) {
    SCOPED_TRACE(std::string("relEps=" + std::to_string(relEps)));
    TStopwatch stopOld;
    TStopwatch stopNew;
    TStopwatch stopRoot;

    RooRealVar& x = dynamic_cast<RooRealVar&>(*observables[0ul]);
    RooArgSet variables(observables);
    variables.add(funcParameters);
    RooRealBinding binding(formula, variables);
    RooIntegrator1D integrator(binding, RooIntegrator1D::Trapezoid, nSteps, relEps);
    OldRooIntegrator1D oldIntegrator(binding, OldRooIntegrator1D::Trapezoid, nSteps, relEps);

    std::vector<double> errors;
    std::vector<double> errorsOld;
    std::vector<double> errorsRootInt;

    RooArgSet initialParameters;
    funcParameters.snapshot(initialParameters);

    constexpr unsigned int nRun = 10000;
    for (unsigned int i=0; i < nRun; ++i) {
      std::vector<double> pars;
      for (const auto p : initialParameters) {
        auto par = static_cast<RooRealVar*>(p);
        pars.push_back(par->getVal() + (par->getMax()-par->getVal())/nRun * i);
      }

      const double target = solution(pars.data(), pars.size());

      stopNew.Start(false);
      const double integral = integrator.integral(pars.data());
      stopNew.Stop();
      errors.push_back(fabs((integral - target)/target));

      stopOld.Start(false);
      const double oldIntegral = oldIntegrator.integral(pars.data());
      stopOld.Stop();
      errorsOld.push_back(fabs((oldIntegral - target)/target));

      {
        std::vector<double> vars(1);
        std::copy(pars.begin(), pars.end(), std::back_inserter(vars));
        auto bindingROOTInt = [&vars,&binding](double x){
          vars[0] = x;
          return binding(vars.data());
        };
        ROOT::Math::IntegratorOneDim rootIntegrator(bindingROOTInt, ROOT::Math::IntegrationOneDim::kADAPTIVE, 1.E-20, relEps, 100);

        stopRoot.Start(false);
        const double rootIntegral = rootIntegrator.Integral(x.getMin(), x.getMax());
        stopRoot.Stop();
        errorsRootInt.push_back(fabs((rootIntegral - target)/target));
      }
    }

    auto mean = [](const std::vector<double>& vec) {
      return TMath::Mean(vec.begin(), vec.end());
    };
    auto stdDev = [&](const std::vector<double>& vec) {
      return TMath::StdDev(vec.begin(), vec.end());
    };
    auto variance = [&](const std::vector<double>& vec) {
      auto tmp = TMath::StdDev(vec.begin(), vec.end());
      return tmp*tmp;
    };
    auto median = [&](const std::vector<double>& vec) {
      return TMath::Median(vec.size(), vec.data());
    };
    auto max = [&](const std::vector<double>& vec) {
      return TMath::MaxElement(vec.size(), vec.data());
    };
    auto q95_99 = [&](const std::vector<double>& vec) -> std::pair<double, double>{
      std::vector<double> q(2);
      std::vector<double> p(1, 0.95);
      p.push_back(0.99);
      TMath::Quantiles(vec.size(), q.size(), const_cast<double*>(vec.data()), q.data(), p.data(), false);
      return {q[0], q[1]};
    };

    const auto q95_99_int1D = q95_99(errors);
    const auto q95_99_intROOT = q95_99(errorsRootInt);

    if (printDetails) {
      std::cout << "Integrating " << name << ", relEps = " << relEps;
      std::cout << "\n\t    \t" << "mean         \tmedian  \tq95     \tq99     \tmax";
      const std::vector<double>* vec = &errorsOld;
      std::cout << "\n\told:\t" << mean(*vec) << "\t" << median(*vec)
          << "\t" << q95_99(*vec).first << "\t" << q95_99(*vec).second << "\t" << max(errorsOld) << "\tt=" << stopOld.CpuTime();
      vec = &errors;
      std::cout << "\n\tnew:\t" << mean(*vec) << "\t" << median(*vec)
          << "\t" << q95_99_int1D.first << "\t" << q95_99_int1D.second << "\t" << max(errorsOld) << "\tt=" << stopNew.CpuTime();
      vec = &errorsRootInt;
      std::cout << "\n\tROOT:\t" << mean(*vec) << "\t" << median(*vec)
          << "\t" << q95_99_intROOT.first << "\t" << q95_99_intROOT.second << "\t" << max(errorsOld) << "\tt=" << stopRoot.CpuTime();
      std::cout << std::endl;

      if (mean(errorsOld) > 10. * relEps) {
        std::cerr << "Old integrator reached poor precision: " << mean(errorsOld) << " +/- " << stdDev(errorsOld)
            << " instead of " << relEps << std::endl;
      }
    }

    EXPECT_LT(mean(errors), relEps) << "RooIntegrator1D should reach target precision.";
    EXPECT_LT(mean(errorsRootInt), relEps) << "ROOT integrator should reach target precision.";

    double mult = 1.;
    if (strstr(::testing::UnitTest::GetInstance()->current_test_info()->name(), "breitWig") && relEps > 1.E-4) {
      // Both integrators have trouble with the Breit-Wigner distribution and high relEpsilon
      // For this one config, we are more forgiving:
      mult = 6.5;
    }
    EXPECT_LT(q95_99_int1D.first, mult * relEps) << "95% quantile of errors exceeds relEpsilon";
    EXPECT_LT(q95_99_int1D.second, 2. * mult * relEps) << "99% quantile of errors exceeds 2.*relEpsilon";

    EXPECT_LT(q95_99_intROOT.first, mult * relEps) << "95% quantile of errors exceeds relEpsilon";
    EXPECT_LT(q95_99_intROOT.second, 2. * mult * relEps) << "95% quantile of errors exceeds relEpsilon";

    if (mean(errorsOld) > 1.E-9) {
      EXPECT_LT(mean(errors), mean(errorsOld)) << "New integrator should be better than old.";
      EXPECT_LT(variance(errors), variance(errorsOld)) << "New integrator should be more stable than old.";
    }
  }
}

TEST(Roo1DIntegrator, ConvergenceSettings_log) {
  RooRealVar x("x", "x", 0.1, 50);
  RooRealVar a("a", "a", 0.2, -100, 1.E5);
  RooFormulaVar formula("formula", "log(a*x)", RooArgList(x, a));
  testConvergenceSettings(formula,
      RooArgSet(x),
      RooArgSet(a),
      [&x](const double* pars, unsigned int nPar){
          const double a = pars[0];
          if (nPar != 1)
            throw std::logic_error("Need npar == 1");

          auto indefInt = [=](double y){
            return 1./a * (y * log(y) - y);
          };
          return indefInt(a*x.getMax()) - indefInt(a*x.getMin());
      },
      "log(a*x)");
}

TEST(Roo1DIntegrator, ConvergenceSettings_pol4) {
  RooRealVar x2("x", "x", -10, 10);
  RooRealVar a2("a", "a",  0.2, -1., 1);
  RooRealVar b2("b", "b",  -0.1, -1., 1);
  RooRealVar c2("c", "c", 0.02, -0.1, 0.1);
  RooRealVar d2("d", "d", -0.,  -0.01, 0.01);
  RooFormulaVar formula("formula", "0.1 + x*(a + x*(b + x*(c + d * x)))", RooArgList(x2, a2, b2, c2, d2));
  testConvergenceSettings(formula,
      RooArgSet(x2),
      RooArgSet(a2, b2, c2, d2),
      [&x2](const double* pars, unsigned int nPar){
          const double a = pars[0];
          const double b = pars[1];
          const double c = pars[2];
          const double d = pars[3];
          if (nPar != 4)
            throw std::logic_error("Need npar == 4");

          auto indefInt = [=](double y){
            return y*(0.1 + y*(1./2.*a + y*(1./3.*b + y*(1./4.*c + 1./5.*d * y))));
          };
          return indefInt(x2.getMax()) - indefInt(x2.getMin());
      },
      "Polynomial 4th order");
}

TEST(Roo1DIntegrator, ConvergenceSettings_breitWig) {
  RooRealVar x3("x", "x", 0., 100.);
  RooRealVar a3("a", "a", 10., 100.);
  RooRealVar b3("b", "b", 3., 10.);
  RooFormulaVar formula("breitwigner", "ROOT::Math::breitwigner_pdf(x, b, a)", RooArgList(x3, a3, b3));
  testConvergenceSettings(formula,
      RooArgSet(x3),
      RooArgSet(a3, b3),
      [&x3](const double* pars, unsigned int nPar) {
        const double a = pars[0];
        const double b = pars[1];
        if (nPar != 2)
          throw std::logic_error("Need npar == 2");

        return ROOT::Math::breitwigner_cdf(x3.getMax(), b, a) - ROOT::Math::breitwigner_cdf(x3.getMin(), b, a);
      },
      "Breit-Wigner distribution");
}

TEST(Roo1DIntegrator, ConvergenceSettings_Erf) {
  RooRealVar x("x", "x", -10., 10.);
  RooRealVar m("m", "m", -5., 5.);
  RooRealVar s("s", "s", 0.5, 10.);
  const std::string funcStr = "ROOT::Math::gaussian_pdf(x, s, m)";
  RooFormulaVar formula("gaus", funcStr.c_str(), RooArgSet(x, m, s));
  testConvergenceSettings(formula,
      RooArgSet(x),
      RooArgSet(m, s),
      [&x](const double* pars, unsigned int nPar) {
        const double m = pars[0];
        const double s = pars[1];
        if (nPar != 2)
          throw std::logic_error("Need npar == 2");

        return ROOT::Math::gaussian_cdf(x.getMax(), s, m) - ROOT::Math::gaussian_cdf(x.getMin(), s, m);
      },
      "Gaussian distribution");
}

TEST(Roo1DIntegrator, RunVsOldIntegrator) {
  constexpr unsigned int nSteps = 25;
  constexpr double relEps = 1.E-50;
  RooRealVar x("x", "x", -100, 100);
  RooRealVar a("a", "a", 0.2, -100, 100);
  RooRealVar b("b", "b", 0.3, -100, 100);

  RooFormulaVar formula("formula", "0.1 + x*(a + b*x)", RooArgList(x, a, b));
  auto solution = [&x](double a, double b){
    auto indefInt = [=](double y){
      return y*(0.1 + y*(1./2.*a + 1./3.*b * y));
    };
    return indefInt(x.getMax()) - indefInt(x.getMin());
  };
  RooRealBinding binding(formula, RooArgSet(x, a, b));

  RooIntegrator1D integrator(binding, RooIntegrator1D::Trapezoid, nSteps, relEps);
  OldRooIntegrator1D old1D(binding, OldRooIntegrator1D::Trapezoid, nSteps, relEps);

  double inputs[2];
  inputs[0] = 0.2; inputs[1] = 0.3;
  a = 0.2;
  b = 0.3;
  double target = solution(0.2, 0.3);
  EXPECT_LE(fabs(integrator.integral(inputs) - target), fabs(old1D.integral(inputs) - target));

  target = solution(4.4, 5.5);
  a = 4.4;
  b = 5.5;
  EXPECT_LE(fabs(integrator.integral() - target), fabs(old1D.integral() - target));
}


TEST(Roo1DIntegrator, RunErf_NStep) {
  RooHelpers::HijackMessageStream hijack(RooFit::WARNING, RooFit::Integration);
  constexpr double sigma = 1.3;
  constexpr double mean  = 0.2;

  for (auto limits : std::initializer_list<std::pair<double, double>>{ {0.1, 1.5}, {-0.3, 2.}, {2.5, 4.5} } ) {
    const double min = limits.first;
    const double max = limits.second;
    RooRealVar theX("theX", "x", min, max);
    std::string funcStr = std::string("ROOT::Math::gaussian_pdf(theX, ") + std::to_string(sigma) + ", " + std::to_string(mean) + ")";
    RooFormulaVar gaus("gaus", funcStr.c_str(), theX);
    RooRealBinding binding(gaus, theX);

    double targetError = 0.;
    for (unsigned int nSteps = 4; nSteps < 24; ++nSteps) {
      RooIntegrator1D integrator(binding, RooIntegrator1D::Trapezoid, nSteps, 1.E-17);
      const double integral = integrator.integral();

      OldRooIntegrator1D oldIntegrator(binding, OldRooIntegrator1D::Trapezoid, nSteps, 1.E-17);
      const double oldIntegral = oldIntegrator.integral();

      const double error = fabs(integral - (ROOT::Math::gaussian_cdf(max, sigma, mean) - ROOT::Math::gaussian_cdf(min, sigma, mean)));
      const double oldError = fabs(oldIntegral - (ROOT::Math::gaussian_cdf(max, sigma, mean) - ROOT::Math::gaussian_cdf(min, sigma, mean)));

      if (oldError != 0.)
        EXPECT_LT(error, oldError);

      if (nSteps == 4) {
        targetError = error;
      } else {
        // Error should go down faster than 0.5^nSteps because the integrator uses series acceleration,
        // so test if it goes down faster than 0.333^nSteps
        targetError /= 3.;
        // But cannot be better than double precision
        EXPECT_LT(error, std::max(targetError, 1.E-16) )    << "For step " << nSteps << " with integral=" << integral;
      }
      if (nSteps > 10)
        EXPECT_LT(error / integral, 1.E-4) << "For step " << nSteps << " with integral=" << integral;
      if (nSteps > 15)
        EXPECT_LT(error / integral, 1.E-6) << "For step " << nSteps << " with integral=" << integral;
      if (nSteps > 21)
        EXPECT_LT(error / integral, 1.E-8) << "For step " << nSteps << " with integral=" << integral;
    }
  }
}


TEST(Roo1DIntegrator, RunErf_Midpoint) {
  const double min=0, max=1;
  RooRealVar theX("theX", "x", min, max);
  RooFormulaVar gaus("gaus", "ROOT::Math::gaussian_pdf(theX, 1, 0)", theX);
  RooRealBinding binding(gaus, theX);
  double targetError = 0.;

  RooHelpers::HijackMessageStream hijack(RooFit::WARNING, RooFit::Integration);

  for (unsigned int nSteps = 4; nSteps < 20; ++nSteps) {
    RooIntegrator1D integrator(binding, RooIntegrator1D::Midpoint, nSteps, 1.E-16);
    const double integral = integrator.integral();
    const double error = fabs(integral - (ROOT::Math::gaussian_cdf(max, 1, 0) - ROOT::Math::gaussian_cdf(min, 1, 0)));
    if (nSteps == 4) {
      targetError = error;
    } else {
      // Error should go down faster than 2^nSteps because of series acceleration.
      targetError /= 3.;
      // But cannot be better than double precision
      EXPECT_LT(error, std::max(targetError, 1.E-16) )    << "For step " << nSteps << " with integral=" << integral;
    }
    if (nSteps > 10)
      EXPECT_LT(error / integral, 1.E-4) << "For step " << nSteps << " with integral=" << integral;
    if (nSteps > 15)
      EXPECT_LT(error / integral, 1.E-6) << "For step " << nSteps << " with integral=" << integral;
  }
}
