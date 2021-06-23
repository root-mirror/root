/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOFIT_MULTIPROCESS_GRADMINIMIZER_H
#define ROOFIT_MULTIPROCESS_GRADMINIMIZER_H

#include <vector>
#include <MultiProcess/Vector.h>
#include <RooGradMinimizerFcn.h>
#include <RooMinimizer.h>

namespace RooFit {
  namespace MultiProcessV1 {
    class GradMinimizerFcn : public RooFit::MultiProcessV1::Vector<RooGradMinimizerFcn> {
     public:
      GradMinimizerFcn(RooAbsReal *funct, RooMinimizerGenericPtr context,
                       std::size_t _N_workers, bool verbose = false);
      GradMinimizerFcn(const GradMinimizerFcn& other);

      ROOT::Math::IMultiGradFunction* Clone() const override;

      void update_state();
      void update_real(std::size_t ix, double val, bool is_constant) override;

        // the const is inherited from ...::evaluate. We are not
      // actually const though, so we use a horrible hack.
//      Double_t evaluate() const override;
//      Double_t evaluate_non_const();

      // --- RESULT LOGISTICS ---
      void send_back_task_result_from_worker(std::size_t task) override;
      void receive_task_result_on_queue(std::size_t task, std::size_t worker_id) override;
      void send_back_results_from_queue_to_master() override;
      void clear_results() override;
      void receive_results_on_master() override;

      // overrides IGradientFunctionMultiDimTempl<double>::Gradient etc from
      // Math/IFunction.h, which are const, but we are not actually const, so
      // we use a const cast hack (the mutable versions):
      void Gradient(const double *x, double *grad) const override;
      void mutable_Gradient(const double *x, double *grad);
      void G2ndDerivative(const double *x, double *g2) const override;
      void mutable_G2ndDerivative(const double *x, double *g2);
      void GStepSize(const double *x, double *gstep) const override;
      void mutable_GStepSize(const double *x, double *gstep);

      void CalculateAll(const double *x);

     private:
      void evaluate_task(std::size_t task) override;
      double get_task_result(std::size_t /*task*/) override;

      // members
      std::size_t N_tasks = 0;
      std::vector<std::size_t> completed_task_ids;
    };

    using GradMinimizer = RooMinimizerTemplate<GradMinimizerFcn, RooFit::MinimizerType::Minuit2, std::size_t>;
  } // namespace MultiProcessV1
} // namespace RooFit

#endif //ROOFIT_MULTIPROCESS_GRADMINIMIZER_H
