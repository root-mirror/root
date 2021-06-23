/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   IP, Inti Pelupessy,  NL eScience Center, i.pelupessy@esciencecenter.nl  *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOFIT_MULTIPROCESS_NLLVAR_H
#define ROOFIT_MULTIPROCESS_NLLVAR_H

#include <map>
#include <MultiProcess/Vector.h>
#include <RooNLLVar.h>

namespace RooFit {
  namespace MultiProcessV1 {

    enum class NLLVarTask {
      all_events,
      single_event,
      bulk_partition,
      interleave
    };

    // for debugging:
    std::ostream& operator<<(std::ostream& out, const NLLVarTask value);


    // --- kahan summation templates ---

    template <typename C>
    typename C::value_type sum_kahan(const C& container) {
      using ValueType = typename C::value_type;
      ValueType sum = 0, carry = 0;
      for (auto element : container) {
        ValueType y = element - carry;
        ValueType t = sum + y;
        carry = (t - sum) - y;
        sum = t;
      }
      return sum;
    }
    
    template <typename IndexType, typename ValueType>
    ValueType sum_kahan(const std::map<IndexType, ValueType>& map) {
      ValueType sum = 0, carry = 0;
      for (auto element : map) {
        ValueType y = element.second - carry;
        ValueType t = sum + y;
        carry = (t - sum) - y;
        sum = t;
      }
      return sum;
    }
    
    template <typename C>
    std::pair<typename C::value_type, typename C::value_type> sum_of_kahan_sums(const C& sum_values, const C& sum_carrys) {
      using ValueType = typename C::value_type;
      ValueType sum = 0, carry = 0;
      for (std::size_t ix = 0; ix < sum_values.size(); ++ix) {
        ValueType y = sum_values[ix];
        carry += sum_carrys[ix];
        y -= carry;
        const ValueType t = sum + y;
        carry = (t - sum) - y;
        sum = t;
      }
      return std::pair<ValueType, ValueType>(sum, carry);
    }
    
    
    
    class NLLVar : public RooFit::MultiProcessV1::Vector<RooNLLVar> {
     public:
      NLLVar(std::size_t NumCPU, NLLVarTask task_mode, const RooNLLVar& nll);
      void init_vars();
      void update_parameters();

      // the const is inherited from RooAbsTestStatistic::evaluate. We are not
      // actually const though, so we use a horrible hack.
      Double_t evaluate() const override;
      Double_t evaluate_non_const();
    
      // --- RESULT LOGISTICS ---
      void send_back_task_result_from_worker(std::size_t task) override;
      void receive_task_result_on_queue(std::size_t task, std::size_t worker_id) override;
      void send_back_results_from_queue_to_master() override;
      void clear_results() override;
      void receive_results_on_master() override;
    
     private:
      void evaluate_task(std::size_t task) override;
      double get_task_result(std::size_t /*task*/) override;

      // members
      std::map<std::size_t, double> carrys;
      double result = 0;
      double carry = 0;
      std::size_t N_tasks = 0;
      NLLVarTask mp_task_mode;
    };

  } // namespace MultiProcessV1
} // namespace RooFit

#endif //ROOFIT_MULTIPROCESS_NLLVAR_H
