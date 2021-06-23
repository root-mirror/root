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

#include <MultiProcess/NLLVar.h>
#include <MultiProcess/messages.h>
#include <RooAbsData.h>
#include <ROOT/RMakeUnique.hxx>  // make_unique
#include <tuple>  // std::tie

namespace RooFit {
  namespace MultiProcessV1 {

    std::ostream& operator<<(std::ostream& out, const NLLVarTask value) {
      const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
      switch(value){
        PROCESS_VAL(NLLVarTask::all_events);
        PROCESS_VAL(NLLVarTask::single_event);
        PROCESS_VAL(NLLVarTask::bulk_partition);
        PROCESS_VAL(NLLVarTask::interleave);
      }
#undef PROCESS_VAL
      return out << s;
    }


    NLLVar::NLLVar(std::size_t NumCPU, NLLVarTask task_mode, const RooNLLVar& nll) :
        RooFit::MultiProcessV1::Vector<RooNLLVar>(NumCPU, nll),  // uses copy constructor for the RooNLLVar part
        mp_task_mode(task_mode)
    {
      if (_gofOpMode == RooAbsTestStatistic::GOFOpMode::MPMaster) {
        TaskManager::remove_job_object(id);
        throw std::logic_error("Cannot create MPRooNLLVar based on a multi-CPU enabled RooNLLVar! The use of the BidirMMapPipe by MPFE in RooNLLVar conflicts with the use of BidirMMapPipe by MultiProcess classes.");
      }

      _vars = RooListProxy("vars", "vars", this);
      init_vars();
      switch (mp_task_mode) {
        case NLLVarTask::all_events: {
          N_tasks = 1;
          break;
        }
        case NLLVarTask::single_event: {
          N_tasks = static_cast<std::size_t>(_data->numEntries());
          break;
        }
        case NLLVarTask::bulk_partition:
        case NLLVarTask::interleave: {
          N_tasks = NumCPU;
          break;
        }
      }

      double_const_methods["getCarry"] = &NLLVar::getCarry;
    }

    void NLLVar::init_vars() {
      // Empty current lists
      _vars.removeAll() ;
      _saveVars.removeAll() ;

      // Retrieve non-constant parameters
      auto vars = std::make_unique<RooArgSet>(*getParameters(RooArgSet()));
      RooArgList varList(*vars);

      // Save in lists
      _vars.add(varList);
      _saveVars.addClone(varList);
    }

    void NLLVar::update_parameters() {
      if (get_manager()->is_master()) {
        for (std::size_t ix = 0u; ix < static_cast<std::size_t>(_vars.getSize()); ++ix) {
          bool valChanged = !_vars[ix].isIdentical(_saveVars[ix], kTRUE);
          bool constChanged = (_vars[ix].isConstant() != _saveVars[ix].isConstant());

          if (valChanged || constChanged) {
            if (constChanged) {
              ((RooRealVar *) &_saveVars[ix])->setConstant(_vars[ix].isConstant());
            }
            // TODO: Check with Wouter why he uses copyCache in MPFE; makes it very difficult to extend, because copyCache is protected (so must be friend). Moved setting value to if-block below.
            //          _saveVars[ix].copyCache(&_vars[ix]);

            // send message to queue (which will relay to workers)
            RooAbsReal * rar_val = dynamic_cast<RooAbsReal *>(&_vars[ix]);
            if (rar_val) {
              Double_t val = rar_val->getVal();
              dynamic_cast<RooRealVar *>(&_saveVars[ix])->setVal(val);
              RooFit::MultiProcessV1::M2Q msg = RooFit::MultiProcessV1::M2Q::update_real;
              Bool_t isC = _vars[ix].isConstant();
              get_manager()->send_from_master_to_queue(msg, id, ix, val, isC);
            }
            // TODO: implement category handling
            //            } else if (dynamic_cast<RooAbsCategory*>(var)) {
            //              M2Q msg = M2Q::update_cat ;
            //              UInt_t cat_ix = ((RooAbsCategory*)var)->getIndex();
            //              *_pipe << msg << ix << cat_ix;
            //            }
          }
        }
      }
    }

    Double_t NLLVar::evaluate() const {
      return const_cast<NLLVar*>(this)->evaluate_non_const();
    }

    Double_t NLLVar::evaluate_non_const() {
      if (get_manager()->is_master()) {
        // update parameters that changed since last calculation (or creation if first time)
        update_parameters();

        // master fills queue with tasks
        for (std::size_t ix = 0; ix < N_tasks; ++ix) {
          JobTask job_task(id, ix);
          get_manager()->to_queue(job_task);
        }
        waiting_for_queued_tasks = true;

        // wait for task results back from workers to master
        gather_worker_results();

        // put the results in vectors for calling sum_of_kahan_sums (TODO: make a map-friendly sum_of_kahan_sums)
        std::vector<double> results_vec, carrys_vec;
        for (auto const &item : results) {
          results_vec.emplace_back(item.second);
          carrys_vec.emplace_back(carrys[item.first]);
        }

        // sum task results
        std::tie(result, carry) = sum_of_kahan_sums(results_vec, carrys_vec);
      }
      return result;
    }

    // --- RESULT LOGISTICS ---

    void NLLVar::send_back_task_result_from_worker(std::size_t task) {
      result = get_task_result(task);
      carry = getCarry();
      get_manager()->send_from_worker_to_queue(id, task, result, carry);
    }

    void NLLVar::receive_task_result_on_queue(std::size_t task, std::size_t worker_id) {
      result = get_manager()->receive_from_worker_on_queue<double>(worker_id);
      carry = get_manager()->receive_from_worker_on_queue<double>(worker_id);
      results[task] = result;
      carrys[task] = carry;
    }

    void NLLVar::send_back_results_from_queue_to_master() {
      get_manager()->send_from_queue_to_master(results.size());
      for (auto const &item : results) {
        get_manager()->send_from_queue_to_master(item.first, item.second, carrys[item.first]);
      }
    }

    void NLLVar::clear_results() {
      // empty results caches
      results.clear();
      carrys.clear();
    }

    void NLLVar::receive_results_on_master() {
      std::size_t N_job_tasks = get_manager()->receive_from_queue_on_master<std::size_t>();
      for (std::size_t task_ix = 0ul; task_ix < N_job_tasks; ++task_ix) {
        std::size_t task_id = get_manager()->receive_from_queue_on_master<std::size_t>();
        results[task_id] = get_manager()->receive_from_queue_on_master<double>();
        carrys[task_id] = get_manager()->receive_from_queue_on_master<double>();
      }
    }

    // --- END OF RESULT LOGISTICS ---

    void NLLVar::evaluate_task(std::size_t task) {
      assert(get_manager()->is_worker());
      std::size_t N_events = static_cast<std::size_t>(_data->numEntries());
      // "default" values (all events in one task)
      std::size_t first = task;
      std::size_t last  = N_events;
      std::size_t step  = 1;
      switch (mp_task_mode) {
        case NLLVarTask::all_events: {
          // default values apply
          break;
        }
        case NLLVarTask::single_event: {
          last = task + 1;
          break;
        }
        case NLLVarTask::bulk_partition: {
          first = N_events * task / N_tasks;
          last  = N_events * (task + 1) / N_tasks;
          break;
        }
        case NLLVarTask::interleave: {
          step = N_tasks;
          break;
        }
      }

      result = evaluatePartition(first, last, step);
    }

    double NLLVar::get_task_result(std::size_t /*task*/) {
      // TODO: this is quite ridiculous, having a get_task_result without task
      // argument. We should have a cache, e.g. a map, that gives the result for
      // a given task. The caller (usually send_back_task_result_from_worker) can
      // then decide whether to erase the value from the cache to keep it clean.
      assert(get_manager()->is_worker());
      return result;
    }

  } // namespace MultiProcess
} // namespace RooFit
