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

#include <MultiProcess/messages.h>

namespace RooFit {
  namespace MultiProcessV1 {

    // for debugging
#define PROCESS_VAL(p) case(p): s = #p; break;

    std::ostream& operator<<(std::ostream& out, const M2Q value){
      const char* s = 0;
      switch(value){
        PROCESS_VAL(M2Q::terminate);
        PROCESS_VAL(M2Q::enqueue);
        PROCESS_VAL(M2Q::retrieve);
        PROCESS_VAL(M2Q::update_real);
        PROCESS_VAL(M2Q::call_double_const_method);
        PROCESS_VAL(M2Q::flush_ostreams);
      }
      return out << s;
    }

    std::ostream& operator<<(std::ostream& out, const Q2M value){
      const char* s = 0;
      switch(value){
        PROCESS_VAL(Q2M::retrieve_rejected);
        PROCESS_VAL(Q2M::retrieve_accepted);
      }
      return out << s;
    }

    std::ostream& operator<<(std::ostream& out, const W2Q value){
      const char* s = 0;
      switch(value){
        PROCESS_VAL(W2Q::dequeue);
        PROCESS_VAL(W2Q::send_result);
      }
      return out << s;
    }

    std::ostream& operator<<(std::ostream& out, const Q2W value){
      const char* s = 0;
      switch(value){
        PROCESS_VAL(Q2W::terminate);
        PROCESS_VAL(Q2W::dequeue_rejected);
        PROCESS_VAL(Q2W::dequeue_accepted);
        PROCESS_VAL(Q2W::update_real);
        PROCESS_VAL(Q2W::result_received);
        PROCESS_VAL(Q2W::call_double_const_method);
        PROCESS_VAL(Q2W::flush_ostreams);
      }
      return out << s;
    }

#undef PROCESS_VAL

  } // namespace MultiProcess
} // namespace RooFit
