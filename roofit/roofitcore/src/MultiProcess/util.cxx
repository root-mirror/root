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

#include <csignal> // kill, SIGKILL
#include <iostream> // cerr, and indirectly WNOHANG, EINTR, W* macros
#include <stdexcept> // runtime_error
#include <sys/wait.h>  // waitpid
#include <string>

#include <RooMsgService.h>

#include <MultiProcess/util.h>

namespace RooFit {
  namespace MultiProcessV1 {

    int wait_for_child(pid_t child_pid, bool may_throw, int retries_before_killing) {
      int status = 0;
      int patience = retries_before_killing;
      pid_t tmp;
      do {
        if (patience-- < 1) {
          ::kill(child_pid, SIGKILL);
        }
        tmp = waitpid(child_pid, &status, WNOHANG);
      } while (
          tmp == 0 // child has not yet changed state, try again
          || (-1 == tmp && EINTR == errno) // retry on interrupted system call
          );

      if (patience < 1) {
        ooccoutD(static_cast<RooAbsArg*>(nullptr),Eval) << "Had to send PID " << child_pid << " " << (-patience+1) << " SIGKILLs";
      }

      if (0 != status) {
        if (WIFEXITED(status)) {
          printf("exited, status=%d\n", WEXITSTATUS(status));
        } else if (WIFSIGNALED(status)) {
          printf("killed by signal %d\n", WTERMSIG(status));
        } else if (WIFSTOPPED(status)) {
          printf("stopped by signal %d\n", WSTOPSIG(status));
        } else if (WIFCONTINUED(status)) {
          printf("continued\n");
        }
      }

      if (-1 == tmp && may_throw) throw std::runtime_error(std::string("waitpid, errno ") + std::to_string(errno));

      return status;
    }

  }
}