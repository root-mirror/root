// Author: Danilo Piparo, Omar Zapata   16/12/2015

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <fcntl.h>
#ifdef _MSC_VER // Visual Studio
#include <winsock2.h>
#include <io.h>
#pragma comment(lib, "Ws2_32.lib")
#define pipe(fds) _pipe(fds, 1048575, _O_BINARY)
#define read _read
#define dup _dup
#define dup2 _dup2
#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define STDERR_FILENO 2
#else
#include <unistd.h>
#endif
#include <string>
#include <iostream>
#include "TInterpreter.h"

bool JupyROOTExecutorImpl(const char *code);
bool JupyROOTDeclarerImpl(const char *code);

class JupyROOTExecutorHandler {
private:
   bool fCapturing = false;
   std::string fStdoutpipe;
   std::string fStderrpipe;
   int fStdout_pipe[2] = {0,0};
   int fStderr_pipe[2] = {0,0};
   int fSaved_stderr = 0;
   int fSaved_stdout = 0;
public:
   JupyROOTExecutorHandler();
   void Poll();
   void InitCapture();
   void EndCapture();
   void Clear();
   std::string &GetStdout();
   std::string &GetStderr();
};


#ifndef F_LINUX_SPECIFIC_BASE
#define F_LINUX_SPECIFIC_BASE       1024
#endif
#ifndef F_SETPIPE_SZ
#define F_SETPIPE_SZ    (F_LINUX_SPECIFIC_BASE + 7)
#endif

constexpr long MAX_PIPE_SIZE = 1048575;

JupyROOTExecutorHandler::JupyROOTExecutorHandler() {}


static void PollImpl(FILE *stdStream, int *pipeHandle, std::string &pipeContent)
{
   int buf_read;
   char ch;
   fflush(stdStream);
   while (true) {
      buf_read = read(pipeHandle[0], &ch, 1);
      if (buf_read == 1) {
         pipeContent += ch;
      } else break;
   }
}

void JupyROOTExecutorHandler::Poll()
{
   PollImpl(stdout, fStdout_pipe, fStdoutpipe);
   PollImpl(stderr, fStderr_pipe, fStderrpipe);
}

static void InitCaptureImpl(int &savedStdStream, int *pipeHandle, int FILENO)
{
   savedStdStream = dup(FILENO);
   if (pipe(pipeHandle) != 0) {
      return;
   }
#ifdef _MSC_VER // Visual Studio
   unsigned long mode = 1;
   ioctlsocket(pipeHandle[0], FIONBIO, &mode);
#else
   long flags_stdout = fcntl(pipeHandle[0], F_GETFL);
   if (flags_stdout == -1) return;
   flags_stdout |= O_NONBLOCK;
   fcntl(pipeHandle[0], F_SETFL, flags_stdout);
   fcntl(pipeHandle[0], F_SETPIPE_SZ, MAX_PIPE_SIZE);
#endif
   dup2(pipeHandle[1], FILENO);
   close(pipeHandle[1]);
}

void JupyROOTExecutorHandler::InitCapture()
{
   if (!fCapturing)  {
      InitCaptureImpl(fSaved_stdout, fStdout_pipe, STDOUT_FILENO);
      InitCaptureImpl(fSaved_stderr, fStderr_pipe, STDERR_FILENO);
      fCapturing = true;
   }
}

void JupyROOTExecutorHandler::EndCapture()
{
   if (fCapturing)  {
      Poll();
      dup2(fSaved_stdout, STDOUT_FILENO);
      dup2(fSaved_stderr, STDERR_FILENO);
      fCapturing = false;
   }
}

void JupyROOTExecutorHandler::Clear()
{
   fStdoutpipe = "";
   fStderrpipe = "";
}

std::string &JupyROOTExecutorHandler::GetStdout()
{
   return fStdoutpipe;
}

std::string &JupyROOTExecutorHandler::GetStderr()
{
   return fStderrpipe;
}

JupyROOTExecutorHandler *JupyROOTExecutorHandler_ptr = nullptr;

////////////////////////////////////////////////////////////////////////////////
bool JupyROOTExecutorImpl(const char *code)
{
   auto status = false;
   try {
      auto err = TInterpreter::kNoError;
      if (gInterpreter->ProcessLine(code, &err)) {
         status = true;
      }

      if (err == TInterpreter::kProcessing) {
         gInterpreter->ProcessLine(".@");
         gInterpreter->ProcessLine("cerr << \"Unbalanced braces. This cell was not processed.\" << endl;");
      }
   } catch (...) {
      status = true;
   }

   return status;
}

bool JupyROOTDeclarerImpl(const char *code)
{
   auto status = false;
   try {
      if (gInterpreter->Declare(code)) {
         status = true;
      }
   } catch (...) {
      status = true;
   }
   return status;
}

extern "C" {

   int JupyROOTExecutor(const char *code)
   {
      return JupyROOTExecutorImpl(code);
   }
   int JupyROOTDeclarer(const char *code)
   {
      return JupyROOTDeclarerImpl(code);
   }

   void JupyROOTExecutorHandler_Clear()
   {
      JupyROOTExecutorHandler_ptr->Clear();
   }

   void JupyROOTExecutorHandler_Ctor()
   {
      if (!JupyROOTExecutorHandler_ptr) {
         JupyROOTExecutorHandler_ptr = new JupyROOTExecutorHandler();
         // Fixes for ROOT-7999
         gInterpreter->ProcessLine("SetErrorHandler((ErrorHandlerFunc_t)&DefaultErrorHandler);");
      }
   }

   void JupyROOTExecutorHandler_Poll()
   {
      JupyROOTExecutorHandler_ptr->Poll();
   }

   void JupyROOTExecutorHandler_EndCapture()
   {
      JupyROOTExecutorHandler_ptr->EndCapture();
   }

   void JupyROOTExecutorHandler_InitCapture()
   {
      JupyROOTExecutorHandler_ptr->InitCapture();
   }

   const char *JupyROOTExecutorHandler_GetStdout()
   {
      return JupyROOTExecutorHandler_ptr->GetStdout().c_str();
   }

   const char *JupyROOTExecutorHandler_GetStderr()
   {
      return JupyROOTExecutorHandler_ptr->GetStderr().c_str();
   }

   void JupyROOTExecutorHandler_Dtor()
   {
      if (!JupyROOTExecutorHandler_ptr) return;
      delete JupyROOTExecutorHandler_ptr;
      JupyROOTExecutorHandler_ptr = nullptr;
   }

}
