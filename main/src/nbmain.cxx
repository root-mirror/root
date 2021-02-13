// @(#)root/main:$Id$
// Author: Enric Tejedor   07/10/15

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// NBMain                                                               //
//                                                                      //
// Main program used to spawn a ROOT notebook                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"

#include "TCollection.h"
#include "TList.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

#include <fstream>
#include <string>
#include <memory>

#define JUPYTER_CMD        "jupyter"
#define NB_OPT             "notebook"
#define JUPYTER_CONF_DIR_V "JUPYTER_CONFIG_DIR"
#define JUPYTER_PATH_V     "JUPYTER_PATH"
#define NB_CONF_DIR        "notebook"
#define ROOTNB_DIR         ".rootnb"
#define COMMIT_FILE        ".rootcommit"
#define JUPYTER_CONFIG     "jupyter_notebook_config.py"

using namespace std;

#ifdef WIN32
#include <process.h>
static string pathsep("\\");
#define execlp _execlp
#else
static string pathsep("/");
#endif

////////////////////////////////////////////////////////////////////////////////
/// Checks whether ROOT notebook files are installed and they are
/// the current version.

static int CheckNbInstallation(string dir)
{
   string commit(gROOT->GetGitCommit());
   string inputfname(dir + pathsep + ROOTNB_DIR + pathsep + COMMIT_FILE);
   ifstream in(inputfname);
   if (in.is_open()) {
      string line;
      in >> line;
      in.close();
      if (line.compare(commit) == 0) return  0; // already installed
      else                           return -1; // install, it's outdated
   }
   else if (gSystem->AccessPathName(inputfname.c_str())) {
      // There is no installation
      return -1;
   }
   else {
      fprintf(stderr,
              "Error checking notebook installation -- cannot open %s\n",
              inputfname.c_str());
      return -2;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Installs ROOT notebook files in the user's home directory.

static bool InstallNbFiles(string source, string dest)
{
   // Create installation directory
   if (gSystem->AccessPathName(dest.c_str())) {
      if (gSystem->mkdir(dest.c_str())) {
         fprintf(stderr,
                 "Error installing notebook configuration files -- cannot create directory %s\n",
                 dest.c_str());
         return false;
      }
   }

   // Copy files in source to dest
   TSystemDirectory dir(source.c_str(), source.c_str());
   std::unique_ptr<TList> files;
   files.reset(dir.GetListOfFiles());
   if (files) {
      TSystemFile *file;
      TListIter it(files.get());
      while ((file = (TSystemFile*)it())) {
         TString s = file->GetName();
         string fname(s.Data());
         string sourcefile = source + pathsep + fname;
         string destfile   = dest   + pathsep + fname;
         if (!file->IsDirectory()) {
            if (gSystem->CopyFile(sourcefile.c_str(), destfile.c_str(), true)) {
               fprintf(stderr,
                       "Error installing notebook configuration files -- cannot copy file %s to %s\n",
                       sourcefile.c_str(), destfile.c_str());
               return false;
            }
         }
         else if (fname.compare(".") && fname.compare("..") && fname.compare("html")) {
            if (!InstallNbFiles(sourcefile, destfile))
               return false;
         }
      }
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates the Jupyter notebook configuration file that sets the
/// necessary environment.

static bool CreateJupyterConfig(string dest, string rootbin, string rootlib, string rootdata)
{
   string jupyconfig = dest + pathsep + JUPYTER_CONFIG;
   ofstream out(jupyconfig, ios::trunc);
   if (out.is_open()) {
      out << "import os" << endl;
      out << "rootbin = '" << rootbin << "'" << endl;
      out << "rootlib = '" << rootlib << "'" << endl;
#ifdef WIN32
      string jsrootsys = rootdata + "\\js\\";
      out << "os.environ['PYTHONPATH']      = '%s' % rootlib + ':' + os.getenv('PYTHONPATH', '')" << endl;
      out << "os.environ['PATH']            = '%s:%s\\bin' % (rootbin,rootbin) + ':' + '%s' % rootlib + ':' + os.getenv('PATH', '')" << endl;
#else
      string jsrootsys = rootdata + "/js/";
      out << "os.environ['PYTHONPATH']      = '%s' % rootlib + ':' + os.getenv('PYTHONPATH', '')" << endl;
      out << "os.environ['PATH']            = '%s:%s/bin' % (rootbin,rootbin) + ':' + os.getenv('PATH', '')" << endl;
      out << "os.environ['LD_LIBRARY_PATH'] = '%s' % rootlib + ':' + os.getenv('LD_LIBRARY_PATH', '')" << endl;
#endif
      out << "c.NotebookApp.extra_static_paths = ['" << jsrootsys << "']" << endl;
      out.close();
      return true;
   }
   else {
      fprintf(stderr,
              "Error installing notebook configuration files -- cannot create IPython config file at %s\n",
              jupyconfig.c_str());
      return false;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a file that stores the current commit id in it.

static bool CreateStamp(string dest)
{
   ofstream out(dest + pathsep + COMMIT_FILE, ios::trunc);
   if (out.is_open()) {
      out << gROOT->GetGitCommit();
      out.close();
      return true;
   }
   else {
      fprintf(stderr,
              "Error installing notebook configuration files -- cannot create %s\n",
              COMMIT_FILE);
      return false;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Spawn a Jupyter notebook customised by ROOT.

int main(int argc, char **argv)
{
   string rootbin(TROOT::GetBinDir().Data());
   string rootlib(TROOT::GetLibDir().Data());
   string rootetc(TROOT::GetEtcDir().Data());
   string rootdata(TROOT::GetDataDir().Data());

   // If needed, install ROOT notebook files in the user's home directory
#ifdef WIN32
   string homedir(getenv("USERPROFILE"));
#else
   string homedir(getenv("HOME"));
#endif
   int inst = CheckNbInstallation(homedir);
   if (inst == -1) {
      // The etc directory contains the ROOT notebook files to install
      string source(rootetc + pathsep + NB_CONF_DIR);
      string dest(homedir + pathsep + ROOTNB_DIR);
      bool res = InstallNbFiles(source, dest) &&
                 CreateJupyterConfig(dest, rootbin, rootlib, rootdata) &&
                 CreateStamp(dest);
      if (!res) return 1;
   }
   else if (inst == -2) return 1;

   // Set IPython directory for the ROOT notebook flavour
   string rootnbpath = homedir + pathsep + ROOTNB_DIR;
   string jupyconfdir(JUPYTER_CONF_DIR_V + ("=" + rootnbpath));
   string jupypathdir(JUPYTER_PATH_V + ("=" + rootnbpath));
   putenv((char *)jupyconfdir.c_str());
   putenv((char *)jupypathdir.c_str());

   char **jargv = new char* [argc + 2];
   jargv[0] = (char *) JUPYTER_CMD;
   jargv[1] = (char *) NB_OPT;
   for (int n=1;n<argc;++n)
      jargv[n+1] = argv[n];
   jargv[argc+1] = nullptr;

   // Execute IPython notebook
   execvp(JUPYTER_CMD, jargv);

   // Exec failed
   fprintf(stderr,
           "Error starting ROOT notebook -- please check that Jupyter is installed\n");

   return 1;
}
