// ************************************************************************* //
// *                                                                       * //
// *                        s t r e s s P r o o f                          * //
// *                                                                       * //
// * This file contains a set of test of PROOF related functionality.      * //
// * The tests can be run as a standalone program or with the interpreter. * //
// * To run as a standalone program:                                       * //
// *                                                                       * //
// *  $ cd $ROOTSYS/test                                                   * //
// *  $ make stressProof                                                   * //
// *                                                                       * //
// * Run stressProof with '-h' to get the list of supported options        * //
// *  $ ./stressProof -h                                                   * //
// *                                                                       * //
// * To run interactively:                                                 * //
// * $ root                                                                * //
// * root[] .include $ROOTSYS/tutorials                                    * //
// * root[] .L stressProof.cxx+                                            * //
// * root[] stressProof(master, tests, wrks, verbose, logfile, dyn, \      * //
// *                    dyn, skipds, h1src, eventsrc, dryrun)              * //
// *                                                                       * //
// * The arguments have the same meaning as above except for               * //
// *     verbose [Int_t]   increasing verbosity (0 == minimal)             * //
// *     dyn     [Bool_t]  if kTRUE run in dynamic startup mode            * //
// *     skipds  [Bool_t]  if kTRUE the dataset related tests are skipped  * //
// *                                                                       * //
// * A certain number of swithces can also be controlled via environment   * //
// * variables: check './stressProof -h'                                   * //
// *                                                                       * //
// * The stressProof function returns 0 on success, 1 on failure.          * //
// *                                                                       * //
// * The successful output looks like this:                                * //
// *                                                                       * //
// *  ******************************************************************   * //
// *  *  Starting  P R O O F - S T R E S S suite                       *   * //
// *  ******************************************************************   * //
// *  *  Log file: /tmp/ProofStress_XrcwBe                                 * //
// *  ******************************************************************   * //
// *   Test  1 : Open a session ................................... OK *   * //
// *   Test  2 : Get session logs ................................. OK *   * //
// *   Test  3 : Simple random number generation .................. OK *   * //
// *   Test  4 : Dataset handling with H1 files ................... OK *   * //
// *   Test  5 : H1: chain processing ............................. OK *   * //
// *   Test  6 : H1: file collection processing ................... OK *   * //
// *   Test  7 : H1: file collection, TPacketizerAdaptive ......... OK *   * //
// *   Test  8 : H1: by-name processing ........................... OK *   * //
// *   Test  9 : H1: multi dataset processing ..................... OK *   * //
// *   Test 10 : H1: multi dataset and entry list ................. OK *   * //
// *   Test 11 : Package management with 'event' .................. OK *   * //
// *   Test 12 : Package argument passing ......................... OK *   * //
// *   Test 13 : Simple 'event' generation ........................ OK *   * //
// *   Test 14 : Input data propagation ........................... OK *   * //
// *   Test 15 : H1, Simple: async mode :.......................... OK *   * //
// *   Test 16 : Admin functionality .............................. OK *   * //
// *   Test 17 : Dynamic sub-mergers functionality ................ OK *   * //
// *   Test 18 : Event range processing ........................... OK *   * //
// *   Test 19 : Event range, TPacketizerAdaptive ................. OK *   * //
// *   Test 20 : File-resident output: merge ...................... OK *   * //
// *   Test 21 : File-resident output: merge w/ submergers ........ OK *   * //
// *   Test 22 : File-resident output: create dataset ............. OK *   * //
// *   Test 23 : File-resident output: multi trees ................ OK *   * //
// *   Test 24 : TTree friends (and TPacketizerFile) .............. OK *   * //
// *   Test 25 : TTree friends, same file ......................... OK *   * //
// *   Test 26 : Handling output via file ......................... OK *   * //
// *   Test 27 : Simple: selector by object ....................... OK *   * //
// *   Test 28 : H1 dataset: selector by object ................... OK *   * //
// *   Test 29 : Chain with TTree in subdirs ...................... OK *   * //
// *  * All registered tests have been passed  :-)                     *   * //
// *  ******************************************************************   * //
// *                                                                       * //
// * The application redirects the processing logs to a log file which is  * //
// * normally deleted at the end of a successful run; if the test fails    * //
// * the caller is asked if she/he wants to keep the log file; if the      * //
// * specifies a log file path of her/his choice, the log file is never    * //
// * deleted.                                                              * //
// *                                                                       * //
// * SKIPPED means that the test cannot be run.                            * //
// *                                                                       * //
// * New tests can be easily added by providing a function performing the  * //
// * test and a name for the test; see examples below.                     * //
// *                                                                       * //
// * It is also possible to trigger the automatic PROOF valgrind setup by  * //
// * means of the env GETPROOF_VALGRIND.                                   * //
// * E.g. to run the master in valgrind do                                 * //
// *                                                                       * //
// *     $ export GETPROOF_VALGRIND="valgrind=master"                      * //
// * or                                                                    * //
// *     $ export GETPROOF_VALGRIND="valgrind=workers"                     * //
// *                                                                       * //
// * before running stressProof. The syntax is the same as for standard    * //
// * PROOF valgrind runs. See                                              * //
// *   http://root.cern.ch/drupal/content/running-proof-query-valgrind     * //
// *                                                                       * //
// ************************************************************************* //

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <memory>
#ifdef WIN32
#include <io.h>
#endif

#include "Getline.h"
#include "TApplication.h"
#include "TChain.h"
#include "TDataMember.h"
#include "TDSet.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THashList.h"
#include "TList.h"
#include "TMacro.h"
#include "TMap.h"
#include "TMath.h"
#include "TMethodCall.h"
#include "TNamed.h"
#include "TNtuple.h"
#include "TParameter.h"
#include "TProof.h"
#include "TProofLog.h"
#include "TProofMgr.h"
#include "TProofOutputFile.h"
#include "TQueryResult.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TSelector.h"
#include "TProofProgressDialog.h"
#include "TProofProgressLog.h"
#include "TProofProgressMemoryPlot.h"
#include "TObjString.h"

#include "proof/getProof.C"

#define PT_NUMTEST 29

static const char *urldef = "proof://localhost:40000";
static TString gtutdir;
static TString gsandbox;
static Int_t gverbose = 1;
static TString gverbproof("kAll");
static TString glogfile;
static Int_t gpoints = 0;
static Bool_t guseprogress = kTRUE;
static Int_t totpoints = 53;
static RedirectHandle_t gRH;
static RedirectHandle_t gRHAdmin;
static Double_t gH1Time = 0;
static Double_t gSimpleTime = 0;
static Double_t gEventTime = 0;
static Int_t gH1Cnt = 0;
static Int_t gSimpleCnt = 0;
static Int_t gEventCnt = 0;
static TStopwatch gTimer;
static Bool_t gTimedOut = kFALSE;
static Bool_t gDynamicStartup = kFALSE;
static Bool_t gLocalCluster = kTRUE;
static Bool_t gSkipDataSetTest = kTRUE;
static Bool_t gUseParallelUnzip = kFALSE;
static Bool_t gClearCache = kFALSE;
static TString gh1src("http://root.cern.ch/files/h1");
static Bool_t gh1ok = kTRUE;
static Bool_t gh1local = kFALSE;
static char gh1sep = '/';
static const char *gh1file[] = { "dstarmb.root", "dstarp1a.root", "dstarp1b.root", "dstarp2.root" };
static TString geventsrc("http://root.cern.ch/files/data");
static Bool_t geventok = kTRUE;
static Bool_t geventlocal = kFALSE;
static Int_t geventnf = 10;
static Long64_t gEventNum = 200000;
static Long64_t gEventFst = 65000;
static Long64_t gEventSiz = 100000;
static TStopwatch gStopwatch;

// The selectors
static TString gTutDir = "$ROOTSYS/tutorials"; // default path
static TList gSelectors;
static TString gH1Sel("/tree/h1analysis.C");
static TString gEventSel("/proof/ProofEvent.C");
static TString gEventProcSel("/proof/ProofEventProc.C");
static TString gSimpleSel("/proof/ProofSimple.C");
static TString gTestsSel("/proof/ProofTests.C");
static TString gNtupleSel("/proof/ProofNtuple.C");
static TString gFriendsSel("/proof/ProofFriends.C");
static TString gAuxSel("/proof/ProofAux.C");

// Special class
static TString gProcFileElem("/proof/ProcFileElements.C");

// Special files
static TString gEmptyInclude("/proof/EmptyInclude.h");
static TString gNtpRndm("/proof/ntprndm.root");

// Special packages
static TString gPackEvent("/proof/event.par");
static TString gPack1("/proof/packtest1.par");
static TString gPack2("/proof/packtest2.par");
static TString gPack3("/proof/packtest3.par");

int stressProof(const char *url = 0,
                const char *tests = 0, Int_t nwrks = -1,
                const char *verbose = "1", const char *logfile = 0,
                Bool_t dyn = kFALSE, Bool_t skipds = kTRUE,
                const char *h1src = 0, const char *eventsrc = 0,
                Bool_t dryrun = kFALSE, Bool_t showcpu = kFALSE,
                Bool_t clearcache = kFALSE, Bool_t useprogress = kTRUE,
                const char *tutdir = 0, Bool_t cleanlog = kFALSE,
                Bool_t keeplog = kTRUE, Bool_t catlog = kFALSE);

//_____________________________batch only_____________________
#ifndef __CINT__
int main(int argc,const char *argv[])
{

   // Request for help?
   if (argc > 1 && !strcmp(argv[1],"-h")) {
      printf(" \n");
      printf(" PROOF test suite\n");
      printf(" \n");
      printf(" Usage:\n");
      printf(" \n");
      printf(" $ ./stressProof [-h] [-n <wrks>] [-d [scope:]level] [-l logfile] [-dyn] [-ds]\n");
      printf("                 [-t tests] [-h1 h1src] [-event src] [-dryrun] [-g] [-cpu]\n");
      printf("                 [-clearcache] [-noprogress] [-tut tutdir] [master]\n");
      printf(" \n");
      printf(" Optional arguments:\n");
      printf("   -h            prints this menu\n");
      printf("   master        entry point of the cluster where to run the test\n");
      printf("                 in the form '[user@]host.domain[:port]'; default 'localhost:40000'\n");
      printf("                 or the env STRESSPROOF_URL\n");
      printf("   -n wrks       number of workers to be started when running on the local host;\n");
      printf("                 default is the nuber of local cores\n");
      printf("   -d [scope:]level\n");
      printf("                 verbosity scope and level; default level is 1; scope refers to PROOF debug\n");
      printf("                 options (see TProofDebug.h) and is enabled only is level > 1 (default kAll).\n");
      printf("   -l logfile    file where to redirect the processing logs; must be writable;\n");
      printf("                 default is a temporary file deleted at the end of the test\n");
      printf("                 in case of success. Some format specifiers can be used in the file name:\n");
      printf("                 - %%p is replaced by the current process ID; \n");
      printf("                 - %%tmp is replaced by the temporary directory; \n");
      printf("                 - %%d is replaced by date and time at test startup in the form YYYYMMDD-HHMM, \n");
      printf("                   e.g. 20120904-1138. \n");
      printf("                 The log file path can be also passed via the env STRESSPROOF_LOGFILE.\n");
      printf("                 In case of failure, the log files of the nodes (master and workers) are saved into\n");
      printf("                 a file called <logfile>.nodes .\n");
      printf("   -c,-cleanlog  delete the logfile specified via '-l' in case of a successful run; by default\n");
      printf("                 the file specified by '-l' is kept in all cases (default log files are deleted\n");
      printf("                 on success); adding this switch allows to keep a user-defined log file only\n");
      printf("                 in case of error.\n");
      printf("   -k,-keeplog   keep all logfiles, including the ones from the PROOF nodes (in one single file)\n");
      printf("                 The paths are printed on the screen.\n");
      printf("   -catlog       prints all the logfiles (also the ones from PROOF nodes) on stdout; useful for\n");
      printf("                 presenting a single aggregated output for automatic tests. If specified in\n");
      printf("                 conjunction with -cleanlog it will only print the logfiles in case of errors\n");
      printf("   -dyn          run the test in dynamicStartup mode\n");
      printf("   -ds           force the dataset test if skipped by default\n");
      printf("   -t tests      run only tests in the comma-separated list and those from which they\n");
      printf("                 depend; ranges can be given with 'first-last' syntax, e.g. '3,6-9,23'\n");
      printf("                 to run tests 3, 6, 7, 8, 9 and 23. Can be also given via the env STRESSPROOF_TESTS.\n");
      printf("   -h1 h1src     specify a location for the H1 files; use h1src=\"download\" to download\n");
      printf("                 to a temporary location (or h1src=\"download=/path/for/local/h1\" to download\n");
      printf("                 to /path/for/local/h1; by default the files are read directly from the\n");
      printf("                 ROOT http server; however this may give failures if the connection is slow.\n");
      printf("                 If h1src ends with '.zip' the program assumes that the 4 files are concatenated\n");
      printf("                 in a single zip archive and are accessible with the standard archive syntax.\n");
      printf("                 Can be also passed via the env STRESSPROOF_H1SRC.\n");
      printf("   -event src\n");
      printf("                 specify a location for the 'event' files; use src=\"download\" to download\n");
      printf("                 to a temporary location (or eventsrc=\"download=/path/for/local/event\" to download\n");
      printf("                 to /path/for/local/event; by default the files are read directly from the\n");
      printf("                 ROOT http server; however this may give failures if the connection is slow\n");
      printf("                 Can be also passed via the env STRESSPROOF_EVENT.\n");
      printf("   -punzip       use parallel unzipping for data-driven processing.\n");
      printf("   -dryrun       only show which tests would be run.\n");
      printf("   -g            enable graphics; default is to run in text mode.\n");
      printf("   -cpu          show CPU times used by each successful test; used for calibration.\n");
      printf("   -clearcache   clear memory cache associated with the files processed when local.\n");
      printf("   -noprogress   do not show progress (escaped chars may confuse some wrapper applications)\n");
      printf("   -tut tutdir   specify alternative location for the ROOT tutorials; allows to run multiple\n");
      printf("                 concurrent instances of stressProof, for tests, for example.\n");
      printf("                 Can be also passed via the env STRESSPROOF_TUTORIALDIR.\n");
      printf(" \n");
      gSystem->Exit(0);
   }

   // Parse options
   const char *url = 0;
   Int_t nWrks = -1;
   const char *verbose = "1";
   Bool_t enablegraphics = kFALSE;
   Bool_t dryrun = kFALSE;
   Bool_t showcpu = kFALSE;
   Bool_t clearcache = kFALSE;
   Bool_t useprogress = kTRUE;
   Bool_t cleanlog = kFALSE;
   Bool_t keeplog = kFALSE;
   Bool_t catlog = kFALSE;
   const char *logfile = 0;
   const char *h1src = 0;
   const char *eventsrc = 0;
   const char *tests = 0;
   const char *tutdir = 0;
   Int_t i = 1;
   while (i < argc) {
      if (!strcmp(argv[i],"-h")) {
         // Ignore if not first argument
         i++;
      } else if (!strcmp(argv[i],"-n")) {
         if (i+1 == argc || argv[i+1][0] == '-') {
            printf(" -n should be followed by the number of workers: ignoring \n");
            i++;
         } else  {
            nWrks = atoi(argv[i+1]);
            i += 2;
         }
      } else if (!strcmp(argv[i],"-d")) {
         if (i+1 == argc || argv[i+1][0] == '-') {
            printf(" -d should be followed by the debug '[scope:]level': ignoring \n");
            i++;
         } else  {
            verbose = argv[i+1];
            i += 2;
         }
      } else if (!strcmp(argv[i],"-l")) {
         if (i+1 == argc || argv[i+1][0] == '-') {
            printf(" -l should be followed by a path: ignoring \n");
            i++;
         } else {
            logfile = argv[i+1];
            i += 2;
         }
      } else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "-cleanlog")) {
         cleanlog = kTRUE;
         i++;
      } else if (!strcmp(argv[i], "-k") || !strcmp(argv[i], "-keeplog")) {
         keeplog = kTRUE;
         i++;
      } else if (!strcmp(argv[i], "-catlog")) {
         catlog = kTRUE;
         i++;
      } else if (!strncmp(argv[i],"-v",2)) {
         // For backward compatibility
         if (!strncmp(argv[i],"-vv",3)) verbose = "2";
         if (!strncmp(argv[i],"-vvv",4)) verbose = "3";
         i++;
      } else if (!strncmp(argv[i],"-dyn",4)) {
         gDynamicStartup = kTRUE;
         i++;
      } else if (!strncmp(argv[i],"-ds",3)) {
         gSkipDataSetTest = kFALSE;
         i++;
      } else if (!strncmp(argv[i],"-punzip",7)) {
         gUseParallelUnzip = kTRUE;
         i++;
      } else if (!strcmp(argv[i],"-t")) {
         if (i+1 == argc || argv[i+1][0] == '-') {
            printf(" -t should be followed by a string or a number: ignoring \n");
            i++;
         } else {
            tests = argv[i+1];
            i += 2;
         }
      } else if (!strcmp(argv[i],"-h1")) {
         if (i+1 == argc || argv[i+1][0] == '-') {
            printf(" -h1 should be followed by a path: ignoring \n");
            i++;
         } else {
            h1src = argv[i+1];
            i += 2;
         }
      } else if (!strcmp(argv[i],"-event")) {
         if (i+1 == argc || argv[i+1][0] == '-') {
            printf(" -event should be followed by a path: ignoring \n");
            i++;
         } else {
            eventsrc = argv[i+1];
            i += 2;
         }
      } else if (!strncmp(argv[i],"-g",2)) {
         enablegraphics = kTRUE;
         i++;
      } else if (!strncmp(argv[i],"-dryrun",7)) {
         dryrun = kTRUE;
         i++;
      } else if (!strncmp(argv[i],"-cpu",4)) {
         showcpu = kTRUE;
         i++;
      } else if (!strncmp(argv[i],"-clearcache",11)) {
         clearcache = kTRUE;
         i++;
      } else if (!strncmp(argv[i],"-noprogress",11)) {
         useprogress = kFALSE;
         i++;
      } else if (!strcmp(argv[i],"-tut")) {
         if (i+1 == argc || argv[i+1][0] == '-') {
            printf(" -tut should be followed by a path: ignoring \n");
            i++;
         } else {
            tutdir = argv[i+1];
            i += 2;
         }
      } else {
         url = argv[i];
         i++;
      }
   }

   // Enable graphics if required
   if (enablegraphics) {
      new TApplication("stressProof", 0, 0);
   } else {
      gROOT->SetBatch(kTRUE);
   }

   int rc = stressProof(url, tests, nWrks, verbose, logfile, gDynamicStartup, gSkipDataSetTest,
                        h1src, eventsrc, dryrun, showcpu, clearcache, useprogress, tutdir, cleanlog,
                        keeplog, catlog);

   gSystem->Exit(rc);
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Print one '.' and count it

Int_t PutPoint()
{
   printf(".");
   return ++gpoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Print some progress information

void PrintStressProgress(Long64_t total, Long64_t processed, Float_t, Long64_t)
{
   gSystem->RedirectOutput(0, 0, &gRH);

   char pc[2] = { '.', ':'};
   static int lstpc = 1;

   int ns = 0;
   if (processed < total) {
      ns += 6;
      fprintf(stderr, "%c %2.0f %%", pc[(lstpc++)%2],
                                    (total ? ((100.0*processed)/total) : 100.0));
   }
   while (ns--) fprintf(stderr, "\b");

   gSystem->RedirectOutput(glogfile, "a", &gRH);
}
////////////////////////////////////////////////////////////////////////////////
/// Dummy PrintProgress

void PrintEmptyProgress(Long64_t, Long64_t, Float_t, Long64_t)
{
   return;
}

// Guard class
class SwitchProgressGuard {
public:
   SwitchProgressGuard(Bool_t force = kFALSE) {
      if (guseprogress || force) {
         gProof->SetPrintProgress(&PrintStressProgress);
      } else {
         gProof->SetPrintProgress(&PrintEmptyProgress);
      }
   }
   ~SwitchProgressGuard() { gProof->SetPrintProgress(0); }
};

////////////////////////////////////////////////////////////////////////////////
/// Remove all non source files associated with seletor at path 'selpath'

void CleanupSelector(const char *selpath)
{
   if (!selpath) return;

   TString dirpath = gSystem->GetDirName(selpath);
   if (gSystem->AccessPathName(dirpath)) return;
   TString selname(gSystem->BaseName(selpath));
   selname.ReplaceAll(".C", "_C");
   void *dirp = gSystem->OpenDirectory(dirpath);
   if (!dirp) return;
   TString fn;
   const char *e = 0;
   while ((e = gSystem->GetDirEntry(dirp))) {
      if (!strncmp(e, selname.Data(), selname.Length())) {
         // Cleanup this entry
         fn.Form("%s/%s", dirpath.Data(), e);
         gSystem->Unlink(fn);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the parallel unzip option

void AssertParallelUnzip()
{
   if (gUseParallelUnzip) {
      gProof->SetParameter("PROOF_UseParallelUnzip", (Int_t)1);
   } else {
      gProof->SetParameter("PROOF_UseParallelUnzip", (Int_t)0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Release the memory cache associated with file 'fn'.

void ReleaseCache(const char *fn)
{
#if defined(R__LINUX)
   TString filename(fn);
   Int_t fd;
   fd = open(filename.Data(), O_RDONLY);
   if (fd > -1) {
      fdatasync(fd);
      posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
      close(fd);
   } else {
      fprintf(stderr, "cannot open file '%s' for cache clean up; errno=%d \n",
                      filename.Data(), errno);
   }
#else
   fprintf(stderr, "ReleaseCache: dummy function: file '%s' untouched ...\n", fn);
#endif
   // Done
   return;
}

//
// Auxilliary classes for testing
//
class RunTimes {
public:
   Double_t fCpu;
   Double_t fReal;
   RunTimes(Double_t c = -1., Double_t r = -1.) : fCpu(c), fReal(r) { }

   void Set(Double_t c = -1., Double_t r = -1.) { if (c > -1.) fCpu = c; if (r > -1.) fReal = r; }
   void Print(const char *tag = "") { printf("%s real: %f s, cpu: %f s\n", tag, fReal, fCpu); }
};
RunTimes operator-(const RunTimes &rt1, const RunTimes &rt2) {
   RunTimes rt(rt1.fCpu - rt2.fCpu, rt1.fReal - rt2.fReal);
   return rt;
}
static RunTimes gProofTimesZero(0.,0.);

typedef Int_t (*ProofTestFun_t)(void *, RunTimes &);
class ProofTest : public TNamed {
private:
   Int_t           fSeq;  // Sequential number for the test
   ProofTestFun_t  fFun;  // Function to be executed for the test
   void           *fArgs; // Arguments to be passed to the function
   TString         fDeps;  // Test dependencies, e.g. "1,3"
   TString         fSels;  // Selectors used, e.g. "h1analysis,ProofSimple"
   Int_t           fDepFrom; // Index for looping over deps
   Int_t           fSelFrom; // Index for looping over selectors
   Bool_t          fEnabled; // kTRUE if this test is enabled
   Double_t        fCpuTime; // CPU time used by the test
   Double_t        fRealTime; // Real time used by the test
   Double_t        fRefReal; // Ref Real time used for PROOF marks
   Double_t        fProofMarks; // PROOF marks
   Bool_t          fUseForMarks; // Use in the calculation of the average PROOF marks

   static Double_t gRefReal[PT_NUMTEST]; // Reference Cpu times

public:
   ProofTest(const char *n, Int_t seq, ProofTestFun_t f, void *a = 0,
             const char *d = "", const char *sel = "", Bool_t useForMarks = kFALSE)
           : TNamed(n,""), fSeq(seq), fFun(f), fArgs(a),
             fDeps(d), fSels(sel), fDepFrom(0), fSelFrom(0), fEnabled(kTRUE),
             fCpuTime(-1.), fRealTime(-1.), fProofMarks(-1.), fUseForMarks(useForMarks)
             { fRefReal = gRefReal[fSeq-1]; }
   virtual ~ProofTest() { }

   void   Disable() { fEnabled = kFALSE; }
   void   Enable() { fEnabled = kTRUE; }
   Bool_t IsEnabled() const { return fEnabled; }

   Int_t  NextDep(Bool_t reset = kFALSE);
   Int_t  NextSel(TString &sel, Bool_t reset = kFALSE);
   Int_t  Num() const { return fSeq; }

   Int_t  Run(Bool_t dryrun = kFALSE, Bool_t showcpu = kFALSE);

   Double_t ProofMarks() const { return fProofMarks; }
   Bool_t UseForMarks() const { return fUseForMarks; }
};

// Reference time measured on a HP DL580 24 core (4 x Intel(R) Xeon(R) CPU X7460
// @ 2.132 GHz, 48GB RAM, 1 Gb/s NIC) with 4 workers.
Double_t ProofTest::gRefReal[PT_NUMTEST] = {
   3.047808,   // #1:  Open a session
   0.021979,   // #2:  Get session logs
   4.779971,   // #3:  Simple random number generation
   0.276155,   // #4:  Dataset handling with H1 files
   5.355514,   // #5:  H1: chain processing
   2.414207,   // #6:  H1: file collection processing
   3.381990,   // #7:  H1: file collection, TPacketizerAdaptive
   3.227942,   // #8:  H1: by-name processing
   3.944204,   // #9:  H1: multi dataset processing
   9.146988,   // #10: H1: multi dataset and entry list
   2.703881,   // #11: Package management with 'event'
   3.814625,   // #12: Package argument passing
   9.028315,   // #13: Simple 'event' generation
   0.123514,   // #14: Input data propagation
   0.129757,   // #15: H1, Simple: async mode
   0.349625,   // #16: Admin functionality
   0.989456,   // #17: Dynamic sub-mergers functionality
   11.23798,   // #18: Event range processing
   6.087582,   // #19: Event range, TPacketizerAdaptive
   2.489555,   // #20: File-resident output: merge
   0.180897,   // #21: File-resident output: merge w/ submergers
   1.417233,   // #22: File-resident output: create dataset
   0.000000,   // #23: File-resident output: multi trees
   7.452465,   // #24: TTree friends (and TPacketizerFile)
   0.259239,   // #25: TTree friends, same file
   6.868858,   // #26: Simple generation: merge-via-file
   6.362017,   // #27: Simple random number generation by TSelector object
   5.519631,   // #28: H1: by-object processing
   7.452465    // #29: Chain with TTree in subdirs
};

//
// Timer to stop asynchronous actions
//
class TTimeOutTimer : public TTimer {
public:
   TTimeOutTimer(Long_t ms);
   Bool_t  Notify();
};

TTimeOutTimer::TTimeOutTimer(Long_t ms)
              : TTimer(ms, kTRUE)
{
   //constructor
   gSystem->AddTimer(this);
}

Bool_t TTimeOutTimer::Notify()
{
   //notifier
   gTimedOut = kTRUE;
   Remove();       // one shot only
   return kTRUE;
}
//------------------------------------------------------------------------------
//_____________________________________________________________________________
Int_t ProofTest::NextDep(Bool_t reset)
{
   // Return index of next dependency or -1 if none (or no more)
   // If reset is kTRUE, reset the internal counter before acting.

   if (reset) fDepFrom = 0;

   TString tkn;
   if (fDeps.Tokenize(tkn, fDepFrom, ",")) {
      if (tkn.IsDigit()) return tkn.Atoi();
   }
   // Not found
   return -1;
}
////////////////////////////////////////////////////////////////////////////////
/// Return index of next dependency or -1 if none (or no more)
/// If reset is kTRUE, reset the internal counter before acting.

Int_t ProofTest::NextSel(TString &sel, Bool_t reset)
{
   if (reset) fSelFrom = 0;
   if (fSels.Tokenize(sel, fSelFrom, ",")) {
      if (!sel.IsNull()) return 0;
   }
   // Not found
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Generic stress steering function; returns 0 on success, -1 on error

Int_t ProofTest::Run(Bool_t dryrun, Bool_t showcpu)
{
   gpoints = 0;
   printf(" Test %2d : %s ", fSeq, GetName());
   PutPoint();
   Int_t rc = 0;
   if (!dryrun) {
      gSystem->RedirectOutput(glogfile, "a", &gRH);
      RunTimes tt;
      gStopwatch.Start();
      rc = (*fFun)(fArgs, tt);
      gStopwatch.Stop();
      fCpuTime = tt.fCpu;
      fRealTime = tt.fReal;
      // Proof marks
      if (fRefReal > 0)
         fProofMarks = (fRealTime > 0) ? 1000 * fRefReal / fRealTime : -1;
      gSystem->RedirectOutput(0, 0, &gRH);
      if (rc == 0) {
         Int_t np = totpoints - strlen(GetName()) - strlen(" OK *");
         while (np--) { printf("."); }
         if (showcpu) {
            printf(" OK * (rt: %f s, cpu: %f s, marks: %.2f)\n",
                   fRealTime, fCpuTime, fProofMarks);
         } else {
            printf(" OK *\n");
         }
      } else if (rc == 1) {
         Int_t np = totpoints - strlen(GetName()) - strlen(" SKIPPED *");
         while (np--) { printf("."); }
         printf(" SKIPPED *\n");
      } else {
         Int_t np = totpoints - strlen(GetName()) - strlen(" FAILED *");
         while (np--) { printf("."); }
         printf(" FAILED *\n");
         gSystem->ShowOutput(&gRH);
      }
   } else {
      if (fEnabled) {
         Int_t np = totpoints - strlen(GetName()) - strlen(" ENABLED *");
         while (np--) { printf("."); }
         printf(" ENABLED *\n");
      } else {
         Int_t np = totpoints - strlen(GetName()) - strlen(" DISABLED *");
         while (np--) { printf("."); }
         printf(" DISABLED *\n");
      }
   }
   // Done
   return rc;
}

// Test functions
Int_t PT_Open(void *, RunTimes &);
Int_t PT_GetLogs(void *, RunTimes &);
Int_t PT_Simple(void *, RunTimes &);
Int_t PT_H1Http(void *, RunTimes &);
Int_t PT_H1FileCollection(void *, RunTimes &);
Int_t PT_H1DataSet(void *, RunTimes &);
Int_t PT_H1MultiDataSet(void *, RunTimes &);
Int_t PT_H1MultiDSetEntryList(void *, RunTimes &);
Int_t PT_DataSets(void *, RunTimes &);
Int_t PT_Packages(void *, RunTimes &);
Int_t PT_Event(void *, RunTimes &);
Int_t PT_InputData(void *, RunTimes &);
Int_t PT_H1SimpleAsync(void *arg, RunTimes &);
Int_t PT_AdminFunc(void *arg, RunTimes &);
Int_t PT_PackageArguments(void *, RunTimes &);
Int_t PT_EventRange(void *, RunTimes &);
Int_t PT_POFNtuple(void *, RunTimes &);
Int_t PT_POFDataset(void *, RunTimes &);
Int_t PT_Friends(void *, RunTimes &);
Int_t PT_TreeSubDirs(void *, RunTimes &);
Int_t PT_SimpleByObj(void *, RunTimes &);
Int_t PT_H1ChainByObj(void *, RunTimes &);
Int_t PT_AssertTutorialDir(const char *tutdir);
Int_t PT_MultiTrees(void *, RunTimes &);
Int_t PT_OutputHandlingViaFile(void *, RunTimes &);

// Auxilliary functions
void PT_GetLastTimes(RunTimes &tt)
{
   // Get in tt the Cpu and Real times used during the run
   // The runtimes
   gStopwatch.Stop();
   tt.fCpu = gStopwatch.CpuTime();
   tt.fReal = gStopwatch.RealTime();
}
void PT_GetLastProofTimes(RunTimes &tt)
{
   // Get in tt the Cpu and Real times used by PROOF since last call
   if (gProof->IsLite()) {
      PT_GetLastTimes(tt);
      return;
   }
   // Update the statistics
   gProof->GetStatistics(kFALSE);
      // The runtimes
   RunTimes proofTimesCurrent(gProof->GetCpuTime(), gProof->GetRealTime());
   tt = proofTimesCurrent - gProofTimesZero;
   gProofTimesZero = proofTimesCurrent;
}

// Arguments structures
typedef struct ptopenargs {            // Open
   const char *url;
   Int_t       nwrks;
} PT_Open_Args_t;

// Packetizer parameters
typedef struct ptpacketizer {
   const char *fName;
   Int_t fType;
} PT_Packetizer_t;

// Options
typedef struct ptoption {
   Int_t fOne;
   Int_t fTwo;
} PT_Option_t;

static PT_Packetizer_t gStd_Old = { "TPacketizerAdaptive", 1 };

////////////////////////////////////////////////////////////////////////////////

int stressProof(const char *url, const char *tests, Int_t nwrks,
                const char *verbose, const char *logfile, Bool_t dyn, Bool_t skipds,
                const char *h1src, const char *eventsrc,
                Bool_t dryrun, Bool_t showcpu, Bool_t clearcache, Bool_t useprogress,
                const char *tutdir, Bool_t cleanlog, Bool_t keeplog, Bool_t catlog)
{
   printf("******************************************************************\n");
   printf("*  Starting  P R O O F - S T R E S S  suite                      *\n");
   printf("******************************************************************\n");

   // Use defaults or environment settings where required
   if (!url) {
      url = getenv("STRESSPROOF_URL");
      if (!url) url = urldef;
   }
   // Set dynamic mode
   gDynamicStartup = (!strcmp(url,"lite://")) ? kFALSE : dyn;

   // Set verbosity
   TString vv(verbose);
   Ssiz_t icol = kNPOS;
   if ((icol = vv.Index(":")) != kNPOS) {
      TString sv = vv(icol+1, vv.Length() - icol -1);
      gverbose = sv.Atoi();
      vv.Remove(icol);
      gverbproof = vv;
   } else {
      gverbose = vv.Atoi();
   }

   // No progress bar if not tty or explicitly not requested (i.e. for ctest)
   guseprogress = useprogress;
   if (isatty(0) == 0 || isatty(1) == 0) guseprogress = kFALSE;
   if (!guseprogress) {
      if (!useprogress) {
         printf("*  Progress not shown (explicit request)                         *\n");
      } else {
         printf("*  Progress not shown (not tty)                                  *\n");
      }
      printf("******************************************************************\n");
   }

   // Notify/warn about the dynamic startup option, if any
   TUrl uu(url), udef(urldef);
   Bool_t extcluster = ((strcmp(uu.GetHost(), udef.GetHost()) ||
                        (uu.GetPort() != udef.GetPort())) && strcmp(url,"lite://"))? kTRUE : kFALSE;
   if (gDynamicStartup && gverbose > 0) {
      // Check url
      if (extcluster) {
         printf("*   WARNING: request to run a test with per-job scheduling on    *\n");
         printf("*            an external cluster: %s .\n", url);
         printf("*            Make sure the dynamic option is set.                *\n");
         printf("*                                                                *\n");
         gDynamicStartup = kFALSE;
      } else {
         printf("*  Runnning in dynamic mode (per-job scheduling)                 *\n");
      }
      printf("*  Tests #15, #22, #23, #26 skipped in dynamic mode             **\n");
      printf("******************************************************************\n");
   }

   if (dryrun) {
      printf("*  Dry-run: only showing what would be done                      *\n");
      printf("******************************************************************\n");
   }

   // Cluster locality
   gLocalCluster = (!extcluster || !strcmp(url, "lite://")) ? kFALSE : kTRUE;

   // Dataset option
   if (!skipds) {
      gSkipDataSetTest = kFALSE;
   } else {
     gSkipDataSetTest = gLocalCluster;
   }

   // Clear cache
   gClearCache = clearcache;

   // Log file path
   Bool_t usedeflog = kTRUE;
   FILE *flog = 0;
   if (!logfile) logfile = getenv("STRESSPROOF_LOGFILE");
   if (logfile && strlen(logfile) > 0 && !dryrun) {
      usedeflog = kFALSE;
      glogfile = logfile;
      if (glogfile.Contains("%tmp"))
         glogfile.ReplaceAll("%tmp", gSystem->TempDirectory());
      if (glogfile.Contains("%p")) {
         TString apid = TString::Format("%d", gSystem->GetPid());
         glogfile.ReplaceAll("%p", apid);
      }
      if (glogfile.Contains("%d")) {
         TString d(TDatime().AsSQLString());
         d.Remove(d.Last(':'));
         d.ReplaceAll("-", "");
         d.ReplaceAll(":", "");
         d.ReplaceAll(" ", "-");
         glogfile.ReplaceAll("%d", d);
      }
      if (!gSystem->AccessPathName(glogfile, kFileExists)) {
         if (gSystem->AccessPathName(glogfile, kWritePermission)) {
            printf(" >>> Cannot write to log file %s - ignore file request\n", logfile);
            usedeflog = kTRUE;
         }
      } else {
         // Create the file
         if (!(flog = fopen(logfile, "w"))) {
            printf(" >>> Cannot create log file %s - ignore file request\n", logfile);
            usedeflog = kTRUE;
         }
      }
   }
   if (usedeflog && !dryrun) {
      glogfile = "ProofStress_";
      if (!(flog = gSystem->TempFileName(glogfile, gSystem->TempDirectory()))) {
         printf(" >>> Cannot create a temporary log file on %s - exit\n", gSystem->TempDirectory());
         return 1;
      }
      fclose(flog);
   }
   if (gverbose > 0) {
      printf("*  Log file: %s\n", glogfile.Data());
      if (cleanlog)
         printf("*  (NB: file will be removed if test is successful)              *\n");
      printf("******************************************************************\n");
   }

   if (gSkipDataSetTest && gverbose > 0) {
      printf("*  Test for dataset handling (#4, #8-10) skipped                **\n");
      printf("******************************************************************\n");
   }
   if (gUseParallelUnzip && gverbose > 0) {
      printf("*  Using parallel unzip where relevant                          **\n");
      printf("******************************************************************\n");
   }
   if (!strcmp(url,"lite://") && gverbose > 0) {
      printf("*  PROOF-Lite session (tests #15 and #16 skipped)               **\n");
      printf("******************************************************************\n");
   }
   if (!h1src) h1src = getenv("STRESSPROOF_H1SRC");
   if (h1src && strlen(h1src)) {
      if (!strcmp(h1src, "download") && extcluster) {
         if (gverbose > 0) {
            printf("*  External clusters: ignoring download request of H1 files\n");
            printf("******************************************************************\n");
         }
      } else if (!gh1src.BeginsWith(h1src)) {
         if (gverbose > 0) {
            printf("*  Taking 'h1' files from: %s\n", h1src);
            printf("******************************************************************\n");
         }
         gh1src = h1src;
         gh1ok = kFALSE;
      }
   }
   if (!eventsrc) eventsrc = getenv("STRESSPROOF_EVENT");
   if (eventsrc && strlen(eventsrc)) {
      if (!strcmp(eventsrc, "download") && extcluster) {
         if (gverbose > 0) {
            printf("*  External clusters: ignoring download request of 'event' files\n");
            printf("******************************************************************\n");
         }
      } else if (!geventsrc.BeginsWith(eventsrc)) {
         if (gverbose > 0) {
            printf("*  Taking 'event' files from: %s\n", eventsrc);
            printf("******************************************************************\n");
         }
         geventsrc = eventsrc;
         geventok = kFALSE;
      }
   }
   if (!tutdir) tutdir = getenv("STRESSPROOF_TUTORIALDIR");
   if (tutdir && strlen(tutdir)) {
      if (!(gTutDir == tutdir)) {
         if (gverbose > 0) {
            printf("*  Taking tutorial files from: %s\n", tutdir);
            printf("******************************************************************\n");
         }
         gTutDir = tutdir;
      }
   }
   if (clearcache && gverbose > 0) {
      printf("*  Clearing cache associated to files, if possible ...          **\n");
      printf("******************************************************************\n");
   }
   if (keeplog && gverbose > 0) {
      printf("*  Keeping logfiles (paths specified at the end)                **\n");
      printf("******************************************************************\n");
   }
   //
   // Reset dataset settings
   gEnv->SetValue("Proof.DataSetManager","");

   //
   // Register tests
   //
   TList *testList = new TList;
   // Simple open
   PT_Open_Args_t PToa = { url, nwrks };
   testList->Add(new ProofTest("Open a session", 1, &PT_Open, (void *)&PToa));
   // Get logs
   testList->Add(new ProofTest("Get session logs", 2, &PT_GetLogs, (void *)&PToa, "1"));
   // Simple histogram generation
   testList->Add(new ProofTest("Simple random number generation", 3, &PT_Simple, 0, "1", "ProofSimple", kTRUE));
   // Test of data set handling with the H1 http files
   testList->Add(new ProofTest("Dataset handling with H1 files", 4, &PT_DataSets, 0, "1"));
   // H1 analysis over HTTP (chain)
   testList->Add(new ProofTest("H1: chain processing", 5, &PT_H1Http, 0, "1", "h1analysis", kTRUE));
   // H1 analysis over HTTP (file collection)
   testList->Add(new ProofTest("H1: file collection processing", 6, &PT_H1FileCollection, 0, "1", "h1analysis", kTRUE));
   // H1 analysis over HTTP: adaptive packetizer
   testList->Add(new ProofTest("H1: file collection, TPacketizerAdaptive", 7,
                               &PT_H1FileCollection, (void *)&gStd_Old, "1", "h1analysis", kTRUE));
   // H1 analysis over HTTP by dataset name
   testList->Add(new ProofTest("H1: by-name processing", 8, &PT_H1DataSet, 0, "1,4", "h1analysis", kTRUE));
   // H1 analysis over HTTP by dataset name split in two
   testList->Add(new ProofTest("H1: multi dataset processing", 9, &PT_H1MultiDataSet, 0, "1,4", "h1analysis", kTRUE));
   // H1 analysis over HTTP by dataset name
   testList->Add(new ProofTest("H1: multi dataset and entry list", 10, &PT_H1MultiDSetEntryList, 0, "1,4", "h1analysis", kTRUE));
   // Test package management with 'event'
   testList->Add(new ProofTest("Package management with 'event'", 11, &PT_Packages, 0, "1"));
   // Test package argument passing
   testList->Add(new ProofTest("Package argument passing", 12, &PT_PackageArguments, 0, "1", "ProofTests"));
   // Simple event analysis
   testList->Add(new ProofTest("Simple 'event' generation", 13, &PT_Event, 0, "1,11", "ProofEvent", kTRUE));
   // Test input data propagation (it only works in the static startup mode)
   testList->Add(new ProofTest("Input data propagation", 14, &PT_InputData, 0, "1", "ProofTests"));
   // Test asynchronous running
   testList->Add(new ProofTest("H1, Simple: async mode", 15, &PT_H1SimpleAsync, 0, "1,3,5", "h1analysis,ProofSimple", kTRUE));
   // Test admin functionality
   testList->Add(new ProofTest("Admin functionality", 16, &PT_AdminFunc, 0, "1"));
   // Test merging via submergers
   PT_Option_t pfoptm = {1, 0};
   testList->Add(new ProofTest("Dynamic sub-mergers functionality", 17,
                                &PT_Simple, (void *)&pfoptm, "1", "ProofSimple", kTRUE));
   // Test range chain and dataset processing EventProc
   testList->Add(new ProofTest("Event range processing", 18,
                               &PT_EventRange, 0, "1,11", "ProofEventProc,ProcFileElements", kTRUE));
   // Test range chain and dataset processing EventProc with TPacketizerAdaptive
   testList->Add(new ProofTest("Event range, TPacketizerAdaptive", 19,
                               &PT_EventRange, (void *)&gStd_Old, "1,11", "ProofEventProc,ProcFileElements", kTRUE));
   // Test TProofOutputFile technology for ntuple creation
   testList->Add(new ProofTest("File-resident output: merge", 20, &PT_POFNtuple, 0, "1", "ProofNtuple", kTRUE));
   // Test TProofOutputFile technology for ntuple creation using submergers
   testList->Add(new ProofTest("File-resident output: merge w/ submergers", 21,
                               &PT_POFNtuple, (void *)&pfoptm, "1", "ProofNtuple", kTRUE));
   // Test TProofOutputFile technology for dataset creation (tests TProofDraw too)
   testList->Add(new ProofTest("File-resident output: create dataset", 22, &PT_POFDataset, 0, "1", "ProofNtuple", kTRUE));
   // Test selecting different TTrees in same files
   testList->Add(new ProofTest("File-resident output: multi trees", 23, &PT_MultiTrees, 0, "1,22", "ProofNtuple", kTRUE));
   // Test TPacketizerFile and TTree friends in separate files
   testList->Add(new ProofTest("TTree friends (and TPacketizerFile)", 24, &PT_Friends, 0, "1", "ProofFriends,ProofAux", kTRUE));
   // Test TPacketizerFile and TTree friends in same file
   Bool_t sameFile = kTRUE;
   testList->Add(new ProofTest("TTree friends, same file", 25,
                               &PT_Friends, (void *)&sameFile, "1", "ProofFriends,ProofAux", kTRUE));

   // Test handling output via file
   testList->Add(new ProofTest("Handling output via file", 26,
                                &PT_OutputHandlingViaFile, 0, "1", "ProofSimple", kTRUE));
   // Simple histogram generation by TSelector object
   testList->Add(new ProofTest("Simple: selector by object", 27, &PT_SimpleByObj, 0, "1", "ProofSimple", kTRUE));
   // H1 analysis over HTTP by TSeletor object
   testList->Add(new ProofTest("H1 chain: selector by object", 28, &PT_H1ChainByObj, 0, "1", "h1analysis", kTRUE));
   // Test TPacketizerFile and TTree friends in separate files
   testList->Add(new ProofTest("Chain with TTree in subdirs", 29, &PT_TreeSubDirs, 0, "1", "ProofFriends,ProofAux", kTRUE));
   // The selectors
   if (PT_AssertTutorialDir(gTutDir) != 0) {
      printf("*  Some of the tutorial files are missing! Stop\n");
      return 1;
   }
   gSelectors.Add(new TNamed("h1analysis", gH1Sel.Data()));
   gSelectors.Add(new TNamed("ProofEvent", gEventSel.Data()));
   gSelectors.Add(new TNamed("ProofEventProc", gEventProcSel.Data()));
   gSelectors.Add(new TNamed("ProofSimple", gSimpleSel.Data()));
   gSelectors.Add(new TNamed("ProofTests", gTestsSel.Data()));
   gSelectors.Add(new TNamed("ProcFileElements", gProcFileElem.Data()));
   gSelectors.Add(new TNamed("ProofNtuple", gNtupleSel.Data()));
   gSelectors.Add(new TNamed("ProofFriends", gFriendsSel.Data()));
   gSelectors.Add(new TNamed("ProofAux", gAuxSel.Data()));
   if (gverbose > 0) {
      if (!dryrun) {
         printf("*  Cleaning-up all non-source files associated to:\n");
      } else {
         printf("*  Would clean-up all non-source files associated to:\n");
      }
   }

   // Check what to run
   ProofTest *t = 0, *treq = 0;
   TIter nxt(testList);
   Bool_t all = kTRUE;
   if (!tests) tests = getenv("STRESSPROOF_TESTS");
   if (tests && strlen(tests)) {
      TString tts(tests), tsg, ts, ten;
      Ssiz_t from = 0;
      while (tts.Tokenize(tsg, from, ",")) {
         if (tsg.CountChar('-') > 1) {
            printf("*                                                               **\r");
            printf("*  Wrong syntax for test range specification: %s\n", tsg.Data());
            continue;
         }
         Ssiz_t fromg = 0;
         Int_t test = -1, last = -1;
         while (tsg.Tokenize(ts, fromg, "-")) {
            if (!ts.IsDigit()) {
               printf("*                                                               **\r");
               printf("*  Test string is not a digit: %s\n", ts.Data());
               continue;
            }
            if (test < 0) {
               test = ts.Atoi();
            } else {
               last = ts.Atoi();
            }
         }
         if (test <= 0) {
            printf("*                                                               **\r");
            printf("*  Non-positive test number: %d\n", test);
            continue;
         }
         const int tmx = PT_NUMTEST;
         if (test > tmx) {
            printf("*                                                               **\r");
            printf("*  Unknown test number: %d\n", test);
            continue;
         }
         if (last > tmx) {
            printf("*                                                               **\r");
            printf("*  Upper test range too large: %d - rescaling\n", last);
            last = tmx;
         } else if (last <= 0) {
            last = test;
         }
         // For final notification
         if (ten != "") ten += ",";
         ten += tsg;
         // Ok, now we can enable
         if (all) {
            all = kFALSE;
            // Disable first all the tests
            while ((t = (ProofTest *)nxt())) { t->Disable(); }
         }
         // Process one by one all enabled tests
         TString cleaned;
         for (Int_t xt = test; xt <= last; xt++) {
            // Locate the ProofTest instance
            nxt.Reset();
            while ((t = (ProofTest *)nxt())) {
               if (t->Num() == xt) treq = t;
            }
            if (!treq) {
               printf("*                                                               **\r");
               printf("*  Test %2d not found among the registered tests - exiting\n", xt);
               printf("******************************************************************\n");
               return 1;
            }
            // Enable the required tests
            Int_t tn = -1;
            while ((tn = treq->NextDep()) > 0) {
               nxt.Reset();
               while ((t = (ProofTest *)nxt())) {
                  if (t->Num() == tn) {
                     t->Enable();
                     break;
                  }
               }
            }
            // Reset associated selectors
            TString sel, tit;
            while ((treq->NextSel(sel)) == 0) {
               TNamed *nm = (TNamed *) gSelectors.FindObject(sel.Data());
               if (nm && cleaned.Index(nm->GetTitle()) == kNPOS) {
                  if (!dryrun) CleanupSelector(nm->GetTitle());
                  if (gverbose > 0) {
                     tit = nm->GetName(); tit.Resize(18);
                     printf("*     %s in %s\n", tit.Data(), gSystem->GetDirName(nm->GetTitle()).Data());
                  }
                  cleaned += TString::Format(":%s:", nm->GetTitle());
               }
            }
            // Enable the required test
            treq->Enable();
         }
      }
      if (!all && !dryrun) {
         // Notify the enabled tests
         printf("*                                                               **\r");
         printf("*  Running only test(s) %s (and related)\n", ten.Data());
         printf("******************************************************************\n");
      }
   }
   if (all) {
      // Clean all the selectors
      TString tit;
      TIter nxs(&gSelectors);
      TNamed *nm = 0;
      while ((nm = (TNamed *)nxs())) {
         if (!dryrun) CleanupSelector(nm->GetTitle());
         if (gverbose > 0) {
            tit = nm->GetName(); tit.Resize(18);
            printf("*     %s in %s\n", tit.Data(), gSystem->GetDirName(nm->GetTitle()).Data());
         }
      }
   }
   if (gverbose > 0)
      printf("******************************************************************\n");

   // If a dry-run show what we would do and exit
   if (dryrun) {
      nxt.Reset();
      while ((t = (ProofTest *)nxt())) { t->Run(kTRUE); }
      printf("******************************************************************\n");
      return 0;
   }

   // Add the ACLiC option to the selector strings
   gH1Sel += "+";
   gEventSel += "+";
   gEventProcSel += "+";
   gSimpleSel += "+";
   gTestsSel += "+";
   gProcFileElem += "+";
   gNtupleSel += "+";
   gFriendsSel += "+";
   gAuxSel += "+";

   //
   // Run the tests
   //
   Bool_t failed = kFALSE;
   nxt.Reset();
   while ((t = (ProofTest *)nxt()))
      if (t->IsEnabled()) {
         if (t->Run(kFALSE, showcpu) < 0) {
            failed = kTRUE;
            break;
         }
      }

   // Done
   Bool_t kept = ((usedeflog || cleanlog) && !keeplog) ? kFALSE : kTRUE;
   if (failed) {
      kept = kTRUE;
      if (usedeflog && !gROOT->IsBatch() && !keeplog) {
         const char *answer = Getline(" Some tests failed: would you like to keep the log file (N,Y)? [Y] ");
         if (answer && (answer[0] == 'N' || answer[0] == 'n')) {
            // Remove log file
            gSystem->Unlink(glogfile);
            kept = kFALSE;
         }
      }
   } else {
      printf("* All registered tests have been passed  :-)                     *\n");
   }

   if (kept) {
      TString logfiles(glogfile);
      // Save also the logs from the workers
      TProofMgr *mgr = gProof ? gProof->GetManager() : 0;
      if (mgr) {
         gSystem->RedirectOutput(glogfile, "a", &gRH);
         TProofLog *pl = mgr->GetSessionLogs();
         if (pl) {
            logfiles += ".nodes";
            pl->Retrieve("*",  TProofLog::kAll, logfiles);
            gSystem->RedirectOutput(0, 0, &gRH);
            SafeDelete(pl);
         } else {
            gSystem->RedirectOutput(0, 0, &gRH);
            printf("+++ Warning: could not get the session logs\n");
         }
      } else {
         printf("+++ Warning: could not attach to manager to get the session logs\n");
      }
      printf("******************************************************************\n");
      printf(" Main log file kept at %s (Proof logs in %s)\n", glogfile.Data(), logfiles.Data());
      if (catlog) {

         // Display all logfiles directly on this terminal. Useful for getting
         // test results without accessing the test machine (i.e. with CDash)
         const size_t readbuf_size = 500;
         char readbuf[readbuf_size];

         printf("******************************************************************\n");
         printf("Content of the main log file: %s\n", glogfile.Data());
         printf("******************************************************************\n");
         std::ifstream glogfile_is( glogfile.Data() );
         if (!glogfile_is) {
            printf("Cannot open %s", glogfile.Data());
         }
         else {
            while ( glogfile_is.good() ) {
               glogfile_is.getline(readbuf, readbuf_size);
               std::cout << readbuf << std::endl;
            }
            glogfile_is.close();
         }

         printf("******************************************************************\n");
         printf("Content of the PROOF servers log files: %s\n", logfiles.Data());
         printf("******************************************************************\n");
         std::ifstream logfiles_is( logfiles.Data() );
         if (!logfiles_is) {
            printf("Cannot open %s", logfiles.Data());
         }
         else {
            while ( logfiles_is.good() ) {
               logfiles_is.getline(readbuf, readbuf_size);
               std::cout << readbuf << std::endl;
            }
            logfiles_is.close();
         }
      }
   } else {
      // Remove log file if not passed by the user
      gSystem->Unlink(glogfile);
   }

   printf("******************************************************************\n");
   if (gProof) {
      // Get average of single PROOF marks
      int navg = 0;
      double avgmarks = -1.;
      nxt.Reset();
      if (gverbose > 1) printf(" PROOFMARKS for single tests (-1 if unmeasured): \n");
      while ((t = (ProofTest *)nxt()))
         if (t->IsEnabled()) {
            if (gverbose > 1) printf(" %d:\t %.2f \n", t->Num(), t->ProofMarks());
            if (t->UseForMarks() && t->ProofMarks() > 0) {
               navg++;
               avgmarks += t->ProofMarks();
            }
         }
      if (navg > 0) avgmarks /= navg;

      gProof->GetStatistics((gverbose > 0));
      // Reference time measured on a HP DL580 24 core (4 x Intel(R) Xeon(R) CPU X7460
      // @ 2.132 GHz, 48GB RAM, 1 Gb/s NIC) with 4 workers.
      const double reftime = 70.169;
      double glbmarks = (gProof->GetRealTime() > 0) ? 1000 * reftime / gProof->GetRealTime() : -1;
      printf(" ROOTMARKS = %.2f (overall: %.2f) ROOT version: %s\t%s@%s\n",
             avgmarks, glbmarks, gROOT->GetVersion(),
             gROOT->GetGitBranch(), gROOT->GetGitCommit());
      // Average from the single tests
      printf("******************************************************************\n");
   }

   // If not PROOF-Lite, stop the daemon used for the test
   if (gProof && !gProof->IsLite() && !extcluster && gROOT->IsBatch()) {
      // Get the manager
      TProofMgr *mgr = gProof->GetManager();
      // Close the instance
      gProof->Close("S");
      delete gProof;
      // Delete the manager
      SafeDelete(mgr);
      // The daemon runs on a port shifted by 1
      if (killXrootdAt(uu.GetPort()+1, "xpdtut") != 0) {
         printf("+++ Warning: test daemon probably still running!\n");
      }
   }

   // Done
   return (failed ? 1 : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Release memory cache associated with the H1 files at 'h1src', if it
/// makes any sense, i.e. are local ...

Int_t PT_H1ReleaseCache(const char *h1src)
{
   if (!h1src || strlen(h1src) <= 0) {
      printf("\n >>> Test failure: src dir undefined\n");
      return -1;
   }

   // If non-local, nothing to do
   if (!gh1local) return 0;

   TString src;
   if (gh1sep == '/') {
      // Loop through the files
      for (Int_t i = 0; i < 4; i++) {
         src.Form("%s/%s", h1src, gh1file[i]);
         ReleaseCache(src.Data());
      }
   } else {
      // Release the zip file ...
      ReleaseCache(h1src);
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that the needed H1 files are available at 'src'
/// If 'src' is "download", the files are download under <tutdir>/h1

Int_t PT_H1AssertFiles(const char *h1src)
{
   if (!h1src || strlen(h1src) <= 0) {
      printf("\n >>> Test failure: src dir undefined\n");
      return -1;
   }

   // Locality
   TUrl u(h1src, kTRUE);
   gh1local = (!strcmp(u.GetProtocol(), "file")) ? kTRUE : kFALSE;

   gh1sep = '/';
   // Special cases
   if (!strncmp(h1src,"download",8)) {
      if (strcmp(h1src,"download")) {
         gh1src = h1src;
         gh1src.ReplaceAll("download=", "");
      }
      if (gh1src.IsNull() || gSystem->AccessPathName(gh1src, kWritePermission))
         gh1src = TString::Format("%s/h1", gtutdir.Data());
      if (gSystem->AccessPathName(gh1src)) {
         if (gSystem->MakeDirectory(gh1src) != 0) {
            printf("\n >>> Test failure: could not create dir %s\n", gh1src.Data());
            return -1;
         }
      }
      // Copy the files now
      Int_t i = 0;
      for (i = 0; i < 4; i++) {
         TString src = TString::Format("http://root.cern.ch/files/h1/%s", gh1file[i]);
         TString dst = TString::Format("%s/%s", gh1src.Data(), gh1file[i]);
         if (!TFile::Cp(src, dst)) {
            printf("\n >>> Test failure: problems retrieving %s\n", src.Data());
            return -1;
         }
         gSystem->RedirectOutput(0, 0, &gRH);
         printf("%d\b", i);
         gSystem->RedirectOutput(glogfile, "a", &gRH);
      }
      gh1local = kTRUE;
      // Done
      gh1ok = kTRUE;
      return 0;
   } else if (TString(h1src).EndsWith(".zip")) {
      // The files are in a zip archive ...
      if (gSystem->AccessPathName(h1src)) {
         printf("\n >>> Test failure: file %s does not exist\n", h1src);
         return -1;
      }
      Int_t i = 0;
      for (i = 0; i < 4; i++) {
         TString src = TString::Format("%s#%s", h1src, gh1file[i]);
         TFile *f = TFile::Open(src);
         if (!f || (f && f->IsZombie())) {
            printf("\n >>> Test failure: file %s not found in archive %s\n", src.Data(), h1src);
            return -1;
         }
         if (guseprogress) {
            gSystem->RedirectOutput(0, 0, &gRH);
            printf("%d\b", i);
            gSystem->RedirectOutput(glogfile, "a", &gRH);
         }
      }
      gh1sep = '#';
   } else {

      // Make sure the files exist at 'src'
      Int_t i = 0;
      for (i = 0; i < 4; i++) {
         TString src = TString::Format("%s/%s", h1src, gh1file[i]);
         if (gSystem->AccessPathName(src)) {
            printf("\n >>> Test failure: file %s does not exist\n", src.Data());
            return -1;
         }
         if (guseprogress) {
            gSystem->RedirectOutput(0, 0, &gRH);
            printf("%d\b", i);
            gSystem->RedirectOutput(glogfile, "a", &gRH);
         }
      }
   }
   gh1src = h1src;

   // Done
   gh1ok = kTRUE;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Release memory cache associated with the event files at 'eventsrc', if it
/// makes any sense, i.e. are local ...

Int_t PT_EventReleaseCache(const char *eventsrc, Int_t nf = 10)
{
   if (!eventsrc || strlen(eventsrc) <= 0) {
      printf("\n >>> Test failure: src dir undefined\n");
      return -1;
   }

   if (nf > 50) {
      printf("\n >>> Test failure: max 50 event files can be checked\n");
      return -1;
   }

   // If non-local, nothing to do
   if (!geventlocal) return 0;

   TString src;
   // Loop through the files
   for (Int_t i = 0; i < nf; i++) {
      src.Form("%s/event_%d.root", eventsrc, i+1);
      ReleaseCache(src.Data());
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that the needed 'event' files are available at 'src'
/// If 'src' is "download", the files are download under <tutdir>/event .
/// By default 10 files are checked; maximum is 50 (idx 1->50 not 0->49).

Int_t PT_EventAssertFiles(const char *eventsrc, Int_t nf = 10)
{
   if (!eventsrc || strlen(eventsrc) <= 0) {
      printf("\n >>> Test failure: src dir undefined\n");
      return -1;
   }

   if (nf > 50) {
      printf("\n >>> Test failure: max 50 event files can be checked\n");
      return -1;
   }

   // Locality
   TUrl u(eventsrc, kTRUE);
   geventlocal = (!strcmp(u.GetProtocol(), "file")) ? kTRUE : kFALSE;

   // Special case
   if (!strncmp(eventsrc,"download",8)) {
      if (strcmp(eventsrc,"download")) {
         geventsrc = eventsrc;
         geventsrc.ReplaceAll("download=", "");
      }
      if (geventsrc.IsNull() || gSystem->AccessPathName(geventsrc, kWritePermission))
         geventsrc = TString::Format("%s/event", gtutdir.Data());
      if (gSystem->AccessPathName(geventsrc)) {
         if (gSystem->MakeDirectory(geventsrc) != 0) {
            printf("\n >>> Test failure: could not create dir %s\n", geventsrc.Data());
            return -1;
         }
      }
      // Copy the files now
      Int_t i = 0;
      for (i = 0; i < nf; i++) {
         TString src = TString::Format("http://root.cern.ch/files/data/event_%d.root", i+1);
         TString dst = TString::Format("%s/event_%d.root", geventsrc.Data(), i+1);
         if (!TFile::Cp(src, dst)) {
            printf("\n >>> Test failure: problems retrieving %s\n", src.Data());
            return -1;
         }
         if (guseprogress) {
            gSystem->RedirectOutput(0, 0, &gRH);
            printf("%d\b", i);
            gSystem->RedirectOutput(glogfile, "a", &gRH);
         }
      }
      geventlocal = kTRUE;
      // Done
      geventok = kTRUE;
      return 0;
   }

   // Make sure the files exist at 'src'
   Int_t i = 0;
   for (i = 0; i < nf; i++) {
      TString src = TString::Format("%s/event_%d.root", eventsrc, i+1);
      if (gSystem->AccessPathName(src)) {
         printf("\n >>> Test failure: file %s does not exist\n", src.Data());
         return -1;
      }
      if (guseprogress) {
         gSystem->RedirectOutput(0, 0, &gRH);
         printf("%d\b", i);
         gSystem->RedirectOutput(glogfile, "a", &gRH);
      }
   }
   geventsrc = eventsrc;

   // Done
   geventok = kTRUE;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that the needed files are available under the specified
/// tutorial directory, setting the relevant variables

Int_t PT_AssertTutorialDir(const char *tutdir)
{
   if (!tutdir || strlen(tutdir) <= 0) {
      printf("\n >>> Test failure: dir undefined\n");
      return -1;
   }

   TString path;
   // Selectors
   gH1Sel.Insert(0, tutdir);
   gSystem->ExpandPathName(gH1Sel);
   if (gSystem->AccessPathName(gH1Sel)) return -1;
   //
   gEventSel.Insert(0, tutdir);
   gSystem->ExpandPathName(gEventSel);
   if (gSystem->AccessPathName(gEventSel)) return -1;
   //
   gEventProcSel.Insert(0, tutdir);
   gSystem->ExpandPathName(gEventProcSel);
   if (gSystem->AccessPathName(gEventProcSel)) return -1;
   //
   gSimpleSel.Insert(0, tutdir);
   gSystem->ExpandPathName(gSimpleSel);
   if (gSystem->AccessPathName(gSimpleSel)) return -1;
   //
   gTestsSel.Insert(0, tutdir);
   gSystem->ExpandPathName(gTestsSel);
   if (gSystem->AccessPathName(gTestsSel)) return -1;
   //
   gNtupleSel.Insert(0, tutdir);
   gSystem->ExpandPathName(gNtupleSel);
   if (gSystem->AccessPathName(gNtupleSel)) return -1;
   //
   gFriendsSel.Insert(0, tutdir);
   gSystem->ExpandPathName(gFriendsSel);
   if (gSystem->AccessPathName(gFriendsSel)) return -1;
   //
   gAuxSel.Insert(0, tutdir);
   gSystem->ExpandPathName(gAuxSel);
   if (gSystem->AccessPathName(gAuxSel)) return -1;

   // Special class
   gProcFileElem.Insert(0, tutdir);
   gSystem->ExpandPathName(gProcFileElem);
   if (gSystem->AccessPathName(gProcFileElem)) return -1;

   // Special files
   gEmptyInclude.Insert(0, tutdir);
   gSystem->ExpandPathName(gEmptyInclude);
   if (gSystem->AccessPathName(gEmptyInclude)) return -1;
   //
   gNtpRndm.Insert(0, tutdir);
   gSystem->ExpandPathName(gNtpRndm);
   if (gSystem->AccessPathName(gNtpRndm)) return -1;

   // Special packages
   gPackEvent.Insert(0, tutdir);
   gSystem->ExpandPathName(gPackEvent);
   if (gSystem->AccessPathName(gPackEvent)) return -1;
   //
   gPack1.Insert(0, tutdir);
   gSystem->ExpandPathName(gPack1);
   if (gSystem->AccessPathName(gPack1)) return -1;
   //
   gPack2.Insert(0, tutdir);
   gSystem->ExpandPathName(gPack2);
   if (gSystem->AccessPathName(gPack2)) return -1;
   //
   gPack3.Insert(0, tutdir);
   gSystem->ExpandPathName(gPack3);
   if (gSystem->AccessPathName(gPack3)) return -1;

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the result of the ProofSimple analysis

Int_t PT_CheckSimple(TQueryResult *qr, Long64_t nevt, Int_t nhist)
{
   if (!qr) {
      printf("\n >>> Test failure: query result not found\n");
      return -1;
   }

   // Make sure the number of processed entries is the one expected
   PutPoint();
   if (qr->GetEntries() != nevt) {
      printf("\n >>> Test failure: wrong number of entries processed: %lld (expected %lld)\n",
             qr->GetEntries(), nevt);
      return -1;
   }

   // Make sure the output list is there
   PutPoint();
   TList *out = qr->GetOutputList();
   if (!out) {
      printf("\n >>> Test failure: output list not found\n");
      return -1;
   }

   // Get the histos
   PutPoint();
   TH1F **hist = new TH1F*[nhist];
   for (Int_t i=0; i < nhist; i++) {
      hist[i] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("h%d",i), out));
      if (!hist[i]) {
         printf("\n >>> Test failure: 'h%d' histo not found\n", i);
         delete[] hist;
         return -1;
      }
   }

   // Check the mean values
   PutPoint();
   for (Int_t i=0; i < nhist; i++) {
      Double_t ave = hist[i]->GetMean();
      Double_t rms = hist[i]->GetRMS();
      if (TMath::Abs(ave) > 5 * rms / TMath::Sqrt(hist[i]->GetEntries())) {
         printf("\n >>> Test failure: 'h%d' histo: mean > 5 * RMS/Sqrt(N)\n", i);
         delete[] hist;
         return -1;
      }
   }

   // Clean up
   delete[] hist;

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the ntuple created by the ProofSimple analysis

Int_t PT_CheckSimpleNtuple(TQueryResult *qr, Long64_t nevt, const char *dsname)
{
   if (!qr) {
      printf("\n >>> Test failure: query result not found\n");
      return -1;
   }

   // Make sure the output list is there
   PutPoint();
   TList *out = qr->GetOutputList();
   if (!out) {
      printf("\n >>> Test failure: output list not found\n");
      return -1;
   }

   // Get the file collection
   PutPoint();
   TFileCollection *fc = dynamic_cast<TFileCollection *>(out->FindObject(dsname));
   if (!fc) {
      printf("\n >>> Test failure: TFileCollection for dataset '%s' not"
             " found in the output list\n", dsname);
      return -1;
   }

   // Check the default tree name
   const char *tname = "/ntuple";
   PutPoint();
   if (!fc->GetDefaultTreeName() || strcmp(fc->GetDefaultTreeName(), tname)) {
      printf("\n >>> Test failure: default tree name does not match (%s != %s)\n",
             fc->GetDefaultTreeName(), tname);
      return -1;
   }

   // Check the number of entries
   PutPoint();
   if (fc->GetTotalEntries(tname) != nevt) {
      printf("\n >>> Test failure: number of entries does not match (%lld != %lld)\n",
             fc->GetTotalEntries(tname), nevt);
      return -1;
   }

   // Check 'pz' histo
   TH1F *hpx = new TH1F("PT_px", "PT_px", 20, -5., 5.);
   PutPoint();
   gProof->DrawSelect(dsname, "px >> PT_px");
   if (TMath::Abs(hpx->GetMean()) > 5 * hpx->GetRMS() / TMath::Sqrt(hpx->GetEntries())) {
      printf("\n >>> Test failure: 'hpx' histo: mean > 5 * RMS/Sqrt(N) (%f,%f)\n",
             hpx->GetMean(), hpx->GetRMS());
      return -1;
   }

   // Check 'pz' histo
   TH1F *hpz = new TH1F("PT_pz", "PT_pz", 20, 0., 20.);
   PutPoint();
   gProof->DrawSelect(dsname, "pz >> PT_pz");
   if (TMath::Abs(hpz->GetMean() - 2.) > 5 * 2. / TMath::Sqrt(hpz->GetEntries())) {
      printf("\n >>> Test failure: 'hpz' histo: (mean - 2) > 5 * RMS/Sqrt(N) (%f,%f)\n",
             hpz->GetMean(), hpz->GetRMS());
      return -1;
   }

   // Check 'random' histo
   TH1F *hpr = new TH1F("PT_rndm", "PT_rndm", 20, 0., 20.);
   PutPoint();
   gProof->DrawSelect(dsname, "random >> PT_rndm");
   if (TMath::Abs(hpr->GetMean() - .5) > 5 * .5 / TMath::Sqrt(hpr->GetEntries())) {
      printf("\n >>> Test failure: 'hpr' histo: (mean - .5) > 5 * RMS/Sqrt(N) (%f,%f)\n",
             hpr->GetMean(), hpr->GetRMS());
      return -1;
   }

   SafeDelete(hpx);
   SafeDelete(hpz);
   SafeDelete(hpr);

   // Clear dsname
   gProof->ClearData(TProof::kDataset |TProof::kForceClear, dsname);

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the result of the H1 analysis

Int_t PT_CheckH1(TQueryResult *qr, Int_t irun = 0)
{
   if (!qr) {
      printf("\n >>> Test failure: output list not found\n");
      return -1;
   }

   // Make sure the number of processed entries is the one expected
   PutPoint();
   Long64_t runEntries[2] = {283813, 7525};
   if (qr->GetEntries() != runEntries[irun]) {
      printf("\n >>> Test failure: wrong number of entries processed: %lld (expected %lld)\n", qr->GetEntries(), runEntries[irun]);
      return -1;
   }

   // Make sure the output list is there
   PutPoint();
   TList *out = qr->GetOutputList();
   if (!out) {
      printf("\n >>> Test failure: output list not found\n");
      return -1;
   }

   // Check the 'hdmd' histo
   PutPoint();
   TH1F *hdmd = dynamic_cast<TH1F*>(out->FindObject("hdmd"));
   if (!hdmd) {
      printf("\n >>> Test failure: 'hdmd' histo not found\n");
      return -1;
   }
   if ((Int_t)(hdmd->GetEntries()) != 7525) {
      printf("\n >>> Test failure: 'hdmd' histo: wrong number"
             " of entries (%d: expected 7525) \n",(Int_t)(hdmd->GetEntries()));
      return -1;
   }
   if (TMath::Abs((hdmd->GetMean() - 0.15512023) / 0.15512023) > 0.001) {
      printf("\n >>> Test failure: 'hdmd' histo: wrong mean"
             " (%f: expected 0.15512023) \n", hdmd->GetMean());
      return -1;
   }

   PutPoint();
   TH2F *h2 = dynamic_cast<TH2F*>(out->FindObject("h2"));
   if (!h2) {
      printf("\n >>> Test failure: 'h2' histo not found\n");
      return -1;
   }
   if ((Int_t)(h2->GetEntries()) != 7525) {
      printf("\n >>> Test failure: 'h2' histo: wrong number"
             " of entries (%d: expected 7525) \n",(Int_t)(h2->GetEntries()));
      return -1;
   }
   if (TMath::Abs((h2->GetMean() - 0.15245688) / 0.15245688) > 0.001) {
      printf("\n >>> Test failure: 'h2' histo: wrong mean"
             " (%f: expected 0.15245688) \n", h2->GetMean());
      return -1;
   }

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the result of the EventProc analysis

Int_t PT_CheckEvent(TQueryResult *qr, const char *pack = "TPacketizer")
{
   if (!qr) {
      printf("\n >>> Test failure: %s: output list not found\n", pack);
      return -1;
   }

   // Make sure the output list is there
   PutPoint();
   TList *out = qr->GetOutputList();
   if (!out) {
      printf("\n >>> Test failure: %s: output list not found\n", pack);
      return -1;
   }

   // Check the 'hdmd' histo
   PutPoint();
   TNamed *nout = dynamic_cast<TNamed *>(out->FindObject("Range_Check"));
   if (!nout) {
      printf("\n >>> Test failure: %s: 'Range_Check' named object not found\n", pack);
      return -1;
   }
   if (strcmp(nout->GetTitle(), "OK")) {
      printf("\n >>> Test failure: %s: 'Range_Check': wrong result: %s \n", pack, nout->GetTitle());
      return -1;
   }

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the result of the ProofNtuple analysis

Int_t PT_CheckNtuple(TQueryResult *qr, Long64_t nevt)
{
   if (!qr) {
      printf("\n >>> Test failure: query result not found\n");
      return -1;
   }

   // Make sure the number of processed entries is the one expected
   PutPoint();
   if (qr->GetEntries() != nevt) {
      printf("\n >>> Test failure: wrong number of entries processed: %lld (expected %lld)\n",
             qr->GetEntries(), nevt);
      return -1;
   }

   // Make sure the output list is there
   PutPoint();
   TList *out = qr->GetOutputList();
   if (!out) {
      printf("\n >>> Test failure: output list not found\n");
      return -1;
   }

   // Get the ntuple form the file
   TProofOutputFile *pof = dynamic_cast<TProofOutputFile*>(out->FindObject("SimpleNtuple.root"));
   if (!pof) {
      printf("\n >>> Test failure: TProofOutputFile not found in the output list\n");
      return -1;
   }

   // Get the file full path
   TString outputFile(pof->GetOutputFileName());
   TString outputName(pof->GetName());
   outputName += ".root";

   // Read the ntuple from the file
   TFile *f = TFile::Open(outputFile);
   if (!f || (f && f->IsZombie())) {
      printf("\n >>> Test failure: could not open file: %s", outputFile.Data());
      return -1;
   }

   // Get the ntuple
   PutPoint();
   TNtuple *ntp = dynamic_cast<TNtuple *>(f->Get("ntuple"));
   if (!ntp) {
      printf("\n >>> Test failure: 'ntuple' not found\n");
      return -1;
   }

   // Check the ntuple content by filling some histos
   TH1F *h1s[3] = {0};
   h1s[0] = new TH1F("h1_1", "3*px+2 with px**2+py**2>1", 50, -15., 15.);
   h1s[1] = new TH1F("h1_2", "2*px+2 with pz>2", 50, -10., 10.);
   h1s[2] = new TH1F("h1_3", "1.3*px+2 with (px^2+py^2>4) && py>0", 50, 0., 8.);
   Float_t px, py, pz;
   Float_t *pp = ntp->GetArgs();
   Long64_t ent = 0;
   while (ent < ntp->GetEntries()) {
      ntp->GetEntry(ent);
      px = pp[0];
      py = pp[1];
      pz = pp[2];
      // Fill the histos
      if (px*px+py*py > 1.) h1s[0]->Fill(3.*px + 2.);
      if (pz > 2.) h1s[1]->Fill(2.*px + 2.);
      if (px*px+py*py > 4. && py > 0.) h1s[2]->Fill(1.3*px + 2.);
      // Go next
      ent++;
   }

   Int_t rch1s = 0;
   TString emsg;
   // Check the histogram entries and mean values
   Int_t hent[3] = { 620, 383, 72};
   Double_t hmea[3] = { 1.992, 2.062 , 3.126};
   for (Int_t i = 0; i < 3; i++) {
      if ((Int_t)(h1s[i]->GetEntries()) != hent[i]) {
         emsg.Form("'%s' histo: wrong number of entries (%d: expected %d)",
                   h1s[i]->GetName(), (Int_t)(h1s[i]->GetEntries()), hent[i]);
         rch1s = -1;
         break;
      }
      if (TMath::Abs((h1s[i]->GetMean() - hmea[i]) / hmea[i]) > 0.001) {
         emsg.Form("'%s' histo: wrong mean (%f: expected %f)",
                   h1s[i]->GetName(), h1s[i]->GetMean(), hmea[i]);
         rch1s = -1;
         break;
      }
   }

   // Cleanup
   for (Int_t i = 0; i < 3; i++) delete h1s[i];
   f->Close();
   delete f;

   // Check the result
   if (rch1s != 0) {
      printf("\n >>> Test failure: %s\n", emsg.Data());
      return -1;
   }

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the result of the ProofNtuple analysis creating a dataset
/// Uses and check also TProofDraw

Int_t PT_CheckDataset(TQueryResult *qr, Long64_t nevt)
{
   if (!qr) {
      printf("\n >>> Test failure: query result not found\n");
      return -1;
   }

   // Make sure the number of processed entries is the one expected
   PutPoint();
   if (qr->GetEntries() != nevt) {
      printf("\n >>> Test failure: wrong number of entries processed: %lld (expected %lld)\n",
             qr->GetEntries(), nevt);
      return -1;
   }

   const char *dsname = "testNtuple";
   // Make sure that the dataset exists
   PutPoint();
   if (!gProof->ExistsDataSet(dsname)) {
      gProof->ShowDataSets();
      printf("\n >>> Test failure: dataset '%s' not found in the repository\n", dsname);
      return -1;
   }
   // ... and that the default tree is 'ntuple'
   gProof->SetDataSetTreeName(dsname, "ntuple");

   // Create the histos
   TH1F *h1s[3] = {0};
   h1s[0] = new TH1F("h1s0", "3*px+2 with px**2+py**2>1", 50, -15., 15.);
   h1s[1] = new TH1F("h1s1", "2*px+2 with pz>2", 50, -10., 10.);
   h1s[2] = new TH1F("h1s2", "1.3*px+2 with (px^2+py^2>4) && py>0", 50, 0., 8.);

   // Fill the histos using TProofDraw
   PutPoint();
   {  SwitchProgressGuard spg;
      gProof->DrawSelect(dsname, "3*px+2 >> h1s0","px**2+py**2>1");
      PutPoint();
      gProof->DrawSelect(dsname, "2*px+2 >> h1s1","pz>2");
      PutPoint();
      gProof->DrawSelect(dsname, "1.3*px+2 >> h1s2","(px^2+py^2>4) && py>0");
   }

   Int_t rch1s = 0;
   TString emsg;
   // Check the histogram entries and mean values
   Float_t hent[3] = { .607700, .364900, .065100};
   Double_t hmea[3] = { 2.022, 2.046 , 3.043};
   Double_t prec = 10. / TMath::Sqrt(nevt);  // ~10 sigma ... conservative
   for (Int_t i = 0; i < 3; i++) {
      Double_t ent = h1s[i]->GetEntries();
      if (TMath::Abs(ent - hent[i] * nevt) / ent > prec) {
         emsg.Form("'%s' histo: wrong number"
               " of entries (%lld: expected %lld)",
                h1s[i]->GetName(), (Long64_t) ent, (Long64_t)(hent[i] *nevt));
         rch1s = -1;
         break;
      }
      Double_t mprec = 5 * h1s[i]->GetRMS() / TMath::Sqrt(h1s[i]->GetEntries()) ;
      if (TMath::Abs((h1s[i]->GetMean() - hmea[i]) / hmea[i]) > mprec ) {
         emsg.Form("'%s' histo: wrong mean (%f: expected %f - RMS: %f)",
                h1s[i]->GetName(), h1s[i]->GetMean(), hmea[i], h1s[i]->GetRMS());
         rch1s = -1;
         break;
      }
   }

   // Cleanup
   for (Int_t i = 0; i < 3; i++) delete h1s[i];

   // Check the result
   if (rch1s != 0) {
      printf("\n >>> Test failure: %s\n", emsg.Data());
      return -1;
   }

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the result of the ProofFriends analysis

Int_t PT_CheckFriends(TQueryResult *qr, Long64_t nevt, bool withfriends)
{
   if (!qr) {
      printf("\n >>> Test failure: query result not found\n");
      return -1;
   }

   // Make sure the number of processed entries is the one expected
   PutPoint();
   if (qr->GetEntries() != nevt) {
      printf("\n >>> Test failure: wrong number of entries processed: %lld (expected %lld)\n",
             qr->GetEntries(), nevt);
      return -1;
   }

   // Make sure the output list is there
   PutPoint();
   TList *out = qr->GetOutputList();
   if (!out) {
      printf("\n >>> Test failure: output list not found\n");
      return -1;
   }

   TString emsg;
   // Check the histogram entries and mean values
   Int_t nchk = (withfriends) ? 4 : 2;
   Int_t rchs = 0;
   const char *hnam[4] = { "histo1", "histo2", "histo3", "histo4" };
   const char *hcls[4] = { "TH2F", "TH1F", "TH1F", "TH2F" };
   Float_t hent[4] = { 1., .227, .0260, .0260};
   TObject *o = 0;
   Double_t ent = -1;
   Double_t prec = 1. / TMath::Sqrt(nevt);
   for (Int_t i = 0; i < nchk; i++) {
      if (!(o = out->FindObject(hnam[i]))) {
         emsg.Form("object '%s' not found", hnam[i]);
         rchs = -1;
         break;
      }
      if (strcmp(o->IsA()->GetName(), hcls[i])) {
         emsg.Form("object '%s' is not '%s'", hnam[i], hcls[i]);
         rchs = -1;
         break;
      }
      if (!strcmp(hcls[i], "TH1F")) {
         ent = ((TH1F *)o)->GetEntries();
      } else {
         ent = ((TH2F *)o)->GetEntries();
      }
      if (TMath::Abs(ent - hent[i] * nevt) / (Double_t)ent > prec) {
         emsg.Form("'%s' histo: wrong number of entries (%lld: expected %lld) \n",
                   o->GetName(), (Long64_t) ent, (Long64_t) (hent[i] * nevt));
         rchs = -1;
         break;
      }
   }

   if (rchs != 0) {
      printf("\n >>> Test failure: %s\n", emsg.Data());
      return -1;
   }

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Test session opening

Int_t PT_Open(void *args, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   PT_Open_Args_t *PToa = (PT_Open_Args_t *)args;
   if (!PToa) {
      printf("\n >>> Test failure: invalid arguments: %p\n", args);
      return -1;
   }

   // Temp dir for PROOF tutorials
   PutPoint();
#if defined(R__MACOSX)
   // Force '/tmp' under macosx, to avoid problems with lengths and symlinks
   TString tmpdir("/tmp"), uspid;
#else
   TString tmpdir(gSystem->TempDirectory()), uspid;
#endif
   UserGroup_t *ug = gSystem->GetUserInfo(gSystem->GetUid());
   if (!ug) {
      printf("\n >>> Test failure: could not get user info");
      return -1;
   }
   if (!tmpdir.EndsWith(ug->fUser.Data())) {
      uspid.Form("/%s/%d", ug->fUser.Data(), gSystem->GetPid());
      delete ug;
   } else {
      uspid.Form("/%d", gSystem->GetPid());
   }
   tmpdir += uspid;
#if !defined(R__MACOSX)
   gtutdir.Form("%s/.proof-tutorial", tmpdir.Data());
#else
   gtutdir.Form("%s/.proof", tmpdir.Data());
#endif
   if (gSystem->AccessPathName(gtutdir)) {
      if (gSystem->mkdir(gtutdir, kTRUE) != 0) {
         printf("\n >>> Test failure: could not assert/create the temporary directory"
                " for the tutorial (%s)", gtutdir.Data());
         return -1;
      }
   }

   // String to initialize the dataset manager
   TString dsetmgrstr;
   dsetmgrstr.Form("file dir:%s/datasets opt:-Cq:As:Sb:", gtutdir.Data());
   gEnv->SetValue("Proof.DataSetManager", dsetmgrstr.Data());

   // String to initialize the package dir
   TString packdir;
   packdir.Form("%s/packages", gtutdir.Data());
   gEnv->SetValue("Proof.PackageDir", packdir.Data());

   // Get the PROOF Session
   PutPoint();
   TProof *p = getProof(PToa->url, PToa->nwrks, gtutdir.Data(), "force", gDynamicStartup, kTRUE);
   if (!p || !(p->IsValid())) {
      printf("\n >>> Test failure: could not start the session\n");
      return -1;
   }

   // Re-check locality: if the logged user name is different from the local one, we may
   // not have all the rights we need, so we go no-local
   if (gLocalCluster) {
      UserGroup_t *pw = gSystem->GetUserInfo();
      if (pw) {
         if (strcmp(pw->fUser, p->GetUser())) gLocalCluster = kFALSE;
         delete pw;
      }
   }

   // Check if it is in dynamic startup mode
   Int_t dyn = 0;
   p->GetRC("Proof.DynamicStartup", dyn);
   if (dyn != 0) gDynamicStartup = kTRUE;

   // Set debug level, if required
   if (gverbose > 1) {
      Int_t debugscope = getDebugEnum(gverbproof.Data());
      p->SetLogLevel(gverbose, debugscope);
   }

   PutPoint();
   if (PToa->nwrks > 0 && p->GetParallel() != PToa->nwrks) {
      printf("\n >>> Test failure: number of workers different from requested\n");
      return -1;
   }

   // Clear Cache
   p->ClearCache();

   // Get some useful info about the cluster (the sandbox dir ...)
   gSystem->RedirectOutput(0, 0, &gRH);
   TString testPrint(TString::Format("%s/testPrint.log", gtutdir.Data()));
   gSystem->RedirectOutput(testPrint, "w", &gRHAdmin);
   gProof->Print();
   gSystem->RedirectOutput(0, 0, &gRHAdmin);
   gSystem->RedirectOutput(glogfile, "a", &gRH);
   TMacro macroPrint(testPrint);
   TObjString *os = macroPrint.GetLineWith("Working directory:");
   if (!os) {
      printf("\n >>> Test failure: problem parsing output from Print()\n");
      return -1;
   }
   Int_t from = strlen("Working directory:") + 1;
   while (os->GetString().Tokenize(gsandbox, from, " ")) {
      if (!gsandbox.IsNull()) break;
   }
   if (gsandbox.IsNull()) {
      printf("\n >>> Test failure: no sandbox dir found\n");
      return -1;
   }
   gsandbox = gSystem->GetDirName(gsandbox);
   gsandbox = gSystem->GetDirName(gsandbox);
   PutPoint();

   // Fill times
   PT_GetLastTimes(tt);

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Test log retrieving

Int_t PT_GetLogs(void *args, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   PT_Open_Args_t *PToa = (PT_Open_Args_t *)args;
   if (!PToa) {
      printf("\n >>> Test failure: invalid arguments: %p\n", args);
      return -1;
   }

   PutPoint();
   TProofLog *pl = TProof::Mgr(PToa->url)->GetSessionLogs();
   if (!pl) {
      printf("\n >>> Test failure: could not get the logs from last session\n");
      return -1;
   }

   PutPoint();
   if (PToa->nwrks > 0 && pl->GetListOfLogs()->GetSize() != (PToa->nwrks + 1)) {
      printf("\n >>> Test failure: number of logs different from workers of workers + 1\n");
      return -1;
   }

   // Fill times
   PT_GetLastTimes(tt);

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Test run for the ProofSimple analysis (see tutorials)

Int_t PT_Simple(void *opts, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   PT_Option_t *ptopt = (PT_Option_t *) opts;

   // Setup submergers if required
   if (ptopt && ptopt->fOne > 0) {
      gProof->SetParameter("PROOF_UseMergers", 0);
   }
   // Setup save-to-file, if required
   TString opt = (ptopt && ptopt->fTwo > 0) ? "stf" : "" ;

   // Define the number of events and histos
   Long64_t nevt = 1000000;
   Int_t nhist = 16;
   // The number of histograms is added as parameter in the input list
   gProof->SetParameter("ProofSimple_NHist", (Long_t)nhist);

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      gProof->Process(gSimpleSel.Data(), nevt, opt);
      gTimer.Stop();
   }

   // Count
   gSimpleCnt++;
   gSimpleTime += gTimer.RealTime();

   // Remove any setting related to submergers
   gProof->DeleteParameters("PROOF_UseMergers");

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckSimple(gProof->GetQueryResult(), nevt, nhist);
}

////////////////////////////////////////////////////////////////////////////////
/// Test output handling via file using ProofSimple (see tutorials)

Int_t PT_OutputHandlingViaFile(void *opts, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }
   // Not yet supported in dynamic mode
   if (gDynamicStartup) {
      return 1;
   }
   PutPoint();

   PT_Option_t *ptopt = (PT_Option_t *) opts;

   // Setup submergers if required
   if (ptopt && ptopt->fOne > 0) {
      gProof->SetParameter("PROOF_UseMergers", 0);
   }
   // Setup save-to-file, if required
   TString opt = (ptopt && ptopt->fTwo > 0) ? "stf" : "" ;

   // Define the number of events and histos
   Long64_t nevt = 1000000 * gProof->GetParallel();
   Int_t nhist = 16;
   // The number of histograms is added as parameter in the input list
   gProof->SetParameter("ProofSimple_NHist", (Long_t)nhist);

   // Merged file pptions to be tested
   const char *testopt[4] = { "stf", "of=proofsimple.root", "of=proofsimple.root;stf",
                                     "of=master:proofsimple.root" };

   for (Int_t i = 0; i < 4; i++) {
      // Clear the list of query results
      if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

      // Save results to file 'proofsimple.root'
      PutPoint();
      {  SwitchProgressGuard spg;
         gTimer.Start();
         gProof->Process(gSimpleSel.Data(), nevt, testopt[i]);
         gTimer.Stop();
      }
      if (PT_CheckSimple(gProof->GetQueryResult(), nevt, nhist) != 0) {
          printf("\n >>> Test failure: output handling via file: option '%s'\n", testopt[i]);
         return -1;
      }
      // Count
      gSimpleCnt++;
      gSimpleTime += gTimer.RealTime();
      // Remove file
      gSystem->Unlink("proofsimple.root");
   }

   // Test dataset creationg with a ntuple
   const char *dsname = "PT_ds_proofsimple";
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();
   if (gProof->ExistsDataSet(dsname)) gProof->RemoveDataSet(dsname);

   // We want the ntuple
   gProof->SetParameter("ProofSimple_Ntuple", "");

   // Save results to file 'proofsimple.root'
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      gProof->Process(gSimpleSel.Data(), nevt, TString::Format("ds=%s|V", dsname));
      gTimer.Stop();
   }
   if (!gProof->ExistsDataSet(dsname)) {
      printf("\n >>> Test failure: output handling via file: dataset '%s' not created\n", dsname);
      return -1;
   }

   // Remove any setting related to submergers
   gProof->DeleteParameters("PROOF_UseMergers");

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckSimpleNtuple(gProof->GetQueryResult(), nevt, dsname);
}

////////////////////////////////////////////////////////////////////////////////
/// Test run for the H1 analysis as a chain reading the data from HTTP

Int_t PT_H1Http(void *, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   // Set/unset the parallel unzip flag
   AssertParallelUnzip();

   // Create the chain
   PutPoint();
   TChain *chain = new TChain("h42");

   // Assert the files, if needed
   if (!gh1ok) {
      if (PT_H1AssertFiles(gh1src.Data()) != 0) {
         gProof->SetPrintProgress(0);
         printf("\n >>> Test failure: could not assert the H1 files\n");
         return -1;
      }
   }
   Int_t i = 0;
   for (i = 0; i < 4; i++) {
      chain->Add(TString::Format("%s%c%s", gh1src.Data(), gh1sep, gh1file[i]));
   }

   // Clear associated memory cache
   if (gClearCache && PT_H1ReleaseCache(gh1src.Data()) != 0) {
      gProof->SetPrintProgress(0);
      printf("\n >>> Test failure: could not clear memory cache for the H1 files\n");
      return -1;
   }

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Process
   PutPoint();
   chain->SetProof();
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      chain->Process(gH1Sel.Data());
      gTimer.Stop();
   }
   gProof->RemoveChain(chain);

   // Count
   gH1Cnt++;
   gH1Time += gTimer.RealTime();

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckH1(gProof->GetQueryResult());
}

////////////////////////////////////////////////////////////////////////////////
/// Test run for the H1 analysis as a file collection reading the data from HTTP

Int_t PT_H1FileCollection(void *arg, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   // Set/unset the parallel unzip flag
   AssertParallelUnzip();

   // Are we asked to change the packetizer strategy?
   if (arg) {
      PT_Packetizer_t *strategy = (PT_Packetizer_t *)arg;
      if (strcmp(strategy->fName, "TPacketizer")) {
         gProof->SetParameter("PROOF_Packetizer", strategy->fName);
      } else {
         if (strategy->fType != 1)
            gProof->SetParameter("PROOF_PacketizerStrategy", strategy->fType);
      }
   }

   // Create the file collection
   PutPoint();
   TFileCollection *fc = new TFileCollection("h42");

   // Assert the files, if needed
   if (!gh1ok) {
      if (PT_H1AssertFiles(gh1src.Data()) != 0) {
         printf("\n >>> Test failure: could not assert the H1 files\n");
         return -1;
      }
   }
   Int_t i = 0;
   for (i = 0; i < 4; i++) {
      fc->Add(new TFileInfo(TString::Format("%s%c%s", gh1src.Data(), gh1sep, gh1file[i])));
   }

   // Clear associated memory cache
   if (gClearCache && PT_H1ReleaseCache(gh1src.Data()) != 0) {
      gProof->SetPrintProgress(0);
      printf("\n >>> Test failure: could not clear memory cache for the H1 files\n");
      return -1;
   }

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      gProof->Process(fc, gH1Sel.Data());
      gTimer.Stop();
   }

   // Restore settings
   gProof->DeleteParameters("PROOF_Packetizer");
   gProof->DeleteParameters("PROOF_PacketizerStrategy");

   // Count
   gH1Cnt++;
   gH1Time += gTimer.RealTime();

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckH1(gProof->GetQueryResult());
}

////////////////////////////////////////////////////////////////////////////////
/// Test run for the H1 analysis as a named dataset reading the data from HTTP

Int_t PT_H1DataSet(void *, RunTimes &tt)
{
   // Checking arguments
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }
   // Not yet supported for PROOF-Lite
   if (gSkipDataSetTest) {
      return 1;
   }
   PutPoint();

   // Set/unset the parallel unzip flag
   AssertParallelUnzip();

   // Name for the target dataset
   const char *dsname = "h1dset";

   // Clear associated memory cache
   if (gClearCache && PT_H1ReleaseCache(gh1src.Data()) != 0) {
      gProof->SetPrintProgress(0);
      printf("\n >>> Test failure: could not clear memory cache for the H1 files\n");
      return -1;
   }

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Process the dataset by name
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      gProof->Process(dsname, gH1Sel.Data());
      gTimer.Stop();
   }

   // Count
   gH1Cnt++;
   gH1Time += gTimer.RealTime();

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckH1(gProof->GetQueryResult());
}

////////////////////////////////////////////////////////////////////////////////
/// Test run for the H1 analysis as a named dataset reading the data from HTTP

Int_t PT_H1MultiDataSet(void *, RunTimes &tt)
{
   // Checking arguments
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }
   // Not yet supported for PROOF-Lite
   if (gSkipDataSetTest) {
      return 1;
   }
   PutPoint();

   // Set/unset the parallel unzip flag
   AssertParallelUnzip();

   // Name for the target dataset
   const char *dsname = "h1dseta h1dsetb";

   // Clear associated memory cache
   if (gClearCache && PT_H1ReleaseCache(gh1src.Data()) != 0) {
      gProof->SetPrintProgress(0);
      printf("\n >>> Test failure: could not clear memory cache for the H1 files\n");
      return -1;
   }

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Process the dataset by name
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      gProof->Process(dsname, gH1Sel.Data());
      gTimer.Stop();
   }

   // Count
   gH1Cnt++;
   gH1Time += gTimer.RealTime();

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckH1(gProof->GetQueryResult());
}

////////////////////////////////////////////////////////////////////////////////
/// Test run using the H1 analysis for the multi-dataset functionality and
/// entry-lists

Int_t PT_H1MultiDSetEntryList(void *, RunTimes &tt)
{
   // Checking arguments
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }
   // Not yet supported for PROOF-Lite
   if (gSkipDataSetTest) {
      return 1;
   }
   PutPoint();

   // Set/unset the parallel unzip flag
   AssertParallelUnzip();

   // Multiple dataset used to create the entry list
   TString dsname("h1dseta|h1dsetb");

   // Clear associated memory cache
   if (gClearCache && PT_H1ReleaseCache(gh1src.Data()) != 0) {
      gProof->SetPrintProgress(0);
      printf("\n >>> Test failure: could not clear memory cache for the H1 files\n");
      return -1;
   }

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Entry-list creation run
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      gProof->Process(dsname, gH1Sel.Data(), "fillList=elist.root");
      gTimer.Stop();
   }

   // Cleanup entry-list from the input list
   TIter nxi(gProof->GetInputList());
   TObject *o = 0;
   while ((o = nxi())) {
      if (!strncmp(o->GetName(), "elist", 6) || !strcmp(o->GetName(), "fillList")) {
         gProof->GetInputList()->Remove(o);
         delete o;
      }
   }

   // Count
   gH1Cnt++;
   gH1Time += gTimer.RealTime();

   // Clear associated memory cache
   if (gClearCache && PT_H1ReleaseCache(gh1src.Data()) != 0) {
      gProof->SetPrintProgress(0);
      printf("\n >>> Test failure: could not clear memory cache for the H1 files\n");
      return -1;
   }

   // Run using the entrylist
   dsname = "h1dseta<<elist.root h1dsetb?enl=elist.root";
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      gProof->Process(dsname, gH1Sel.Data());
      gTimer.Stop();
   }

   // Unlink the entry list file
   gSystem->Unlink("elist.root");

   // Cleanup entry-list from the input list
   nxi.Reset();
   while ((o = nxi())) {
      if (!strncmp(o->GetName(), "elist", 6)) {
         gProof->GetInputList()->Remove(o);
         delete o;
      }
   }

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckH1(gProof->GetQueryResult(), 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Test dataset registration, verification, usage, removal.
/// Use H1 analysis files on HTTP as example

Int_t PT_DataSets(void *, RunTimes &tt)
{
   // Checking arguments
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }
   // Not yet supported for PROOF-Lite
   if (gSkipDataSetTest) {
      return 1;
   }
   PutPoint();

   // Cleanup the area
   PutPoint();
   TMap *dsm = gProof->GetDataSets();
   if (!dsm) {
      printf("\n >>> Test failure: could not retrieve map of datasets (even empty)!\n");
      return -1;
   }
   if (dsm->GetSize() > 0) {
      // Remove the datasets already registered
      TIter nxd(dsm);
      TObjString *os = 0;
      while ((os = (TObjString *)nxd())) {
         gProof->RemoveDataSet(os->GetName());
      }
      // Check the result
      delete dsm;
      dsm = gProof->GetDataSets();
      if (!dsm || dsm->GetSize() > 0) {
         printf("\n >>> Test failure: could not cleanup the dataset area! (%p)\n", dsm);
         delete dsm;
         return -1;
      }
   }
   delete dsm;

   // Create the file collection
   PutPoint();
   TFileCollection *fc = new TFileCollection();
   TFileCollection *fca = new TFileCollection();
   TFileCollection *fcb = new TFileCollection();

   // Assert the files, if needed
   if (!gh1ok) {
      if (PT_H1AssertFiles(gh1src.Data()) != 0) {
         printf("\n >>> Test failure: could not assert the H1 files\n");
         return -1;
      }
   }
   Int_t i = 0;
   for (i = 0; i < 4; i++) {
      fc->Add(new TFileInfo(TString::Format("%s%c%s", gh1src.Data(), gh1sep, gh1file[i])));
      if (i < 2) {
         fca->Add(new TFileInfo(TString::Format("%s%c%s", gh1src.Data(), gh1sep, gh1file[i])));
      } else {
         fcb->Add(new TFileInfo(TString::Format("%s%c%s", gh1src.Data(), gh1sep, gh1file[i])));
      }
   }
   fc->Update();
   fca->Update();
   fcb->Update();

   // Name for this dataset
   const char *dsname = "h1dset";
   const char *dsnamea = "h1dseta";
   const char *dsnameb = "h1dsetb";

   // Register the dataset
   PutPoint();
   gProof->RegisterDataSet(dsname, fc);
   gProof->RegisterDataSet(dsnamea, fca);
   gProof->RegisterDataSet(dsnameb, fcb);
   // Check the result
   dsm = gProof->GetDataSets();
   if (!dsm || dsm->GetSize() != 3) {
      printf("\n >>> Test failure: could not register '%s,%s,%s' (%p)\n",
             dsname, dsnamea, dsnameb, dsm);
      delete dsm;
      return -1;
   }
   delete dsm;

   // Test removal
   PutPoint();
   gProof->RemoveDataSet(dsname);
   gProof->RemoveDataSet(dsnamea);
   gProof->RemoveDataSet(dsnameb);
   // Check the result
   dsm = gProof->GetDataSets();
   if (!dsm || dsm->GetSize() != 0) {
      printf("\n >>> Test failure: could not cleanup '%s,%s,%s' (%p)\n",
             dsname, dsnamea, dsnameb, dsm);
      delete dsm;
      return -1;
   }
   delete dsm;

   // Re-register the dataset
   PutPoint();
   gProof->RegisterDataSet(dsname, fc);
   gProof->RegisterDataSet(dsnamea, fca);
   gProof->RegisterDataSet(dsnameb, fcb);
   // Check the result
   dsm = gProof->GetDataSets();
   if (!dsm || dsm->GetSize() != 3) {
      printf("\n >>> Test failure: could not re-register '%s,%s,%s' (%p)\n",
             dsname, dsnamea, dsnameb, dsm);
      delete dsm;
      return -1;
   }
   delete dsm;

   // Verify the dataset
   PutPoint();
   if (gProof->VerifyDataSet(dsname) != 0) {
      printf("\n >>> Test failure: could not verify '%s'!\n", dsname);
      return -1;
   }
   if (gProof->VerifyDataSet(dsnamea) != 0) {
      printf("\n >>> Test failure: could not verify '%s'!\n", dsnamea);
      return -1;
   }
   if (gProof->VerifyDataSet(dsnameb) != 0) {
      printf("\n >>> Test failure: could not verify '%s'!\n", dsnameb);
      return -1;
   }
   gProof->ShowDataSets();

   // Remove the file collection
   delete fc;
   delete fca;
   delete fcb;

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Test package clearing, uploading, enabling, removal.
/// Use event.par as example.

Int_t PT_Packages(void *, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   // Cleanup the area
   PutPoint();
   TList *packs = gProof->GetListOfPackages();
   if (!packs) {
      printf("\n >>> Test failure: could not retrieve list of packages (even empty)!\n");
      return -1;
   }
   if (packs->GetSize() > 0) {
      // Remove the packages already available
      gProof->ClearPackages();
      // Check the result
      packs = gProof->GetListOfPackages();
      if (!packs || packs->GetSize() > 0) {
         printf("\n >>> Test failure: could not cleanup the package area!\n");
         return -1;
      }
   }

   // Name and location for this package
   const char *pack = "event";

   // Upload the package
   PutPoint();
   gProof->UploadPackage(gPackEvent);
   // Check the result
   packs = gProof->GetListOfPackages();
   if (!packs || packs->GetSize() != 1) {
      printf("\n >>> Test failure: could not upload '%s'!\n", gPackEvent.Data());
      return -1;
   }

   // Test cleanup
   PutPoint();
   gProof->ClearPackage(pack);
   // Check the result
   packs = gProof->GetListOfPackages();
   if (!packs || packs->GetSize() != 0) {
      printf("\n >>> Test failure: could not cleanup '%s'!\n", pack);
      return -1;
   }

   // Re-upload the package
   PutPoint();
   gProof->UploadPackage(gPackEvent);
   // Check the result
   packs = gProof->GetListOfPackages();
   if (!packs || packs->GetSize() != 1) {
      printf("\n >>> Test failure: could not re-upload '%s'!\n", gPackEvent.Data());
      return -1;
   }

   // Enable the package
   PutPoint();
   gProof->EnablePackage(pack);
   // Check the result
   packs = gProof->GetListOfEnabledPackages();
   if (!packs || packs->GetSize() != 1) {
      printf("\n >>> Test failure: could not enable '%s'!\n", pack);
      return -1;
   }

   // Fill times
   PT_GetLastTimes(tt);

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Test run for the ProofEvent analysis (see tutorials)

Int_t PT_Event(void *, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   // Define the number of events
   Long64_t nevt = 100000;

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gProof->Process(gEventSel.Data(), nevt);
   }

   // Make sure the query result is there
   PutPoint();
   TQueryResult *qr = 0;
   if (!(qr = gProof->GetQueryResult())) {
      printf("\n >>> Test failure: query result not found\n");
      return -1;
   }

   // Make sure the number of processed entries is the one expected
   PutPoint();
   if (qr->GetEntries() != nevt) {
      printf("\n >>> Test failure: wrong number of entries processed: %lld (expected %lld)\n",
             qr->GetEntries(), nevt);
      return -1;
   }

   // Make sure the output list is there
   PutPoint();
   if (!(gProof->GetOutputList())) {
      printf("\n >>> Test failure: output list not found\n");
      return -1;
   }

   // Check the 'histo'
   PutPoint();
   TH1F *histo = dynamic_cast<TH1F*>(gProof->GetOutputList()->FindObject("histo"));
   if (!histo) {
      printf("\n >>> Test failure: 'histo' not found\n");
      return -1;
   }

   // Check the mean values
   Double_t ave = histo->GetMean();
   Double_t rms = histo->GetRMS();
   if (TMath::Abs(ave - 50) > 10 * rms / TMath::Sqrt(histo->GetEntries())) {
      printf("\n >>> Test failure: 'histo': mean > 5 * RMS/Sqrt(N)\n");
      return -1;
   }

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Test input data functionality

Int_t PT_InputData(void *, RunTimes &tt)
{
   // Checking arguments
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }
   PutPoint();

   // Create the test information to be send via input and retrieved
   TH1F *h1 = new TH1F("h1data","Input data from file",100,-5.,5.);
   h1->FillRandom("gaus", 1000);
   TList *h1list = new TList;
   h1list->SetName("h1list");
   h1list->SetOwner(kTRUE);
   h1list->Add(h1);
   h1list->Add(new TParameter<Double_t>("h1avg", h1->GetMean()));
   h1list->Add(new TParameter<Double_t>("h1rms", h1->GetRMS()));
   TString datafile = glogfile;
   datafile += ("_h1data.root");
   TFile *f = TFile::Open(datafile, "RECREATE");
   if (!f) {
      printf("\n >>> Test failure: could not open file for input data\n");
      return -1;
   }
   f->cd();
   h1list->Write(0, TObject::kSingleKey, 0);
   f->Close();
   gProof->SetInputDataFile(datafile.Data());

   // Histo to be sent from memory
   TH1F *h2 = new TH1F("h2data","Input data from memory",100,-5.,5.);
   h2->FillRandom("gaus", 1000);
   TList *h2list = new TList;
   h2list->SetName("h2list");
   h2list->SetOwner(kTRUE);
   h2list->Add(h2);
   h2list->Add(new TParameter<Double_t>("h2avg", h2->GetMean()));
   h2list->Add(new TParameter<Double_t>("h2rms", h2->GetRMS()));
   gProof->AddInputData(h2list);

   // Normal input parameter
   gProof->AddInput(new TNamed("InputObject", glogfile.Data()));

   // Type of test
   gProof->AddInput(new TNamed("ProofTests_Type", "InputData"));

   // Define the number of events
   Long64_t nevt = 1;

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gProof->Process(gTestsSel.Data(), nevt);
   }

   // Cleanup
   gSystem->Unlink(datafile.Data());
   gProof->ClearInputData(h1list);
   gProof->ClearInputData(h2list);
   delete h1list;
   delete h2list;

   // Make sure the query result is there
   PutPoint();
   TQueryResult *qr = 0;
   if (!(qr = gProof->GetQueryResult())) {
      printf("\n >>> Test failure: query result not found\n");
      return -1;
   }

   // Make sure the output list is there
   PutPoint();
   if (!(gProof->GetOutputList())) {
      printf("\n >>> Test failure: output list not found\n");
      return -1;
   }

   // Check the 'histo's
   PutPoint();
   TH1I *stat = dynamic_cast<TH1I*>(gProof->GetOutputList()->FindObject("TestStat"));
   if (!stat) {
      printf("\n >>> Test failure: 'TestStat' histo not found\n");
      return -1;
   }

   // Test how many workers got everything successfully
   Int_t nw = (Int_t) stat->GetBinContent(1);
   PutPoint();

   if (TMath::Abs(stat->GetBinContent(2) - nw) > .1) {
      printf("\n >>> Test failure: histo 'h1' not correctly received on all workers (%.0f/%d)\n",
             stat->GetBinContent(2), nw);
      return -1;
   }
   if (TMath::Abs(stat->GetBinContent(3) - nw) > .1) {
      printf("\n >>> Test failure: histo 'h2' not correctly received on all workers (%.0f/%d)\n",
             stat->GetBinContent(3), nw);
      return -1;
   }
   if (TMath::Abs(stat->GetBinContent(4) - nw) > .1) {
      printf("\n >>> Test failure: test input object not correctly received on all workers (%.0f/%d)\n",
             stat->GetBinContent(4), nw);
      return -1;
   }

   // Fill times
   PT_GetLastTimes(tt);

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Testing passing arguments to packages

Int_t PT_PackageArguments(void *, RunTimes &tt)
{
   // Checking arguments
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }
   PutPoint();

   // Passing a 'const char *': upload packtest1
   const char *pack1 = "packtest1";
   PutPoint();
   gProof->UploadPackage(gPack1);
   // Check the result
   TList *packs = gProof->GetListOfPackages();
   if (!packs || !packs->FindObject(pack1)) {
      printf("\n >>> Test failure: could not upload '%s'!\n", gPack1.Data());
      return -1;
   }
   // Enable the package now passing a 'const char *' argument
   TString arg("ProofTest.ConstChar");
   if (gProof->EnablePackage(pack1, arg) != 0) {
      printf("\n >>> Test failure: could not enable '%s' with argument: '%s'!\n", gPack1.Data(), arg.Data());
      return -1;
   }

   // Type of test
   gProof->SetParameter("ProofTests_Type", "PackTest1");

   // Define the number of events
   Long64_t nevt = 1;

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Variable to check
   gProof->SetParameter("testenv", arg.Data());

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gProof->Process(gTestsSel.Data(), nevt);
   }

   // Some cleanup
   gProof->ClearPackage(pack1);
   gProof->DeleteParameters("ProofTests_Type");
   gProof->DeleteParameters("testenv");

   // Make sure the query result is there
   PutPoint();
   TQueryResult *qr = 0;
   if (!(qr = gProof->GetQueryResult())) {
      printf("\n >>> Test failure: query result not found\n");
      return -1;
   }

   // Make sure the output list is there
   PutPoint();
   if (!(gProof->GetOutputList())) {
      printf("\n >>> Test failure: output list not found\n");
      return -1;
   }

   // Check the 'histo's
   PutPoint();
   TH1I *stat = dynamic_cast<TH1I*>(gProof->GetOutputList()->FindObject("TestStat"));
   if (!stat) {
      printf("\n >>> Test failure: 'TestStat' histo not found\n");
      return -1;
   }

   // Test how many workers got everything successfully
   Int_t nw = (Int_t) stat->GetBinContent(1);
   PutPoint();

   if (TMath::Abs(stat->GetBinContent(2) - nw) > .1) {
      printf("\n >>> Test failure: var '%s' not correctly set on all workers (%.0f/%d)\n",
             arg.Data(), stat->GetBinContent(2), nw);
      return -1;
   }

   // Passing a 'TList *': upload packtest2
   const char *pack2 = "packtest2";
   PutPoint();
   gProof->UploadPackage(gPack2);
   // Check the result
   packs = gProof->GetListOfPackages();
   if (!packs || !packs->FindObject(pack2)) {
      printf("\n >>> Test failure: could not upload '%s'!\n", gPack2.Data());
      return -1;
   }
   // Testing recursive enabling via dependencies: upload packtest3
   const char *pack3 = "packtest3";
   PutPoint();
   gProof->UploadPackage(gPack3);
   // Check the result
   packs = gProof->GetListOfPackages();
   if (!packs || !packs->FindObject(pack3)) {
      printf("\n >>> Test failure: could not upload '%s'!\n", gPack3.Data());
      return -1;
   }
   // Create the argument list
   TList *argls = new TList;
   argls->Add(new TNamed("ProofTest.ArgOne", "2."));
   argls->Add(new TNamed("ProofTest.ArgTwo", "3."));
   argls->Add(new TNamed("ProofTest.ArgThree", "4."));
   // Enable the package now passing the 'TList *' argument
   if (gProof->EnablePackage(pack3, argls) != 0) {
      printf("\n >>> Test failure: could not enable '%s' with argument: '%s'!\n", gPack3.Data(), arg.Data());
      return -1;
   }
   // Check list of enabled packages
   TList *enpkg = gProof->GetListOfEnabledPackages();
   if (!enpkg || enpkg->GetSize() < 3) {
      printf("\n >>> Test failure: not all requested packages enabled\n");
      if (!enpkg->FindObject("packtest1")) printf("\n >>> Test failure: 'packtest1' not enabled\n");
      if (!enpkg->FindObject("packtest2")) printf("\n >>> Test failure: 'packtest2' not enabled\n");
      if (!enpkg->FindObject("packtest3")) printf("\n >>> Test failure: 'packtest3' not enabled\n");
      return -1;
   }

   // Type of test
   gProof->SetParameter("ProofTests_Type", "PackTest2");

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Variable to check
   TString envs("ProofTest.ArgOne,ProofTest.ArgTwo,ProofTest.ArgThree");
   gProof->SetParameter("testenv", envs.Data());

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gProof->Process(gTestsSel.Data(), nevt);
   }

   // Some cleanup
   gProof->ClearPackage(pack2);
   gProof->DeleteParameters("ProofTests_Type");
   gProof->DeleteParameters("testenv");

   // Make sure the query result is there
   PutPoint();
   qr = 0;
   if (!(qr = gProof->GetQueryResult())) {
      printf("\n >>> Test failure: query result not found\n");
      return -1;
   }

   // Make sure the output list is there
   PutPoint();
   if (!(gProof->GetOutputList())) {
      printf("\n >>> Test failure: output list not found\n");
      return -1;
   }

   // Check the 'histo's
   PutPoint();
   if (!(stat = dynamic_cast<TH1I*>(gProof->GetOutputList()->FindObject("TestStat")))) {
      printf("\n >>> Test failure: 'TestStat' histo not found\n");
      return -1;
   }

   // Test how many workers got everything successfully
   nw = (Int_t) stat->GetBinContent(1);
   PutPoint();

   if (TMath::Abs(stat->GetBinContent(2) - nw) > .1) {
      printf("\n >>> Test failure: var 'ProofTest.ArgOne' not correctly set on all workers (%.0f/%d)\n",
             stat->GetBinContent(2), nw);
      return -1;
   }

   if (TMath::Abs(stat->GetBinContent(3) - nw) > .1) {
      printf("\n >>> Test failure: var 'ProofTest.ArgTwo' not correctly set on all workers (%.0f/%d)\n",
             stat->GetBinContent(3), nw);
      return -1;
   }

   if (TMath::Abs(stat->GetBinContent(4) - nw) > .1) {
      printf("\n >>> Test failure: var 'ProofTest.ArgThree' not correctly set on all workers (%.0f/%d)\n",
             stat->GetBinContent(4), nw);
      return -1;
   }

   // Fill times
   PT_GetLastTimes(tt);

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Test run for the H1 and Simple analysis in asynchronous mode

Int_t PT_H1SimpleAsync(void *arg, RunTimes &tt)
{
   // Checking arguments
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }
   // Not yet supported for PROOF-Lite
   if (gProof->IsLite()) {
      return 1;
   }
   // Not supported in dynamic mode
   if (gDynamicStartup) {
      return 1;
   }
   PutPoint();

   // Set/unset the parallel unzip flag
   AssertParallelUnzip();

   // Are we asked to change the packetizer strategy?
   if (arg) {
      PT_Packetizer_t *strategy = (PT_Packetizer_t *)arg;
      if (strcmp(strategy->fName, "TPacketizer")) {
         gProof->SetParameter("PROOF_Packetizer", strategy->fName);
      } else {
         if (strategy->fType != 1)
            gProof->SetParameter("PROOF_PacketizerStrategy", strategy->fType);
      }
   }

   // Create the file collection
   PutPoint();
   TFileCollection *fc = new TFileCollection("h42");

   // Assert the files, if needed
   if (!gh1ok) {
      if (PT_H1AssertFiles(gh1src.Data()) != 0) {
         printf("\n >>> Test failure: could not assert the H1 files\n");
         return -1;
      }
   }
   Int_t i = 0;
   for (i = 0; i < 4; i++) {
      fc->Add(new TFileInfo(TString::Format("%s%c%s", gh1src.Data(), gh1sep, gh1file[i])));
   }

   // Clear associated memory cache
   if (gClearCache && PT_H1ReleaseCache(gh1src.Data()) != 0) {
      gProof->SetPrintProgress(0);
      printf("\n >>> Test failure: could not clear memory cache for the H1 files\n");
      return -1;
   }

   // Clear the list of query results
   PutPoint();
   if (gProof->GetQueryResults()) {
      gProof->GetQueryResults()->Clear();
      gProof->Remove("cleanupdir");
   }

   // Define the number of events and histos
   Long64_t nevt = 1000000;
   Int_t nhist = 16;
   // Submit the processing requests
   PutPoint();
   {  SwitchProgressGuard spg;
      gProof->Process(fc, gH1Sel.Data(), "ASYN");

      // The number of histograms is added as parameter in the input list
      gProof->SetParameter("ProofSimple_NHist", (Long_t)nhist);
      gProof->Process(gSimpleSel.Data(), nevt, "ASYN");

      // Wait a bit as a function of previous runnings
      Double_t dtw = 10;
      if (gH1Cnt > 0 && gSimpleTime > 0) {
         dtw = gH1Time / gH1Cnt + gSimpleTime / gSimpleCnt + 1;
         if (dtw < 10) dtw = 10;
      }
      Int_t tw = (Int_t) (5 * dtw);

      gTimedOut = kFALSE;
      TTimeOutTimer t(tw*1000);

      // Wait for the processing
      while (!gProof->IsIdle() && !gTimedOut)
         gSystem->InnerLoop();

   }
   PutPoint();

   // Restore settings
   gProof->DeleteParameters("PROOF_Packetizer");
   gProof->DeleteParameters("PROOF_PacketizerStrategy");

   // Retrieve the list of available query results
   TList *ql = gProof->GetQueryResults();
   if (ql && ql->GetSize() > 0) {
      ql->Print();
   }

   TString ref;
   TIter nxq(ql, kIterBackward);
   Int_t nd = 2;
   TQueryResult *qr = 0;
   while ((qr = (TQueryResult *)nxq()) && nd > 0) {
      ref.Form("%s:%s", qr->GetTitle(), qr->GetName());
      gProof->Retrieve(ref);
      qr = gProof->GetQueryResult(ref);
      if (qr && qr->GetSelecImp()) {
         if (!strcmp(qr->GetSelecImp()->GetTitle(), "h1analysis")) {
            PutPoint();
            if (PT_CheckH1(qr) != 0) return -1;
            nd--;
         } else if (!strcmp(qr->GetSelecImp()->GetTitle(), "ProofSimple")) {
            PutPoint();
            if (PT_CheckSimple(qr, nevt, nhist) != 0) return -1;
            nd--;
         } else {
            printf("\n >>> Test failure: query with unexpected selector '%s'\n", qr->GetSelecImp()->GetTitle());
            return -1;
         }
      } else {
         printf("\n >>> Test failure: query undefined (%p) or with empty selector macro ('%s')\n", qr, ref.Data());
         return -1;
      }
   }

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Done
   PutPoint();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Test run for the admin functionality

Int_t PT_AdminFunc(void *, RunTimes &tt)
{
   // Checking arguments
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }
   // Not yet supported for PROOF-Lite
   if (gProof->IsLite()) {
      return 1;
   }
   PutPoint();
   // Attach to the manager
   TProofMgr *mgr = gProof->GetManager();
   if (!mgr) {
      printf("\n >>> Test failure: no PROOF manager found\n");
      return -1;
   }
   PutPoint();
   // Directory for this test
   TString testDir(TString::Format("%s/stressProof-Admin", gtutdir.Data()));
   if (gSystem->AccessPathName(testDir)) {
      // Create the directory
      if (gSystem->MakeDirectory(testDir)) {
         printf("\n >>> Test failure: cannot create %s\n", testDir.Data());
         return -1;
      }
   }
   // Create a small test file
   TMacro testMacro;
   testMacro.AddLine("// Test macro");
   testMacro.AddLine("#include \"TSystem.h\"");
   testMacro.AddLine("void testMacro(Int_t opt = 1)");
   testMacro.AddLine("{");
   testMacro.AddLine("   if (opt == 1) {");
   testMacro.AddLine("      Printf(\"Pid: \", gSystem->GetPid());");
   testMacro.AddLine("   }");
   testMacro.AddLine("}");
   // Save the file in the temporary area
   TString testFile(TString::Format("%s/testMacro.C", testDir.Data()));
   if (!gSystem->AccessPathName(testFile)) {
      // The file exists: remove it first
      if (gSystem->Unlink(testFile)) {
         printf("\n >>> Test failure: cannot unlink %s\n", testFile.Data());
         return -1;
      }
   }
   testMacro.SaveSource(testFile);
   FileStat_t stloc;
   if (gSystem->GetPathInfo(testFile, stloc) != 0) {
      // The file was not created
      printf("\n >>> Test failure: file %s was not created\n", testFile.Data());
      return -1;
   }
   // Reference checksum
   std::unique_ptr<TMD5> testMacroMd5(testMacro.Checksum());
   if (!testMacroMd5.get()) {
      // MD5 sum not calculated
      printf("\n >>> Test failure: could not calculate the md5 sum of the test macro\n");
      return -1;
   }
   PutPoint();

   // Send the file to the sandbox
   if (mgr->PutFile(testFile, "~/", "force") != 0) {
      // The file was not sent over correctly
      printf("\n >>> Test failure: problems sending file to master sandbox\n");
      return -1;
   }
   PutPoint();
   // Retrieve the file
   TString testFile1(TString::Format("%s/testMacro.cxx", testDir.Data()));
   if (mgr->GetFile("~/testMacro.C", testFile1, "force") != 0) {
      // The file was not retrieved correctly
      printf("\n >>> Test failure: problems retrieving file from the master sandbox\n");
      return -1;
   }
   PutPoint();

   // Test 'ls'
   gSystem->RedirectOutput(0, 0, &gRH);
   TString testLs(TString::Format("%s/testLs.log", testDir.Data()));
   gSystem->RedirectOutput(testLs, "w", &gRHAdmin);
   mgr->Ls("~/testMacro.C");
   gSystem->RedirectOutput(0, 0, &gRHAdmin);
   gSystem->RedirectOutput(glogfile, "a", &gRH);
   TMacro macroLs(testLs);
   TString testLsLine = TString::Format("%s/testMacro.C", gsandbox.Data());
   // The first part of <tmp> maybe sligthly different
#if defined(R__MACOSX)
   if (testLsLine.Index(".proof") != kNPOS)
      testLsLine.Remove(0, testLsLine.Index(".proof"));
#else
   if (testLsLine.Index(".proof-tutorial") != kNPOS)
      testLsLine.Remove(0, testLsLine.Index(".proof-tutorial"));
#endif
   if (!macroLs.GetLineWith(testLsLine)) {
      printf("\n >>> Test failure: Ls: output not consistent (line: '%s')\n", testLsLine.Data());
      printf(" >>> Log file: '%s'\n", testLs.Data());
      printf("+++ BOF +++\n");
      macroLs.Print();
      printf("+++ EOF +++\n");
      return -1;
   }
   PutPoint();

   // Test 'more'
   gSystem->RedirectOutput(0, 0, &gRH);
   TString testMore(TString::Format("%s/testMore.log", testDir.Data()));
   gSystem->RedirectOutput(testMore, "w", &gRHAdmin);
   mgr->More("~/testMacro.C");
   gSystem->RedirectOutput(0, 0, &gRHAdmin);
   gSystem->RedirectOutput(glogfile, "a", &gRH);
   TMacro macroMore(testMore);
   if (macroMore.GetListOfLines()->GetSize() < 2) {
      printf("\n >>> Test failure: More output not too short: %d lines\n",
                                   macroMore.GetListOfLines()->GetSize());
      return -1;
   }
   TObjString *os = (TObjString *) macroMore.GetListOfLines()->First();
   while (!os->GetString().BeginsWith("// Test macro")) {
      macroMore.GetListOfLines()->Remove(macroMore.GetListOfLines()->First());
      os = (TObjString *) macroMore.GetListOfLines()->First();
   }
   std::unique_ptr<TMD5> testMoreMd5(macroMore.Checksum());
   if (!testMoreMd5.get()) {
      // MD5 sum not calculated
      printf("\n >>> Test failure: could not calculate the md5 sum of the 'more' result\n");
      return -1;
   }
   if (strcmp(testMacroMd5->AsString(), testMoreMd5->AsString())) {
      printf("\n >>> Test failure: More: result not consistent with reference: {%s, %s}\n",
                                         testMacroMd5->AsString(), testMoreMd5->AsString());
      return -1;
   }
   PutPoint();

   // Test 'stat'
   FileStat_t strem;
   if (mgr->Stat("~/testMacro.C", strem) != 0) {
      // Stat failure
      printf("\n >>> Test failure: could not stat remotely the test file\n");
      return -1;
   }
   if (strem.fSize != stloc.fSize) {
      // Stat failure
      printf("\n >>> Test failure: stat sizes inconsistent: %lld vs %lld (bytes)\n", strem.fSize, stloc.fSize);
      return -1;
   }
   PutPoint();

   // Test 'cp' and 'md5sum;
   if (mgr->Cp("http://root.cern.ch/files/h1/dstarmb.root", "~/") != 0) {
      // Cp failure
      printf("\n >>> Test failure: could not retrieve remotely dstarmb.root from the Root Web server\n");
      return -1;
   }
   TString sum;
   if (mgr->Md5sum("~/dstarmb.root", sum) != 0) {
      // MD5
      printf("\n >>> Test failure: calculating the md5sum of dstarmb.root\n");
      return -1;
   }
   if (sum != "0a60055370e16d954f90fb50c2d1a801") {
      // MD5 wrong
      printf("\n >>> Test failure: wrong value for md5sum of dstarmb.root: %s\n", sum.Data());
      return -1;
   }
   PutPoint();

   // Fill times
   PT_GetLastTimes(tt);

   // Done
   PutPoint();
   return 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Test processing of sub-samples (entries-from-first) from files with the
/// 'event' structures

Int_t PT_EventRange(void *arg, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   // Set/unset the parallel unzip flag
   AssertParallelUnzip();

   // Are we asked to change the packetizer strategy?
   const char *pack = "TPacketizer";
   if (arg) {
      PT_Packetizer_t *strategy = (PT_Packetizer_t *)arg;
      if (strcmp(strategy->fName, "TPacketizer")) {
         gProof->SetParameter("PROOF_Packetizer", strategy->fName);
         pack = strategy->fName;
      } else {
         if (strategy->fType != 1)
            gProof->SetParameter("PROOF_PacketizerStrategy", strategy->fType);
      }
   }

   // Test first with a chain
   PutPoint();
   TChain *chain = new TChain("EventTree");

   // Assert the files, if needed
   if (!geventok) {
      if (PT_EventAssertFiles(geventsrc.Data(), geventnf) != 0) {
         gProof->SetPrintProgress(0);
         printf("\n >>> Test failure: could not assert the event files\n");
         return -1;
      }
   }
   Int_t i = 0;
   for (i = 0; i < geventnf; i++) {
      chain->AddFile(TString::Format("%s/event_%d.root", geventsrc.Data(), i+1));
   }

   // Clear associated memory cache
   if (gClearCache && PT_EventReleaseCache(geventsrc.Data(), geventnf) != 0) {
      gProof->SetPrintProgress(0);
      printf("\n >>> Test failure: could not clear memory cache for the event files\n");
      return -1;
   }

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Load special class for event ranges checks
   if (gProof->Load(TString::Format("%s,%s", gProcFileElem.Data(), gEmptyInclude.Data())) != 0) {
      gProof->SetPrintProgress(0);
      printf("\n >>> Test failure: could not load auxilliary files %s and %s\n",
             gProcFileElem.Data(), gEmptyInclude.Data());
      return -1;
   }

   // Add some parameters for later checking in Terminate
   Int_t ifst = gEventFst / gEventSiz + 1;
   Long64_t efst = gEventFst - (ifst - 1) * gEventSiz;
   TString ffst = TString::Format("%s/event_%d.root?fst=%lld", geventsrc.Data(), ifst, efst);
   gProof->SetParameter("Range_First_File", ffst.Data());

   Int_t ilst = (gEventFst + gEventNum) / gEventSiz + 1;
   Long64_t elst = (gEventFst + gEventNum) - (ilst - 1) * gEventSiz - 1;
   TString flst = TString::Format("%s/event_%d.root?lst=%lld", geventsrc.Data(), ilst, elst);
   gProof->SetParameter("Range_Last_File", flst.Data());
   gProof->SetParameter("Range_Num_Files", (Int_t) (ilst - ifst + 1));

   // Process
   PutPoint();
   chain->SetProof();
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      chain->Process(gEventProcSel.Data(), "", gEventNum, gEventFst);
      gTimer.Stop();
   }
   gProof->RemoveChain(chain);

   // Count
   gEventCnt++;
   gEventTime += gTimer.RealTime();

   // Check the result
   Int_t rcch = PT_CheckEvent(gProof->GetQueryResult(), pack);

   // We are done if not dataset test possible
   if (gSkipDataSetTest) {
      // Restore settings
      gProof->DeleteParameters("PROOF_Packetizer");
      gProof->DeleteParameters("PROOF_PacketizerStrategy");
      return rcch;
   }

   // Create the dataset
   TFileCollection *fc = new TFileCollection("dsevent", "", "");
   for (i = 0; i < geventnf; i++) {
      fc->Add(new TFileInfo(TString::Format("%s/event_%d.root", geventsrc.Data(), i+1)));
   }

   // Register the dataset
   PutPoint();
   gProof->RegisterDataSet("dsevent", fc);

   // Check the result
   if (!gProof->ExistsDataSet("dsevent")) {
      printf("\n >>> Test failure: could not register 'dsevent'\n");
      return -1;
   }

   // Verify the dataset
   PutPoint();
   if (gProof->VerifyDataSet("dsevent") < 0) {
      printf("\n >>> Test failure: could not verify 'dsevent'!\n");
      return -1;
   }

   // Clear associated memory cache
   if (gClearCache && PT_EventReleaseCache(geventsrc.Data(), geventnf) != 0) {
      gProof->SetPrintProgress(0);
      printf("\n >>> Test failure: could not clear memory cache for the event files\n");
      return -1;
   }

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      gProof->Process("dsevent", gEventProcSel.Data(), "", gEventNum, gEventFst);
      gTimer.Stop();
   }
   gProof->RemoveDataSet("dsevent");

   // Restore settings
   gProof->DeleteParameters("PROOF_Packetizer");
   gProof->DeleteParameters("PROOF_PacketizerStrategy");

   // Count
   gEventCnt++;
   gEventTime += gTimer.RealTime();

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckEvent(gProof->GetQueryResult(), pack);
}

////////////////////////////////////////////////////////////////////////////////
/// Test TProofOutputFile technology to create a ntuple, with or without
/// submergers

Int_t PT_POFNtuple(void *opts, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   PT_Option_t *ptopt = (PT_Option_t *) opts;

   // Setup submergers if required
   if (ptopt && ptopt->fTwo > 0) {
      gProof->SetParameter("PROOF_UseMergers", 0);
   }

   // Output file
   TString fout("<datadir>/ProofNtuple.root");
   gProof->AddInput(new TNamed("PROOF_OUTPUTFILE", fout.Data()));

   // We use the 'NtpRndm' for a fixed values of randoms; we need to send over the file
   gProof->SetInputDataFile(gNtpRndm);
   // Set the related parameter
   gProof->SetParameter("PROOF_USE_NTP_RNDM","yes");

   // Define the number of events and histos
   Long64_t nevt = 1000;

   // We do not plot the ntuple (we are in batch mode)
   gProof->SetParameter("PROOF_NTUPLE_DONT_PLOT", "yes");

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      gProof->Process(gNtupleSel.Data(), nevt);
      gTimer.Stop();
   }

   // Remove any setting related to submergers
   gProof->DeleteParameters("PROOF_UseMergers");
   gProof->DeleteParameters("PROOF_NTUPLE_DONT_PLOT");
   gProof->DeleteParameters("PROOF_USE_NTP_RNDM");
   gProof->SetInputDataFile(0);

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckNtuple(gProof->GetQueryResult(), nevt);
}

////////////////////////////////////////////////////////////////////////////////
/// Test TProofOutputFile technology to create a dataset

Int_t PT_POFDataset(void *, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   const char *dsname = "testNtuple";
   // Clean-up any existing dataset with that name
   if (gProof->ExistsDataSet(dsname)) gProof->RemoveDataSet(dsname);

   // Ask for registration of the dataset (the default is the the TFileCollection is return
   // without registration; the name of the TFileCollection is the name of the dataset
   gProof->SetParameter("SimpleNtuple.root", dsname);

   // We use the 'NtpRndm' for a fixed values of randoms; we need to send over the file
   gProof->SetInputDataFile(gNtpRndm);
   // Set the related parameter
   gProof->SetParameter("PROOF_USE_NTP_RNDM","yes");

   // Define the number of events and histos
   Long64_t nevt = 1000000;

   // We do not plot the ntuple (we are in batch mode)
   gProof->SetParameter("PROOF_NTUPLE_DONT_PLOT", "yes");

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      gProof->Process(gNtupleSel.Data(), nevt);
      gTimer.Stop();
   }

   // Remove any setting related to submergers
   gProof->DeleteParameters("PROOF_NTUPLE_DONT_PLOT");
   gProof->DeleteParameters("SimpleNtuple.root");
   gProof->DeleteParameters("PROOF_USE_NTP_RNDM");
   gProof->SetInputDataFile(0);

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckDataset(gProof->GetQueryResult(), nevt);
}

////////////////////////////////////////////////////////////////////////////////
/// Test processing of multiple trees in the same files

Int_t PT_MultiTrees(void *, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   const char *dsname = "testNtuple";
   // There must be a dataset 'testNtuple' already registered and validated
   if (!gProof->ExistsDataSet(dsname)) {
      printf("\n >>> Test failure: dataset '%s' does not exist\n", dsname);
      return -1;
   }

   // Get the associated TFileCollection
   TFileCollection *fc = gProof->GetDataSet(dsname);
   if (!fc) {
      printf("\n >>> Test failure: unable to get TFileCollection for dataset '%s'\n", dsname);
      return -1;
   }

   // Now create a TDSet out of the TFileCollection
   TDSet *dset = new TDSet("testntps", "ntuple", "/", "TTree");
   TChain *ch1 = new TChain("ntuple");
   TChain *ch2 = new TChain("ntuple2");
   TIter nxf(fc->GetList());
   TFileInfo *fi = 0;
   while ((fi = (TFileInfo *) nxf())) {
      dset->Add(fi->GetCurrentUrl()->GetUrl());
      ch1->Add(fi->GetCurrentUrl()->GetUrl());
      ch2->Add(fi->GetCurrentUrl()->GetUrl());
   }

   // Check the ntuple content by filling some histos
   TH1F *h1s[2] = {0};
   h1s[0] = new TH1F("h1_1", "3*px+2 with px**2+py**2>1", 50, -15., 15.);
   h1s[1] = new TH1F("h1_2", "vx**2+vy**2 with abs(vz)<.1", 50, 0., 10.);

   Int_t rch1s = 0;
   TString emsg;
   const char *type[3] = { "dsname", "TDSet", "TChain" };
   for (Int_t j = 0; j < 3; j++) {

      PutPoint();

      if (j == 0) {
         // Fill the first histo from the first ntuple
         gProof->SetDataSetTreeName(dsname, "ntuple");
         {  SwitchProgressGuard spg;
            gProof->DrawSelect(dsname, "3*px+2>>h1_1", "px*px+py*py>1");
         }
         // Fill the second histo from the second ntuple
         gProof->SetDataSetTreeName(dsname, "ntuple2");
         {  SwitchProgressGuard spg;
            gProof->DrawSelect(dsname, "vx*vx+vy*vy>>h1_2", "vz>-0.1&&vz<0.1");
         }
      } else if (j == 1) {
         // Fill the first histo from the first ntuple
         {  SwitchProgressGuard spg;
            gProof->DrawSelect(dset, "3*px+2>>h1_1", "px*px+py*py>1");
         }
         // Fill the second histo from the second ntuple
         dset->SetObjName("ntuple2");
         {  SwitchProgressGuard spg;
            gProof->DrawSelect(dset, "vx*vx+vy*vy>>h1_2", "vz>-0.1&&vz<0.1");
         }
      } else {
         // Fill the first histo from the first ntuple
         {  SwitchProgressGuard spg;
            ch1->Draw("3*px+2>>h1_1", "px*px+py*py>1");
         }
         // Fill the second histo from the second ntuple
         {  SwitchProgressGuard spg;
            ch2->Draw("vx*vx+vy*vy>>h1_2", "vz>-0.1&&vz<0.1");
         }
      }

      rch1s = 0;
      // Check the histogram entries and mean values
      Int_t hent[2] = { 607700, 96100};
      Double_t hmea[2] = { 2.022, 1.859};
      for (Int_t i = 0; i < 2; i++) {
         if ((Int_t)(h1s[i]->GetEntries()) != hent[i]) {
            emsg.Form("%s: '%s' histo: wrong number of entries (%d: expected %d)",
                      type[j], h1s[i]->GetName(), (Int_t)(h1s[i]->GetEntries()), hent[i]);
            rch1s = -1;
            break;
         }
         if (TMath::Abs((h1s[i]->GetMean() - hmea[i]) / hmea[i]) > 0.001) {
            emsg.Form("%s: '%s' histo: wrong mean (%f: expected %f)",
                      type[j], h1s[i]->GetName(), h1s[i]->GetMean(), hmea[i]);
            rch1s = -1;
            break;
         }
      }
   }

   // Cleanup
   for (Int_t i = 0; i < 2; i++) delete h1s[i];

   // Check the result
   if (rch1s != 0) {
      printf("\n >>> Test failure: %s\n", emsg.Data());
      return -1;
   }

   // Clean-up
   gProof->RemoveDataSet(dsname);

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return rch1s;
}

////////////////////////////////////////////////////////////////////////////////
/// Test processing of TTree friends in PROOF

Int_t PT_Friends(void *sf, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   // Not supported in dynamic mode
   if (gDynamicStartup) {
      return 1;
   }
   PutPoint();

   // Separate or same file ?
   Bool_t sameFile = (sf) ? kTRUE : kFALSE;

   // File generation: we use TPacketizerFile in here to create two files per node
   TList *wrks = gProof->GetListOfSlaveInfos();
   if (!wrks) {
      printf("\n >>> Test failure: could not get the list of information about the workers\n");
      return -1;
   }

   // Create the map
   TString fntree;
   TMap *files = new TMap;
   files->SetName("PROOF_FilesToProcess");
   TIter nxwi(wrks);
   TSlaveInfo *wi = 0;
   while ((wi = (TSlaveInfo *) nxwi())) {
      fntree.Form("tree_%s.root", wi->GetOrdinal());
      THashList *wrklist = (THashList *) files->GetValue(wi->GetName());
      if (!wrklist) {
         wrklist = new THashList;
         wrklist->SetName(wi->GetName());
         files->Add(new TObjString(wi->GetName()), wrklist);
      }
      wrklist->Add(new TObjString(fntree));
   }
   Int_t nwrk = wrks->GetSize();

   // Generate the files
   gProof->AddInput(files);
   if (sameFile) {
      Printf("runProof: friend tree stored in the same file as the main tree");
      gProof->SetParameter("ProofAux_Action", "GenerateTreesSameFile");
   } else {
      gProof->SetParameter("ProofAux_Action", "GenerateTrees");
   }

   // File generation: define the number of events per worker
   Long64_t nevt = 1000;
   gProof->SetParameter("ProofAux_NEvents", (Long64_t)nevt);
   // Special Packetizer
   gProof->SetParameter("PROOF_Packetizer", "TPacketizerFile");
   // Now process
   gProof->Process(gAuxSel.Data(), 1);
   // Remove the packetizer specifications
   gProof->DeleteParameters("PROOF_Packetizer");

   // Check that we got some output
   if (!gProof->GetOutputList()) {
      printf("\n >>> Test failure: output list not found!\n");
      return -1;
   }

   // Create the TDSet objects
   TDSet *dset = new TDSet("Tmain", "Tmain");
   TDSet *dsetf = new TDSet("Tfrnd", "Tfrnd");
   // Fill them with the information found in the output list
   Bool_t foundMain = kFALSE, foundFriend = kFALSE;
   TIter nxo(gProof->GetOutputList());
   TObject *o = 0;
   TObjString *os = 0;
   while ((o = nxo())) {
      TList *l = dynamic_cast<TList *> (o);
      if (l && !strncmp(l->GetName(), "MainList-", 9)) {
         foundMain = kTRUE;
         TIter nxf(l);
         while ((os = (TObjString *) nxf()))
            dset->Add(os->GetName());
      }
   }
   nxo.Reset();
   while ((o = nxo())) {
      TList *l = dynamic_cast<TList *> (o);
      if (l && !strncmp(l->GetName(), "FriendList-", 11)) {
         foundFriend = kTRUE;
         TIter nxf(l);
         while ((os = (TObjString *) nxf()))
            dsetf->Add(os->GetName());
      }
   }

   // If we did not found the main or the friend meta info we fail
   if (!foundMain || !foundFriend) {
      printf("\n >>> Test failure: 'main' or 'friend' meta info missing!\n");
      return -1;
   }

   // Connect the two datasets for processing
   dset->AddFriend(dsetf, "friend");

   // We do not plot the ntuple (we are in batch mode)
   gProof->SetParameter("PROOF_DONT_PLOT", "yes");

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      dset->Process(gFriendsSel.Data());
      gTimer.Stop();
   }

   // Remove any setting
   gProof->DeleteParameters("PROOF_DONT_PLOT");
   gProof->GetInputList()->Remove(files);
   files->SetOwner(kTRUE);
   SafeDelete(files);
   // Clear the files created by this run
   gProof->ClearData(TProof::kUnregistered | TProof::kForceClear);

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckFriends(gProof->GetQueryResult(), nevt * nwrk, 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Test processing of TTree in subdirectories

Int_t PT_TreeSubDirs(void*, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   // Not supported in dynamic mode
   if (gDynamicStartup) {
      return 1;
   }
   PutPoint();

   // File generation: we use TPacketizerFile in here to create two files per node
   TList *wrks = gProof->GetListOfSlaveInfos();
   if (!wrks) {
      printf("\n >>> Test failure: could not get the list of information about the workers\n");
      return -1;
   }

   // Create the map
   TString fntree;
   TMap *files = new TMap;
   files->SetName("PROOF_FilesToProcess");
   TIter nxwi(wrks);
   TSlaveInfo *wi = 0;
   while ((wi = (TSlaveInfo *) nxwi())) {
      fntree.Form("tree_%s.root", wi->GetOrdinal());
      THashList *wrklist = (THashList *) files->GetValue(wi->GetName());
      if (!wrklist) {
         wrklist = new THashList;
         wrklist->SetName(wi->GetName());
         files->Add(new TObjString(wi->GetName()), wrklist);
      }
      wrklist->Add(new TObjString(fntree));
   }
   Int_t nwrk = wrks->GetSize();

   // Generate the files
   gProof->AddInput(files);
   gProof->SetParameter("ProofAux_Action", "GenerateTrees:dir1/dir2/dir3");

   // File generation: define the number of events per worker
   Long64_t nevt = 1000;
   gProof->SetParameter("ProofAux_NEvents", (Long64_t)nevt);
   // Special Packetizer
   gProof->SetParameter("PROOF_Packetizer", "TPacketizerFile");
   // Now process
   gProof->Process(gAuxSel.Data(), 1);
   // Remove the packetizer specifications
   gProof->DeleteParameters("PROOF_Packetizer");

   // Check that we got some output
   if (!gProof->GetOutputList()) {
      printf("\n >>> Test failure: output list not found!\n");
      return -1;
   }

   // Create the TChain objects
   TChain *dset = new TChain("dir1/dir2/dir3/Tmain");
   // Fill them with the information found in the output list
   Bool_t foundMain = kFALSE;
   TIter nxo(gProof->GetOutputList());
   TObject *o = 0;
   TObjString *os = 0;
   while ((o = nxo())) {
      TList *l = dynamic_cast<TList *> (o);
      if (l && !strncmp(l->GetName(), "MainList-", 9)) {
         foundMain = kTRUE;
         TIter nxf(l);
         while ((os = (TObjString *) nxf()))
            dset->Add(os->GetName());
      }
   }
   dset->SetProof();

   // If we did not found the main or the friend meta info we fail
   if (!foundMain) {
      printf("\n >>> Test failure: 'main' meta info missing!\n");
      return -1;
   }

   // We do not plot the ntuple (we are in batch mode)
   gProof->SetParameter("PROOF_DONT_PLOT", "yes");

   // We do use friends
   gProof->SetParameter("PROOF_NO_FRIENDS", "yes");

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      dset->Process(gFriendsSel.Data());
      gTimer.Stop();
   }

   // Remove any setting
   gProof->DeleteParameters("PROOF_DONT_PLOT");
   gProof->GetInputList()->Remove(files);
   files->SetOwner(kTRUE);
   SafeDelete(files);
   // Clear the files created by this run
   gProof->ClearData(TProof::kUnregistered | TProof::kForceClear);

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckFriends(gProof->GetQueryResult(), nevt * nwrk, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Test run for the ProofSimple analysis (see tutorials) passing the
/// selector by object

Int_t PT_SimpleByObj(void *submergers, RunTimes &tt)
{
   // Checking arguments
   PutPoint();
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }

   // Setup submergers if required
   if (submergers) {
      gProof->SetParameter("PROOF_UseMergers", 0);
   }

   // Define the number of events and histos
   Long64_t nevt = 1000000;
   Int_t nhist = 16;
   // The number of histograms is set inside the selector object; make sure
   // that it is not passed in the input list
   gProof->DeleteParameters("ProofSimple_NHist");

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Define TSelector object. We use  reflection to avoid including the header,
   // so being able to change the tutorial directory
   TString emsg;
   gProof->Load(gSimpleSel);
   TSelector *sel = TSelector::GetSelector(gSimpleSel);
   if (sel) {
      TClass *cl = sel->IsA();
      if (cl) {
         TDataMember *dm = cl->GetDataMember("fNhist");
         if (dm) {
            TMethodCall *setter = dm->SetterMethod(cl);
            if (setter) {
               setter->Execute(sel, TString::Format("%d", nhist).Data());
            } else {
               emsg = "no SetterMethod for fNhist: check version of ProofSimple";
            }
         } else {
            emsg = "fNhist not found";
         }
      } else {
         emsg = "IsA() failed";
      }
   } else {
      emsg = "GetSelector failed";
   }
   if (!emsg.IsNull()) {
      printf("\n >>> Test failure: initializing ProofSimple selector: %s\n", emsg.Data());
      return -1;
   }

   // Process
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      gProof->Process(sel, nevt);
      gTimer.Stop();
   }

   // Count
   gSimpleCnt++;
   gSimpleTime += gTimer.RealTime();

   // Remove any setting related to submergers
   gProof->DeleteParameters("PROOF_UseMergers");

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckSimple(gProof->GetQueryResult(), nevt, nhist);
}

////////////////////////////////////////////////////////////////////////////////
/// Test run for the H1 analysis as a chain reading the data from HTTP and
/// passing the selector by object

Int_t PT_H1ChainByObj(void *, RunTimes &tt)
{
   // Checking arguments
   if (!gProof) {
      printf("\n >>> Test failure: no PROOF session found\n");
      return -1;
   }
   // Not yet supported for PROOF-Lite
   if (gSkipDataSetTest) {
      return 1;
   }
   PutPoint();

   // Set/unset the parallel unzip flag
   AssertParallelUnzip();

   // Create the chain
   PutPoint();
   TChain *chain = new TChain("h42");

   // Assert the files, if needed
   if (!gh1ok) {
      if (PT_H1AssertFiles(gh1src.Data()) != 0) {
         gProof->SetPrintProgress(0);
         printf("\n >>> Test failure: could not assert the H1 files\n");
         return -1;
      }
   }
   Int_t i = 0;
   for (i = 0; i < 4; i++) {
      chain->Add(TString::Format("%s%c%s", gh1src.Data(), gh1sep, gh1file[i]));
   }

   // Clear associated memory cache
   if (gClearCache && PT_H1ReleaseCache(gh1src.Data()) != 0) {
      gProof->SetPrintProgress(0);
      printf("\n >>> Test failure: could not clear memory cache for the H1 files\n");
      return -1;
   }

   // Clear the list of query results
   if (gProof->GetQueryResults()) gProof->GetQueryResults()->Clear();

   // Load TSelector
   gProof->Load(gH1Sel);
   TSelector *sel = TSelector::GetSelector(gH1Sel);

   // Process
   PutPoint();
   chain->SetProof();
   PutPoint();
   {  SwitchProgressGuard spg;
      gTimer.Start();
      chain->Process(sel);
      gTimer.Stop();
   }
   gProof->RemoveChain(chain);

   // Count
   gH1Cnt++;
   gH1Time += gTimer.RealTime();

   // The runtimes
   PT_GetLastProofTimes(tt);

   // Check the results
   PutPoint();
   return PT_CheckH1(gProof->GetQueryResult());
}
