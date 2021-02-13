/**
  \file hadd.cxx
  \brief This program will add histograms (see note) and Trees from a list of root files and write them to a target root file.
  The target file is newly created and must not be
  identical to one of the source files.

  Syntax:
  ```{.cpp}
       hadd targetfile source1 source2 ...
  ```
  or
  ```{.cpp}
       hadd -f targetfile source1 source2 ...
  ```
  (targetfile is overwritten if it exists)

  \param -a   Append to the output
  \param -f   Force overwriting of output file.
  \param -f[0-9] Set target compression level. 0 = uncompressed, 6 = highly compressed.
  \param -fk  Sets the target file to contain the baskets with the same compression
              as the input files (unless -O is specified). Compresses the meta data
              using the compression level specified in the first input or the
              compression setting after fk (for example 206 when using -fk206)
  \param -ff  The compression level used is the one specified in the first input

  \param -k   Skip corrupt or non-existent files, do not exit
  \param -O   Re-optimize basket size when merging TTree
  \param -v   Explicitly set the verbosity level: 0 request no output, 99 is the default
  \param -j   Parallelise the execution in multiple processes
  \param -dbg  Parallelise the execution in multiple processes in debug mode (Does not delete  partial  files  stored
              inside working directory)
  \param -d   Carry out the partial multiprocess execution in the specified directory
  \param -n   Open at most `n` at once (use 0 to request to use the system maximum)
  \param -experimental-io-features `<feature>` Enables the corresponding experimental feature for output trees
  \return hadd returns a status code: 0 if OK, -1 otherwise

  When the -f option is specified, one can also specify the compression
  level of the target file. By default the compression level is 1 (kDefaultZLIB), but
  if "-f0" is specified, the target file will not be compressed.
  if "-f6" is specified, the compression level 6 will be used.

  For example assume 3 files f1, f2, f3 containing histograms hn and Trees Tn
   - f1 with h1 h2 h3 T1
   - f2 with h1 h4 T1 T2
   - f3 with h5
  the result of
  ```
    hadd -f x.root f1.root f2.root f3.root
  ```
  will be a file x.root with h1 h2 h3 h4 h5 T1 T2
  where
   - h1 will be the sum of the 2 histograms in f1 and f2
   - T1 will be the merge of the Trees in f1 and f2

  The files may contain sub-directories.

  If the source files contains histograms and Trees, one can skip
  the Trees with
  ```
       hadd -T targetfile source1 source2 ...
  ```

  Wildcarding and indirect files are also supported
  ```
      hadd result.root  myfil*.root
  ```
  will merge all files in myfil*.root
  ```
      hadd result.root file1.root @list.txt file2. root myfil*.root
  ```
  will merge file1.root, file2.root, all files in myfil*.root
  and all files in the indirect text file list.txt ("@" as the first
  character of the file indicates an indirect file. An indirect file
  is a text file containing a list of other files, including other
  indirect files, one line per file).

  If the sources and and target compression levels are identical (default),
  the program uses the TChain::Merge function with option "fast", ie
  the merge will be done without  unzipping or unstreaming the baskets
  (i.e. direct copy of the raw byte on disk). The "fast" mode is typically
  5 times faster than the mode unzipping and unstreaming the baskets.

  If the option -cachesize is used, hadd will resize (or disable if 0) the
  prefetching cache use to speed up I/O operations.

  For options that take a size as argument, a decimal number of bytes is expected.
  If the number ends with a `k`, `m`, `g`, etc., the number is multiplied
  by 1000 (1K), 1000000 (1MB), 1000000000 (1G), etc.
  If this prefix is followed by `i`, the number is multiplied by the traditional
  1024 (1KiB), 1048576 (1MiB), 1073741824 (1GiB), etc.
  The prefix can be optionally followed by B whose casing is ignored,
  eg. 1k, 1K, 1Kb and 1KB are the same.

  \note By default histograms are added. However hadd does not support the case where
         histograms have their bit TH1::kIsAverage set.

  \authors Rene Brun, Dirk Geppert, Sven A. Schmidt, Toby Burnett
*/
#include "Compression.h"
#include <ROOT/RConfig.hxx>
#include "ROOT/TIOFeatures.hxx"
#include "TFile.h"
#include "THashList.h"
#include "TKey.h"
#include "TClass.h"
#include "TSystem.h"
#include "TUUID.h"
#include "ROOT/StringConv.hxx"
#include "snprintf.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <climits>
#include <sstream>
#include "haddCommandLineOptionsHelp.h"

#include "TFileMerger.h"
#ifndef R__WIN32
#include "ROOT/TProcessExecutor.hxx"
#endif

////////////////////////////////////////////////////////////////////////////////

int main( int argc, char **argv )
{
   if ( argc < 3 || "-h" == std::string(argv[1]) || "--help" == std::string(argv[1]) ) {
         fprintf(stderr, kCommandLineOptionsHelp);
         return 1;
   }

   ROOT::TIOFeatures features;
   Bool_t append = kFALSE;
   Bool_t force = kFALSE;
   Bool_t skip_errors = kFALSE;
   Bool_t reoptimize = kFALSE;
   Bool_t noTrees = kFALSE;
   Bool_t keepCompressionAsIs = kFALSE;
   Bool_t useFirstInputCompression = kFALSE;
   Bool_t multiproc = kFALSE;
   Bool_t debug = kFALSE;
   Int_t maxopenedfiles = 0;
   Int_t verbosity = 99;
   TString cacheSize;
   SysInfo_t s;
   gSystem->GetSysInfo(&s);
   auto nProcesses = s.fCpus;
   auto workingDir = gSystem->TempDirectory();
   int outputPlace = 0;
   int ffirst = 2;
   Int_t newcomp = -1;
   for( int a = 1; a < argc; ++a ) {
      if ( strcmp(argv[a],"-T") == 0 ) {
         noTrees = kTRUE;
         ++ffirst;
      } else if ( strcmp(argv[a],"-a") == 0 ) {
         append = kTRUE;
         ++ffirst;
      } else if ( strcmp(argv[a],"-f") == 0 ) {
         force = kTRUE;
         ++ffirst;
      } else if ( strcmp(argv[a],"-k") == 0 ) {
         skip_errors = kTRUE;
         ++ffirst;
      } else if ( strcmp(argv[a],"-O") == 0 ) {
         reoptimize = kTRUE;
         ++ffirst;
      } else if (strcmp(argv[a], "-dbg") == 0) {
         debug = kTRUE;
         verbosity = kTRUE;
         ++ffirst;
      } else if (strcmp(argv[a], "-d") == 0) {
         if (a + 1 != argc && argv[a + 1][0] != '-') {
            if (gSystem->AccessPathName(argv[a + 1])) {
               std::cerr << "Error: could not access the directory specified: " << argv[a + 1]
                         << ". We will use the system's temporal directory.\n";
            } else {
               workingDir = argv[a + 1];
            }
            ++a;
            ++ffirst;
         } else {
            std::cout << "-d: no directory specified.  We will use the system's temporal directory.\n";
         }
         ++ffirst;
      } else if (strcmp(argv[a], "-j") == 0) {
         // If the number of processes is not specified, use the default.
         if (a + 1 != argc && argv[a + 1][0] != '-') {
            // number of processes specified
            Long_t request = 1;
            for (char *c = argv[a + 1]; *c != '\0'; ++c) {
               if (!isdigit(*c)) {
                  // Wrong number of Processes. Use the default:
                  std::cerr << "Error: could not parse the number of processes to run in parallel passed after -j: "
                            << argv[a + 1] << ". We will use the system maximum.\n";
                  request = 0;
                  break;
               }
            }
            if (request == 1) {
               request = strtol(argv[a + 1], 0, 10);
               if (request < kMaxLong && request >= 0) {
                  nProcesses = (Int_t)request;
                  ++a;
                  ++ffirst;
                  std::cout << "Parallelizing  with " << nProcesses << " processes.\n";
               } else {
                  std::cerr << "Error: could not parse the number of processes to use passed after -j: " << argv[a + 1]
                            << ". We will use the default value (number of logical cores).\n";
               }
            }
         }
         multiproc = kTRUE;
         ++ffirst;
      } else if ( strcmp(argv[a],"-cachesize=") == 0 ) {
         int size;
         static const size_t arglen = strlen("-cachesize=");
         auto parseResult = ROOT::FromHumanReadableSize(argv[a]+arglen,size);
         if (parseResult == ROOT::EFromHumanReadableSize::kParseFail) {
            std::cerr << "Error: could not parse the cache size passed after -cachesize: "
                      << argv[a + 1] << ". We will use the default value.\n";
         } else if (parseResult == ROOT::EFromHumanReadableSize::kOverflow) {
            double m;
            const char *munit = nullptr;
            ROOT::ToHumanReadableSize(INT_MAX,false,&m,&munit);
            std::cerr << "Error: the cache size passed after -cachesize is too large: "
                      << argv[a + 1] << " is greater than " << m << munit
                      << ". We will use the default value.\n";
         } else {
            cacheSize = "cachesize=";
            cacheSize.Append(argv[a]+1);
         }
         ++ffirst;
      } else if ( strcmp(argv[a],"-cachesize") == 0 ) {
         if (a+1 >= argc) {
            std::cerr << "Error: no cache size number was provided after -cachesize.\n";
         } else {
            int size;
            auto parseResult = ROOT::FromHumanReadableSize(argv[a+1],size);
            if (parseResult == ROOT::EFromHumanReadableSize::kParseFail) {
               std::cerr << "Error: could not parse the cache size passed after -cachesize: "
                         << argv[a + 1] << ". We will use the default value.\n";
            } else if (parseResult == ROOT::EFromHumanReadableSize::kOverflow) {
               double m;
               const char *munit = nullptr;
               ROOT::ToHumanReadableSize(INT_MAX,false,&m,&munit);
               std::cerr << "Error: the cache size passed after -cachesize is too large: "
                         << argv[a + 1] << " is greater than " << m << munit
                         << ". We will use the default value.\n";
               ++a;
               ++ffirst;
            } else {
               cacheSize = "cachesize=";
               cacheSize.Append(argv[a+1]);
               ++a;
               ++ffirst;
            }
         }
         ++ffirst;
      } else if (!strcmp(argv[a], "-experimental-io-features")) {
         if (a+1 >= argc) {
            std::cerr << "Error: no IO feature was specified after -experimental-io-features; ignoring\n";
         } else {
            std::stringstream ss;
            ss.str(argv[++a]);
            ++ffirst;
            std::string item;
            while (std::getline(ss, item, ',')) {
               if (!features.Set(item)) {
                  std::cerr << "Ignoring unknown feature request: " << item << std::endl;
               }
            }
         }
         ++ffirst;
      } else if ( strcmp(argv[a],"-n") == 0 ) {
         if (a+1 >= argc) {
            std::cerr << "Error: no maximum number of opened was provided after -n.\n";
         } else {
            Long_t request = strtol(argv[a+1], 0, 10);
            if (request < kMaxLong && request >= 0) {
               maxopenedfiles = (Int_t)request;
               ++a;
               ++ffirst;
            } else {
               std::cerr << "Error: could not parse the max number of opened file passed after -n: " << argv[a+1] << ". We will use the system maximum.\n";
            }
         }
         ++ffirst;
      } else if ( strcmp(argv[a],"-v") == 0 ) {
         if (a+1 == argc || argv[a+1][0] == '-') {
            // Verbosity level was not specified use the default:
            verbosity = 99;
//         if (a+1 >= argc) {
//            std::cerr << "Error: no verbosity level was provided after -v.\n";
         } else {
            Bool_t hasFollowupNumber = kTRUE;
            for (char *c = argv[a+1]; *c != '\0'; ++c) {
               if (!isdigit(*c)) {
                  // Verbosity level was not specified use the default:
                  hasFollowupNumber = kFALSE;
                  break;
               }
            }
            if (hasFollowupNumber) {
               Long_t request = strtol(argv[a+1], 0, 10);
               if (request < kMaxLong && request >= 0) {
                  verbosity = (Int_t)request;
                  ++a;
                  ++ffirst;
               } else {
                  verbosity = 99;
                  std::cerr << "Error: could not parse the verbosity level passed after -v: " << argv[a+1] << ". We will use the default value (99).\n";
               }
            }
         }
         ++ffirst;
      } else if ( argv[a][0] == '-' ) {
         bool farg = false;
         if (force && argv[a][1] == 'f') {
            // Bad argument
            std::cerr << "Error: Using option " << argv[a] << " more than once is not supported.\n";
            ++ffirst;
            farg = true;
         }
         const char *prefix = "";
         if (argv[a][1] == 'f' && argv[a][2] == 'k') {
            farg = true;
            force = kTRUE;
            keepCompressionAsIs = kTRUE;
            prefix = "k";
         }
         if (argv[a][1] == 'f' && argv[a][2] == 'f') {
            farg = true;
            force = kTRUE;
            useFirstInputCompression = kTRUE;
            if (argv[a][3] != '\0') {
               std::cerr << "Error: option -ff should not have any suffix: " << argv[a] << " (suffix has been ignored)\n";
            }
         }
         char ft[7];
         for (int alg = 0; !useFirstInputCompression && alg <= 5; ++alg) {
            for( int j=0; j<=9; ++j ) {
               const int comp = (alg*100)+j;
               snprintf(ft,7,"-f%s%d",prefix,comp);
               if (!strcmp(argv[a],ft)) {
                  farg = true;
                  force = kTRUE;
                  newcomp = comp;
                  break;
               }
            }
         }
         if (!farg) {
            // Bad argument
            std::cerr << "Error: option " << argv[a] << " is not a supported option.\n";
         }
         ++ffirst;
      } else if (!outputPlace) {
         outputPlace = a;
      }
   }

   gSystem->Load("libTreePlayer");

   const char *targetname = 0;
   if (outputPlace) {
      targetname = argv[outputPlace];
   } else {
      targetname = argv[ffirst-1];
   }

   if (verbosity > 1) {
      std::cout << "hadd Target file: " << targetname << std::endl;
   }

   TFileMerger fileMerger(kFALSE, kFALSE);
   fileMerger.SetMsgPrefix("hadd");
   fileMerger.SetPrintLevel(verbosity - 1);
   if (maxopenedfiles > 0) {
      fileMerger.SetMaxOpenedFiles(maxopenedfiles);
   }
   if (newcomp == -1) {
      if (useFirstInputCompression || keepCompressionAsIs) {
         // grab from the first file.
         TFile *firstInput = nullptr;
         if (argv[ffirst] && argv[ffirst][0]=='@') {
            std::ifstream indirect_file(argv[ffirst]+1);
            if( ! indirect_file.is_open() ) {
               std::cerr<< "hadd could not open indirect file " << (argv[ffirst]+1) << std::endl;
               return 1;
            }
            std::string line;
            while( indirect_file ){
               if( std::getline(indirect_file, line) && line.length() ) {
                  firstInput = TFile::Open(line.c_str());
                  break;
               }
            }
         } else {
            firstInput = TFile::Open(argv[ffirst]);
         }
         if (firstInput && !firstInput->IsZombie())
            newcomp = firstInput->GetCompressionSettings();
         else
            newcomp = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault % 100;
         delete firstInput;
      } else newcomp = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault % 100; // default compression level.
   }
   if (verbosity > 1) {
      if (keepCompressionAsIs && !reoptimize)
         std::cout << "hadd compression setting for meta data: " << newcomp << '\n';
      else
         std::cout << "hadd compression setting for all output: " << newcomp << '\n';
   }
   if (append) {
      if (!fileMerger.OutputFile(targetname, "UPDATE", newcomp)) {
         std::cerr << "hadd error opening target file for update :" << argv[ffirst-1] << "." << std::endl;
         exit(2);
      }
   } else if (!fileMerger.OutputFile(targetname, force, newcomp)) {
      std::cerr << "hadd error opening target file (does " << argv[ffirst-1] << " exist?)." << std::endl;
      if (!force) std::cerr << "Pass \"-f\" argument to force re-creation of output file." << std::endl;
      exit(1);
   }

   auto filesToProcess = argc - ffirst;
   auto step = (filesToProcess + nProcesses - 1) / nProcesses;
   if (multiproc && step < 3) {
      // At least 3 files per process
      step = 3;
      nProcesses = (filesToProcess + step - 1) / step;
      std::cout << "Each process should handle at least 3 files for efficiency.";
      std::cout << " Setting the number of processes to: " << nProcesses << std::endl;
   }
   if (nProcesses == 1)
      multiproc = kFALSE;

   std::vector<std::string> partialFiles;

#ifndef R__WIN32
   // this is commented out only to try to prevent false positive detection
   // from several anti-virus engines on Windows, and multiproc is not
   // supported on Windows anyway
   if (multiproc) {
      auto uuid = TUUID();
      auto partialTail = uuid.AsString();
      for (auto i = 0; (i * step) < filesToProcess; i++) {
         std::stringstream buffer;
         buffer << workingDir << "/partial" << i << "_" << partialTail << ".root";
         partialFiles.emplace_back(buffer.str());
      }
   }
#endif

   auto mergeFiles = [&](TFileMerger &merger) {
      if (reoptimize) {
         merger.SetFastMethod(kFALSE);
      } else {
         if (!keepCompressionAsIs && merger.HasCompressionChange()) {
            // Don't warn if the user any request re-optimization.
            std::cout << "hadd Sources and Target have different compression levels" << std::endl;
            std::cout << "hadd merging will be slower" << std::endl;
         }
      }
      merger.SetNotrees(noTrees);
      merger.SetMergeOptions(cacheSize);
      merger.SetIOFeatures(features);
      Bool_t status;
      if (append)
         status = merger.PartialMerge(TFileMerger::kIncremental | TFileMerger::kAll);
      else
         status = merger.Merge();
      return status;
   };

   auto sequentialMerge = [&](TFileMerger &merger, int start, int nFiles) {

      for (auto i = start; i < (start + nFiles) && i < argc; i++) {
         if (argv[i] && argv[i][0] == '@') {
            std::ifstream indirect_file(argv[i] + 1);
            if (!indirect_file.is_open()) {
               std::cerr << "hadd could not open indirect file " << (argv[i] + 1) << std::endl;
               return kFALSE;
            }
            while (indirect_file) {
               std::string line;
               if (std::getline(indirect_file, line) && line.length() && !merger.AddFile(line.c_str())) {
                  return kFALSE;
               }
            }
         } else if (!merger.AddFile(argv[i])) {
            if (skip_errors) {
               std::cerr << "hadd skipping file with error: " << argv[i] << std::endl;
            } else {
               std::cerr << "hadd exiting due to error in " << argv[i] << std::endl;
               return kFALSE;
            }
         }
      }
      return mergeFiles(merger);
   };

   auto parallelMerge = [&](int start) {
      TFileMerger mergerP(kFALSE, kFALSE);
      mergerP.SetMsgPrefix("hadd");
      mergerP.SetPrintLevel(verbosity - 1);
      if (maxopenedfiles > 0) {
         mergerP.SetMaxOpenedFiles(maxopenedfiles / nProcesses);
      }
      if (!mergerP.OutputFile(partialFiles[(start - ffirst) / step].c_str(), newcomp)) {
         std::cerr << "hadd error opening target partial file" << std::endl;
         exit(1);
      }
      return sequentialMerge(mergerP, start, step);
   };

   auto reductionFunc = [&]() {
      for (const auto &pf : partialFiles) {
         fileMerger.AddFile(pf.c_str());
      }
      return mergeFiles(fileMerger);
   };

   Bool_t status;

#ifndef R__WIN32
   if (multiproc) {
      ROOT::TProcessExecutor p(nProcesses);
      auto res = p.Map(parallelMerge, ROOT::TSeqI(ffirst, argc, step));
      status = std::accumulate(res.begin(), res.end(), 0U) == partialFiles.size();
      if (status) {
         status = reductionFunc();
      } else {
         std::cout << "hadd failed at the parallel stage" << std::endl;
      }
      if (!debug) {
         for (const auto &pf : partialFiles) {
            gSystem->Unlink(pf.c_str());
         }
      }
   } else {
      status = sequentialMerge(fileMerger, ffirst, filesToProcess);
   }
#else
   status = sequentialMerge(fileMerger, ffirst, filesToProcess);
#endif

   if (status) {
      if (verbosity == 1) {
         std::cout << "hadd merged " << fileMerger.GetMergeList()->GetEntries() << " input files in " << targetname
                   << ".\n";
      }
      return 0;
   } else {
      if (verbosity == 1) {
         std::cout << "hadd failure during the merge of " << fileMerger.GetMergeList()->GetEntries()
                   << " input files in " << targetname << ".\n";
      }
      return 1;
   }
}
