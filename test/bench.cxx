// this test program compares the I/O performance obtained with
// all STL collections of objects or pointers to objects and also
// Root collection class TClonesArray.
// Trees in compression and non compression mode are created for each
// of the following cases:
//  -STLcollection<THit>
//  -STLcollection<THit*>
//  -TClonesArray(TObjHit) in no split mode
//  -TClonesArray(TObjHit) in split mode
// where:
//  THit is a class not derived from TObject
//  TObjHit derives from TObject and THit
//
//  run with
//     bench
//   or
//     bench -m   to stream objects memberwise
//
// The test prints a summary table comparing performances for all above cases
// (CPU, file size, compression factors).
// Reference numbers on a Pentium IV 2.4 Ghz machine are given as reference.
//      Authors:  Rene Brun, Markus Frank

#include "TROOT.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TStreamerInfo.h"
#include "snprintf.h"

#include "TBench.h"

using namespace std;

struct TBenchData {
   TBenchData() : cp1(0), nbytes1(0), cp2w(0), cp2r(0), cx3(0), nbytes3(0), cp3w(0), cp3r(0)  {}
   TBenchData(const char *name, Double_t i_cp1, Float_t i_cx3, Long64_t i_nbytes1, Long64_t i_nbytes3, Double_t i_cp2w, Double_t i_cp3w, Double_t i_cp2r, Double_t i_cp3r)
     :  fName(name), cp1(i_cp1), nbytes1(i_nbytes1), cp2w(i_cp2w), cp2r(i_cp2r), cx3(i_cx3), nbytes3(i_nbytes3), cp3w(i_cp3w), cp3r(i_cp3r) {}
   TString fName;
   Double_t rt1;
   Double_t cp1;
   Float_t cx;
   Long64_t nbytes1;
   Double_t rt2w;
   Double_t cp2w;
   Double_t rt2r;
   Double_t cp2r;
   Float_t cx3;
   Long64_t nbytes3;
   Double_t rt3w;
   Double_t cp3w;
   Double_t rt3r;
   Double_t cp3r;

   Double_t cptot() { return cp1 + cp2w + cp2r + cp3w + cp3r; }
};

template <class TGen> TBenchData runTest(const char *name, int nevents, int nhits, int splitlevel)
{
   static TStopwatch timer;

   TBenchData data;
   data.fName = name;

   timer.Start();
   TGen *STLhit = new TGen(nhits);
   STLhit->MakeTree(0,nevents,0,0,data.cx);
   timer.Stop();
   data.rt1 = timer.RealTime();
   data.cp1 = timer.CpuTime();
   printf("1 %-26s  : RT=%6.2f s  Cpu=%6.2f s\n",data.fName.Data(),data.rt1,data.cp1);

   timer.Start(kTRUE);
   data.nbytes1 = STLhit->MakeTree(1,nevents,0,splitlevel,data.cx);
   timer.Stop();
   data.rt2w = timer.RealTime();
   data.cp2w = timer.CpuTime();
   printf("2 %-26s w: RT=%6.2f s  Cpu=%6.2f s, size= %8lld bytes, cx=%5.2f\n",data.fName.Data(),data.rt2w,data.cp2w,data.nbytes1,data.cx);

   timer.Start(kTRUE);
   STLhit->ReadTree();
   timer.Stop();
   data.rt2r = timer.RealTime();
   data.cp2r = timer.CpuTime();
   printf("3 %-26s r: RT=%6.2f s  Cpu=%6.2f s\n",data.fName.Data(),data.rt2r,data.cp2r);

   timer.Start(kTRUE);
   data.nbytes3 = STLhit->MakeTree(1,nevents,1,splitlevel,data.cx3);
   timer.Stop();
   data.rt3w = timer.RealTime();
   data.cp3w = timer.CpuTime();
   printf("4 %-26s w: RT=%6.2f s  Cpu=%6.2f s, size= %8lld bytes, cx=%5.2f\n",data.fName.Data(),data.rt3w,data.cp3w,data.nbytes3,data.cx3);

   timer.Start(kTRUE);
   STLhit->ReadTree();
   timer.Stop();
   data.rt3r = timer.RealTime();
   data.cp3r = timer.CpuTime();
   printf("5 %-26s r: RT=%6.2f s  Cpu=%6.2f s\n",data.fName.Data(),data.rt3r,data.cp3r);

   delete STLhit;
   return data;
}

template <class TGen> TBenchData runTest(const char *name, int nevents, int nhits, int splitlevel, Double_t &cptot, vector<TBenchData> &results)
{
   TBenchData data = runTest<TGen>( name, nevents, nhits, splitlevel);
   cptot += data.cptot();
   results.push_back( data );
   return data;
}

int main(int argc, char **argv)
{
   std::string inclRootSys = ("-I" + TROOT::GetRootSys() + "/test").Data();
   TROOT::AddExtraInterpreterArgs({inclRootSys});

   bool writereferences = false;
   bool memberwise = false;
   bool shortrun = false;

   TVirtualStreamerInfo::SetStreamMemberWise(kFALSE);
   // by default stream objects objectwise
   // if program option "-m" is specified, stream memberwise
   for(int a=1; a<argc; ++a) {
      if (strstr(argv[a],"-m")) {
         TVirtualStreamerInfo::SetStreamMemberWise(kTRUE);
         printf("bench option -m specified. Streaming objects memberwise\n");
         memberwise = true;
      } else if (strstr(argv[a],"-r")) {
         writereferences = true;
      } else if (strstr(argv[a],"-s")) {
         shortrun = true;
      }
   }
   int nhits       = shortrun ? 10 : 1000;
   int nevents     = shortrun ? 40 : 400;

   Double_t cptot = 0;

   //delete temp file used for the benchmark
   gSystem->Unlink("bench.root");

   vector<TBenchData> results;
   vector<TBenchData> references;
   references.push_back( TBenchData( "vector<THit> level=99",   0.44,   5.37, 39724266, 7392775,   0.96,   2.15,   0.32,   0.73 ) );
   references.push_back( TBenchData( "vector<THit> level= 0",   0.29,   3.76, 42053455, 11177799,   0.85,   3.22,   0.68,   0.93 ) );
   references.push_back( TBenchData( "vector<THit> level= 0 MW",   0.28,   4.96, 39654255, 7997901,   0.63,   1.91,   0.33,   0.52 ) );
   references.push_back( TBenchData( "list<THit> level=99",   0.36,   5.37, 39724301, 7393378,   0.75,   2.33,   0.50,  0.64 ) );
   references.push_back( TBenchData( "list<THit> level= 0",   0.37,   3.76, 42053474, 11177716,   0.94,   2.96,   0.64,   1.02 ) );
   references.push_back( TBenchData( "list<THit> level= 0 MW",   0.37,   4.96, 39654274, 7998366,   0.76,   2.26,   0.37,   0.63 ) );
   references.push_back( TBenchData( "deque<THit> level=99",   0.38,   5.37, 39724311, 7392970,   0.70,   2.01,   0.30,   0.79 ) );
   references.push_back( TBenchData( "deque<THit> level= 0",   0.38,   3.76, 42053480, 11177240,   0.95,   2.70,   0.81,   0.97 ) );
   references.push_back( TBenchData( "deque<THit> level= 0 MW",   0.37,   4.96, 39654280, 7998569,   0.76,   2.07,   0.42,   0.88 ) );
   references.push_back( TBenchData( "set<THit> level=99",   0.41,   6.21, 39724292, 6396267,   0.81,   2.09,   0.62, 0.71 ) );
   references.push_back( TBenchData( "set<THit> level= 0",   0.41,   4.39, 42053469, 9573788,   0.99,   2.65,   0.89, 1.10 ) );
   references.push_back( TBenchData( "set<THit> level= 0 MW",   0.41,   5.79, 39654269, 6851326,   0.85,   1.62,   0.45,   0.69 ) );
   references.push_back( TBenchData( "multiset<THit> level=99",   0.44,   6.21, 39724346, 6396195,   0.83,   2.25,   0.46,   0.71 ) );
   references.push_back( TBenchData( "multiset<THit> level= 0",   0.43,   4.39, 42053503, 9573512,   1.04,   2.83,   0.76,   1.09 ) );
   references.push_back( TBenchData( "multiset<THit> level= 0 MW",   0.44,   5.79, 39654303, 6851637,   0.88,   2.03, 0.62,   0.69 ) );
   references.push_back( TBenchData( "map<THit> level=99",   0.48,   5.26, 41335814, 7861427,   0.91,   2.39,   0.67, 0.77 ) );
   references.push_back( TBenchData( "map<THit> level= 0",   0.48,   3.84, 43653478, 11357815,   1.08,   2.69,   1.14,  1.19 ) );
   references.push_back( TBenchData( "map<THit> level= 0 MW",   0.48,   3.68, 43655963, 11866052,   1.18,   2.94,   1.07,   1.17 ) );
   references.push_back( TBenchData( "multimap<THit> level=99",   0.49,   5.26, 41335869, 7861144,   0.90,   2.52,   0.48,   0.78 ) );
   references.push_back( TBenchData( "multimap<THit> level= 0",   0.49,   3.84, 43653510, 11358500,   1.09,   2.72,   1.13,   1.18 ) );
   references.push_back( TBenchData( "multimap<THit> level= 0 MW",   0.48,   3.68, 43655998, 11866816,   1.15,   2.95,  1.04,   1.17 ) );
   references.push_back( TBenchData( "vector<THit*> level=25599",   0.25,   5.14, 40137890, 7806511,   0.57,   1.78,0.33,   0.87 ) );
   references.push_back( TBenchData( "list<THit*> level=25599",   0.33,   5.14, 40137936, 7806412,   0.63,   2.13,   0.42,   0.69 ) );
   references.push_back( TBenchData( "deque<THit*> level=25599",   0.31,   5.14, 40137945, 7806520,   0.67,   1.89,   0.35,   0.88 ) );
   references.push_back( TBenchData( "set<THit*> level=25599",   0.38,   5.74, 40137925, 6990746,   0.72,   2.20,   0.70,   0.95 ) );
   references.push_back( TBenchData( "multiset<THit*> level=25599",   0.36,   5.85, 40137981, 6865586,   0.73,   2.17,  0.68,   0.86 ) );
   references.push_back( TBenchData( "map<THit*> level=99",   0.37,   3.71, 46864802, 12643706,   1.28,   3.24,   1.34,   2.11 ) );
   references.push_back( TBenchData( "multimap<THit*> level=99",   0.36,   3.71, 46864850, 12642782,   1.28,   3.20,1.38,   1.58 ) );
   references.push_back( TBenchData( "vector<THit*> level=99 (NS)",   0.24,   3.72, 45257044, 12170238,   1.08,   3.05,   1.43,   1.44 ) );
   references.push_back( TBenchData( "list<THit*> level=99 (NS)",   0.30,   3.72, 45257088, 12169539,   1.14,   2.88, 1.45,   1.49 ) );
   references.push_back( TBenchData( "deque<THit*> level=99 (NS)",   0.32,   3.72, 45257098, 12167790,   1.16,   2.69,  1.41,   1.44 ) );
   references.push_back( TBenchData( "set<THit*> level=99 (NS)",   0.38,   4.29, 45257081, 10550593,   1.23,   3.01,1.28,   2.14 ) );
   references.push_back( TBenchData( "multiset<THit*> level=99 (NS)",   0.45,   4.28, 45257136, 10572288,   1.27,   3.03,   1.61,   1.50 ) );
   references.push_back( TBenchData( "TClonesArray(TObjHit) level= 0",   0.28,   4.96, 39666751, 8006156,   0.55,   1.76,   0.16,   0.64 ) );
   references.push_back( TBenchData( "TClonesArray(TObjHit) level=99",   0.37,   5.37, 39722797, 7392577,   0.55,   2.07,   0.17,   0.46 ) );

   /// STL VECTOR
   runTest<TSTLhit>( "vector<THit> level=99", nevents, nhits, 99, cptot, results );
   /// STL VECTOR not split.
   runTest<TSTLhit>( "vector<THit> level= 0", nevents, nhits,  0, cptot, results );

   /// STL VECTOR not split, member wise mode
   memberwise = TVirtualStreamerInfo::SetStreamMemberWise(true);
   runTest<TSTLhit>( "vector<THit> level= 0 MW", nevents, nhits,  0, cptot, results );
   TVirtualStreamerInfo::SetStreamMemberWise(memberwise);

   // STL list
   runTest<TSTLhitList>( "list<THit> level=99", nevents, nhits, 99, cptot, results );
   runTest<TSTLhitList>( "list<THit> level= 0", nevents, nhits,  0, cptot, results );
   memberwise = TVirtualStreamerInfo::SetStreamMemberWise(true);
   runTest<TSTLhitList>( "list<THit> level= 0 MW", nevents, nhits,  0, cptot, results );
   TVirtualStreamerInfo::SetStreamMemberWise(memberwise);

   // STL DEQUE
   runTest<TSTLhitDeque>( "deque<THit> level=99", nevents, nhits, 99, cptot, results );
   runTest<TSTLhitDeque>( "deque<THit> level= 0", nevents, nhits,  0, cptot, results );
   memberwise = TVirtualStreamerInfo::SetStreamMemberWise(true);
   runTest<TSTLhitDeque>( "deque<THit> level= 0 MW", nevents, nhits,  0, cptot, results );
   TVirtualStreamerInfo::SetStreamMemberWise(memberwise);

   // STL SET
   runTest<TSTLhitSet>( "set<THit> level=99", nevents, nhits, 99, cptot, results );
   runTest<TSTLhitSet>( "set<THit> level= 0", nevents, nhits,  0, cptot, results );
   memberwise = TVirtualStreamerInfo::SetStreamMemberWise(true);
   runTest<TSTLhitSet>( "set<THit> level= 0 MW", nevents, nhits,  0, cptot, results );
   TVirtualStreamerInfo::SetStreamMemberWise(memberwise);

   // STL MULTI SET
   runTest<TSTLhitMultiset>( "multiset<THit> level=99", nevents, nhits, 99, cptot, results );
   runTest<TSTLhitMultiset>( "multiset<THit> level= 0", nevents, nhits,  0, cptot, results );
   memberwise = TVirtualStreamerInfo::SetStreamMemberWise(true);
   runTest<TSTLhitMultiset>( "multiset<THit> level= 0 MW", nevents, nhits,  0, cptot, results );
   TVirtualStreamerInfo::SetStreamMemberWise(memberwise);

   // STL map
   runTest<TSTLhitMap>( "map<THit> level=99", nevents, nhits, 99, cptot, results );
   runTest<TSTLhitMap>( "map<THit> level= 0", nevents, nhits,  0, cptot, results );
   memberwise = TVirtualStreamerInfo::SetStreamMemberWise(true);
   runTest<TSTLhitMap>( "map<THit> level= 0 MW", nevents, nhits,  0, cptot, results );
   TVirtualStreamerInfo::SetStreamMemberWise(memberwise);

   // STL multimap
   runTest<TSTLhitMultiMap>( "multimap<THit> level=99", nevents, nhits, 99, cptot, results );
   runTest<TSTLhitMultiMap>( "multimap<THit> level= 0", nevents, nhits,  0, cptot, results );
   memberwise = TVirtualStreamerInfo::SetStreamMemberWise(true);
   runTest<TSTLhitMultiMap>( "multimap<THit> level= 0 MW", nevents, nhits,  0, cptot, results );
   TVirtualStreamerInfo::SetStreamMemberWise(memberwise);

   //__________________________________________________________________________
   //
   //testing STL vector of pointers to THit
   runTest<TSTLhitStar>( "vector<THit*> level=25599", nevents, nhits, 25599, cptot, results );

   // STL list*
   runTest<TSTLhitStarList>( "list<THit*> level=25599", nevents, nhits, 25599, cptot, results );

   // STL DEQUE*
   runTest<TSTLhitStarDeque>( "deque<THit*> level=25599", nevents, nhits, 25599, cptot, results );

   // STL SET*
   runTest<TSTLhitStarSet>( "set<THit*> level=25599", nevents, nhits, 25599, cptot, results );

   // STL MULTI SET*
   runTest<TSTLhitStarMultiSet>( "multiset<THit*> level=25599", nevents, nhits, 25599, cptot, results );

   // STL MAP*
   runTest<TSTLhitStarMap>( "map<THit*> level=99", nevents, nhits, 99, cptot, results );

   // STL MULTIMAP*
   runTest<TSTLhitStarMultiMap>( "multimap<THit*> level=99", nevents, nhits, 99, cptot, results );

   //__________________________________________________________________________
   //
   //testing STL vector of pointers to THit (NOSPLIT)
   runTest<TSTLhitStar>( "vector<THit*> level=99 (NS)", nevents, nhits, 99, cptot, results );

   // STL list* (NOSPLIT)
   runTest<TSTLhitStarList>( "list<THit*> level=99 (NS)", nevents, nhits, 99, cptot, results );

   // STL DEQUE* (NOSPLIT)
   runTest<TSTLhitStarDeque>( "deque<THit*> level=99 (NS)", nevents, nhits, 99, cptot, results );

   // STL SET* (NOSPLIT)
   runTest<TSTLhitStarSet>( "set<THit*> level=99 (NS)", nevents, nhits, 99, cptot, results );

   // STL MULTI SET* (NOSPLIT)
   runTest<TSTLhitStarMultiSet>( "multiset<THit*> level=99 (NS)", nevents, nhits, 99, cptot, results );

   //___________________________________________________________________________
   //
   //testing TClonesArray of TObjHit deriving from THit
   runTest<TCloneshit>( "TClonesArray(TObjHit) level= 0", nevents, nhits,  0, cptot, results );
   runTest<TCloneshit>( "TClonesArray(TObjHit) level=99", nevents, nhits, 99, cptot, results );

   Double_t cpref = 195.43;
   Double_t rootmarks = cpref*800/cptot;

   for(unsigned int t=references.size(); t<results.size(); ++t) {
      references.push_back(TBenchData());
   }

   //print all results
   char line1[100], line2[100];
   printf("\n");
   printf("*******************************************************************************\n");
   snprintf(line1,100,"Comparing STL vector with TClonesArray: Root %-8s",gROOT->GetVersion());
   printf("*       %s                 *\n",line1);
   Bool_t UNIX = strcmp(gSystem->GetName(), "Unix") == 0;
   if (UNIX) {
      TString sp = gSystem->GetFromPipe("uname -a");
      sp.Resize(60);
      printf("*  SYS: %s\n",sp.Data());
      if (strstr(gSystem->GetBuildNode(),"Linux")) {
         sp = gSystem->GetFromPipe("lsb_release -d -s");
         printf("*  SYS: %s\n",sp.Data());
      }
      if (strstr(gSystem->GetBuildNode(),"Darwin")) {
         sp  = gSystem->GetFromPipe("sw_vers -productVersion");
         sp += " Mac OS X ";
         printf("*  SYS: %s\n",sp.Data());
      }
   } else {
      const char *os = gSystem->Getenv("OS");
      snprintf(line2,100,"Windows");
      if (!os) printf("*  Windows 95\n");
      else     printf("*  %s %s \n",os,gSystem->Getenv("PROCESSOR_IDENTIFIER"));
   }
   printf("*     Reference machine pcbrun4.cern.ch  RedHat Linux 7.3                     *\n");
   printf("*         (Pentium IV 2.4 Ghz 512 Mbytes RAM, IDE disk)                       *\n");
   printf("*           (send your results to rootdev@cern.ch)                            *\n");
   printf("*******************************************************************************\n");
   printf("* Time to fill the structures  (seconds)   Reference      cx      Reference   *\n");
   printf("*******************************************************************************\n");
   for(unsigned int t=0; t<results.size() && t<references.size(); ++t) {
      printf("* %-30s %6.2f       %5.2f     %5.2f       %5.2f     *\n",results[t].fName.Data(),results[t].cp1,references[t].cp1,results[t].cx3,references[t].cx3);
   }
   printf("*******************************************************************************\n");
   printf("* Size of file in bytes          comp 0    Reference    comp 1    Reference   *\n");
   printf("*******************************************************************************\n");
   for(unsigned int t=0; t<results.size() && t<references.size(); ++t) {
      printf("* %-30s %8lld   %8lld   %8lld  %8lld     *\n",results[t].fName.Data(),results[t].nbytes1,references[t].nbytes1,results[t].nbytes3,references[t].nbytes3);
   }
   printf("*******************************************************************************\n");
   printf("* Time to write in seconds       comp 0    Reference    comp 1    Reference   *\n");
   printf("*******************************************************************************\n");
   for(unsigned int t=0; t<results.size() && t<references.size(); ++t) {
      printf("* %-30s %6.2f      %6.2f    %6.2f      %6.2f     *\n",results[t].fName.Data(),results[t].cp2w,references[t].cp2w, results[t].cp3w, references[t].cp3w);
   }
   printf("*******************************************************************************\n");
   printf("* Time to read in seconds        comp 0    Reference    comp 1    Reference   *\n");
   printf("*******************************************************************************\n");
   for(unsigned int t=0; t<results.size() && t<references.size(); ++t) {
      printf("* %-30s %6.2f      %6.2f    %6.2f      %6.2f     *\n",results[t].fName.Data(),results[t].cp2r,references[t].cp2r,results[t].cp3r,references[t].cp3r);
   }
   printf("*******************************************************************************\n");
   printf("* Total CPU time              %8.2f    %8.2f                            *\n",cptot,cpref);
   printf("* Estimated ROOTMARKS         %8.2f      800.00                            *\n",rootmarks);
   printf("******************************************************************************\n");

   if (writereferences) {
      for(unsigned int t=0; t<results.size() && t<references.size(); ++t) {
         printf("references.push_back( TBenchData( \"%s\", %6.2f, %6.2f, %lld, %lld, %6.2f, %6.2f, %6.2f, %6.2f ) );\n",
                results[t].fName.Data(),results[t].cp1,results[t].cx3,results[t].nbytes1,results[t].nbytes3,results[t].cp2w,results[t].cp3w,results[t].cp2r,results[t].cp3r);
      }
   }
   return 0;
}
