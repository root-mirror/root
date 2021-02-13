// @(#)root/test:$name:  $:$id: stressLinear.cxx,v 1.15 2002/10/25 10:47:51 rdm exp $
// Authors: Fons Rademakers, Eddy Offermann  Jan 2004

/////////////////////////////////////////////////////////////////
//
//    R O O T   LINEAR ALGEBRA T E S T  S U I T E  and  B E N C H M A R K S
//    =====================================================================
//
// The suite of programs below tests many elements of the vector, matrix
// and matrix decomposition classes.
// To run in batch, do
//   stressLinear      : with no parameters, run standard test with matrices upto 100x100
//   stressLinear 30   : run test with matrices upto 30x30
//   stressLinear 30 1 : run test with matrices upto 30x30 in verbose mode
//
// To run interactively, do
// root -b
//  Root > .x stressLinear.cxx(100)   run standard test with matrices upto 100x100
//  Root > .x stressLinear.cxx(30)    run with matrices upto 30x30
//  Root > .x stressLinear.cxx(30,1)  run with matrices upto 30x30 in verbose mode
//
// Several tests are run sequentially. Each test will produce one line (Test OK or Test FAILED) .
// At the end of the test a table is printed showing the global results
// Real Time and Cpu Time.
// One single number (ROOTMARKS) is also calculated showing the relative
// performance of your machine compared to a reference machine
// a Pentium IV 2.4 Ghz) with 512 MBytes of memory
// and 120 GBytes IDE disk.
//
// An example of output when all the tests run OK is shown below:

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Linear Algebra Package -- Matrix Verifications.                      //
//                                                                      //
// This file implements a large set of TMat operation tests.            //
// *******************************************************************  //
// *  Linear Algebra - S T R E S S suite                             *  //
// *******************************************************************  //
// *******************************************************************  //
// *  Starting  Matrix - S T R E S S                                 *  //
// *******************************************************************  //
// Test  1 : Allocation, Resizing.................................. OK  //
// Test  2 : Filling, Inserting, Using............................. OK  //
// Test  3 : Uniform matrix operations............................. OK  //
// Test  4 : Binary Matrix element-by-element operations............OK  //
// Test  5 : Matrix transposition...................................OK  //
// Test  6 : Haar/Hilbert Matrix....................................OK  //
// Test  7 : Matrix promises........................................OK  //
// Test  8 : Matrix Norms...........................................OK  //
// Test  9 : Matrix Determinant.....................................OK  //
// Test 10 : General Matrix Multiplications.........................OK  //
// Test 11 : Symmetric Matrix Multiplications.......................OK  //
// Test 12 : Matrix Vector Multiplications..........................OK  //
// Test 13 : Matrix Inversion.......................................OK  //
// Test 14 : Matrix Persistence.....................................OK  //
// ******************************************************************   //
// *  Starting  Sparse Matrix - S T R E S S                         *   //
// ******************************************************************   //
// Test  1 : Allocation, Resizing...................................OK  //
// Test  2 : Filling, Inserting, Using..............................OK  //
// Test  3 : Uniform matrix operations..............................OK  //
// Test  4 : Binary Matrix element-by-element operations............OK  //
// Test  5 : Matrix transposition...................................OK  //
// Test  6 : Matrix Norms...........................................OK  //
// Test  7 : General Matrix Multiplications.........................OK  //
// Test  8 : Matrix Vector Multiplications..........................OK  //
// Test  9 : Matrix Slices to Vectors...............................OK  //
// Test 10 : Matrix Persistence.....................................OK  //
// *******************************************************************  //
// *  Starting  Vector - S T R E S S                                 *  //
// *******************************************************************  //
// Test  1 : Allocation, Filling, Resizing......................... OK  //
// Test  2 : Uniform vector operations............................. OK  //
// Test  3 : Binary vector element-by-element operations............OK  //
// Test  4 : Vector Norms...........................................OK  //
// Test  5 : Matrix Slices to Vectors...............................OK  //
// Test  6 : Vector Persistence.....................................OK  //
// *******************************************************************  //
// *  Starting  Linear Algebra - S T R E S S                         *  //
// *******************************************************************  //
// Test  1 : Decomposition / Reconstruction........................ OK  //
// Test  2 : Linear Equations...................................... OK  //
// Test  3 : Pseudo-Inverse, Moore-Penrose......................... OK  //
// Test  4 : Eigen - Values/Vectors.................................OK  //
// Test  5 : Decomposition Persistence..............................OK  //
// *******************************************************************  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <snprintf.h>
#include <TSystem.h>
#include <TFile.h>
#include <TBenchmark.h>
#include <TArrayD.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>

#include "TMatrixF.h"
#include "TMatrixFSym.h"
#include "TMatrixFSparse.h"
#include "TMatrixFLazy.h"
#include "TVectorF.h"

#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSparse.h"
#include "TMatrixDLazy.h"
#include "TMatrixDUtils.h"
#include "TVectorD.h"

#include "TDecompLU.h"
#include "TDecompChol.h"
#include "TDecompQRH.h"
#include "TDecompSVD.h"
#include "TDecompBK.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"

void stressLinear                  (Int_t maxSizeReq=100,Int_t verbose=0);
void StatusPrint                   (Int_t id,const TString &title,Bool_t status);

void mstress_allocation            (Int_t msize);
void mstress_matrix_fill           (Int_t rsize,Int_t csize);
void mstress_element_op            (Int_t rsize,Int_t csize);
void mstress_binary_ebe_op         (Int_t rsize, Int_t csize);
void mstress_transposition         (Int_t msize);
void mstress_special_creation      (Int_t dim);
void mstress_matrix_fill           (Int_t rsize,Int_t csize);
void mstress_matrix_promises       (Int_t dim);
void mstress_norms                 (Int_t rsize,Int_t csize);
void mstress_determinant           (Int_t msize);
void mstress_mm_multiplications    ();
void mstress_sym_mm_multiplications(Int_t msize);
void mstress_vm_multiplications    ();
void mstress_inversion             ();
void mstress_matrix_io             ();

void spstress_allocation           (Int_t msize);
void spstress_matrix_fill          (Int_t rsize,Int_t csize);
void spstress_element_op           (Int_t rsize,Int_t csize);
void spstress_binary_ebe_op        (Int_t rsize, Int_t csize);
void spstress_transposition        (Int_t msize);
void spstress_norms                (Int_t rsize,Int_t csize);
void spstress_mm_multiplications   ();
void spstress_vm_multiplications   ();
void spstress_matrix_slices        (Int_t vsize);
void spstress_matrix_io            ();

void vstress_allocation            (Int_t msize);
void vstress_element_op            (Int_t vsize);
void vstress_binary_op             (Int_t vsize);
void vstress_norms                 (Int_t vsize);
void vstress_matrix_slices         (Int_t vsize);
void vstress_vector_io             ();

Bool_t test_svd_expansion          (const TMatrixD &A);
void   astress_decomp              ();
void   astress_lineqn              ();
void   astress_pseudo              ();
void   astress_eigen               (Int_t msize);
void   astress_decomp_io           (Int_t msize);

void   stress_backward_io          ();

void   cleanup                     ();

//_____________________________batch only_____________________
#ifndef __CINT__

int main(int argc,const char *argv[])
{
  gROOT->SetBatch();
  Int_t maxSizeReq = 100;
  if (argc > 1)  maxSizeReq = atoi(argv[1]);
  Int_t verbose = 0;
  if (argc > 2)  verbose = (atoi(argv[2]) > 0 ? 1 : 0);

  gBenchmark = new TBenchmark();
  stressLinear(maxSizeReq,verbose);

  cleanup();
  return 0;
}

#endif

#define EPSILON 1.0e-14

Int_t gVerbose      = 0;
Int_t gNrLoop;
const Int_t nrSize  = 20;
const Int_t gSizeA[] = {5,6,7,8,9,10,15,20,30,40,50,60,70,80,90,100,300,500,700,1000};

//________________________________common part_________________________

void stressLinear(Int_t maxSizeReq,Int_t verbose)
{
  std::cout << "******************************************************************" <<std::endl;
  std::cout << "*  Starting  Linear Algebra - S T R E S S suite                  *" <<std::endl;
  std::cout << "******************************************************************" <<std::endl;
  std::cout << "******************************************************************" <<std::endl;

  gVerbose = verbose;

  gBenchmark->Start("stressLinear");

  gNrLoop = nrSize-1;
  while (gNrLoop > 0 && maxSizeReq < gSizeA[gNrLoop])
    gNrLoop--;

  const Int_t maxSize = gSizeA[gNrLoop];

  // Matrix
  {
    std::cout << "*  Starting  Matrix - S T R E S S                                *" <<std::endl;
    std::cout << "******************************************************************" <<std::endl;
    mstress_allocation(maxSize);
    mstress_matrix_fill(maxSize,maxSize/2);
    mstress_element_op(maxSize,maxSize/2);
    mstress_binary_ebe_op(maxSize/2,maxSize);
    mstress_transposition(maxSize);
    mstress_special_creation(maxSize);
#ifndef __CINT__
    mstress_matrix_promises(maxSize);
#endif
    mstress_norms(maxSize,maxSize/2);
    mstress_determinant(maxSize);
    mstress_mm_multiplications();
    mstress_sym_mm_multiplications(maxSize);
    mstress_vm_multiplications();
    mstress_inversion();

    mstress_matrix_io();
    std::cout << "******************************************************************" <<std::endl;
  }

  // Sparse Matrix
  {
    std::cout << "*  Starting  Sparse Matrix - S T R E S S                         *" <<std::endl;
    std::cout << "******************************************************************" <<std::endl;
    spstress_allocation(maxSize);
    spstress_matrix_fill(maxSize,maxSize/2);
    spstress_element_op(maxSize,maxSize/2);
    spstress_binary_ebe_op(maxSize/2,maxSize);
    spstress_transposition(maxSize);
    spstress_norms(maxSize,maxSize/2);
    spstress_mm_multiplications();
    spstress_vm_multiplications();
    spstress_matrix_slices(maxSize);
    spstress_matrix_io();
    std::cout << "******************************************************************" <<std::endl;
  }

  {
    std::cout << "*  Starting  Vector - S T R E S S                                *" <<std::endl;
    std::cout << "******************************************************************" <<std::endl;
    vstress_allocation(maxSize);
    vstress_element_op(maxSize);
    vstress_binary_op(maxSize);
    vstress_norms(maxSize);
    vstress_matrix_slices(maxSize);
    vstress_vector_io();
    std::cout << "******************************************************************" <<std::endl;
  }

  // Linear Algebra
  {
    std::cout << "*  Starting  Linear Algebra - S T R E S S                        *" <<std::endl;
    std::cout << "******************************************************************" <<std::endl;
    astress_decomp();
    astress_lineqn();
    astress_pseudo();
    astress_eigen(5);
    astress_decomp_io(10);
    std::cout << "******************************************************************" <<std::endl;
  }

  //Backward Compatibility of Streamers
  {
    std::cout << "*  Starting  Backward IO compatibility - S T R E S S             *" <<std::endl;
    std::cout << "******************************************************************" <<std::endl;
    stress_backward_io();
    std::cout << "******************************************************************" <<std::endl;
  }

  gBenchmark->Stop("stressLinear");

  //Print table with results
  Bool_t UNIX = strcmp(gSystem->GetName(), "Unix") == 0;
  printf("******************************************************************\n");
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
    const Char_t *os = gSystem->Getenv("OS");
    if (!os) printf("*  SYS: Windows 95\n");
    else     printf("*  SYS: %s %s \n",os,gSystem->Getenv("PROCESSOR_IDENTIFIER"));
  }

  printf("******************************************************************\n");
  gBenchmark->Print("stressLinear");
  const Int_t nr = 7;
  const Double_t x_b12[] = { 10.,   30.,   50.,   100.,  300.,  500.,    700.};
  const Double_t y_b12[] = {10.74, 15.72, 20.00, 35.79, 98.77, 415.34, 1390.33};

  TGraph gr(nr,x_b12,y_b12);
  Double_t ct = gBenchmark->GetCpuTime("stressLinear");
  const Double_t rootmarks = 600*gr.Eval(maxSize)/ct;

  printf("******************************************************************\n");
  printf("*  ROOTMARKS =%6.1f   *  Root%-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),
         gROOT->GetVersionDate(),gROOT->GetVersionTime());
  printf("******************************************************************\n");
}

////////////////////////////////////////////////////////////////////////////////
/// Print test program number and its title

void StatusPrint(Int_t id,const TString &title,Bool_t status)
{
  const Int_t kMAX = 65;
  Char_t number[4];
  snprintf(number,4,"%2d",id);
  TString header = TString("Test ")+number+" : "+title;
  const Int_t nch = header.Length();
  for (Int_t i = nch; i < kMAX; i++) header += '.';
  std::cout << header << (status ? "OK" : "FAILED") << std::endl;
}

//------------------------------------------------------------------------
//          Test allocation functions and compatibility check
//
void mstress_allocation(Int_t msize)
{
  if (gVerbose)
    std::cout << "\n\n---> Test allocation and compatibility check" << std::endl;

  Int_t i,j;
  Bool_t ok = kTRUE;

  TMatrixD m1(4,msize);
  for (i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++)
    for (j = m1.GetColLwb(); j <= m1.GetColUpb(); j++)
      m1(i,j) = TMath::Pi()*i+TMath::E()*j;

  TMatrixD m2(0,3,0,msize-1);
  TMatrixD m3(1,4,0,msize-1);
  TMatrixD m4(m1);

  if (gVerbose) {
    std::cout << "\nStatus information reported for matrix m3:" << std::endl;
    std::cout << "  Row lower bound ... " << m3.GetRowLwb() << std::endl;
    std::cout << "  Row upper bound ... " << m3.GetRowUpb() << std::endl;
    std::cout << "  Col lower bound ... " << m3.GetColLwb() << std::endl;
    std::cout << "  Col upper bound ... " << m3.GetColUpb() << std::endl;
    std::cout << "  No. rows ..........." << m3.GetNrows()  << std::endl;
    std::cout << "  No. cols ..........." << m3.GetNcols()  << std::endl;
    std::cout << "  No. of elements ...." << m3.GetNoElements() << std::endl;
  }

  if (gVerbose)
    std::cout << "\nCheck matrices 1 & 2 for compatibility" << std::endl;
  ok &= AreCompatible(m1,m2,gVerbose);

  if (gVerbose)
    std::cout << "Check matrices 1 & 4 for compatibility" << std::endl;
  ok &= AreCompatible(m1,m4,gVerbose);

  if (gVerbose)
    std::cout << "m2 has to be compatible with m3 after resizing to m3" << std::endl;
  m2.ResizeTo(m3);
  ok &= AreCompatible(m2,m3,gVerbose);

  TMatrixD m5(m1.GetNrows()+1,m1.GetNcols()+5);
  for (i = m5.GetRowLwb(); i <= m5.GetRowUpb(); i++)
    for (j = m5.GetColLwb(); j <= m5.GetColUpb(); j++)
      m5(i,j) = TMath::Pi()*i+TMath::E()*j;

  if (gVerbose)
    std::cout << "m1 has to be compatible with m5 after resizing to m5" << std::endl;
  m1.ResizeTo(m5.GetNrows(),m5.GetNcols());
  ok &= AreCompatible(m1,m5,gVerbose);

  if (gVerbose)
    std::cout << "m1 has to be equal to m4 after stretching and shrinking" << std::endl;
  m1.ResizeTo(m4.GetNrows(),m4.GetNcols());
  ok &= VerifyMatrixIdentity(m1,m4,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "m5 has to be equal to m1 after shrinking" << std::endl;
  m5.ResizeTo(m1.GetNrows(),m1.GetNcols());
  ok &= VerifyMatrixIdentity(m1,m5,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "stretching and shrinking for small matrices (stack)" << std::endl;
  if (gVerbose)
    std::cout << "m8 has to be equal to m7 after stretching and shrinking" << std::endl;
  TMatrixD m6(4,4);
  for (i = m6.GetRowLwb(); i <= m6.GetRowUpb(); i++)
    for (j = m6.GetColLwb(); j <= m6.GetColUpb(); j++)
      m6(i,j) = TMath::Pi()*i+TMath::E()*j;
  TMatrixD m8(3,3);
  for (i = m8.GetRowLwb(); i <= m8.GetRowUpb(); i++)
    for (j = m8.GetColLwb(); j <= m8.GetColUpb(); j++)
      m8(i,j) = TMath::Pi()*i+TMath::E()*j;
  TMatrixD m7(m8);

  m8.ResizeTo(4,4);
  m8.ResizeTo(3,3);
  ok &= VerifyMatrixIdentity(m7,m8,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "m6 has to be equal to m8 after shrinking" << std::endl;
  m6.ResizeTo(3,3);
  ok &= VerifyMatrixIdentity(m6,m8,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(1,"Allocation, Resizing",ok);
}

class FillMatrix : public TElementPosActionD {
   Int_t no_elems,no_cols;
   void Operation(Double_t &element) const
      { element = 4*TMath::Pi()/no_elems * (fI*no_cols+fJ); }
public:
   FillMatrix() {}
   FillMatrix(const TMatrixD &m) :
         no_elems(m.GetNoElements()),no_cols(m.GetNcols()) { }
};

//
//------------------------------------------------------------------------
//          Test Filling of matrix
//
void mstress_matrix_fill(Int_t rsize,Int_t csize)
{
  if (gVerbose)
    std::cout << "\n\n---> Test different matrix filling methods\n" << std::endl;

  Bool_t ok = kTRUE;
  if (gVerbose)
    std::cout << "Creating m  with Apply function..." << std::endl;
  TMatrixD m(-1,rsize-2,1,csize);
#ifndef __CINT__
  FillMatrix f(m);
  m.Apply(f);
#else
  for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
    for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
      m(i,j) = 4*TMath::Pi()/m.GetNoElements() * (i*m.GetNcols()+j);
#endif

  {
    if (gVerbose)
      std::cout << "Check identity between m and matrix filled through (i,j)" << std::endl;

    TMatrixD m_overload1(-1,rsize-2,1,csize);
    TMatrixD m_overload2(-1,rsize-2,1,csize);

    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
    {
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
      {
        const Double_t val = 4*TMath::Pi()/rsize/csize*(i*csize+j);
        m_overload1(i,j)  = val;
        m_overload2[i][j] = val;
      }
    }

    ok &= VerifyMatrixIdentity(m,m_overload1,gVerbose,EPSILON);
    if (gVerbose)
      std::cout << "Check identity between m and matrix filled through [i][j]" << std::endl;
    ok &= VerifyMatrixIdentity(m,m_overload2,gVerbose,EPSILON);
    if (gVerbose)
      std::cout << "Check identity between matrix filled through [i][j] and (i,j)" << std::endl;
    ok &= VerifyMatrixIdentity(m_overload1,m_overload2,gVerbose,EPSILON);
  }

  {
    TArrayD a_fortran(rsize*csize);
    TArrayD a_c      (rsize*csize);
    for (Int_t i = 0; i < rsize; i++)
    {
      for (Int_t j = 0; j < csize; j++)
      {
        a_c[i*csize+j]       = 4*TMath::Pi()/rsize/csize*((i-1)*csize+j+1);
        a_fortran[i+rsize*j] = a_c[i*csize+j];
      }
    }

    if (gVerbose)
      std::cout << "Creating m_fortran by filling with fortran stored matrix" << std::endl;
    TMatrixD m_fortran(-1,rsize-2,1,csize,a_fortran.GetArray(),"F");
    if (gVerbose)
      std::cout << "Check identity between m and m_fortran" << std::endl;
    ok &= VerifyMatrixIdentity(m,m_fortran,gVerbose,EPSILON);

    if (gVerbose)
      std::cout << "Creating m_c by filling with c stored matrix" << std::endl;
    TMatrixD m_c(-1,rsize-2,1,csize,a_c.GetArray());
    if (gVerbose)
      std::cout << "Check identity between m and m_c" << std::endl;
    ok &= VerifyMatrixIdentity(m,m_c,gVerbose,EPSILON);
  }

  {
    if (gVerbose)
      std::cout << "Check insertion/extraction of sub-matrices" << std::endl;
    {
      TMatrixD m_sub1 = m;
      m_sub1.ResizeTo(0,rsize-2,2,csize);
      TMatrixD m_sub2 = m.GetSub(0,rsize-2,2,csize,"");
      ok &= VerifyMatrixIdentity(m_sub1,m_sub2,gVerbose,EPSILON);
    }

    {
      TMatrixD m2(-1,rsize-2,1,csize);
      TMatrixD m_part1 = m.GetSub(0,rsize-2,2,csize,"");
      TMatrixD m_part2 = m.GetSub(0,rsize-2,1,1,"");
      TMatrixD m_part3 = m.GetSub(-1,-1,2,csize,"");
      TMatrixD m_part4 = m.GetSub(-1,-1,1,1,"");
      m2.SetSub(0,2,m_part1);
      m2.SetSub(0,1,m_part2);
      m2.SetSub(-1,2,m_part3);
      m2.SetSub(-1,1,m_part4);
      ok &= VerifyMatrixIdentity(m,m2,gVerbose,EPSILON);
    }

    {
      TMatrixD m2(-1,rsize-2,1,csize);
      TMatrixD m_part1 = m.GetSub(0,rsize-2,2,csize,"S");
      TMatrixD m_part2 = m.GetSub(0,rsize-2,1,1,"S");
      TMatrixD m_part3 = m.GetSub(-1,-1,2,csize,"S");
      TMatrixD m_part4 = m.GetSub(-1,-1,1,1,"S");
      m2.SetSub(0,2,m_part1);
      m2.SetSub(0,1,m_part2);
      m2.SetSub(-1,2,m_part3);
      m2.SetSub(-1,1,m_part4);
      ok &= VerifyMatrixIdentity(m,m2,gVerbose,EPSILON);
    }
  }

  {
    if (gVerbose)
      std::cout << "Check sub-matrix views" << std::endl;
    {
      TMatrixD m3(-1,rsize-2,1,csize);
      TMatrixDSub(m3,0,rsize-2,2,csize) = TMatrixDSub(m,0,rsize-2,2,csize);
      TMatrixDSub(m3,0,rsize-2,1,1)     = TMatrixDSub(m,0,rsize-2,1,1);
      TMatrixDSub(m3,-1,-1,2,csize)     = TMatrixDSub(m,-1,-1,2,csize);
      TMatrixDSub(m3,-1,-1,1,1)         = TMatrixDSub(m,-1,-1,1,1);
      ok &= VerifyMatrixIdentity(m,m3,gVerbose,EPSILON);

      TMatrixD unit(3,3);
      TMatrixDSub(m3,1,3,1,3)  = unit.UnitMatrix();
      TMatrixDSub(m3,1,3,1,3) *= m.GetSub(1,3,1,3);
      ok &= VerifyMatrixIdentity(m,m3,gVerbose,EPSILON);

      TMatrixDSub(m3,0,rsize-2,2,csize) = 1.0;
      TMatrixDSub(m3,0,rsize-2,1,1)     = 1.0;
      TMatrixDSub(m3,-1,-1,2,csize)     = 1.0;
      TMatrixDSub(m3,-1,-1,1,1)         = 1.0;
      ok &= (m3 == 1.0);

    }
  }

  {
    if (gVerbose)
      std::cout << "Check array Use" << std::endl;
    {
      TMatrixD *m1a = new TMatrixD(m);
      TMatrixD *m2a = new TMatrixD();
      m2a->Use(m1a->GetRowLwb(),m1a->GetRowUpb(),m1a->GetColLwb(),m1a->GetColUpb(),m1a->GetMatrixArray());
      ok &= VerifyMatrixIdentity(m,*m2a,gVerbose,EPSILON);
      m2a->Sqr();
      TMatrixD m4 = m; m4.Sqr();
      ok &= VerifyMatrixIdentity(m4,*m1a,gVerbose,EPSILON);
      delete m1a;
      delete m2a;
    }
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(2,"Filling, Inserting, Using",ok);
}

//
//------------------------------------------------------------------------
//                Test uniform element operations
//
typedef  Double_t (*dfunc)(Double_t);
class ApplyFunction : public TElementActionD {
   dfunc fFunc;
   void Operation(Double_t &element) const { element = fFunc(Double_t(element)); }
public:
   ApplyFunction(dfunc func) : fFunc(func) { }
};

void mstress_element_op(Int_t rsize,Int_t csize)
{
  Bool_t ok = kTRUE;
  const Double_t pattern = 8.625;

  TMatrixD m(-1,rsize-2,1,csize);

  if (gVerbose)
    std::cout << "\nWriting zeros to m..." << std::endl;
  {
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for(Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        m(i,j) = 0;
    ok &= VerifyMatrixValue(m,0.,gVerbose,EPSILON);
  }

  if (gVerbose)
    std::cout << "Creating zero m1 ..." << std::endl;
  TMatrixD m1(TMatrixD::kZero, m);
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "Comparing m1 with 0 ..." << std::endl;
  R__ASSERT(m1 == 0);
  R__ASSERT(!(m1 != 0));

  if (gVerbose)
    std::cout << "Writing a pattern " << pattern << " by assigning to m(i,j)..." << std::endl;
  {
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        m(i,j) = pattern;
    ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  }

  if (gVerbose)
    std::cout << "Writing the pattern by assigning to m1 as a whole ..."  << std::endl;
  m1 = pattern;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "Comparing m and m1 ..." << std::endl;
  R__ASSERT(m == m1);
  if (gVerbose)
    std::cout << "Comparing (m=0) and m1 ..." << std::endl;
  R__ASSERT(!(m.Zero() == m1));

  if (gVerbose)
    std::cout << "Clearing m1 ..." << std::endl;
  m1.Zero();
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nClear m and add the pattern" << std::endl;
  m.Zero();
  m += pattern;
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "   add the doubled pattern with the negative sign" << std::endl;
  m += -2*pattern;
  ok &= VerifyMatrixValue(m,-pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "   subtract the trippled pattern with the negative sign" << std::endl;
  m -= -3*pattern;
  ok &= VerifyMatrixValue(m,2*pattern,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nVerify comparison operations when all elems are the same" << std::endl;
  m = pattern;
  R__ASSERT( m == pattern && !(m != pattern) );
  R__ASSERT( m > 0 && m >= pattern && m <= pattern );
  R__ASSERT( m > -pattern && m >= -pattern );
  R__ASSERT( m <= pattern && !(m < pattern) );
  m -= 2*pattern;
  R__ASSERT( m  < -pattern/2 && m <= -pattern/2 );
  R__ASSERT( m  >= -pattern && !(m > -pattern) );

  if (gVerbose)
    std::cout << "\nVerify comparison operations when not all elems are the same" << std::endl;
  m = pattern; m(m.GetRowUpb(),m.GetColUpb()) = pattern-1;
  R__ASSERT( !(m == pattern) && !(m != pattern) );
  R__ASSERT( m != 0 );                   // none of elements are 0
  R__ASSERT( !(m >= pattern) && m <= pattern && !(m<pattern) );
  R__ASSERT( !(m <= pattern-1) && m >= pattern-1 && !(m>pattern-1) );

  if (gVerbose)
    std::cout << "\nAssign 2*pattern to m by repeating additions" << std::endl;
  m = 0; m += pattern; m += pattern;
  if (gVerbose)
    std::cout << "Assign 2*pattern to m1 by multiplying by two " << std::endl;
  m1 = pattern; m1 *= 2;
  ok &= VerifyMatrixValue(m1,2*pattern,gVerbose,EPSILON);
  R__ASSERT( m == m1 );
  if (gVerbose)
    std::cout << "Multiply m1 by one half returning it to the 1*pattern" << std::endl;
  m1 *= 1/2.;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nAssign -pattern to m and m1" << std::endl;
  m.Zero(); m -= pattern; m1 = -pattern;
  ok &= VerifyMatrixValue(m,-pattern,gVerbose,EPSILON);
  R__ASSERT( m == m1 );
  if (gVerbose)
    std::cout << "m = sqrt(sqr(m)); m1 = abs(m1); Now m and m1 have to be the same" << std::endl;
  m.Sqr();
  ok &= VerifyMatrixValue(m,pattern*pattern,gVerbose,EPSILON);
  m.Sqrt();
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  m1.Abs();
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  ok &= VerifyMatrixIdentity(m1,m,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nCheck out to see that sin^2(x) + cos^2(x) = 1" << std::endl;
  {
#ifndef __CINT__
    FillMatrix f(m);
    m.Apply(f);
#else
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        m(i,j) = 4*TMath::Pi()/m.GetNoElements() * (i*m.GetNcols()+j);
#endif
  }
  m1 = m;
  {
#ifndef __CINT__
    ApplyFunction s(&TMath::Sin);
    ApplyFunction c(&TMath::Cos);
    m.Apply(s);
    m1.Apply(c);
#else
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++) {
        m(i,j)  = TMath::Sin(m(i,j));
        m1(i,j) = TMath::Cos(m1(i,j));
      }
    }
#endif
  }
  m.Sqr();
  m1.Sqr();
  m += m1;
  ok &= VerifyMatrixValue(m,1.,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(3,"Uniform matrix operations",ok);
}

//
//------------------------------------------------------------------------
//        Test binary matrix element-by-element operations
//
void mstress_binary_ebe_op(Int_t rsize, Int_t csize)
{
  if (gVerbose)
    std::cout << "\n---> Test Binary Matrix element-by-element operations" << std::endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 4.25;

  TMatrixD m(2,rsize+1,0,csize-1);
  TMatrixD m1(TMatrixD::kZero,m);
  TMatrixD mp(TMatrixD::kZero,m);

  {
    for (Int_t i = mp.GetRowLwb(); i <= mp.GetRowUpb(); i++)
      for (Int_t j = mp.GetColLwb(); j <= mp.GetColUpb(); j++)
        mp(i,j) = (i-m.GetNrows()/2.)*j*pattern;
  }

  if (gVerbose)
    std::cout << "\nVerify assignment of a matrix to the matrix" << std::endl;
  m = pattern;
  m1.Zero();
  m1 = m;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  R__ASSERT( m1 == m );

  if (gVerbose)
    std::cout << "\nAdding the matrix to itself, uniform pattern " << pattern << std::endl;
  m.Zero(); m = pattern;
  m1 = m; m1 += m1;
  ok &= VerifyMatrixValue(m1,2*pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "  subtracting two matrices ..." << std::endl;
  m1 -= m;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "  subtracting the matrix from itself" << std::endl;
  // Hiding the self-subtraction from the compiler, causing warning by -Wself-assign-overloaded in clang 7.0
  m1 -= static_cast<TMatrixD&>(m1);
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "  adding two matrices together" << std::endl;
  m1 += m;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);

  if (gVerbose) {
    std::cout << "\nArithmetic operations on matrices with not the same elements" << std::endl;
    std::cout << "   adding mp to the zero matrix..." << std::endl;
  }
  m.Zero(); m += mp;
  ok &= VerifyMatrixIdentity(m,mp,gVerbose,EPSILON);
  m1 = m;
  if (gVerbose)
    std::cout << "   making m = 3*mp and m1 = 3*mp, via add() and succesive mult" << std::endl;
  Add(m,2.,mp);
  m1 += m1; m1 += mp;
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "   clear both m and m1, by subtracting from itself and via add()" << std::endl;
  // Hiding the self-subtraction from the compiler, causing warning by -Wself-assign-overloaded in clang 7.0
  m1 -= static_cast<TMatrixD&>(m1);
  Add(m,-3.,mp);
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);

  if (gVerbose) {
    std::cout << "\nTesting element-by-element multiplications and divisions" << std::endl;
    std::cout << "   squaring each element with sqr() and via multiplication" << std::endl;
  }
  m = mp; m1 = mp;
  m.Sqr();
  ElementMult(m1,m1);
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "   compare (m = pattern^2)/pattern with pattern" << std::endl;
  m = pattern; m1 = pattern;
  m.Sqr();
  ElementDiv(m,m1);
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  if (gVerbose)
    Compare(m1,m);

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(4,"Binary Matrix element-by-element operations",ok);
}

//
//------------------------------------------------------------------------
//              Verify matrix transposition
//
void mstress_transposition(Int_t msize)
{
  if (gVerbose) {
    std::cout << "\n---> Verify matrix transpose "
            "for matrices of a characteristic size " << msize << std::endl;
  }

  Bool_t ok = kTRUE;
  {
    if (gVerbose)
      std::cout << "\nCheck to see that a square UnitMatrix stays the same";
    TMatrixD m(msize,msize);
    m.UnitMatrix();
    TMatrixD mt(TMatrixD::kTransposed,m);
    ok &= ( m == mt ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      std::cout << "\nTest a non-square UnitMatrix";
    TMatrixD m(msize,msize+1);
    m.UnitMatrix();
    TMatrixD mt(TMatrixD::kTransposed,m);
    R__ASSERT(m.GetNrows() == mt.GetNcols() && m.GetNcols() == mt.GetNrows() );
    for (Int_t i = m.GetRowLwb(); i <= TMath::Min(m.GetRowUpb(),m.GetColUpb()); i++)
      for (Int_t j = m.GetColLwb(); j <= TMath::Min(m.GetRowUpb(),m.GetColUpb()); j++)
        ok &= ( m(i,j) == mt(i,j) ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      std::cout << "\nCheck to see that a symmetric (Hilbert)Matrix stays the same";
    TMatrixD m = THilbertMatrixD(msize,msize);
    TMatrixD mt(TMatrixD::kTransposed,m);
    ok &= ( m == mt ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      std::cout << "\nCheck transposing a non-symmetric matrix";
    TMatrixD m = THilbertMatrixD(msize+1,msize);
    m(1,2) = TMath::Pi();
    TMatrixD mt(TMatrixD::kTransposed,m);
    R__ASSERT(m.GetNrows() == mt.GetNcols() && m.GetNcols() == mt.GetNrows());
    R__ASSERT(mt(2,1)  == (Double_t)TMath::Pi() && mt(1,2)  != (Double_t)TMath::Pi());
    R__ASSERT(mt[2][1] == (Double_t)TMath::Pi() && mt[1][2] != (Double_t)TMath::Pi());

    if (gVerbose)
      std::cout << "\nCheck double transposing a non-symmetric matrix" << std::endl;
    TMatrixD mtt(TMatrixD::kTransposed,mt);
    ok &= ( m == mtt ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(5,"Matrix transposition",ok);
}

//
//------------------------------------------------------------------------
//           Test special matrix creation
//
class MakeHilbert : public TElementPosActionD {
  void Operation(Double_t &element) const { element = 1./(fI+fJ+1); }
public:
  MakeHilbert() { }
};

#if !defined (__CINT__) || defined (__MAKECINT__)
class TestUnit : public TElementPosActionD {
  mutable Int_t fIsUnit;
  void Operation(Double_t &element) const
      { if (fIsUnit)
          fIsUnit = ((fI==fJ) ? (element == 1.0) : (element == 0)); }
public:
  TestUnit() : fIsUnit(0==0) { }
  Int_t is_indeed_unit() const { return fIsUnit; }
};
#else
  Bool_t is_indeed_unit(TMatrixD &m) {
    Bool_t isUnit = kTRUE;
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++) {
        if (isUnit)
          isUnit = ((i==j) ? (m(i,j) == 1.0) : (m(i,j) == 0));
      }
    return isUnit;
  }
#endif

void mstress_special_creation(Int_t dim)
{
  if (gVerbose)
    std::cout << "\n---> Check creating some special matrices of dimension " << dim << std::endl;

  Bool_t ok = kTRUE;
  {
    if (gVerbose)
      std::cout << "\ntest creating Hilbert matrices" << std::endl;
    TMatrixD m = THilbertMatrixD(dim+1,dim);
    TMatrixD m1(TMatrixD::kZero,m);
    ok &= ( !(m == m1) ) ? kTRUE : kFALSE;
    ok &= ( m != 0 ) ? kTRUE : kFALSE;
#ifndef __CINT__
    MakeHilbert mh;
    m1.Apply(mh);
#else
    for (Int_t i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++)
      for (Int_t j = m1.GetColLwb(); j <= m1.GetColUpb(); j++)
        m1(i,j) = 1./(i+j+1);
#endif
    ok &= ( m1 != 0 ) ? kTRUE : kFALSE;
    ok &= ( m == m1 ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      std::cout << "test creating zero matrix and copy constructor" << std::endl;
    TMatrixD m = THilbertMatrixD(dim,dim+1);
    ok &= ( m != 0 ) ? kTRUE : kFALSE;
    TMatrixD m1(m);               // Applying the copy constructor
    ok &= ( m1 == m ) ? kTRUE : kFALSE;
    TMatrixD m2(TMatrixD::kZero,m);
    ok &= ( m2 == 0 ) ? kTRUE : kFALSE;
    ok &= ( m != 0 ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      std::cout << "test creating unit matrices" << std::endl;
    TMatrixD m(dim,dim);
#ifndef __CINT__
    {
      TestUnit test_unit;
      m.Apply(test_unit);
      ok &= ( !test_unit.is_indeed_unit() ) ? kTRUE : kFALSE;
    }
#else
    ok &= ( !is_indeed_unit(m) ) ? kTRUE : kFALSE;
#endif
    m.UnitMatrix();
#ifndef __CINT__
    {
      TestUnit test_unit;
       m.Apply(test_unit);
       ok &= ( test_unit.is_indeed_unit() ) ? kTRUE : kFALSE;
    }
#else
    ok &= ( is_indeed_unit(m) ) ? kTRUE : kFALSE;
#endif
    m.ResizeTo(dim-1,dim);
    TMatrixD m2(TMatrixD::kUnit,m);
#ifndef __CINT__
    {
      TestUnit test_unit;
      m2.Apply(test_unit);
      ok &= ( test_unit.is_indeed_unit() ) ? kTRUE : kFALSE;
    }
#else
    ok &= ( is_indeed_unit(m2) ) ? kTRUE : kFALSE;
#endif
    m.ResizeTo(dim,dim-2);
    m.UnitMatrix();
#ifndef __CINT__
    {
      TestUnit test_unit;
      m.Apply(test_unit);
      ok &= ( test_unit.is_indeed_unit() ) ? kTRUE : kFALSE;
    }
#else
    ok &= ( is_indeed_unit(m) ) ? kTRUE : kFALSE;
#endif
  }

  {
    if (gVerbose)
      std::cout << "check to see that Haar matrix has *exactly* orthogonal columns" << std::endl;
    Int_t j;
    const Int_t order = 5;
    const TMatrixD haar = THaarMatrixD(order);
    ok &= ( haar.GetNrows() == (1<<order) &&
               haar.GetNrows() == haar.GetNcols() ) ? kTRUE : kFALSE;
    TVectorD colj(1<<order);
    TVectorD coll(1<<order);
    for (j = haar.GetColLwb(); j <= haar.GetColUpb(); j++) {
      colj = TMatrixDColumn_const(haar,j);
      ok &= (TMath::Abs(colj*colj-1.0) <= 1.0e-15 ) ? kTRUE : kFALSE;
      for (Int_t l = j+1; l <= haar.GetColUpb(); l++) {
        coll = TMatrixDColumn_const(haar,l);
        const Double_t val = colj*coll;
        ok &= ( TMath::Abs(val) <= 1.0e-15 ) ? kTRUE : kFALSE;
      }
    }

    if (gVerbose)
      std::cout << "make Haar (sub)matrix and test it *is* a submatrix" << std::endl;
    const Int_t no_sub_cols = (1<<order) - 3;
    const TMatrixD haar_sub = THaarMatrixD(order,no_sub_cols);
    ok &= ( haar_sub.GetNrows() == (1<<order) &&
               haar_sub.GetNcols() == no_sub_cols ) ? kTRUE : kFALSE;
    for (j = haar_sub.GetColLwb(); j <= haar_sub.GetColUpb(); j++) {
      colj = TMatrixDColumn_const(haar,j);
      coll = TMatrixDColumn_const(haar_sub,j);
      ok &= VerifyVectorIdentity(colj,coll,gVerbose);
    }
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(6,"Haar/Hilbert Matrix",ok);
}

#ifndef __CINT__
//
//------------------------------------------------------------------------
//           Test matrix promises
//
class hilbert_matrix_promise : public TMatrixDLazy {
  void FillIn(TMatrixD &m) const { m = THilbertMatrixD(m.GetRowLwb(),m.GetRowUpb(),
                                                   m.GetColLwb(),m.GetColUpb()); }

public:
  hilbert_matrix_promise(Int_t nrows,Int_t ncols)
     : TMatrixDLazy(nrows,ncols) {}
  hilbert_matrix_promise(Int_t row_lwb,Int_t row_upb,
                         Int_t col_lwb,Int_t col_upb)
     : TMatrixDLazy(row_lwb,row_upb,col_lwb,col_upb) { }
};

void mstress_matrix_promises(Int_t dim)
{
  if (gVerbose)
    std::cout << "\n---> Check making/forcing promises, (lazy)matrices of dimension " << dim << std::endl;

  Bool_t ok = kTRUE;
  {
    if (gVerbose)
      std::cout << "\nmake a promise and force it by a constructor" << std::endl;
    TMatrixD m  = hilbert_matrix_promise(dim,dim+1);
    TMatrixD m1 = THilbertMatrixD(dim,dim+1);
    ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  }

  {
    if (gVerbose)
      std::cout << "make a promise and force it by an assignment" << std::endl;
    TMatrixD m(-1,dim,0,dim);
    m = hilbert_matrix_promise(-1,dim,0,dim);
    TMatrixD m1 = THilbertMatrixD(-1,dim,0,dim);
    ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(7,"Matrix promises",ok);
}
#endif

//
//------------------------------------------------------------------------
//             Verify the norm calculation
//
void mstress_norms(Int_t rsize_req,Int_t csize_req)
{
  if (gVerbose)
    std::cout << "\n---> Verify norm calculations" << std::endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 10.25;

  Int_t rsize = rsize_req;
  if (rsize%2 != 0)  rsize--;
  Int_t csize = csize_req;
  if (csize%2 != 0)  csize--;
  if (rsize%2 == 1 || csize%2 == 1) {
    std::cout << "rsize: " << rsize <<std::endl;
    std::cout << "csize: " << csize <<std::endl;
    Fatal("mstress_norms","Sorry, size of the matrix to test must be even for this test\n");
  }

  TMatrixD m(rsize,csize);

  if (gVerbose)
    std::cout << "\nAssign " << pattern << " to all the elements and check norms" << std::endl;
  m = pattern;
  if (gVerbose)
    std::cout << "  1. (col) norm should be pattern*nrows" << std::endl;
  ok &= ( m.Norm1() == pattern*m.GetNrows() ) ? kTRUE : kFALSE;
  ok &= ( m.Norm1() == m.ColNorm() ) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "  Inf (row) norm should be pattern*ncols" << std::endl;
  ok &= ( m.NormInf() == pattern*m.GetNcols() ) ? kTRUE : kFALSE;
  ok &= ( m.NormInf() == m.RowNorm() ) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "  Square of the Eucl norm has got to be pattern^2 * no_elems" << std::endl;
  ok &= ( m.E2Norm() == (pattern*pattern)*m.GetNoElements() ) ? kTRUE : kFALSE;
  TMatrixD m1(TMatrixD::kZero,m);
  ok &= ( m.E2Norm() == E2Norm(m,m1) ) ? kTRUE : kFALSE;

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(8,"Matrix Norms",ok);
}

//
//------------------------------------------------------------------------
//              Verify the determinant evaluation
//
void mstress_determinant(Int_t msize)
{
  if (gVerbose)
    std::cout << "\n---> Verify determinant evaluation for a square matrix of size " << msize << std::endl;

  Bool_t ok = kTRUE;
  TMatrixD m(msize,msize);
  const Double_t pattern = 2.5;

  {
    if (gVerbose)
      std::cout << "\nCheck to see that the determinant of the unit matrix is one" <<std::endl;
    m.UnitMatrix();
    Double_t d1,d2;
    m.Determinant(d1,d2);
    const Double_t det = d1*TMath::Power(2.,d2);
    if (gVerbose) {
      std::cout << "det = " << det << " deviation= " << TMath::Abs(det-1);
      std::cout << ( (TMath::Abs(det-1) < EPSILON) ? " OK" : " too large") <<std::endl;
    }
    ok &= (TMath::Abs(det-1) < EPSILON) ? kTRUE : kFALSE;
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "\nCheck the determinant for the matrix with " << pattern << " at the diagonal" << std::endl;
    {
      for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
        for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
          m(i,j) = ( i==j ? pattern : 0 );
    }
    Double_t d1,d2;
    m.Determinant(d1,d2);
    const Double_t d1_abs = TMath::Abs(d1);
    const Double_t log_det1 = d2*TMath::Log(2.)+TMath::Log(d1_abs);
    const Double_t log_det2 = ((Double_t)m.GetNrows())*TMath::Log(pattern);
    if (gVerbose) {
      std::cout << "log_det1 = " << log_det1 << " deviation= " << TMath::Abs(log_det1-log_det2);
      std::cout << ( (TMath::Abs(log_det1-log_det2) < msize*EPSILON) ? " OK" : " too large") <<std::endl;
    }
    ok &= ( TMath::Abs(log_det1-log_det2) < msize*EPSILON  ) ? kTRUE : kFALSE;
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "\nCheck the determinant of the transposed matrix" << std::endl;
    m.UnitMatrix();
    m(1,2) = pattern;
    TMatrixD m_tran(TMatrixD::kTransposed,m);
    ok &= ( !(m == m_tran) ) ? kTRUE : kFALSE;
    const Double_t det1 = m.Determinant();
    const Double_t det2 = m_tran.Determinant();
    if (gVerbose) {
      std::cout << "det1 = " << det1 << " deviation= " << TMath::Abs(det1-det2);
      std::cout << ( (TMath::Abs(det1-det2) < msize*EPSILON) ? " OK" : " too large") <<std::endl;
    }
    ok &= ( TMath::Abs(det1-det2) < msize*EPSILON ) ? kTRUE : kFALSE;
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "\nswap two rows/cols of a matrix through method 1 and watch det's sign" << std::endl;
    m.UnitMatrix();
    TMatrixDRow(m,3).Assign(pattern);
    Double_t d1,d2;
    m.Determinant(d1,d2);
    TMatrixDRow row1(m,1);
    TVectorD vrow1(m.GetRowLwb(),m.GetRowUpb()); vrow1 = row1;
    TVectorD vrow3(m.GetRowLwb(),m.GetRowUpb()); vrow3 = TMatrixDRow(m,3);
    row1 = vrow3; TMatrixDRow(m,3) = vrow1;
    Double_t d1_p,d2_p;
    m.Determinant(d1_p,d2_p);
    if (gVerbose) {
      std::cout << "d1 = " << d1 << " deviation= " << TMath::Abs(d1+d1_p);
      std::cout << ( ( d1 == -d1_p ) ? " OK" : " too large") <<std::endl;
    }
    ok &= ( d1 == -d1_p ) ? kTRUE : kFALSE;
    TMatrixDColumn col2(m,2);
    TVectorD vcol2(m.GetRowLwb(),m.GetRowUpb()); vcol2 = col2;
    TVectorD vcol4(m.GetRowLwb(),m.GetRowUpb()); vcol4 = TMatrixDColumn(m,4);
    col2 = vcol4; TMatrixDColumn(m,4) = vcol2;
    m.Determinant(d1_p,d2_p);
    if (gVerbose) {
      std::cout << "d1 = " << d1 << " deviation= " << TMath::Abs(d1-d1_p);
      std::cout << ( ( d1 == d1_p ) ? " OK" : " too large") <<std::endl;
    }
    ok &= ( d1 == d1_p ) ? kTRUE : kFALSE;
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "\nswap two rows/cols of a matrix through method 2 and watch det's sign" << std::endl;
    m.UnitMatrix();
    TMatrixDRow(m,3).Assign(pattern);
    Double_t d1,d2;
    m.Determinant(d1,d2);

    TMatrixD m_save( m);
    TMatrixDRow(m,1) = TMatrixDRow(m_save,3);
    TMatrixDRow(m,3) = TMatrixDRow(m_save,1);
    Double_t d1_p,d2_p;
    m.Determinant(d1_p,d2_p);
    if (gVerbose) {
      std::cout << "d1 = " << d1 << " deviation= " << TMath::Abs(d1+d1_p);
      std::cout << ( ( d1 == -d1_p ) ? " OK" : " too large") <<std::endl;
    }
    ok &= ( d1 == -d1_p ) ? kTRUE : kFALSE;

    m_save = m;
    TMatrixDColumn(m,2) = TMatrixDColumn(m_save,4);
    TMatrixDColumn(m,4) = TMatrixDColumn(m_save,2);
    m.Determinant(d1_p,d2_p);
    if (gVerbose) {
      std::cout << "d1 = " << d1 << " deviation= " << TMath::Abs(d1-d1_p);
      std::cout << ( ( d1 == d1_p ) ? " OK" : " too large") <<std::endl;
    }
    ok &= ( d1 == d1_p ) ? kTRUE : kFALSE;
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "\nCheck the determinant for the matrix with " << pattern << " at the anti-diagonal" << std::endl;
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        m(i,j) = ( i==(m.GetColUpb()+m.GetColLwb()-j) ? pattern : 0 );

    Double_t d1,d2;
    m.Determinant(d1,d2);
    const Double_t d1_abs = TMath::Abs(d1);
    const Double_t log_det1 = d2*TMath::Log(2.)+TMath::Log(d1_abs);
    const Double_t log_det2 = ((Double_t)m.GetNrows())*TMath::Log(pattern);
    const Double_t sign     = ( m.GetNrows()*(m.GetNrows()-1)/2 & 1 ? -1 : 1 );
    if (gVerbose) {
      std::cout << "log_det1 = " << log_det1 << " deviation= " << TMath::Abs(log_det1-log_det2);
      std::cout << ( (TMath::Abs(log_det1-log_det2) < msize*EPSILON) ? " OK" : " too large") <<std::endl;
      std::cout << ( sign * d1 > 0. ? " OK sign" : " wrong sign") <<std::endl;
    }
    ok &= ( TMath::Abs(log_det1-log_det2) < msize*EPSILON  ) ? kTRUE : kFALSE;
    ok &= ( sign * d1 > 0. ) ? kTRUE : kFALSE;
  }

  // Next test disabled because it produces (of course) a Warning
  if (0) {
    if (gVerbose)
       std::cout << "\nCheck the determinant for the singular matrix"
                    "\n\tdefined as above with zero first row" << std::endl;
    m.Zero();
    {
      for (Int_t i = m.GetRowLwb()+1; i <= m.GetRowUpb(); i++)
        for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
          m(i,j) = ( i==(m.GetColUpb()+m.GetColLwb()-j) ? pattern : 0 );
    }
    std::cout << "\n\tdeterminant is " << m.Determinant();
    ok &= ( m.Determinant() == 0 ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    std::cout << "\nCheck out the determinant of the Hilbert matrix";
  TMatrixDSym H = THilbertMatrixDSym(3);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    std::cout << "\n    3x3 Hilbert matrix: exact determinant 1/2160 ";
    std::cout << "\n                              computed    1/"<< 1/H.Determinant();
  }

  H.ResizeTo(4,4);
  H = THilbertMatrixDSym(4);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    std::cout << "\n    4x4 Hilbert matrix: exact determinant 1/6048000 ";
    std::cout << "\n                              computed    1/"<< 1/H.Determinant();
  }

  H.ResizeTo(5,5);
  H = THilbertMatrixDSym(5);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    std::cout << "\n    5x5 Hilbert matrix: exact determinant 3.749295e-12";
    std::cout << "\n                              computed    "<< H.Determinant();
  }

  if (gVerbose) {
    TDecompQRH qrh(H);
    Double_t d1,d2;
    qrh.Det(d1,d2);
    std::cout  << "\n qrh det = " << d1*TMath::Power(2.0,d2) <<std::endl;
  }

  if (gVerbose) {
    TDecompChol chol(H);
    Double_t d1,d2;
    chol.Det(d1,d2);
    std::cout  << "\n chol det = " << d1*TMath::Power(2.0,d2) <<std::endl;
  }

  if (gVerbose) {
    TDecompSVD svd(H);
    Double_t d1,d2;
    svd.Det(d1,d2);
    std::cout  << "\n svd det = " << d1*TMath::Power(2.0,d2) <<std::endl;
  }

  H.ResizeTo(7,7);
  H = THilbertMatrixDSym(7);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    std::cout << "\n    7x7 Hilbert matrix: exact determinant 4.8358e-25";
    std::cout << "\n                              computed    "<< H.Determinant();
  }

  H.ResizeTo(9,9);
  H = THilbertMatrixDSym(9);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    std::cout << "\n    9x9 Hilbert matrix: exact determinant 9.72023e-43";
    std::cout << "\n                              computed    "<< H.Determinant();
  }

  H.ResizeTo(10,10);
  H = THilbertMatrixDSym(10);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    std::cout << "\n    10x10 Hilbert matrix: exact determinant 2.16418e-53";
    std::cout << "\n                              computed    "<< H.Determinant();
  }

  if (gVerbose)
  std::cout << "\nDone\n" << std::endl;

  StatusPrint(9,"Matrix Determinant",ok);
}

//
//------------------------------------------------------------------------
//               Verify matrix multiplications
//
void mstress_mm_multiplications()
{
  Bool_t ok = kTRUE;

  Int_t iloop = 0;
  Int_t nr    = 0;
  while (iloop <= gNrLoop) {
    const Int_t msize = gSizeA[iloop];
    const Double_t epsilon = EPSILON*msize/100;

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    Int_t i,j;
    if (msize <= 5) {
       iloop++;
       continue;  // some references to m(3,..)
    }

    if (verbose)
      std::cout << "\n---> Verify matrix multiplications "
              "for matrices of the characteristic size " << msize << std::endl;

    {
      if (verbose)
        std::cout << "\nTest inline multiplications of the UnitMatrix" << std::endl;
      TMatrixD m = THilbertMatrixD(-1,msize,-1,msize);
      TMatrixD u(TMatrixD::kUnit,m);
      m(3,1) = TMath::Pi();
      u *= m;
      ok &= VerifyMatrixIdentity(u,m,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Test inline multiplications by a DiagMat" << std::endl;
      TMatrixD m = THilbertMatrixD(msize+3,msize);
      m(1,3) = TMath::Pi();
      TVectorD v(msize);
      for (i = v.GetLwb(); i <= v.GetUpb(); i++)
        v(i) = 1+i;
      TMatrixD diag(msize,msize);
      TMatrixDDiag d(diag);
      d = v;
      TMatrixD eth = m;
      for (i = eth.GetRowLwb(); i <= eth.GetRowUpb(); i++)
        for (j = eth.GetColLwb(); j <= eth.GetColUpb(); j++)
          eth(i,j) *= v(j);
      m *= diag;
      ok &= VerifyMatrixIdentity(m,eth,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Test XPP = X where P is a permutation matrix" << std::endl;
      TMatrixD m = THilbertMatrixD(msize-1,msize);
      m(2,3) = TMath::Pi();
      TMatrixD eth = m;
      TMatrixD p(msize,msize);
      for (i = p.GetRowLwb(); i <= p.GetRowUpb(); i++)
        p(p.GetRowUpb()+p.GetRowLwb()-i,i) = 1;
      m *= p;
      m *= p;
      ok &= VerifyMatrixIdentity(m,eth,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Test general matrix multiplication through inline mult" << std::endl;
      TMatrixD m = THilbertMatrixD(msize-2,msize);
      m(3,3) = TMath::Pi();
      TMatrixD mt(TMatrixD::kTransposed,m);
      TMatrixD p = THilbertMatrixD(msize,msize);
      TMatrixDDiag(p) += 1;
      TMatrixD mp(m,TMatrixD::kMult,p);
      TMatrixD m1 = m;
      m *= p;
      ok &= VerifyMatrixIdentity(m,mp,verbose,epsilon);
      TMatrixD mp1(mt,TMatrixD::kTransposeMult,p);
      VerifyMatrixIdentity(m,mp1,verbose,epsilon);
      ok &= ( !(m1 == m) );
      TMatrixD mp2(TMatrixD::kZero,m1);
      ok &= ( mp2 == 0 );
      mp2.Mult(m1,p);
      ok &= VerifyMatrixIdentity(m,mp2,verbose,epsilon);

      if (verbose)
        std::cout << "Test XP=X*P  vs XP=X;XP*=P" << std::endl;
      TMatrixD mp3 = m1*p;
      ok &= VerifyMatrixIdentity(m,mp3,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Check n * m  == n * m_sym; m_sym * n == m * n; m_sym * m_sym == m * m" <<std::endl;

      const TMatrixD     n     = THilbertMatrixD(0,msize-1,0,msize-1);
      const TMatrixD     m     = n;
      const TMatrixDSym  m_sym = THilbertMatrixDSym(0,msize-1);

      const TMatrixD nm1 = n * m_sym;
      const TMatrixD nm2 = n * m;
      const TMatrixD mn1 = m_sym * n;
      const TMatrixD mn2 = m * n;
      const TMatrixD mm1 = m_sym * m_sym;
      const TMatrixD mm2 = m * m;

      ok &= VerifyMatrixIdentity(nm1,nm2,verbose,epsilon);
      ok &= VerifyMatrixIdentity(mn1,mn2,verbose,epsilon);
      ok &= VerifyMatrixIdentity(mm1,mm2,verbose,epsilon);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop++;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        std::cout << "General Matrix Multiplications failed for size " << msize << std::endl;
      break;
    }
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "Check to see UU' = U'U = E when U is the Haar matrix" << std::endl;
    const Int_t order = 5;
    const Int_t no_sub_cols = (1<<order)-5;
    TMatrixD haar_sub = THaarMatrixD(5,no_sub_cols);
    TMatrixD haar_sub_t(TMatrixD::kTransposed,haar_sub);
    TMatrixD hsths(haar_sub_t,TMatrixD::kMult,haar_sub);
    TMatrixD hsths1(TMatrixD::kZero,hsths); hsths1.Mult(haar_sub_t,haar_sub);
    TMatrixD hsths_eth(TMatrixD::kUnit,hsths);
    ok &= ( hsths.GetNrows() == no_sub_cols && hsths.GetNcols() == no_sub_cols );
    ok &= VerifyMatrixIdentity(hsths,hsths_eth,gVerbose,EPSILON);
    ok &= VerifyMatrixIdentity(hsths1,hsths_eth,gVerbose,EPSILON);
    TMatrixD haar = THaarMatrixD(order);
    TMatrixD unit(TMatrixD::kUnit,haar);
    TMatrixD haar_t(TMatrixD::kTransposed,haar);
    TMatrixD hth(haar,TMatrixD::kTransposeMult,haar);
    TMatrixD hht(haar,TMatrixD::kMult,haar_t);
    TMatrixD hht1 = haar; hht1 *= haar_t;
    TMatrixD hht2(TMatrixD::kZero,haar); hht2.Mult(haar,haar_t);
    ok &= VerifyMatrixIdentity(unit,hth,gVerbose,EPSILON);
    ok &= VerifyMatrixIdentity(unit,hht,gVerbose,EPSILON);
    ok &= VerifyMatrixIdentity(unit,hht1,gVerbose,EPSILON);
    ok &= VerifyMatrixIdentity(unit,hht2,gVerbose,EPSILON);
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(10,"General Matrix Multiplications",ok);
}

//
//------------------------------------------------------------------------
//               Verify symmetric matrix multiplications
//
void mstress_sym_mm_multiplications(Int_t msize)
{
  if (gVerbose)
    std::cout << "\n---> Verify symmetric matrix multiplications "
            "for matrices of the characteristic size " << msize << std::endl;

  Int_t i,j;
  Bool_t ok = kTRUE;

  const Double_t epsilon = EPSILON*msize/100;

  {
    if (gVerbose)
      std::cout << "\nTest inline multiplications of the UnitMatrix" << std::endl;
    TMatrixD m = THilbertMatrixD(-1,msize,-1,msize);
    TMatrixDSym m_sym(-1,msize,m.GetMatrixArray());
    TMatrixDSym u(TMatrixDSym::kUnit,m_sym);
    TMatrixD u2 = u * m_sym;
    ok &= VerifyMatrixIdentity(u2,m_sym,gVerbose,epsilon);
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "\nTest symmetric multiplications" << std::endl;
    {
      if (gVerbose)
        std::cout << "\n  Test m * m_sym == m_sym * m == m_sym * m_sym  multiplications" << std::endl;
      TMatrixD m = THilbertMatrixD(-1,msize,-1,msize);
      TMatrixDSym m_sym(-1,msize,m.GetMatrixArray());
      TMatrixD mm      = m * m;
      TMatrixD mm_sym1 = m_sym * m_sym;
      TMatrixD mm_sym2 = m     * m_sym;
      TMatrixD mm_sym3 = m_sym * m;
      ok &= VerifyMatrixIdentity(mm,mm_sym1,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mm,mm_sym2,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mm,mm_sym3,gVerbose,epsilon);
    }

    {
      if (gVerbose)
        std::cout << "\n  Test m^T * m_sym == m_sym^T * m == m_sym^T * m_sym  multiplications" << std::endl;
      TMatrixD m = THilbertMatrixD(-1,msize,-1,msize);
      TMatrixDSym m_sym(-1,msize,m.GetMatrixArray());
      TMatrixD mtm      = TMatrixD(m    ,TMatrixD::kTransposeMult,m);
      TMatrixD mtm_sym1 = TMatrixD(m_sym,TMatrixD::kTransposeMult,m_sym);
      TMatrixD mtm_sym2 = TMatrixD(m    ,TMatrixD::kTransposeMult,m_sym);
      TMatrixD mtm_sym3 = TMatrixD(m_sym,TMatrixD::kTransposeMult,m);
      ok &= VerifyMatrixIdentity(mtm,mtm_sym1,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mtm,mtm_sym2,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mtm,mtm_sym3,gVerbose,epsilon);
    }

    {
      if (gVerbose)
        std::cout << "\n  Test m * m_sym^T == m_sym * m^T == m_sym * m_sym^T  multiplications" << std::endl;
      TMatrixD m = THilbertMatrixD(-1,msize,-1,msize);
      TMatrixDSym m_sym(-1,msize,m.GetMatrixArray());
      TMatrixD mmt      = TMatrixD(m    ,TMatrixD::kMultTranspose,m);
      TMatrixD mmt_sym1 = TMatrixD(m_sym,TMatrixD::kMultTranspose,m_sym);
      TMatrixD mmt_sym2 = TMatrixD(m    ,TMatrixD::kMultTranspose,m_sym);
      TMatrixD mmt_sym3 = TMatrixD(m_sym,TMatrixD::kMultTranspose,m);
      ok &= VerifyMatrixIdentity(mmt,mmt_sym1,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mmt,mmt_sym2,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mmt,mmt_sym3,gVerbose,epsilon);
    }

    {
      if (gVerbose)
        std::cout << "\n  Test n * m_sym == n * m multiplications" << std::endl;
      TMatrixD n = THilbertMatrixD(-1,msize,-1,msize);
      TMatrixD m = n;
      n(1,3) = TMath::Pi();
      n(3,1) = TMath::Pi();
      TMatrixDSym m_sym(-1,msize,m.GetMatrixArray());
      TMatrixD nm1 = n * m_sym;
      TMatrixD nm2 = n * m;
      ok &= VerifyMatrixIdentity(nm1,nm2,gVerbose,epsilon);
    }
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "Test inline multiplications by a DiagMatrix" << std::endl;
    TMatrixDSym ms = THilbertMatrixDSym(msize);
    ms(1,3) = TMath::Pi();
    ms(3,1) = TMath::Pi();
    TVectorD v(msize);
    for (i = v.GetLwb(); i <= v.GetUpb(); i++)
      v(i) = 1+i;
    TMatrixDSym diag(msize);
    TMatrixDDiag d(diag); d = v;
    TMatrixDSym eth = ms;
    for (i = eth.GetRowLwb(); i <= eth.GetRowUpb(); i++)
      for (j = eth.GetColLwb(); j <= eth.GetColUpb(); j++)
        eth(i,j) *= v(j);
    TMatrixD m2 = ms * diag;
    ok &= VerifyMatrixIdentity(m2,eth,gVerbose,epsilon);
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "Test XPP = X where P is a permutation matrix" << std::endl;
    TMatrixDSym ms = THilbertMatrixDSym(msize);
    ms(2,3) = TMath::Pi();
    ms(3,2) = TMath::Pi();
    TMatrixDSym eth = ms;
    TMatrixDSym p(msize);
    for (i = p.GetRowLwb(); i <= p.GetRowUpb(); i++)
      p(p.GetRowUpb()+p.GetRowLwb()-i,i) = 1;
    TMatrixD m2 = ms * p;
    m2 *= p;
    ok &= VerifyMatrixIdentity(m2,eth,gVerbose,epsilon);
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "Test general matrix multiplication through inline mult" << std::endl;
    TMatrixDSym ms = THilbertMatrixDSym(msize);
    ms(2,3) = TMath::Pi();
    ms(3,2) = TMath::Pi();
    TMatrixDSym mt(TMatrixDSym::kTransposed,ms);
    TMatrixDSym p = THilbertMatrixDSym(msize);
    TMatrixDDiag(p) += 1;
    TMatrixD mp(ms,TMatrixD::kMult,p);
    TMatrixDSym m1 = ms;
    TMatrixD m3(ms,TMatrixD::kMult,p);
    memcpy(ms.GetMatrixArray(),m3.GetMatrixArray(),msize*msize*sizeof(Double_t));
    ok &= VerifyMatrixIdentity(ms,mp,gVerbose,epsilon);
    TMatrixD mp1(mt,TMatrixD::kTransposeMult,p);
    ok &= VerifyMatrixIdentity(ms,mp1,gVerbose,epsilon);
    ok &= ( !(m1 == ms) ) ? kTRUE : kFALSE;
    TMatrixDSym mp2(TMatrixDSym::kZero,ms);
    ok &= ( mp2 == 0 ) ? kTRUE : kFALSE;

    if (gVerbose)
      std::cout << "Test XP=X*P  vs XP=X;XP*=P" << std::endl;
    TMatrixD mp3 = m1*p;
    ok &= VerifyMatrixIdentity(ms,mp3,gVerbose,epsilon);
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "Check to see UU' = U'U = E when U is the Haar matrix" << std::endl;
    {
      const Int_t order = 5;
      const Int_t no_sub_cols = (1<<order)-5;
      TMatrixD haarb = THaarMatrixD(5,no_sub_cols);
      TMatrixD haarb_t(TMatrixD::kTransposed,haarb);
      TMatrixD hth(haarb_t,TMatrixD::kMult,haarb);
      TMatrixDSym  hth1(TMatrixDSym::kAtA,haarb);
      ok &= VerifyMatrixIdentity(hth,hth1,gVerbose,epsilon);
    }

    {
      TMatrixD haar = THaarMatrixD(5);
      TMatrixD unit(TMatrixD::kUnit,haar);
      TMatrixD haar_t(TMatrixD::kTransposed,haar);
      TMatrixDSym hths(TMatrixDSym::kAtA,haar);
      TMatrixD hht(haar,TMatrixD::kMult,haar_t);
      TMatrixD hht1 = haar; hht1 *= haar_t;
      ok &= VerifyMatrixIdentity(unit,hths,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(unit,hht,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(unit,hht1,gVerbose,epsilon);
    }
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "Check to see A.Similarity(Haar) = Haar * A * Haar^T" <<
              " and A.SimilarityT(Haar) = Haar^T * A * Haar" << std::endl;
    {
      TMatrixD    h  = THaarMatrixD(5);
      TMatrixDSym ms = THilbertMatrixDSym(1<<5);
      TMatrixD    hmht = h*TMatrixD(ms,TMatrixD::kMultTranspose,h);
      ok &= VerifyMatrixIdentity(ms.Similarity(h),hmht,gVerbose,epsilon);
      TMatrixD    htmh = TMatrixD(h,TMatrixD::kTransposeMult,ms)*h;
      ok &= VerifyMatrixIdentity(ms.SimilarityT(h),htmh,gVerbose,epsilon);
    }
    if (gVerbose)
      std::cout << "Check to see A.Similarity(B_sym) = A.Similarity(B)" << std::endl;
    {
      TMatrixDSym nsym = THilbertMatrixDSym(5);
      TMatrixD    n    = THilbertMatrixD(5,5);
      TMatrixDSym ms   = THilbertMatrixDSym(5);
      ok &= VerifyMatrixIdentity(ms.Similarity(nsym),ms.Similarity(n),gVerbose,epsilon);
    }
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(11,"Symmetric Matrix Multiplications",ok);
}

//
//------------------------------------------------------------------------
//               Verify vector-matrix multiplications
//
void mstress_vm_multiplications()
{
  Bool_t ok = kTRUE;

  Int_t iloop = gNrLoop;
  Int_t nr    = 0;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Double_t epsilon = EPSILON*msize;

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    if (verbose)
      std::cout << "\n---> Verify vector-matrix multiplications "
             "for matrices of the characteristic size " << msize << std::endl;

    {
      if (verbose)
        std::cout << "\nCheck shrinking a vector by multiplying by a non-sq unit matrix" << std::endl;
      TVectorD vb(-2,msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()-i;
      ok &= ( vb != 0 ) ? kTRUE : kFALSE;
      TMatrixD mc(1,msize-2,-2,msize);       // contracting matrix
      mc.UnitMatrix();
      TVectorD v1 = vb;
      TVectorD v2 = vb;
      v1 *= mc;
      v2.ResizeTo(1,msize-2);
      ok &= VerifyVectorIdentity(v1,v2,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Check expanding a vector by multiplying by a non-sq unit matrix" << std::endl;
      TVectorD vb(msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()+i;
      ok &= ( vb != 0 ) ? kTRUE : kFALSE;
      TMatrixD me(2,msize+5,0,msize-1);    // expanding matrix
      me.UnitMatrix();
      TVectorD v1 = vb;
      TVectorD v2 = vb;
      v1 *= me;
      v2.ResizeTo(v1);
      ok &= VerifyVectorIdentity(v1,v2,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Check general matrix-vector multiplication" << std::endl;
      TVectorD vb(msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()+i;
      TMatrixD vm(msize,1);
      TMatrixDColumn(vm,0) = vb;
      TMatrixD m = THilbertMatrixD(0,msize,0,msize-1);
      vb *= m;
      ok &= ( vb.GetLwb() == 0 ) ? kTRUE : kFALSE;
      TMatrixD mvm(m,TMatrixD::kMult,vm);
      TMatrixD mvb(msize+1,1);
      TMatrixDColumn(mvb,0) = vb;
      ok &= VerifyMatrixIdentity(mvb,mvm,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Check symmetric matrix-vector multiplication" << std::endl;
      TVectorD vb(msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()+i;
      TMatrixD vm(msize,1);
      TMatrixDColumn(vm,0) = vb;
      TMatrixDSym ms = THilbertMatrixDSym(0,msize-1);
      vb *= ms;
      ok &= ( vb.GetLwb() == 0 ) ? kTRUE : kFALSE;
      TMatrixD mvm(ms,TMatrixD::kMult,vm);
      TMatrixD mvb(msize,1);
      TMatrixDColumn(mvb,0) = vb;
      ok &= VerifyMatrixIdentity(mvb,mvm,verbose,epsilon);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop--;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        std::cout << "Vector Matrix Multiplications failed for size " << msize << std::endl;
      break;
    }
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(12,"Matrix Vector Multiplications",ok);
}

//
//------------------------------------------------------------------------
//               Verify matrix inversion
//
void mstress_inversion()
{
  Bool_t ok = kTRUE;

  Int_t iloop = gNrLoop;
  Int_t nr    = 0;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Double_t epsilon = EPSILON*msize/10;

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    if (verbose)
      std::cout << "\n---> Verify matrix inversion for square matrices of size " << msize << std::endl;
    {
      if (verbose)
        std::cout << "\nTest inversion of a diagonal matrix" << std::endl;
      TMatrixD m(-1,msize,-1,msize);
      TMatrixD mi(TMatrixD::kZero,m);
      for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
        m(i,i)=i-m.GetRowLwb()+1;
        mi(i,i) = 1/m(i,i);
      }
      TMatrixD mi1(TMatrixD::kInverted,m);
      m.Invert();
      ok &= VerifyMatrixIdentity(m,mi,verbose,epsilon);
      ok &= VerifyMatrixIdentity(mi1,mi,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Test inversion of an orthonormal (Haar) matrix" << std::endl;
      const Int_t order = Int_t(TMath::Log(msize)/TMath::Log(2));
      TMatrixD mr = THaarMatrixD(order);
      TMatrixD morig = mr;
      TMatrixD mt(TMatrixD::kTransposed,mr);
      Double_t det = -1;         // init to a wrong val to see if it's changed
      mr.Invert(&det);
      ok &= VerifyMatrixIdentity(mr,mt,verbose,epsilon);
      ok &= ( TMath::Abs(det-1) <= msize*epsilon ) ? kTRUE : kFALSE;
      if (verbose) {
        std::cout << "det = " << det << " deviation= " << TMath::Abs(det-1);
        std::cout << ( (TMath::Abs(det-1) < msize*epsilon) ? " OK" : " too large") <<std::endl;
      }
      TMatrixD mti(TMatrixD::kInverted,mt);
      ok &= VerifyMatrixIdentity(mti,morig,verbose,msize*epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Test inversion of a good matrix with diagonal dominance" << std::endl;
      TMatrixD mh = THilbertMatrixD(msize,msize);
      TMatrixDDiag(mh) += 1;
      TMatrixD morig = mh;
      Double_t det_inv = 0;
      const Double_t det_comp = mh.Determinant();
      mh.Invert(&det_inv);
      if (verbose) {
        std::cout << "\tcomputed determinant             " << det_comp << std::endl;
        std::cout << "\tdeterminant returned by Invert() " << det_inv << std::endl;
      }

      if (verbose)
        std::cout << "\tcheck to see M^(-1) * M is E" << std::endl;
      TMatrixD mim(mh,TMatrixD::kMult,morig);
      TMatrixD unit(TMatrixD::kUnit,mh);
      ok &= VerifyMatrixIdentity(mim,unit,verbose,epsilon);

      if (verbose)
        std::cout << "Test inversion through the matrix decompositions" << std::endl;

      TMatrixDSym ms = THilbertMatrixDSym(msize);
      TMatrixDDiag(ms) += 1;
      if (verbose)
        std::cout << "Test inversion through SVD" << std::endl;
      TMatrixD inv_svd (msize,msize); TDecompSVD  svd (ms); svd.Invert(inv_svd);
      ok &= VerifyMatrixIdentity(inv_svd,mh,verbose,epsilon);
      if (verbose)
        std::cout << "Test inversion through LU" << std::endl;
      TMatrixD inv_lu  (msize,msize); TDecompLU   lu  (ms); lu.Invert(inv_lu);
      ok &= VerifyMatrixIdentity(inv_lu,mh,verbose,epsilon);
      if (verbose)
        std::cout << "Test inversion through Cholesky" << std::endl;
      TMatrixDSym inv_chol(msize); TDecompChol chol(ms); chol.Invert(inv_chol);
      ok &= VerifyMatrixIdentity(inv_chol,mh,verbose,epsilon);
      if (verbose)
        std::cout << "Test inversion through QRH" << std::endl;
      TMatrixD inv_qrh (msize,msize); TDecompQRH  qrh (ms); qrh.Invert(inv_qrh);
      ok &= VerifyMatrixIdentity(inv_qrh,mh,verbose,epsilon);
      if (verbose)
        std::cout << "Test inversion through Bunch-Kaufman" << std::endl;
      TMatrixDSym inv_bk(msize); TDecompBK bk(ms); chol.Invert(inv_bk);
      ok &= VerifyMatrixIdentity(inv_bk,mh,verbose,epsilon);

      if (verbose)
        std::cout << "\tcheck to see M * M^(-1) is E" << std::endl;
      TMatrixD mmi = morig; mmi *= mh;
      ok &= VerifyMatrixIdentity(mmi,unit,verbose,epsilon);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop--;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        std::cout << "Matrix Inversion failed for size " << msize << std::endl;
      break;
    }
  }

  {
    if (gVerbose) {
      std::cout << "Check to see that Invert() and InvertFast() give identical results" <<std::endl;
      std::cout << " for size < (7x7)" <<std::endl;
    }
    Int_t size;
    for (size = 2; size < 7; size++) {
      TMatrixD m1 = THilbertMatrixD(size,size);
      m1(0,1) = TMath::Pi();
      TMatrixDDiag(m1) += 1;
      TMatrixD m2 = m1;
      Double_t det1 = 0.0;
      Double_t det2 = 0.0;
      m1.Invert(&det1);
      m2.InvertFast(&det2);
      ok &= VerifyMatrixIdentity(m1,m2,gVerbose,EPSILON);
      ok &= (TMath::Abs(det1-det2) < EPSILON);
      if (gVerbose) {
        std::cout << "det(Invert)= " << det1 << "  det(InvertFast)= " << det2 <<std::endl;
        std::cout << " deviation= " << TMath::Abs(det1-det2);
        std::cout << ( (TMath::Abs(det1-det2) <  EPSILON) ? " OK" : " too large") <<std::endl;
      }
    }
    for (size = 2; size < 7; size++) {
      TMatrixDSym ms1 = THilbertMatrixDSym(size);
      TMatrixDDiag(ms1) += 1;
      TMatrixDSym ms2 = ms1;
      Double_t det1 = 0.0;
      Double_t det2 = 0.0;
      ms1.Invert(&det1);
      ms2.InvertFast(&det2);
      ok &= VerifyMatrixIdentity(ms1,ms2,gVerbose,EPSILON);
      ok &= (TMath::Abs(det1-det2) < EPSILON);
      if (gVerbose) {
        std::cout << "det(Invert)= " << det1 << "  det(InvertFast)= " << det2 <<std::endl;
        std::cout << " deviation= " << TMath::Abs(det1-det2);
        std::cout << ( (TMath::Abs(det1-det2) <  EPSILON) ? " OK" : " too large") <<std::endl;
      }
    }
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(13,"Matrix Inversion",ok);
}

//
//------------------------------------------------------------------------
//           Test matrix I/O
//
void mstress_matrix_io()
{
  if (gVerbose)
    std::cout << "\n---> Test matrix I/O" << std::endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = TMath::Pi();

  TFile *f = new TFile("stress-vmatrix.root", "RECREATE");

  Char_t name[80];
  Int_t iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TMatrixD m(msize,msize);
    m = pattern;

    Double_t *pattern_array = new Double_t[msize*msize];
    for (Int_t i = 0; i < msize*msize; i++)
      pattern_array[i] = pattern;
    TMatrixD ma;
    ma.Use(msize,msize,pattern_array);

    TMatrixDSym ms(msize);
    ms = pattern;

    if (verbose)
      std::cout << "\nWrite matrix m to database" << std::endl;
    snprintf(name,80,"m_%d",msize);
    m.Write(name);

    if (verbose)
      std::cout << "\nWrite matrix ma which adopts to database" << std::endl;
    snprintf(name,80,"ma_%d",msize);
    ma.Write(name);

    if (verbose)
      std::cout << "\nWrite symmetric matrix ms to database" << std::endl;
    snprintf(name,80,"ms_%d",msize);
    ms.Write(name);

    delete [] pattern_array;

    iloop--;
  }

  if (gVerbose)
    std::cout << "\nClose database" << std::endl;
  delete f;

  if (gVerbose)
    std::cout << "\nOpen database in read-only mode and read matrix" << std::endl;
  TFile *f1 = new TFile("stress-vmatrix.root");

  iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TMatrixD m(msize,msize);
    m = pattern;
    snprintf(name,80,"m_%d",msize);
    TMatrixD *mr  = (TMatrixD*) f1->Get(name);
    snprintf(name,80,"ma_%d",msize);
    TMatrixD *mar = (TMatrixD*) f1->Get(name);
    snprintf(name,80,"ms_%d",msize);
    TMatrixDSym *msr = (TMatrixDSym*) f1->Get(name);

    if (verbose)
      std::cout << "\nRead matrix should be same as original" << std::endl;
    ok &= ((*mr)  == m) ? kTRUE : kFALSE;
    ok &= ((*mar) == m) ? kTRUE : kFALSE;
    ok &= ((*msr) == m) ? kTRUE : kFALSE;

    iloop--;
  }

  delete f1;

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(14,"Matrix Persistence",ok);
}

//------------------------------------------------------------------------
//          Test allocation functions and compatibility check
//
void spstress_allocation(Int_t msize)
{
  if (gVerbose)
    std::cout << "\n\n---> Test allocation and compatibility check" << std::endl;

  Int_t i,j;
  Bool_t ok = kTRUE;

  TMatrixDSparse m1(0,3,0,msize-1);
  {
    Int_t nr = 4*msize;
    Int_t    *irow = new Int_t[nr];
    Int_t    *icol = new Int_t[nr];
    Double_t *val  = new Double_t[nr];

    Int_t n = 0;
    for (i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++) {
      for (j = m1.GetColLwb(); j <= m1.GetColUpb(); j++) {
        irow[n] = i;
        icol[n] = j;
        val[n] = TMath::Pi()*i+TMath::E()*j;
        n++;
      }
    }
    m1.SetMatrixArray(nr,irow,icol,val);
    delete [] irow;
    delete [] icol;
    delete [] val;
  }

  TMatrixDSparse m2(0,3,0,msize-1);
  TMatrixDSparse m3(1,4,0,msize-1);
  TMatrixDSparse m4(m1);

  if (gVerbose) {
    std::cout << "\nStatus information reported for matrix m3:" << std::endl;
    std::cout << "  Row lower bound ... " << m3.GetRowLwb() << std::endl;
    std::cout << "  Row upper bound ... " << m3.GetRowUpb() << std::endl;
    std::cout << "  Col lower bound ... " << m3.GetColLwb() << std::endl;
    std::cout << "  Col upper bound ... " << m3.GetColUpb() << std::endl;
    std::cout << "  No. rows ..........." << m3.GetNrows()  << std::endl;
    std::cout << "  No. cols ..........." << m3.GetNcols()  << std::endl;
    std::cout << "  No. of elements ...." << m3.GetNoElements() << std::endl;
  }

  if (gVerbose)
    std::cout << "Check matrices 1 & 4 for compatibility" << std::endl;
  ok &= AreCompatible(m1,m4,gVerbose);

  if (gVerbose)
    std::cout << "m2 has to be compatible with m3 after resizing to m3" << std::endl;
  m2.ResizeTo(m3);
  ok &= AreCompatible(m2,m3,gVerbose);

  TMatrixD m5_d(m1.GetNrows()+1,m1.GetNcols()+5);
  for (i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++)
    for (j = m1.GetColLwb(); j <= m1.GetColUpb(); j++)
      m5_d(i,j) = TMath::Pi()*i+TMath::E()*j;
  TMatrixDSparse m5(m5_d);

  if (gVerbose)
    std::cout << "m1 has to be compatible with m5 after resizing to m5" << std::endl;
  m1.ResizeTo(m5.GetNrows(),m5.GetNcols());
  ok &= AreCompatible(m1,m5,gVerbose);

  if (gVerbose)
    std::cout << "m1 has to be equal to m4 after stretching and shrinking" << std::endl;
  m1.ResizeTo(m4.GetNrows(),m4.GetNcols());
  ok &= VerifyMatrixIdentity(m1,m4,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "m5 has to be equal to m1 after shrinking" << std::endl;
  m5.ResizeTo(m1.GetNrows(),m1.GetNcols());
  ok &= VerifyMatrixIdentity(m1,m5,gVerbose,msize*EPSILON);

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(1,"Allocation, Resizing",ok);
}

//
//------------------------------------------------------------------------
//          Test Filling of matrix
//
void spstress_matrix_fill(Int_t rsize,Int_t csize)
{
  Bool_t ok = kTRUE;

  if (rsize < 4) {
     Error("spstress_matrix_fill", "rsize should be >= 4");
     ok = kFALSE;
     StatusPrint(2, "Filling, Inserting, Using", ok);
     return;
  }

  if (csize < 4) {
    Error("spstress_matrix_fill","csize should be >= 4");
    ok = kFALSE;
    StatusPrint(2,"Filling, Inserting, Using",ok);
    return;
  }

  TMatrixD m_d(-1,rsize-2,1,csize);
  for (Int_t i = m_d.GetRowLwb(); i <= m_d.GetRowUpb(); i++)
    for (Int_t j = m_d.GetColLwb(); j <= m_d.GetColUpb(); j++)
      // Keep column 2 on zero
      if (j != 2)
        m_d(i,j) = TMath::Pi()*i+TMath::E()*j;
  TMatrixDSparse m(m_d);

  {
    if (gVerbose)
      std::cout << "Check filling through operator(i,j) without setting sparse index" << std::endl;
    TMatrixDSparse m1(-1,rsize-2,1,csize);

    for (Int_t i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++)
      for (Int_t j = m1.GetColLwb(); j <= m1.GetColUpb(); j++)
        if (j != 2)
          m1(i,j) = TMath::Pi()*i+TMath::E()*j;
    ok &= VerifyMatrixIdentity(m1,m,gVerbose,EPSILON);
  }

  {
    if (gVerbose)
      std::cout << "Check filling through operator(i,j)" << std::endl;
    TMatrixDSparse m2(-1,rsize-2,1,csize);
    m2.SetSparseIndex(m);

    for (Int_t i = m2.GetRowLwb(); i <= m2.GetRowUpb(); i++)
      for (Int_t j = m2.GetColLwb(); j <= m2.GetColUpb(); j++)
        if (j != 2)
          m2(i,j) = TMath::Pi()*i+TMath::E()*j;
    ok &= VerifyMatrixIdentity(m2,m,gVerbose,EPSILON);
  }

  {
    if (gVerbose)
      std::cout << "Check insertion/extraction of sub-matrices" << std::endl;
    {
      TMatrixDSparse m_sub1 = m;
      m_sub1.ResizeTo(0,rsize-2,2,csize);
      TMatrixDSparse m_sub2 = m.GetSub(0,rsize-2,2,csize,"");
      ok &= VerifyMatrixIdentity(m_sub1,m_sub2,gVerbose,EPSILON);
    }

    {
      TMatrixDSparse m3(-1,rsize-2,1,csize);
      TMatrixDSparse m_part1 = m.GetSub(-1,rsize-2,1,csize,"");
      m3.SetSub(-1,1,m_part1);
      ok &= VerifyMatrixIdentity(m,m3,gVerbose,EPSILON);
    }

    {
      TMatrixDSparse m4(-1,rsize-2,1,csize);
      TMatrixDSparse m_part1 = m.GetSub(0,rsize-2,2,csize,"");
      TMatrixDSparse m_part2 = m.GetSub(0,rsize-2,1,1,"");
      TMatrixDSparse m_part3 = m.GetSub(-1,-1,2,csize,"");
      TMatrixDSparse m_part4 = m.GetSub(-1,-1,1,1,"");
      m4.SetSub(0,2,m_part1);
      m4.SetSub(0,1,m_part2);
      m4.SetSub(-1,2,m_part3);
      m4.SetSub(-1,1,m_part4);
      ok &= VerifyMatrixIdentity(m,m4,gVerbose,EPSILON);
    }

    {
      // change the insertion order
      TMatrixDSparse m5(-1,rsize-2,1,csize);
      TMatrixDSparse m_part1 = m.GetSub(0,rsize-2,2,csize,"");
      TMatrixDSparse m_part2 = m.GetSub(0,rsize-2,1,1,"");
      TMatrixDSparse m_part3 = m.GetSub(-1,-1,2,csize,"");
      TMatrixDSparse m_part4 = m.GetSub(-1,-1,1,1,"");
      m5.SetSub(-1,1,m_part4);
      m5.SetSub(-1,2,m_part3);
      m5.SetSub(0,1,m_part2);
      m5.SetSub(0,2,m_part1);
      ok &= VerifyMatrixIdentity(m,m5,gVerbose,EPSILON);
    }

    {
      TMatrixDSparse m6(-1,rsize-2,1,csize);
      TMatrixDSparse m_part1 = m.GetSub(0,rsize-2,2,csize,"S");
      TMatrixDSparse m_part2 = m.GetSub(0,rsize-2,1,1,"S");
      TMatrixDSparse m_part3 = m.GetSub(-1,-1,2,csize,"S");
      TMatrixDSparse m_part4 = m.GetSub(-1,-1,1,1,"S");
      m6.SetSub(0,2,m_part1);
      m6.SetSub(0,1,m_part2);
      m6.SetSub(-1,2,m_part3);
      m6.SetSub(-1,1,m_part4);
      ok &= VerifyMatrixIdentity(m,m6,gVerbose,EPSILON);
    }
  }

  {
    if (gVerbose)
      std::cout << "Check insertion/extraction of rows" << std::endl;

    TMatrixDSparse m1 = m;
    TVectorD v1(1,csize);
    TVectorD v2(1,csize);
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
      v1 = TMatrixDSparseRow(m,i);
      m1.InsertRow(i,1,v1.GetMatrixArray());
      ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
      m1.InsertRow(i,3,v1.GetMatrixArray()+2,v1.GetNrows()-2);
      ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);

      m1.ExtractRow(i,1,v2.GetMatrixArray());
      ok &= VerifyVectorIdentity(v1,v2,gVerbose,EPSILON);
      m1.ExtractRow(i,3,v2.GetMatrixArray()+2,v1.GetNrows()-2);
      ok &= VerifyVectorIdentity(v1,v2,gVerbose,EPSILON);
    }
  }

  {
    if (gVerbose)
      std::cout << "Check array Use" << std::endl;
    {
      TMatrixDSparse *m1a = new TMatrixDSparse(m);
      TMatrixDSparse *m2a = new TMatrixDSparse();
      m2a->Use(*m1a);
      m2a->Sqr();
      TMatrixDSparse m7 = m; m7.Sqr();
      ok &= VerifyMatrixIdentity(m7,*m1a,gVerbose,EPSILON);
      delete m1a;
      delete m2a;
    }
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(2,"Filling, Inserting, Using",ok);
}

//
//------------------------------------------------------------------------
//                Test uniform element operations
//
void spstress_element_op(Int_t rsize,Int_t csize)
{
  Bool_t ok = kTRUE;
  const Double_t pattern = 8.625;

  TMatrixDSparse m(-1,rsize-2,1,csize);

  if (gVerbose)
    std::cout << "Creating zero m1 ..." << std::endl;
  TMatrixDSparse m1(TMatrixDSparse::kZero, m);
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "Comparing m1 with 0 ..." << std::endl;
  R__ASSERT(m1 == 0);
  R__ASSERT(!(m1 != 0));

  if (gVerbose)
    std::cout << "Writing a pattern " << pattern << " by assigning through SetMatrixArray..." << std::endl;
  {
    const Int_t nr = rsize*csize;
    Int_t    *irow = new Int_t[nr];
    Int_t    *icol = new Int_t[nr];
    Double_t *val  = new Double_t[nr];

    Int_t n = 0;
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++) {
        irow[n] = i;
        icol[n] = j;
        val[n] = pattern;
        n++;
      }
    }
    m.SetMatrixArray(nr,irow,icol,val);
    delete [] irow;
    delete [] icol;
    delete [] val;
    ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  }

  if (gVerbose)
    std::cout << "Writing the pattern by assigning to m1 as a whole ..."  << std::endl;
  m1.SetSparseIndex(m);
  m1 = pattern;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "Comparing m and m1 ..." << std::endl;
  R__ASSERT(m == m1);
  if (gVerbose)
    std::cout << "Comparing (m=0) and m1 ..." << std::endl;
  R__ASSERT(!((m=0) == m1));

  if (gVerbose)
    std::cout << "Clearing m1 ..." << std::endl;
  m1.Zero();
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nSet m = pattern" << std::endl;
  m = pattern;
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "   add the doubled pattern with the negative sign" << std::endl;
  m += -2*pattern;
  ok &= VerifyMatrixValue(m,-pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "   subtract the trippled pattern with the negative sign" << std::endl;
  m -= -3*pattern;
  ok &= VerifyMatrixValue(m,2*pattern,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nVerify comparison operations when all elems are the same" << std::endl;
  m = pattern;
  R__ASSERT( m == pattern && !(m != pattern) );
  R__ASSERT( m > 0 && m >= pattern && m <= pattern );
  R__ASSERT( m > -pattern && m >= -pattern );
  R__ASSERT( m <= pattern && !(m < pattern) );
  m -= 2*pattern;
  R__ASSERT( m  < -pattern/2 && m <= -pattern/2 );
  R__ASSERT( m  >= -pattern && !(m > -pattern) );

  if (gVerbose)
    std::cout << "\nVerify comparison operations when not all elems are the same" << std::endl;
  {
    Int_t nr = rsize*csize;
    Int_t    *irow = new Int_t[nr];
    Int_t    *icol = new Int_t[nr];
    Double_t *val  = new Double_t[nr];

    Int_t n = 0;
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++) {
        irow[n] = i;
        icol[n] = j;
        val[n] = pattern;
        n++;
      }
    }
    val[n-1] = pattern-1;
    m.SetMatrixArray(nr,irow,icol,val);
    delete [] irow;
    delete [] icol;
    delete [] val;
  }

  R__ASSERT( !(m == pattern) && !(m != pattern) );
  R__ASSERT( m != 0 );                   // none of elements are 0
  R__ASSERT( !(m >= pattern) && m <= pattern && !(m<pattern) );
  R__ASSERT( !(m <= pattern-1) && m >= pattern-1 && !(m>pattern-1) );
  if (gVerbose)
    std::cout << "\nAssign 2*pattern to m by repeating additions" << std::endl;
  m = 0; m += pattern; m += pattern;
  if (gVerbose)
    std::cout << "Assign 2*pattern to m1 by multiplying by two " << std::endl;
  m1.SetSparseIndex(m);
  m1 = pattern; m1 *= 2;
  ok &= VerifyMatrixValue(m1,2*pattern,gVerbose,EPSILON);
  R__ASSERT( m == m1 );
  if (gVerbose)
    std::cout << "Multiply m1 by one half returning it to the 1*pattern" << std::endl;
  m1 *= 1/2.;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nAssign -pattern to m and m1" << std::endl;
  m = 0; m -= pattern; m1 = -pattern;
  ok &= VerifyMatrixValue(m,-pattern,gVerbose,EPSILON);
  R__ASSERT( m == m1 );
  if (gVerbose)
    std::cout << "m = sqrt(sqr(m)); m1 = abs(m1); Now m and m1 have to be the same" << std::endl;
  m.Sqr();
  ok &= VerifyMatrixValue(m,pattern*pattern,gVerbose,EPSILON);
  m.Sqrt();
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  m1.Abs();
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  ok &= VerifyMatrixIdentity(m1,m,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(3,"Uniform matrix operations",ok);
}

//
//------------------------------------------------------------------------
//        Test binary matrix element-by-element operations
//
void spstress_binary_ebe_op(Int_t rsize, Int_t csize)
{
  if (gVerbose)
    std::cout << "\n---> Test Binary Matrix element-by-element operations" << std::endl;

  Bool_t ok = kTRUE;

  const Double_t pattern = 4.25;

  TMatrixD m_d(2,rsize+1,0,csize-1); m_d = 1;
  TMatrixDSparse m (TMatrixDSparse::kZero,m_d); m.SetSparseIndex (m_d);
  TMatrixDSparse m1(TMatrixDSparse::kZero,m);   m1.SetSparseIndex(m_d);

  TMatrixDSparse mp(TMatrixDSparse::kZero,m1);  mp.SetSparseIndex(m_d);
  {
    for (Int_t i = mp.GetRowLwb(); i <= mp.GetRowUpb(); i++)
      for (Int_t j = mp.GetColLwb(); j <= mp.GetColUpb(); j++)
        mp(i,j) = TMath::Pi()*i+TMath::E()*(j+1);
  }

  if (gVerbose)
    std::cout << "\nVerify assignment of a matrix to the matrix" << std::endl;
  m  = pattern;
  m1 = m;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  R__ASSERT( m1 == m );

  if (gVerbose)
    std::cout << "\nAdding the matrix to itself, uniform pattern " << pattern << std::endl;
  m.Zero(); m.SetSparseIndex(m_d); m = pattern;

  m1 = m; m1 += m1;
  ok &= VerifyMatrixValue(m1,2*pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "  subtracting two matrices ..." << std::endl;
  m1 -= m;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "  subtracting the matrix from itself" << std::endl;
  // Hiding the self-subtraction from the compiler, causing warning by -Wself-assign-overloaded in clang 7.0
  m1 -= static_cast<TMatrixDSparse&>(m1);
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);
  m1.SetSparseIndex(m_d);

  if (gVerbose) {
    std::cout << "\nArithmetic operations on matrices with not the same elements" << std::endl;
    std::cout << "   adding mp to the zero matrix..." << std::endl;
  }
  m.Zero(); m.SetSparseIndex(m_d); m += mp;
  ok &= VerifyMatrixIdentity(m,mp,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "   making m = 3*mp and m1 = 3*mp, via add() and succesive mult" << std::endl;
  m1 = m;
  Add(m,2.,mp);
  m1 += m1; m1 += mp;
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "   clear both m and m1, by subtracting from itself and via add()" << std::endl;
  // Hiding the self-subtraction from the compiler, causing warning by -Wself-assign-overloaded in clang 7.0
  m1 -= static_cast<TMatrixDSparse&>(m1);
  Add(m,-3.,mp);
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);

  if (gVerbose) {
    std::cout << "\nTesting element-by-element multiplications and divisions" << std::endl;
    std::cout << "   squaring each element with sqr() and via multiplication" << std::endl;
  }
  m.SetSparseIndex(m_d);  m = mp;
  m1.SetSparseIndex(m_d); m1 = mp;
  m.Sqr();
  ElementMult(m1,m1);
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "   compare (m = pattern^2)/pattern with pattern" << std::endl;
  m = pattern; m1 = pattern;
  m.Sqr();
  ElementDiv(m,m1);
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  if (gVerbose)
    Compare(m1,m);

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(4,"Binary Matrix element-by-element operations",ok);
}

//
//------------------------------------------------------------------------
//              Verify matrix transposition
//
void spstress_transposition(Int_t msize)
{
  if (gVerbose) {
    std::cout << "\n---> Verify matrix transpose "
            "for matrices of a characteristic size " << msize << std::endl;
  }

  Bool_t ok = kTRUE;
  {
    if (gVerbose)
      std::cout << "\nCheck to see that a square UnitMatrix stays the same";
    TMatrixDSparse m(msize,msize);
    m.UnitMatrix();
    TMatrixDSparse mt(TMatrixDSparse::kTransposed,m);
    ok &= ( m == mt ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      std::cout << "\nTest a non-square UnitMatrix";
    TMatrixDSparse m(msize,msize+1);
    m.UnitMatrix();
    TMatrixDSparse mt(TMatrixDSparse::kTransposed,m);
    R__ASSERT(m.GetNrows() == mt.GetNcols() && m.GetNcols() == mt.GetNrows() );
    const Int_t rowlwb = m.GetRowLwb();
    const Int_t collwb = m.GetColLwb();
    const Int_t upb = TMath::Min(m.GetRowUpb(),m.GetColUpb());
    TMatrixDSparse m_part  = m.GetSub(rowlwb,upb,collwb,upb);
    TMatrixDSparse mt_part = mt.GetSub(rowlwb,upb,collwb,upb);
    ok &= ( m_part == mt_part ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      std::cout << "\nCheck to see that a symmetric (Hilbert)Matrix stays the same";
    TMatrixDSparse m = TMatrixD(THilbertMatrixD(msize,msize));
    TMatrixDSparse mt(TMatrixDSparse::kTransposed,m);
    ok &= ( m == mt ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      std::cout << "\nCheck transposing a non-symmetric matrix";
    TMatrixDSparse m = TMatrixD(THilbertMatrixD(msize+1,msize));
    m(1,2) = TMath::Pi();
    TMatrixDSparse mt(TMatrixDSparse::kTransposed,m);
    R__ASSERT(m.GetNrows() == mt.GetNcols() && m.GetNcols() == mt.GetNrows());
    R__ASSERT(mt(2,1)  == (Double_t)TMath::Pi() && mt(1,2)  != (Double_t)TMath::Pi());
    R__ASSERT(mt[2][1] == (Double_t)TMath::Pi() && mt[1][2] != (Double_t)TMath::Pi());

    if (gVerbose)
      std::cout << "\nCheck double transposing a non-symmetric matrix" << std::endl;
    TMatrixDSparse mtt(TMatrixDSparse::kTransposed,mt);
    ok &= ( m == mtt ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(5,"Matrix transposition",ok);
}

//
//------------------------------------------------------------------------
//             Verify the norm calculation
//
void spstress_norms(Int_t rsize_req,Int_t csize_req)
{
  if (gVerbose)
    std::cout << "\n---> Verify norm calculations" << std::endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 10.25;

  Int_t rsize = rsize_req;
  if (rsize%2 != 0)  rsize--;
  Int_t csize = csize_req;
  if (csize%2 != 0)  csize--;
  if (rsize%2 == 1 || csize%2 == 1) {
    std::cout << "rsize: " << rsize <<std::endl;
    std::cout << "csize: " << csize <<std::endl;
    Fatal("spstress_norms","Sorry, size of the matrix to test must be even for this test\n");
  }

  TMatrixD m_d(rsize,csize); m_d = 1;
  TMatrixDSparse m(rsize,csize); m.SetSparseIndex(m_d);

  if (gVerbose)
    std::cout << "\nAssign " << pattern << " to all the elements and check norms" << std::endl;
  m = pattern;
  if (gVerbose)
    std::cout << "  1. (col) norm should be pattern*nrows" << std::endl;
  ok &= ( m.Norm1() == pattern*m.GetNrows() ) ? kTRUE : kFALSE;
  ok &= ( m.Norm1() == m.ColNorm() ) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "  Inf (row) norm should be pattern*ncols" << std::endl;
  ok &= ( m.NormInf() == pattern*m.GetNcols() ) ? kTRUE : kFALSE;
  ok &= ( m.NormInf() == m.RowNorm() ) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "  Square of the Eucl norm has got to be pattern^2 * no_elems" << std::endl;
  ok &= ( m.E2Norm() == (pattern*pattern)*m.GetNoElements() ) ? kTRUE : kFALSE;
  TMatrixDSparse m1(m); m1 = 1;
  ok &= ( m.E2Norm() == E2Norm(m+1.,m1) ) ? kTRUE : kFALSE;

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(6,"Matrix Norms",ok);
}

//
//------------------------------------------------------------------------
//               Verify matrix multiplications
//
void spstress_mm_multiplications()
{
  Bool_t ok = kTRUE;
  Int_t i;

  Int_t iloop = 0;
  Int_t nr    = 0;
  while (iloop <= gNrLoop) {
    const Int_t msize = gSizeA[iloop];
    const Double_t epsilon = EPSILON*msize/100;

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    if (msize <= 5) {
       iloop++;
       continue;  // some references to m(3,..)
    }

    if (verbose)
      std::cout << "\n---> Verify matrix multiplications "
              "for matrices of the characteristic size " << msize << std::endl;

    {
      if (verbose)
        std::cout << "\nTest inline multiplications of the UnitMatrix" << std::endl;
      if (ok)
      {
        TMatrixDSparse m = TMatrixD(THilbertMatrixD(-1,msize,-1,msize));
        TMatrixDSparse ur(TMatrixDSparse::kUnit,m);
        TMatrixDSparse ul(TMatrixDSparse::kUnit,m);
        m(3,1) = TMath::Pi();
        ul *= m;
        m  *= ur;
        ok &= VerifyMatrixIdentity(ul,m,verbose,epsilon);
      }

      if (ok)
      {
        TMatrixD m_d = THilbertMatrixD(-1,msize,-1,msize);
        TMatrixDSparse ur(TMatrixDSparse::kUnit,m_d);
        TMatrixDSparse ul(TMatrixDSparse::kUnit,m_d);
        m_d(3,1) = TMath::Pi();
        ul *= m_d;
        m_d  *= TMatrixD(ur);
        ok &= VerifyMatrixIdentity(ul,m_d,verbose,epsilon);
      }
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Test XPP = X where P is a permutation matrix" << std::endl;
      TMatrixDSparse m = TMatrixD(THilbertMatrixD(msize-1,msize));
      m(2,3) = TMath::Pi();
      TMatrixDSparse eth = m;
      TMatrixDSparse p(msize,msize);
      {
        Int_t    *irow = new Int_t[msize];
        Int_t    *icol = new Int_t[msize];
        Double_t *val  = new Double_t[msize];

        Int_t n = 0;
        for (i = p.GetRowLwb(); i <= p.GetRowUpb(); i++) {
          irow[n] = p.GetRowUpb()+p.GetRowLwb()-i;
          icol[n] = i;
          val[n] = 1;
          n++;
        }
        p.SetMatrixArray(msize,irow,icol,val);
        delete [] irow;
        delete [] icol;
        delete [] val;
      }
      m *= p;
      m *= p;
      ok &= VerifyMatrixIdentity(m,eth,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Test general matrix multiplication through inline mult" << std::endl;
      TMatrixDSparse m = TMatrixD(THilbertMatrixD(msize-2,msize));
      m(3,3) = TMath::Pi();
      TMatrixD p_d = THilbertMatrixD(msize,msize);
      TMatrixDDiag(p_d) += 1;
      TMatrixDSparse mp(m,TMatrixDSparse::kMult,p_d);
      TMatrixDSparse m1 = m;
      m *= p_d;
      ok &= VerifyMatrixIdentity(m,mp,verbose,epsilon);
      TMatrixDSparse pt_d(TMatrixDSparse::kTransposed,p_d);
      TMatrixDSparse mp1(m1,TMatrixDSparse::kMultTranspose,pt_d);
      VerifyMatrixIdentity(m,mp1,verbose,epsilon);

      ok &= ( !(m1 == m) );
      TMatrixDSparse mp2(TMatrixDSparse::kZero,m1);
      ok &= ( mp2 == 0 );
      mp2.SetSparseIndex(m1);
      mp2.Mult(m1,p_d);
      ok &= VerifyMatrixIdentity(m,mp2,verbose,epsilon);

      if (verbose)
        std::cout << "Test XP=X*P  vs XP=X;XP*=P" << std::endl;
      TMatrixDSparse mp3 = m1*p_d;
      ok &= VerifyMatrixIdentity(m,mp3,verbose,epsilon);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop++;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        std::cout << "General Matrix Multiplications failed for size " << msize << std::endl;
      break;
    }
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "Check to see UU' = U'U = E when U is the Haar matrix" << std::endl;
    const Int_t order = 5;
    const Int_t no_sub_cols = (1<<order)-order;
    TMatrixDSparse haar_sub = TMatrixD(THaarMatrixD(order,no_sub_cols));
    TMatrixDSparse haar_sub_t(TMatrixDSparse::kTransposed,haar_sub);
    TMatrixDSparse hsths(haar_sub_t,TMatrixDSparse::kMult,haar_sub);

    TMatrixDSparse square(no_sub_cols,no_sub_cols);
    TMatrixDSparse units(TMatrixDSparse::kUnit,square);
    ok &= ( hsths.GetNrows() == no_sub_cols && hsths.GetNcols() == no_sub_cols );
    ok &= VerifyMatrixIdentity(hsths,units,gVerbose,EPSILON);

    TMatrixDSparse haar = TMatrixD(THaarMatrixD(order));
    TMatrixDSparse haar_t(TMatrixDSparse::kTransposed,haar);
    TMatrixDSparse unit(TMatrixDSparse::kUnit,haar);

    TMatrixDSparse hth(haar_t,TMatrixDSparse::kMult,haar);
    ok &= VerifyMatrixIdentity(hth,unit,gVerbose,EPSILON);

    TMatrixDSparse hht(haar,TMatrixDSparse::kMultTranspose,haar);
    ok &= VerifyMatrixIdentity(hht,unit,gVerbose,EPSILON);

    TMatrixDSparse hht1 = haar; hht1 *= haar_t;
    ok &= VerifyMatrixIdentity(hht1,unit,gVerbose,EPSILON);

    TMatrixDSparse hht2(TMatrixDSparse::kZero,haar);
    hht2.SetSparseIndex(hht1);
    hht2.Mult(haar,haar_t);
    ok &= VerifyMatrixIdentity(hht2,unit,gVerbose,EPSILON);
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(7,"General Matrix Multiplications",ok);
}

//
//------------------------------------------------------------------------
//               Verify vector-matrix multiplications
//
void spstress_vm_multiplications()
{
  Bool_t ok = kTRUE;

  Int_t iloop = gNrLoop;
  Int_t nr    = 0;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Double_t epsilon = EPSILON*msize;

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    if (verbose)
      std::cout << "\n---> Verify vector-matrix multiplications "
             "for matrices of the characteristic size " << msize << std::endl;

    {
      if (verbose)
        std::cout << "\nCheck shrinking a vector by multiplying by a non-sq unit matrix" << std::endl;
      TVectorD vb(-2,msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()-i;
      ok &= ( vb != 0 ) ? kTRUE : kFALSE;
      TMatrixDSparse mc(1,msize-2,-2,msize);       // contracting matrix
      mc.UnitMatrix();
      TVectorD v1 = vb;
      TVectorD v2 = vb;
      v1 *= mc;
      v2.ResizeTo(1,msize-2);
      ok &= VerifyVectorIdentity(v1,v2,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Check expanding a vector by multiplying by a non-sq unit matrix" << std::endl;
      TVectorD vb(msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()+i;
      ok &= ( vb != 0 ) ? kTRUE : kFALSE;
      TMatrixDSparse me(2,msize+5,0,msize-1);    // expanding matrix
      me.UnitMatrix();
      TVectorD v1 = vb;
      TVectorD v2 = vb;
      v1 *= me;
      v2.ResizeTo(v1);
      ok &= VerifyVectorIdentity(v1,v2,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        std::cout << "Check general matrix-vector multiplication" << std::endl;
      TVectorD vb(msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()+i;
      TMatrixD vm(msize,1);
      TMatrixDColumn(vm,0) = vb;
      TMatrixD hilbert_with_zeros = THilbertMatrixD(0,msize,0,msize-1);
      TMatrixDRow   (hilbert_with_zeros,3).Assign(0.0);
      TMatrixDColumn(hilbert_with_zeros,3) = 0.0;
      const TMatrixDSparse m = hilbert_with_zeros;
      vb *= m;
      ok &= ( vb.GetLwb() == 0 ) ? kTRUE : kFALSE;
      TMatrixDSparse mvm(m,TMatrixDSparse::kMult,vm);
      TMatrixD mvb(msize+1,1);
      TMatrixDColumn(mvb,0) = vb;
      ok &= VerifyMatrixIdentity(mvb,mvm,verbose,epsilon);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop--;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        std::cout << "Vector Matrix Multiplications failed for size " << msize << std::endl;
      break;
    }
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(8,"Matrix Vector Multiplications",ok);
}

//
//------------------------------------------------------------------------
//           Test operations with vectors and sparse matrix slices
//
void spstress_matrix_slices(Int_t vsize)
{
  if (gVerbose)
    std::cout << "\n---> Test operations with vectors and sparse matrix slices" << std::endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 8.625;

  TVectorD vc(0,vsize);
  TVectorD vr(0,vsize+1);
  TMatrixD       m_d(0,vsize,0,vsize+1); m_d = pattern;
  TMatrixDSparse m(m_d);

  Int_t i,j;
  if (gVerbose)
    std::cout << "\nCheck modifying the matrix row-by-row" << std::endl;
  m = pattern;
  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
    TMatrixDSparseRow(m,i) = pattern+2;
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    vr = TMatrixDSparseRow(m,i);
    ok &= VerifyVectorValue(vr,pattern+2,gVerbose,EPSILON);
    vr = TMatrixDSparseRow(m,i+1 > m.GetRowUpb() ? m.GetRowLwb() : i+1);
    ok &= VerifyVectorValue(vr,pattern,gVerbose,EPSILON);
    TMatrixDSparseRow(m,i) += -2;
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
    vr = pattern-2;
    TMatrixDSparseRow(m,i) = vr;
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    {
      TMatrixDSparseRow mr(m,i);
      for (j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        mr(j) *= 8;
    }
    vr = TMatrixDSparseRow(m,i);
    ok &= VerifyVectorValue(vr,8*(pattern-2),gVerbose,EPSILON);
    TMatrixDSparseRow(m,i) *= 1./8;
    TMatrixDSparseRow(m,i) += 2;
    vr = TMatrixDSparseRow(m,i);
    ok &= VerifyVectorValue(vr,pattern,gVerbose,EPSILON);
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    std::cout << "\nCheck modifying the matrix diagonal" << std::endl;
  m = pattern;
  TMatrixDSparseDiag td = m;
  td = pattern-3;
  ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
  vc = TMatrixDSparseDiag(m);
  ok &= VerifyVectorValue(vc,pattern-3,gVerbose,EPSILON);
  td += 3;
  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  vc = pattern+3;
  td = vc;
  ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
  {
    TMatrixDSparseDiag md(m);
    for (j = 0; j < md.GetNdiags(); j++)
      md(j) /= 1.5;
  }
  vc = TMatrixDSparseDiag(m);
  ok &= VerifyVectorValue(vc,(pattern+3)/1.5,gVerbose,EPSILON);
  TMatrixDSparseDiag(m) *= 1.5;
  TMatrixDSparseDiag(m) += -3;
  vc = TMatrixDSparseDiag(m);
  ok &= VerifyVectorValue(vc,pattern,gVerbose,EPSILON);
  ok &= ( m == pattern ) ? kTRUE : kFALSE;

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;
  StatusPrint(9,"Matrix Slices to Vectors",ok);
}

//
//------------------------------------------------------------------------
//           Test matrix I/O
//
void spstress_matrix_io()
{
  if (gVerbose)
    std::cout << "\n---> Test matrix I/O" << std::endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = TMath::Pi();

  TFile *f = new TFile("stress-vmatrix.root", "RECREATE");

  Char_t name[80];
  Int_t iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TMatrixD m_d(msize,msize); m_d = pattern;
    TMatrixDSparse m(m_d);
    TMatrixDSparse ma;
    ma.Use(m);

    if (verbose)
      std::cout << "\nWrite matrix m to database" << std::endl;
    snprintf(name,80,"m_%d",msize);
    m.Write(name);

    if (verbose)
      std::cout << "\nWrite matrix ma which adopts to database" << std::endl;
    snprintf(name,80,"ma_%d",msize);
    ma.Write(name);

    iloop--;
  }

  if (gVerbose)
    std::cout << "\nClose database" << std::endl;
  delete f;

  if (gVerbose)
    std::cout << "\nOpen database in read-only mode and read matrix" << std::endl;
  TFile *f1 = new TFile("stress-vmatrix.root");

  iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TMatrixD m_d(msize,msize); m_d = pattern;
    TMatrixDSparse m(m_d);
    snprintf(name,80,"m_%d",msize);
    TMatrixDSparse *mr  = (TMatrixDSparse*) f1->Get(name);
    snprintf(name,80,"ma_%d",msize);
    TMatrixDSparse *mar = (TMatrixDSparse*) f1->Get(name);
    snprintf(name,80,"ms_%d",msize);

    if (verbose)
      std::cout << "\nRead matrix should be same as original" << std::endl;
    ok &= ((*mr)  == m) ? kTRUE : kFALSE;
    ok &= ((*mar) == m) ? kTRUE : kFALSE;

    iloop--;
  }

  delete f1;

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(10,"Matrix Persistence",ok);
}

//------------------------------------------------------------------------
//          Test allocation functions and compatibility check
//
void vstress_allocation(Int_t msize)
{
  if (gVerbose)
    std::cout << "\n\n---> Test allocation and compatibility check" << std::endl;

  Int_t i;
  Bool_t ok = kTRUE;
  TVectorD v1(msize);
  TVectorD v2(0,msize-1);
  TVectorD v3(1,msize);
  TVectorD v4(v1);

  if (gVerbose) {
    std::cout << "\nStatus information reported for vector v3:" << std::endl;
    std::cout << "  Lower bound ... " << v3.GetLwb() << std::endl;
    std::cout << "  Upper bound ... " << v3.GetUpb() << std::endl;
    std::cout << "  No. of elements " << v3.GetNoElements() << std::endl;
  }

  if (gVerbose)
    std::cout << "\nCheck vectors 1 & 2 for compatibility" << std::endl;
  ok &= AreCompatible(v1,v2,gVerbose);

  if (gVerbose)
    std::cout << "Check vectors 1 & 4 for compatibility" << std::endl;
  ok &= AreCompatible(v1,v4,gVerbose);

  if (gVerbose)
    std::cout << "v2 has to be compatible with v3 after resizing to v3" << std::endl;
  v2.ResizeTo(v3);
  ok &= AreCompatible(v2,v3,gVerbose);

  TVectorD v5(v1.GetUpb()+5);
  if (gVerbose)
    std::cout << "v1 has to be compatible with v5 after resizing to v5.upb" << std::endl;
  v1.ResizeTo(v5.GetNoElements());
  ok &= AreCompatible(v1,v5,gVerbose);

  {
    if (gVerbose)
      std::cout << "Check that shrinking does not change remaining elements" << std::endl;
    TVectorD vb(-1,msize);
    for (i = vb.GetLwb(); i <= vb.GetUpb(); i++)
      vb(i) = i+TMath::Pi();
    TVectorD v = vb;
    ok &= ( v == vb ) ? kTRUE : kFALSE;
    ok &= ( v != 0 ) ? kTRUE : kFALSE;
    v.ResizeTo(0,msize/2);
    for (i = v.GetLwb(); i <= v.GetUpb(); i++)
      ok &= ( v(i) == vb(i) ) ? kTRUE : kFALSE;
    if (gVerbose)
      std::cout << "Check that expansion expands by zeros" << std::endl;
    const Int_t old_nelems = v.GetNoElements();
    const Int_t old_lwb    = v.GetLwb();
    v.ResizeTo(vb);
    ok &= ( !(v == vb) ) ? kTRUE : kFALSE;
    for (i = old_lwb; i < old_lwb+old_nelems; i++)
      ok &= ( v(i) == vb(i) ) ? kTRUE : kFALSE;
    for (i = v.GetLwb(); i < old_lwb; i++)
      ok &= ( v(i) == 0 ) ? kTRUE : kFALSE;
    for (i = old_lwb+old_nelems; i <= v.GetUpb(); i++)
      ok &= ( v(i) == 0 ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;
  StatusPrint(1,"Allocation, Filling, Resizing",ok);
}

//
//------------------------------------------------------------------------
//                Test uniform element operations
//
class SinAction : public TElementActionD {
  void Operation(Double_t &element) const { element = TMath::Sin(element); }
  public:
    SinAction() { }
};

class CosAction : public TElementPosActionD {
  Double_t factor;
  void Operation(Double_t &element) const { element = TMath::Cos(factor*fI); }
  public:
    CosAction(Int_t no_elems): factor(2*TMath::Pi()/no_elems) { }
};

void vstress_element_op(Int_t vsize)
{
  if (gVerbose)
    std::cout << "\n---> Test operations that treat each element uniformly" << std::endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = TMath::Pi();

  TVectorD v(-1,vsize-2);
  TVectorD v1(v);

  if (gVerbose)
    std::cout << "\nWriting zeros to v..." << std::endl;
  for (Int_t i = v.GetLwb(); i <= v.GetUpb(); i++)
    v(i) = 0;
  ok &= VerifyVectorValue(v,0.0,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "Clearing v1 ..." << std::endl;
  v1.Zero();
  ok &= VerifyVectorValue(v1,0.0,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "Comparing v1 with 0 ..." << std::endl;
  ok &= (v1 == 0) ? kTRUE : kFALSE;

  if (gVerbose)
    std::cout << "Writing a pattern " << pattern << " by assigning to v(i)..." << std::endl;
  {
    for (Int_t i = v.GetLwb(); i <= v.GetUpb(); i++)
      v(i) = pattern;
    ok &= VerifyVectorValue(v,pattern,gVerbose,EPSILON);
  }

  if (gVerbose)
    std::cout << "Writing the pattern by assigning to v1 as a whole ..." << std::endl;
  v1 = pattern;
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "Comparing v and v1 ..." << std::endl;
  ok &= (v == v1) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "Comparing (v=0) and v1 ..." << std::endl;
  ok &= (!(v.Zero() == v1)) ? kTRUE : kFALSE;

  if (gVerbose)
    std::cout << "\nClear v and add the pattern" << std::endl;
  v.Zero();
  v += pattern;
  ok &= VerifyVectorValue(v,pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "   add the doubled pattern with the negative sign" << std::endl;
  v += -2*pattern;
  ok &= VerifyVectorValue(v,-pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "   subtract the trippled pattern with the negative sign" << std::endl;
  v -= -3*pattern;
  ok &= VerifyVectorValue(v,2*pattern,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nVerify comparison operations" << std::endl;
  v = pattern;
  ok &= ( v == pattern && !(v != pattern) && v >= pattern && v <= pattern ) ? kTRUE : kFALSE;
  ok &= ( v > 0 && v >= 0 ) ? kTRUE : kFALSE;
  ok &= ( v > -pattern && v >= -pattern ) ? kTRUE : kFALSE;
  ok &= ( v < pattern+1 && v <= pattern+1 ) ? kTRUE : kFALSE;
  v(v.GetUpb()) += 1;
  ok &= ( !(v==pattern)      && !(v != pattern)  && v != pattern-1 ) ? kTRUE : kFALSE;
  ok &= ( v >= pattern       && !(v > pattern)   && !(v >= pattern+1) ) ? kTRUE : kFALSE;
  ok &= ( v <= pattern+1.001 && !(v < pattern+1) && !(v <= pattern) ) ? kTRUE : kFALSE;

  if (gVerbose)
    std::cout << "\nAssign 2*pattern to v by repeating additions" << std::endl;
  v = 0; v += pattern; v += pattern;
  if (gVerbose)
    std::cout << "Assign 2*pattern to v1 by multiplying by two" << std::endl;
  v1 = pattern; v1 *= 2;
  ok &= VerifyVectorValue(v1,2*pattern,gVerbose,EPSILON);
  ok &= ( v == v1 ) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "Multiply v1 by one half returning it to the 1*pattern" << std::endl;
  v1 *= 1/2.;
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nAssign -pattern to v and v1" << std::endl;
  v.Zero(); v -= pattern; v1 = -pattern;
  ok &= VerifyVectorValue(v,-pattern,gVerbose,EPSILON);
  ok &= ( v == v1 ) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "v = sqrt(sqr(v)); v1 = abs(v1); Now v and v1 have to be the same" << std::endl;
  v.Sqr();
  ok &= VerifyVectorValue(v,pattern*pattern,gVerbose,EPSILON);
  v.Sqrt();
  ok &= VerifyVectorValue(v,pattern,gVerbose,EPSILON);
  v1.Abs();
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);
  ok &= ( v == v1 ) ? kTRUE : kFALSE;

  {
    if (gVerbose)
      std::cout << "\nCheck out to see that sin^2(x) + cos^2(x) = 1" << std::endl;
    for (Int_t i = v.GetLwb(); i <= v.GetUpb(); i++)
      v(i) = 2*TMath::Pi()/v.GetNoElements()*i;
#ifndef __CINT__
    SinAction s;
    v.Apply(s);
    CosAction c(v.GetNoElements());
    v1.Apply(c);
#else
    for (Int_t i = v.GetLwb(); i <= v.GetUpb(); i++) {
      v(i)  = TMath::Sin(v(i));
      v1(i) = TMath::Cos(2*TMath::Pi()/v1.GetNrows()*i);
    }
#endif
    TVectorD v2 = v;
    TVectorD v3 = v1;
    v.Sqr();
    v1.Sqr();
    v += v1;
    ok &= VerifyVectorValue(v,1.,gVerbose,EPSILON);
  }

  if (gVerbose)
    std::cout << "\nVerify constructor with initialization" << std::endl;
#ifndef __CINT__
  TVectorD vi(0,4,0.0,1.0,2.0,3.0,4.0,"END");
#else
  Double_t vval[] = {0.0,1.0,2.0,3.0,4.0};
  TVectorD vi(5,vval);
#endif
  TVectorD vit(5);
  {
    for (Int_t i = vit.GetLwb(); i <= vit.GetUpb(); i++)
      vit(i) = Double_t(i);
    ok &= VerifyVectorIdentity(vi,vit,gVerbose,EPSILON);
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;
  StatusPrint(2,"Uniform vector operations",ok);
}

//
//------------------------------------------------------------------------
//                 Test binary vector operations
//
void vstress_binary_op(Int_t vsize)
{
  if (gVerbose)
    std::cout << "\n---> Test Binary Vector operations" << std::endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = TMath::Pi();

  const Double_t epsilon = EPSILON*vsize/10;

  TVectorD v(2,vsize+1);
  TVectorD v1(v);

  if (gVerbose)
    std::cout << "\nVerify assignment of a vector to the vector" << std::endl;
  v = pattern;
  v1.Zero();
  v1 = v;
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);
  ok &= ( v1 == v ) ? kTRUE : kFALSE;

  if (gVerbose)
    std::cout << "\nAdding one vector to itself, uniform pattern " << pattern << std::endl;
  v.Zero(); v = pattern;
  v1 = v; v1 += v1;
  ok &= VerifyVectorValue(v1,2*pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "  subtracting two vectors ..." << std::endl;
  v1 -= v;
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "  subtracting the vector from itself" << std::endl;
  // Hiding the self-subtraction from the compiler, causing warning by -Wself-assign-overloaded in clang 7.0
  v1 -= static_cast<TVectorD&>(v1);
  ok &= VerifyVectorValue(v1,0.,gVerbose,EPSILON);
  if (gVerbose)
    std::cout << "  adding two vectors together" << std::endl;
  v1 += v;
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);

  TVectorD vp(2,vsize+1);
  {
    for (Int_t i = vp.GetLwb(); i <= vp.GetUpb(); i++)
      vp(i) = (i-vp.GetNoElements()/2.)*pattern;
  }

  if (gVerbose) {
    std::cout << "\nArithmetic operations on vectors with not the same elements" << std::endl;
    std::cout << "   adding vp to the zero vector..." << std::endl;
  }
  v.Zero();
  ok &= ( v == 0.0 ) ? kTRUE : kFALSE;
  v += vp;
  ok &= VerifyVectorIdentity(v,vp,gVerbose,epsilon);
  v1 = v;
  if (gVerbose)
    std::cout << "   making v = 3*vp and v1 = 3*vp, via add() and succesive mult" << std::endl;
  Add(v,2.,vp);
  v1 += v1; v1 += vp;
  ok &= VerifyVectorIdentity(v,v1,gVerbose,epsilon);
  if (gVerbose)
    std::cout << "   clear both v and v1, by subtracting from itself and via add()" << std::endl;
  // Hiding the self-subtraction from the compiler, causing warning by -Wself-assign-overloaded in clang 7.0
  v1 -= static_cast<TVectorD&>(v1);
  Add(v,-3.,vp);
  ok &= VerifyVectorIdentity(v,v1,gVerbose,epsilon);

  if (gVerbose) {
    std::cout << "\nTesting element-by-element multiplications and divisions" << std::endl;
    std::cout << "   squaring each element with sqr() and via multiplication" << std::endl;
  }
  v = vp; v1 = vp;
  v.Sqr();
  ElementMult(v1,v1);
  ok &= VerifyVectorIdentity(v,v1,gVerbose,epsilon);
  if (gVerbose)
    std::cout << "   compare (v = pattern^2)/pattern with pattern" << std::endl;
  v = pattern; v1 = pattern;
  v.Sqr();
  ElementDiv(v,v1);
  ok &= VerifyVectorValue(v,pattern,gVerbose,epsilon);
  if (gVerbose)
   Compare(v1,v);

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;
  StatusPrint(3,"Binary vector element-by-element operations",ok);
}

//
//------------------------------------------------------------------------
//               Verify the norm calculation
//
void vstress_norms(Int_t vsize)
{
  if (gVerbose)
    std::cout << "\n---> Verify norm calculations" << std::endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 10.25;

  if ( vsize % 2 == 1 )
    Fatal("vstress_norms", "size of the vector to test must be even for this test\n");

  TVectorD v(vsize);
  TVectorD v1(v);

  if (gVerbose)
    std::cout << "\nAssign " << pattern << " to all the elements and check norms" << std::endl;
  v = pattern;
  if (gVerbose)
    std::cout << "  1. norm should be pattern*no_elems" << std::endl;
  ok &= ( v.Norm1() == pattern*v.GetNoElements() ) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "  Square of the 2. norm has got to be pattern^2 * no_elems" << std::endl;
  ok &= ( v.Norm2Sqr() == (pattern*pattern)*v.GetNoElements() ) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "  Inf norm should be pattern itself" << std::endl;
  ok &= ( v.NormInf() == pattern ) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "  Scalar product of vector by itself is the sqr(2. vector norm)" << std::endl;
  ok &= ( v.Norm2Sqr() == v*v ) ? kTRUE : kFALSE;

  Double_t ap_step = 1;
  Double_t ap_a0   = -pattern;
  Int_t n = v.GetNoElements();
  if (gVerbose) {
    std::cout << "\nAssign the arithm progression with 1. term " << ap_a0 <<
            "\nand the difference " << ap_step << std::endl;
  }
  {
    for (Int_t i = v.GetLwb(); i <= v.GetUpb(); i++)
      v(i) = (i-v.GetLwb())*ap_step + ap_a0;
  }
  Int_t l = TMath::Min(TMath::Max((int)TMath::Ceil(-ap_a0/ap_step),0),n);
  Double_t norm = (2*ap_a0+(l+n-1)*ap_step)/2*(n-l) +
                  (-2*ap_a0-(l-1)*ap_step)/2*l;
  if (gVerbose)
    std::cout << "  1. norm should be " << norm << std::endl;
  ok &= ( v.Norm1() == norm ) ? kTRUE : kFALSE;
  norm = n*( (ap_a0*ap_a0)+ap_a0*ap_step*(n-1)+(ap_step*ap_step)*(n-1)*(2*n-1)/6);
  if (gVerbose) {
    std::cout << "  Square of the 2. norm has got to be "
            "n*[ a0^2 + a0*q*(n-1) + q^2/6*(n-1)*(2n-1) ], or " << norm << std::endl;
  }
  ok &= ( TMath::Abs( (v.Norm2Sqr()-norm)/norm ) < EPSILON ) ? kTRUE : kFALSE;

  norm = TMath::Max(TMath::Abs(v(v.GetLwb())),TMath::Abs(v(v.GetUpb())));
  if (gVerbose)
    std::cout << "  Inf norm should be max(abs(a0),abs(a0+(n-1)*q)), ie " << norm << std::endl;
  ok &= ( v.NormInf() == norm ) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "  Scalar product of vector by itself is the sqr(2. vector norm)" << std::endl;
  ok &= ( v.Norm2Sqr() == v*v ) ? kTRUE : kFALSE;

#if 0
  v1.Zero();
  Compare(v,v1);  // they are not equal (of course)
#endif

  if (gVerbose)
    std::cout << "\nConstruct v1 to be orthogonal to v as v(n), -v(n-1), v(n-2)..." << std::endl;
  {
    for (Int_t i = 0; i < v1.GetNoElements(); i++)
      v1(i+v1.GetLwb()) = v(v.GetUpb()-i) * ( i % 2 == 1 ? -1 : 1 );
  }
  if (gVerbose)
    std::cout << "||v1|| has got to be equal ||v|| regardless of the norm def" << std::endl;
  ok &= ( v1.Norm1()    == v.Norm1() ) ? kTRUE : kFALSE;
  ok &= ( v1.Norm2Sqr() == v.Norm2Sqr() ) ? kTRUE : kFALSE;
  ok &= ( v1.NormInf()  == v.NormInf() ) ? kTRUE : kFALSE;
  if (gVerbose)
    std::cout << "But the scalar product has to be zero" << std::endl;
  ok &= ( v1 * v == 0 ) ? kTRUE : kFALSE;

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;
  StatusPrint(4,"Vector Norms",ok);
}

//
//------------------------------------------------------------------------
//           Test operations with vectors and matrix slices
//
void vstress_matrix_slices(Int_t vsize)
{
  if (gVerbose)
    std::cout << "\n---> Test operations with vectors and matrix slices" << std::endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 8.625;

  TVectorD vc(0,vsize);
  TVectorD vr(0,vsize+1);
  TMatrixD m(0,vsize,0,vsize+1);

  Int_t i,j;
  if (gVerbose)
    std::cout << "\nCheck modifying the matrix column-by-column" << std::endl;
  m = pattern;
  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetColLwb(); i <= m.GetColUpb(); i++) {
    TMatrixDColumn(m,i) = pattern-1;
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    TMatrixDColumn(m,i) *= 2;
    vc = TMatrixDColumn(m,i);
    ok &= VerifyVectorValue(vc,2*(pattern-1),gVerbose);
    vc = TMatrixDColumn(m,i+1 > m.GetColUpb() ? m.GetColLwb() : i+1);
    ok &= VerifyVectorValue(vc,pattern,gVerbose,EPSILON);
    TMatrixDColumn(m,i) *= 0.5;
    TMatrixDColumn(m,i) += 1;
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetColLwb(); i <= m.GetColUpb(); i++) {
    vc = pattern+1;
    TMatrixDColumn(m,i) = vc;
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    {
      TMatrixDColumn mc(m,i);
      for (j = m.GetRowLwb(); j <= m.GetRowUpb(); j++)
        mc(j) *= 4;
    }
    vc = TMatrixDColumn(m,i);
    ok &= VerifyVectorValue(vc,4*(pattern+1),gVerbose,EPSILON);
    TMatrixDColumn(m,i) *= 0.25;
    TMatrixDColumn(m,i) += -1;
    vc = TMatrixDColumn(m,i);
    ok &= VerifyVectorValue(vc,pattern,gVerbose,EPSILON);
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    std::cout << "\nCheck modifying the matrix row-by-row" << std::endl;
  m = pattern;
  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
     TMatrixDRow(m,i).Assign(pattern+2);
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    vr = TMatrixDRow(m,i);
    ok &= VerifyVectorValue(vr,pattern+2,gVerbose,EPSILON);
    vr = TMatrixDRow(m,i+1 > m.GetRowUpb() ? m.GetRowLwb() : i+1);
    ok &= VerifyVectorValue(vr,pattern,gVerbose,EPSILON);
    TMatrixDRow(m,i) += -2;
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
    vr = pattern-2;
    TMatrixDRow(m,i) = vr;
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    {
      TMatrixDRow mr(m,i);
      for (j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        mr(j) *= 8;
    }
    vr = TMatrixDRow(m,i);
    ok &= VerifyVectorValue(vr,8*(pattern-2),gVerbose,EPSILON);
    TMatrixDRow(m,i) *= 1./8;
    TMatrixDRow(m,i) += 2;
    vr = TMatrixDRow(m,i);
    ok &= VerifyVectorValue(vr,pattern,gVerbose,EPSILON);
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    std::cout << "\nCheck modifying the matrix diagonal" << std::endl;
  m = pattern;
  TMatrixDDiag td = m;
  td = pattern-3;
  ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
  vc = TMatrixDDiag(m);
  ok &= VerifyVectorValue(vc,pattern-3,gVerbose,EPSILON);
  td += 3;
  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  vc = pattern+3;
  td = vc;
  ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
  {
    TMatrixDDiag md(m);
    for (j = 0; j < md.GetNdiags(); j++)
      md(j) /= 1.5;
  }
  vc = TMatrixDDiag(m);
  ok &= VerifyVectorValue(vc,(pattern+3)/1.5,gVerbose,EPSILON);
  TMatrixDDiag(m) *= 1.5;
  TMatrixDDiag(m) += -3;
  vc = TMatrixDDiag(m);
  ok &= VerifyVectorValue(vc,pattern,gVerbose,EPSILON);
  ok &= ( m == pattern ) ? kTRUE : kFALSE;

  if (gVerbose) {
    std::cout << "\nCheck out to see that multiplying by diagonal is column-wise"
            "\nmatrix multiplication" << std::endl;
  }
  TMatrixD mm(m);
  TMatrixD m1(m.GetRowLwb(),TMath::Max(m.GetRowUpb(),m.GetColUpb()),
              m.GetColLwb(),TMath::Max(m.GetRowUpb(),m.GetColUpb()));
  TVectorD vc1(vc),vc2(vc);
  for (i = m.GetRowLwb(); i < m.GetRowUpb(); i++)
     TMatrixDRow(m,i).Assign(pattern+i);      // Make a multiplicand
  mm = m;                          // Save it

  m1 = pattern+10;
  for (i = vr.GetLwb(); i <= vr.GetUpb(); i++)
    vr(i) = i+2;
  TMatrixDDiag td2 = m1;
  td2 = vr;
  ok &= ( !(m1 == pattern+10) ) ? kTRUE : kFALSE;

  m *= TMatrixDDiag(m1);
  for (i = m.GetColLwb(); i <= m.GetColUpb(); i++) {
    vc1 = TMatrixDColumn(mm,i);
    vc1 *= vr(i);                    // Do a column-wise multiplication
    vc2 = TMatrixDColumn(m,i);
    ok &= VerifyVectorIdentity(vc1,vc2,gVerbose,EPSILON);
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;
  StatusPrint(5,"Matrix Slices to Vectors",ok);
}

//
//------------------------------------------------------------------------
//           Test vector I/O
//
void vstress_vector_io()
{
  if (gVerbose)
    std::cout << "\n---> Test vector I/O" << std::endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = TMath::Pi();

  TFile *f = new TFile("stress-vvector.root","RECREATE");

  Char_t name[80];
  Int_t iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TVectorD v(msize);
    v = pattern;

    Double_t *pattern_array = new Double_t[msize];
    for (Int_t i = 0; i < msize; i++)
      pattern_array[i] = pattern;
    TVectorD va;
    va.Use(msize,pattern_array);

    if (verbose)
      std::cout << "\nWrite vector v to database" << std::endl;

    snprintf(name,80,"v_%d",msize);
    v.Write(name);
    snprintf(name,80,"va_%d",msize);
    va.Write(name);

    delete [] pattern_array;

    iloop--;
  }

  if (gVerbose)
    std::cout << "\nClose database" << std::endl;
  delete f;

  if (gVerbose)
    std::cout << "\nOpen database in read-only mode and read vector" << std::endl;
  TFile *f1 = new TFile("stress-vvector.root");

  iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TVectorD v(msize);
    v = pattern;

    snprintf(name,80,"v_%d",msize);
    TVectorD *vr  = (TVectorD*) f1->Get(name);
    snprintf(name,80,"va_%d",msize);
    TVectorD *var = (TVectorD*) f1->Get(name);

    if (verbose)
      std::cout << "\nRead vector should be same as original still in memory" << std::endl;
    ok &= ((*vr) == v)  ? kTRUE : kFALSE;
    ok &= ((*var) == v) ? kTRUE : kFALSE;

    iloop--;
  }

  delete f1;

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;
  StatusPrint(6,"Vector Persistence",ok);
}

Bool_t test_svd_expansion(const TMatrixD &A)
{
  if (gVerbose)
    std::cout << "\n\nSVD-decompose matrix A and check if we can compose it back\n" <<std::endl;

  Bool_t ok = kTRUE;

  TDecompSVD svd(A);
  if (gVerbose) {
    std::cout << "left factor U" <<std::endl;
    svd.GetU().Print();
    std::cout << "Vector of Singular values" <<std::endl;
    svd.GetSig().Print();
    std::cout << "right factor V" <<std::endl;
    svd.GetV().Print();
  }

  {
    if (gVerbose)
      std::cout << "\tchecking that U is orthogonal indeed, i.e., U'U=E and UU'=E" <<std::endl;
    const Int_t nRows = svd.GetU().GetNrows();
    const Int_t nCols = svd.GetU().GetNcols();
    TMatrixD E1(nRows,nRows); E1.UnitMatrix();
    TMatrixD E2(nCols,nCols); E2.UnitMatrix();
    TMatrixD ut(TMatrixD::kTransposed,svd.GetU());
    ok &= VerifyMatrixIdentity(ut * svd.GetU(),E2,gVerbose,100*EPSILON);
    ok &= VerifyMatrixIdentity(svd.GetU() * ut,E1,gVerbose,100*EPSILON);
  }

  {
    if (gVerbose)
      std::cout << "\tchecking that V is orthogonal indeed, i.e., V'V=E and VV'=E" <<std::endl;
    const Int_t nRows = svd.GetV().GetNrows();
    const Int_t nCols = svd.GetV().GetNcols();
    TMatrixD E1(nRows,nRows); E1.UnitMatrix();
    TMatrixD E2(nCols,nCols); E2.UnitMatrix();
    TMatrixD vt(TMatrixD::kTransposed,svd.GetV());
    ok &= VerifyMatrixIdentity(vt * svd.GetV(),E2,gVerbose,100*EPSILON);
    ok &= VerifyMatrixIdentity(svd.GetV() * vt,E1,gVerbose,100*EPSILON);
  }

  {
    if (gVerbose)
      std::cout << "\tchecking that U*Sig*V' is indeed A" <<std::endl;
    const Int_t nRows = svd.GetU().GetNrows();
    const Int_t nCols = svd.GetV().GetNcols();
    TMatrixD s(nRows,nCols);
    TMatrixDDiag diag(s); diag = svd.GetSig();
    TMatrixD vt(TMatrixD::kTransposed,svd.GetV());
    TMatrixD tmp = s * vt;
    ok &= VerifyMatrixIdentity(A,svd.GetU() * tmp,gVerbose,100*EPSILON);
    if (gVerbose) {
      std::cout << "U*Sig*V'" <<std::endl;
      (svd.GetU()*tmp).Print();
    }
  }

  return ok;
}

#ifndef __CINT__
// Make a matrix from an array (read it row-by-row)
class MakeMatrix : public TMatrixDLazy {
  const Double_t *array;
        Int_t     no_elems;
  void FillIn(TMatrixD& m) const {
    R__ASSERT( m.GetNrows() * m.GetNcols() == no_elems );
    const Double_t *ap = array;
          Double_t *mp = m.GetMatrixArray();
    for (Int_t i = 0; i < no_elems; i++)
      *mp++ = *ap++;
  }

public:
  MakeMatrix(Int_t nrows,Int_t ncols,
             const Double_t *_array,Int_t _no_elems)
    :TMatrixDLazy(nrows,ncols), array(_array), no_elems(_no_elems) {}
  MakeMatrix(Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb,
             const Double_t *_array,Int_t _no_elems)
    : TMatrixDLazy(row_lwb,row_upb,col_lwb,col_upb),
      array(_array), no_elems(_no_elems) {}
};
#else
TMatrixD MakeMatrix(Int_t nrows,Int_t ncols,
                    const Double_t *_array,Int_t _no_elems)
{
  TMatrixD m(nrows,ncols,_array);
  return m;
}
#endif

void astress_decomp()
{
  Bool_t ok = kTRUE;

  {
    if (gVerbose)
      std::cout << "\nBrandt example page 503\n" <<std::endl;

    Double_t array0[] = { -2,1,0,0, 2,1,0,0, 0,0,0,0, 0,0,0,0 };
    ok &= test_svd_expansion(MakeMatrix(4,4,array0,sizeof(array0)/sizeof(array0[0])));
  }

  {
    if (gVerbose)
      std::cout << "\nRotated by PI/2 Matrix Diag(1,4,9)\n" <<std::endl;

    Double_t array1[] = {0,-4,0,  1,0,0,  0,0,9 };
    ok &= test_svd_expansion(MakeMatrix(3,3,array1,sizeof(array1)/sizeof(array1[0])));
  }

  {
    if (gVerbose)
      std::cout << "\nExample from the Forsythe, Malcolm, Moler's book\n" <<std::endl;

    Double_t array2[] =
         { 1,6,11, 2,7,12, 3,8,13, 4,9,14, 5,10,15};
    ok &= test_svd_expansion(MakeMatrix(5,3,array2,sizeof(array2)/sizeof(array2[0])));
  }

  {
    if (gVerbose)
      std::cout << "\nExample from the Wilkinson, Reinsch's book\n" <<
              "Singular numbers are 0, 19.5959, 20, 0, 35.3270\n" <<std::endl;

    Double_t array3[] =
        { 22, 10,  2,   3,  7,    14,  7, 10,  0,  8,
          -1, 13, -1, -11,  3,    -3, -2, 13, -2,  4,
           9,  8,  1,  -2,  4,     9,  1, -7,  5, -1,
           2, -6,  6,   5,  1,     4,  5,  0, -2,  2 };
    ok &= test_svd_expansion(MakeMatrix(8,5,array3,sizeof(array3)/sizeof(array3[0])));
  }

  {
    if (gVerbose)
      std::cout << "\nExample from the Wilkinson, Reinsch's book\n" <<
              "Ordered singular numbers are Sig[21-k] = sqrt(k*(k-1))\n" <<std::endl;
    TMatrixD A(21,20);
    for (Int_t irow = A.GetRowLwb(); irow <= A.GetRowUpb(); irow++)
      for (Int_t icol = A.GetColLwb(); icol <= A.GetColUpb(); icol++)
        A(irow,icol) = (irow>icol ? 0 : ((irow==icol) ? 20-icol : -1));

    ok &= test_svd_expansion(A);
  }

  if (0)
  {
    if (gVerbose) {
      std::cout << "\nTest by W. Meier <wmeier@manu.com> to catch an obscure "
           << "bug in QR\n" <<std::endl;
      std::cout << "expect singular values to be\n"
           << "1.4666e-024   1.828427   3.828427   4.366725  7.932951\n" <<std::endl;
    }

    Double_t array4[] =
        {  1,  2,  0,  0,  0,
           0,  2,  3,  0,  0,
           0,  0,  0,  4,  0,
           0,  0,  0,  4,  5,
           0,  0,  0,  0,  5 };
    ok &= test_svd_expansion(MakeMatrix(5,5,array4,sizeof(array4)/sizeof(array4[0])));
  }

  {
    const TMatrixD m1 = THilbertMatrixD(5,5);
    TDecompLU lu(m1);
    ok &= VerifyMatrixIdentity(lu.GetMatrix(),m1,gVerbose,100*EPSILON);
  }

  {
    const TMatrixD m2 = THilbertMatrixD(5,5);
    const TMatrixDSym mtm(TMatrixDSym::kAtA,m2);
    TDecompChol chol(mtm);
    ok &= VerifyMatrixIdentity(chol.GetMatrix(),mtm,gVerbose,100*EPSILON);
  }

  if (gVerbose)
    std::cout << "\nDone" <<std::endl;

  StatusPrint(1,"Decomposition / Reconstruction",ok);
}

void astress_lineqn()
{
  if (gVerbose)
    std::cout << "\nSolve Ax=b where A is a Hilbert matrix and b(i) = sum_j Aij\n" <<std::endl;

  Bool_t ok = kTRUE;

  Int_t iloop = gNrLoop;
  Int_t nr    = 0;

  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    if (verbose)
      std::cout << "\nSolve Ax=b for size = " << msize <<std::endl;

    // Since The Hilbert matrix is accuracy "challenged", I will use a diagonaly
    // dominant one fore sizes > 100, otherwise the verification might fail

    TMatrixDSym m1 = THilbertMatrixDSym(-1,msize-2);
    TMatrixDDiag diag1(m1);
    diag1 += 1.;

    TVectorD rowsum1(-1,msize-2); rowsum1.Zero();
    TVectorD colsum1(-1,msize-2); colsum1.Zero();
    for (Int_t irow = m1.GetRowLwb(); irow <= m1.GetColUpb(); irow++) {
      for (Int_t icol = m1.GetColLwb(); icol <= m1.GetColUpb(); icol++) {
        rowsum1(irow) += m1(irow,icol);
        colsum1(icol) += m1(irow,icol);
      }
    }

    TMatrixDSym m2 = THilbertMatrixDSym(msize);
    TMatrixDDiag diag2(m2);
    diag2 += 1.;

    TVectorD rowsum2(msize); rowsum2.Zero();
    TVectorD colsum2(msize); colsum2.Zero();
    for (Int_t irow = m2.GetRowLwb(); irow <= m2.GetColUpb(); irow++) {
      for (Int_t icol = m2.GetColLwb(); icol <= m2.GetColUpb(); icol++) {
        rowsum2(irow) += m2(irow,icol);
        colsum2(icol) += m2(irow,icol);
      }
    }

    TVectorD b1(-1,msize-2);
    TVectorD b2(msize);
    {
      TDecompLU lu(m1,1.0e-20);
      b1 = rowsum1;
      lu.Solve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
      b1 = colsum1;
      lu.TransSolve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
    }

    {
      TDecompChol chol(m1,1.0e-20);
      b1 = rowsum1;
      chol.Solve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
      b1 = colsum1;
      chol.TransSolve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
    }

    {
      TDecompQRH qrh1(m1,1.0e-20);
      b1 = rowsum1;
      qrh1.Solve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);

      TDecompQRH qrh2(m2,1.0e-20);
      b2 = colsum2;
      qrh2.TransSolve(b2);
      if (msize < 10)
        ok &= VerifyVectorValue(b2,1.0,verbose,msize*EPSILON);
    }

    {
      TDecompSVD svd1(m1);
      b1 = rowsum1;
      svd1.Solve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);

      TDecompSVD svd2(m2);
      b2 = colsum2;
      svd2.TransSolve(b2);
      if (msize < 10)
        ok &= VerifyVectorValue(b2,1.0,verbose,msize*EPSILON);
    }

    {
      TDecompBK bk(m1,1.0e-20);
      b1 = rowsum1;
      bk.Solve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
      b1 = colsum1;
      bk.TransSolve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop--;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        std::cout << "Linear Equations failed for size " << msize << std::endl;
      break;
    }
  }

  if (gVerbose)
    std::cout << "\nDone" <<std::endl;

  StatusPrint(2,"Linear Equations",ok);
}

void astress_pseudo()
{
// The pseudo-inverse of A is found by "inverting" the SVD of A.
// To be more precise, we use SVD to solve the equation
// AX=B where B is a unit matrix.
//
// After we compute the pseudo-inverse, we verify the Moore-Penrose
// conditions: Matrix X is a pseudo-inverse of A iff
//      AXA = A
//      XAX = X
//      XA = (XA)' (i.e., XA is symmetric)
//      AX = (AX)' (i.e., AX is symmetric)

  Bool_t ok = kTRUE;

  // Allocate and fill matrix A
  enum {nrows = 4, ncols = 3};
#ifndef __CINT__
  const Double_t A_data [nrows*ncols] =
#else
  const Double_t A_data [12] =
#endif
   {0, 0, 0,
    0, 0, 0,
    1, 1, 0,
    4, 3, 0};
  TMatrixD A(nrows,ncols,A_data);

  // Perform the SVD decomposition of the transposed matrix.

  TDecompSVD svd(A);
  if (gVerbose)
    std::cout << "\ncondition number " << svd.Condition() << std::endl;

  // Compute the inverse as a solution of A*Ainv = E
  TMatrixD Ainv(nrows,nrows); Ainv.UnitMatrix();
  svd.MultiSolve(Ainv);
  Ainv.ResizeTo(ncols,nrows);
  TMatrixD Ainv2(nrows,nrows);
  svd.Invert(Ainv2);
  Ainv2.ResizeTo(ncols,nrows);
  ok &= VerifyMatrixIdentity(Ainv,Ainv2,gVerbose,EPSILON);

  if (gVerbose) {
    std::cout << "\nChecking the Moore-Penrose conditions for the Ainv" << std::endl;
    std::cout << "\nVerify that Ainv * A * Ainv is Ainv" << std::endl;
  }
  ok &= VerifyMatrixIdentity(Ainv * A * Ainv, Ainv,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nVerify that A * Ainv * A is A" << std::endl;
  ok &= VerifyMatrixIdentity(A * Ainv * A, A,gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nVerify that Ainv * A is symmetric" << std::endl;
  ok &= VerifyMatrixIdentity(Ainv * A, (Ainv*A).T(),gVerbose,EPSILON);

  if (gVerbose)
    std::cout << "\nVerify that A * Ainv is symmetric" << std::endl;
  ok &= VerifyMatrixIdentity(A * Ainv, (A*Ainv).T(),gVerbose,EPSILON);

  StatusPrint(3,"Pseudo-Inverse, Moore-Penrose",ok);
}

void astress_eigen(Int_t msize)
{
  Bool_t ok = kTRUE;

  // Check :
  // 1.  the eigen values of M' * M are the same as the singular values
  //     squared of the SVD of M .
  // 2.  M' * M  x_i =  lambda_i x_i where x_i and lambda_i are the ith
  //     eigen - vector/value .

  const TMatrixD m = THilbertMatrixD(msize,msize);

  TDecompSVD svd(m);
  TVectorD sig = svd.GetSig(); sig.Sqr();

  // Symmetric matrix EigenVector algorithm
  TMatrixDSym mtm1(TMatrixDSym::kAtA,m);
  const TMatrixDSymEigen eigen1(mtm1);
  const TVectorD eigenVal1 = eigen1.GetEigenValues();

  ok &= VerifyVectorIdentity(sig,eigenVal1,gVerbose,EPSILON);

  // Check that the eigen vectors by comparing M'M x to lambda x
  const TMatrixD mtm1x = mtm1 * eigen1.GetEigenVectors();
  TMatrixD lam_x1 = eigen1.GetEigenVectors();
  lam_x1.NormByRow(eigen1.GetEigenValues(),"");

  ok &= VerifyMatrixIdentity(mtm1x,lam_x1,gVerbose,EPSILON);

  // General matrix EigenVector algorithm tested on symmetric matrix
  TMatrixD mtm2(m,TMatrixD::kTransposeMult,m);
  const TMatrixDEigen eigen2(mtm2);
  const TVectorD eigenVal2 = eigen2.GetEigenValuesRe();

  ok &= VerifyVectorIdentity(sig,eigenVal2,gVerbose,EPSILON);

  // Check that the eigen vectors by comparing M'M x to lambda x
  const TMatrixD mtm2x = mtm2 * eigen2.GetEigenVectors();
  TMatrixD lam_x2 = eigen2.GetEigenVectors();
  lam_x2.NormByRow(eigen2.GetEigenValuesRe(),"");

  ok &= VerifyMatrixIdentity(mtm2x,lam_x2,gVerbose,EPSILON);

  // The Imaginary part of the eigenvalues should be zero
  const TVectorD eigenValIm = eigen2.GetEigenValuesIm();
  TVectorD epsVec(msize); epsVec = EPSILON;
  ok &= VerifyVectorIdentity(epsVec,eigenValIm,gVerbose,EPSILON);

  StatusPrint(4,"Eigen - Values/Vectors",ok);
}

void astress_decomp_io(Int_t msize)
{
//
//------------------------------------------------------------------------
//           Test decomposition I/O
//
  if (gVerbose)
    std::cout << "\n---> Test decomp I/O" << std::endl;

  Bool_t ok = kTRUE;

  const TMatrixDSym m = THilbertMatrixDSym(msize);
  TVectorD rowsum(msize); rowsum.Zero();
  TVectorD colsum(msize); colsum.Zero();

  for (Int_t irow = 0; irow < m.GetNrows(); irow++) {
    for (Int_t icol = 0; icol < m.GetNcols(); icol++) {
      rowsum(irow) += m(irow,icol);
      colsum(icol) += m(irow,icol);
    }
  }

  if (gVerbose)
    std::cout << "\nWrite decomp m to database" << std::endl;

  TFile *f = new TFile("stress-vdecomp.root", "RECREATE");

  TDecompLU   lu(m,1.0e-20);
  TDecompQRH  qrh(m,1.0e-20);
  TDecompChol chol(m,1.0e-20);
  TDecompSVD  svd(m);
  TDecompBK   bk(m,1.0e-20);

  lu.Write("lu");
  qrh.Write("qrh");
  chol.Write("chol");
  svd.Write("svd");
  bk.Write("bk");

  if (gVerbose)
    std::cout << "\nClose database" << std::endl;
  delete f;

  if (gVerbose)
    std::cout << "\nOpen database in read-only mode and read matrix" << std::endl;
  TFile *f1 = new TFile("stress-vdecomp.root");

  if (gVerbose)
    std::cout << "\nRead decompositions should create same solutions" << std::endl;
  {
    TDecompLU *rlu = (TDecompLU*)  f1->Get("lu");
    TVectorD b1(rowsum);
    lu.Solve(b1);
    TVectorD b2(rowsum);
    rlu->Solve(b2);
    ok &= (b1 == b2);
    b1 = colsum;
    lu.TransSolve(b1);
    b2 = colsum;
    rlu->TransSolve(b2);
    ok &= (b1 == b2);
  }

  {
    TDecompChol *rchol = (TDecompChol*) f1->Get("chol");
    TVectorD b1(rowsum);
    chol.Solve(b1);
    TVectorD b2(rowsum);
    rchol->Solve(b2);
    ok &= (b1 == b2);
    b1 = colsum;
    chol.TransSolve(b1);
    b2 = colsum;
    rchol->TransSolve(b2);
    ok &= (b1 == b2);
  }

  {
    TDecompQRH *rqrh = (TDecompQRH*) f1->Get("qrh");
    TVectorD b1(rowsum);
    qrh.Solve(b1);
    TVectorD b2(rowsum);
    rqrh->Solve(b2);
    ok &= (b1 == b2);
    b1 = colsum;
    qrh.TransSolve(b1);
    b2 = colsum;
    rqrh->TransSolve(b2);
    ok &= (b1 == b2);
  }

  {
    TDecompSVD *rsvd = (TDecompSVD*) f1->Get("svd");
    TVectorD b1(rowsum);
    svd.Solve(b1);
    TVectorD b2(rowsum);
    rsvd->Solve(b2);
    ok &= (b1 == b2);
    b1 = colsum;
    svd.TransSolve(b1);
    b2 = colsum;
    rsvd->TransSolve(b2);
    ok &= (b1 == b2);
  }

  {
    TDecompBK *rbk = (TDecompBK*) f1->Get("bk");
    TVectorD b1(rowsum);
    bk.Solve(b1);
    TVectorD b2(rowsum);
    rbk->Solve(b2);
    ok &= (b1 == b2);
    b1 = colsum;
    bk.TransSolve(b1);
    b2 = colsum;
    rbk->TransSolve(b2);
    ok &= (b1 == b2);
  }

  delete f1;

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(5,"Decomposition Persistence",ok);
}

void stress_backward_io()
{
  TFile::SetCacheFileDir(".");
  TFile *f = TFile::Open("http://root.cern.ch/files/linearIO.root","CACHEREAD");

  TMatrixF mf1 = THilbertMatrixF(-5,5,-5,5);
  mf1[1][2] = TMath::Pi();
  TMatrixFSym mf2 = THilbertMatrixFSym(-5,5);
  TVectorF vf_row(mf1.GetRowLwb(),mf1.GetRowUpb()); vf_row = TMatrixFRow(mf1,3);

  TMatrixF    *mf1_r    = (TMatrixF*)    f->Get("mf1");
  TMatrixFSym *mf2_r    = (TMatrixFSym*) f->Get("mf2");
  TVectorF    *vf_row_r = (TVectorF*)    f->Get("vf_row");

  Bool_t ok = kTRUE;

  ok &= ((*mf1_r) == mf1) ? kTRUE : kFALSE;
  ok &= ((*mf2_r) == mf2) ? kTRUE : kFALSE;
  ok &= ((*vf_row_r) == vf_row) ? kTRUE : kFALSE;

  TMatrixD md1 = THilbertMatrixD(-5,5,-5,5);
  md1[1][2] = TMath::Pi();
  TMatrixDSym md2 = THilbertMatrixDSym(-5,5);
  TMatrixDSparse md3 = md1;
  TVectorD vd_row(md1.GetRowLwb(),md1.GetRowUpb()); vd_row = TMatrixDRow(md1,3);

  TMatrixD       *md1_r    = (TMatrixD*)       f->Get("md1");
  TMatrixDSym    *md2_r    = (TMatrixDSym*)    f->Get("md2");
  TMatrixDSparse *md3_r    = (TMatrixDSparse*) f->Get("md3");
  TVectorD       *vd_row_r = (TVectorD*)       f->Get("vd_row");

  ok &= ((*md1_r) == md1) ? kTRUE : kFALSE;
  ok &= ((*md2_r) == md2) ? kTRUE : kFALSE;
  ok &= ((*md3_r) == md3) ? kTRUE : kFALSE;
  ok &= ((*vd_row_r) == vd_row) ? kTRUE : kFALSE;

  StatusPrint(1,"Streamers",ok);
}

void cleanup()
{
  gSystem->Unlink("stress-vmatrix.root");
  gSystem->Unlink("stress-vvector.root");
  gSystem->Unlink("stress-vdecomp.root");
}
