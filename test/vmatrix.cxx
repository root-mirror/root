// @(#)root/test:$Id$
// Author: Fons Rademakers and Eddy Offermann  Nov 2003

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Linear Algebra Package -- Matrix Verifications.                      //
//                                                                      //
// This file implements a large set of TMat operation tests.            //
// *******************************************************************  //
// *  Starting  Matrix - S T R E S S suite                              //
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
// *******************************************************************  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//_____________________________batch only_____________________
#ifndef __CINT__

#include <iostream>
#include "snprintf.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TArrayD.h"
#include "TMath.h"

#include "TDecompLU.h"
#include "TDecompQRH.h"
#include "TDecompSVD.h"

void stress_matrix                (Int_t verbose);
void StatusPrint                  (Int_t id,const TString &title,Int_t status);

void stress_allocation            ();
void stress_matrix_fill           (Int_t rsize,Int_t csize);
void stress_element_op            (Int_t rsize,Int_t csize);
void stress_binary_ebe_op         (Int_t rsize, Int_t csize);
void stress_transposition         (Int_t msize);
void stress_special_creation      (Int_t dim);
void stress_matrix_fill           (Int_t rsize,Int_t csize);
void stress_matrix_promises       (Int_t dim);
void stress_norms                 (Int_t rsize,Int_t csize);
void stress_determinant           (Int_t msize);
void stress_mm_multiplications    (Int_t msize);
void stress_sym_mm_multiplications(Int_t msize);
void stress_vm_multiplications    (Int_t msize);
void stress_inversion             (Int_t msize);
void stress_matrix_io             ();

int main(int argc,char **argv)
{
  Int_t verbose = 0;
  Char_t c;
  while (argc > 1 && argv[1][0] == '-' && argv[1][1] != 0) {
    for (Int_t i = 1; (c = argv[1][i]) != 0; i++) {
      switch (c) {
        case 'v':
          verbose = 1;
          break;
        default:
          Error("vmatrix", "unknown flag -%c",c);
          break;
      }
    }
    argc--;
    argv++;
  }
  stress_matrix(verbose);
  return 0;
}
#endif

#define EPSILON 1.0e-14

Int_t gVerbose = 0;

//________________________________common part_________________________

void stress_matrix(Int_t verbose=0)
{
  std::cout << "******************************************************************" <<std::endl;
  std::cout << "*  Starting  Matrix - S T R E S S suite                          *" <<std::endl;
  std::cout << "******************************************************************" <<std::endl;

  gVerbose = verbose;
  stress_allocation();
  stress_matrix_fill(20,10);
  stress_element_op(20,10);
  stress_binary_ebe_op(10,20);
  stress_transposition(20);
  stress_special_creation(20);
#ifndef __CINT__
  stress_matrix_promises(20);
#endif
  stress_norms(40,20);
  stress_determinant(20);
  stress_mm_multiplications(20);
  stress_sym_mm_multiplications(20);
  stress_vm_multiplications(20);
  stress_inversion(20);

  stress_matrix_io();
  std::cout << "******************************************************************" <<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Print test program number and its title
///   const Int_t kMAX = 65;
///   TString header = TString("Test ")+Form("%2d",id)+" : "+title;
///   const Int_t nch = header.Length();
///   for (Int_t i = nch; i < kMAX; i++) header += '.';
///   std::cout << header << (status ? "OK" : "FAILED") << std::endl;
/// Print test program number and its title

void StatusPrint(Int_t id,const Char_t *title,Bool_t status)
{
   const Int_t kMAX = 65;
   char header[80];
   snprintf(header,80,"Test %2d : %s",id,title);
   Int_t nch = strlen(header);
   for (Int_t i=nch;i<kMAX;i++) header[i] = '.';
   header[kMAX] = 0;
   header[kMAX-1] = ' ';
   std::cout << header << (status ? "OK" : "FAILED") << std::endl;
}

//------------------------------------------------------------------------
//          Test allocation functions and compatibility check
//
void stress_allocation()
{
  if (gVerbose)
    std::cout << "\n\n---> Test allocation and compatibility check" << std::endl;

  Bool_t ok = kTRUE;

  Int_t i,j;
  TMatrixD m1(4,20);
  for (i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++)
    for (j = m1.GetColLwb(); j <= m1.GetColUpb(); j++)
      m1(i,j) = TMath::Pi()*i+TMath::E()*j;

  TMatrixD m2(0,3,0,19);
  TMatrixD m3(1,4,0,19);
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
void stress_matrix_fill(Int_t rsize,Int_t csize)
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
      std::cout << "Check array Use" << std::endl;
    {
      TMatrixD *m1 = new TMatrixD(m);
      TMatrixD *m2 = new TMatrixD();
      m2->Use(m1->GetRowLwb(),m1->GetRowUpb(),m1->GetColLwb(),m1->GetColUpb(),m1->GetMatrixArray());
      ok &= VerifyMatrixIdentity(m,*m2,gVerbose,EPSILON);
      m2->Sqr();
      TMatrixD m3 = m; m3.Sqr();
      ok &= VerifyMatrixIdentity(m3,*m1,gVerbose,EPSILON);
      delete m1;
      delete m2;
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
typedef  double (*dfunc)(double);
class ApplyFunction : public TElementActionD {
   dfunc fFunc;
   void Operation(Double_t &element) const { element = fFunc(double(element)); }
public:
   ApplyFunction(dfunc func) : fFunc(func) { }
};

void stress_element_op(Int_t rsize,Int_t csize)
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
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);

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
void stress_binary_ebe_op(Int_t rsize, Int_t csize)
{
  if (gVerbose)
    std::cout << "\n---> Test Binary Matrix element-by-element operations" << std::endl;

  Bool_t ok = kTRUE;
  const double pattern = 4.25;

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
void stress_transposition(Int_t msize)
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

#ifndef __CINT__
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

void stress_special_creation(Int_t dim)
{
  if (gVerbose)
    std::cout << "\n---> Check creating some special matrices of dimension " << dim << std::endl;

  Int_t j;
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
      for (j = m1.GetColLwb(); j <= m1.GetColUpb(); j++)
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

void stress_matrix_promises(Int_t dim)
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

//
//------------------------------------------------------------------------
//             Verify the norm calculation
//
void stress_norms(Int_t rsize,Int_t csize)
{
  if (gVerbose)
    std::cout << "\n---> Verify norm calculations" << std::endl;

  Bool_t ok = kTRUE;
  const double pattern = 10.25;

  if (rsize % 2 == 1 || csize %2 == 1)
    Fatal("stress_norms","Sorry, size of the matrix to test must be even for this test\n");

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
void stress_determinant(Int_t msize)
{
  if (gVerbose)
    std::cout << "\n---> Verify determinant evaluation for a square matrix of size " << msize << std::endl;

  Bool_t ok = kTRUE;
  TMatrixD m(msize,msize);
  const double pattern = 2.5;

  if (gVerbose)
    std::cout << "\nCheck to see that the determinant of the unit matrix is one";
  m.UnitMatrix();
  if (gVerbose)
     std::cout << "\n\tdeterminant is " << m.Determinant();
  ok &= ( m.Determinant() == 1 ) ? kTRUE : kFALSE;

  if (gVerbose)
    std::cout << "\nCheck the determinant for the matrix with " << pattern << " at the diagonal";
  {
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        m(i,j) = ( i==j ? pattern : 0 );
  }
  if (gVerbose)
     std::cout << "\n\tdeterminant is " << m.Determinant() << " should be " << TMath::Power(pattern,(double)m.GetNrows()) <<std::endl;
  ok &= ( TMath::Abs(m.Determinant()-TMath::Power(pattern,(double)m.GetNrows())) < DBL_EPSILON  ) ? kTRUE : kFALSE;

  if (gVerbose)
    std::cout << "\nCheck the determinant of the transposed matrix";
  m.UnitMatrix();
  m(1,2) = pattern;
  TMatrixD m_tran(TMatrixD::kTransposed,m);
  ok &= ( !(m == m_tran) ) ? kTRUE : kFALSE;
  ok &= ( m.Determinant() == m_tran.Determinant() ) ? kTRUE : kFALSE;

  {
    if (gVerbose)
      std::cout << "\nswap two rows/cols of a matrix through method 1 and watch det's sign";
    m.UnitMatrix();
    TMatrixDRow(m,3).Assign(pattern);
    const double det1 = m.Determinant();
    TMatrixDRow row1(m,1);
    TVectorD vrow1(m.GetRowLwb(),m.GetRowUpb()); vrow1 = row1;
    TVectorD vrow3(m.GetRowLwb(),m.GetRowUpb()); vrow3 = TMatrixDRow(m,3);
    row1 = vrow3; TMatrixDRow(m,3) = vrow1;
    ok &= ( m.Determinant() == -det1 ) ? kTRUE : kFALSE;
    TMatrixDColumn col2(m,2);
    TVectorD vcol2(m.GetRowLwb(),m.GetRowUpb()); vcol2 = col2;
    TVectorD vcol4(m.GetRowLwb(),m.GetRowUpb()); vcol4 = TMatrixDColumn(m,4);
    col2 = vcol4; TMatrixDColumn(m,4) = vcol2;
    ok &= ( m.Determinant() == det1 ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      std::cout << "\nswap two rows/cols of a matrix through method 2 and watch det's sign";
    m.UnitMatrix();
    TMatrixDRow(m,3).Assign(pattern);
    const double det1 = m.Determinant();

    TMatrixD m_save( m);
    TMatrixDRow(m,1) = TMatrixDRow(m_save,3);
    TMatrixDRow(m,3) = TMatrixDRow(m_save,1);
    ok &= ( m.Determinant() == -det1 ) ? kTRUE : kFALSE;

    m_save = m;
    TMatrixDColumn(m,2) = TMatrixDColumn(m_save,4);
    TMatrixDColumn(m,4) = TMatrixDColumn(m_save,2);
    ok &= ( m.Determinant() == det1 ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    std::cout << "\nCheck the determinant for the matrix with " << pattern << " at the anti-diagonal";
  {
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        m(i,j) = ( i==(m.GetColUpb()+m.GetColLwb()-j) ? pattern : 0 );
    ok &= ( m.Determinant() == TMath::Power(pattern,(double)m.GetNrows()) *
               ( m.GetNrows()*(m.GetNrows()-1)/2 & 1 ? -1 : 1 ) ) ? kTRUE : kFALSE;
  }

  if (0)
  {
    if (gVerbose)
       std::cout << "\nCheck the determinant for the singular matrix"
                    "\n\tdefined as above with zero first row";
    m.Zero();
    {
      for (Int_t i = m.GetRowLwb()+1; i <= m.GetRowUpb(); i++)
        for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
          m(i,j) = ( i==(m.GetColUpb()+m.GetColLwb()-j) ? pattern : 0 );
    }
    if (gVerbose)
       std::cout << "\n\tdeterminant is " << m.Determinant();
    ok &= ( m.Determinant() == 0 ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    std::cout << "\nCheck out the determinant of the Hilbert matrix";
  TMatrixD H = THilbertMatrixD(3,3);
  if (gVerbose) {
    std::cout << "\n    3x3 Hilbert matrix: exact determinant 1/2160 ";
    std::cout << "\n                              computed    1/"<< 1/H.Determinant();
  }

  H.ResizeTo(4,4);
  H = THilbertMatrixD(4,4);
  if (gVerbose) {
    std::cout << "\n    4x4 Hilbert matrix: exact determinant 1/6048000 ";
    std::cout << "\n                              computed    1/"<< 1/H.Determinant();
  }

  H.ResizeTo(5,5);
  H = THilbertMatrixD(5,5);
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
    TDecompSVD svd(H);
    Double_t d1,d2;
    svd.Det(d1,d2);
    std::cout  << "\n svd det = " << d1*TMath::Power(2.0,d2) <<std::endl;
  }

  H.ResizeTo(7,7);
  H = THilbertMatrixD(7,7);
  if (gVerbose) {
    std::cout << "\n    7x7 Hilbert matrix: exact determinant 4.8358e-25";
    std::cout << "\n                              computed    "<< H.Determinant();
  }

  H.ResizeTo(9,9);
  H = THilbertMatrixD(9,9);
  if (gVerbose) {
    std::cout << "\n    9x9 Hilbert matrix: exact determinant 9.72023e-43";
    std::cout << "\n                              computed    "<< H.Determinant();
  }

  H.ResizeTo(10,10);
  H = THilbertMatrixD(10,10);
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
void stress_mm_multiplications(Int_t msize)
{
  if (gVerbose)
    std::cout << "\n---> Verify matrix multiplications "
            "for matrices of the characteristic size " << msize << std::endl;

  const Double_t epsilon = EPSILON*msize/100;

  Int_t i,j;
  Bool_t ok = kTRUE;
  {
    if (gVerbose)
      std::cout << "\nTest inline multiplications of the UnitMatrix" << std::endl;
    TMatrixD m = THilbertMatrixD(-1,msize,-1,msize);
    TMatrixD u(TMatrixD::kUnit,m);
    m(3,1) = TMath::Pi();
    u *= m;
    ok &= VerifyMatrixIdentity(u,m,gVerbose,epsilon);
  }

  {
    if (gVerbose)
      std::cout << "Test inline multiplications by a DiagMat" << std::endl;
    TMatrixD m = THilbertMatrixD(msize+3,msize);
    m(1,3) = TMath::Pi();
    TVectorD v(msize);
    for (i = v.GetLwb(); i <= v.GetUpb(); i++)
      v(i) = 1+i;
    TMatrixD diag(msize,msize);
    TMatrixDDiag d = TMatrixDDiag(diag);
    d = v;
    TMatrixD eth = m;
    for (i = eth.GetRowLwb(); i <= eth.GetRowUpb(); i++)
      for (j = eth.GetColLwb(); j <= eth.GetColUpb(); j++)
        eth(i,j) *= v(j);
    m *= diag;
    ok &= VerifyMatrixIdentity(m,eth,gVerbose,epsilon);
  }

  {
    if (gVerbose)
      std::cout << "Test XPP = X where P is a permutation matrix" << std::endl;
    TMatrixD m = THilbertMatrixD(msize-1,msize);
    m(2,3) = TMath::Pi();
    TMatrixD eth = m;
    TMatrixD p(msize,msize);
    for (i = p.GetRowLwb(); i <= p.GetRowUpb(); i++)
      p(p.GetRowUpb()+p.GetRowLwb()-i,i) = 1;
    m *= p;
    m *= p;
    ok &= VerifyMatrixIdentity(m,eth,gVerbose,epsilon);
  }

  {
    if (gVerbose)
      std::cout << "Test general matrix multiplication through inline mult" << std::endl;
    TMatrixD m = THilbertMatrixD(msize-2,msize);
    m(3,3) = TMath::Pi();
    TMatrixD mt(TMatrixD::kTransposed,m);
    TMatrixD p = THilbertMatrixD(msize,msize);
    TMatrixDDiag(p) += 1;
    TMatrixD mp(m,TMatrixD::kMult,p);
    TMatrixD m1 = m;
    m *= p;
    ok &= VerifyMatrixIdentity(m,mp,gVerbose,epsilon);
    TMatrixD mp1(mt,TMatrixD::kTransposeMult,p);
    VerifyMatrixIdentity(m,mp1,gVerbose,epsilon);
    ok &= ( !(m1 == m) );
    TMatrixD mp2(TMatrixD::kZero,m1);
    ok &= ( mp2 == 0 );
    mp2.Mult(m1,p);
    ok &= VerifyMatrixIdentity(m,mp2,gVerbose,epsilon);

    if (gVerbose)
      std::cout << "Test XP=X*P  vs XP=X;XP*=P" << std::endl;
    TMatrixD mp3 = m1*p;
    ok &= VerifyMatrixIdentity(m,mp3,gVerbose,epsilon);
  }

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
    TMatrixD haar = THaarMatrixD(5);
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
void stress_sym_mm_multiplications(Int_t msize)
{
  if (gVerbose)
    std::cout << "\n---> Verify symmetric matrix multiplications "
            "for matrices of the characteristic size " << msize << std::endl;

  Bool_t ok = kTRUE;

  Int_t i,j;
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
    TMatrixDSym m = THilbertMatrixDSym(msize);
    m(1,3) = TMath::Pi();
    m(3,1) = TMath::Pi();
    TVectorD v(msize);
    for (i = v.GetLwb(); i <= v.GetUpb(); i++)
      v(i) = 1+i;
    TMatrixDSym diag(msize);
    TMatrixDDiag d(diag); d = v;
    TMatrixDSym eth = m;
    for (i = eth.GetRowLwb(); i <= eth.GetRowUpb(); i++)
      for (j = eth.GetColLwb(); j <= eth.GetColUpb(); j++)
        eth(i,j) *= v(j);
    TMatrixD m2 = m * diag;
    ok &= VerifyMatrixIdentity(m2,eth,gVerbose,epsilon);
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "Test XPP = X where P is a permutation matrix" << std::endl;
    TMatrixDSym m = THilbertMatrixDSym(msize);
    m(2,3) = TMath::Pi();
    m(3,2) = TMath::Pi();
    TMatrixDSym eth = m;
    TMatrixDSym p(msize);
    for (i = p.GetRowLwb(); i <= p.GetRowUpb(); i++)
      p(p.GetRowUpb()+p.GetRowLwb()-i,i) = 1;
    TMatrixD m2 = m * p;
    m2 *= p;
    ok &= VerifyMatrixIdentity(m2,eth,gVerbose,epsilon);
  }

  if (ok)
  {
    if (gVerbose)
      std::cout << "Test general matrix multiplication through inline mult" << std::endl;
    TMatrixDSym m = THilbertMatrixDSym(msize);
    m(2,3) = TMath::Pi();
    m(3,2) = TMath::Pi();
    TMatrixDSym mt(TMatrixDSym::kTransposed,m);
    TMatrixDSym p = THilbertMatrixDSym(msize);
    TMatrixDDiag(p) += 1;
    TMatrixD mp(m,TMatrixD::kMult,p);
    TMatrixDSym m1 = m;
    TMatrixD m3(m,TMatrixD::kMult,p);
    memcpy(m.GetMatrixArray(),m3.GetMatrixArray(),msize*msize*sizeof(Double_t));
    ok &= VerifyMatrixIdentity(m,mp,gVerbose,epsilon);
    TMatrixD mp1(mt,TMatrixD::kTransposeMult,p);
    ok &= VerifyMatrixIdentity(m,mp1,gVerbose,epsilon);
    ok &= ( !(m1 == m) ) ? kTRUE : kFALSE;
    TMatrixDSym mp2(TMatrixDSym::kZero,m);
    ok &= ( mp2 == 0 ) ? kTRUE : kFALSE;

    if (gVerbose)
      std::cout << "Test XP=X*P  vs XP=X;XP*=P" << std::endl;
    TMatrixD mp3 = m1*p;
    ok &= VerifyMatrixIdentity(m,mp3,gVerbose,epsilon);
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
      TMatrixDSym  hth(TMatrixDSym::kAtA,haar);
      TMatrixD hht(haar,TMatrixD::kMult,haar_t);
      TMatrixD hht1 = haar; hht1 *= haar_t;
      ok &= VerifyMatrixIdentity(unit,hth,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(unit,hht,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(unit,hht1,gVerbose,epsilon);
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
void stress_vm_multiplications(Int_t msize)
{
  if (gVerbose)
    std::cout << "\n---> Verify vector-matrix multiplications "
           "for matrices of the characteristic size " << msize << std::endl;

  const Double_t epsilon = EPSILON*msize/100;

  Bool_t ok = kTRUE;
  {
    if (gVerbose)
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
    ok &= VerifyVectorIdentity(v1,v2,gVerbose,epsilon);
  }

  {
    if (gVerbose)
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
    ok &= VerifyVectorIdentity(v1,v2,gVerbose,epsilon);
  }

  {
    if (gVerbose)
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
    ok &= VerifyMatrixIdentity(mvb,mvm,gVerbose,epsilon);
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(12,"Matrix Vector Multiplications",ok);
}

//
//------------------------------------------------------------------------
//               Verify matrix inversion
//
void stress_inversion(Int_t msize)
{
  if (gVerbose)
    std::cout << "\n---> Verify matrix inversion for square matrices of size " << msize << std::endl;

  const Double_t epsilon = EPSILON*msize/10;

  Bool_t ok = kTRUE;
  {
    if (gVerbose)
      std::cout << "\nTest inversion of a diagonal matrix" << std::endl;
    TMatrixD m(-1,msize,-1,msize);
    TMatrixD mi(TMatrixD::kZero,m);
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      mi(i,i) = 1/(m(i,i)=i-m.GetRowLwb()+1);
    TMatrixD mi1(TMatrixD::kInverted,m);
    m.Invert();
    ok &= VerifyMatrixIdentity(m,mi,gVerbose,epsilon);
    ok &= VerifyMatrixIdentity(mi1,mi,gVerbose,epsilon);
  }

  {
    if (gVerbose)
      std::cout << "Test inversion of an orthonormal (Haar) matrix" << std::endl;
    TMatrixD m = THaarMatrixD(3);
    TMatrixD morig = m;
    TMatrixD mt(TMatrixD::kTransposed,m);
    double det = -1;         // init to a wrong val to see if it's changed
    m.Invert(&det);
    ok &= ( TMath::Abs(det-1) <= msize*epsilon ) ? kTRUE : kFALSE;
    ok &= VerifyMatrixIdentity(m,mt,gVerbose,epsilon);
    TMatrixD mti(TMatrixD::kInverted,mt);
    ok &= VerifyMatrixIdentity(mti,morig,gVerbose,msize*epsilon);
  }

  {
    if (gVerbose)
      std::cout << "Test inversion of a good matrix with diagonal dominance" << std::endl;
    TMatrixD m = THilbertMatrixD(msize,msize);
    TMatrixDDiag(m) += 1;
    TMatrixD morig = m;
    Double_t det_inv = 0;
    const Double_t det_comp = m.Determinant();
    m.Invert(&det_inv);
    if (gVerbose) {
      std::cout << "\tcomputed determinant             " << det_comp << std::endl;
      std::cout << "\tdeterminant returned by Invert() " << det_inv << std::endl;
    }

    if (gVerbose)
      std::cout << "\tcheck to see M^(-1) * M is E" << std::endl;
    TMatrixD mim(m,TMatrixD::kMult,morig);
    TMatrixD unit(TMatrixD::kUnit,m);
    ok &= VerifyMatrixIdentity(mim,unit,gVerbose,epsilon);

    if (gVerbose)
      std::cout << "\tcheck to see M * M^(-1) is E" << std::endl;
    TMatrixD mmi = morig; mmi *= m;
    ok &= VerifyMatrixIdentity(mmi,unit,gVerbose,epsilon);
  }

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(13,"Matrix Inversion",ok);
}

//
//------------------------------------------------------------------------
//           Test matrix I/O
//
void stress_matrix_io()
{
  if (gVerbose)
    std::cout << "\n---> Test matrix I/O" << std::endl;

  Bool_t ok = kTRUE;
  const double pattern = TMath::Pi();

  TMatrixD m(40,40);
  m = pattern;

  if (gVerbose)
    std::cout << "\nWrite matrix m to database" << std::endl;

  TFile *f = new TFile("vmatrix.root", "RECREATE");

  m.Write("m");

  if (gVerbose)
    std::cout << "\nClose database" << std::endl;
  delete f;

  if (gVerbose)
    std::cout << "\nOpen database in read-only mode and read matrix" << std::endl;
  TFile *f1 = new TFile("vmatrix.root");

  TMatrixD *mr = (TMatrixD*) f1->Get("m");

  if (gVerbose)
    std::cout << "\nRead matrix should be same as original still in memory" << std::endl;
  ok &= ((*mr) == m) ? kTRUE : kFALSE;

  delete f1;

  if (gVerbose)
    std::cout << "\nDone\n" << std::endl;

  StatusPrint(14,"Matrix Persistence",ok);
}
