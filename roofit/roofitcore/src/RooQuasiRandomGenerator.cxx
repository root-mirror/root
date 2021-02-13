/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooQuasiRandomGenerator.cxx
\class RooQuasiRandomGenerator
\ingroup Roofitcore

This class generates the quasi-random (aka "low discrepancy")
sequence for dimensions up to 12 using the Niederreiter base 2
algorithm described in Bratley, Fox, Niederreiter, ACM Trans.
Model. Comp. Sim. 2, 195 (1992). This implementation was adapted
from the 0.9 beta release of the GNU scientific library.
Quasi-random number sequences are useful for improving the
convergence of a Monte Carlo integration.
**/

#include "RooFit.h"

#include "RooQuasiRandomGenerator.h"
#include "RooMsgService.h"

#include <iostream>
#include <cassert>

using namespace std;

ClassImp(RooQuasiRandomGenerator);


////////////////////////////////////////////////////////////////////////////////
/// Perform one-time initialization of our static coefficient array if necessary
/// and initialize our workspace.

RooQuasiRandomGenerator::RooQuasiRandomGenerator() 
{
  if(!_coefsCalculated) {
    calculateCoefs(MaxDimension);
    _coefsCalculated= kTRUE;
  }
  // allocate workspace memory
  _nextq= new Int_t[MaxDimension];
  reset();
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooQuasiRandomGenerator::~RooQuasiRandomGenerator() 
{
  delete[] _nextq;
}


////////////////////////////////////////////////////////////////////////////////
/// Reset the workspace to its initial state.

void RooQuasiRandomGenerator::reset() 
{
  _sequenceCount= 0;
  for(Int_t dim= 0; dim < MaxDimension; dim++) _nextq[dim]= 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Generate the next number in the sequence for the specified dimension.
/// The maximum dimension supported is 12.

Bool_t RooQuasiRandomGenerator::generate(UInt_t dimension, Double_t vector[]) 
{
  /* Load the result from the saved state. */
  static const Double_t recip = 1.0/(double)(1U << NBits); /* 2^(-nbits) */

  UInt_t dim;
  for(dim=0; dim < dimension; dim++) {
    vector[dim] = _nextq[dim] * recip;
  }

  /* Find the position of the least-significant zero in sequence_count.
   * This is the bit that changes in the Gray-code representation as
   * the count is advanced.
   */
  Int_t r(0),c(_sequenceCount);
  while(1) {
    if((c % 2) == 1) {
      ++r;
      c /= 2;
    }
    else break;
  }
  if(r >= NBits) {
    oocoutE((TObject*)0,Integration) << "RooQuasiRandomGenerator::generate: internal error!" << endl;
    return kFALSE;
  }

  /* Calculate the next state. */
  for(dim=0; dim<dimension; dim++) {
    _nextq[dim] ^= _cj[r][dim];
  }
  _sequenceCount++;

  return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate the coefficients for the given number of dimensions

void RooQuasiRandomGenerator::calculateCoefs(UInt_t dimension) 
{
  int ci[NBits][NBits];
  int v[NBits+MaxDegree+1];
  int r;
  unsigned int i_dim;

  for(i_dim=0; i_dim<dimension; i_dim++) {

    const int poly_index = i_dim + 1;
    int j, k;

    /* Niederreiter (page 56, after equation (7), defines two
     * variables Q and U.  We do not need Q explicitly, but we
     * do need U.
     */
    int u = 0;

    /* For each dimension, we need to calculate powers of an
     * appropriate irreducible polynomial, see Niederreiter
     * page 65, just below equation (19).
     * Copy the appropriate irreducible polynomial into PX,
     * and its degree into E.  Set polynomial B = PX ** 0 = 1.
     * M is the degree of B.  Subsequently B will hold higher
     * powers of PX.
     */
    int pb[MaxDegree+1];
    int px[MaxDegree+1];
    int px_degree = _polyDegree[poly_index];
    int pb_degree = 0;

    for(k=0; k<=px_degree; k++) {
      px[k] = _primitivePoly[poly_index][k];
      pb[k] = 0;
    }
    pb[0] = 1;

    for(j=0; j<NBits; j++) {

      /* If U = 0, we need to set B to the next power of PX
       * and recalculate V.
       */
      if(u == 0) calculateV(px, px_degree, pb, &pb_degree, v, NBits+MaxDegree);

      /* Now C is obtained from V.  Niederreiter
       * obtains A from V (page 65, near the bottom), and then gets
       * C from A (page 56, equation (7)).  However this can be done
       * in one step.  Here CI(J,R) corresponds to
       * Niederreiter's C(I,J,R).
       */
      for(r=0; r<NBits; r++) {
        ci[r][j] = v[r+u];
      }

      /* Advance Niederreiter's state variables. */
      ++u;
      if(u == px_degree) u = 0;
    }

    /* The array CI now holds the values of C(I,J,R) for this value
     * of I.  We pack them into array CJ so that CJ(I,R) holds all
     * the values of C(I,J,R) for J from 1 to NBITS.
     */
    for(r=0; r<NBits; r++) {
      int term = 0;
      for(j=0; j<NBits; j++) {
        term = 2*term + ci[r][j];
      }
      _cj[r][i_dim] = term;
    }

  }
}


////////////////////////////////////////////////////////////////////////////////
/// Internal function

void RooQuasiRandomGenerator::calculateV(const int px[], int px_degree,
					 int pb[], int * pb_degree, int v[], int maxv) 
{
  const int nonzero_element = 1;    /* nonzero element of Z_2  */
  const int arbitrary_element = 1;  /* arbitray element of Z_2 */

  /* The polynomial ph is px**(J-1), which is the value of B on arrival.
   * In section 3.3, the values of Hi are defined with a minus sign:
   * don't forget this if you use them later !
   */
  int ph[MaxDegree+1];
  /* int ph_degree = *pb_degree; */
  int bigm = *pb_degree;      /* m from section 3.3 */
  int m;                      /* m from section 2.3 */
  int r, k, kj;

  for(k=0; k<=MaxDegree; k++) {
    ph[k] = pb[k];
  }

  /* Now multiply B by PX so B becomes PX**J.
   * In section 2.3, the values of Bi are defined with a minus sign :
   * don't forget this if you use them later !
   */
   polyMultiply(px, px_degree, pb, *pb_degree, pb, pb_degree);
   m = *pb_degree;

  /* Now choose a value of Kj as defined in section 3.3.
   * We must have 0 <= Kj < E*J = M.
   * The limit condition on Kj does not seem very relevant
   * in this program.
   */
  /* Quoting from BFN: "Our program currently sets each K_q
   * equal to eq. This has the effect of setting all unrestricted
   * values of v to 1."
   * Actually, it sets them to the arbitrary chosen value.
   * Whatever.
   */
  kj = bigm;

  /* Now choose values of V in accordance with
   * the conditions in section 3.3.
   */
  for(r=0; r<kj; r++) {
    v[r] = 0;
  }
  v[kj] = 1;


  if(kj >= bigm) {
    for(r=kj+1; r<m; r++) {
      v[r] = arbitrary_element;
    }
  }
  else {
    /* This block is never reached. */

    int term = sub(0, ph[kj]);

    for(r=kj+1; r<bigm; r++) {
      v[r] = arbitrary_element;

      /* Check the condition of section 3.3,
       * remembering that the H's have the opposite sign.  [????????]
       */
      term = sub(term, mul(ph[r], v[r]));
    }

    /* Now v[bigm] != term. */
    v[bigm] = add(nonzero_element, term);

    for(r=bigm+1; r<m; r++) {
      v[r] = arbitrary_element;
    }
  }

  /* Calculate the remaining V's using the recursion of section 2.3,
   * remembering that the B's have the opposite sign.
   */
  for(r=0; r<=maxv-m; r++) {
    int term = 0;
    for(k=0; k<m; k++) {
      term = sub(term, mul(pb[k], v[r+k]));
    }
    v[r+m] = term;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Internal function

void RooQuasiRandomGenerator::polyMultiply(const int pa[], int pa_degree, const int pb[],
					   int pb_degree, int pc[], int  * pc_degree) 
{
  int j, k;
  int pt[MaxDegree+1];
  int pt_degree = pa_degree + pb_degree;

  for(k=0; k<=pt_degree; k++) {
    int term = 0;
    for(j=0; j<=k; j++) {
      const int conv_term = mul(pa[k-j], pb[j]);
      term = add(term, conv_term);
    }
    pt[k] = term;
  }

  for(k=0; k<=pt_degree; k++) {
    pc[k] = pt[k];
  }
  for(k=pt_degree+1; k<=MaxDegree; k++) {
    pc[k] = 0;
  }

  *pc_degree = pt_degree;
}


////////////////////////////////////////////////////////////////////////////////

Int_t RooQuasiRandomGenerator::_cj[RooQuasiRandomGenerator::NBits]
[RooQuasiRandomGenerator::MaxDimension];

/* primitive polynomials in binary encoding */

////////////////////////////////////////////////////////////////////////////////

const Int_t RooQuasiRandomGenerator::_primitivePoly[RooQuasiRandomGenerator::MaxDimension+1]

////////////////////////////////////////////////////////////////////////////////

[RooQuasiRandomGenerator::MaxPrimitiveDegree+1] =
{
  { 1, 0, 0, 0, 0, 0 },  /*  1               */
  { 0, 1, 0, 0, 0, 0 },  /*  x               */
  { 1, 1, 0, 0, 0, 0 },  /*  1 + x           */
  { 1, 1, 1, 0, 0, 0 },  /*  1 + x + x^2     */
  { 1, 1, 0, 1, 0, 0 },  /*  1 + x + x^3     */
  { 1, 0, 1, 1, 0, 0 },  /*  1 + x^2 + x^3   */
  { 1, 1, 0, 0, 1, 0 },  /*  1 + x + x^4     */
  { 1, 0, 0, 1, 1, 0 },  /*  1 + x^3 + x^4   */
  { 1, 1, 1, 1, 1, 0 },  /*  1 + x + x^2 + x^3 + x^4   */
  { 1, 0, 1, 0, 0, 1 },  /*  1 + x^2 + x^5             */
  { 1, 0, 0, 1, 0, 1 },  /*  1 + x^3 + x^5             */
  { 1, 1, 1, 1, 0, 1 },  /*  1 + x + x^2 + x^3 + x^5   */
  { 1, 1, 1, 0, 1, 1 }   /*  1 + x + x^2 + x^4 + x^5   */
};

/* degrees of primitive polynomials */

////////////////////////////////////////////////////////////////////////////////

const Int_t RooQuasiRandomGenerator::_polyDegree[RooQuasiRandomGenerator::MaxDimension+1] =
{
  0, 1, 1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5
};

Bool_t RooQuasiRandomGenerator::_coefsCalculated= kFALSE;
