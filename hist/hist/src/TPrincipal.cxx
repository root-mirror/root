// @(#)root/hist:$Id$
// Author: Christian Holm Christensen    1/8/2000

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TPrincipal
    \ingroup Hist
Principal Components Analysis (PCA)

The current implementation is based on the LINTRA package from CERNLIB
by R. Brun, H. Hansroul, and J. Kubler.
The class has been implemented by Christian Holm Christensen in August 2000.

## Introduction

In many applications of various fields of research, the treatment of
large amounts of data requires powerful techniques capable of rapid
data reduction and analysis. Usually, the quantities most
conveniently measured by the experimentalist, are not necessarily the
most significant for classification and analysis of the data. It is
then useful to have a way of selecting an optimal set of variables
necessary for the recognition process and reducing the dimensionality
of the problem, resulting in an easier classification procedure.

This paper describes the implementation of one such method of
feature selection, namely the principal components analysis. This
multidimensional technique is well known in the field of pattern
recognition and and its use in Particle Physics has been documented
elsewhere (cf. H. Wind, <I>Function Parameterization</I>, CERN
72-21).

## Overview
Suppose we have prototypes which are trajectories of particles,
passing through a spectrometer. If one measures the passage of the
particle at say 8 fixed planes, the trajectory is described by an
8-component vector:
\f[
  \mathbf{x} = \left(x_0, x_1, \ldots, x_7\right)
\f]
in 8-dimensional pattern space.

One proceeds by generating a a representative tracks sample and
building up the covariance matrix \f$\mathsf{C}\f$. Its eigenvectors and
eigenvalues are computed by standard methods, and thus a new basis is
obtained for the original 8-dimensional space the expansion of the
prototypes,
\f[
  \mathbf{x}_m = \sum^7_{i=0} a_{m_i} \mathbf{e}_i
  \quad
  \mbox{where}
  \quad
  a_{m_i} = \mathbf{x}^T\bullet\mathbf{e}_i
\f]
allows the study of the behavior of the coefficients \f$a_{m_i}\f$ for all
the tracks of the sample. The eigenvectors which are insignificant for
the trajectory description in the expansion will have their
corresponding coefficients \f$a_{m_i}\f$ close to zero for all the
prototypes.

On one hand, a reduction of the dimensionality is then obtained by
omitting these least significant vectors in the subsequent analysis.

On the other hand, in the analysis of real data, these least
significant variables(?) can be used for the pattern
recognition problem of extracting the valid combinations of
coordinates describing a true trajectory from the set of all possible
wrong combinations.

The program described here performs this principal components analysis
on a sample of data provided by the user. It computes the covariance
matrix, its eigenvalues ands corresponding eigenvectors and exhibits
the behavior of the principal components \f$a_{m_i}\f$, thus providing
to the user all the means of understanding their data.

## Principal Components Method
Let's consider a sample of \f$M\f$ prototypes each being characterized by
\f$P\f$ variables \f$x_0, x_1, \ldots, x_{P-1}\f$. Each prototype is a point, or a
column vector, in a \f$P\f$-dimensional *Pattern space*.
\f[
  \mathbf{x} = \left[\begin{array}{c}
    x_0\\x_1\\\vdots\\x_{P-1}\end{array}\right]\,,
\f]
where each \f$x_n\f$ represents the particular value associated with the
\f$n\f$-dimension.

Those \f$P\f$ variables are the quantities accessible to the
experimentalist, but are not necessarily the most significant for the
classification purpose.

The *Principal Components Method* consists of applying a
*linear* transformation to the original variables. This
transformation is described by an orthogonal matrix and is equivalent
to a rotation of the original pattern space into a new set of
coordinate vectors, which hopefully provide easier feature
identification and dimensionality reduction.

Let's define the covariance matrix:
\f[
  \mathsf{C} = \left\langle\mathbf{y}\mathbf{y}^T\right\rangle
  \quad\mbox{where}\quad
  \mathbf{y} = \mathbf{x} - \left\langle\mathbf{x}\right\rangle\,,
\f]
and the brackets indicate mean value over the sample of \f$M\f$
prototypes.

This matrix \f$\mathsf{C}\f$ is real, positive definite, symmetric, and will
have all its eigenvalues greater then zero. It will now be show that
among the family of all the complete orthonormal bases of the pattern
space, the base formed by the eigenvectors of the covariance matrix
and belonging to the largest eigenvalues, corresponds to the most
significant features of the description of the original prototypes.

let the prototypes be expanded on into a set of \f$N\f$ basis vectors
\f$\mathbf{e}_n, n=0,\ldots,N,N+1, \ldots, P-1\f$
\f[
  \mathbf{y}_i = \sum^N_{i=0} a_{i_n} \mathbf{e}_n,
  \quad
  i = 1, \ldots, M,
  \quad
  N < P-1
\f]
The `best' feature coordinates \f$\mathbf{e}_n\f$, spanning a *feature
space*,  will be obtained by minimizing the error due to this
truncated expansion, i.e.,
\f[
  \min\left(E_N\right) =
  \min\left[\left\langle\left(\mathbf{y}_i - \sum^N_{i=0} a_{i_n} \mathbf{e}_n\right)^2\right\rangle\right]
\f]
with the conditions:
\f[
  \mathbf{e}_k\bullet\mathbf{e}_j = \delta_{jk} =
  \left\{\begin{array}{rcl}
    1 & \mbox{for} & k = j\\
    0 & \mbox{for} & k \neq j
  \end{array}\right.
\f]
Multiplying (3) by \f$\mathbf{e}^T_n\f$ using (5),
we get
\f[
  a_{i_n} = \mathbf{y}_i^T\bullet\mathbf{e}_n\,,
\f]
so the error becomes
\f{eqnarray*}{
  E_N &=&
  \left\langle\left[\sum_{n=N+1}^{P-1}  a_{i_n}\mathbf{e}_n\right]^2\right\rangle\nonumber\\
  &=&
  \left\langle\left[\sum_{n=N+1}^{P-1}  \mathbf{y}_i^T\bullet\mathbf{e}_n\mathbf{e}_n\right]^2\right\rangle\nonumber\\
  &=&
  \left\langle\sum_{n=N+1}^{P-1}  \mathbf{e}_n^T\mathbf{y}_i\mathbf{y}_i^T\mathbf{e}_n\right\rangle\nonumber\\
  &=&
  \sum_{n=N+1}^{P-1}  \mathbf{e}_n^T\mathsf{C}\mathbf{e}_n
\f}
The minimization of the sum in (7) is obtained when each
term \f$\mathbf{e}_n^\mathsf{C}\mathbf{e}_n\f$ is minimum, since \f$\mathsf{C}\f$ is
positive definite. By the method of Lagrange multipliers, and the
condition (5), we get
\f[
  E_N = \sum^{P-1}_{n=N+1} \left(\mathbf{e}_n^T\mathsf{C}\mathbf{e}_n -
    l_n\mathbf{e}_n^T\bullet\mathbf{e}_n + l_n\right)
\f]
The minimum condition \f$\frac{dE_N}{d\mathbf{e}^T_n} = 0\f$ leads to the
equation
\f[
  \mathsf{C}\mathbf{e}_n = l_n\mathbf{e}_n\,,
\f]
which shows that \f$\mathbf{e}_n\f$ is an eigenvector of the covariance
matrix \f$\mathsf{C}\f$ with eigenvalue \f$l_n\f$. The estimated minimum error is
then given by
\f[
  E_N \sim \sum^{P-1}_{n=N+1} \mathbf{e}_n^T\bullet l_n\mathbf{e}_n
      = \sum^{P-1}_{n=N+1}  l_n\,,
\f]
where \f$l_n,\,n=N+1,\ldots,P\f$ \f$l_n,\,n=N+1,\ldots,P-1\f$ are the eigenvalues associated with the
omitted eigenvectors in the expansion (3). Thus, by choosing
the \f$N\f$ largest eigenvalues, and their associated eigenvectors, the
error \f$E_N\f$ is minimized.

The transformation matrix to go from the pattern space to the feature
space consists of the ordered eigenvectors \f$\mathbf{e}_1,\ldots,\mathbf{e}_P\f$
\f$\mathbf{e}_0,\ldots,\mathbf{e}_{P-1}\f$ for its columns
\f[
  \mathsf{T} = \left[
    \begin{array}{cccc}
      \mathbf{e}_0 &
      \mathbf{e}_1 &
      \vdots &
      \mathbf{e}_{P-1}
    \end{array}\right]
  = \left[
    \begin{array}{cccc}
      \mathbf{e}_{0_0} &  \mathbf{e}_{1_0} & \cdots &  \mathbf{e}_{{P-1}_0}\\
      \mathbf{e}_{0_1} &  \mathbf{e}_{1_1} & \cdots &  \mathbf{e}_{{P-1}_1}\\
      \vdots        &  \vdots        & \ddots &  \vdots \\
      \mathbf{e}_{0_{P-1}} &  \mathbf{e}_{1_{P-1}} & \cdots &  \mathbf{e}_{{P-1}_{P-1}}\\
    \end{array}\right]
\f]
This is an orthogonal transformation, or rotation, of the pattern
space and feature selection results in ignoring certain coordinates
in the transformed space.

Christian Holm August 2000, CERN
*/

#include "TPrincipal.h"

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"
#include "TMath.h"
#include "TList.h"
#include "TH2.h"
#include "TDatime.h"
#include "TBrowser.h"
#include "TROOT.h"
#include "Riostream.h"


ClassImp(TPrincipal);

////////////////////////////////////////////////////////////////////////////////
/// Empty constructor. Do not use.

TPrincipal::TPrincipal()
  : fMeanValues(0),
    fSigmas(0),
    fCovarianceMatrix(1,1),
    fEigenVectors(1,1),
    fEigenValues(0),
    fOffDiagonal(0),
    fStoreData(kFALSE)
{
   fTrace              = 0;
   fHistograms         = 0;
   fIsNormalised       = kFALSE;
   fNumberOfDataPoints = 0;
   fNumberOfVariables  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor. Argument is number of variables in the sample of data
/// Options are:
///  - N       Normalize the covariance matrix (default)
///  - D       Store input data (default)
///
/// The created object is  named "principal" by default.

TPrincipal::TPrincipal(Int_t nVariables, Option_t *opt)
  : fMeanValues(nVariables),
    fSigmas(nVariables),
    fCovarianceMatrix(nVariables,nVariables),
    fEigenVectors(nVariables,nVariables),
    fEigenValues(nVariables),
    fOffDiagonal(nVariables),
    fStoreData(kFALSE)
{
   if (nVariables <= 1) {
      Error("TPrincipal", "You can't be serious - nVariables == 1!!!");
      return;
   }

   SetName("principal");

   fTrace              = 0;
   fHistograms         = 0;
   fIsNormalised       = kFALSE;
   fNumberOfDataPoints = 0;
   fNumberOfVariables  = nVariables;
   while (strlen(opt) > 0) {
      switch(*opt++) {
         case 'N':
         case 'n':
            fIsNormalised = kTRUE;
            break;
         case 'D':
         case 'd':
            fStoreData    = kTRUE;
            break;
         default:
            break;
      }
   }

   if (!fMeanValues.IsValid())
      Error("TPrincipal","Couldn't create vector mean values");
   if (!fSigmas.IsValid())
      Error("TPrincipal","Couldn't create vector sigmas");
   if (!fCovarianceMatrix.IsValid())
      Error("TPrincipal","Couldn't create covariance matrix");
   if (!fEigenVectors.IsValid())
      Error("TPrincipal","Couldn't create eigenvector matrix");
   if (!fEigenValues.IsValid())
      Error("TPrincipal","Couldn't create eigenvalue vector");
   if (!fOffDiagonal.IsValid())
      Error("TPrincipal","Couldn't create offdiagonal vector");
   if (fStoreData) {
      fUserData.ResizeTo(nVariables*1000);
      fUserData.Zero();
      if (!fUserData.IsValid())
         Error("TPrincipal","Couldn't create user data vector");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TPrincipal::TPrincipal(const TPrincipal& pr) :
  TNamed(pr),
  fNumberOfDataPoints(pr.fNumberOfDataPoints),
  fNumberOfVariables(pr.fNumberOfVariables),
  fMeanValues(pr.fMeanValues),
  fSigmas(pr.fSigmas),
  fCovarianceMatrix(pr.fCovarianceMatrix),
  fEigenVectors(pr.fEigenVectors),
  fEigenValues(pr.fEigenValues),
  fOffDiagonal(pr.fOffDiagonal),
  fUserData(pr.fUserData),
  fTrace(pr.fTrace),
  fHistograms(pr.fHistograms),
  fIsNormalised(pr.fIsNormalised),
  fStoreData(pr.fStoreData)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TPrincipal& TPrincipal::operator=(const TPrincipal& pr)
{
   if(this!=&pr) {
      TNamed::operator=(pr);
      fNumberOfDataPoints=pr.fNumberOfDataPoints;
      fNumberOfVariables=pr.fNumberOfVariables;
      fMeanValues=pr.fMeanValues;
      fSigmas=pr.fSigmas;
      fCovarianceMatrix=pr.fCovarianceMatrix;
      fEigenVectors=pr.fEigenVectors;
      fEigenValues=pr.fEigenValues;
      fOffDiagonal=pr.fOffDiagonal;
      fUserData=pr.fUserData;
      fTrace=pr.fTrace;
      fHistograms=pr.fHistograms;
      fIsNormalised=pr.fIsNormalised;
      fStoreData=pr.fStoreData;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TPrincipal::~TPrincipal()
{
   if (fHistograms) {
      fHistograms->Delete();
      delete fHistograms;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add a data point and update the covariance matrix. The input
/// array must be <TT>fNumberOfVariables</TT> long.
///
///
/// The Covariance matrix and mean values of the input data is calculated
/// on the fly by the following equations:
///
/// \f[
/// \left<x_i\right>^{(0)}  = x_{i0}
/// \f]
///
///
/// \f[
/// \left<x_i\right>^{(n)} = \left<x_i\right>^{(n-1)}
/// + \frac1n \left(x_{in} - \left<x_i\right>^{(n-1)}\right)
/// \f]
///
/// \f[
/// C_{ij}^{(0)} = 0
/// \f]
///
///
///
/// \f[
/// C_{ij}^{(n)} = C_{ij}^{(n-1)}
/// + \frac1{n-1}\left[\left(x_{in} - \left<x_i\right>^{(n)}\right)
///   \left(x_{jn} - \left<x_j\right>^{(n)}\right)\right]
/// - \frac1n C_{ij}^{(n-1)}
/// \f]
///
/// since this is a really fast method, with no rounding errors (please
/// refer to CERN 72-21 pp. 54-106).
///
///
/// The data is stored internally in a <TT>TVectorD</TT>, in the following
/// way:
///
/// \f[
/// \mathbf{x} = \left[\left(x_{0_0},\ldots,x_{{P-1}_0}\right),\ldots,
///     \left(x_{0_i},\ldots,x_{{P-1}_i}\right), \ldots\right]
/// \f]
///
/// With \f$P\f$ as defined in the class description.

void TPrincipal::AddRow(const Double_t *p)
{
   if (!p)
      return;

   // Increment the data point counter
   Int_t i,j;
   if (++fNumberOfDataPoints == 1) {
      for (i = 0; i < fNumberOfVariables; i++)
         fMeanValues(i) = p[i];
   }
   else {

      const Double_t invnp = 1. / Double_t(fNumberOfDataPoints);
      const Double_t invnpM1 = 1. /(Double_t(fNumberOfDataPoints - 1));
      const Double_t cor = 1. - invnp;
      // use directly vector array for faster element access
      Double_t * meanValues = fMeanValues.GetMatrixArray();
      Double_t * covMatrix =  fCovarianceMatrix.GetMatrixArray();
      for (i = 0; i < fNumberOfVariables; i++) {

         meanValues[i] *= cor;
         meanValues[i] += p[i] * invnp;
         const Double_t t1 = (p[i] - meanValues[i]) * invnpM1;

         // Setting Matrix (lower triangle) elements
         for (j = 0; j < i + 1; j++) {
            const Int_t index = i * fNumberOfVariables + j;
            covMatrix[index] *= cor;
            covMatrix[index] += t1 * (p[j] - meanValues[j]);
         }
      }
   }

   // Store data point in internal vector
   // If the vector isn't big enough to hold the new data, then
   // expand the vector by half it's size.
   if (!fStoreData)
      return;
   Int_t size = fUserData.GetNrows();
   if (fNumberOfDataPoints * fNumberOfVariables > size)
      fUserData.ResizeTo(size + size/2);

   for (i = 0; i < fNumberOfVariables; i++) {
      j = (fNumberOfDataPoints-1) * fNumberOfVariables + i;
      fUserData.GetMatrixArray()[j] = p[i];
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Browse the TPrincipal object in the TBrowser.

void TPrincipal::Browse(TBrowser *b)
{
   if (fHistograms) {
      TIter next(fHistograms);
      TH1* h = 0;
      while ((h = (TH1*)next()))
         b->Add(h,h->GetName());
   }

   if (fStoreData)
      b->Add(&fUserData,"User Data");
   b->Add(&fCovarianceMatrix,"Covariance Matrix");
   b->Add(&fMeanValues,"Mean value vector");
   b->Add(&fSigmas,"Sigma value vector");
   b->Add(&fEigenValues,"Eigenvalue vector");
   b->Add(&fEigenVectors,"Eigenvector Matrix");

}

////////////////////////////////////////////////////////////////////////////////
/// Clear the data in Object. Notice, that's not possible to change
/// the dimension of the original data.

void TPrincipal::Clear(Option_t *opt)
{
   if (fHistograms) {
      fHistograms->Delete(opt);
   }

   fNumberOfDataPoints = 0;
   fTrace              = 0;
   fCovarianceMatrix.Zero();
   fEigenVectors.Zero();
   fEigenValues.Zero();
   fMeanValues.Zero();
   fSigmas.Zero();
   fOffDiagonal.Zero();

   if (fStoreData) {
      fUserData.ResizeTo(fNumberOfVariables * 1000);
      fUserData.Zero();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return a row of the user supplied data.
/// If row is out of bounds, 0 is returned.
/// It's up to the user to delete the returned array.
/// Row 0 is the first row;

const Double_t *TPrincipal::GetRow(Int_t row)
{
   if (row >= fNumberOfDataPoints)
      return 0;

   if (!fStoreData)
      return 0;

   Int_t index   = row  * fNumberOfVariables;
   return &fUserData(index);
}


////////////////////////////////////////////////////////////////////////////////
/// Generates the file `<filename>`, with `.C` appended if it does
/// argument doesn't end in .cxx or .C.
///
/// The file contains the implementation of two functions
/// ~~~ {.cpp}
///    void X2P(Double_t *x, Double *p)
///    void P2X(Double_t *p, Double *x, Int_t nTest)
/// ~~~
/// which does the same as  `TPrincipal::X2P` and `TPrincipal::P2X`
/// respectively. Please refer to these methods.
///
/// Further, the static variables:
/// ~~~ {.cpp}
///    Int_t    gNVariables
///    Double_t gEigenValues[]
///    Double_t gEigenVectors[]
///    Double_t gMeanValues[]
///    Double_t gSigmaValues[]
/// ~~~
/// are initialized. The only ROOT header file needed is Rtypes.h
///
/// See TPrincipal::MakeRealCode for a list of options

void TPrincipal::MakeCode(const char *filename, Option_t *opt)
{
   TString outName(filename);
   if (!outName.EndsWith(".C") && !outName.EndsWith(".cxx"))
      outName += ".C";

   MakeRealCode(outName.Data(),"",opt);
}

////////////////////////////////////////////////////////////////////////////////
/// Make histograms of the result of the analysis.
/// The option string say which histograms to create
///     - X         Histogram original data
///     - P         Histogram principal components corresponding to
///                 original data
///     - D         Histogram the difference between the original data
///                 and the projection of principal unto a lower
///                 dimensional subspace (2D histograms)
///     - E         Histogram the eigenvalues
///     - S         Histogram the square of the residues
///                 (see `TPrincipal::SumOfSquareResiduals`)
/// The histograms will be named `<name>_<type><number>`, where `<name>`
/// is the first argument, `<type>` is one of X,P,D,E,S, and `<number>`
/// is the variable.

void TPrincipal::MakeHistograms(const char *name, Option_t *opt)
{
   Bool_t makeX  = kFALSE;
   Bool_t makeD  = kFALSE;
   Bool_t makeP  = kFALSE;
   Bool_t makeE  = kFALSE;
   Bool_t makeS  = kFALSE;

   Int_t len     = strlen(opt);
   Int_t i,j,k;
   for (i = 0; i < len; i++) {
      switch (opt[i]) {
         case 'X':
         case 'x':
            if (fStoreData)
               makeX = kTRUE;
            break;
         case 'd':
         case 'D':
            if (fStoreData)
               makeD = kTRUE;
            break;
         case 'P':
         case 'p':
            if (fStoreData)
               makeP = kTRUE;
            break;
         case 'E':
         case 'e':
            makeE = kTRUE;
            break;
         case 's':
         case 'S':
            if (fStoreData)
               makeS = kTRUE;
            break;
         default:
            Warning("MakeHistograms","Unknown option: %c",opt[i]);
      }
   }

   // If no option was given, then exit gracefully
   if (!makeX && !makeD && !makeP && !makeE && !makeS)
      return;

   // If the list of histograms doesn't exist, create it.
   if (!fHistograms)
      fHistograms = new TList;

   // Don't create the histograms if they are already in the TList.
   if (makeX && fHistograms->FindObject(Form("%s_x000",name)))
      makeX = kFALSE;
   if (makeD && fHistograms->FindObject(Form("%s_d000",name)))
      makeD = kFALSE;
   if (makeP && fHistograms->FindObject(Form("%s_p000",name)))
      makeP = kFALSE;
   if (makeE && fHistograms->FindObject(Form("%s_e",name)))
      makeE = kFALSE;
   if (makeS && fHistograms->FindObject(Form("%s_s",name)))
      makeS = kFALSE;

   TH1F **hX  = 0;
   TH2F **hD  = 0;
   TH1F **hP  = 0;
   TH1F *hE   = 0;
   TH1F *hS   = 0;

   // Initialize the arrays of histograms needed
   if (makeX)
      hX = new TH1F * [fNumberOfVariables];

   if (makeD)
      hD = new TH2F * [fNumberOfVariables];

   if (makeP)
      hP = new TH1F * [fNumberOfVariables];

   if (makeE){
      hE = new TH1F(Form("%s_e",name), "Eigenvalues of Covariance matrix",
         fNumberOfVariables,0,fNumberOfVariables);
      hE->SetXTitle("Eigenvalue");
      fHistograms->Add(hE);
   }

   if (makeS) {
      hS = new TH1F(Form("%s_s",name),"E_{N}",
         fNumberOfVariables-1,1,fNumberOfVariables);
      hS->SetXTitle("N");
      hS->SetYTitle("#sum_{i=1}^{M} (x_{i} - x'_{N,i})^{2}");
      fHistograms->Add(hS);
   }

   // Initialize sub elements of the histogram arrays
   for (i = 0; i < fNumberOfVariables; i++) {
      if (makeX) {
         // We allow 4 sigma spread in the original data in our
         // histogram.
         Double_t xlowb  = fMeanValues(i) - 4 * fSigmas(i);
         Double_t xhighb = fMeanValues(i) + 4 * fSigmas(i);
         Int_t    xbins  = fNumberOfDataPoints/100;
         hX[i]           = new TH1F(Form("%s_x%03d", name, i),
            Form("Pattern space, variable %d", i),
            xbins,xlowb,xhighb);
         hX[i]->SetXTitle(Form("x_{%d}",i));
         fHistograms->Add(hX[i]);
      }

      if(makeD) {
         // The upper limit below is arbitrary!!!
         Double_t dlowb  = 0;
         Double_t dhighb = 20;
         Int_t    dbins  = fNumberOfDataPoints/100;
         hD[i]           = new TH2F(Form("%s_d%03d", name, i),
            Form("Distance from pattern to "
            "feature space, variable %d", i),
            dbins,dlowb,dhighb,
            fNumberOfVariables-1,
            1,
            fNumberOfVariables);
         hD[i]->SetXTitle(Form("|x_{%d} - x'_{%d,N}|/#sigma_{%d}",i,i,i));
         hD[i]->SetYTitle("N");
         fHistograms->Add(hD[i]);
      }

      if(makeP) {
         // For some reason, the trace of the none-scaled matrix
         // (see TPrincipal::MakeNormalised) should enter here. Taken
         // from LINTRA code.
         Double_t et = TMath::Abs(fEigenValues(i) * fTrace);
         Double_t plowb   = -10 * TMath::Sqrt(et);
         Double_t phighb  = -plowb;
         Int_t    pbins   = 100;
         hP[i]            = new TH1F(Form("%s_p%03d", name, i),
            Form("Feature space, variable %d", i),
            pbins,plowb,phighb);
         hP[i]->SetXTitle(Form("p_{%d}",i));
         fHistograms->Add(hP[i]);
      }

      if (makeE)
         // The Eigenvector histogram is easy
         hE->Fill(i,fEigenValues(i));

   }
   if (!makeX && !makeP && !makeD && !makeS) {
      if (hX)
         delete[] hX;
      if (hD)
         delete[] hD;
      if (hP)
         delete[] hP;
      return;
   }

   Double_t *x = 0;
   Double_t *p = new Double_t[fNumberOfVariables];
   Double_t *d = new Double_t[fNumberOfVariables];
   for (i = 0; i < fNumberOfDataPoints; i++) {

      // Zero arrays
      for (j = 0; j < fNumberOfVariables; j++)
         p[j] = d[j] = 0;

      // update the original data histogram
      x  = (Double_t*)(GetRow(i));
      R__ASSERT(x);

      if (makeP||makeD||makeS)
         // calculate the corresponding principal component
         X2P(x,p);

      if (makeD || makeS) {
         // Calculate the difference between the original data, and the
         // same project onto principal components, and then to a lower
         // dimensional sub-space
         for (j = fNumberOfVariables; j > 0; j--) {
            P2X(p,d,j);

            for (k = 0; k < fNumberOfVariables; k++) {
               // We use the absolute value of the difference!
               d[k] = x[k] - d[k];

               if (makeS)
                  hS->Fill(j,d[k]*d[k]);

               if (makeD) {
                  d[k] = TMath::Abs(d[k]) / (fIsNormalised ? fSigmas(k) : 1);
                  (hD[k])->Fill(d[k],j);
               }
            }
         }
      }

      if (makeX||makeP) {
         // If we are asked to make any of these histograms, we have to
         // go here
         for (j = 0; j < fNumberOfVariables; j++) {
            if (makeX)
               (hX[j])->Fill(x[j]);

            if (makeP)
               (hP[j])->Fill(p[j]);
         }
      }
   }
   // Clean up
   if (hX)
      delete [] hX;
   if (hD)
      delete [] hD;
   if (hP)
      delete [] hP;
   if (d)
      delete [] d;
   if (p)
      delete [] p;

   // Normalize the residues
   if (makeS)
      hS->Scale(Double_t(1.)/fNumberOfDataPoints);
}

////////////////////////////////////////////////////////////////////////////////
/// Normalize the covariance matrix

void TPrincipal::MakeNormalised()
{
   Int_t i,j;
   for (i = 0; i < fNumberOfVariables; i++) {
      fSigmas(i) = TMath::Sqrt(fCovarianceMatrix(i,i));
      if (fIsNormalised)
         for (j = 0; j <= i; j++)
            fCovarianceMatrix(i,j) /= (fSigmas(i) * fSigmas(j));

      fTrace += fCovarianceMatrix(i,i);
   }

   // Fill remaining parts of matrix, and scale.
   for (i = 0; i < fNumberOfVariables; i++)
      for (j = 0; j <= i; j++) {
         fCovarianceMatrix(i,j) /= fTrace;
         fCovarianceMatrix(j,i) = fCovarianceMatrix(i,j);
      }

}

////////////////////////////////////////////////////////////////////////////////
/// Generate the file <classname>PCA.cxx which contains the
/// implementation of two methods:
/// ~~~ {.cpp}
///    void <classname>::X2P(Double_t *x, Double *p)
///    void <classname>::P2X(Double_t *p, Double *x, Int_t nTest)
/// ~~~
/// which does the same as  TPrincipal::X2P and TPrincipal::P2X
/// respectively. Please refer to these methods.
///
/// Further, the public static members:
/// ~~~ {.cpp}
///    Int_t    <classname>::fgNVariables
///    Double_t <classname>::fgEigenValues[]
///    Double_t <classname>::fgEigenVectors[]
///    Double_t <classname>::fgMeanValues[]
///    Double_t <classname>::fgSigmaValues[]
/// ~~~
/// are initialized, and assumed to exist. The class declaration is
/// assumed to be in <classname>.h and assumed to be provided by the
/// user.
///
/// See TPrincipal::MakeRealCode for a list of options
///
/// The minimal class definition is:
/// ~~~ {.cpp}
///   class <classname> {
///   public:
///     static Int_t    fgNVariables;
///     static Double_t fgEigenVectors[];
///     static Double_t fgEigenValues[];
///     static Double_t fgMeanValues[];
///     static Double_t fgSigmaValues[];
///
///     void X2P(Double_t *x, Double_t *p);
///     void P2X(Double_t *p, Double_t *x, Int_t nTest);
///   };
/// ~~~
/// Whether the methods <classname>::X2P and <classname>::P2X should
/// be static or not, is up to the user.

void TPrincipal::MakeMethods(const char *classname, Option_t *opt)
{

   MakeRealCode(Form("%sPCA.cxx", classname), classname, opt);
}


////////////////////////////////////////////////////////////////////////////////
/// Perform the principal components analysis.
/// This is done in several stages in the TMatrix::EigenVectors method:
///  - Transform the covariance matrix into a tridiagonal matrix.
///  - Find the eigenvalues and vectors of the tridiagonal matrix.

void TPrincipal::MakePrincipals()
{
   // Normalize covariance matrix
   MakeNormalised();

   TMatrixDSym sym; sym.Use(fCovarianceMatrix.GetNrows(),fCovarianceMatrix.GetMatrixArray());
   TMatrixDSymEigen eigen(sym);
   fEigenVectors = eigen.GetEigenVectors();
   fEigenValues  = eigen.GetEigenValues();
   //make sure that eigenvalues are positive
   for (Int_t i = 0; i < fNumberOfVariables; i++) {
      if (fEigenValues[i] < 0) fEigenValues[i] = -fEigenValues[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This is the method that actually generates the code for the
/// transformations to and from feature space and pattern space
/// It's called by TPrincipal::MakeCode and TPrincipal::MakeMethods.
///
/// The options are: NONE so far

void TPrincipal::MakeRealCode(const char *filename, const char *classname,
                              Option_t *)
{
   Bool_t  isMethod = (classname[0] == '\0' ? kFALSE : kTRUE);
   const char *prefix   = (isMethod ? Form("%s::", classname) : "");
   const char *cv_qual  = (isMethod ? "" : "static ");

   std::ofstream outFile(filename,std::ios::out|std::ios::trunc);
   if (!outFile) {
      Error("MakeRealCode","couldn't open output file '%s'",filename);
      return;
   }

   std::cout << "Writing on file \"" << filename << "\" ... " << std::flush;
   //
   // Write header of file
   //
   // Emacs mode line ;-)
   outFile << "// -*- mode: c++ -*-" << std::endl;
   // Info about creator
   outFile << "// " << std::endl
      << "// File " << filename
      << " generated by TPrincipal::MakeCode" << std::endl;
   // Time stamp
   TDatime date;
   outFile << "// on " << date.AsString() << std::endl;
   // ROOT version info
   outFile << "// ROOT version " << gROOT->GetVersion()
      << std::endl << "//" << std::endl;
   // General information on the code
   outFile << "// This file contains the functions " << std::endl
      << "//" << std::endl
      << "//    void  " << prefix
      << "X2P(Double_t *x, Double_t *p); " << std::endl
      << "//    void  " << prefix
      << "P2X(Double_t *p, Double_t *x, Int_t nTest);"
      << std::endl << "//" << std::endl
      << "// The first for transforming original data x in " << std::endl
      << "// pattern space, to principal components p in " << std::endl
      << "// feature space. The second function is for the" << std::endl
      << "// inverse transformation, but using only nTest" << std::endl
      << "// of the principal components in the expansion" << std::endl
      << "// " << std::endl
      << "// See TPrincipal class documentation for more "
      << "information " << std::endl << "// " << std::endl;
   // Header files
   outFile << "#ifndef __CINT__" << std::endl;
   if (isMethod)
      // If these are methods, we need the class header
      outFile << "#include \"" << classname << ".h\"" << std::endl;
   else
      // otherwise, we need the typedefs of Int_t and Double_t
      outFile << "#include <Rtypes.h> // needed for Double_t etc" << std::endl;
   // Finish the preprocessor block
   outFile << "#endif" << std::endl << std::endl;

   //
   // Now for the data
   //
   // We make the Eigenvector matrix, Eigenvalue vector, Sigma vector,
   // and Mean value vector static, since all are needed in both
   // functions. Also ,the number of variables are stored in a static
   // variable.
   outFile << "//" << std::endl
      << "// Static data variables"  << std::endl
      << "//" << std::endl;
   outFile << cv_qual << "Int_t    " << prefix << "gNVariables = "
      << fNumberOfVariables << ";" << std::endl;

   // Assign the values to the Eigenvector matrix. The elements are
   // stored row-wise, that is element
   //    M[i][j] = e[i * nVariables + j]
   // where i and j are zero-based.
   outFile << std::endl << "// Assignment of eigenvector matrix." << std::endl
      << "// Elements are stored row-wise, that is" << std::endl
      << "//    M[i][j] = e[i * nVariables + j] " << std::endl
      << "// where i and j are zero-based" << std::endl;
   outFile << cv_qual << "Double_t " << prefix
      << "gEigenVectors[] = {" << std::flush;
   Int_t i,j;
   for (i = 0; i < fNumberOfVariables; i++) {
      for (j = 0; j < fNumberOfVariables; j++) {
         Int_t index = i * fNumberOfVariables + j;
         outFile << (index != 0 ? "," : "" ) << std::endl
            << "  "  << fEigenVectors(i,j) << std::flush;
      }
   }
   outFile << "};" << std::endl << std::endl;

   // Assignment to eigenvalue vector. Zero-based.
   outFile << "// Assignment to eigen value vector. Zero-based." << std::endl;
   outFile << cv_qual << "Double_t " << prefix
      << "gEigenValues[] = {" << std::flush;
   for (i = 0; i < fNumberOfVariables; i++)
      outFile << (i != 0 ? "," : "") << std::endl
      << "  " << fEigenValues(i) << std::flush;
   outFile << std::endl << "};" << std::endl << std::endl;

   // Assignment to mean Values vector. Zero-based.
   outFile << "// Assignment to mean value vector. Zero-based." << std::endl;
   outFile << cv_qual << "Double_t " << prefix
      << "gMeanValues[] = {" << std::flush;
   for (i = 0; i < fNumberOfVariables; i++)
      outFile << (i != 0 ? "," : "") << std::endl
      << "  " << fMeanValues(i) << std::flush;
   outFile << std::endl << "};" << std::endl << std::endl;

   // Assignment to mean Values vector. Zero-based.
   outFile << "// Assignment to sigma value vector. Zero-based." << std::endl;
   outFile << cv_qual << "Double_t " << prefix
      << "gSigmaValues[] = {" << std::flush;
   for (i = 0; i < fNumberOfVariables; i++)
      outFile << (i != 0 ? "," : "") << std::endl
      << "  " << (fIsNormalised ? fSigmas(i) : 1) << std::flush;
   //    << "  " << fSigmas(i) << std::flush;
   outFile << std::endl << "};" << std::endl << std::endl;

   //
   // Finally we reach the functions themselves
   //
   // First: void x2p(Double_t *x, Double_t *p);
   //
   outFile << "// " << std::endl
      << "// The "
      << (isMethod ? "method " : "function ")
      << "  void " << prefix
      << "X2P(Double_t *x, Double_t *p)"
      << std::endl << "// " << std::endl;
   outFile << "void " << prefix
      << "X2P(Double_t *x, Double_t *p) {" << std::endl
      << "  for (Int_t i = 0; i < gNVariables; i++) {" << std::endl
      << "    p[i] = 0;" << std::endl
      << "    for (Int_t j = 0; j < gNVariables; j++)" << std::endl
      << "      p[i] += (x[j] - gMeanValues[j]) " << std::endl
      << "        * gEigenVectors[j *  gNVariables + i] "
      << "/ gSigmaValues[j];" << std::endl << std::endl << "  }"
      << std::endl << "}" << std::endl << std::endl;
   //
   // Now: void p2x(Double_t *p, Double_t *x, Int_t nTest);
   //
   outFile << "// " << std::endl << "// The "
      << (isMethod ? "method " : "function ")
      << "  void " << prefix
      << "P2X(Double_t *p, Double_t *x, Int_t nTest)"
      << std::endl << "// " << std::endl;
   outFile << "void " << prefix
      << "P2X(Double_t *p, Double_t *x, Int_t nTest) {" << std::endl
      << "  for (Int_t i = 0; i < gNVariables; i++) {" << std::endl
      << "    x[i] = gMeanValues[i];" << std::endl
      << "    for (Int_t j = 0; j < nTest; j++)" << std::endl
      << "      x[i] += p[j] * gSigmaValues[i] " << std::endl
      << "        * gEigenVectors[i *  gNVariables + j];" << std::endl
      << "  }" << std::endl << "}" << std::endl << std::endl;

   // EOF
   outFile << "// EOF for " << filename << std::endl;

   // Close the file
   outFile.close();

   std::cout << "done" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate x as a function of nTest of the most significant
/// principal components p, and return it in x.
/// It's the users responsibility to make sure that both x and p are
/// of the right size (i.e., memory must be allocated for x).

void TPrincipal::P2X(const Double_t *p, Double_t *x, Int_t nTest)
{
   for (Int_t i = 0; i < fNumberOfVariables; i++){
      x[i] = fMeanValues(i);
      for (Int_t j = 0; j < nTest; j++)
         x[i] += p[j] * (fIsNormalised ? fSigmas(i) : 1)
         * fEigenVectors(i,j);
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Print the statistics
/// Options are
///  - M            Print mean values of original data
///  - S            Print sigma values of original data
///  - E            Print eigenvalues of covariance matrix
///  - V            Print eigenvectors of covariance matrix
/// Default is MSE

void TPrincipal::Print(Option_t *opt) const
{
   Bool_t printV = kFALSE;
   Bool_t printM = kFALSE;
   Bool_t printS = kFALSE;
   Bool_t printE = kFALSE;

   Int_t len     = strlen(opt);
   for (Int_t i = 0; i < len; i++) {
      switch (opt[i]) {
         case 'V':
         case 'v':
            printV = kTRUE;
            break;
         case 'M':
         case 'm':
            printM = kTRUE;
            break;
         case 'S':
         case 's':
            printS = kTRUE;
            break;
         case 'E':
         case 'e':
            printE = kTRUE;
            break;
         default:
            Warning("Print", "Unknown option '%c'",opt[i]);
            break;
      }
   }

   if (printM||printS||printE) {
      std::cout << " Variable #  " << std::flush;
      if (printM)
         std::cout << "| Mean Value " << std::flush;
      if (printS)
         std::cout << "|   Sigma    " << std::flush;
      if (printE)
         std::cout << "| Eigenvalue" << std::flush;
      std::cout << std::endl;

      std::cout << "-------------" << std::flush;
      if (printM)
         std::cout << "+------------" << std::flush;
      if (printS)
         std::cout << "+------------" << std::flush;
      if (printE)
         std::cout << "+------------" << std::flush;
      std::cout << std::endl;

      for (Int_t i = 0; i < fNumberOfVariables; i++) {
         std::cout << std::setw(12) << i << " " << std::flush;
         if (printM)
            std::cout << "| " << std::setw(10) << std::setprecision(4)
            << fMeanValues(i) << " " << std::flush;
         if (printS)
            std::cout << "| " << std::setw(10) << std::setprecision(4)
            << fSigmas(i) << " " << std::flush;
         if (printE)
            std::cout << "| " << std::setw(10) << std::setprecision(4)
            << fEigenValues(i) << " " << std::flush;
         std::cout << std::endl;
      }
      std::cout << std::endl;
   }

   if(printV) {
      for (Int_t i = 0; i < fNumberOfVariables; i++) {
         std::cout << "Eigenvector # " << i << std::flush;
         TVectorD v(fNumberOfVariables);
         v = TMatrixDColumn_const(fEigenVectors,i);
         v.Print();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the sum of the square residuals, that is
///
/// \f[
/// E_N = \sum_{i=0}^{P-1} \left(x_i - x^\prime_i\right)^2
/// \f]
///
/// where \f$x^\prime_i = \sum_{j=i}^N p_i e_{n_j}\f$
/// is the \f$i^{\mbox{th}}\f$ component of the principal vector, corresponding to
/// \f$x_i\f$, the original data; I.e., the square distance to the space
/// spanned by \f$N\f$ eigenvectors.

void TPrincipal::SumOfSquareResiduals(const Double_t *x, Double_t *s)
{

   if (!x)
      return;

   Double_t p[100];
   Double_t xp[100];

   X2P(x,p);
   for (Int_t i = fNumberOfVariables-1; i >= 0; i--) {
      P2X(p,xp,i);
      for (Int_t j = 0; j < fNumberOfVariables; j++) {
         s[i] += (x[j] - xp[j])*(x[j] - xp[j]);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Test the PCA, bye calculating the sum square of residuals
/// (see method SumOfSquareResiduals), and display the histogram

void TPrincipal::Test(Option_t *)
{
   MakeHistograms("pca","S");

   if (!fStoreData)
      return;

   TH1 *pca_s = 0;
   if (fHistograms) pca_s = (TH1*)fHistograms->FindObject("pca_s");
   if (!pca_s) {
      Warning("Test", "Couldn't get histogram of square residuals");
      return;
   }

   pca_s->Draw();
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the principal components from the original data vector
/// x, and return it in p.
///
/// It's the users responsibility to make sure that both x and p are
/// of the right size (i.e., memory must be allocated for p).

void TPrincipal::X2P(const Double_t *x, Double_t *p)
{
   for (Int_t i = 0; i < fNumberOfVariables; i++){
      p[i] = 0;
      for (Int_t j = 0; j < fNumberOfVariables; j++)
         p[i] += (x[j] - fMeanValues(j)) * fEigenVectors(j,i) /
         (fIsNormalised ? fSigmas(j) : 1);
   }

}
