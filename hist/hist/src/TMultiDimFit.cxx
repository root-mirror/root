// @(#)root/hist:$Id$
// Author: Christian Holm Christensen 07/11/2000

/** \class TMultiDimFit
    \ingroup Hist

 Multidimensional Fits in ROOT.
 ## Overview
 A common problem encountered in different fields of applied science is
 to find an expression for one physical quantity in terms of several
 others, which are directly measurable.

 An example in high energy physics is the evaluation of the momentum of
 a charged particle from the observation of its trajectory in a magnetic
 field.  The problem is to relate the momentum of the particle to the
 observations, which may consists of positional measurements at
 intervals along the particle trajectory.

 The exact functional relationship between the measured quantities
 (e.g., the space-points) and the dependent quantity (e.g., the
 momentum) is in general not known, but one possible way of solving the
 problem, is to find an expression which reliably approximates the
 dependence of the momentum on the observations.

 This explicit function of the observations can be obtained by a
 <I>least squares</I> fitting procedure applied to a representative
 sample of the data, for which the dependent quantity (e.g., momentum)
 and the independent observations are known. The function can then be
 used to compute the quantity of interest for new observations of the
 independent variables.

 This class <TT>TMultiDimFit</TT> implements such a procedure in
 ROOT. It is largely based on the CERNLIB MUDIFI package [2].
 Though the basic concepts are still sound, and
 therefore kept, a few implementation details have changed, and this
 class can take advantage of MINUIT [4] to improve the errors
 of the fitting, thanks to the class TMinuit.

 In [5] and [6] H. Wind demonstrates the utility
 of this procedure in the context of tracking, magnetic field
 parameterisation, and so on. The outline of the method used in this
 class is based on Winds discussion, and I refer these two excellents
 text for more information.

 And example of usage is given in multidimfit.C.

 ## The Method
 Let \f$ D \f$ by the dependent quantity of interest, which depends smoothly
 on the observable quantities \f$ x_1, \ldots, x_N \f$ which we'll denote by
 \f$\mathbf{x}\f$. Given a training sample of \f$ M\f$ tuples of the form, (TMultiDimFit::AddRow)

 \f[
     \left(\mathbf{x}_j, D_j, E_j\right)\quad,
 \f]
 where \f$\mathbf{x}_j = (x_{1,j},\ldots,x_{N,j})\f$ are \f$ N\f$ independent
 variables, \f$ D_j\f$ is the known, quantity dependent at \f$\mathbf{x}_j\f$ and \f$ E_j\f$ is
 the square error in \f$ D_j\f$, the class will try to find the parameterization
 \f[
     D_p(\mathbf{x}) = \sum_{l=1}^{L} c_l \prod_{i=1}^{N} p_{li}\left(x_i\right)
     = \sum_{l=1}^{L} c_l F_l(\mathbf{x})
 \f]
 such that

 \f[
     S \equiv \sum_{j=1}^{M} \left(D_j - D_p\left(\mathbf{x}_j\right)\right)^2
 \f]
 is minimal. Here \f$p_{li}(x_i)\f$ are monomials, or Chebyshev or Legendre
 polynomials, labelled \f$l = 1, \ldots, L\f$, in each variable \f$ x_i\f$,\f$ i=1, \ldots, N\f$.

 So what TMultiDimFit does, is to determine the number of terms \f$ L\f$, and then \f$ L\f$ terms
 (or functions) \f$ F_l\f$, and the \f$ L\f$ coefficients \f$ c_l\f$, so that \f$ S\f$ is minimal
 (TMultiDimFit::FindParameterization).

 Of course it's more than a little unlikely that \f$ S\f$ will ever become
 exact zero as a result of the procedure outlined below. Therefore, the
 user is asked to provide a minimum relative error \f$ \epsilon\f$ (TMultiDimFit::SetMinRelativeError),
 and \f$ S\f$ will be considered minimized when

 \f[
   R = \frac{S}{\sum_{j=1}^M D_j^2} < \epsilon
 \f]
 Optionally, the user may impose a functional expression by specifying
 the powers of each variable in \f$ L\f$ specified functions \f$ F_1, \ldots,F_L\f$ (TMultiDimFit::SetPowers).
 In that case, only the coefficients \f$ c_l\f$ is calculated by the class.

 ## Limiting the Number of Terms
 As always when dealing with fits, there's a real chance of *over fitting*. As is well-known, it's
 always possible to fit an \f$ N-1\f$ polynomial in \f$ x\f$ to \f$ N\f$ points \f$ (x,y)\f$ with
 \f$\chi^2 = 0\f$, but the polynomial is not likely to fit new data at all [1].
 Therefore, the user is asked to provide an upper limit, \f$ L_{max}\f$ to the number of terms in
 \f$ D_p\f$ (TMultiDimFit::SetMaxTerms).

 However, since there's an infinite number of \f$ F_l\f$ to choose from, the
 user is asked to give the maximum power. \f$ P_{max,i}\f$, of each variable
 \f$ x_i\f$ to be considered in the minimization of \f$ S\f$ (TMultiDimFit::SetMaxPowers).

 One way of obtaining values for the maximum power in variable \f$ i\f$, is
 to perform a regular fit to the dependent quantity \f$ D\f$, using a
 polynomial only in \f$ x_i\f$. The maximum power is \f$ P_{max,i}\f$ is then the
 power that does not significantly improve the one-dimensional
 least-square fit over \f$ x_i\f$ to \f$ D\f$ [5].

 There are still a huge amount of possible choices for \f$ F_l\f$; in fact
 there are \f$\prod_{i=1}^{N} (P_{max,i} + 1)\f$ possible
 choices. Obviously we need to limit this. To this end, the user is
 asked to set a *power control limit*, \f$ Q\f$ (TMultiDimFit::SetPowerLimit), and a function
 \f$ F_l\f$ is only accepted if
 \f[
   Q_l = \sum_{i=1}^{N} \frac{P_{li}}{P_{max,i}} < Q
 \f]
 where \f$ P_{li}\f$ is the leading power of variable \f$ x_i\f$ in function \f$ F_l\f$ (TMultiDimFit::MakeCandidates).
 So the number of functions increase with \f$ Q\f$ (1, 2 is fine, 5 is way out).

 ## Gram-Schmidt Orthogonalisation</A>
 To further reduce the number of functions in the final expression,
 only those functions that significantly reduce \f$ S\f$ is chosen. What
 `significant' means, is chosen by the user, and will be
 discussed below (see [2.3](TMultiFimFit.html#sec:selectiondetail)).

 The functions \f$ F_l\f$ are generally not orthogonal, which means one will
 have to evaluate all possible \f$ F_l\f$'s over all data-points before
 finding the most significant [1]. We can, however, do
 better then that. By applying the *modified Gram-Schmidt
 orthogonalisation* algorithm [5] [3] to the
 functions \f$ F_l\f$, we can evaluate the contribution to the reduction of
 \f$ S\f$ from each function in turn, and we may delay the actual inversion
 of the curvature-matrix (TMultiDimFit::MakeGramSchmidt).

 So we are let to consider an \f$ M\times L\f$ matrix \f$\mathsf{F}\f$, an
 element of which is given by
 \f[
   f_{jl} = F_j\left(x_{1j} , x_{2j}, \ldots, x_{Nj}\right)
   = F_l(\mathbf{x}_j)\,  \quad\mbox{with}~j=1,2,\ldots,M,
 \f]
 where \f$ j\f$ labels the \f$ M\f$ rows in the training sample and \f$ l\f$ labels
 \f$ L\f$ functions of \f$ N\f$ variables, and \f$ L \leq M\f$. That is, \f$ f_{jl}\f$ is
 the term (or function) numbered \f$ l\f$ evaluated at the data point
 \f$ j\f$. We have to normalise \f$\mathbf{x}_j\f$ to \f$ [-1,1]\f$ for this to
 succeed [5] (TMultiDimFit::MakeNormalized). We then define a
 matrix \f$\mathsf{W}\f$ of which the columns \f$\mathbf{w}_j\f$ are given by
 \f{eqnarray*}{
   \mathbf{w}_1 &=& \mathbf{f}_1 = F_1\left(\mathbf x_1\right)\\
   \mathbf{w}_l &=& \mathbf{f}_l - \sum^{l-1}_{k=1} \frac{\mathbf{f}_l \bullet
   \mathbf{w}_k}{\mathbf{w}_k^2}\mathbf{w}_k\,.
 \f}
 and \f$\mathbf{w}_{l}\f$ is the component of \f$\mathbf{f}_{l} \f$ orthogonal
 to \f$\mathbf{w}_{1}, \ldots, \mathbf{w}_{l-1}\f$. Hence we obtain [3],
 \f[
   \mathbf{w}_k\bullet\mathbf{w}_l = 0\quad\mbox{if}~k \neq l\quad.
 \f]
 We now take as a new model \f$\mathsf{W}\mathbf{a}\f$. We thus want to
 minimize
 \f[
   S\equiv \left(\mathbf{D} - \mathsf{W}\mathbf{a}\right)^2\quad,
 \f]
 where \f$\mathbf{D} = \left(D_1,\ldots,D_M\right)\f$ is a vector of the
 dependent quantity in the sample. Differentiation with respect to
 \f$ a_j\f$ gives, using [6], <a name="eq:dS2"></a>
 \f[
   \mathbf{D}\bullet\mathbf{w}_l - a_l\mathbf{w}_l^2 = 0
 \f]
 or
 \f[
   a_l = \frac{\mathbf{D}_l\bullet\mathbf{w}_l}{\mathbf{w}_l^2}
 \f]
 Let \f$ S_j\f$ be the sum of squares of residuals when taking \f$ j\f$ functions
 into account. Then
 \f[
   S_l = \left[\mathbf{D} - \sum^l_{k=1} a_k\mathbf{w}_k\right]^2
   = \mathbf{D}^2 - 2\mathbf{D} \sum^l_{k=1} a_k\mathbf{w}_k
   + \sum^l_{k=1} a_k^2\mathbf{w}_k^2
 \f]
 Using [9], we see that
 \f{eqnarray*}{
   S_l &=& \mathbf{D}^2 - 2 \sum^l_{k=1} a_k^2\mathbf{w}_k^2 +
   \sum^j_{k=1} a_k^2\mathbf{w}_k^2\nonumber\\
   &=& \mathbf{D}^2 - \sum^l_{k=1} a_k^2\mathbf{w}_k^2\nonumber\\
   &=& \mathbf{D}^2 - \sum^l_{k=1} \frac{\left(\mathbf D\bullet \mathbf
   w_k\right)}{\mathbf w_k^2}
 \f}
 So for each new function \f$ F_l\f$ included in the model, we get a
 reduction of the sum of squares of residuals of \f$a_l^2\mathbf{w}_l^2\f$,
 where \f$\mathbf{w}_l\f$ is given by [4] and \f$ a_l\f$ by [9]. Thus, using
 the Gram-Schmidt orthogonalisation, we
 can decide if we want to include this function in the final model,
 *before* the matrix inversion.

 ## Function Selection Based on Residual
 Supposing that \f$ L-1\f$ steps of the procedure have been performed, the
 problem now is to consider the \f$L^{\mbox{th}}\f$ function.

 The sum of squares of residuals can be written as
 \f[
   S_L = \textbf{D}^T\bullet\textbf{D} -
   \sum^L_{l=1}a^2_l\left(\textbf{w}_l^T\bullet\textbf{w}_l\right)
 \f]
 where the relation [9] have been taken into account. The
 contribution of the \f$L^{\mbox{th}}\f$ function to the reduction of S, is
 given by
 \f[
   \Delta S_L = a^2_L\left(\textbf{w}_L^T\bullet\textbf{w}_L\right)
 \f]
 Two test are now applied to decide whether this \f$L^{\mbox{th}}\f$
 function is to be included in the final expression, or not.

 ## Test 1
 Denoting by \f$ H_{L-1}\f$ the subspace spanned by \f$\textbf{w}_1,\ldots,\textbf{w}_{L-1}\f$
 the function \d$\textbf{w}_L\d$ is by construction (see 4) the projection of the function
 \f$ F_L\f$ onto the direction perpendicular to \f$ H_{L-1}\f$. Now, if the
 length of \f$\textbf{w}_L\f$ (given by \f$\textbf{w}_L\bullet\textbf{w}_L\f$)
 is very small compared to the length of \f$\textbf{f}_L\f$ this new
 function can not contribute much to the reduction of the sum of
 squares of residuals. The test consists then in calculating the angle
 \f$ \theta \f$ between the two vectors \f$\textbf{w}_L\f$ \f$ \textbf {f}_L\f$
 (see also figure 1) and requiring that it's
 *greater* then a threshold value which the user must set (TMultiDimFit::SetMinAngle).

 \image html multidimfit_img86.gif "Figure 1: (a) angle \f$\theta\f$ between \f$\textbf{w}_l\f$ and \f$\textbf{f}_L\f$, (b) angle \f$ \phi \f$ between \f$\textbf{w}_L\f$ and \f$\textbf{D}\f$"

 ## Test 2
 Let \f$\textbf{D}\f$ be the data vector to be fitted. As illustrated in
 figure 1, the \f$L^{\mbox{th}}\f$ function \f$\textbf{w}_L\f$
 will contribute significantly to the reduction of \f$ S\f$, if the angle
 \f$\phi^\prime\f$ between \f$\textbf{w}_L\f$ and \f$\textbf{D}\f$ is smaller than
 an upper limit \f$ \phi \f$, defined by the user (MultiDimFit::SetMaxAngle)

 However, the method automatically readjusts the value of this angle
 while fitting is in progress, in order to make the selection criteria
 less and less difficult to be fulfilled. The result is that the
 functions contributing most to the reduction of \f$ S\f$ are chosen first
 (TMultiDimFit::TestFunction).

 In case \f$ \phi \f$ isn't defined, an alternative method of
 performing this second test is used: The \f$L^{\mbox{th}}\f$
 function \f$\textbf{f}_L\f$ is accepted if (refer also to equation (13))
 \f[
   \Delta S_L > \frac{S_{L-1}}{L_{max}-L}
 \f]
 where  \f$ S_{L-1}\f$ is the sum of the \f$ L-1\f$ first residuals from the
 \f$ L-1\f$ functions previously accepted; and \f$ L_{max}\f$ is the total number
 of functions allowed in the final expression of the fit (defined by
 user).

 From this we see, that by restricting \f$ L_{max}\f$ -- the number of
 terms in the final model -- the fit is more difficult to perform,
 since the above selection criteria is more limiting.

 The more coefficients we evaluate, the more the sum of squares of
 residuals \f$ S\f$ will be reduced. We can evaluate \f$ S\f$ before inverting
 \f$\mathsf{B}\f$ as shown below.

 ## Coefficients and Coefficient Errors
 Having found a parameterization, that is the \f$ F_l\f$'s and \f$ L\f$, that
 minimizes \f$ S\f$, we still need to determine the coefficients
 \f$ c_l\f$. However, it's a feature of how we choose the significant
 functions, that the evaluation of the \f$ c_l\f$'s becomes trivial [5]. To derive
 \f$\mathbf{c}\f$, we first note that
 equation (4) can be written as
 \f[
   \mathsf{F} = \mathsf{W}\mathsf{B}
 \f]
 where
 \f{eqnarray*}{
   b_{ij} = \frac{\mathbf{f}_j \bullet \mathbf{w}_i}{\mathbf{w}_i^2}
     & \mbox{if} & i < j\\
   1 & \mbox{if} & i = j\\
   0 & \mbox{if} & i > j
 \f}
 Consequently, \f$\mathsf{B}\f$ is an upper triangle matrix, which can be
 readily inverted. So we now evaluate
 \f[
   \mathsf{F}\mathsf{B}^{-1} = \mathsf{W}
 \f]
 The model \f$\mathsf{W}\mathbf{a}\f$ can therefore be written as
 \f$(\mathsf{F}\mathsf{B}^{-1})\mathbf{a} = \mathsf{F}(\mathsf{B}^{-1}\mathbf{a})\,.\f$

 The original model \f$\mathsf{F}\mathbf{c}\f$ is therefore identical with
 this if
 \f[
   \mathbf{c} = \left(\mathsf{B}^{-1}\mathbf{a}\right) =
   \left[\mathbf{a}^T\left(\mathsf{B}^{-1}\right)^T\right]^T\,.
 \f]
 The reason we use \f$\left(\mathsf{B}^{-1}\right)^T\f$ rather then
 \f$\mathsf{B}^{-1}\f$ is to save storage, since \f$\left(\mathsf{B}^{-1}\right)^T\f$
 can be stored in the same matrix as \f$\mathsf{B}\f$ (TMultiDimFit::MakeCoefficients).
 The errors in the coefficients is calculated by inverting the curvature matrix
 of the non-orthogonal functions \f$ f_{lj}\f$ [1] (TMultiDimFit::MakeCoefficientErrors).

 ## Considerations
 It's important to realize that the training sample should be
 representative of the problem at hand, in particular along the borders
 of the region of interest. This is because the algorithm presented
 here, is a *interpolation*, rather then a *extrapolation* [5].

 Also, the independent variables \f$ x_{i}\f$ need to be linear
 independent, since the procedure will perform poorly if they are not
 [5]. One can find an linear transformation from ones
 original variables \f$ \xi_{i}\f$ to a set of linear independent variables
 \f$ x_{i}\f$, using a *Principal Components Analysis* (see TPrincipal), and
 then use the transformed variable as input to this class [5] [6].

 H. Wind also outlines a method for parameterising a multidimensional
 dependence over a multidimensional set of variables. An example
 of the method from [5], is a follows (please refer to
 [5] for a full discussion):

 1. Define \f$\mathbf{P} = (P_1, \ldots, P_5)\f$ are the 5 dependent
 quantities that define a track.
 2. Compute, for \f$ M\f$ different values of \f$\mathbf{P}\f$, the tracks
 through the magnetic field, and determine the corresponding
 \f$\mathbf{x} = (x_1, \ldots, x_N)\f$.
 3. Use the simulated observations to determine, with a simple
 approximation, the values of \f$\mathbf{P}_j\f$. We call these values
 \f$\mathbf{P}^\prime_j, j = 1, \ldots, M\f$.
 4. Determine from \f$\mathbf{x}\f$ a set of at least five relevant
 coordinates \f$\mathbf{x}^\prime\f$, using contrains, *or
 alternative:*
 5. Perform a Principal Component Analysis (using TPrincipal), and use
 to get a linear transformation \f$\mathbf{x} \rightarrow \mathbf{x}^\prime\f$, so that
 \f$\mathbf{x}^\prime\f$ are constrained and linear independent.
 6. Perform a Principal Component Analysis on
 \f$Q_i = P_i / P^\prime_i\, i = 1, \ldots, 5\f$, to get linear
 indenpendent (among themselves, but not independent of \f$\mathbf{x}\f$) quantities
 \f$\mathbf{Q}^\prime\f$
 7. For each component \f$Q^\prime_i\f$ make a multidimensional fit,
 using \f$\mathbf{x}^\prime\f$ as the variables, thus determining a set of
 coefficients \f$\mathbf{c}_i\f$.

 To process data, using this parameterisation, do
 1. Test wether the observation \f$\mathbf{x}\f$ within the domain of
 the parameterization, using the result from the Principal Component
 Analysis.
 2. Determine \f$\mathbf{P}^\prime\f$ as before.
 3. Determine \f$\mathbf{x}^\prime\f$ as before.
 4. Use the result of the fit to determine \f$\mathbf{Q}^\prime\f$.
 5. Transform back to \f$\mathbf{P}\f$ from \f$\mathbf{Q}^\prime\f$, using
 the result from the Principal Component Analysis.

 ## Testing the parameterization
 The class also provides functionality for testing the, over the
 training sample, found parameterization (TMultiDimFit::Fit). This is done by passing
 the class a test sample of \f$ M_t\f$ tuples of the form
 \f$(\mathbf{x}_{t,j},D_{t,j}, E_{t,j})\f$, where \f$\mathbf{x}_{t,j}\f$ are the independent
 variables, \f$ D_{t,j}\f$ the known, dependent quantity, and \f$ E_{t,j}\f$ is
 the square error in \f$ D_{t,j}\f$ (TMultiDimFit::AddTestRow).

 The parameterization is then evaluated at every \f$\mathbf{x}_t\f$ in the
 test sample, and
 \f[
   S_t \equiv \sum_{j=1}^{M_t} \left(D_{t,j} -
   D_p\left(\mathbf{x}_{t,j}\right)\right)^2
 \f]
 is evaluated. The relative error over the test sample
 \f[
   R_t = \frac{S_t}{\sum_{j=1}^{M_t} D_{t,j}^2}
 \f]
 should not be to low or high compared to \f$ R\f$ from the training
 sample. Also, multiple correlation coefficient from both samples should
 be fairly close, otherwise one of the samples is not representative of
 the problem. A large difference in the reduced \f$ \chi^2\f$ over the two
 samples indicate an over fit, and the maximum number of terms in the
 parameterisation should be reduced.

 It's possible to use [4] to further improve the fit, using the test sample.

 Christian Holm

 ## Bibliography
 - <a name="bevington"></a> Philip R. Bevington and D. Keith Robinson. *Data Reduction and Error Analysis for
   the Physical Sciences*. McGraw-Hill, 2 edition, 1992.
 - <a name="mudifi"></a> R. Brun et al. *Long writeup DD/75-23*, CERN, 1980.
 - Gene H. Golub and Charles F. van Loan. *Matrix Computations*.
   John Hopkins University Press, Baltimore, 3 edition, 1996.
 - <a name="minuit"></a>F. James. *Minuit*. Long writeup D506, CERN, 1998.
 - <a name="wind72"></a>H. Wind. *Function parameterization*. Proceedings of the 1972 CERN Computing and Data Processing
   School, volume 72-21 of Yellow report. CERN, 1972.
 - <a name="wind81"></a>H. Wind. 1. principal component analysis, 2. pattern recognition for track
   finding, 3. interpolation and functional representation. Yellow report EP/81-12, CERN, 1981.

[1]: classTMultiDimFit.html#bevington
[2]: classTMultiDimFit.html#mudifi
[4]: classTMultiDimFit.html#minuit
[5]: classTMultiDimFit.html#wind72
[6]: classTMultiDimFit.html#wind81
[9]: classTMultiDimFit.html#eq:dS2
*/


#include "Riostream.h"
#include "TMultiDimFit.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TList.h"
#include "TBrowser.h"
#include "TDecompChol.h"
#include "TDatime.h"


#define RADDEG (180. / TMath::Pi())
#define DEGRAD (TMath::Pi() / 180.)
#define HIST_XORIG     0
#define HIST_DORIG     1
#define HIST_XNORM     2
#define HIST_DSHIF     3
#define HIST_RX        4
#define HIST_RD        5
#define HIST_RTRAI     6
#define HIST_RTEST     7
#define PARAM_MAXSTUDY 1
#define PARAM_SEVERAL  2
#define PARAM_RELERR   3
#define PARAM_MAXTERMS 4


////////////////////////////////////////////////////////////////////////////////

static void mdfHelper(int&, double*, double&, double*, int);

////////////////////////////////////////////////////////////////////////////////

ClassImp(TMultiDimFit);

//____________________________________________________________________
// Static instance. Used with mdfHelper and TMinuit
TMultiDimFit* TMultiDimFit::fgInstance = 0;


////////////////////////////////////////////////////////////////////////////////
/// Empty CTOR. Do not use

TMultiDimFit::TMultiDimFit()
{
   fMeanQuantity           = 0;
   fMaxQuantity            = 0;
   fMinQuantity            = 0;
   fSumSqQuantity          = 0;
   fSumSqAvgQuantity       = 0;

   fNVariables             = 0;
   fSampleSize             = 0;
   fTestSampleSize         = 0;

   fMinAngle               = 1;
   fMaxAngle               = 0;
   fMaxTerms               = 0;
   fMinRelativeError       = 0;
   fMaxPowers              = 0;
   fPowerLimit             = 0;

   fMaxFunctions           = 0;
   fFunctionCodes          = 0;
   fMaxStudy               = 0;
   fMaxFuncNV              = 0;

   fMaxPowersFinal         = 0;
   fPowers                 = 0;
   fPowerIndex             = 0;

   fMaxResidual            = 0;
   fMinResidual            = 0;
   fMaxResidualRow         = 0;
   fMinResidualRow         = 0;
   fSumSqResidual          = 0;

   fNCoefficients          = 0;
   fRMS                    = 0;
   fChi2                   = 0;
   fParameterisationCode   = 0;

   fError                  = 0;
   fTestError              = 0;
   fPrecision              = 0;
   fTestPrecision          = 0;
   fCorrelationCoeff       = 0;
   fTestCorrelationCoeff   = 0;

   fHistograms             = 0;
   fHistogramMask          = 0;
   fBinVarX                = 100;
   fBinVarY                = 100;

   fFitter                 = 0;
   fPolyType               = kMonomials;
   fShowCorrelation        = kFALSE;
   fIsUserFunction         = kFALSE;
   fIsVerbose              = kFALSE;

}


////////////////////////////////////////////////////////////////////////////////
/// Constructor
/// Second argument is the type of polynomials to use in
/// parameterisation, one of:
///      TMultiDimFit::kMonomials
///      TMultiDimFit::kChebyshev
///      TMultiDimFit::kLegendre
///
/// Options:
///   K      Compute (k)correlation matrix
///   V      Be verbose
///
/// Default is no options.
///

TMultiDimFit::TMultiDimFit(Int_t dimension,
                           EMDFPolyType type,
                           Option_t *option)
: TNamed("multidimfit","Multi-dimensional fit object"),
fQuantity(dimension),
fSqError(dimension),
fVariables(dimension*100),
fMeanVariables(dimension),
fMaxVariables(dimension),
fMinVariables(dimension)
{
   fgInstance = this;

   fMeanQuantity           = 0;
   fMaxQuantity            = 0;
   fMinQuantity            = 0;
   fSumSqQuantity          = 0;
   fSumSqAvgQuantity       = 0;

   fNVariables             = dimension;
   fSampleSize             = 0;
   fTestSampleSize         = 0;

   fMinAngle               = 1;
   fMaxAngle               = 0;
   fMaxTerms               = 0;
   fMinRelativeError       = 0.01;
   fMaxPowers              = new Int_t[dimension];
   fPowerLimit             = 1;

   fMaxFunctions           = 0;
   fFunctionCodes          = 0;
   fMaxStudy               = 0;
   fMaxFuncNV              = 0;

   fMaxPowersFinal         = new Int_t[dimension];
   fPowers                 = 0;
   fPowerIndex             = 0;

   fMaxResidual            = 0;
   fMinResidual            = 0;
   fMaxResidualRow         = 0;
   fMinResidualRow         = 0;
   fSumSqResidual          = 0;

   fNCoefficients          = 0;
   fRMS                    = 0;
   fChi2                   = 0;
   fParameterisationCode   = 0;

   fError                  = 0;
   fTestError              = 0;
   fPrecision              = 0;
   fTestPrecision          = 0;
   fCorrelationCoeff       = 0;
   fTestCorrelationCoeff   = 0;

   fHistograms             = 0;
   fHistogramMask          = 0;
   fBinVarX                = 100;
   fBinVarY                = 100;

   fFitter                 = 0;
   fPolyType               = type;
   fShowCorrelation        = kFALSE;
   fIsUserFunction         = kFALSE;
   fIsVerbose              = kFALSE;
   TString opt             = option;
   opt.ToLower();

   if (opt.Contains("k")) fShowCorrelation = kTRUE;
   if (opt.Contains("v")) fIsVerbose       = kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

TMultiDimFit::~TMultiDimFit()
{
   delete [] fPowers;
   delete [] fMaxPowers;
   delete [] fMaxPowersFinal;
   delete [] fPowerIndex;
   delete [] fFunctionCodes;
   if (fHistograms) fHistograms->Clear("nodelete");
   delete fHistograms;
}


////////////////////////////////////////////////////////////////////////////////
/// Add a row consisting of fNVariables independent variables, the
/// known, dependent quantity, and optionally, the square error in
/// the dependent quantity, to the training sample to be used for the
/// parameterization.
/// The mean of the variables and quantity is calculated on the fly,
/// as outlined in TPrincipal::AddRow.
/// This sample should be representative of the problem at hand.
/// Please note, that if no error is given Poisson statistics is
/// assumed and the square error is set to the value of dependent
/// quantity.  See also the
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::AddRow(const Double_t *x, Double_t D, Double_t E)
{
   if (!x)
      return;

   if (++fSampleSize == 1) {
      fMeanQuantity  = D;
      fMaxQuantity   = D;
      fMinQuantity   = D;
      fSumSqQuantity = D * D;// G.Q. erratum on August 15th, 2008
   }
   else {
      fMeanQuantity  *= 1 - 1./Double_t(fSampleSize);
      fMeanQuantity  += D / Double_t(fSampleSize);
      fSumSqQuantity += D * D;

      if (D >= fMaxQuantity) fMaxQuantity = D;
      if (D <= fMinQuantity) fMinQuantity = D;
   }


   // If the vector isn't big enough to hold the new data, then
   // expand the vector by half it's size.
   Int_t size = fQuantity.GetNrows();
   if (fSampleSize > size) {
      fQuantity.ResizeTo(size + size/2);
      fSqError.ResizeTo(size + size/2);
   }

   // Store the value
   fQuantity(fSampleSize-1) = D;
   fSqError(fSampleSize-1) = (E == 0 ? D : E);

   // Store data point in internal vector
   // If the vector isn't big enough to hold the new data, then
   // expand the vector by half it's size
   size = fVariables.GetNrows();
   if (fSampleSize * fNVariables > size)
      fVariables.ResizeTo(size + size/2);


   // Increment the data point counter
   Int_t i,j;
   for (i = 0; i < fNVariables; i++) {
      if (fSampleSize == 1) {
         fMeanVariables(i) = x[i];
         fMaxVariables(i)  = x[i];
         fMinVariables(i)  = x[i];
      }
      else {
         fMeanVariables(i) *= 1 - 1./Double_t(fSampleSize);
         fMeanVariables(i) += x[i] / Double_t(fSampleSize);

         // Update the maximum value for this component
         if (x[i] >= fMaxVariables(i)) fMaxVariables(i)  = x[i];

         // Update the minimum value for this component
         if (x[i] <= fMinVariables(i)) fMinVariables(i)  = x[i];

      }

      // Store the data.
      j = (fSampleSize-1) * fNVariables + i;
      fVariables(j) = x[i];
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Add a row consisting of fNVariables independent variables, the
/// known, dependent quantity, and optionally, the square error in
/// the dependent quantity, to the test sample to be used for the
/// test of the parameterization.
/// This sample needn't be representative of the problem at hand.
/// Please note, that if no error is given Poisson statistics is
/// assumed and the square error is set to the value of dependent
/// quantity.  See also the
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::AddTestRow(const Double_t *x, Double_t D, Double_t E)
{
   if (fTestSampleSize++ == 0) {
      fTestQuantity.ResizeTo(fNVariables);
      fTestSqError.ResizeTo(fNVariables);
      fTestVariables.ResizeTo(fNVariables * 100);
   }

   // If the vector isn't big enough to hold the new data, then
   // expand the vector by half it's size.
   Int_t size = fTestQuantity.GetNrows();
   if (fTestSampleSize > size) {
      fTestQuantity.ResizeTo(size + size/2);
      fTestSqError.ResizeTo(size + size/2);
   }

   // Store the value
   fTestQuantity(fTestSampleSize-1) = D;
   fTestSqError(fTestSampleSize-1) = (E == 0 ? D : E);

   // Store data point in internal vector
   // If the vector isn't big enough to hold the new data, then
   // expand the vector by half it's size
   size = fTestVariables.GetNrows();
   if (fTestSampleSize * fNVariables > size)
      fTestVariables.ResizeTo(size + size/2);


   // Increment the data point counter
   Int_t i,j;
   for (i = 0; i < fNVariables; i++) {
      j = fNVariables * (fTestSampleSize - 1) + i;
      fTestVariables(j) = x[i];

      if (x[i] > fMaxVariables(i))
         Warning("AddTestRow", "variable %d (row: %d) too large: %f > %f",
                 i, fTestSampleSize, x[i], fMaxVariables(i));
      if (x[i] < fMinVariables(i))
         Warning("AddTestRow", "variable %d (row: %d) too small: %f < %f",
                 i, fTestSampleSize, x[i], fMinVariables(i));
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Browse the TMultiDimFit object in the TBrowser.

void TMultiDimFit::Browse(TBrowser* b)
{
   if (fHistograms) {
      TIter next(fHistograms);
      TH1* h = 0;
      while ((h = (TH1*)next()))
         b->Add(h,h->GetName());
   }
   if (fVariables.IsValid())
      b->Add(&fVariables, "Variables (Training)");
   if (fQuantity.IsValid())
      b->Add(&fQuantity, "Quantity (Training)");
   if (fSqError.IsValid())
      b->Add(&fSqError, "Error (Training)");
   if (fMeanVariables.IsValid())
      b->Add(&fMeanVariables, "Mean of Variables (Training)");
   if (fMaxVariables.IsValid())
      b->Add(&fMaxVariables, "Mean of Variables (Training)");
   if (fMinVariables.IsValid())
      b->Add(&fMinVariables, "Min of Variables (Training)");
   if (fTestVariables.IsValid())
      b->Add(&fTestVariables, "Variables (Test)");
   if (fTestQuantity.IsValid())
      b->Add(&fTestQuantity, "Quantity (Test)");
   if (fTestSqError.IsValid())
      b->Add(&fTestSqError, "Error (Test)");
   if (fFunctions.IsValid())
      b->Add(&fFunctions, "Functions");
   if(fCoefficients.IsValid())
      b->Add(&fCoefficients,"Coefficients");
   if(fCoefficientsRMS.IsValid())
      b->Add(&fCoefficientsRMS,"Coefficients Errors");
   if (fOrthFunctions.IsValid())
      b->Add(&fOrthFunctions, "Orthogonal Functions");
   if (fOrthFunctionNorms.IsValid())
      b->Add(&fOrthFunctionNorms, "Orthogonal Functions Norms");
   if (fResiduals.IsValid())
      b->Add(&fResiduals, "Residuals");
   if(fOrthCoefficients.IsValid())
      b->Add(&fOrthCoefficients,"Orthogonal Coefficients");
   if (fOrthCurvatureMatrix.IsValid())
      b->Add(&fOrthCurvatureMatrix,"Orthogonal curvature matrix");
   if(fCorrelationMatrix.IsValid())
      b->Add(&fCorrelationMatrix,"Correlation Matrix");
   if (fFitter)
      b->Add(fFitter, fFitter->GetName());
}


////////////////////////////////////////////////////////////////////////////////
/// Clear internal structures and variables

void TMultiDimFit::Clear(Option_t *option)
{
   Int_t i, j, n = fNVariables, m = fMaxFunctions;

   // Training sample, dependent quantity
   fQuantity.Zero();
   fSqError.Zero();
   fMeanQuantity                 = 0;
   fMaxQuantity                  = 0;
   fMinQuantity                  = 0;
   fSumSqQuantity                = 0;
   fSumSqAvgQuantity             = 0;

   // Training sample, independent variables
   fVariables.Zero();
   fNVariables                   = 0;
   fSampleSize                   = 0;
   fMeanVariables.Zero();
   fMaxVariables.Zero();
   fMinVariables.Zero();

   // Test sample
   fTestQuantity.Zero();
   fTestSqError.Zero();
   fTestVariables.Zero();
   fTestSampleSize               = 0;

   // Functions
   fFunctions.Zero();
   //for (i = 0; i < fMaxTerms; i++)  fPowerIndex[i]    = 0;
   //for (i = 0; i < fMaxTerms; i++)  fFunctionCodes[i] = 0;
   fMaxFunctions                 = 0;
   fMaxStudy                     = 0;
   fOrthFunctions.Zero();
   fOrthFunctionNorms.Zero();

   // Control parameters
   fMinRelativeError             = 0;
   fMinAngle                     = 0;
   fMaxAngle                     = 0;
   fMaxTerms                     = 0;

   // Powers
   for (i = 0; i < n; i++) {
      fMaxPowers[i]               = 0;
      fMaxPowersFinal[i]          = 0;
      for (j = 0; j < m; j++)
         fPowers[i * n + j]        = 0;
   }
   fPowerLimit                   = 0;

   // Residuals
   fMaxResidual                  = 0;
   fMinResidual                  = 0;
   fMaxResidualRow               = 0;
   fMinResidualRow               = 0;
   fSumSqResidual                = 0;

   // Fit
   fNCoefficients                = 0;
   fOrthCoefficients             = 0;
   fOrthCurvatureMatrix          = 0;
   fRMS                          = 0;
   fCorrelationMatrix.Zero();
   fError                        = 0;
   fTestError                    = 0;
   fPrecision                    = 0;
   fTestPrecision                = 0;

   // Coefficients
   fCoefficients.Zero();
   fCoefficientsRMS.Zero();
   fResiduals.Zero();
   fHistograms->Clear(option);

   // Options
   fPolyType                     = kMonomials;
   fShowCorrelation              = kFALSE;
   fIsUserFunction               = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Evaluate parameterization at point x. Optional argument coeff is
/// a vector of coefficients for the parameterisation, fNCoefficients
/// elements long.

Double_t TMultiDimFit::Eval(const Double_t *x, const Double_t* coeff) const
{
   Double_t returnValue = fMeanQuantity;
   Double_t term        = 0;
   Int_t    i, j;

   for (i = 0; i < fNCoefficients; i++) {
      // Evaluate the ith term in the expansion
      term = (coeff ? coeff[i] : fCoefficients(i));
      for (j = 0; j < fNVariables; j++) {
         // Evaluate the factor (polynomial) in the j-th variable.
         Int_t    p  =  fPowers[fPowerIndex[i] * fNVariables + j];
         Double_t y  =  1 + 2. / (fMaxVariables(j) - fMinVariables(j))
         * (x[j] - fMaxVariables(j));
         term        *= EvalFactor(p,y);
      }
      // Add this term to the final result
      returnValue += term;
   }
   return returnValue;
}


////////////////////////////////////////////////////////////////////////////////
/// Evaluate parameterization error at point x. Optional argument coeff is
/// a vector of coefficients for the parameterisation, fNCoefficients
/// elements long.

Double_t TMultiDimFit::EvalError(const Double_t *x, const Double_t* coeff) const
{
   Double_t returnValue = 0;
   Double_t term        = 0;
   Int_t    i, j;

   for (i = 0; i < fNCoefficients; i++) {
      //     std::cout << "Error coef " << i << " -> " << fCoefficientsRMS(i) << std::endl;
   }
   for (i = 0; i < fNCoefficients; i++) {
      // Evaluate the ith term in the expansion
      term = (coeff ? coeff[i] : fCoefficientsRMS(i));
      for (j = 0; j < fNVariables; j++) {
         // Evaluate the factor (polynomial) in the j-th variable.
         Int_t    p  =  fPowers[fPowerIndex[i] * fNVariables + j];
         Double_t y  =  1 + 2. / (fMaxVariables(j) - fMinVariables(j))
         * (x[j] - fMaxVariables(j));
         term        *= EvalFactor(p,y);
         //   std::cout << "i,j " << i << ", " << j << "  "  << p << "  " << y << "  " << EvalFactor(p,y) << "  " << term << std::endl;
      }
      // Add this term to the final result
      returnValue += term*term;
      //      std::cout << " i = " << i << " value = " << returnValue << std::endl;
   }
   returnValue = sqrt(returnValue);
   return returnValue;
}


////////////////////////////////////////////////////////////////////////////////
/// PRIVATE METHOD:
/// Calculate the control parameter from the passed powers

Double_t TMultiDimFit::EvalControl(const Int_t *iv) const
{
   Double_t s = 0;
   Double_t epsilon = 1e-6; // a small number
   for (Int_t i = 0; i < fNVariables; i++) {
      if (fMaxPowers[i] != 1)
         s += (epsilon + iv[i] - 1) / (epsilon + fMaxPowers[i] - 1);
   }
   return s;
}

////////////////////////////////////////////////////////////////////////////////
/// PRIVATE METHOD:
/// Evaluate function with power p at variable value x

Double_t TMultiDimFit::EvalFactor(Int_t p, Double_t x) const
{
   Int_t    i   = 0;
   Double_t p1  = 1;
   Double_t p2  = 0;
   Double_t p3  = 0;
   Double_t r   = 0;

   switch(p) {
      case 1:
         r = 1;
         break;
      case 2:
         r =  x;
         break;
      default:
         p2 = x;
         for (i = 3; i <= p; i++) {
            p3 = p2 * x;
            if (fPolyType == kLegendre)
               p3 = ((2 * i - 3) * p2 * x - (i - 2) * p1) / (i - 1);
            else if (fPolyType == kChebyshev)
               p3 = 2 * x * p2 - p1;
            p1 = p2;
            p2 = p3;
         }
         r = p3;
   }

   return r;
}


////////////////////////////////////////////////////////////////////////////////
/// Find the parameterization
///
/// Options:
///     None so far
///
/// For detailed description of what this entails, please refer to the
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::FindParameterization(Option_t *)
{
   MakeNormalized();
   MakeCandidates();
   MakeParameterization();
   MakeCoefficients();
   MakeCoefficientErrors();
   MakeCorrelation();
}

////////////////////////////////////////////////////////////////////////////////
/// Try to fit the found parameterisation to the test sample.
///
/// Options
///     M     use Minuit to improve coefficients
///
/// Also, refer to
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::Fit(Option_t *option)
{
   Int_t i, j;
   Double_t*      x    = new Double_t[fNVariables];
   Double_t  sumSqD    = 0;
   Double_t    sumD    = 0;
   Double_t  sumSqR    = 0;
   Double_t    sumR    = 0;

   // Calculate the residuals over the test sample
   for (i = 0; i < fTestSampleSize; i++) {
      for (j = 0; j < fNVariables; j++)
         x[j] = fTestVariables(i * fNVariables + j);
      Double_t res =  fTestQuantity(i) - Eval(x);
      sumD         += fTestQuantity(i);
      sumSqD       += fTestQuantity(i) * fTestQuantity(i);
      sumR         += res;
      sumSqR       += res * res;
      if (TESTBIT(fHistogramMask,HIST_RTEST))
         ((TH1D*)fHistograms->FindObject("res_test"))->Fill(res);
   }
   Double_t dAvg         = sumSqD - (sumD * sumD) / fTestSampleSize;
   Double_t rAvg         = sumSqR - (sumR * sumR) / fTestSampleSize;
   fTestCorrelationCoeff = (dAvg - rAvg) / dAvg;
   fTestError            = sumSqR;
   fTestPrecision        = sumSqR / sumSqD;

   TString opt(option);
   opt.ToLower();

   if (!opt.Contains("m"))
      MakeChi2();

   if (fNCoefficients * 50 > fTestSampleSize)
      Warning("Fit", "test sample is very small");

   if (!opt.Contains("m")) {
      Error("Fit", "invalid option");
      delete [] x;
      return;
   }

   fFitter = TVirtualFitter::Fitter(0,fNCoefficients);
   if (!fFitter) {
      Error("Fit", "Cannot create Fitter");
      delete [] x;
      return;
   }
   fFitter->SetFCN(mdfHelper);

   const Int_t  maxArgs = 16;
   Int_t           args = 1;
   Double_t*   arglist  = new Double_t[maxArgs];
   arglist[0]           = -1;
   fFitter->ExecuteCommand("SET PRINT",arglist,args);

   for (i = 0; i < fNCoefficients; i++) {
      Double_t startVal = fCoefficients(i);
      Double_t startErr = fCoefficientsRMS(i);
      fFitter->SetParameter(i, Form("coeff%02d",i),
                            startVal, startErr, 0, 0);
   }

   // arglist[0]           = 0;
   args                 = 1;
   // fFitter->ExecuteCommand("SET PRINT",arglist,args);
   fFitter->ExecuteCommand("MIGRAD",arglist,args);

   for (i = 0; i < fNCoefficients; i++) {
      Double_t val = 0, err = 0, low = 0, high = 0;
      fFitter->GetParameter(i, Form("coeff%02d",i),
                            val, err, low, high);
      fCoefficients(i)    = val;
      fCoefficientsRMS(i) = err;
   }
   delete [] x;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the static instance.

TMultiDimFit* TMultiDimFit::Instance()
{
   return fgInstance;
}

////////////////////////////////////////////////////////////////////////////////
/// PRIVATE METHOD:
/// Create list of candidate functions for the parameterisation. See
/// also
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::MakeCandidates()
{
   Int_t i = 0;
   Int_t j = 0;
   Int_t k = 0;

   // The temporary array to store the powers in. We don't need to
   // initialize this array however.
   fMaxFuncNV = fNVariables * fMaxFunctions;
   Int_t *powers = new Int_t[fMaxFuncNV];

   // store of `control variables'
   Double_t* control  = new Double_t[fMaxFunctions];

   // We've better initialize the variables
   Int_t *iv = new Int_t[fNVariables];
   for (i = 0; i < fNVariables; i++)
      iv[i] = 1;

   if (!fIsUserFunction) {

      // Number of funcs selected
      Int_t     numberFunctions = 0;

      // Absolute max number of functions
      Int_t maxNumberFunctions = 1;
      for (i = 0; i < fNVariables; i++)
         maxNumberFunctions *= fMaxPowers[i];

      while (kTRUE) {
         // Get the control value for this function
         Double_t s = EvalControl(iv);

         if (s <= fPowerLimit) {

            // Call over-loadable method Select, as to allow the user to
            // interfere with the selection of functions.
            if (Select(iv)) {
               numberFunctions++;

               // If we've reached the user defined limit of how many
               // functions we can consider, break out of the loop
               if (numberFunctions > fMaxFunctions)
                  break;

               // Store the control value, so we can sort array of powers
               // later on
               control[numberFunctions-1] = Int_t(1.0e+6*s);

               // Store the powers in powers array.
               for (i = 0; i < fNVariables; i++) {
                  j = (numberFunctions - 1) * fNVariables + i;
                  powers[j] = iv[i];
               }
            } // if (Select())
         } // if (s <= fPowerLimit)

         for (i = 0; i < fNVariables; i++)
            if (iv[i] < fMaxPowers[i])
               break;

         // If all variables have reached their maximum power, then we
         // break out of the loop
         if (i == fNVariables) {
            fMaxFunctions = numberFunctions;
            break;
         }

         // Next power in variable i
         if (i < fNVariables) iv[i]++;

         for (j = 0; j < i; j++)
            iv[j] = 1;
      } // while (kTRUE)
   }
   else {
      // In case the user gave an explicit function
      for (i = 0; i < fMaxFunctions; i++) {
         // Copy the powers to working arrays
         for (j = 0; j < fNVariables; j++) {
            powers[i * fNVariables + j] = fPowers[i * fNVariables + j];
            iv[j]                 = fPowers[i * fNVariables + j];
         }

         control[i] = Int_t(1.0e+6*EvalControl(iv));
      }
   }

   // Now we need to sort the powers according to least `control
   // variable'
   Int_t *order = new Int_t[fMaxFunctions];
   for (i = 0; i < fMaxFunctions; i++)
      order[i] = i;
   fMaxFuncNV = fMaxFunctions * fNVariables;
   fPowers = new Int_t[fMaxFuncNV];

   for (i = 0; i < fMaxFunctions; i++) {
      Double_t x = control[i];
      Int_t    l = order[i];
      k = i;

      for (j = i; j < fMaxFunctions; j++) {
         if (control[j] <= x) {
            x = control[j];
            l = order[j];
            k = j;
         }
      }

      if (k != i) {
         control[k] = control[i];
         control[i] = x;
         order[k]   = order[i];
         order[i]   = l;
      }
   }

   for (i = 0; i < fMaxFunctions; i++)
      for (j = 0; j < fNVariables; j++)
         fPowers[i * fNVariables + j] = powers[order[i] * fNVariables + j];

   delete [] control;
   delete [] powers;
   delete [] order;
   delete [] iv;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate Chi square over either the test sample. The optional
/// argument coeff is a vector of coefficients to use in the
/// evaluation of the parameterisation. If coeff == 0, then the found
/// coefficients is used.
/// Used my MINUIT for fit (see TMultDimFit::Fit)

Double_t TMultiDimFit::MakeChi2(const Double_t* coeff)
{
   fChi2 = 0;
   Int_t i, j;
   Double_t* x = new Double_t[fNVariables];
   for (i = 0; i < fTestSampleSize; i++) {
      // Get the stored point
      for (j = 0; j < fNVariables; j++)
         x[j] = fTestVariables(i * fNVariables + j);

      // Evaluate function. Scale to shifted values
      Double_t f = Eval(x,coeff);

      // Calculate contribution to Chic square
      fChi2 += 1. / TMath::Max(fTestSqError(i),1e-20)
      * (fTestQuantity(i) - f) * (fTestQuantity(i) - f);
   }

   // Clean up
   delete [] x;

   return fChi2;
}


////////////////////////////////////////////////////////////////////////////////
/// Generate the file <filename> with .C appended if argument doesn't
/// end in .cxx or .C. The contains the implementation of the
/// function:
///
///   Double_t <funcname>(Double_t *x)
///
/// which does the same as TMultiDimFit::Eval. Please refer to this
/// method.
///
/// Further, the static variables:
///
///     Int_t    gNVariables
///     Int_t    gNCoefficients
///     Double_t gDMean
///     Double_t gXMean[]
///     Double_t gXMin[]
///     Double_t gXMax[]
///     Double_t gCoefficient[]
///     Int_t    gPower[]
///
/// are initialized. The only ROOT header file needed is Rtypes.h
///
/// See TMultiDimFit::MakeRealCode for a list of options

void TMultiDimFit::MakeCode(const char* filename, Option_t *option)
{

   TString outName(filename);
   if (!outName.EndsWith(".C") && !outName.EndsWith(".cxx"))
      outName += ".C";

   MakeRealCode(outName.Data(),"",option);
}



////////////////////////////////////////////////////////////////////////////////
/// PRIVATE METHOD:
/// Compute the errors on the coefficients. For this to be done, the
/// curvature matrix of the non-orthogonal functions, is computed.

void TMultiDimFit::MakeCoefficientErrors()
{
   Int_t    i = 0;
   Int_t    j = 0;
   Int_t    k = 0;
   TVectorD iF(fSampleSize);
   TVectorD jF(fSampleSize);
   fCoefficientsRMS.ResizeTo(fNCoefficients);

   TMatrixDSym curvatureMatrix(fNCoefficients);

   // Build the curvature matrix
   for (i = 0; i < fNCoefficients; i++) {
      iF = TMatrixDRow(fFunctions,i);
      for (j = 0; j <= i; j++) {
         jF = TMatrixDRow(fFunctions,j);
         for (k = 0; k < fSampleSize; k++)
            curvatureMatrix(i,j) +=
            1 / TMath::Max(fSqError(k), 1e-20) * iF(k) * jF(k);
         curvatureMatrix(j,i) = curvatureMatrix(i,j);
      }
   }

   // Calculate Chi Square
   fChi2 = 0;
   for (i = 0; i < fSampleSize; i++) {
      Double_t f = 0;
      for (j = 0; j < fNCoefficients; j++)
         f += fCoefficients(j) * fFunctions(j,i);
      fChi2 += 1. / TMath::Max(fSqError(i),1e-20) * (fQuantity(i) - f)
      * (fQuantity(i) - f);
   }

   // Invert the curvature matrix
   const TVectorD diag = TMatrixDDiag_const(curvatureMatrix);
   curvatureMatrix.NormByDiag(diag);

   TDecompChol chol(curvatureMatrix);
   if (!chol.Decompose())
      Error("MakeCoefficientErrors", "curvature matrix is singular");
   chol.Invert(curvatureMatrix);

   curvatureMatrix.NormByDiag(diag);

   for (i = 0; i < fNCoefficients; i++)
      fCoefficientsRMS(i) = TMath::Sqrt(curvatureMatrix(i,i));
}


////////////////////////////////////////////////////////////////////////////////
/// PRIVATE METHOD:
/// Invert the model matrix B, and compute final coefficients. For a
/// more thorough discussion of what this means, please refer to the
/// <a href="#TMultiDimFit:description">class description</a>
///
/// First we invert the lower triangle matrix fOrthCurvatureMatrix
/// and store the inverted matrix in the upper triangle.

void TMultiDimFit::MakeCoefficients()
{
   Int_t i = 0, j = 0;
   Int_t col = 0, row = 0;

   // Invert the B matrix
   for (col = 1; col < fNCoefficients; col++) {
      for (row = col - 1; row > -1; row--) {
         fOrthCurvatureMatrix(row,col) = 0;
         for (i = row; i <= col ; i++)
            fOrthCurvatureMatrix(row,col) -=
            fOrthCurvatureMatrix(i,row)
            * fOrthCurvatureMatrix(i,col);
      }
   }

   // Compute the final coefficients
   fCoefficients.ResizeTo(fNCoefficients);

   for (i = 0; i < fNCoefficients; i++) {
      Double_t sum = 0;
      for (j = i; j < fNCoefficients; j++)
         sum += fOrthCurvatureMatrix(i,j) * fOrthCoefficients(j);
      fCoefficients(i) = sum;
   }

   // Compute the final residuals
   fResiduals.ResizeTo(fSampleSize);
   for (i = 0; i < fSampleSize; i++)
      fResiduals(i) = fQuantity(i);

   for (i = 0; i < fNCoefficients; i++)
      for (j = 0; j < fSampleSize; j++)
         fResiduals(j) -= fCoefficients(i) * fFunctions(i,j);

   // Compute the max and minimum, and squared sum of the evaluated
   // residuals
   fMinResidual = 10e10;
   fMaxResidual = -10e10;
   Double_t sqRes  = 0;
   for (i = 0; i < fSampleSize; i++){
      sqRes += fResiduals(i) * fResiduals(i);
      if (fResiduals(i) <= fMinResidual) {
         fMinResidual     = fResiduals(i);
         fMinResidualRow  = i;
      }
      if (fResiduals(i) >= fMaxResidual) {
         fMaxResidual     = fResiduals(i);
         fMaxResidualRow  = i;
      }
   }

   fCorrelationCoeff = fSumSqResidual / fSumSqAvgQuantity;
   fPrecision        = TMath::Sqrt(sqRes / fSumSqQuantity);

   // If we use histograms, fill some more
   if (TESTBIT(fHistogramMask,HIST_RD) ||
       TESTBIT(fHistogramMask,HIST_RTRAI) ||
       TESTBIT(fHistogramMask,HIST_RX)) {
      for (i = 0; i < fSampleSize; i++) {
         if (TESTBIT(fHistogramMask,HIST_RD))
            ((TH2D*)fHistograms->FindObject("res_d"))->Fill(fQuantity(i),
                                                            fResiduals(i));
         if (TESTBIT(fHistogramMask,HIST_RTRAI))
            ((TH1D*)fHistograms->FindObject("res_train"))->Fill(fResiduals(i));

         if (TESTBIT(fHistogramMask,HIST_RX))
            for (j = 0; j < fNVariables; j++)
               ((TH2D*)fHistograms->FindObject(Form("res_x_%d",j)))
               ->Fill(fVariables(i * fNVariables + j),fResiduals(i));
      }
   } // If histograms

}


////////////////////////////////////////////////////////////////////////////////
/// PRIVATE METHOD:
/// Compute the correlation matrix

void TMultiDimFit::MakeCorrelation()
{
   if (!fShowCorrelation)
      return;

   fCorrelationMatrix.ResizeTo(fNVariables,fNVariables+1);

   Double_t d2      = 0;
   Double_t ddotXi  = 0; // G.Q. needs to be reinitialized in the loop over i fNVariables
   Double_t xiNorm  = 0; // G.Q. needs to be reinitialized in the loop over i fNVariables
   Double_t xidotXj = 0; // G.Q. needs to be reinitialized in the loop over j fNVariables
   Double_t xjNorm  = 0; // G.Q. needs to be reinitialized in the loop over j fNVariables

   Int_t i, j, k, l, m;  // G.Q. added m variable
   for (i = 0; i < fSampleSize; i++)
      d2 += fQuantity(i) * fQuantity(i);

   for (i = 0; i < fNVariables; i++) {
      ddotXi = 0.; // G.Q. reinitialisation
      xiNorm = 0.; // G.Q. reinitialisation
      for (j = 0; j< fSampleSize; j++) {
         // Index of sample j of variable i
         k =  j * fNVariables + i;
         ddotXi += fQuantity(j) * (fVariables(k) - fMeanVariables(i));
         xiNorm += (fVariables(k) - fMeanVariables(i))
         * (fVariables(k) - fMeanVariables(i));
      }
      fCorrelationMatrix(i,0) = ddotXi / TMath::Sqrt(d2 * xiNorm);

      for (j = 0; j < i; j++) {
         xidotXj = 0.; // G.Q. reinitialisation
         xjNorm = 0.; // G.Q. reinitialisation
         for (k = 0; k < fSampleSize; k++) {
            // Index of sample j of variable i
            // l =  j * fNVariables + k;  // G.Q.
            l =  k * fNVariables + j; // G.Q.
            m =  k * fNVariables + i; // G.Q.
                                      // G.Q.        xidotXj += (fVariables(i) - fMeanVariables(i))
                                      // G.Q.          * (fVariables(l) - fMeanVariables(j));
            xidotXj += (fVariables(m) - fMeanVariables(i))
            * (fVariables(l) - fMeanVariables(j));  // G.Q. modified index for Xi
            xjNorm  += (fVariables(l) - fMeanVariables(j))
            * (fVariables(l) - fMeanVariables(j));
         }
         //fCorrelationMatrix(i+1,j) = xidotXj / TMath::Sqrt(xiNorm * xjNorm);
         fCorrelationMatrix(i,j+1) = xidotXj / TMath::Sqrt(xiNorm * xjNorm);
      }
   }
}



////////////////////////////////////////////////////////////////////////////////
/// PRIVATE METHOD:
/// Make Gram-Schmidt orthogonalisation. The class description gives
/// a thorough account of this algorithm, as well as
/// references. Please refer to the
/// <a href="#TMultiDimFit:description">class description</a>

Double_t TMultiDimFit::MakeGramSchmidt(Int_t function)
{

   // calculate w_i, that is, evaluate the current function at data
   // point i
   Double_t f2                        = 0;
   fOrthCoefficients(fNCoefficients)      = 0;
   fOrthFunctionNorms(fNCoefficients)  = 0;
   Int_t j        = 0;
   Int_t k        = 0;

   for (j = 0; j < fSampleSize; j++) {
      fFunctions(fNCoefficients, j) = 1;
      fOrthFunctions(fNCoefficients, j) = 0;
      // First, however, we need to calculate f_fNCoefficients
      for (k = 0; k < fNVariables; k++) {
         Int_t    p   =  fPowers[function * fNVariables + k];
         Double_t x   =  fVariables(j * fNVariables + k);
         fFunctions(fNCoefficients, j) *= EvalFactor(p,x);
      }

      // Calculate f dot f in f2
      f2 += fFunctions(fNCoefficients,j) *  fFunctions(fNCoefficients,j);
      // Assign to w_fNCoefficients f_fNCoefficients
      fOrthFunctions(fNCoefficients, j) = fFunctions(fNCoefficients, j);
   }

   // the first column of w is equal to f
   for (j = 0; j < fNCoefficients; j++) {
      Double_t fdw = 0;
      // Calculate (f_fNCoefficients dot w_j) / w_j^2
      for (k = 0; k < fSampleSize; k++) {
         fdw += fFunctions(fNCoefficients, k) * fOrthFunctions(j,k)
         / fOrthFunctionNorms(j);
      }

      fOrthCurvatureMatrix(fNCoefficients,j) = fdw;
      // and subtract it from the current value of w_ij
      for (k = 0; k < fSampleSize; k++)
         fOrthFunctions(fNCoefficients,k) -= fdw * fOrthFunctions(j,k);
   }

   for (j = 0; j < fSampleSize; j++) {
      // calculate squared length of w_fNCoefficients
      fOrthFunctionNorms(fNCoefficients) +=
      fOrthFunctions(fNCoefficients,j)
      * fOrthFunctions(fNCoefficients,j);

      // calculate D dot w_fNCoefficients in A
      fOrthCoefficients(fNCoefficients) += fQuantity(j)
      * fOrthFunctions(fNCoefficients, j);
   }

   // First test, but only if didn't user specify
   if (!fIsUserFunction)
      if (TMath::Sqrt(fOrthFunctionNorms(fNCoefficients) / (f2 + 1e-10))
          < TMath::Sin(fMinAngle*DEGRAD))
         return 0;

   // The result found by this code for the first residual is always
   // much less then the one found be MUDIFI. That's because it's
   // supposed to be. The cause is the improved precision of Double_t
   // over DOUBLE PRECISION!
   fOrthCurvatureMatrix(fNCoefficients,fNCoefficients) = 1;
   Double_t b = fOrthCoefficients(fNCoefficients);
   fOrthCoefficients(fNCoefficients) /= fOrthFunctionNorms(fNCoefficients);

   // Calculate the residual from including this fNCoefficients.
   Double_t dResidur = fOrthCoefficients(fNCoefficients) * b;

   return dResidur;
}


////////////////////////////////////////////////////////////////////////////////
/// Make histograms of the result of the analysis. This message
/// should be sent after having read all data points, but before
/// finding the parameterization
///
/// Options:
///     A         All the below
///     X         Original independent variables
///     D         Original dependent variables
///     N         Normalised independent variables
///     S         Shifted dependent variables
///     R1        Residuals versus normalised independent variables
///     R2        Residuals versus dependent variable
///     R3        Residuals computed on training sample
///     R4        Residuals computed on test sample
///
/// For a description of these quantities, refer to
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::MakeHistograms(Option_t *option)
{
   TString opt(option);
   opt.ToLower();

   if (opt.Length() < 1)
      return;

   if (!fHistograms)
      fHistograms = new TList;

   // Counter variable
   Int_t i = 0;

   // Histogram of original variables
   if (opt.Contains("x") || opt.Contains("a")) {
      SETBIT(fHistogramMask,HIST_XORIG);
      for (i = 0; i < fNVariables; i++)
         if (!fHistograms->FindObject(Form("x_%d_orig",i)))
            fHistograms->Add(new TH1D(Form("x_%d_orig",i),
                                      Form("Original variable # %d",i),
                                      fBinVarX, fMinVariables(i),
                                      fMaxVariables(i)));
   }

   // Histogram of original dependent variable
   if (opt.Contains("d") || opt.Contains("a")) {
      SETBIT(fHistogramMask,HIST_DORIG);
      if (!fHistograms->FindObject("d_orig"))
         fHistograms->Add(new TH1D("d_orig", "Original Quantity",
                                   fBinVarX, fMinQuantity, fMaxQuantity));
   }

   // Histograms of normalized variables
   if (opt.Contains("n") || opt.Contains("a")) {
      SETBIT(fHistogramMask,HIST_XNORM);
      for (i = 0; i < fNVariables; i++)
         if (!fHistograms->FindObject(Form("x_%d_norm",i)))
            fHistograms->Add(new TH1D(Form("x_%d_norm",i),
                                      Form("Normalized variable # %d",i),
                                      fBinVarX, -1,1));
   }

   // Histogram of shifted dependent variable
   if (opt.Contains("s") || opt.Contains("a")) {
      SETBIT(fHistogramMask,HIST_DSHIF);
      if (!fHistograms->FindObject("d_shifted"))
         fHistograms->Add(new TH1D("d_shifted", "Shifted Quantity",
                                   fBinVarX, fMinQuantity - fMeanQuantity,
                                   fMaxQuantity - fMeanQuantity));
   }

   // Residual from training sample versus independent variables
   if (opt.Contains("r1") || opt.Contains("a")) {
      SETBIT(fHistogramMask,HIST_RX);
      for (i = 0; i < fNVariables; i++)
         if (!fHistograms->FindObject(Form("res_x_%d",i)))
            fHistograms->Add(new TH2D(Form("res_x_%d",i),
                                      Form("Computed residual versus x_%d", i),
                                      fBinVarX, -1,    1,
                                      fBinVarY,
                                      fMinQuantity - fMeanQuantity,
                                      fMaxQuantity - fMeanQuantity));
   }

   // Residual from training sample versus. dependent variable
   if (opt.Contains("r2") || opt.Contains("a")) {
      SETBIT(fHistogramMask,HIST_RD);
      if (!fHistograms->FindObject("res_d"))
         fHistograms->Add(new TH2D("res_d",
                                   "Computed residuals vs Quantity",
                                   fBinVarX,
                                   fMinQuantity - fMeanQuantity,
                                   fMaxQuantity - fMeanQuantity,
                                   fBinVarY,
                                   fMinQuantity - fMeanQuantity,
                                   fMaxQuantity - fMeanQuantity));
   }

   // Residual from training sample
   if (opt.Contains("r3") || opt.Contains("a")) {
      SETBIT(fHistogramMask,HIST_RTRAI);
      if (!fHistograms->FindObject("res_train"))
         fHistograms->Add(new TH1D("res_train",
                                   "Computed residuals over training sample",
                                   fBinVarX, fMinQuantity - fMeanQuantity,
                                   fMaxQuantity - fMeanQuantity));

   }
   if (opt.Contains("r4") || opt.Contains("a")) {
      SETBIT(fHistogramMask,HIST_RTEST);
      if (!fHistograms->FindObject("res_test"))
         fHistograms->Add(new TH1D("res_test",
                                   "Distribution of residuals from test",
                                   fBinVarX,fMinQuantity - fMeanQuantity,
                                   fMaxQuantity - fMeanQuantity));
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Generate the file <classname>MDF.cxx which contains the
/// implementation of the method:
///
///   Double_t <classname>::MDF(Double_t *x)
///
/// which does the same as  TMultiDimFit::Eval. Please refer to this
/// method.
///
/// Further, the public static members:
///
///   Int_t    <classname>::fgNVariables
///   Int_t    <classname>::fgNCoefficients
///   Double_t <classname>::fgDMean
///   Double_t <classname>::fgXMean[]       //[fgNVariables]
///   Double_t <classname>::fgXMin[]        //[fgNVariables]
///   Double_t <classname>::fgXMax[]        //[fgNVariables]
///   Double_t <classname>::fgCoefficient[] //[fgNCoeffficents]
///   Int_t    <classname>::fgPower[]       //[fgNCoeffficents*fgNVariables]
///
/// are initialized, and assumed to exist. The class declaration is
/// assumed to be in <classname>.h and assumed to be provided by the
/// user.
///
/// See TMultiDimFit::MakeRealCode for a list of options
///
/// The minimal class definition is:
///
///   class <classname> {
///   public:
///     Int_t    <classname>::fgNVariables;     // Number of variables
///     Int_t    <classname>::fgNCoefficients;  // Number of terms
///     Double_t <classname>::fgDMean;          // Mean from training sample
///     Double_t <classname>::fgXMean[];        // Mean from training sample
///     Double_t <classname>::fgXMin[];         // Min from training sample
///     Double_t <classname>::fgXMax[];         // Max from training sample
///     Double_t <classname>::fgCoefficient[];  // Coefficients
///     Int_t    <classname>::fgPower[];        // Function powers
///
///     Double_t Eval(Double_t *x);
///   };
///
/// Whether the method <classname>::Eval should be static or not, is
/// up to the user.

void TMultiDimFit::MakeMethod(const Char_t* classname, Option_t* option)
{
   MakeRealCode(Form("%sMDF.cxx", classname), classname, option);
}



////////////////////////////////////////////////////////////////////////////////
/// PRIVATE METHOD:
/// Normalize data to the interval [-1;1]. This is needed for the
/// classes method to work.

void TMultiDimFit::MakeNormalized()
{
   Int_t i = 0;
   Int_t j = 0;
   Int_t k = 0;

   for (i = 0; i < fSampleSize; i++) {
      if (TESTBIT(fHistogramMask,HIST_DORIG))
         ((TH1D*)fHistograms->FindObject("d_orig"))->Fill(fQuantity(i));

      fQuantity(i) -= fMeanQuantity;
      fSumSqAvgQuantity  += fQuantity(i) * fQuantity(i);

      if (TESTBIT(fHistogramMask,HIST_DSHIF))
         ((TH1D*)fHistograms->FindObject("d_shifted"))->Fill(fQuantity(i));

      for (j = 0; j < fNVariables; j++) {
         Double_t range = 1. / (fMaxVariables(j) - fMinVariables(j));
         k              = i * fNVariables + j;

         // Fill histograms of original independent variables
         if (TESTBIT(fHistogramMask,HIST_XORIG))
            ((TH1D*)fHistograms->FindObject(Form("x_%d_orig",j)))
            ->Fill(fVariables(k));

         // Normalise independent variables
         fVariables(k) = 1 + 2 * range * (fVariables(k) - fMaxVariables(j));

         // Fill histograms of normalised independent variables
         if (TESTBIT(fHistogramMask,HIST_XNORM))
            ((TH1D*)fHistograms->FindObject(Form("x_%d_norm",j)))
            ->Fill(fVariables(k));

      }
   }
   // Shift min and max of dependent variable
   fMaxQuantity -= fMeanQuantity;
   fMinQuantity -= fMeanQuantity;

   // Shift mean of independent variables
   for (i = 0; i < fNVariables; i++) {
      Double_t range = 1. / (fMaxVariables(i) - fMinVariables(i));
      fMeanVariables(i) = 1 + 2 * range * (fMeanVariables(i)
                                           - fMaxVariables(i));
   }
}


////////////////////////////////////////////////////////////////////////////////
/// PRIVATE METHOD:
/// Find the parameterization over the training sample. A full account
/// of the algorithm is given in the
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::MakeParameterization()
{
   Int_t     i              = -1;
   Int_t     j              = 0;
   Int_t     k              = 0;
   Int_t     maxPass        = 3;
   Int_t     studied        = 0;
   Double_t  squareResidual = fSumSqAvgQuantity;
   fNCoefficients            = 0;
   fSumSqResidual           = fSumSqAvgQuantity;
   fFunctions.ResizeTo(fMaxTerms,fSampleSize);
   fOrthFunctions.ResizeTo(fMaxTerms,fSampleSize);
   fOrthFunctionNorms.ResizeTo(fMaxTerms);
   fOrthCoefficients.ResizeTo(fMaxTerms);
   fOrthCurvatureMatrix.ResizeTo(fMaxTerms,fMaxTerms);
   fFunctions = 1;

   fFunctionCodes = new Int_t[fMaxFunctions];
   fPowerIndex    = new Int_t[fMaxTerms];
   Int_t l;
   for (l=0;l<fMaxFunctions;l++) fFunctionCodes[l] = 0;
   for (l=0;l<fMaxTerms;l++)     fPowerIndex[l]    = 0;

   if (fMaxAngle != 0)  maxPass = 100;
   if (fIsUserFunction) maxPass = 1;

   // Loop over the number of functions we want to study.
   // increment inspection counter
   while(kTRUE) {

      // Reach user defined limit of studies
      if (studied++ >= fMaxStudy) {
         fParameterisationCode = PARAM_MAXSTUDY;
         break;
      }

      // Considered all functions several times
      if (k >= maxPass) {
         fParameterisationCode = PARAM_SEVERAL;
         break;
      }

      // increment function counter
      i++;

      // If we've reached the end of the functions, restart pass
      if (i == fMaxFunctions) {
         if (fMaxAngle != 0)
            fMaxAngle += (90 - fMaxAngle) / 2;
         i = 0;
         studied--;
         k++;
         continue;
      }
      if (studied == 1)
         fFunctionCodes[i] = 0;
      else if (fFunctionCodes[i] >= 2)
         continue;

      // Print a happy message
      if (fIsVerbose && studied == 1)
         std::cout << "Coeff   SumSqRes    Contrib   Angle      QM   Func"
         << "     Value        W^2  Powers" << std::endl;

      // Make the Gram-Schmidt
      Double_t dResidur = MakeGramSchmidt(i);

      if (dResidur == 0) {
         // This function is no good!
         // First test is in MakeGramSchmidt
         fFunctionCodes[i] = 1;
         continue;
      }

      // If user specified function, assume they know what they are doing
      if (!fIsUserFunction) {
         // Flag this function as considered
         fFunctionCodes[i] = 2;

         // Test if this function contributes to the fit
         if (!TestFunction(squareResidual, dResidur)) {
            fFunctionCodes[i] = 1;
            continue;
         }
      }

      // If we get to here, the function currently considered is
      // fNCoefficients, so we increment the counter
      // Flag this function as OK, and store and the number in the
      // index.
      fFunctionCodes[i]          = 3;
      fPowerIndex[fNCoefficients] = i;
      fNCoefficients++;

      // We add the current contribution to the sum of square of
      // residuals;
      squareResidual -= dResidur;


      // Calculate control parameter from this function
      for (j = 0; j < fNVariables; j++) {
         if (fNCoefficients == 1
             || fMaxPowersFinal[j] <= fPowers[i * fNVariables + j] - 1)
            fMaxPowersFinal[j] = fPowers[i * fNVariables + j] - 1;
      }
      Double_t s = EvalControl(&fPowers[i * fNVariables]);

      // Print the statistics about this function
      if (fIsVerbose) {
         std::cout << std::setw(5)  << fNCoefficients << " "
         << std::setw(10) << std::setprecision(4) << squareResidual << " "
         << std::setw(10) << std::setprecision(4) << dResidur << " "
         << std::setw(7)  << std::setprecision(3) << fMaxAngle << " "
         << std::setw(7)  << std::setprecision(3) << s << " "
         << std::setw(5)  << i << " "
         << std::setw(10) << std::setprecision(4)
         << fOrthCoefficients(fNCoefficients-1) << " "
         << std::setw(10) << std::setprecision(4)
         << fOrthFunctionNorms(fNCoefficients-1) << " "
         << std::flush;
         for (j = 0; j < fNVariables; j++)
            std::cout << " " << fPowers[i * fNVariables + j] - 1 << std::flush;
         std::cout << std::endl;
      }

      if (fNCoefficients >= fMaxTerms /* && fIsVerbose */) {
         fParameterisationCode = PARAM_MAXTERMS;
         break;
      }

      Double_t err  = TMath::Sqrt(TMath::Max(1e-20,squareResidual) /
                                  fSumSqAvgQuantity);
      if (err < fMinRelativeError) {
         fParameterisationCode = PARAM_RELERR;
         break;
      }

   }

   fError          = TMath::Max(1e-20,squareResidual);
   fSumSqResidual -= fError;
   fRMS = TMath::Sqrt(fError / fSampleSize);
}


////////////////////////////////////////////////////////////////////////////////
/// PRIVATE METHOD:
/// This is the method that actually generates the code for the
/// evaluation the parameterization on some point.
/// It's called by TMultiDimFit::MakeCode and TMultiDimFit::MakeMethod.
///
/// The options are: NONE so far

void TMultiDimFit::MakeRealCode(const char *filename,
                                const char *classname,
                                Option_t *)
{
   Int_t i, j;

   Bool_t  isMethod     = (classname[0] == '\0' ? kFALSE : kTRUE);
   const char *prefix   = (isMethod ? Form("%s::", classname) : "");
   const char *cv_qual  = (isMethod ? "" : "static ");

   std::ofstream outFile(filename,std::ios::out|std::ios::trunc);
   if (!outFile) {
      Error("MakeRealCode","couldn't open output file '%s'",filename);
      return;
   }

   if (fIsVerbose)
      std::cout << "Writing on file \"" << filename << "\" ... " << std::flush;
   //
   // Write header of file
   //
   // Emacs mode line ;-)
   outFile << "// -*- mode: c++ -*-" << std::endl;
   // Info about creator
   outFile << "// " << std::endl
   << "// File " << filename
   << " generated by TMultiDimFit::MakeRealCode" << std::endl;
   // Time stamp
   TDatime date;
   outFile << "// on " << date.AsString() << std::endl;
   // ROOT version info
   outFile << "// ROOT version " << gROOT->GetVersion()
   << std::endl << "//" << std::endl;
   // General information on the code
   outFile << "// This file contains the function " << std::endl
   << "//" << std::endl
   << "//    double  " << prefix << "MDF(double *x); " << std::endl
   << "//" << std::endl
   << "// For evaluating the parameterization obtained" << std::endl
   << "// from TMultiDimFit and the point x" << std::endl
   << "// " << std::endl
   << "// See TMultiDimFit class documentation for more "
   << "information " << std::endl << "// " << std::endl;
   // Header files
   if (isMethod)
      // If these are methods, we need the class header
      outFile << "#include \"" << classname << ".h\"" << std::endl;

   //
   // Now for the data
   //
   outFile << "//" << std::endl
   << "// Static data variables"  << std::endl
   << "//" << std::endl;
   outFile << cv_qual << "int    " << prefix << "gNVariables    = "
   << fNVariables << ";" << std::endl;
   outFile << cv_qual << "int    " << prefix << "gNCoefficients = "
   << fNCoefficients << ";" << std::endl;
   outFile << cv_qual << "double " << prefix << "gDMean         = "
   << fMeanQuantity << ";" << std::endl;

   // Assignment to mean vector.
   outFile << "// Assignment to mean vector." << std::endl;
   outFile << cv_qual << "double " << prefix
   << "gXMean[] = {" << std::endl;
   for (i = 0; i < fNVariables; i++)
      outFile << (i != 0 ? ", " : "  ") << fMeanVariables(i) << std::flush;
   outFile << " };" << std::endl << std::endl;

   // Assignment to minimum vector.
   outFile << "// Assignment to minimum vector." << std::endl;
   outFile << cv_qual << "double " << prefix
   << "gXMin[] = {" << std::endl;
   for (i = 0; i < fNVariables; i++)
      outFile << (i != 0 ? ", " : "  ") << fMinVariables(i) << std::flush;
   outFile << " };" << std::endl << std::endl;

   // Assignment to maximum vector.
   outFile << "// Assignment to maximum vector." << std::endl;
   outFile << cv_qual << "double " << prefix
   << "gXMax[] = {" << std::endl;
   for (i = 0; i < fNVariables; i++)
      outFile << (i != 0 ? ", " : "  ") << fMaxVariables(i) << std::flush;
   outFile << " };" << std::endl << std::endl;

   // Assignment to coefficients vector.
   outFile << "// Assignment to coefficients vector." << std::endl;
   outFile << cv_qual << "double " << prefix
   << "gCoefficient[] = {" << std::flush;
   for (i = 0; i < fNCoefficients; i++)
      outFile << (i != 0 ? "," : "") << std::endl
      << "  " << fCoefficients(i) << std::flush;
   outFile << std::endl << " };" << std::endl << std::endl;

   // Assignment to error coefficients vector.
   outFile << "// Assignment to error coefficients vector." << std::endl;
   outFile << cv_qual << "double " << prefix
   << "gCoefficientRMS[] = {" << std::flush;
   for (i = 0; i < fNCoefficients; i++)
      outFile << (i != 0 ? "," : "") << std::endl
      << "  " << fCoefficientsRMS(i) << std::flush;
   outFile << std::endl << " };" << std::endl << std::endl;

   // Assignment to powers vector.
   outFile << "// Assignment to powers vector." << std::endl
   << "// The powers are stored row-wise, that is" << std::endl
   << "//  p_ij = " << prefix
   << "gPower[i * NVariables + j];" << std::endl;
   outFile << cv_qual << "int    " << prefix
   << "gPower[] = {" << std::flush;
   for (i = 0; i < fNCoefficients; i++) {
      for (j = 0; j < fNVariables; j++) {
         if (j != 0) outFile << std::flush << "  ";
         else        outFile << std::endl << "  ";
         outFile << fPowers[fPowerIndex[i] * fNVariables + j]
         << (i == fNCoefficients - 1 && j == fNVariables - 1 ? "" : ",")
         << std::flush;
      }
   }
   outFile << std::endl << "};" << std::endl << std::endl;


   //
   // Finally we reach the function itself
   //
   outFile << "// " << std::endl
   << "// The "
   << (isMethod ? "method " : "function ")
   << "  double " << prefix
   << "MDF(double *x)"
   << std::endl << "// " << std::endl;
   outFile << "double " << prefix
   << "MDF(double *x) {" << std::endl
   << "  double returnValue = " << prefix << "gDMean;" << std::endl
   << "  int    i = 0, j = 0, k = 0;" << std::endl
   << "  for (i = 0; i < " << prefix << "gNCoefficients ; i++) {"
   << std::endl
   << "    // Evaluate the ith term in the expansion" << std::endl
   << "    double term = " << prefix << "gCoefficient[i];"
   << std::endl
   << "    for (j = 0; j < " << prefix << "gNVariables; j++) {"
   << std::endl
   << "      // Evaluate the polynomial in the jth variable." << std::endl
   << "      int power = "<< prefix << "gPower["
   << prefix << "gNVariables * i + j]; " << std::endl
   << "      double p1 = 1, p2 = 0, p3 = 0, r = 0;" << std::endl
   << "      double v =  1 + 2. / ("
   << prefix << "gXMax[j] - " << prefix
   << "gXMin[j]) * (x[j] - " << prefix << "gXMax[j]);" << std::endl
   << "      // what is the power to use!" << std::endl
   << "      switch(power) {" << std::endl
   << "      case 1: r = 1; break; " << std::endl
   << "      case 2: r = v; break; " << std::endl
   << "      default: " << std::endl
   << "        p2 = v; " << std::endl
   << "        for (k = 3; k <= power; k++) { " << std::endl
   << "          p3 = p2 * v;" << std::endl;
   if (fPolyType == kLegendre)
      outFile << "          p3 = ((2 * i - 3) * p2 * v - (i - 2) * p1)"
      << " / (i - 1);" << std::endl;
   if (fPolyType == kChebyshev)
      outFile << "          p3 = 2 * v * p2 - p1; " << std::endl;
   outFile << "          p1 = p2; p2 = p3; " << std::endl << "        }" << std::endl
   << "        r = p3;" << std::endl << "      }" << std::endl
   << "      // multiply this term by the poly in the jth var" << std::endl
   << "      term *= r; " << std::endl << "    }" << std::endl
   << "    // Add this term to the final result" << std::endl
   << "    returnValue += term;" << std::endl << "  }" << std::endl
   << "  return returnValue;" << std::endl << "}" << std::endl << std::endl;

   // EOF
   outFile << "// EOF for " << filename << std::endl;

   // Close the file
   outFile.close();

   if (fIsVerbose)
      std::cout << "done" << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Print statistics etc.
/// Options are
///   P        Parameters
///   S        Statistics
///   C        Coefficients
///   R        Result of parameterisation
///   F        Result of fit
///   K        Correlation Matrix
///   M        Pretty print formula
///

void TMultiDimFit::Print(Option_t *option) const
{
   Int_t i = 0;
   Int_t j = 0;

   TString opt(option);
   opt.ToLower();

   if (opt.Contains("p")) {
      // Print basic parameters for this object
      std::cout << "User parameters:" << std::endl
      << "----------------" << std::endl
      << " Variables:                    " << fNVariables << std::endl
      << " Data points:                  " << fSampleSize << std::endl
      << " Max Terms:                    " << fMaxTerms << std::endl
      << " Power Limit Parameter:        " << fPowerLimit << std::endl
      << " Max functions:                " << fMaxFunctions << std::endl
      << " Max functions to study:       " << fMaxStudy << std::endl
      << " Max angle (optional):         " << fMaxAngle << std::endl
      << " Min angle:                    " << fMinAngle << std::endl
      << " Relative Error accepted:      " << fMinRelativeError << std::endl
      << " Maximum Powers:               " << std::flush;
      for (i = 0; i < fNVariables; i++)
         std::cout << " " << fMaxPowers[i] - 1 << std::flush;
      std::cout << std::endl << std::endl
      << " Parameterisation will be done using " << std::flush;
      if (fPolyType == kChebyshev)
         std::cout << "Chebyshev polynomials" << std::endl;
      else if (fPolyType == kLegendre)
         std::cout << "Legendre polynomials" << std::endl;
      else
         std::cout << "Monomials" << std::endl;
      std::cout << std::endl;
   }

   if (opt.Contains("s")) {
      // Print statistics for read data
      std::cout << "Sample statistics:" << std::endl
      << "------------------" << std::endl
      << "                 D"  << std::flush;
      for (i = 0; i < fNVariables; i++)
         std::cout << " " << std::setw(10) << i+1 << std::flush;
      std::cout << std::endl << " Max:   " << std::setw(10) << std::setprecision(7)
      << fMaxQuantity << std::flush;
      for (i = 0; i < fNVariables; i++)
         std::cout << " " << std::setw(10) << std::setprecision(4)
         << fMaxVariables(i) << std::flush;
      std::cout << std::endl << " Min:   " << std::setw(10) << std::setprecision(7)
      << fMinQuantity << std::flush;
      for (i = 0; i < fNVariables; i++)
         std::cout << " " << std::setw(10) << std::setprecision(4)
         << fMinVariables(i) << std::flush;
      std::cout << std::endl << " Mean:  " << std::setw(10) << std::setprecision(7)
      << fMeanQuantity << std::flush;
      for (i = 0; i < fNVariables; i++)
         std::cout << " " << std::setw(10) << std::setprecision(4)
         << fMeanVariables(i) << std::flush;
      std::cout << std::endl << " Function Sum Squares:         " << fSumSqQuantity
      << std::endl << std::endl;
   }

   if (opt.Contains("r")) {
      std::cout << "Results of Parameterisation:" << std::endl
      << "----------------------------" << std::endl
      << " Total reduction of square residuals    "
      << fSumSqResidual << std::endl
      << " Relative precision obtained:           "
      << fPrecision   << std::endl
      << " Error obtained:                        "
      << fError << std::endl
      << " Multiple correlation coefficient:      "
      << fCorrelationCoeff   << std::endl
      << " Reduced Chi square over sample:        "
      << fChi2 / (fSampleSize - fNCoefficients) << std::endl
      << " Maximum residual value:                "
      << fMaxResidual << std::endl
      << " Minimum residual value:                "
      << fMinResidual << std::endl
      << " Estimated root mean square:            "
      << fRMS << std::endl
      << " Maximum powers used:                   " << std::flush;
      for (j = 0; j < fNVariables; j++)
         std::cout << fMaxPowersFinal[j] << " " << std::flush;
      std::cout << std::endl
      << " Function codes of candidate functions." << std::endl
      << "  1: considered,"
      << "  2: too little contribution,"
      << "  3: accepted." << std::flush;
      for (i = 0; i < fMaxFunctions; i++) {
         if (i % 60 == 0)
            std::cout << std::endl << " " << std::flush;
         else if (i % 10 == 0)
            std::cout << " " << std::flush;
         std::cout << fFunctionCodes[i];
      }
      std::cout << std::endl << " Loop over candidates stopped because " << std::flush;
      switch(fParameterisationCode){
         case PARAM_MAXSTUDY:
            std::cout << "max allowed studies reached" << std::endl; break;
         case PARAM_SEVERAL:
            std::cout << "all candidates considered several times" << std::endl; break;
         case PARAM_RELERR:
            std::cout << "wanted relative error obtained" << std::endl; break;
         case PARAM_MAXTERMS:
            std::cout << "max number of terms reached" << std::endl; break;
         default:
            std::cout << "some unknown reason" << std::endl;
            break;
      }
      std::cout << std::endl;
   }

   if (opt.Contains("f")) {
      std::cout << "Results of Fit:" << std::endl
      << "---------------" << std::endl
      << " Test sample size:                      "
      << fTestSampleSize << std::endl
      << " Multiple correlation coefficient:      "
      << fTestCorrelationCoeff << std::endl
      << " Relative precision obtained:           "
      << fTestPrecision   << std::endl
      << " Error obtained:                        "
      << fTestError << std::endl
      << " Reduced Chi square over sample:        "
      << fChi2 / (fSampleSize - fNCoefficients) << std::endl
      << std::endl;
      if (fFitter) {
         fFitter->PrintResults(1,1);
         std::cout << std::endl;
      }
   }

   if (opt.Contains("c")){
      std::cout << "Coefficients:" << std::endl
      << "-------------" << std::endl
      << "   #         Value        Error   Powers" << std::endl
      << " ---------------------------------------" << std::endl;
      for (i = 0; i < fNCoefficients; i++) {
         std::cout << " " << std::setw(3) << i << "  "
         << std::setw(12) << fCoefficients(i) << "  "
         << std::setw(12) << fCoefficientsRMS(i) << "  " << std::flush;
         for (j = 0; j < fNVariables; j++)
            std::cout << " " << std::setw(3)
            << fPowers[fPowerIndex[i] * fNVariables + j] - 1 << std::flush;
         std::cout << std::endl;
      }
      std::cout << std::endl;
   }
   if (opt.Contains("k") && fCorrelationMatrix.IsValid()) {
      std::cout << "Correlation Matrix:" << std::endl
      << "-------------------";
      fCorrelationMatrix.Print();
   }

   if (opt.Contains("m")) {
      std::cout << "Parameterization:" << std::endl
      << "-----------------" << std::endl
      << "  Normalised variables: " << std::endl;
      for (i = 0; i < fNVariables; i++)
         std::cout << "\ty_" << i << "\t= 1 + 2 * (x_" << i << " - "
         << fMaxVariables(i) << ") / ("
         << fMaxVariables(i) << " - " << fMinVariables(i) << ")"
         << std::endl;
      std::cout << std::endl
      << "  f(";
      for (i = 0; i < fNVariables; i++) {
         std::cout << "y_" << i;
         if (i != fNVariables-1) std::cout << ", ";
      }
      std::cout << ") = ";
      for (i = 0; i < fNCoefficients; i++) {
         if (i != 0)
            std::cout << std::endl << "\t" << (fCoefficients(i) < 0 ? "- " : "+ ")
            << TMath::Abs(fCoefficients(i));
         else
            std::cout << fCoefficients(i);
         for (j = 0; j < fNVariables; j++) {
            Int_t p = fPowers[fPowerIndex[i] * fNVariables + j];
            switch (p) {
               case 1: break;
               case 2: std::cout << " * y_" << j; break;
               default:
                  switch(fPolyType) {
                     case kLegendre:  std::cout << " * L_" << p-1 << "(y_" << j << ")"; break;
                     case kChebyshev: std::cout << " * C_" << p-1 << "(y_" << j << ")"; break;
                     default:         std::cout << " * y_" << j << "^" << p-1; break;
                  }
            }

         }
      }
      std::cout << std::endl;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Selection method. User can override this method for specialized
/// selection of acceptable functions in fit. Default is to select
/// all. This message is sent during the build-up of the function
/// candidates table once for each set of powers in
/// variables. Notice, that the argument array contains the powers
/// PLUS ONE. For example, to De select the function
///     f = x1^2 * x2^4 * x3^5,
/// this method should return kFALSE if given the argument
///     { 3, 4, 6 }

Bool_t TMultiDimFit::Select(const Int_t *)
{
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the max angle (in degrees) between the initial data vector to
/// be fitted, and the new candidate function to be included in the
/// fit.  By default it is 0, which automatically chooses another
/// selection criteria. See also
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::SetMaxAngle(Double_t ang)
{
   if (ang >= 90 || ang < 0) {
      Warning("SetMaxAngle", "angle must be in [0,90)");
      return;
   }

   fMaxAngle = ang;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the min angle (in degrees) between a new candidate function
/// and the subspace spanned by the previously accepted
/// functions. See also
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::SetMinAngle(Double_t ang)
{
   if (ang > 90 || ang <= 0) {
      Warning("SetMinAngle", "angle must be in [0,90)");
      return;
   }

   fMinAngle = ang;

}


////////////////////////////////////////////////////////////////////////////////
/// Define a user function. The input array must be of the form
/// (p11, ..., p1N, ... ,pL1, ..., pLN)
/// Where N is the dimension of the data sample, L is the number of
/// terms (given in terms) and the first number, labels the term, the
/// second the variable.  More information is given in the
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::SetPowers(const Int_t* powers, Int_t terms)
{
   fIsUserFunction = kTRUE;
   fMaxFunctions   = terms;
   fMaxTerms       = terms;
   fMaxStudy       = terms;
   fMaxFuncNV      = fMaxFunctions * fNVariables;
   fPowers         = new Int_t[fMaxFuncNV];
   Int_t i, j;
   for (i = 0; i < fMaxFunctions; i++)
      for(j = 0; j < fNVariables; j++)
         fPowers[i * fNVariables + j] = powers[i * fNVariables + j]  + 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the user parameter for the function selection. The bigger the
/// limit, the more functions are used. The meaning of this variable
/// is defined in the
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::SetPowerLimit(Double_t limit)
{
   fPowerLimit = limit;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the maximum power to be considered in the fit for each
/// variable. See also
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::SetMaxPowers(const Int_t* powers)
{
   if (!powers)
      return;

   for (Int_t i = 0; i < fNVariables; i++)
      fMaxPowers[i] = powers[i]+1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the acceptable relative error for when sum of square
/// residuals is considered minimized. For a full account, refer to
/// the
/// <a href="#TMultiDimFit:description">class description</a>

void TMultiDimFit::SetMinRelativeError(Double_t error)
{
   fMinRelativeError = error;
}


////////////////////////////////////////////////////////////////////////////////
/// PRIVATE METHOD:
/// Test whether the currently considered function contributes to the
/// fit. See also
/// <a href="#TMultiDimFit:description">class description</a>

Bool_t TMultiDimFit::TestFunction(Double_t squareResidual,
                                  Double_t dResidur)
{
   if (fNCoefficients != 0) {
      // Now for the second test:
      if (fMaxAngle == 0) {
         // If the user hasn't supplied a max angle do the test as,
         if (dResidur <
             squareResidual / (fMaxTerms - fNCoefficients + 1 + 1E-10)) {
            return kFALSE;
         }
      }
      else {
         // If the user has provided a max angle, test if the calculated
         // angle is less then the max angle.
         if (TMath::Sqrt(dResidur/fSumSqAvgQuantity) <
             TMath::Cos(fMaxAngle*DEGRAD)) {
            return kFALSE;
         }
      }
   }
   // If we get here, the function is OK
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Helper function for doing the minimisation of Chi2 using Minuit

void mdfHelper(int& /*npar*/, double* /*divs*/, double& chi2,
               double* coeffs, int /*flag*/)
{
   // Get pointer  to current TMultiDimFit object.
   TMultiDimFit* mdf = TMultiDimFit::Instance();
   chi2     = mdf->MakeChi2(coeffs);
}
