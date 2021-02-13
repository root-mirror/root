//////////////////////////////////////////////////////////////////////////
//
// 'BASIC FUNCTIONALITY' RooFit tutorial macro #101
//
// Fitting, plotting, toy data generation on one-dimensional p.d.f
//
// pdf = gauss(x,m,s)
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooFormulaVar.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooUnitTest.h"

using namespace RooFit ;


// Elementary operations on a gaussian PDF
class TestBasic101 : public RooUnitTest
{
public:
  TestBasic101(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Fitting,plotting & event generation of basic p.d.f",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // S e t u p   m o d e l
    // ---------------------

    // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
    RooRealVar x("x","x",-10,10) ;
    RooRealVar mean("mean","mean of gaussian",1,-10,10) ;
    RooRealVar sigma("sigma","width of gaussian",1,0.1,10) ;

    // Build gaussian p.d.f in terms of x,mean and sigma
    RooGaussian gauss("gauss","gaussian PDF",x,mean,sigma) ;

    // Construct plot frame in 'x'
    RooPlot* xframe = x.frame(Title("Gaussian p.d.f.")) ;


    // P l o t   m o d e l   a n d   c h a n g e   p a r a m e t e r   v a l u e s
    // ---------------------------------------------------------------------------

    // Plot gauss in frame (i.e. in x)
    gauss.plotOn(xframe) ;

    // Change the value of sigma to 3
    sigma.setVal(3) ;

    // Plot gauss in frame (i.e. in x) and draw frame on canvas
    gauss.plotOn(xframe,LineColor(kRed),Name("another")) ;


    // G e n e r a t e   e v e n t s
    // -----------------------------

    // Generate a dataset of 1000 events in x from gauss
    RooDataSet* data = gauss.generate(x,10000) ;

    // Make a second plot frame in x and draw both the
    // data and the p.d.f in the frame
    RooPlot* xframe2 = x.frame(Title("Gaussian p.d.f. with data")) ;
    data->plotOn(xframe2) ;
    gauss.plotOn(xframe2) ;


    // F i t   m o d e l   t o   d a t a
    // -----------------------------

    // Fit pdf to data
    gauss.fitTo(*data) ;


    // --- Post processing for stressRooFit ---
    regPlot(xframe ,"rf101_plot1") ;
    regPlot(xframe2,"rf101_plot2") ;

    delete data ;

    return kTRUE ;
  }
} ;


//////////////////////////////////////////////////////////////////////////
//
// 'BASIC FUNCTIONALITY' RooFit tutorial macro #102
//
// Importing data from ROOT TTrees and THx histograms
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit ;


class TestBasic102 : public RooUnitTest
{
public:
  TestBasic102(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Data import methods",refFile,writeRef,verbose) {} ;

  TH1* makeTH1()
  {
    // Create ROOT TH1 filled with a Gaussian distribution

    TH1D* hh = new TH1D("hh","hh",25,-10,10) ;
    for (int i=0 ; i<100 ; i++) {
      hh->Fill(gRandom->Gaus(0,3)) ;
    }
    return hh ;
  }


  TTree* makeTTree()
  {
    // Create ROOT TTree filled with a Gaussian distribution in x and a uniform distribution in y

    TTree* tree = new TTree("tree","tree") ;
    Double_t* px = new Double_t ;
    Double_t* py = new Double_t ;
    tree->Branch("x",px,"x/D") ;
    tree->Branch("y",py,"y/D") ;
    for (int i=0 ; i<100 ; i++) {
      *px = gRandom->Gaus(0,3) ;
      *py = gRandom->Uniform()*30 - 15 ;
      tree->Fill() ;
    }

    //delete px ;
    //delete py ;

    return tree ;
  }

  Bool_t testCode() {

    ////////////////////////////////////////////////////////
    // I m p o r t i n g   R O O T   h i s t o g r a m s  //
    ////////////////////////////////////////////////////////

    // I m p o r t   T H 1   i n t o   a   R o o D a t a H i s t
    // ---------------------------------------------------------

    // Create a ROOT TH1 histogram
    TH1* hh = makeTH1() ;

    // Declare observable x
    RooRealVar x("x","x",-10,10) ;

    // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
    RooDataHist dh("dh","dh",x,Import(*hh)) ;


    // P l o t   a n d   f i t   a   R o o D a t a H i s t
    // ---------------------------------------------------

    // Make plot of binned dataset showing Poisson error bars (RooFit default)
    RooPlot* frame = x.frame(Title("Imported TH1 with Poisson error bars")) ;
    dh.plotOn(frame) ;

    // Fit a Gaussian p.d.f to the data
    RooRealVar mean("mean","mean",0,-10,10) ;
    RooRealVar sigma("sigma","sigma",3,0.1,10) ;
    RooGaussian gauss("gauss","gauss",x,mean,sigma) ;
    gauss.fitTo(dh) ;
    gauss.plotOn(frame) ;

    // P l o t   a n d   f i t   a   R o o D a t a H i s t   w i t h   i n t e r n a l   e r r o r s
    // ---------------------------------------------------------------------------------------------

    // If histogram has custom error (i.e. its contents is does not originate from a Poisson process
    // but e.g. is a sum of weighted events) you can data with symmetric 'sum-of-weights' error instead
    // (same error bars as shown by ROOT)
    RooPlot* frame2 = x.frame(Title("Imported TH1 with internal errors")) ;
    dh.plotOn(frame2,DataError(RooAbsData::SumW2)) ;
    gauss.plotOn(frame2) ;

    // Please note that error bars shown (Poisson or SumW2) are for visualization only, the are NOT used
    // in a maximum likelihood fit
    //
    // A (binned) ML fit will ALWAYS assume the Poisson error interpretation of data (the mathematical definition
    // of likelihood does not take any external definition of errors). Data with non-unit weights can only be correctly
    // fitted with a chi^2 fit (see rf602_chi2fit.C)


    ////////////////////////////////////////////////
    // I m p o r t i n g   R O O T  T T r e e s   //
    ////////////////////////////////////////////////


    // I m p o r t   T T r e e   i n t o   a   R o o D a t a S e t
    // -----------------------------------------------------------

    TTree* tree = makeTTree() ;

    // Define 2nd observable y
    RooRealVar y("y","y",-10,10) ;

    // Construct unbinned dataset importing tree branches x and y matching between branches and RooRealVars
    // is done by name of the branch/RRV
    //
    // Note that ONLY entries for which x,y have values within their allowed ranges as defined in
    // RooRealVar x and y are imported. Since the y values in the import tree are in the range [-15,15]
    // and RRV y defines a range [-10,10] this means that the RooDataSet below will have less entries than the TTree 'tree'

    RooDataSet ds("ds","ds",RooArgSet(x,y),Import(*tree)) ;


    // P l o t   d a t a s e t   w i t h   m u l t i p l e   b i n n i n g   c h o i c e s
    // ------------------------------------------------------------------------------------

    // Print unbinned dataset with default frame binning (100 bins)
    RooPlot* frame3 = y.frame(Title("Unbinned data shown in default frame binning")) ;
    ds.plotOn(frame3) ;

    // Print unbinned dataset with custom binning choice (20 bins)
    RooPlot* frame4 = y.frame(Title("Unbinned data shown with custom binning")) ;
    ds.plotOn(frame4,Binning(20)) ;

    // Draw all frames on a canvas
    regPlot(frame ,"rf102_plot1") ;
    regPlot(frame2,"rf102_plot2") ;
    regPlot(frame3,"rf102_plot3") ;
    regPlot(frame4,"rf102_plot4") ;

    delete hh ;
    delete tree ;

    return kTRUE ;
  }
} ;







//////////////////////////////////////////////////////////////////////////
//
// 'BASIC FUNCTIONALITY' RooFit tutorial macro #103
//
// Interpreted functions and p.d.f.s
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooConstVar.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"

using namespace RooFit ;


class TestBasic103 : public RooUnitTest
{
public:
  TestBasic103(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Interpreted expression p.d.f.",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    /////////////////////////////////////////////////////////
    // G e n e r i c   i n t e r p r e t e d   p . d . f . //
    /////////////////////////////////////////////////////////

    // Declare observable x
    RooRealVar x("x","x",-20,20) ;

    // C o n s t r u c t   g e n e r i c   p d f   f r o m   i n t e r p r e t e d   e x p r e s s i o n
    // -------------------------------------------------------------------------------------------------

    // To construct a proper p.d.f, the formula expression is explicitly normalized internally by dividing
    // it by a numeric integral of the expresssion over x in the range [-20,20]
    //
    RooRealVar alpha("alpha","alpha",5,0.1,10) ;
    RooGenericPdf genpdf("genpdf","genpdf","(1+0.1*abs(x)+sin(sqrt(abs(x*alpha+0.1))))",RooArgSet(x,alpha)) ;


    // S a m p l e ,   f i t   a n d   p l o t   g e n e r i c   p d f
    // ---------------------------------------------------------------

    // Generate a toy dataset from the interpreted p.d.f
    RooDataSet* data = genpdf.generate(x,10000) ;

    // Fit the interpreted p.d.f to the generated data
    genpdf.fitTo(*data) ;

    // Make a plot of the data and the p.d.f overlaid
    RooPlot* xframe = x.frame(Title("Interpreted expression pdf")) ;
    data->plotOn(xframe) ;
    genpdf.plotOn(xframe) ;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // S t a n d a r d   p . d . f   a d j u s t   w i t h   i n t e r p r e t e d   h e l p e r   f u n c t i o n //
    //                                                                                                             //
    // Make a gauss(x,sqrt(mean2),sigma) from a standard RooGaussian                                               //
    //                                                                                                             //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // C o n s t r u c t   s t a n d a r d   p d f  w i t h   f o r m u l a   r e p l a c i n g   p a r a m e t e r
    // ------------------------------------------------------------------------------------------------------------

    // Construct parameter mean2 and sigma
    RooRealVar mean2("mean2","mean^2",10,0,200) ;
    RooRealVar sigma("sigma","sigma",3,0.1,10) ;

    // Construct interpreted function mean = sqrt(mean^2)
    RooFormulaVar mean("mean","mean","sqrt(mean2)",mean2) ;

    // Construct a gaussian g2(x,sqrt(mean2),sigma) ;
    RooGaussian g2("g2","h2",x,mean,sigma) ;


    // G e n e r a t e   t o y   d a t a
    // ---------------------------------

    // Construct a separate gaussian g1(x,10,3) to generate a toy Gaussian dataset with mean 10 and width 3
    RooGaussian g1("g1","g1",x,RooConst(10),RooConst(3)) ;
    RooDataSet* data2 = g1.generate(x,1000) ;


    // F i t   a n d   p l o t   t a i l o r e d   s t a n d a r d   p d f
    // -------------------------------------------------------------------

    // Fit g2 to data from g1
    RooFitResult* r = g2.fitTo(*data2,Save()) ;

    // Plot data on frame and overlay projection of g2
    RooPlot* xframe2 = x.frame(Title("Tailored Gaussian pdf")) ;
    data2->plotOn(xframe2) ;
    g2.plotOn(xframe2) ;

    regPlot(xframe,"rf103_plot1") ;
    regPlot(xframe2,"rf103_plot2") ;
    regResult(r,"rf103_fit1") ;

    delete data ;
    delete data2 ;

    return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'BASIC FUNCTIONALITY' RooFit tutorial macro #105
//
//  Demonstration of binding ROOT Math functions as RooFit functions
//  and pdfs
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TMath.h"
#include "TF1.h"
#include "Math/DistFunc.h"
#include "RooCFunction1Binding.h"
#include "RooCFunction3Binding.h"
#include "RooTFnBinding.h"

using namespace RooFit ;

class TestBasic105 : public RooUnitTest
{
public:
  TestBasic105(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("C++ function binding operator p.d.f",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // B i n d   T M a t h : : E r f   C   f u n c t i o n
    // ---------------------------------------------------

    // Bind one-dimensional TMath::Erf function as RooAbsReal function
    RooRealVar x("x","x",-3,3) ;
    RooAbsReal* erf = bindFunction("erf",TMath::Erf,x) ;

    // Plot erf on frame
    RooPlot* frame1 = x.frame(Title("TMath::Erf bound as RooFit function")) ;
    erf->plotOn(frame1) ;


    // B i n d   R O O T : : M a t h : : b e t a _ p d f   C   f u n c t i o n
    // -----------------------------------------------------------------------

    // Bind pdf ROOT::Math::Beta with three variables as RooAbsPdf function
    // exclude x=0 and x=1 points from range since beta_pdf diverges at x=0 and x=1
    // when a < 1 and/or b <1
    RooRealVar x2("x2","x2",0.001,0.999) ;
    RooRealVar a("a","a",5,0,10) ;
    RooRealVar b("b","b",2,0,10) ;
    RooAbsPdf* beta = bindPdf("beta",ROOT::Math::beta_pdf,x2,a,b) ;

    // Generate some events and fit
    RooDataSet* data = beta->generate(x2,10000) ;
    beta->fitTo(*data) ;

    // Plot data and pdf on frame
    RooPlot* frame2 = x2.frame(Title("ROOT::Math::Beta bound as RooFit pdf")) ;
    data->plotOn(frame2) ;
    beta->plotOn(frame2) ;



    // B i n d   R O O T   T F 1   a s   R o o F i t   f u n c t i o n
    // ---------------------------------------------------------------

    // Create a ROOT TF1 function
    TF1 *fa1 = new TF1("fa1","sin(x)/x",0,10);

    // Create an observable
    RooRealVar x3("x3","x3",0.01,20) ;

    // Create binding of TF1 object to above observable
    RooAbsReal* rfa1 = bindFunction(fa1,x3) ;

    // Make plot frame in observable, plot TF1 binding function
    RooPlot* frame3 = x3.frame(Title("TF1 bound as RooFit function")) ;
    rfa1->plotOn(frame3) ;


    regPlot(frame1,"rf105_plot1") ;
    regPlot(frame2,"rf105_plot2") ;
    regPlot(frame3,"rf105_plot3") ;

    delete erf ;
    delete beta ;
    delete fa1 ;
    delete rfa1 ;
    delete data ;

    return kTRUE ;
  }
} ;
/////////////////////////////////////////////////////////////////////////
//
// 'BASIC FUNCTIONALITY' RooFit tutorial macro #108
//
// Plotting unbinned data with alternate and variable binnings
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooBMixDecay.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
using namespace RooFit ;


class TestBasic108 : public RooUnitTest
{
public:
  TestBasic108(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Non-standard binning in counting and asymmetry plots",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // S e t u p   m o d e l
    // ---------------------

    // Build a B decay p.d.f with mixing
    RooRealVar dt("dt","dt",-20,20) ;
    RooRealVar dm("dm","dm",0.472) ;
    RooRealVar tau("tau","tau",1.547) ;
    RooRealVar w("w","mistag rate",0.1) ;
    RooRealVar dw("dw","delta mistag rate",0.) ;

    RooCategory mixState("mixState","B0/B0bar mixing state") ;
    mixState.defineType("mixed",-1) ;
    mixState.defineType("unmixed",1) ;
    RooCategory tagFlav("tagFlav","Flavour of the tagged B0") ;
    tagFlav.defineType("B0",1) ;
    tagFlav.defineType("B0bar",-1) ;

    // Build a gaussian resolution model
    RooRealVar dterr("dterr","dterr",0.1,1.0) ;
    RooRealVar bias1("bias1","bias1",0) ;
    RooRealVar sigma1("sigma1","sigma1",0.1) ;
    RooGaussModel gm1("gm1","gauss model 1",dt,bias1,sigma1) ;

    // Construct Bdecay (x) gauss
    RooBMixDecay bmix("bmix","decay",dt,mixState,tagFlav,tau,dm,w,dw,gm1,RooBMixDecay::DoubleSided) ;


    // S a m p l e   d a t a   f r o m   m o d e l
    // --------------------------------------------

    // Sample 2000 events in (dt,mixState,tagFlav) from bmix
    RooDataSet *data = bmix.generate(RooArgSet(dt,mixState,tagFlav),2000) ;



    // S h o w   d t   d i s t r i b u t i o n   w i t h   c u s t o m   b i n n i n g
    // -------------------------------------------------------------------------------

    // Make plot of dt distribution of data in range (-15,15) with fine binning for dt>0 and coarse binning for dt<0

    // Create binning object with range (-15,15)
    RooBinning tbins(-15,15) ;

    // Add 60 bins with uniform spacing in range (-15,0)
    tbins.addUniform(60,-15,0) ;

    // Add 15 bins with uniform spacing in range (0,15)
    tbins.addUniform(15,0,15) ;

    // Make plot with specified binning
    RooPlot* dtframe = dt.frame(Range(-15,15),Title("dt distribution with custom binning")) ;
    data->plotOn(dtframe,Binning(tbins)) ;
    bmix.plotOn(dtframe) ;

    // NB: Note that bin density for each bin is adjusted to that of default frame binning as shown
    // in Y axis label (100 bins --> Events/0.4*Xaxis-dim) so that all bins represent a consistent density distribution


    // S h o w   m i x s t a t e   a s y m m e t r y  w i t h   c u s t o m   b i n n i n g
    // ------------------------------------------------------------------------------------

    // Make plot of dt distribution of data asymmetry in 'mixState' with variable binning

    // Create binning object with range (-10,10)
    RooBinning abins(-10,10) ;

    // Add boundaries at 0, (-1,1), (-2,2), (-3,3), (-4,4) and (-6,6)
    abins.addBoundary(0) ;
    abins.addBoundaryPair(1) ;
    abins.addBoundaryPair(2) ;
    abins.addBoundaryPair(3) ;
    abins.addBoundaryPair(4) ;
    abins.addBoundaryPair(6) ;

    // Create plot frame in dt
    RooPlot* aframe = dt.frame(Range(-10,10),Title("mixState asymmetry distribution with custom binning")) ;

    // Plot mixState asymmetry of data with specified customg binning
    data->plotOn(aframe,Asymmetry(mixState),Binning(abins)) ;

    // Plot corresponding property of p.d.f
    bmix.plotOn(aframe,Asymmetry(mixState)) ;

    // Adjust vertical range of plot to sensible values for an asymmetry
    aframe->SetMinimum(-1.1) ;
    aframe->SetMaximum(1.1) ;

    // NB: For asymmetry distributions no density corrects are needed (and are thus not applied)


    regPlot(dtframe,"rf108_plot1") ;
    regPlot(aframe,"rf108_plot2") ;

    delete data ;

    return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'BASIC FUNCTIONALITY' RooFit tutorial macro #109
//
// Calculating chi^2 from histograms and curves in RooPlots,
// making histogram of residual and pull distributions
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHist.h"
using namespace RooFit ;

class TestBasic109 : public RooUnitTest
{
public:
  TestBasic109(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Calculation of chi^2 and residuals in plots",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // S e t u p   m o d e l
    // ---------------------

    // Create observables
    RooRealVar x("x","x",-10,10) ;

    // Create Gaussian
    RooRealVar sigma("sigma","sigma",3,0.1,10) ;
    RooRealVar mean("mean","mean",0,-10,10) ;
    RooGaussian gauss("gauss","gauss",x,RooConst(0),sigma) ;

    // Generate a sample of 1000 events with sigma=3
    RooDataSet* data = gauss.generate(x,10000) ;

    // Change sigma to 3.15
    sigma=3.15 ;


    // P l o t   d a t a   a n d   s l i g h t l y   d i s t o r t e d   m o d e l
    // ---------------------------------------------------------------------------

    // Overlay projection of gauss with sigma=3.15 on data with sigma=3.0
    RooPlot* frame1 = x.frame(Title("Data with distorted Gaussian pdf"),Bins(40)) ;
    data->plotOn(frame1,DataError(RooAbsData::SumW2)) ;
    gauss.plotOn(frame1) ;


    // C a l c u l a t e   c h i ^ 2
    // ------------------------------

    // Show the chi^2 of the curve w.r.t. the histogram
    // If multiple curves or datasets live in the frame you can specify
    // the name of the relevant curve and/or dataset in chiSquare()
    regValue(frame1->chiSquare(),"rf109_chi2") ;


    // S h o w   r e s i d u a l   a n d   p u l l   d i s t s
    // -------------------------------------------------------

    // Construct a histogram with the residuals of the data w.r.t. the curve
    RooHist* hresid = frame1->residHist() ;

    // Construct a histogram with the pulls of the data w.r.t the curve
    RooHist* hpull = frame1->pullHist() ;

    // Create a new frame to draw the residual distribution and add the distribution to the frame
    RooPlot* frame2 = x.frame(Title("Residual Distribution")) ;
    frame2->addPlotable(hresid,"P") ;

    // Create a new frame to draw the pull distribution and add the distribution to the frame
    RooPlot* frame3 = x.frame(Title("Pull Distribution")) ;
    frame3->addPlotable(hpull,"P") ;

    regPlot(frame1,"rf109_plot1") ;
    regPlot(frame2,"rf109_plot2") ;
    regPlot(frame3,"rf109_plot3") ;

    delete data ;
    //delete hresid ;
    //delete hpull ;

    return kTRUE ;
  }
} ;
/////////////////////////////////////////////////////////////////////////
//
// 'BASIC FUNCTIONALITY' RooFit tutorial macro #110
//
// Examples on normalization of p.d.f.s,
// integration of p.d.fs, construction
// of cumulative distribution functions from p.d.f.s
// in one dimension
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
#include "TCanvas.h"
using namespace RooFit ;

class TestBasic110 : public RooUnitTest
{
public:
  TestBasic110(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Normalization of p.d.f.s in 1D",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // S e t u p   m o d e l
    // ---------------------

    // Create observables x,y
    RooRealVar x("x","x",-10,10) ;

    // Create p.d.f. gaussx(x,-2,3)
    RooGaussian gx("gx","gx",x,RooConst(-2),RooConst(3)) ;


    // R e t r i e v e   r a w  &   n o r m a l i z e d   v a l u e s   o f   R o o F i t   p . d . f . s
    // --------------------------------------------------------------------------------------------------

    // Return 'raw' unnormalized value of gx
    regValue(gx.getVal(),"rf110_gx") ;

    // Return value of gx normalized over x in range [-10,10]
    RooArgSet nset(x) ;

    regValue(gx.getVal(&nset),"rf110_gx_Norm[x]") ;

    // Create object representing integral over gx
    // which is used to calculate  gx_Norm[x] == gx / gx_Int[x]
    RooAbsReal* igx = gx.createIntegral(x) ;
    regValue(igx->getVal(),"rf110_gx_Int[x]") ;


    // I n t e g r a t e   n o r m a l i z e d   p d f   o v e r   s u b r a n g e
    // ----------------------------------------------------------------------------

    // Define a range named "signal" in x from -5,5
    x.setRange("signal",-5,5) ;

    // Create an integral of gx_Norm[x] over x in range "signal"
    // This is the fraction of of p.d.f. gx_Norm[x] which is in the
    // range named "signal"
    RooAbsReal* igx_sig = gx.createIntegral(x,NormSet(x),Range("signal")) ;
    regValue(igx_sig->getVal(),"rf110_gx_Int[x|signal]_Norm[x]") ;



    // C o n s t r u c t   c u m u l a t i v e   d i s t r i b u t i o n   f u n c t i o n   f r o m   p d f
    // -----------------------------------------------------------------------------------------------------

    // Create the cumulative distribution function of gx
    // i.e. calculate Int[-10,x] gx(x') dx'
    RooAbsReal* gx_cdf = gx.createCdf(x) ;

    // Plot cdf of gx versus x
    RooPlot* frame = x.frame(Title("c.d.f of Gaussian p.d.f")) ;
    gx_cdf->plotOn(frame) ;


    regPlot(frame,"rf110_plot1") ;

    delete igx ;
    delete igx_sig ;
    delete gx_cdf ;

    return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'BASIC FUNCTIONALITY' RooFit tutorial macro #111
//
// Configuration and customization of how numeric (partial) integrals
// are executed
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooLandau.h"
#include "RooArgSet.h"
#include <iomanip>
using namespace RooFit ;

class TestBasic111 : public RooUnitTest
{
public:
  TestBasic111(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Numeric integration configuration",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // A d j u s t   g l o b a l   1 D   i n t e g r a t i o n   p r e c i s i o n
    // ----------------------------------------------------------------------------

    // Example: Change global precision for 1D integrals from 1e-7 to 1e-6
    //
    // The relative epsilon (change as fraction of current best integral estimate) and
    // absolute epsilon (absolute change w.r.t last best integral estimate) can be specified
    // separately. For most p.d.f integrals the relative change criterium is the most important,
    // however for certain non-p.d.f functions that integrate out to zero a separate absolute
    // change criterium is necessary to declare convergence of the integral
    //
    // NB: This change is for illustration only. In general the precision should be at least 1e-7
    // for normalization integrals for MINUIT to succeed.
    //
    RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-6) ;
    RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-6) ;


    // N u m e r i c   i n t e g r a t i o n   o f   l a n d a u   p d f
    // ------------------------------------------------------------------

    // Construct p.d.f without support for analytical integrator for demonstration purposes
    RooRealVar x("x","x",-10,10) ;
    RooLandau landau("landau","landau",x,RooConst(0),RooConst(0.1)) ;


    // Calculate integral over landau with default choice of numeric integrator
    RooAbsReal* intLandau = landau.createIntegral(x) ;
    Double_t val = intLandau->getVal() ;
    regValue(val,"rf111_val1") ;


    // S a m e   w i t h   c u s t o m   c o n f i g u r a t i o n
    // -----------------------------------------------------------


    // Construct a custom configuration which uses the adaptive Gauss-Kronrod technique
    // for closed 1D integrals
    RooNumIntConfig customConfig(*RooAbsReal::defaultIntegratorConfig()) ;
#ifdef R__HAS_MATHMORE
    customConfig.method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D") ;
#endif

    // Calculate integral over landau with custom integral specification
    RooAbsReal* intLandau2 = landau.createIntegral(x,NumIntConfig(customConfig)) ;
    Double_t val2 = intLandau2->getVal() ;
    regValue(val2,"rf111_val2") ;


    // A d j u s t i n g   d e f a u l t   c o n f i g   f o r   a   s p e c i f i c   p d f
    // -------------------------------------------------------------------------------------


    // Another possibility: associate custom numeric integration configuration as default for object 'landau'
    landau.setIntegratorConfig(customConfig) ;


    // Calculate integral over landau custom numeric integrator specified as object default
    RooAbsReal* intLandau3 = landau.createIntegral(x) ;
    Double_t val3 = intLandau3->getVal() ;
    regValue(val3,"rf111_val3") ;


    delete intLandau ;
    delete intLandau2 ;
    delete intLandau3 ;

    return kTRUE ;
  }
} ;
/////////////////////////////////////////////////////////////////////////
//
// 'ADDITION AND CONVOLUTION' RooFit tutorial macro #201
//
// Composite p.d.f with signal and background component
//
// pdf = f_bkg * bkg(x,a0,a1) + (1-fbkg) * (f_sig1 * sig1(x,m,s1 + (1-f_sig1) * sig2(x,m,s2)))
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic201 : public RooUnitTest
{
public:
  TestBasic201(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Addition operator p.d.f.",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // S e t u p   c o m p o n e n t   p d f s
    // ---------------------------------------

    // Declare observable x
    RooRealVar x("x","x",0,10) ;

    // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
    RooRealVar mean("mean","mean of gaussians",5) ;
    RooRealVar sigma1("sigma1","width of gaussians",0.5) ;
    RooRealVar sigma2("sigma2","width of gaussians",1) ;

    RooGaussian sig1("sig1","Signal component 1",x,mean,sigma1) ;
    RooGaussian sig2("sig2","Signal component 2",x,mean,sigma2) ;

    // Build Chebychev polynomial p.d.f.
    RooRealVar a0("a0","a0",0.5,0.,1.) ;
    RooRealVar a1("a1","a1",-0.2,-1.,1.) ;
    RooChebychev bkg("bkg","Background",x,RooArgSet(a0,a1)) ;


    ////////////////////////////////////////////////////
    // M E T H O D   1 - T w o   R o o A d d P d f s  //
    ////////////////////////////////////////////////////


    // A d d   s i g n a l   c o m p o n e n t s
    // ------------------------------------------

    // Sum the signal components into a composite signal p.d.f.
    RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
    RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),sig1frac) ;


    // A d d  s i g n a l   a n d   b a c k g r o u n d
    // ------------------------------------------------

    // Sum the composite signal and background
    RooRealVar bkgfrac("bkgfrac","fraction of background",0.5,0.,1.) ;
    RooAddPdf  model("model","g1+g2+a",RooArgList(bkg,sig),bkgfrac) ;


    // S a m p l e ,   f i t   a n d   p l o t   m o d e l
    // ---------------------------------------------------

    // Generate a data sample of 1000 events in x from model
    RooDataSet *data = model.generate(x,1000) ;

    // Fit model to data
    model.fitTo(*data) ;

    // Plot data and PDF overlaid
    RooPlot* xframe = x.frame(Title("Example of composite pdf=(sig1+sig2)+bkg")) ;
    data->plotOn(xframe) ;
    model.plotOn(xframe) ;

    // Overlay the background component of model with a dashed line
    model.plotOn(xframe,Components(bkg),LineStyle(kDashed)) ;

    // Overlay the background+sig2 components of model with a dotted line
    model.plotOn(xframe,Components(RooArgSet(bkg,sig2)),LineStyle(kDotted)) ;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // M E T H O D   2 - O n e   R o o A d d P d f   w i t h   r e c u r s i v e   f r a c t i o n s  //
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Construct sum of models on one go using recursive fraction interpretations
    //
    //   model2 = bkg + (sig1 + sig2)
    //
    RooAddPdf  model2("model","g1+g2+a",RooArgList(bkg,sig1,sig2),RooArgList(bkgfrac,sig1frac),kTRUE) ;

    // NB: Each coefficient is interpreted as the fraction of the
    // left-hand component of the i-th recursive sum, i.e.
    //
    //   sum4 = A + ( B + ( C + D)  with fraction fA, fB and fC expands to
    //
    //   sum4 = fA*A + (1-fA)*(fB*B + (1-fB)*(fC*C + (1-fC)*D))


    // P l o t   r e c u r s i v e   a d d i t i o n   m o d e l
    // ---------------------------------------------------------
    model2.plotOn(xframe,LineColor(kRed),LineStyle(kDashed)) ;
    model2.plotOn(xframe,Components(RooArgSet(bkg,sig2)),LineColor(kRed),LineStyle(kDashed)) ;


    regPlot(xframe,"rf201_plot1") ;

    delete data ;
    return kTRUE ;

  }

} ;

//////////////////////////////////////////////////////////////////////////
//
// 'ADDITION AND CONVOLUTION' RooFit tutorial macro #202
//
// Setting up an extended maximum likelihood fit
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic202 : public RooUnitTest
{
public:
  TestBasic202(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Extended ML fits to addition operator p.d.f.s",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // S e t u p   c o m p o n e n t   p d f s
    // ---------------------------------------

    // Declare observable x
    RooRealVar x("x","x",0,10) ;

    // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
    RooRealVar mean("mean","mean of gaussians",5) ;
    RooRealVar sigma1("sigma1","width of gaussians",0.5) ;
    RooRealVar sigma2("sigma2","width of gaussians",1) ;

    RooGaussian sig1("sig1","Signal component 1",x,mean,sigma1) ;
    RooGaussian sig2("sig2","Signal component 2",x,mean,sigma2) ;

    // Build Chebychev polynomial p.d.f.
    RooRealVar a0("a0","a0",0.5,0.,1.) ;
    RooRealVar a1("a1","a1",-0.2,-1.,1.) ;
    RooChebychev bkg("bkg","Background",x,RooArgSet(a0,a1)) ;

    // Sum the signal components into a composite signal p.d.f.
    RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
    RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),sig1frac) ;

    /////////////////////
    // M E T H O D   1 //
    /////////////////////


    // C o n s t r u c t   e x t e n d e d   c o m p o s i t e   m o d e l
    // -------------------------------------------------------------------

    // Sum the composite signal and background into an extended pdf nsig*sig+nbkg*bkg
    RooRealVar nsig("nsig","number of signal events",500,0.,10000) ;
    RooRealVar nbkg("nbkg","number of background events",500,0,10000) ;
    RooAddPdf  model("model","(g1+g2)+a",RooArgList(bkg,sig),RooArgList(nbkg,nsig)) ;



    // S a m p l e ,   f i t   a n d   p l o t   e x t e n d e d   m o d e l
    // ---------------------------------------------------------------------

    // Generate a data sample of expected number events in x from model
    // = model.expectedEvents() = nsig+nbkg
    RooDataSet *data = model.generate(x) ;

    // Fit model to data, extended ML term automatically included
    model.fitTo(*data) ;

    // Plot data and PDF overlaid, use expected number of events for p.d.f projection normalization
    // rather than observed number of events (==data->numEntries())
    RooPlot* xframe = x.frame(Title("extended ML fit example")) ;
    data->plotOn(xframe) ;
    model.plotOn(xframe,Normalization(1.0,RooAbsReal::RelativeExpected)) ;

    // Overlay the background component of model with a dashed line
    model.plotOn(xframe,Components(bkg),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected)) ;

    // Overlay the background+sig2 components of model with a dotted line
    model.plotOn(xframe,Components(RooArgSet(bkg,sig2)),LineStyle(kDotted),Normalization(1.0,RooAbsReal::RelativeExpected)) ;


    /////////////////////
    // M E T H O D   2 //
    /////////////////////

    // C o n s t r u c t   e x t e n d e d   c o m p o n e n t s   f i r s t
    // ---------------------------------------------------------------------

    // Associated nsig/nbkg as expected number of events with sig/bkg
    RooExtendPdf esig("esig","extended signal p.d.f",sig,nsig) ;
    RooExtendPdf ebkg("ebkg","extended background p.d.f",bkg,nbkg) ;


    // S u m   e x t e n d e d   c o m p o n e n t s   w i t h o u t   c o e f s
    // -------------------------------------------------------------------------

    // Construct sum of two extended p.d.f. (no coefficients required)
    RooAddPdf  model2("model2","(g1+g2)+a",RooArgList(ebkg,esig)) ;


    regPlot(xframe,"rf202_plot1") ;

    delete data ;
    return kTRUE ;

  }

} ;
/////////////////////////////////////////////////////////////////////////
//
// 'ADDITION AND CONVOLUTION' RooFit tutorial macro #203
//
// Fitting and plotting in sub ranges
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
using namespace RooFit ;


class TestBasic203 : public RooUnitTest
{
public:
  TestBasic203(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Basic fitting and plotting in ranges",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // S e t u p   m o d e l
    // ---------------------

    // Construct observables x
    RooRealVar x("x","x",-10,10) ;

    // Construct gaussx(x,mx,1)
    RooRealVar mx("mx","mx",0,-10,10) ;
    RooGaussian gx("gx","gx",x,mx,RooConst(1)) ;

    // Construct px = 1 (flat in x)
    RooPolynomial px("px","px",x) ;

    // Construct model = f*gx + (1-f)px
    RooRealVar f("f","f",0.,1.) ;
    RooAddPdf model("model","model",RooArgList(gx,px),f) ;

    // Generated 10000 events in (x,y) from p.d.f. model
    RooDataSet* modelData = model.generate(x,10000) ;

    // F i t   f u l l   r a n g e
    // ---------------------------

    // Fit p.d.f to all data
    RooFitResult* r_full = model.fitTo(*modelData,Save(kTRUE)) ;


    // F i t   p a r t i a l   r a n g e
    // ----------------------------------

    // Define "signal" range in x as [-3,3]
    x.setRange("signal",-3,3) ;

    // Fit p.d.f only to data in "signal" range
    RooFitResult* r_sig = model.fitTo(*modelData,Save(kTRUE),Range("signal")) ;


    // P l o t   /   p r i n t   r e s u l t s
    // ---------------------------------------

    // Make plot frame in x and add data and fitted model
    RooPlot* frame = x.frame(Title("Fitting a sub range")) ;
    modelData->plotOn(frame) ;
    model.plotOn(frame,LineStyle(kDashed),LineColor(kRed)) ; // Add shape in full ranged dashed
    model.plotOn(frame) ; // By default only fitted range is shown

    regPlot(frame,"rf203_plot") ;
    regResult(r_full,"rf203_r_full") ;
    regResult(r_sig,"rf203_r_sig") ;

    delete modelData ;
    return kTRUE;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'ADDITION AND CONVOLUTION' RooFit tutorial macro #204
//
// Extended maximum likelihood fit with alternate range definition
// for observed number of events.
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic204 : public RooUnitTest
{
public:
  TestBasic204(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Extended ML fit in sub range",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // S e t u p   c o m p o n e n t   p d f s
    // ---------------------------------------

    // Declare observable x
    RooRealVar x("x","x",0,10) ;

    // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
    RooRealVar mean("mean","mean of gaussians",5) ;
    RooRealVar sigma1("sigma1","width of gaussians",0.5) ;
    RooRealVar sigma2("sigma2","width of gaussians",1) ;

    RooGaussian sig1("sig1","Signal component 1",x,mean,sigma1) ;
    RooGaussian sig2("sig2","Signal component 2",x,mean,sigma2) ;

    // Build Chebychev polynomial p.d.f.
    RooRealVar a0("a0","a0",0.5,0.,1.) ;
    RooRealVar a1("a1","a1",-0.2,-1.,1.) ;
    RooChebychev bkg("bkg","Background",x,RooArgSet(a0,a1)) ;

    // Sum the signal components into a composite signal p.d.f.
    RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
    RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),sig1frac) ;


    // C o n s t r u c t   e x t e n d e d   c o m p s   wi t h   r a n g e   s p e c
    // ------------------------------------------------------------------------------

    // Define signal range in which events counts are to be defined
    x.setRange("signalRange",4,6) ;

    // Associated nsig/nbkg as expected number of events with sig/bkg _in_the_range_ "signalRange"
    RooRealVar nsig("nsig","number of signal events in signalRange",500,0.,10000) ;
    RooRealVar nbkg("nbkg","number of background events in signalRange",500,0,10000) ;
    RooExtendPdf esig("esig","extended signal p.d.f",sig,nsig,"signalRange") ;
    RooExtendPdf ebkg("ebkg","extended background p.d.f",bkg,nbkg,"signalRange") ;


    // S u m   e x t e n d e d   c o m p o n e n t s
    // ---------------------------------------------

    // Construct sum of two extended p.d.f. (no coefficients required)
    RooAddPdf  model("model","(g1+g2)+a",RooArgList(ebkg,esig)) ;


    // S a m p l e   d a t a ,   f i t   m o d e l
    // -------------------------------------------

    // Generate 1000 events from model so that nsig,nbkg come out to numbers <<500 in fit
    RooDataSet *data = model.generate(x,1000) ;


    // Perform unbinned extended ML fit to data
    RooFitResult* r = model.fitTo(*data,Extended(kTRUE),Save()) ;


    regResult(r,"rf204_result") ;

    delete data ;
    return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'ADDITION AND CONVOLUTION' RooFit tutorial macro #205
//
// Options for plotting components of composite p.d.f.s.
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic205 : public RooUnitTest
{
public:
  TestBasic205(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Component plotting variations",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // S e t u p   c o m p o s i t e    p d f
    // --------------------------------------

    // Declare observable x
    RooRealVar x("x","x",0,10) ;

    // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
    RooRealVar mean("mean","mean of gaussians",5) ;
    RooRealVar sigma1("sigma1","width of gaussians",0.5) ;
    RooRealVar sigma2("sigma2","width of gaussians",1) ;
    RooGaussian sig1("sig1","Signal component 1",x,mean,sigma1) ;
    RooGaussian sig2("sig2","Signal component 2",x,mean,sigma2) ;

    // Sum the signal components into a composite signal p.d.f.
    RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
    RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),sig1frac) ;

    // Build Chebychev polynomial p.d.f.
    RooRealVar a0("a0","a0",0.5,0.,1.) ;
    RooRealVar a1("a1","a1",-0.2,-1.,1.) ;
    RooChebychev bkg1("bkg1","Background 1",x,RooArgSet(a0,a1)) ;

    // Build expontential pdf
    RooRealVar alpha("alpha","alpha",-1) ;
    RooExponential bkg2("bkg2","Background 2",x,alpha) ;

    // Sum the background components into a composite background p.d.f.
    RooRealVar bkg1frac("sig1frac","fraction of component 1 in background",0.2,0.,1.) ;
    RooAddPdf bkg("bkg","Signal",RooArgList(bkg1,bkg2),sig1frac) ;

    // Sum the composite signal and background
    RooRealVar bkgfrac("bkgfrac","fraction of background",0.5,0.,1.) ;
    RooAddPdf  model("model","g1+g2+a",RooArgList(bkg,sig),bkgfrac) ;



    // S e t u p   b a s i c   p l o t   w i t h   d a t a   a n d   f u l l   p d f
    // ------------------------------------------------------------------------------

    // Generate a data sample of 1000 events in x from model
    RooDataSet *data = model.generate(x,1000) ;

    // Plot data and complete PDF overlaid
    RooPlot* xframe  = x.frame(Title("Component plotting of pdf=(sig1+sig2)+(bkg1+bkg2)")) ;
    data->plotOn(xframe) ;
    model.plotOn(xframe) ;

    // Clone xframe for use below
    RooPlot* xframe2 = x.frame(Title("Component plotting of pdf=(sig1+sig2)+(bkg1+bkg2)"),Name("xframe2")) ;


    // M a k e   c o m p o n e n t   b y   o b j e c t   r e f e r e n c e
    // --------------------------------------------------------------------

    // Plot single background component specified by object reference
    model.plotOn(xframe,Components(bkg),LineColor(kRed)) ;

    // Plot single background component specified by object reference
    model.plotOn(xframe,Components(bkg2),LineStyle(kDashed),LineColor(kRed)) ;

    // Plot multiple background components specified by object reference
    // Note that specified components may occur at any level in object tree
    // (e.g bkg is component of 'model' and 'sig2' is component 'sig')
    model.plotOn(xframe,Components(RooArgSet(bkg,sig2)),LineStyle(kDotted)) ;



    // M a k e   c o m p o n e n t   b y   n a m e  /   r e g e x p
    // ------------------------------------------------------------

    // Plot single background component specified by name
    model.plotOn(xframe2,Components("bkg"),LineColor(kCyan)) ;

    // Plot multiple background components specified by name
    model.plotOn(xframe2,Components("bkg1,sig2"),LineStyle(kDotted),LineColor(kCyan)) ;

    // Plot multiple background components specified by regular expression on name
    model.plotOn(xframe2,Components("sig*"),LineStyle(kDashed),LineColor(kCyan)) ;

    // Plot multiple background components specified by multiple regular expressions on name
    model.plotOn(xframe2,Components("bkg1,sig*"),LineStyle(kDashed),LineColor(kYellow),Invisible()) ;


    regPlot(xframe,"rf205_plot1") ;
    regPlot(xframe2,"rf205_plot2") ;

    delete data ;
    return kTRUE ;

  }

} ;
/////////////////////////////////////////////////////////////////////////
//
// 'ADDITION AND CONVOLUTION' RooFit tutorial macro #208
//
// One-dimensional numeric convolution
// (require ROOT to be compiled with --enable-fftw3)
//
// pdf = landau(t) (x) gauss(t)
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TPluginManager.h"
#include "TROOT.h"

using namespace RooFit ;



class TestBasic208 : public RooUnitTest
{
public:
  TestBasic208(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("FFT Convolution operator p.d.f.",refFile,writeRef,verbose) {} ;

  Bool_t isTestAvailable() {
     // only if ROOT was build with fftw3 enabled
     TString conffeatures = gROOT->GetConfigFeatures();
     if(conffeatures.Contains("fftw3")) {
        TPluginHandler *h;
        if ((h = gROOT->GetPluginManager()->FindHandler("TVirtualFFT"))) {
           if (h->LoadPlugin() == -1) {
              gROOT->ProcessLine("new TNamed ;") ;
              return kFALSE;
           } else {
              return kTRUE ;
           }
        }
     }
     return kFALSE ;
  }

  Double_t ctol() { return 1e-2 ; } // Account for difficult shape of Landau distribution

  Bool_t testCode() {

    // S e t u p   c o m p o n e n t   p d f s
    // ---------------------------------------

    // Construct observable
    RooRealVar t("t","t",-10,30) ;

    // Construct landau(t,ml,sl) ;
    RooRealVar ml("ml","mean landau",5.,-20,20) ;
    RooRealVar sl("sl","sigma landau",1,0.1,10) ;
    RooLandau landau("lx","lx",t,ml,sl) ;

    // Construct gauss(t,mg,sg)
    RooRealVar mg("mg","mg",0) ;
    RooRealVar sg("sg","sg",2,0.1,10) ;
    RooGaussian gauss("gauss","gauss",t,mg,sg) ;


    // C o n s t r u c t   c o n v o l u t i o n   p d f
    // ---------------------------------------

    // Set #bins to be used for FFT sampling to 10000
    t.setBins(10000,"cache") ;

    // Construct landau (x) gauss
    RooFFTConvPdf lxg("lxg","landau (X) gauss",t,landau,gauss) ;


    // S a m p l e ,   f i t   a n d   p l o t   c o n v o l u t e d   p d f
    // ----------------------------------------------------------------------

    // Sample 1000 events in x from gxlx
    RooDataSet* data = lxg.generate(t,10000) ;

    // Fit gxlx to data
    lxg.fitTo(*data) ;

    // Plot data, landau pdf, landau (X) gauss pdf
    RooPlot* frame = t.frame(Title("landau (x) gauss convolution")) ;
    data->plotOn(frame) ;
    lxg.plotOn(frame) ;
    landau.plotOn(frame,LineStyle(kDashed)) ;

    regPlot(frame,"rf208_plot1") ;

    delete data ;
    return kTRUE ;

  }
} ;



/////////////////////////////////////////////////////////////////////////
//
// 'ADDITION AND CONVOLUTION' RooFit tutorial macro #209
//
// Decay function p.d.fs with optional B physics
// effects (mixing and CP violation) that can be
// analytically convolved with e.g. Gaussian resolution
// functions
//
// pdf1 = decay(t,tau) (x) delta(t)
// pdf2 = decay(t,tau) (x) gauss(t,m,s)
// pdf3 = decay(t,tau) (x) (f*gauss1(t,m1,s1) + (1-f)*gauss2(t,m1,s1))
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussModel.h"
#include "RooAddModel.h"
#include "RooTruthModel.h"
#include "RooDecay.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
using namespace RooFit ;



class TestBasic209 : public RooUnitTest
{
public:
  TestBasic209(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Analytical convolution operator",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // B - p h y s i c s   p d f   w i t h   t r u t h   r e s o l u t i o n
    // ---------------------------------------------------------------------

    // Variables of decay p.d.f.
    RooRealVar dt("dt","dt",-10,10) ;
    RooRealVar tau("tau","tau",1.548) ;

    // Build a truth resolution model (delta function)
    RooTruthModel tm("tm","truth model",dt) ;

    // Construct decay(t) (x) delta(t)
    RooDecay decay_tm("decay_tm","decay",dt,tau,tm,RooDecay::DoubleSided) ;

    // Plot p.d.f. (dashed)
    RooPlot* frame = dt.frame(Title("Bdecay (x) resolution")) ;
    decay_tm.plotOn(frame,LineStyle(kDashed)) ;


    // B - p h y s i c s   p d f   w i t h   G a u s s i a n   r e s o l u t i o n
    // ----------------------------------------------------------------------------

    // Build a gaussian resolution model
    RooRealVar bias1("bias1","bias1",0) ;
    RooRealVar sigma1("sigma1","sigma1",1) ;
    RooGaussModel gm1("gm1","gauss model 1",dt,bias1,sigma1) ;

    // Construct decay(t) (x) gauss1(t)
    RooDecay decay_gm1("decay_gm1","decay",dt,tau,gm1,RooDecay::DoubleSided) ;

    // Plot p.d.f.
    decay_gm1.plotOn(frame) ;


    // B - p h y s i c s   p d f   w i t h   d o u b l e   G a u s s i a n   r e s o l u t i o n
    // ------------------------------------------------------------------------------------------

    // Build another gaussian resolution model
    RooRealVar bias2("bias2","bias2",0) ;
    RooRealVar sigma2("sigma2","sigma2",5) ;
    RooGaussModel gm2("gm2","gauss model 2",dt,bias2,sigma2) ;

    // Build a composite resolution model f*gm1+(1-f)*gm2
    RooRealVar gm1frac("gm1frac","fraction of gm1",0.5) ;
    RooAddModel gmsum("gmsum","sum of gm1 and gm2",RooArgList(gm1,gm2),gm1frac) ;

    // Construct decay(t) (x) (f*gm1 + (1-f)*gm2)
    RooDecay decay_gmsum("decay_gmsum","decay",dt,tau,gmsum,RooDecay::DoubleSided) ;

    // Plot p.d.f. (red)
    decay_gmsum.plotOn(frame,LineColor(kRed)) ;

    regPlot(frame,"rf209_plot1") ;

    return kTRUE ;

  }
} ;
/////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #301
//
// Multi-dimensional p.d.f.s through composition, e.g. substituting a
// p.d.f parameter with a function that depends on other observables
//
// pdf = gauss(x,f(y),s) with f(y) = a0 + a1*y
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolyVar.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
using namespace RooFit ;



class TestBasic301 : public RooUnitTest
{
public:
  TestBasic301(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Composition extension of basic p.d.f",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // S e t u p   c o m p o s e d   m o d e l   g a u s s ( x , m ( y ) , s )
  // -----------------------------------------------------------------------

  // Create observables
  RooRealVar x("x","x",-5,5) ;
  RooRealVar y("y","y",-5,5) ;

  // Create function f(y) = a0 + a1*y
  RooRealVar a0("a0","a0",-0.5,-5,5) ;
  RooRealVar a1("a1","a1",-0.5,-1,1) ;
  RooPolyVar fy("fy","fy",y,RooArgSet(a0,a1)) ;

  // Creat gauss(x,f(y),s)
  RooRealVar sigma("sigma","width of gaussian",0.5) ;
  RooGaussian model("model","Gaussian with shifting mean",x,fy,sigma) ;


  // S a m p l e   d a t a ,   p l o t   d a t a   a n d   p d f   o n   x   a n d   y
  // ---------------------------------------------------------------------------------

  // Generate 10000 events in x and y from model
  RooDataSet *data = model.generate(RooArgSet(x,y),10000) ;

  // Plot x distribution of data and projection of model on x = Int(dy) model(x,y)
  RooPlot* xframe = x.frame() ;
  data->plotOn(xframe) ;
  model.plotOn(xframe) ;

  // Plot x distribution of data and projection of model on y = Int(dx) model(x,y)
  RooPlot* yframe = y.frame() ;
  data->plotOn(yframe) ;
  model.plotOn(yframe) ;

  // Make two-dimensional plot in x vs y
  TH1* hh_model = model.createHistogram("hh_model",x,Binning(50),YVar(y,Binning(50))) ;
  hh_model->SetLineColor(kBlue) ;


  regPlot(xframe,"rf301_plot1") ;
  regPlot(yframe,"rf302_plot2") ;
  regTH(hh_model,"rf302_model2d") ;

  delete data ;
  return kTRUE ;
  }
} ;


//////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #302
//
//  Utility functions classes available for use in tailoring
//  of composite (multidimensional) pdfs
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooFormulaVar.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooPolyVar.h"
#include "TCanvas.h"
#include "TH1.h"

using namespace RooFit ;


class TestBasic302 : public RooUnitTest
{
public:
  TestBasic302(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Sum and product utility functions",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   o b s e r v a b l e s ,   p a r a m e t e r s
  // -----------------------------------------------------------

  // Create observables
  RooRealVar x("x","x",-5,5) ;
  RooRealVar y("y","y",-5,5) ;

  // Create parameters
  RooRealVar a0("a0","a0",-1.5,-5,5) ;
  RooRealVar a1("a1","a1",-0.5,-1,1) ;
  RooRealVar sigma("sigma","width of gaussian",0.5) ;


  // U s i n g   R o o F o r m u l a V a r   t o   t a i l o r   p d f
  // -----------------------------------------------------------------------

  // Create interpreted function f(y) = a0 - a1*sqrt(10*abs(y))
  RooFormulaVar fy_1("fy_1","a0-a1*sqrt(10*abs(y))",RooArgSet(y,a0,a1)) ;

  // Create gauss(x,f(y),s)
  RooGaussian model_1("model_1","Gaussian with shifting mean",x,fy_1,sigma) ;



  // U s i n g   R o o P o l y V a r   t o   t a i l o r   p d f
  // -----------------------------------------------------------------------

  // Create polynomial function f(y) = a0 + a1*y
  RooPolyVar fy_2("fy_2","fy_2",y,RooArgSet(a0,a1)) ;

  // Create gauss(x,f(y),s)
  RooGaussian model_2("model_2","Gaussian with shifting mean",x,fy_2,sigma) ;



  // U s i n g   R o o A d d i t i o n   t o   t a i l o r   p d f
  // -----------------------------------------------------------------------

  // Create sum function f(y) = a0 + y
  RooAddition fy_3("fy_3","a0+y",RooArgSet(a0,y)) ;

  // Create gauss(x,f(y),s)
  RooGaussian model_3("model_3","Gaussian with shifting mean",x,fy_3,sigma) ;



  // U s i n g   R o o P r o d u c t   t o   t a i l o r   p d f
  // -----------------------------------------------------------------------

  // Create product function f(y) = a1*y
  RooProduct fy_4("fy_4","a1*y",RooArgSet(a1,y)) ;

  // Create gauss(x,f(y),s)
  RooGaussian model_4("model_4","Gaussian with shifting mean",x,fy_4,sigma) ;



  // P l o t   a l l   p d f s
  // ----------------------------

  // Make two-dimensional plots in x vs y
  TH1* hh_model_1 = model_1.createHistogram("hh_model_1",x,Binning(50),YVar(y,Binning(50))) ;
  TH1* hh_model_2 = model_2.createHistogram("hh_model_2",x,Binning(50),YVar(y,Binning(50))) ;
  TH1* hh_model_3 = model_3.createHistogram("hh_model_3",x,Binning(50),YVar(y,Binning(50))) ;
  TH1* hh_model_4 = model_4.createHistogram("hh_model_4",x,Binning(50),YVar(y,Binning(50))) ;
  hh_model_1->SetLineColor(kBlue) ;
  hh_model_2->SetLineColor(kBlue) ;
  hh_model_3->SetLineColor(kBlue) ;
  hh_model_4->SetLineColor(kBlue) ;

  regTH(hh_model_1,"rf202_model2d_1") ;
  regTH(hh_model_2,"rf202_model2d_2") ;
  regTH(hh_model_3,"rf202_model2d_3") ;
  regTH(hh_model_4,"rf202_model2d_4") ;

  return kTRUE ;
  }
} ;
/////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #303
//
// Use of tailored p.d.f as conditional p.d.fs.s
//
// pdf = gauss(x,f(y),sx | y ) with f(y) = a0 + a1*y
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolyVar.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TH1.h"
using namespace RooFit ;


class TestBasic303 : public RooUnitTest
{
public:

RooDataSet* makeFakeDataXY()
{
  RooRealVar x("x","x",-10,10) ;
  RooRealVar y("y","y",-10,10) ;
  RooArgSet coord(x,y) ;

  RooDataSet* d = new RooDataSet("d","d",RooArgSet(x,y)) ;

  for (int i=0 ; i<10000 ; i++) {
    Double_t tmpy = gRandom->Gaus(0,10) ;
    Double_t tmpx = gRandom->Gaus(0.5*tmpy,1) ;
    if (fabs(tmpy)<10 && fabs(tmpx)<10) {
      x = tmpx ;
      y = tmpy ;
      d->add(coord) ;
    }

  }

  return d ;
}



  TestBasic303(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Conditional use of F(x|y)",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // S e t u p   c o m p o s e d   m o d e l   g a u s s ( x , m ( y ) , s )
  // -----------------------------------------------------------------------

  // Create observables
  RooRealVar x("x","x",-10,10) ;
  RooRealVar y("y","y",-10,10) ;

  // Create function f(y) = a0 + a1*y
  RooRealVar a0("a0","a0",-0.5,-5,5) ;
  RooRealVar a1("a1","a1",-0.5,-1,1) ;
  RooPolyVar fy("fy","fy",y,RooArgSet(a0,a1)) ;

  // Creat gauss(x,f(y),s)
  RooRealVar sigma("sigma","width of gaussian",0.5,0.1,2.0) ;
  RooGaussian model("model","Gaussian with shifting mean",x,fy,sigma) ;


  // Obtain fake external experimental dataset with values for x and y
  RooDataSet* expDataXY = makeFakeDataXY() ;



  // G e n e r a t e   d a t a   f r o m   c o n d i t i o n a l   p . d . f   m o d e l ( x | y )
  // ---------------------------------------------------------------------------------------------

  // Make subset of experimental data with only y values
  RooDataSet* expDataY= (RooDataSet*) expDataXY->reduce(y) ;

  // Generate 10000 events in x obtained from _conditional_ model(x|y) with y values taken from experimental data
  RooDataSet *data = model.generate(x,ProtoData(*expDataY)) ;



  // F i t   c o n d i t i o n a l   p . d . f   m o d e l ( x | y )   t o   d a t a
  // ---------------------------------------------------------------------------------------------

  model.fitTo(*expDataXY,ConditionalObservables(y)) ;



  // P r o j e c t   c o n d i t i o n a l   p . d . f   o n   x   a n d   y   d i m e n s i o n s
  // ---------------------------------------------------------------------------------------------

  // Plot x distribution of data and projection of model on x = 1/Ndata sum(data(y_i)) model(x;y_i)
  RooPlot* xframe = x.frame() ;
  expDataXY->plotOn(xframe) ;
  model.plotOn(xframe,ProjWData(*expDataY)) ;


  // Speed up (and approximate) projection by using binned clone of data for projection
  RooAbsData* binnedDataY = expDataY->binnedClone() ;
  model.plotOn(xframe,ProjWData(*binnedDataY),LineColor(kCyan),LineStyle(kDotted),Name("Alt1")) ;


  // Show effect of projection with too coarse binning
  ((RooRealVar*)expDataY->get()->find("y"))->setBins(5) ;
  RooAbsData* binnedDataY2 = expDataY->binnedClone() ;
  model.plotOn(xframe,ProjWData(*binnedDataY2),LineColor(kRed),Name("Alt2")) ;


  regPlot(xframe,"rf303_plot1") ;

  delete binnedDataY ;
  delete binnedDataY2 ;
  delete expDataXY ;
  delete expDataY ;
  delete data ;

  return kTRUE ;
  }
} ;




/////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #304
//
// Simple uncorrelated multi-dimensional p.d.f.s
//
// pdf = gauss(x,mx,sx) * gauss(y,my,sy)
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;



class TestBasic304 : public RooUnitTest
{
public:
  TestBasic304(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Product operator p.d.f. with uncorrelated terms",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   c o m p o n e n t   p d f s   i n   x   a n d   y
  // ----------------------------------------------------------------

  // Create two p.d.f.s gaussx(x,meanx,sigmax) gaussy(y,meany,sigmay) and its variables
  RooRealVar x("x","x",-5,5) ;
  RooRealVar y("y","y",-5,5) ;

  RooRealVar meanx("mean1","mean of gaussian x",2) ;
  RooRealVar meany("mean2","mean of gaussian y",-2) ;
  RooRealVar sigmax("sigmax","width of gaussian x",1) ;
  RooRealVar sigmay("sigmay","width of gaussian y",5) ;

  RooGaussian gaussx("gaussx","gaussian PDF",x,meanx,sigmax) ;
  RooGaussian gaussy("gaussy","gaussian PDF",y,meany,sigmay) ;



  // C o n s t r u c t   u n c o r r e l a t e d   p r o d u c t   p d f
  // -------------------------------------------------------------------

  // Multiply gaussx and gaussy into a two-dimensional p.d.f. gaussxy
  RooProdPdf  gaussxy("gaussxy","gaussx*gaussy",RooArgList(gaussx,gaussy)) ;



  // S a m p l e   p d f ,   p l o t   p r o j e c t i o n   o n   x   a n d   y
  // ---------------------------------------------------------------------------

  // Generate 10000 events in x and y from gaussxy
  RooDataSet *data = gaussxy.generate(RooArgSet(x,y),10000) ;

  // Plot x distribution of data and projection of gaussxy on x = Int(dy) gaussxy(x,y)
  RooPlot* xframe = x.frame(Title("X projection of gauss(x)*gauss(y)")) ;
  data->plotOn(xframe) ;
  gaussxy.plotOn(xframe) ;

  // Plot x distribution of data and projection of gaussxy on y = Int(dx) gaussxy(x,y)
  RooPlot* yframe = y.frame(Title("Y projection of gauss(x)*gauss(y)")) ;
  data->plotOn(yframe) ;
  gaussxy.plotOn(yframe) ;

  regPlot(xframe,"rf304_plot1") ;
  regPlot(yframe,"rf304_plot2") ;

  delete data ;

  return kTRUE ;
  }
} ;



/////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #305
//
// Multi-dimensional p.d.f.s with conditional p.d.fs in product
//
// pdf = gauss(x,f(y),sx | y ) * gauss(y,ms,sx)    with f(y) = a0 + a1*y
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolyVar.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
using namespace RooFit ;



class TestBasic305 : public RooUnitTest
{
public:
  TestBasic305(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Product operator p.d.f. with conditional term",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   c o n d i t i o n a l   p d f   g x ( x | y )
  // -----------------------------------------------------------

  // Create observables
  RooRealVar x("x","x",-5,5) ;
  RooRealVar y("y","y",-5,5) ;

  // Create function f(y) = a0 + a1*y
  RooRealVar a0("a0","a0",-0.5,-5,5) ;
  RooRealVar a1("a1","a1",-0.5,-1,1) ;
  RooPolyVar fy("fy","fy",y,RooArgSet(a0,a1)) ;

  // Create gaussx(x,f(y),sx)
  RooRealVar sigmax("sigma","width of gaussian",0.5) ;
  RooGaussian gaussx("gaussx","Gaussian in x with shifting mean in y",x,fy,sigmax) ;



  // C r e a t e   p d f   g y ( y )
  // -----------------------------------------------------------

  // Create gaussy(y,0,5)
  RooGaussian gaussy("gaussy","Gaussian in y",y,RooConst(0),RooConst(3)) ;



  // C r e a t e   p r o d u c t   g x ( x | y ) * g y ( y )
  // -------------------------------------------------------

  // Create gaussx(x,sx|y) * gaussy(y)
  RooProdPdf model("model","gaussx(x|y)*gaussy(y)",gaussy,Conditional(gaussx,x)) ;



  // S a m p l e ,   f i t   a n d   p l o t   p r o d u c t   p d f
  // ---------------------------------------------------------------

  // Generate 1000 events in x and y from model
  RooDataSet *data = model.generate(RooArgSet(x,y),10000) ;

  // Plot x distribution of data and projection of model on x = Int(dy) model(x,y)
  RooPlot* xframe = x.frame() ;
  data->plotOn(xframe) ;
  model.plotOn(xframe) ;

  // Plot x distribution of data and projection of model on y = Int(dx) model(x,y)
  RooPlot* yframe = y.frame() ;
  data->plotOn(yframe) ;
  model.plotOn(yframe) ;

  // Make two-dimensional plot in x vs y
  TH1* hh_model = model.createHistogram("hh_model_rf305",x,Binning(50),YVar(y,Binning(50))) ;
  hh_model->SetLineColor(kBlue) ;

  regTH(hh_model,"rf305_model2d") ;
  regPlot(xframe,"rf305_plot1") ;
  regPlot(yframe,"rf305_plot2") ;

  delete data ;

  return kTRUE ;

  }
} ;



//////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #306
//
// Complete example with use of conditional p.d.f. with per-event errors
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooLandau.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH2D.h"
using namespace RooFit ;



class TestBasic306 : public RooUnitTest
{
public:
  TestBasic306(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Conditional use of per-event error p.d.f. F(t|dt)",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // B - p h y s i c s   p d f   w i t h   p e r - e v e n t  G a u s s i a n   r e s o l u t i o n
  // ----------------------------------------------------------------------------------------------

  // Observables
  RooRealVar dt("dt","dt",-10,10) ;
  RooRealVar dterr("dterr","per-event error on dt",0.01,10) ;

  // Build a gaussian resolution model scaled by the per-event error = gauss(dt,bias,sigma*dterr)
  RooRealVar bias("bias","bias",0,-10,10) ;
  RooRealVar sigma("sigma","per-event error scale factor",1,0.1,10) ;
  RooGaussModel gm("gm1","gauss model scaled bt per-event error",dt,bias,sigma,dterr) ;

  // Construct decay(dt) (x) gauss1(dt|dterr)
  RooRealVar tau("tau","tau",1.548) ;
  RooDecay decay_gm("decay_gm","decay",dt,tau,gm,RooDecay::DoubleSided) ;



  // C o n s t r u c t   f a k e   ' e x t e r n a l '   d a t a    w i t h   p e r - e v e n t   e r r o r
  // ------------------------------------------------------------------------------------------------------

  // Use landau p.d.f to get somewhat realistic distribution with long tail
  RooLandau pdfDtErr("pdfDtErr","pdfDtErr",dterr,RooConst(1),RooConst(0.25)) ;
  RooDataSet* expDataDterr = pdfDtErr.generate(dterr,10000) ;



  // S a m p l e   d a t a   f r o m   c o n d i t i o n a l   d e c a y _ g m ( d t | d t e r r )
  // ---------------------------------------------------------------------------------------------

  // Specify external dataset with dterr values to use decay_dm as conditional p.d.f.
  RooDataSet* data = decay_gm.generate(dt,ProtoData(*expDataDterr)) ;



  // F i t   c o n d i t i o n a l   d e c a y _ d m ( d t | d t e r r )
  // ---------------------------------------------------------------------

  // Specify dterr as conditional observable
  decay_gm.fitTo(*data,ConditionalObservables(dterr)) ;



  // P l o t   c o n d i t i o n a l   d e c a y _ d m ( d t | d t e r r )
  // ---------------------------------------------------------------------


  // Make two-dimensional plot of conditional p.d.f in (dt,dterr)
  TH1* hh_decay = decay_gm.createHistogram("hh_decay",dt,Binning(50),YVar(dterr,Binning(50))) ;
  hh_decay->SetLineColor(kBlue) ;


  // Plot decay_gm(dt|dterr) at various values of dterr
  RooPlot* frame = dt.frame(Title("Slices of decay(dt|dterr) at various dterr")) ;
  for (Int_t ibin=0 ; ibin<100 ; ibin+=20) {
    dterr.setBin(ibin) ;
    decay_gm.plotOn(frame,Normalization(5.),Name(Form("curve_slice_%d",ibin))) ;
  }


  // Make projection of data an dt
  RooPlot* frame2 = dt.frame(Title("Projection of decay(dt|dterr) on dt")) ;
  data->plotOn(frame2) ;

  // Make projection of decay(dt|dterr) on dt.
  //
  // Instead of integrating out dterr, make a weighted average of curves
  // at values dterr_i as given in the external dataset.
  // (The kTRUE argument bins the data before projection to speed up the process)
  decay_gm.plotOn(frame2,ProjWData(*expDataDterr,kTRUE)) ;


  regTH(hh_decay,"rf306_model2d") ;
  regPlot(frame,"rf306_plot1") ;
  regPlot(frame2,"rf306_plot2") ;

  delete expDataDterr ;
  delete data ;

  return kTRUE ;
  }
} ;

//////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #307
//
// Complete example with use of full p.d.f. with per-event errors
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooLandau.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
using namespace RooFit ;


class TestBasic307 : public RooUnitTest
{
public:
  TestBasic307(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Full per-event error p.d.f. F(t|dt)G(dt)",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // B - p h y s i c s   p d f   w i t h   p e r - e v e n t  G a u s s i a n   r e s o l u t i o n
  // ----------------------------------------------------------------------------------------------

  // Observables
  RooRealVar dt("dt","dt",-10,10) ;
  RooRealVar dterr("dterr","per-event error on dt",0.1,10) ;

  // Build a gaussian resolution model scaled by the per-event error = gauss(dt,bias,sigma*dterr)
  RooRealVar bias("bias","bias",0,-10,10) ;
  RooRealVar sigma("sigma","per-event error scale factor",1,0.1,10) ;
  RooGaussModel gm("gm1","gauss model scaled bt per-event error",dt,bias,sigma,dterr) ;

  // Construct decay(dt) (x) gauss1(dt|dterr)
  RooRealVar tau("tau","tau",1.548) ;
  RooDecay decay_gm("decay_gm","decay",dt,tau,gm,RooDecay::DoubleSided) ;



  // C o n s t r u c t   e m p i r i c a l   p d f   f o r   p e r - e v e n t   e r r o r
  // -----------------------------------------------------------------

  // Use landau p.d.f to get empirical distribution with long tail
  RooLandau pdfDtErr("pdfDtErr","pdfDtErr",dterr,RooConst(1),RooConst(0.25)) ;
  RooDataSet* expDataDterr = pdfDtErr.generate(dterr,10000) ;

  // Construct a histogram pdf to describe the shape of the dtErr distribution
  RooDataHist* expHistDterr = expDataDterr->binnedClone() ;
  RooHistPdf pdfErr("pdfErr","pdfErr",dterr,*expHistDterr) ;


  // C o n s t r u c t   c o n d i t i o n a l   p r o d u c t   d e c a y _ d m ( d t | d t e r r ) * p d f ( d t e r r )
  // ----------------------------------------------------------------------------------------------------------------------

  // Construct production of conditional decay_dm(dt|dterr) with empirical pdfErr(dterr)
  RooProdPdf model("model","model",pdfErr,Conditional(decay_gm,dt)) ;

  // (Alternatively you could also use the landau shape pdfDtErr)
  //RooProdPdf model("model","model",pdfDtErr,Conditional(decay_gm,dt)) ;



  // S a m p l e,   f i t   a n d   p l o t   p r o d u c t   m o d e l
  // ------------------------------------------------------------------

  // Specify external dataset with dterr values to use model_dm as conditional p.d.f.
  RooDataSet* data = model.generate(RooArgSet(dt,dterr),10000) ;



  // F i t   c o n d i t i o n a l   d e c a y _ d m ( d t | d t e r r )
  // ---------------------------------------------------------------------

  // Specify dterr as conditional observable
  model.fitTo(*data) ;



  // P l o t   c o n d i t i o n a l   d e c a y _ d m ( d t | d t e r r )
  // ---------------------------------------------------------------------


  // Make projection of data an dt
  RooPlot* frame = dt.frame(Title("Projection of model(dt|dterr) on dt")) ;
  data->plotOn(frame) ;
  model.plotOn(frame) ;


  regPlot(frame,"rf307_plot1") ;

  delete expDataDterr ;
  delete expHistDterr ;
  delete data ;

  return kTRUE ;

  }
} ;
/////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #308
//
// Examples on normalization of p.d.f.s,
// integration of p.d.fs, construction
// of cumulative distribution functions from p.d.f.s
// in two dimensions
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
using namespace RooFit ;


class TestBasic308 : public RooUnitTest
{
public:
  TestBasic308(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Normalization of p.d.f.s in 2D",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // S e t u p   m o d e l
  // ---------------------

  // Create observables x,y
  RooRealVar x("x","x",-10,10) ;
  RooRealVar y("y","y",-10,10) ;

  // Create p.d.f. gaussx(x,-2,3), gaussy(y,2,2)
  RooGaussian gx("gx","gx",x,RooConst(-2),RooConst(3)) ;
  RooGaussian gy("gy","gy",y,RooConst(+2),RooConst(2)) ;

  // Create gxy = gx(x)*gy(y)
  RooProdPdf gxy("gxy","gxy",RooArgSet(gx,gy)) ;



  // R e t r i e v e   r a w  &   n o r m a l i z e d   v a l u e s   o f   R o o F i t   p . d . f . s
  // --------------------------------------------------------------------------------------------------

  // Return 'raw' unnormalized value of gx
  regValue(gxy.getVal(),"rf308_gxy") ;

  // Return value of gxy normalized over x _and_ y in range [-10,10]
  RooArgSet nset_xy(x,y) ;
  regValue(gxy.getVal(&nset_xy),"rf308_gx_Norm[x,y]") ;

  // Create object representing integral over gx
  // which is used to calculate  gx_Norm[x,y] == gx / gx_Int[x,y]
  RooAbsReal* igxy = gxy.createIntegral(RooArgSet(x,y)) ;
  regValue(igxy->getVal(),"rf308_gx_Int[x,y]") ;

  // NB: it is also possible to do the following

  // Return value of gxy normalized over x in range [-10,10] (i.e. treating y as parameter)
  RooArgSet nset_x(x) ;
  regValue(gxy.getVal(&nset_x),"rf308_gx_Norm[x]") ;

  // Return value of gxy normalized over y in range [-10,10] (i.e. treating x as parameter)
  RooArgSet nset_y(y) ;
  regValue(gxy.getVal(&nset_y),"rf308_gx_Norm[y]") ;



  // I n t e g r a t e   n o r m a l i z e d   p d f   o v e r   s u b r a n g e
  // ----------------------------------------------------------------------------

  // Define a range named "signal" in x from -5,5
  x.setRange("signal",-5,5) ;
  y.setRange("signal",-3,3) ;

  // Create an integral of gxy_Norm[x,y] over x and y in range "signal"
  // This is the fraction of of p.d.f. gxy_Norm[x,y] which is in the
  // range named "signal"
  RooAbsReal* igxy_sig = gxy.createIntegral(RooArgSet(x,y),NormSet(RooArgSet(x,y)),Range("signal")) ;
  regValue(igxy_sig->getVal(),"rf308_gx_Int[x,y|signal]_Norm[x,y]") ;



  // C o n s t r u c t   c u m u l a t i v e   d i s t r i b u t i o n   f u n c t i o n   f r o m   p d f
  // -----------------------------------------------------------------------------------------------------

  // Create the cumulative distribution function of gx
  // i.e. calculate Int[-10,x] gx(x') dx'
  RooAbsReal* gxy_cdf = gxy.createCdf(RooArgSet(x,y)) ;

  // Plot cdf of gx versus x
  TH1* hh_cdf = gxy_cdf->createHistogram("hh_cdf",x,Binning(40),YVar(y,Binning(40))) ;

  regTH(hh_cdf,"rf308_cdf") ;

  delete igxy_sig ;
  delete igxy ;
  delete gxy_cdf ;

  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #309
//
// Projecting p.d.f and data slices in discrete observables
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooBMixDecay.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic310 : public RooUnitTest
{
public:
  TestBasic310(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Data and p.d.f projection in category slice",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   B   d e c a y   p d f   w it h   m i x i n g
  // ----------------------------------------------------------

  // Decay time observables
  RooRealVar dt("dt","dt",-20,20) ;

  // Discrete observables mixState (B0tag==B0reco?) and tagFlav (B0tag==B0(bar)?)
  RooCategory mixState("mixState","B0/B0bar mixing state") ;
  RooCategory tagFlav("tagFlav","Flavour of the tagged B0") ;

  // Define state labels of discrete observables
  mixState.defineType("mixed",-1) ;
  mixState.defineType("unmixed",1) ;
  tagFlav.defineType("B0",1) ;
  tagFlav.defineType("B0bar",-1) ;

  // Model parameters
  RooRealVar dm("dm","delta m(B)",0.472,0.,1.0) ;
  RooRealVar tau("tau","B0 decay time",1.547,1.0,2.0) ;
  RooRealVar w("w","Flavor Mistag rate",0.03,0.0,1.0) ;
  RooRealVar dw("dw","Flavor Mistag rate difference between B0 and B0bar",0.01) ;

  // Build a gaussian resolution model
  RooRealVar bias1("bias1","bias1",0) ;
  RooRealVar sigma1("sigma1","sigma1",0.01) ;
  RooGaussModel gm1("gm1","gauss model 1",dt,bias1,sigma1) ;

  // Construct a decay pdf, smeared with single gaussian resolution model
  RooBMixDecay bmix_gm1("bmix","decay",dt,mixState,tagFlav,tau,dm,w,dw,gm1,RooBMixDecay::DoubleSided) ;

  // Generate BMixing data with above set of event errors
  RooDataSet *data = bmix_gm1.generate(RooArgSet(dt,tagFlav,mixState),20000) ;



  // P l o t   f u l l   d e c a y   d i s t r i b u t i o n
  // ----------------------------------------------------------

  // Create frame, plot data and pdf projection (integrated over tagFlav and mixState)
  RooPlot* frame = dt.frame(Title("Inclusive decay distribution")) ;
  data->plotOn(frame) ;
  bmix_gm1.plotOn(frame) ;



  // P l o t   d e c a y   d i s t r .   f o r   m i x e d   a n d   u n m i x e d   s l i c e   o f   m i x S t a t e
  // ------------------------------------------------------------------------------------------------------------------

  // Create frame, plot data (mixed only)
  RooPlot* frame2 = dt.frame(Title("Decay distribution of mixed events")) ;
  data->plotOn(frame2,Cut("mixState==mixState::mixed")) ;

  // Position slice in mixState at "mixed" and plot slice of pdf in mixstate over data (integrated over tagFlav)
  bmix_gm1.plotOn(frame2,Slice(mixState,"mixed")) ;

  // Create frame, plot data (unmixed only)
  RooPlot* frame3 = dt.frame(Title("Decay distribution of unmixed events")) ;
  data->plotOn(frame3,Cut("mixState==mixState::unmixed")) ;

  // Position slice in mixState at "unmixed" and plot slice of pdf in mixstate over data (integrated over tagFlav)
  bmix_gm1.plotOn(frame3,Slice(mixState,"unmixed")) ;


  regPlot(frame,"rf310_plot1") ;
  regPlot(frame2,"rf310_plot2") ;
  regPlot(frame3,"rf310_plot3") ;

  delete data ;

  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #310
//
// Projecting p.d.f and data ranges in continuous observables
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic311 : public RooUnitTest
{
public:
  TestBasic311(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Data and p.d.f projection in sub range",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   3 D   p d f   a n d   d a t a
  // -------------------------------------------

  // Create observables
  RooRealVar x("x","x",-5,5) ;
  RooRealVar y("y","y",-5,5) ;
  RooRealVar z("z","z",-5,5) ;

  // Create signal pdf gauss(x)*gauss(y)*gauss(z)
  RooGaussian gx("gx","gx",x,RooConst(0),RooConst(1)) ;
  RooGaussian gy("gy","gy",y,RooConst(0),RooConst(1)) ;
  RooGaussian gz("gz","gz",z,RooConst(0),RooConst(1)) ;
  RooProdPdf sig("sig","sig",RooArgSet(gx,gy,gz)) ;

  // Create background pdf poly(x)*poly(y)*poly(z)
  RooPolynomial px("px","px",x,RooArgSet(RooConst(-0.1),RooConst(0.004))) ;
  RooPolynomial py("py","py",y,RooArgSet(RooConst(0.1),RooConst(-0.004))) ;
  RooPolynomial pz("pz","pz",z) ;
  RooProdPdf bkg("bkg","bkg",RooArgSet(px,py,pz)) ;

  // Create composite pdf sig+bkg
  RooRealVar fsig("fsig","signal fraction",0.1,0.,1.) ;
  RooAddPdf model("model","model",RooArgList(sig,bkg),fsig) ;

  RooDataSet* data = model.generate(RooArgSet(x,y,z),20000) ;



  // P r o j e c t   p d f   a n d   d a t a   o n   x
  // -------------------------------------------------

  // Make plain projection of data and pdf on x observable
  RooPlot* frame = x.frame(Title("Projection of 3D data and pdf on X"),Bins(40)) ;
  data->plotOn(frame) ;
  model.plotOn(frame) ;



  // P r o j e c t   p d f   a n d   d a t a   o n   x   i n   s i g n a l   r a n g e
  // ----------------------------------------------------------------------------------

  // Define signal region in y and z observables
  y.setRange("sigRegion",-1,1) ;
  z.setRange("sigRegion",-1,1) ;

  // Make plot frame
  RooPlot* frame2 = x.frame(Title("Same projection on X in signal range of (Y,Z)"),Bins(40)) ;

  // Plot subset of data in which all observables are inside "sigRegion"
  // For observables that do not have an explicit "sigRegion" range defined (e.g. observable)
  // an implicit definition is used that is identical to the full range (i.e. [-5,5] for x)
  data->plotOn(frame2,CutRange("sigRegion")) ;

  // Project model on x, integrating projected observables (y,z) only in "sigRegion"
  model.plotOn(frame2,ProjectionRange("sigRegion")) ;


  regPlot(frame,"rf311_plot1") ;
  regPlot(frame2,"rf312_plot2") ;

  delete data ;

  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #312
//
// Performing fits in multiple (disjoint) ranges in one or more dimensions
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooFitResult.h"
using namespace RooFit ;


class TestBasic312 : public RooUnitTest
{
public:
  TestBasic312(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Fit in multiple rectangular ranges",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   2 D   p d f   a n d   d a t a
  // -------------------------------------------

  // Define observables x,y
  RooRealVar x("x","x",-10,10) ;
  RooRealVar y("y","y",-10,10) ;

  // Construct the signal pdf gauss(x)*gauss(y)
  RooRealVar mx("mx","mx",1,-10,10) ;
  RooRealVar my("my","my",1,-10,10) ;

  RooGaussian gx("gx","gx",x,mx,RooConst(1)) ;
  RooGaussian gy("gy","gy",y,my,RooConst(1)) ;

  RooProdPdf sig("sig","sig",gx,gy) ;

  // Construct the background pdf (flat in x,y)
  RooPolynomial px("px","px",x) ;
  RooPolynomial py("py","py",y) ;
  RooProdPdf bkg("bkg","bkg",px,py) ;

  // Construct the composite model sig+bkg
  RooRealVar f("f","f",0.,1.) ;
  RooAddPdf model("model","model",RooArgList(sig,bkg),f) ;

  // Sample 10000 events in (x,y) from the model
  RooDataSet* modelData = model.generate(RooArgSet(x,y),10000) ;



  // D e f i n e   s i g n a l   a n d   s i d e b a n d   r e g i o n s
  // -------------------------------------------------------------------

  // Construct the SideBand1,SideBand2,Signal regions
  //
  //                    |
  //      +-------------+-----------+
  //      |             |           |
  //      |    Side     |   Sig     |
  //      |    Band1    |   nal     |
  //      |             |           |
  //    --+-------------+-----------+--
  //      |                         |
  //      |           Side          |
  //      |           Band2         |
  //      |                         |
  //      +-------------+-----------+
  //                    |

  x.setRange("SB1",-10,+10) ;
  y.setRange("SB1",-10,0) ;

  x.setRange("SB2",-10,0) ;
  y.setRange("SB2",0,+10) ;

  x.setRange("SIG",0,+10) ;
  y.setRange("SIG",0,+10) ;

  x.setRange("FULL",-10,+10) ;
  y.setRange("FULL",-10,+10) ;


  // P e r f o r m   f i t s   i n   i n d i v i d u a l   s i d e b a n d   r e g i o n s
  // -------------------------------------------------------------------------------------

  // Perform fit in SideBand1 region (RooAddPdf coefficients will be interpreted in full range)
  RooFitResult* r_sb1 = model.fitTo(*modelData,Range("SB1"),Save()) ;

  // Perform fit in SideBand2 region (RooAddPdf coefficients will be interpreted in full range)
  RooFitResult* r_sb2 = model.fitTo(*modelData,Range("SB2"),Save()) ;



  // P e r f o r m   f i t s   i n   j o i n t    s i d e b a n d   r e g i o n s
  // -----------------------------------------------------------------------------

  // Now perform fit to joint 'L-shaped' sideband region 'SB1|SB2'
  // (RooAddPdf coefficients will be interpreted in full range)
  RooFitResult* r_sb12 = model.fitTo(*modelData,Range("SB1,SB2"),Save()) ;


  regResult(r_sb1,"rf312_fit_sb1") ;
  regResult(r_sb2,"rf312_fit_sb2") ;
  regResult(r_sb12,"rf312_fit_sb12") ;

  delete modelData ;

  return kTRUE ;

  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #313
//
// Working with parameterized ranges to define non-rectangular regions
// for fitting and integration
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic313 : public RooUnitTest
{
public:
  TestBasic313(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Integration over non-rectangular regions",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   3 D   p d f
  // -------------------------

  // Define observable (x,y,z)
  RooRealVar x("x","x",0,10) ;
  RooRealVar y("y","y",0,10) ;
  RooRealVar z("z","z",0,10) ;

  // Define 3 dimensional pdf
  RooRealVar z0("z0","z0",-0.1,1) ;
  RooPolynomial px("px","px",x,RooConst(0)) ;
  RooPolynomial py("py","py",y,RooConst(0)) ;
  RooPolynomial pz("pz","pz",z,z0) ;
  RooProdPdf pxyz("pxyz","pxyz",RooArgSet(px,py,pz)) ;



  // D e f i n e d   n o n - r e c t a n g u l a r   r e g i o n   R   i n   ( x , y , z )
  // -------------------------------------------------------------------------------------

  //
  // R = Z[0 - 0.1*Y^2] * Y[0.1*X - 0.9*X] * X[0 - 10]
  //

  // Construct range parameterized in "R" in y [ 0.1*x, 0.9*x ]
  RooFormulaVar ylo("ylo","0.1*x",x) ;
  RooFormulaVar yhi("yhi","0.9*x",x) ;
  y.setRange("R",ylo,yhi) ;

  // Construct parameterized ranged "R" in z [ 0, 0.1*y^2 ]
  RooFormulaVar zlo("zlo","0.0*y",y) ;
  RooFormulaVar zhi("zhi","0.1*y*y",y) ;
  z.setRange("R",zlo,zhi) ;



  // C a l c u l a t e   i n t e g r a l   o f   n o r m a l i z e d   p d f   i n   R
  // ----------------------------------------------------------------------------------

  // Create integral over normalized pdf model over x,y,z in "R" region
  RooAbsReal* intPdf = pxyz.createIntegral(RooArgSet(x,y,z),RooArgSet(x,y,z),"R") ;

  // Plot value of integral as function of pdf parameter z0
  RooPlot* frame = z0.frame(Title("Integral of pxyz over x,y,z in region R")) ;
  intPdf->plotOn(frame) ;


  regPlot(frame,"rf313_plot1") ;

  delete intPdf ;

  return kTRUE;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #314
//
// Working with parameterized ranges in a fit. This an example of a
// fit with an acceptance that changes per-event
//
//  pdf = exp(-t/tau) with t[tmin,5]
//
//  where t and tmin are both observables in the dataset
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooFitResult.h"

using namespace RooFit ;


class TestBasic314 : public RooUnitTest
{
public:
  TestBasic314(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Fit with non-rectangular observable boundaries",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // D e f i n e   o b s e r v a b l e s   a n d   d e c a y   p d f
  // ---------------------------------------------------------------

  // Declare observables
  RooRealVar t("t","t",0,5) ;
  RooRealVar tmin("tmin","tmin",0,0,5) ;

  // Make parameterized range in t : [tmin,5]
  t.setRange(tmin,RooConst(t.getMax())) ;

  // Make pdf
  RooRealVar tau("tau","tau",-1.54,-10,-0.1) ;
  RooExponential model("model","model",t,tau) ;



  // C r e a t e   i n p u t   d a t a
  // ------------------------------------

  // Generate complete dataset without acceptance cuts (for reference)
  RooDataSet* dall = model.generate(t,10000) ;

  // Generate a (fake) prototype dataset for acceptance limit values
  RooDataSet* tmp = RooGaussian("gmin","gmin",tmin,RooConst(0),RooConst(0.5)).generate(tmin,5000) ;

  // Generate dataset with t values that observe (t>tmin)
  RooDataSet* dacc = model.generate(t,ProtoData(*tmp)) ;



  // F i t   p d f   t o   d a t a   i n   a c c e p t a n c e   r e g i o n
  // -----------------------------------------------------------------------

  RooFitResult* r = model.fitTo(*dacc,Save()) ;



  // P l o t   f i t t e d   p d f   o n   f u l l   a n d   a c c e p t e d   d a t a
  // ---------------------------------------------------------------------------------

  // Make plot frame, add datasets and overlay model
  RooPlot* frame = t.frame(Title("Fit to data with per-event acceptance")) ;
  dall->plotOn(frame,MarkerColor(kRed),LineColor(kRed)) ;
  model.plotOn(frame) ;
  dacc->plotOn(frame,Name("dacc")) ;

  // Print fit results to demonstrate absence of bias
  regResult(r,"rf314_fit") ;
  regPlot(frame,"rf314_plot1") ;

  delete tmp ;
  delete dacc ;
  delete dall ;

  return kTRUE;
  }
} ;


//////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #315
//
// Marginizalization of multi-dimensional p.d.f.s through integration
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPolyVar.h"
#include "TH1.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooConstVar.h"
using namespace RooFit ;



class TestBasic315 : public RooUnitTest
{
public:
  TestBasic315(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("P.d.f. marginalization through integration",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    // C r e a t e   p d f   m ( x , y )  =  g x ( x | y ) * g ( y )
    // --------------------------------------------------------------

    // Increase default precision of numeric integration
    // as this exercise has high sensitivity to numeric integration precision
    RooAbsPdf::defaultIntegratorConfig()->setEpsRel(1e-8) ;
    RooAbsPdf::defaultIntegratorConfig()->setEpsAbs(1e-8) ;

    // Create observables
    RooRealVar x("x","x",-5,5) ;
    RooRealVar y("y","y",-2,2) ;

    // Create function f(y) = a0 + a1*y
    RooRealVar a0("a0","a0",0) ;
    RooRealVar a1("a1","a1",-1.5,-3,1) ;
    RooPolyVar fy("fy","fy",y,RooArgSet(a0,a1)) ;

    // Create gaussx(x,f(y),sx)
    RooRealVar sigmax("sigmax","width of gaussian",0.5) ;
    RooGaussian gaussx("gaussx","Gaussian in x with shifting mean in y",x,fy,sigmax) ;

    // Create gaussy(y,0,2)
    RooGaussian gaussy("gaussy","Gaussian in y",y,RooConst(0),RooConst(2)) ;

    // Create gaussx(x,sx|y) * gaussy(y)
    RooProdPdf model("model","gaussx(x|y)*gaussy(y)",gaussy,Conditional(gaussx,x)) ;



    // M a r g i n a l i z e   m ( x , y )   t o   m ( x )
    // ----------------------------------------------------

    // modelx(x) = Int model(x,y) dy
    RooAbsPdf* modelx = model.createProjection(y) ;



    // U s e   m a r g i n a l i z e d   p . d . f .   a s   r e g u l a r   1 - D   p . d . f .
    // ------------------------------------------------------------------------------------------

    // Sample 1000 events from modelx
    RooAbsData* data = modelx->generateBinned(x,1000) ;

    // Fit modelx to toy data
    modelx->fitTo(*data) ;

    // Plot modelx over data
    RooPlot* frame = x.frame(40) ;
    data->plotOn(frame) ;
    modelx->plotOn(frame) ;

    regPlot(frame,"rf315_frame") ;

    delete data ;
    delete modelx ;

    return kTRUE ;
  }
} ;





//////////////////////////////////////////////////////////////////////////
//
// 'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #316
//
// Using the likelihood ratio techique to construct a signal enhanced
// one-dimensional projection of a multi-dimensional p.d.f.
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic316 : public RooUnitTest
{
public:
  TestBasic316(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Likelihood ratio projection plot",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   3 D   p d f   a n d   d a t a
  // -------------------------------------------

  // Create observables
  RooRealVar x("x","x",-5,5) ;
  RooRealVar y("y","y",-5,5) ;
  RooRealVar z("z","z",-5,5) ;

  // Create signal pdf gauss(x)*gauss(y)*gauss(z)
  RooGaussian gx("gx","gx",x,RooConst(0),RooConst(1)) ;
  RooGaussian gy("gy","gy",y,RooConst(0),RooConst(1)) ;
  RooGaussian gz("gz","gz",z,RooConst(0),RooConst(1)) ;
  RooProdPdf sig("sig","sig",RooArgSet(gx,gy,gz)) ;

  // Create background pdf poly(x)*poly(y)*poly(z)
  RooPolynomial px("px","px",x,RooArgSet(RooConst(-0.1),RooConst(0.004))) ;
  RooPolynomial py("py","py",y,RooArgSet(RooConst(0.1),RooConst(-0.004))) ;
  RooPolynomial pz("pz","pz",z) ;
  RooProdPdf bkg("bkg","bkg",RooArgSet(px,py,pz)) ;

  // Create composite pdf sig+bkg
  RooRealVar fsig("fsig","signal fraction",0.1,0.,1.) ;
  RooAddPdf model("model","model",RooArgList(sig,bkg),fsig) ;

  RooDataSet* data = model.generate(RooArgSet(x,y,z),20000) ;



  // P r o j e c t   p d f   a n d   d a t a   o n   x
  // -------------------------------------------------

  // Make plain projection of data and pdf on x observable
  RooPlot* frame = x.frame(Title("Projection of 3D data and pdf on X"),Bins(40)) ;
  data->plotOn(frame) ;
  model.plotOn(frame) ;



  // D e f i n e   p r o j e c t e d   s i g n a l   l i k e l i h o o d   r a t i o
  // ----------------------------------------------------------------------------------

  // Calculate projection of signal and total likelihood on (y,z) observables
  // i.e. integrate signal and composite model over x
  RooAbsPdf* sigyz = sig.createProjection(x) ;
  RooAbsPdf* totyz = model.createProjection(x) ;

  // Construct the log of the signal / signal+background probability
  RooFormulaVar llratio_func("llratio","log10(@0)-log10(@1)",RooArgList(*sigyz,*totyz)) ;



  // P l o t   d a t a   w i t h   a   L L r a t i o   c u t
  // -------------------------------------------------------

  // Calculate the llratio value for each event in the dataset
  data->addColumn(llratio_func) ;

  // Extract the subset of data with large signal likelihood
  RooDataSet* dataSel = (RooDataSet*) data->reduce(Cut("llratio>0.7")) ;

  // Make plot frame
  RooPlot* frame2 = x.frame(Title("Same projection on X with LLratio(y,z)>0.7"),Bins(40)) ;

  // Plot select data on frame
  dataSel->plotOn(frame2) ;



  // M a k e   M C   p r o j e c t i o n   o f   p d f   w i t h   s a m e   L L r a t i o   c u t
  // ---------------------------------------------------------------------------------------------

  // Generate large number of events for MC integration of pdf projection
  RooDataSet* mcprojData = model.generate(RooArgSet(x,y,z),10000) ;

  // Calculate LL ratio for each generated event and select MC events with llratio)0.7
  mcprojData->addColumn(llratio_func) ;
  RooDataSet* mcprojDataSel = (RooDataSet*) mcprojData->reduce(Cut("llratio>0.7")) ;

  // Project model on x, integrating projected observables (y,z) with Monte Carlo technique
  // on set of events with the same llratio cut as was applied to data
  model.plotOn(frame2,ProjWData(*mcprojDataSel)) ;


  regPlot(frame,"rf316_plot1") ;
  regPlot(frame2,"rf316_plot2") ;

  delete data ;
  delete dataSel ;
  delete mcprojData ;
  delete mcprojDataSel ;
  delete sigyz ;
  delete totyz ;

  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'DATA AND CATEGORIES' RooFit tutorial macro #402
//
// Tools for manipulation of (un)binned datasets
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TFile.h"
using namespace RooFit ;


class TestBasic402 : public RooUnitTest
{
public:
  TestBasic402(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Basic operations on datasets",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // Binned (RooDataHist) and unbinned datasets (RooDataSet) share
  // many properties and inherit from a common abstract base class
  // (RooAbsData), that provides an interface for all operations
  // that can be performed regardless of the data format

  RooRealVar  x("x","x",-10,10) ;
  RooRealVar  y("y","y", 0, 40) ;
  RooCategory c("c","c") ;
  c.defineType("Plus",+1) ;
  c.defineType("Minus",-1) ;



  // B a s i c   O p e r a t i o n s   o n   u n b i n n e d   d a t a s e t s
  // --------------------------------------------------------------

  // RooDataSet is an unbinned dataset (a collection of points in N-dimensional space)
  RooDataSet d("d","d",RooArgSet(x,y,c)) ;

  // Unlike RooAbsArgs (RooAbsPdf,RooFormulaVar,....) datasets are not attached to
  // the variables they are constructed from. Instead they are attached to an internal
  // clone of the supplied set of arguments

  // Fill d with dummy values
  Int_t i ;
  for (i=0 ; i<1000 ; i++) {
    x = i/50 - 10 ;
    y = sqrt(1.0*i) ;
    c.setLabel((i%2)?"Plus":"Minus") ;

    // We must explicitly refer to x,y,c here to pass the values because
    // d is not linked to them (as explained above)
    d.add(RooArgSet(x,y,c)) ;
  }


  // R e d u c i n g ,   A p p e n d i n g   a n d   M e r g i n g
  // -------------------------------------------------------------

  // The reduce() function returns a new dataset which is a subset of the original
  RooDataSet* d1 = (RooDataSet*) d.reduce(RooArgSet(x,c)) ;
  RooDataSet* d2 = (RooDataSet*) d.reduce(RooArgSet(y)) ;
  RooDataSet* d3 = (RooDataSet*) d.reduce("y>5.17") ;
  RooDataSet* d4 = (RooDataSet*) d.reduce(RooArgSet(x,c),"y>5.17") ;

  regValue(d3->numEntries(),"rf403_nd3") ;
  regValue(d4->numEntries(),"rf403_nd4") ;

  // The merge() function adds two data set column-wise
  d1->merge(d2) ;

  // The append() function addes two datasets row-wise
  d1->append(*d3) ;

  regValue(d1->numEntries(),"rf403_nd1") ;




  // O p e r a t i o n s   o n   b i n n e d   d a t a s e t s
  // ---------------------------------------------------------

  // A binned dataset can be constructed empty, from an unbinned dataset, or
  // from a ROOT native histogram (TH1,2,3)

  // The binning of real variables (like x,y) is done using their fit range
  // 'get/setRange()' and number of specified fit bins 'get/setBins()'.
  // Category dimensions of binned datasets get one bin per defined category state
  x.setBins(10) ;
  y.setBins(10) ;
  RooDataHist dh("dh","binned version of d",RooArgSet(x,y),d) ;

  RooPlot* yframe = y.frame(Bins(10),Title("Operations on binned datasets")) ;
  dh.plotOn(yframe) ; // plot projection of 2D binned data on y

  // Reduce the 2-dimensional binned dataset to a 1-dimensional binned dataset
  //
  // All reduce() methods are interfaced in RooAbsData. All reduction techniques
  // demonstrated on unbinned datasets can be applied to binned datasets as well.
  RooDataHist* dh2 = (RooDataHist*) dh.reduce(y,"x>0") ;

  // Add dh2 to yframe and redraw
  dh2->plotOn(yframe,LineColor(kRed),MarkerColor(kRed),Name("dh2")) ;

  regPlot(yframe,"rf402_plot1") ;

  delete d1 ;
  delete d2 ;
  delete d3 ;
  delete d4 ;
  delete dh2 ;
  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'DATA AND CATEGORIES' RooFit tutorial macro #403
//
// Using weights in unbinned datasets
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooFitResult.h"
using namespace RooFit ;


class TestBasic403 : public RooUnitTest
{
public:
  TestBasic403(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Fits with weighted datasets",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   o b s e r v a b l e   a n d   u n w e i g h t e d   d a t a s e t
  // -------------------------------------------------------------------------------

  // Declare observable
  RooRealVar x("x","x",-10,10) ;
  x.setBins(40) ;

  // Construction a uniform pdf
  RooPolynomial p0("px","px",x) ;

  // Sample 1000 events from pdf
  RooDataSet* data = p0.generate(x,1000) ;



  // C a l c u l a t e   w e i g h t   a n d   m a k e   d a t a s e t   w e i g h t e d
  // -----------------------------------------------------------------------------------

  // Construct formula to calculate (fake) weight for events
  RooFormulaVar wFunc("w","event weight","(x*x+10)",x) ;

  // Add column with variable w to previously generated dataset
  RooRealVar* w = (RooRealVar*) data->addColumn(wFunc) ;

  // Instruct dataset d in interpret w as event weight rather than as observable
  RooDataSet dataw(data->GetName(),data->GetTitle(),data,*data->get(),0,w->GetName()) ;
  //data->setWeightVar(*w) ;


  // U n b i n n e d   M L   f i t   t o   w e i g h t e d   d a t a
  // ---------------------------------------------------------------

  // Construction quadratic polynomial pdf for fitting
  RooRealVar a0("a0","a0",1) ;
  RooRealVar a1("a1","a1",0,-1,1) ;
  RooRealVar a2("a2","a2",1,0,10) ;
  RooPolynomial p2("p2","p2",x,RooArgList(a0,a1,a2),0) ;

  // Fit quadratic polynomial to weighted data

  // NOTE: Maximum likelihood fit to weighted data does in general
  //       NOT result in correct error estimates, unless individual
  //       event weights represent Poisson statistics themselves.
  //       In general, parameter error reflect precision of SumOfWeights
  //       events rather than NumEvents events. See comparisons below
  RooFitResult* r_ml_wgt = p2.fitTo(dataw,Save()) ;



  // P l o t   w e i g h e d   d a t a   a n d   f i t   r e s u l t
  // ---------------------------------------------------------------

  // Construct plot frame
  RooPlot* frame = x.frame(Title("Unbinned ML fit, binned chi^2 fit to weighted data")) ;

  // Plot data using sum-of-weights-squared error rather than Poisson errors
  dataw.plotOn(frame,DataError(RooAbsData::SumW2)) ;

  // Overlay result of 2nd order polynomial fit to weighted data
  p2.plotOn(frame) ;



  // M L  F i t   o f   p d f   t o   e q u i v a l e n t  u n w e i g h t e d   d a t a s e t
  // -----------------------------------------------------------------------------------------

  // Construct a pdf with the same shape as p0 after weighting
  RooGenericPdf genPdf("genPdf","x*x+10",x) ;

  // Sample a dataset with the same number of events as data
  RooDataSet* data2 = genPdf.generate(x,1000) ;

  // Sample a dataset with the same number of weights as data
  RooDataSet* data3 = genPdf.generate(x,43000) ;

  // Fit the 2nd order polynomial to both unweighted datasets and save the results for comparison
  RooFitResult* r_ml_unw10 = p2.fitTo(*data2,Save()) ;
  RooFitResult* r_ml_unw43 = p2.fitTo(*data3,Save()) ;


  // C h i 2   f i t   o f   p d f   t o   b i n n e d   w e i g h t e d   d a t a s e t
  // ------------------------------------------------------------------------------------

  // Construct binned clone of unbinned weighted dataset
  RooDataHist* binnedData = dataw.binnedClone() ;

  // Perform chi2 fit to binned weighted dataset using sum-of-weights errors
  //
  // NB: Within the usual approximations of a chi2 fit, a chi2 fit to weighted
  // data using sum-of-weights-squared errors does give correct error
  // estimates
  RooChi2Var chi2("chi2","chi2",p2,*binnedData) ;
  RooMinuit m(chi2) ;
  m.migrad() ;
  m.hesse() ;

  // Plot chi^2 fit result on frame as well
  RooFitResult* r_chi2_wgt = m.save() ;
  p2.plotOn(frame,LineStyle(kDashed),LineColor(kRed),Name("p2_alt")) ;



  // C o m p a r e   f i t   r e s u l t s   o f   c h i 2 , M L   f i t s   t o   ( u n ) w e i g h t e d   d a t a
  // ---------------------------------------------------------------------------------------------------------------

  // Note that ML fit on 1Kevt of weighted data is closer to result of ML fit on 43Kevt of unweighted data
  // than to 1Kevt of unweighted data, whereas the reference chi^2 fit with SumW2 error gives a result closer to
  // that of an unbinned ML fit to 1Kevt of unweighted data.

  regResult(r_ml_unw10,"rf403_ml_unw10") ;
  regResult(r_ml_unw43,"rf403_ml_unw43") ;
  regResult(r_ml_wgt  ,"rf403_ml_wgt") ;
  regResult(r_chi2_wgt ,"rf403_ml_chi2") ;
  regPlot(frame,"rf403_plot1") ;

  delete binnedData ;
  delete data ;
  delete data2 ;
  delete data3 ;

  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'DATA AND CATEGORIES' RooFit tutorial macro #404
//
// Working with RooCategory objects to describe discrete variables
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPolynomial.h"
#include "RooCategory.h"
#include "Roo1DTable.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic404 : public RooUnitTest
{
public:
  TestBasic404(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Categories basic functionality",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C o n s t r u c t    a   c a t e g o r y   w i t h   l a b e l s
  // ----------------------------------------------------------------

  // Define a category with labels only
  RooCategory tagCat("tagCat","Tagging category") ;
  tagCat.defineType("Lepton") ;
  tagCat.defineType("Kaon") ;
  tagCat.defineType("NetTagger-1") ;
  tagCat.defineType("NetTagger-2") ;



  // C o n s t r u c t    a   c a t e g o r y   w i t h   l a b e l s   a n d   i n d e c e s
  // ----------------------------------------------------------------------------------------

  // Define a category with explicitly numbered states
  RooCategory b0flav("b0flav","B0 flavour eigenstate") ;
  b0flav.defineType("B0",-1) ;
  b0flav.defineType("B0bar",1) ;



  // G e n e r a t e   d u m m y   d a t a  f o r   t a b u l a t i o n   d e m o
  // ----------------------------------------------------------------------------

  // Generate a dummy dataset
  RooRealVar x("x","x",0,10) ;
  RooDataSet *data = RooPolynomial("p","p",x).generate(RooArgSet(x,b0flav,tagCat),10000) ;



  // P r i n t   t a b l e s   o f   c a t e g o r y   c o n t e n t s   o f   d a t a s e t s
  // ------------------------------------------------------------------------------------------

  // Tables are equivalent of plots for categories
  Roo1DTable* btable = data->table(b0flav) ;

  // Create table for subset of events matching cut expression
  Roo1DTable* ttable = data->table(tagCat,"x>8.23") ;

  // Create table for all (tagCat x b0flav) state combinations
  Roo1DTable* bttable = data->table(RooArgSet(tagCat,b0flav)) ;

  // Retrieve number of events from table
  // Number can be non-integer if source dataset has weighed events
  Double_t nb0 = btable->get("B0") ;
  regValue(nb0,"rf404_nb0") ;

  // Retrieve fraction of events with "Lepton" tag
  Double_t fracLep = ttable->getFrac("Lepton") ;
  regValue(fracLep,"rf404_fracLep") ;


  // D e f i n i n g   r a n g e s   f o r   p l o t t i n g ,   f i t t i n g   o n   c a t e g o r i e s
  // ------------------------------------------------------------------------------------------------------

  // Define named range as comma separated list of labels
  tagCat.setRange("good","Lepton,Kaon") ;

  // Or add state names one by one
  tagCat.addToRange("soso","NetTagger-1") ;
  tagCat.addToRange("soso","NetTagger-2") ;

  // Use category range in dataset reduction specification
  RooDataSet* goodData = (RooDataSet*) data->reduce(CutRange("good")) ;
  Roo1DTable* gtable = goodData->table(tagCat) ;


  regTable(btable,"rf404_btable") ;
  regTable(ttable,"rf404_ttable") ;
  regTable(bttable,"rf404_bttable") ;
  regTable(gtable,"rf404_gtable") ;


  delete goodData ;
  delete data ;
  return kTRUE ;

  }

} ;
//////////////////////////////////////////////////////////////////////////
//
// 'DATA AND CATEGORIES' RooFit tutorial macro #405
//
// Demonstration of real-->discrete mapping functions
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCategory.h"
#include "RooThresholdCategory.h"
#include "RooBinningCategory.h"
#include "Roo1DTable.h"
#include "RooArgusBG.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooRealConstant.h"
using namespace RooFit ;


class TestBasic405 : public RooUnitTest
{
public:
  TestBasic405(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Real-to-category functions",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {


  // D e f i n e   p d f   i n   x ,   s a m p l e   d a t a s e t   i n   x
  // ------------------------------------------------------------------------


  // Define a dummy PDF in x
  RooRealVar x("x","x",0,10) ;
  RooArgusBG a("a","argus(x)",x,RooRealConstant::value(10),RooRealConstant::value(-1)) ;

  // Generate a dummy dataset
  RooDataSet *data = a.generate(x,10000) ;



  // C r e a t e   a   t h r e s h o l d   r e a l - > c a t   f u n c t i o n
  // --------------------------------------------------------------------------

  // A RooThresholdCategory is a category function that maps regions in a real-valued
  // input observable observables to state names. At construction time a 'default'
  // state name must be specified to which all values of x are mapped that are not
  // otherwise assigned
  RooThresholdCategory xRegion("xRegion","region of x",x,"Background") ;

  // Specify thresholds and state assignments one-by-one.
  // Each statement specifies that all values _below_ the given value
  // (and above any lower specified threshold) are mapped to the
  // category state with the given name
  //
  // Background | SideBand | Signal | SideBand | Background
  //           4.23       5.23     8.23       9.23
  xRegion.addThreshold(4.23,"Background") ;
  xRegion.addThreshold(5.23,"SideBand") ;
  xRegion.addThreshold(8.23,"Signal") ;
  xRegion.addThreshold(9.23,"SideBand") ;



  // U s e   t h r e s h o l d   f u n c t i o n   t o   p l o t   d a t a   r e g i o n s
  // -------------------------------------------------------------------------------------

  // Add values of threshold function to dataset so that it can be used as observable
  data->addColumn(xRegion) ;

  // Make plot of data in x
  RooPlot* xframe = x.frame(Title("Demo of threshold and binning mapping functions")) ;
  data->plotOn(xframe) ;

  // Use calculated category to select sideband data
  data->plotOn(xframe,Cut("xRegion==xRegion::SideBand"),MarkerColor(kRed),LineColor(kRed),Name("data_cut")) ;



  // C r e a t e   a   b i n n i n g    r e a l - > c a t   f u n c t i o n
  // ----------------------------------------------------------------------

  // A RooBinningCategory is a category function that maps bins of a (named) binning definition
  // in a real-valued input observable observables to state names. The state names are automatically
  // constructed from the variable name, the binning name and the bin number. If no binning name
  // is specified the default binning is mapped

  x.setBins(10,"coarse") ;
  RooBinningCategory xBins("xBins","coarse bins in x",x,"coarse") ;



  // U s e   b i n n i n g   f u n c t i o n   f o r   t a b u l a t i o n   a n d   p l o t t i n g
  // -----------------------------------------------------------------------------------------------

  // Print table of xBins state multiplicity. Note that xBins does not need to be an observable in data
  // it can be a function of observables in data as well
  Roo1DTable* xbtable = data->table(xBins) ;

  // Add values of xBins function to dataset so that it can be used as observable
  RooCategory* xb = (RooCategory*) data->addColumn(xBins) ;

  // Define range "alt" as including bins 1,3,5,7,9
  xb->setRange("alt","x_coarse_bin1,x_coarse_bin3,x_coarse_bin5,x_coarse_bin7,x_coarse_bin9") ;

  // Construct subset of data matching range "alt" but only for the first 5000 events and plot it on the fram
  RooDataSet* dataSel = (RooDataSet*) data->reduce(CutRange("alt"),EventRange(0,5000)) ;
//   dataSel->plotOn(xframe,MarkerColor(kGreen),LineColor(kGreen),Name("data_sel")) ;


  regTable(xbtable,"rf405_xbtable") ;
  regPlot(xframe,"rf405_plot1") ;

  delete data ;
  delete dataSel ;

  return kTRUE ;

  }

} ;
//////////////////////////////////////////////////////////////////////////
//
// 'DATA AND CATEGORIES' RooFit tutorial macro #406
//
// Demonstration of discrete-->discrete (invertable) functions
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPolynomial.h"
#include "RooCategory.h"
#include "RooMappedCategory.h"
#include "RooMultiCategory.h"
#include "RooSuperCategory.h"
#include "Roo1DTable.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic406 : public RooUnitTest
{
public:
  TestBasic406(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Category-to-category functions",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C o n s t r u c t  t w o   c a t e g o r i e s
  // ----------------------------------------------

  // Define a category with labels only
  RooCategory tagCat("tagCat","Tagging category") ;
  tagCat.defineType("Lepton") ;
  tagCat.defineType("Kaon") ;
  tagCat.defineType("NetTagger-1") ;
  tagCat.defineType("NetTagger-2") ;

  // Define a category with explicitly numbered states
  RooCategory b0flav("b0flav","B0 flavour eigenstate") ;
  b0flav.defineType("B0",-1) ;
  b0flav.defineType("B0bar",1) ;

  // Construct a dummy dataset with random values of tagCat and b0flav
  RooRealVar x("x","x",0,10) ;
  RooPolynomial p("p","p",x) ;
  RooDataSet* data = p.generate(RooArgSet(x,b0flav,tagCat),10000) ;



  // C r e a t e   a   c a t - > c a t   m  a p p i n g   c a t e g o r y
  // ---------------------------------------------------------------------

  // A RooMappedCategory is category->category mapping function based on string expression
  // The constructor takes an input category an a default state name to which unassigned
  // states are mapped
  RooMappedCategory tcatType("tcatType","tagCat type",tagCat,"Cut based") ;

  // Enter fully specified state mappings
  tcatType.map("Lepton","Cut based") ;
  tcatType.map("Kaon","Cut based") ;

  // Enter a wilcard expression mapping
  tcatType.map("NetTagger*","Neural Network") ;

  // Make a table of the mapped category state multiplicit in data
  Roo1DTable* mtable = data->table(tcatType) ;



  // C r e a t e   a   c a t   X   c a t   p r o d u c t   c a t e g o r y
  // ----------------------------------------------------------------------

  // A SUPER-category is 'product' of _lvalue_ categories. The state names of a super
  // category is a composite of the state labels of the input categories
  RooSuperCategory b0Xtcat("b0Xtcat","b0flav X tagCat",RooArgSet(b0flav,tagCat)) ;

  // Make a table of the product category state multiplicity in data
  Roo1DTable* stable = data->table(b0Xtcat) ;

  // Since the super category is an lvalue, assignment is explicitly possible
  b0Xtcat.setLabel("{B0bar;Lepton}") ;



  // A MULTI-category is a 'product' of any category (function). The state names of a super
  // category is a composite of the state labels of the input categories
  RooMultiCategory b0Xttype("b0Xttype","b0flav X tagType",RooArgSet(b0flav,tcatType)) ;

  // Make a table of the product category state multiplicity in data
  Roo1DTable* xtable = data->table(b0Xttype) ;

  regTable(mtable,"rf406_mtable") ;
  regTable(stable,"rf406_stable") ;
  regTable(xtable,"rf406_xtable") ;

  delete data ;

  return kTRUE ;
  }

} ;
//////////////////////////////////////////////////////////////////////////
//
// 'ORGANIZATION AND SIMULTANEOUS FITS' RooFit tutorial macro #501
//
// Using simultaneous p.d.f.s to describe simultaneous fits to multiple
// datasets
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic501 : public RooUnitTest
{
public:
  TestBasic501(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Simultaneous p.d.f. operator",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   m o d e l   f o r   p h y s i c s   s a m p l e
  // -------------------------------------------------------------

  // Create observables
  RooRealVar x("x","x",-8,8) ;

  // Construct signal pdf
  RooRealVar mean("mean","mean",0,-8,8) ;
  RooRealVar sigma("sigma","sigma",0.3,0.1,10) ;
  RooGaussian gx("gx","gx",x,mean,sigma) ;

  // Construct background pdf
  RooRealVar a0("a0","a0",-0.1,-1,1) ;
  RooRealVar a1("a1","a1",0.004,-1,1) ;
  RooChebychev px("px","px",x,RooArgSet(a0,a1)) ;

  // Construct composite pdf
  RooRealVar f("f","f",0.2,0.,1.) ;
  RooAddPdf model("model","model",RooArgList(gx,px),f) ;



  // C r e a t e   m o d e l   f o r   c o n t r o l   s a m p l e
  // --------------------------------------------------------------

  // Construct signal pdf.
  // NOTE that sigma is shared with the signal sample model
  RooRealVar mean_ctl("mean_ctl","mean_ctl",-3,-8,8) ;
  RooGaussian gx_ctl("gx_ctl","gx_ctl",x,mean_ctl,sigma) ;

  // Construct the background pdf
  RooRealVar a0_ctl("a0_ctl","a0_ctl",-0.1,-1,1) ;
  RooRealVar a1_ctl("a1_ctl","a1_ctl",0.5,-0.1,1) ;
  RooChebychev px_ctl("px_ctl","px_ctl",x,RooArgSet(a0_ctl,a1_ctl)) ;

  // Construct the composite model
  RooRealVar f_ctl("f_ctl","f_ctl",0.5,0.,1.) ;
  RooAddPdf model_ctl("model_ctl","model_ctl",RooArgList(gx_ctl,px_ctl),f_ctl) ;



  // G e n e r a t e   e v e n t s   f o r   b o t h   s a m p l e s
  // ---------------------------------------------------------------

  // Generate 1000 events in x and y from model
  RooDataSet *data = model.generate(RooArgSet(x),100) ;
  RooDataSet *data_ctl = model_ctl.generate(RooArgSet(x),2000) ;



  // C r e a t e   i n d e x   c a t e g o r y   a n d   j o i n   s a m p l e s
  // ---------------------------------------------------------------------------

  // Define category to distinguish physics and control samples events
  RooCategory sample("sample","sample") ;
  sample.defineType("physics") ;
  sample.defineType("control") ;

  // Construct combined dataset in (x,sample)
  RooDataSet combData("combData","combined data",x,Index(sample),Import("physics",*data),Import("control",*data_ctl)) ;



  // C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
  // -----------------------------------------------------------------------------------

  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;

  // Associate model with the physics state and model_ctl with the control state
  simPdf.addPdf(model,"physics") ;
  simPdf.addPdf(model_ctl,"control") ;



  // P e r f o r m   a   s i m u l t a n e o u s   f i t
  // ---------------------------------------------------

  // Perform simultaneous fit of model to data and model_ctl to data_ctl
  simPdf.fitTo(combData) ;



  // P l o t   m o d e l   s l i c e s   o n   d a t a    s l i c e s
  // ----------------------------------------------------------------

  // Make a frame for the physics sample
  RooPlot* frame1 = x.frame(Bins(30),Title("Physics sample")) ;

  // Plot all data tagged as physics sample
  combData.plotOn(frame1,Cut("sample==sample::physics")) ;

  // Plot "physics" slice of simultaneous pdf.
  // NBL You _must_ project the sample index category with data using ProjWData
  // as a RooSimultaneous makes no prediction on the shape in the index category
  // and can thus not be integrated
  simPdf.plotOn(frame1,Slice(sample,"physics"),ProjWData(sample,combData)) ;
  simPdf.plotOn(frame1,Slice(sample,"physics"),Components("px"),ProjWData(sample,combData),LineStyle(kDashed)) ;

  // The same plot for the control sample slice
  RooPlot* frame2 = x.frame(Bins(30),Title("Control sample")) ;
  combData.plotOn(frame2,Cut("sample==sample::control")) ;
  simPdf.plotOn(frame2,Slice(sample,"control"),ProjWData(sample,combData)) ;
  simPdf.plotOn(frame2,Slice(sample,"control"),Components("px_ctl"),ProjWData(sample,combData),LineStyle(kDashed)) ;


  regPlot(frame1,"rf501_plot1") ;
  regPlot(frame2,"rf501_plot2") ;

  delete data ;
  delete data_ctl ;

  return kTRUE ;

  }

} ;
//////////////////////////////////////////////////////////////////////////
//
// 'ORGANIZATION AND SIMULTANEOUS FITS' RooFit tutorial macro #501
//
// Using simultaneous p.d.f.s to describe simultaneous fits to multiple
// datasets
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooWorkspace.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooAddModel.h"
#include "RooDecay.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic599 : public RooUnitTest
{
public:
  TestBasic599(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Workspace and p.d.f. persistence",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

    if (_write) {

      RooWorkspace *w = new RooWorkspace("TestBasic11_ws") ;

      regWS(w,"Basic11_ws") ;

      // Build Gaussian PDF in X
      RooRealVar x("x","x",-10,10) ;
      RooRealVar meanx("meanx","mean of gaussian",-1) ;
      RooRealVar sigmax("sigmax","width of gaussian",3) ;
      RooGaussian gaussx("gaussx","gaussian PDF",x,meanx,sigmax) ;

      // Build Gaussian PDF in Y
      RooRealVar y("y","y",-10,10) ;
      RooRealVar meany("meany","mean of gaussian",-1) ;
      RooRealVar sigmay("sigmay","width of gaussian",3) ;
      RooGaussian gaussy("gaussy","gaussian PDF",y,meany,sigmay) ;

      // Make product of X and Y
      RooProdPdf gaussxy("gaussxy","gaussx*gaussy",RooArgSet(gaussx,gaussy)) ;

      // Make flat bkg in X and Y
      RooPolynomial flatx("flatx","flatx",x) ;
      RooPolynomial flaty("flaty","flaty",x) ;
      RooProdPdf flatxy("flatxy","flatx*flaty",RooArgSet(flatx,flaty)) ;

      // Make sum of gaussxy and flatxy
      RooRealVar frac("frac","frac",0.5,0.,1.) ;
      RooAddPdf sumxy("sumxy","sumxy",RooArgList(gaussxy,flatxy),frac) ;

      // Store p.d.f in workspace
      w->import(gaussx) ;
      w->import(gaussxy,RenameConflictNodes("set2")) ;
      w->import(sumxy,RenameConflictNodes("set3")) ;

      // Make reference plot of GaussX
      RooPlot* frame1 = x.frame() ;
      gaussx.plotOn(frame1) ;
      regPlot(frame1,"Basic11_gaussx_framex") ;

      // Make reference plots for GaussXY
      RooPlot* frame2 = x.frame() ;
      gaussxy.plotOn(frame2) ;
      regPlot(frame2,"Basic11_gaussxy_framex") ;

      RooPlot* frame3 = y.frame() ;
      gaussxy.plotOn(frame3) ;
      regPlot(frame3,"Basic11_gaussxy_framey") ;

      // Make reference plots for SumXY
      RooPlot* frame4 = x.frame() ;
      sumxy.plotOn(frame4) ;
      regPlot(frame4,"Basic11_sumxy_framex") ;

      RooPlot* frame5 = y.frame() ;
      sumxy.plotOn(frame5) ;
      regPlot(frame5,"Basic11_sumxy_framey") ;

      // Analytically convolved p.d.f.s

      // Build a simple decay PDF
      RooRealVar dt("dt","dt",-20,20) ;
      RooRealVar tau("tau","tau",1.548) ;

      // Build a gaussian resolution model
      RooRealVar bias1("bias1","bias1",0) ;
      RooRealVar sigma1("sigma1","sigma1",1) ;
      RooGaussModel gm1("gm1","gauss model 1",dt,bias1,sigma1) ;

      // Construct a decay PDF, smeared with single gaussian resolution model
      RooDecay decay_gm1("decay_gm1","decay",dt,tau,gm1,RooDecay::DoubleSided) ;

      // Build another gaussian resolution model
      RooRealVar bias2("bias2","bias2",0) ;
      RooRealVar sigma2("sigma2","sigma2",5) ;
      RooGaussModel gm2("gm2","gauss model 2",dt,bias2,sigma2) ;

      // Build a composite resolution model
      RooRealVar gm1frac("gm1frac","fraction of gm1",0.5) ;
      RooAddModel gmsum("gmsum","sum of gm1 and gm2",RooArgList(gm1,gm2),gm1frac) ;

      // Construct a decay PDF, smeared with double gaussian resolution model
      RooDecay decay_gmsum("decay_gmsum","decay",dt,tau,gmsum,RooDecay::DoubleSided) ;

      w->import(decay_gm1) ;
      w->import(decay_gmsum,RenameConflictNodes("set3")) ;

      RooPlot* frame6 = dt.frame() ;
      decay_gm1.plotOn(frame6) ;
      regPlot(frame6,"Basic11_decay_gm1_framedt") ;

      RooPlot* frame7 = dt.frame() ;
      decay_gmsum.plotOn(frame7) ;
      regPlot(frame7,"Basic11_decay_gmsum_framedt") ;

      // Construct simultaneous p.d.f
      RooCategory cat("cat","cat") ;
      cat.defineType("A") ;
      cat.defineType("B") ;
      RooSimultaneous sim("sim","sim",cat) ;
      sim.addPdf(gaussxy,"A") ;
      sim.addPdf(flatxy,"B") ;

      w->import(sim,RenameConflictNodes("set4")) ;

      // Make plot with dummy dataset for index projection
      RooDataHist dh("dh","dh",cat) ;
      cat.setLabel("A") ;
      dh.add(cat) ;
      cat.setLabel("B") ;
      dh.add(cat) ;

      RooPlot* frame8 = x.frame() ;
      sim.plotOn(frame8,ProjWData(cat,dh),Project(cat)) ;

      regPlot(frame8,"Basic11_sim_framex") ;


    } else {

      RooWorkspace* w = getWS("Basic11_ws") ;
      if (!w) return kFALSE ;

      // Retrieve p.d.f from workspace
      RooAbsPdf* gaussx = w->pdf("gaussx") ;

      // Make test plot and offer for comparison against ref plot
      RooPlot* frame1 = w->var("x")->frame() ;
      gaussx->plotOn(frame1) ;
      regPlot(frame1,"Basic11_gaussx_framex") ;

      // Retrieve p.d.f from workspace
      RooAbsPdf* gaussxy = w->pdf("gaussxy") ;

      // Make test plot and offer for comparison against ref plot
      RooPlot* frame2 = w->var("x")->frame() ;
      gaussxy->plotOn(frame2) ;
      regPlot(frame2,"Basic11_gaussxy_framex") ;

      // Make test plot and offer for comparison against ref plot
      RooPlot* frame3 = w->var("y")->frame() ;
      gaussxy->plotOn(frame3) ;
      regPlot(frame3,"Basic11_gaussxy_framey") ;

      // Retrieve p.d.f from workspace
      RooAbsPdf* sumxy = w->pdf("sumxy") ;

      // Make test plot and offer for comparison against ref plot
      RooPlot* frame4 = w->var("x")->frame() ;
      sumxy->plotOn(frame4) ;
      regPlot(frame4,"Basic11_sumxy_framex") ;

      // Make test plot and offer for comparison against ref plot
      RooPlot* frame5 = w->var("y")->frame() ;
      sumxy->plotOn(frame5) ;
      regPlot(frame5,"Basic11_sumxy_framey") ;

      // Retrieve p.d.f from workspace
      RooAbsPdf* decay_gm1 = w->pdf("decay_gm1") ;

      // Make test plot and offer for comparison against ref plot
      RooPlot* frame6 = w->var("dt")->frame() ;
      decay_gm1->plotOn(frame6) ;
      regPlot(frame6,"Basic11_decay_gm1_framedt") ;

      // Retrieve p.d.f from workspace
      RooAbsPdf* decay_gmsum = w->pdf("decay_gmsum") ;

      // Make test plot and offer for comparison against ref plot
      RooPlot* frame7 = w->var("dt")->frame() ;
      decay_gmsum->plotOn(frame7) ;
      regPlot(frame7,"Basic11_decay_gmsum_framedt") ;

      // Retrieve p.d.f. from workspace
      RooAbsPdf* sim = w->pdf("sim") ;
      RooCategory* cat = w->cat("cat") ;

      // Make plot with dummy dataset for index projection
      RooPlot* frame8 = w->var("x")->frame() ;

      RooDataHist dh("dh","dh",*cat) ;
      cat->setLabel("A") ;
      dh.add(*cat) ;
      cat->setLabel("B") ;
      dh.add(*cat) ;

      sim->plotOn(frame8,ProjWData(*cat,dh),Project(*cat)) ;

      regPlot(frame8,"Basic11_sim_framex") ;

    }

    // "Workspace persistence"
    return kTRUE ;
  }
} ;
/////////////////////////////////////////////////////////////////////////
//
// 'LIKELIHOOD AND MINIMIZATION' RooFit tutorial macro #601
//
// Interactive minimization with MINUIT
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
using namespace RooFit ;


class TestBasic601 : public RooUnitTest
{
public:
  TestBasic601(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Interactive Minuit",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // S e t u p   p d f   a n d   l i k e l i h o o d
  // -----------------------------------------------

  // Observable
  RooRealVar x("x","x",-20,20) ;

  // Model (intentional strong correlations)
  RooRealVar mean("mean","mean of g1 and g2",0) ;
  RooRealVar sigma_g1("sigma_g1","width of g1",3) ;
  RooGaussian g1("g1","g1",x,mean,sigma_g1) ;

  RooRealVar sigma_g2("sigma_g2","width of g2",4,3.0,6.0) ;
  RooGaussian g2("g2","g2",x,mean,sigma_g2) ;

  RooRealVar frac("frac","frac",0.5,0.0,1.0) ;
  RooAddPdf model("model","model",RooArgList(g1,g2),frac) ;

  // Generate 1000 events
  RooDataSet* data = model.generate(x,1000) ;

  // Construct unbinned likelihood
  RooNLLVar nll("nll","nll",model,*data) ;


  // I n t e r a c t i v e   m i n i m i z a t i o n ,   e r r o r   a n a l y s i s
  // -------------------------------------------------------------------------------

  // Create MINUIT interface object
  RooMinuit m(nll) ;

  // Call MIGRAD to minimize the likelihood
  m.migrad() ;

  // Run HESSE to calculate errors from d2L/dp2
  m.hesse() ;

  // Run MINOS on sigma_g2 parameter only
  m.minos(sigma_g2) ;


  // S a v i n g   r e s u l t s ,   c o n t o u r   p l o t s
  // ---------------------------------------------------------

  // Save a snapshot of the fit result. This object contains the initial
  // fit parameters, the final fit parameters, the complete correlation
  // matrix, the EDM, the minimized FCN , the last MINUIT status code and
  // the number of times the RooFit function object has indicated evaluation
  // problems (e.g. zero probabilities during likelihood evaluation)
  RooFitResult* r = m.save() ;


  // C h a n g e   p a r a m e t e r   v a l u e s ,   f l o a t i n g
  // -----------------------------------------------------------------

  // At any moment you can manually change the value of a (constant)
  // parameter
  mean = 0.3 ;

  // Rerun MIGRAD,HESSE
  m.migrad() ;
  m.hesse() ;

  // Now fix sigma_g2
  sigma_g2.setConstant(kTRUE) ;

  // Rerun MIGRAD,HESSE
  m.migrad() ;
  m.hesse() ;

  RooFitResult* r2 = m.save() ;

  regResult(r,"rf601_r") ;
  regResult(r2,"rf601_r2") ;

  delete data ;

  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'LIKELIHOOD AND MINIMIZATION' RooFit tutorial macro #602
//
// Setting up a binning chi^2 fit
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic602 : public RooUnitTest
{
public:
  TestBasic602(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Chi2 minimization",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // S e t u p   m o d e l
  // ---------------------

  // Declare observable x
  RooRealVar x("x","x",0,10) ;

  // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
  RooRealVar mean("mean","mean of gaussians",5) ;
  RooRealVar sigma1("sigma1","width of gaussians",0.5) ;
  RooRealVar sigma2("sigma2","width of gaussians",1) ;

  RooGaussian sig1("sig1","Signal component 1",x,mean,sigma1) ;
  RooGaussian sig2("sig2","Signal component 2",x,mean,sigma2) ;

  // Build Chebychev polynomial p.d.f.
  RooRealVar a0("a0","a0",0.5,0.,1.) ;
  RooRealVar a1("a1","a1",-0.2,-1.,1.) ;
  RooChebychev bkg("bkg","Background",x,RooArgSet(a0,a1)) ;

  // Sum the signal components into a composite signal p.d.f.
  RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
  RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),sig1frac) ;

  // Sum the composite signal and background
  RooRealVar bkgfrac("bkgfrac","fraction of background",0.5,0.,1.) ;
  RooAddPdf  model("model","g1+g2+a",RooArgList(bkg,sig),bkgfrac) ;


  // C r e a t e   b i n n e d   d a t a s e t
  // -----------------------------------------

  RooDataSet* d = model.generate(x,10000) ;
  RooDataHist* dh = d->binnedClone() ;


  // Construct a chi^2 of the data and the model,
  // which is the input probability density scaled
  // by the number of events in the dataset
  RooChi2Var chi2("chi2","chi2",model,*dh) ;

  // Use RooMinuit interface to minimize chi^2
  RooMinuit m(chi2) ;
  m.migrad() ;
  m.hesse() ;

  RooFitResult* r = m.save() ;

  regResult(r,"rf602_r") ;

  delete d ;
  delete dh ;

  return kTRUE ;
  }
} ;
/////////////////////////////////////////////////////////////////////////
//
// 'LIKELIHOOD AND MINIMIZATION' RooFit tutorial macro #604
//
// Fitting with constraints
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
using namespace RooFit ;


class TestBasic604 : public RooUnitTest
{
public:
  TestBasic604(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Auxiliary observable constraints",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   m o d e l  a n d   d a t a s e t
  // ----------------------------------------------

  // Construct a Gaussian p.d.f
  RooRealVar x("x","x",-10,10) ;

  RooRealVar m("m","m",0,-10,10) ;
  RooRealVar s("s","s",2,0.1,10) ;
  RooGaussian gauss("gauss","gauss(x,m,s)",x,m,s) ;

  // Construct a flat p.d.f (polynomial of 0th order)
  RooPolynomial poly("poly","poly(x)",x) ;

  // Construct model = f*gauss + (1-f)*poly
  RooRealVar f("f","f",0.5,0.,1.) ;
  RooAddPdf model("model","model",RooArgSet(gauss,poly),f) ;

  // Generate small dataset for use in fitting below
  RooDataSet* d = model.generate(x,50) ;



  // C r e a t e   c o n s t r a i n t   p d f
  // -----------------------------------------

  // Construct Gaussian constraint p.d.f on parameter f at 0.8 with resolution of 0.1
  RooGaussian fconstraint("fconstraint","fconstraint",f,RooConst(0.8),RooConst(0.1)) ;



  // M E T H O D   1   -   A d d   i n t e r n a l   c o n s t r a i n t   t o   m o d e l
  // -------------------------------------------------------------------------------------

  // Multiply constraint term with regular p.d.f using RooProdPdf
  // Specify in fitTo() that internal constraints on parameter f should be used

  // Multiply constraint with p.d.f
  RooProdPdf modelc("modelc","model with constraint",RooArgSet(model,fconstraint)) ;

  // Fit modelc without use of constraint term
  RooFitResult* r1 = modelc.fitTo(*d,Save()) ;

  // Fit modelc with constraint term on parameter f
  RooFitResult* r2 = modelc.fitTo(*d,Constrain(f),Save()) ;



  // M E T H O D   2   -     S p e c i f y   e x t e r n a l   c o n s t r a i n t   w h e n   f i t t i n g
  // -------------------------------------------------------------------------------------------------------

  // Construct another Gaussian constraint p.d.f on parameter f at 0.8 with resolution of 0.1
  RooGaussian fconstext("fconstext","fconstext",f,RooConst(0.2),RooConst(0.1)) ;

  // Fit with external constraint
  RooFitResult* r3 = model.fitTo(*d,ExternalConstraints(fconstext),Save()) ;


  regResult(r1,"rf604_r1") ;
  regResult(r2,"rf604_r2") ;
  regResult(r3,"rf604_r3") ;

  delete d ;

  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'LIKELIHOOD AND MINIMIZATION' RooFit tutorial macro #605
//
// Working with the profile likelihood estimator
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooNLLVar.h"
#include "RooProfileLL.h"
#include "RooMinuit.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic605 : public RooUnitTest
{
public:
  TestBasic605(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Profile Likelihood operator",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   m o d e l   a n d   d a t a s e t
  // -----------------------------------------------

  // Observable
  RooRealVar x("x","x",-20,20) ;

  // Model (intentional strong correlations)
  RooRealVar mean("mean","mean of g1 and g2",0,-10,10) ;
  RooRealVar sigma_g1("sigma_g1","width of g1",3) ;
  RooGaussian g1("g1","g1",x,mean,sigma_g1) ;

  RooRealVar sigma_g2("sigma_g2","width of g2",4,3.0,6.0) ;
  RooGaussian g2("g2","g2",x,mean,sigma_g2) ;

  RooRealVar frac("frac","frac",0.5,0.0,1.0) ;
  RooAddPdf model("model","model",RooArgList(g1,g2),frac) ;

  // Generate 1000 events
  RooDataSet* data = model.generate(x,1000) ;



  // C o n s t r u c t   p l a i n   l i k e l i h o o d
  // ---------------------------------------------------

  // Construct unbinned likelihood
  RooNLLVar nll("nll","nll",model,*data) ;

  // Minimize likelihood w.r.t all parameters before making plots
  RooMinuit(nll).migrad() ;

  // Plot likelihood scan frac
  RooPlot* frame1 = frac.frame(Bins(10),Range(0.01,0.95),Title("LL and profileLL in frac")) ;
  nll.plotOn(frame1,ShiftToZero()) ;

  // Plot likelihood scan in sigma_g2
  RooPlot* frame2 = sigma_g2.frame(Bins(10),Range(3.3,5.0),Title("LL and profileLL in sigma_g2")) ;
  nll.plotOn(frame2,ShiftToZero()) ;



  // C o n s t r u c t   p r o f i l e   l i k e l i h o o d   i n   f r a c
  // -----------------------------------------------------------------------

  // The profile likelihood estimator on nll for frac will minimize nll w.r.t
  // all floating parameters except frac for each evaluation
  RooProfileLL pll_frac("pll_frac","pll_frac",nll,frac) ;

  // Plot the profile likelihood in frac
  pll_frac.plotOn(frame1,LineColor(kRed)) ;

  // Adjust frame maximum for visual clarity
  frame1->SetMinimum(0) ;
  frame1->SetMaximum(3) ;



  // C o n s t r u c t   p r o f i l e   l i k e l i h o o d   i n   s i g m a _ g 2
  // -------------------------------------------------------------------------------

  // The profile likelihood estimator on nll for sigma_g2 will minimize nll
  // w.r.t all floating parameters except sigma_g2 for each evaluation
  RooProfileLL pll_sigmag2("pll_sigmag2","pll_sigmag2",nll,sigma_g2) ;

  // Plot the profile likelihood in sigma_g2
  pll_sigmag2.plotOn(frame2,LineColor(kRed)) ;

  // Adjust frame maximum for visual clarity
  frame2->SetMinimum(0) ;
  frame2->SetMaximum(3) ;


  regPlot(frame1,"rf605_plot1") ;
  regPlot(frame2,"rf605_plot2") ;

  delete data ;

  return kTRUE ;
  }
} ;

//////////////////////////////////////////////////////////////////////////
//
// 'LIKELIHOOD AND MINIMIZATION' RooFit tutorial macro #606
//
// Understanding and customizing error handling in likelihood evaluations
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgusBG.h"
#include "RooNLLVar.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


class TestBasic606 : public RooUnitTest
{
public:
  TestBasic606(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("NLL error handling",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   m o d e l  a n d   d a t a s e t
  // ----------------------------------------------

  // Observable
  RooRealVar m("m","m",5.20,5.30) ;

  // Parameters
  RooRealVar m0("m0","m0",5.291,5.20,5.30) ;
  RooRealVar k("k","k",-30,-50,-10) ;

  // Pdf
  RooArgusBG argus("argus","argus",m,m0,k) ;

  // Sample 1000 events in m from argus
  RooDataSet* data = argus.generate(m,1000) ;



  // P l o t   m o d e l   a n d   d a t a
  // --------------------------------------

  RooPlot* frame1 = m.frame(Bins(40),Title("Argus model and data")) ;
  data->plotOn(frame1) ;
  argus.plotOn(frame1) ;



  // F i t   m o d e l   t o   d a t a
  // ---------------------------------

  argus.fitTo(*data,PrintEvalErrors(10),Warnings(kFALSE)) ;
  m0.setError(0.1) ;
  argus.fitTo(*data,PrintEvalErrors(0),EvalErrorWall(kFALSE),Warnings(kFALSE)) ;



  // P l o t   l i k e l i h o o d   a s   f u n c t i o n   o f   m 0
  // ------------------------------------------------------------------

  // Construct likelihood function of model and data
  RooNLLVar nll("nll","nll",argus,*data) ;

  // Plot likelihood in m0 in range that includes problematic values
  // In this configuration no messages are printed for likelihood evaluation errors,
  // but if an likelihood value evaluates with error, the corresponding value
  // on the curve will be set to the value given in EvalErrorValue().

  RooPlot* frame2 = m0.frame(Range(5.288,5.293),Title("-log(L) scan vs m0, problematic regions masked")) ;
  nll.plotOn(frame2,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;


  regPlot(frame1,"rf606_plot1") ;
  regPlot(frame2,"rf606_plot3") ; // 3 is the reference of the plot

  delete data ;
  return kTRUE ;

  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'LIKELIHOOD AND MINIMIZATION' RooFit tutorial macro #607
//
// Demonstration of options of the RooFitResult class
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"

using namespace RooFit ;


class TestBasic607 : public RooUnitTest
{
public:
  TestBasic607(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Fit Result functionality",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   p d f ,   d a t a
  // --------------------------------

  // Declare observable x
  RooRealVar x("x","x",0,10) ;

  // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
  RooRealVar mean("mean","mean of gaussians",5,-10,10) ;
  RooRealVar sigma1("sigma1","width of gaussians",0.5,0.1,10) ;
  RooRealVar sigma2("sigma2","width of gaussians",1,0.1,10) ;

  RooGaussian sig1("sig1","Signal component 1",x,mean,sigma1) ;
  RooGaussian sig2("sig2","Signal component 2",x,mean,sigma2) ;

  // Build Chebychev polynomial p.d.f.
  RooRealVar a0("a0","a0",0.5,0.,1.) ;
  RooRealVar a1("a1","a1",-0.2,-1.,1.) ;
  RooChebychev bkg("bkg","Background",x,RooArgSet(a0,a1)) ;

  // Sum the signal components into a composite signal p.d.f.
  RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
  RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),sig1frac) ;

  // Sum the composite signal and background
  RooRealVar bkgfrac("bkgfrac","fraction of background",0.5,0.,1.) ;
  RooAddPdf  model("model","g1+g2+a",RooArgList(bkg,sig),bkgfrac) ;

  // Generate 1000 events
  RooDataSet* data = model.generate(x,1000) ;



  // F i t   p d f   t o   d a t a ,   s a v e   f i t r e s u l t
  // -------------------------------------------------------------

  // Perform fit and save result
  RooFitResult* r = model.fitTo(*data,Save()) ;


  // V i s u a l i z e   c o r r e l a t i o n   m a t r i x
  // -------------------------------------------------------

  // Construct 2D color plot of correlation matrix
  gStyle->SetOptStat(0) ;
  gStyle->SetPalette(1) ;
  TH2* hcorr = r->correlationHist() ;


  // Sample dataset with parameter values according to distribution
  // of covariance matrix of fit result
  RooDataSet randPars("randPars","randPars",r->floatParsFinal()) ;
  for (Int_t i=0 ; i<10000 ; i++) {
    randPars.add(r->randomizePars()) ;
  }

  // make histogram of 2D distribution in sigma1 vs sig1frac
  TH1* hhrand = randPars.createHistogram("hhrand",sigma1,Binning(35,0.25,0.65),YVar(sig1frac,Binning(35,0.3,1.1))) ;

  regTH(hcorr,"rf607_hcorr") ;
  regTH(hhrand,"rf607_hhand") ;

  delete data ;
  delete r ;

  return kTRUE ;

  }
} ;


//////////////////////////////////////////////////////////////////////////
//
// 'LIKELIHOOD AND MINIMIZATION' RooFit tutorial macro #609
//
// Setting up a chi^2 fit to an unbinned dataset with X,Y,err(Y)
// values (and optionally err(X) values)
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPolyVar.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TRandom.h"

using namespace RooFit ;

class TestBasic609 : public RooUnitTest
{
public:
  TestBasic609(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Chi^2 fit to X-Y dataset",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   d a t a s e t   w i t h   X   a n d   Y   v a l u e s
  // -------------------------------------------------------------------

  // Make weighted XY dataset with asymmetric errors stored
  // The StoreError() argument is essential as it makes
  // the dataset store the error in addition to the values
  // of the observables. If errors on one or more observables
  // are asymmetric, one can store the asymmetric error
  // using the StoreAsymError() argument

  RooRealVar x("x","x",-11,11) ;
  RooRealVar y("y","y",-10,200) ;
  RooDataSet dxy("dxy","dxy",RooArgSet(x,y),StoreError(RooArgSet(x,y))) ;

  // Fill an example dataset with X,err(X),Y,err(Y) values
  for (int i=0 ; i<=10 ; i++) {

    // Set X value and error
    x = -10 + 2*i;
    x.setError( i<5 ? 0.5/1. : 1.0/1. ) ;

    // Set Y value and error
    y = x.getVal() * x.getVal() + 4*fabs(RooRandom::randomGenerator()->Gaus()) ;
    y.setError(sqrt(y.getVal())) ;

    dxy.add(RooArgSet(x,y)) ;
  }



  // P e r f o r m   c h i 2   f i t   t o   X + / - d x   a n d   Y + / - d Y   v a l u e s
  // ---------------------------------------------------------------------------------------

  // Make fit function
  RooRealVar a("a","a",0.0,-10,10) ;
  RooRealVar b("b","b",0.0,-100,100) ;
  RooPolyVar f("f","f",x,RooArgList(b,a,RooConst(1))) ;

  // Plot dataset in X-Y interpretation
  RooPlot* frame = x.frame(Title("Chi^2 fit of function set of (X#pmdX,Y#pmdY) values")) ;
  dxy.plotOnXY(frame,YVar(y)) ;

  // Fit chi^2 using X and Y errors
  f.chi2FitTo(dxy,YVar(y)) ;

  // Overlay fitted function
  f.plotOn(frame) ;

  // Alternative: fit chi^2 integrating f(x) over ranges defined by X errors, rather
  // than taking point at center of bin
  f.chi2FitTo(dxy,YVar(y),Integrate(kTRUE)) ;

  // Overlay alternate fit result
  f.plotOn(frame,LineStyle(kDashed),LineColor(kRed),Name("alternate")) ;

  regPlot(frame,"rf609_frame") ;

  return kTRUE ;
  }
} ;



//////////////////////////////////////////////////////////////////////////
//
// 'SPECIAL PDFS' RooFit tutorial macro #701
//
// Unbinned maximum likelihood fit of an efficiency eff(x) function to
// a dataset D(x,cut), where cut is a category encoding a selection, of which
// the efficiency as function of x should be described by eff(x)
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooEfficiency.h"
#include "RooPolynomial.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


// Elementary operations on a gaussian PDF
class TestBasic701 : public RooUnitTest
{
public:
  TestBasic701(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Efficiency operator p.d.f. 1D",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C o n s t r u c t   e f f i c i e n c y   f u n c t i o n   e ( x )
  // -------------------------------------------------------------------

  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
  RooRealVar x("x","x",-10,10) ;

  // Efficiency function eff(x;a,b)
  RooRealVar a("a","a",0.4,0,1) ;
  RooRealVar b("b","b",5) ;
  RooRealVar c("c","c",-1,-10,10) ;
  RooFormulaVar effFunc("effFunc","(1-a)+a*cos((x-c)/b)",RooArgList(a,b,c,x)) ;



  // C o n s t r u c t   c o n d i t i o n a l    e f f i c i e n c y   p d f   E ( c u t | x )
  // ------------------------------------------------------------------------------------------

  // Acceptance state cut (1 or 0)
  RooCategory cut("cut","cutr") ;
  cut.defineType("accept",1) ;
  cut.defineType("reject",0) ;

  // Construct efficiency p.d.f eff(cut|x)
  RooEfficiency effPdf("effPdf","effPdf",effFunc,cut,"accept") ;



  // G e n e r a t e   d a t a   ( x ,   c u t )   f r o m   a   t o y   m o d e l
  // -----------------------------------------------------------------------------

  // Construct global shape p.d.f shape(x) and product model(x,cut) = eff(cut|x)*shape(x)
  // (These are _only_ needed to generate some toy MC here to be used later)
  RooPolynomial shapePdf("shapePdf","shapePdf",x,RooConst(-0.095)) ;
  RooProdPdf model("model","model",shapePdf,Conditional(effPdf,cut)) ;

  // Generate some toy data from model
  RooDataSet* data = model.generate(RooArgSet(x,cut),10000) ;



  // F i t   c o n d i t i o n a l   e f f i c i e n c y   p d f   t o   d a t a
  // --------------------------------------------------------------------------

  // Fit conditional efficiency p.d.f to data
  effPdf.fitTo(*data,ConditionalObservables(x)) ;



  // P l o t   f i t t e d ,   d a t a   e f f i c i e n c y
  // --------------------------------------------------------

  // Plot distribution of all events and accepted fraction of events on frame
  RooPlot* frame1 = x.frame(Bins(20),Title("Data (all, accepted)")) ;
  data->plotOn(frame1) ;
  data->plotOn(frame1,Cut("cut==cut::accept"),MarkerColor(kRed),LineColor(kRed)) ;

  // Plot accept/reject efficiency on data overlay fitted efficiency curve
  RooPlot* frame2 = x.frame(Bins(20),Title("Fitted efficiency")) ;
  data->plotOn(frame2,Efficiency(cut)) ; // needs ROOT version >= 5.21
  effFunc.plotOn(frame2,LineColor(kRed)) ;


  regPlot(frame1,"rf701_plot1") ;
  regPlot(frame2,"rf701_plot2") ;


  delete data ;

  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'SPECIAL PDFS' RooFit tutorial macro #702
//
// Unbinned maximum likelihood fit of an efficiency eff(x) function to
// a dataset D(x,cut), where cut is a category encoding a selection whose
// efficiency as function of x should be described by eff(x)
//
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCategory.h"
#include "RooEfficiency.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooFormulaVar.h"
#include "TCanvas.h"
#include "TH1.h"
#include "RooPlot.h"
using namespace RooFit ;


// Elementary operations on a gaussian PDF
class TestBasic702 : public RooUnitTest
{
public:
  TestBasic702(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Efficiency operator p.d.f. 2D",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  Bool_t flat=kFALSE ;

  // C o n s t r u c t   e f f i c i e n c y   f u n c t i o n   e ( x , y )
  // -----------------------------------------------------------------------

  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
  RooRealVar x("x","x",-10,10) ;
  RooRealVar y("y","y",-10,10) ;

  // Efficiency function eff(x;a,b)
  RooRealVar ax("ax","ay",0.6,0,1) ;
  RooRealVar bx("bx","by",5) ;
  RooRealVar cx("cx","cy",-1,-10,10) ;

  RooRealVar ay("ay","ay",0.2,0,1) ;
  RooRealVar by("by","by",5) ;
  RooRealVar cy("cy","cy",-1,-10,10) ;

  RooFormulaVar effFunc("effFunc","((1-ax)+ax*cos((x-cx)/bx))*((1-ay)+ay*cos((y-cy)/by))",RooArgList(ax,bx,cx,x,ay,by,cy,y)) ;

  // Acceptance state cut (1 or 0)
  RooCategory cut("cut","cutr") ;
  cut.defineType("accept",1) ;
  cut.defineType("reject",0) ;



  // C o n s t r u c t   c o n d i t i o n a l    e f f i c i e n c y   p d f   E ( c u t | x , y )
  // ---------------------------------------------------------------------------------------------

  // Construct efficiency p.d.f eff(cut|x)
  RooEfficiency effPdf("effPdf","effPdf",effFunc,cut,"accept") ;



  // G e n e r a t e   d a t a   ( x , y , c u t )   f r o m   a   t o y   m o d e l
  // -------------------------------------------------------------------------------

  // Construct global shape p.d.f shape(x) and product model(x,cut) = eff(cut|x)*shape(x)
  // (These are _only_ needed to generate some toy MC here to be used later)
  RooPolynomial shapePdfX("shapePdfX","shapePdfX",x,RooConst(flat?0:-0.095)) ;
  RooPolynomial shapePdfY("shapePdfY","shapePdfY",y,RooConst(flat?0:+0.095)) ;
  RooProdPdf shapePdf("shapePdf","shapePdf",RooArgSet(shapePdfX,shapePdfY)) ;
  RooProdPdf model("model","model",shapePdf,Conditional(effPdf,cut)) ;

  // Generate some toy data from model
  RooDataSet* data = model.generate(RooArgSet(x,y,cut),10000) ;



  // F i t   c o n d i t i o n a l   e f f i c i e n c y   p d f   t o   d a t a
  // --------------------------------------------------------------------------

  // Fit conditional efficiency p.d.f to data
  effPdf.fitTo(*data,ConditionalObservables(RooArgSet(x,y))) ;



  // P l o t   f i t t e d ,   d a t a   e f f i c i e n c y
  // --------------------------------------------------------

  // Make 2D histograms of all data, selected data and efficiency function
  TH1* hh_data_all = data->createHistogram("hh_data_all",x,Binning(8),YVar(y,Binning(8))) ;
  TH1* hh_data_sel = data->createHistogram("hh_data_sel",x,Binning(8),YVar(y,Binning(8)),Cut("cut==cut::accept")) ;
  TH1* hh_eff      = effFunc.createHistogram("hh_eff",x,Binning(50),YVar(y,Binning(50))) ;

  // Some adjustsment for good visualization
  hh_data_all->SetMinimum(0) ;
  hh_data_sel->SetMinimum(0) ;
  hh_eff->SetMinimum(0) ;
  hh_eff->SetLineColor(kBlue) ;

  regTH(hh_data_all,"rf702_hh_data_all") ;
  regTH(hh_data_sel,"rf702_hh_data_sel") ;
  regTH(hh_eff,"rf702_hh_eff") ;

  delete data ;

  return kTRUE;

  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'SPECIAL PDFS' RooFit tutorial macro #703
//
// Using a product of an (acceptance) efficiency and a p.d.f as p.d.f.
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooEffProd.h"
#include "RooFormulaVar.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


// Elementary operations on a gaussian PDF
class TestBasic703 : public RooUnitTest
{
public:
  TestBasic703(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Efficiency product operator p.d.f",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // D e f i n e   o b s e r v a b l e s   a n d   d e c a y   p d f
  // ---------------------------------------------------------------

  // Declare observables
  RooRealVar t("t","t",0,5) ;

  // Make pdf
  RooRealVar tau("tau","tau",-1.54,-4,-0.1) ;
  RooExponential model("model","model",t,tau) ;



  // D e f i n e   e f f i c i e n c y   f u n c t i o n
  // ---------------------------------------------------

  // Use error function to simulate turn-on slope
  RooFormulaVar eff("eff","0.5*(TMath::Erf((t-1)/0.5)+1)",t) ;



  // D e f i n e   d e c a y   p d f   w i t h   e f f i c i e n c y
  // ---------------------------------------------------------------

  // Multiply pdf(t) with efficiency in t
  RooEffProd modelEff("modelEff","model with efficiency",model,eff) ;



  // P l o t   e f f i c i e n c y ,   p d f
  // ----------------------------------------

  RooPlot* frame1 = t.frame(Title("Efficiency")) ;
  eff.plotOn(frame1,LineColor(kRed)) ;

  RooPlot* frame2 = t.frame(Title("Pdf with and without efficiency")) ;

  model.plotOn(frame2,LineStyle(kDashed)) ;
  modelEff.plotOn(frame2) ;



  // G e n e r a t e   t o y   d a t a ,   f i t   m o d e l E f f   t o   d a t a
  // ------------------------------------------------------------------------------

  // Generate events. If the input pdf has an internal generator, the internal generator
  // is used and an accept/reject sampling on the efficiency is applied.
  RooDataSet* data = modelEff.generate(t,10000) ;

  // Fit pdf. The normalization integral is calculated numerically.
  modelEff.fitTo(*data) ;

  // Plot generated data and overlay fitted pdf
  RooPlot* frame3 = t.frame(Title("Fitted pdf with efficiency")) ;
  data->plotOn(frame3) ;
  modelEff.plotOn(frame3) ;


  regPlot(frame1,"rf703_plot1") ;
  regPlot(frame2,"rf703_plot2") ;
  regPlot(frame3,"rf703_plot3") ;


  delete data ;
  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'SPECIAL PDFS' RooFit tutorial macro #704
//
// Using a p.d.f defined by a sum of real-valued amplitude components
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooTruthModel.h"
#include "RooFormulaVar.h"
#include "RooRealSumPdf.h"
#include "RooPolyVar.h"
#include "RooProduct.h"
#include "TH1.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


// Elementary operations on a gaussian PDF
class TestBasic704 : public RooUnitTest
{
public:
  TestBasic704(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Amplitude sum operator p.d.f",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // S e t u p   2 D   a m p l i t u d e   f u n c t i o n s
  // -------------------------------------------------------

  // Observables
  RooRealVar t("t","time",-1.,15.);
  RooRealVar cosa("cosa","cos(alpha)",-1.,1.);

  // Use RooTruthModel to obtain compiled implementation of sinh/cosh modulated decay functions
  RooRealVar tau("tau","#tau",1.5);
  RooRealVar deltaGamma("deltaGamma","deltaGamma", 0.3);
  RooTruthModel tm("tm","tm",t) ;
  RooFormulaVar coshGBasis("coshGBasis","exp(-@0/ @1)*cosh(@0*@2/2)",RooArgList(t,tau,deltaGamma));
  RooFormulaVar sinhGBasis("sinhGBasis","exp(-@0/ @1)*sinh(@0*@2/2)",RooArgList(t,tau,deltaGamma));
  RooAbsReal* coshGConv = tm.convolution(&coshGBasis,&t);
  RooAbsReal* sinhGConv = tm.convolution(&sinhGBasis,&t);

  // Construct polynomial amplitudes in cos(a)
  RooPolyVar poly1("poly1","poly1",cosa,RooArgList(RooConst(0.5),RooConst(0.2),RooConst(0.2)),0);
  RooPolyVar poly2("poly2","poly2",cosa,RooArgList(RooConst(1),RooConst(-0.2),RooConst(3)),0);

  // Construct 2D amplitude as uncorrelated product of amp(t)*amp(cosa)
  RooProduct  ampl1("ampl1","amplitude 1",RooArgSet(poly1,*coshGConv));
  RooProduct  ampl2("ampl2","amplitude 2",RooArgSet(poly2,*sinhGConv));



  // C o n s t r u c t   a m p l i t u d e   s u m   p d f
  // -----------------------------------------------------

  // Amplitude strengths
  RooRealVar f1("f1","f1",1,0,2) ;
  RooRealVar f2("f2","f2",0.5,0,2) ;

  // Construct pdf
  RooRealSumPdf pdf("pdf","pdf",RooArgList(ampl1,ampl2),RooArgList(f1,f2)) ;

  // Generate some toy data from pdf
  RooDataSet* data = pdf.generate(RooArgSet(t,cosa),10000);

  // Fit pdf to toy data with only amplitude strength floating
  pdf.fitTo(*data) ;



  // P l o t   a m p l i t u d e   s u m   p d f
  // -------------------------------------------

  // Make 2D plots of amplitudes
  TH1* hh_cos = ampl1.createHistogram("hh_cos",t,Binning(50),YVar(cosa,Binning(50))) ;
  TH1* hh_sin = ampl2.createHistogram("hh_sin",t,Binning(50),YVar(cosa,Binning(50))) ;
  hh_cos->SetLineColor(kBlue) ;
  hh_sin->SetLineColor(kBlue) ;


  // Make projection on t, plot data, pdf and its components
  // Note component projections may be larger than sum because amplitudes can be negative
  RooPlot* frame1 = t.frame();
  data->plotOn(frame1);
  pdf.plotOn(frame1);
  pdf.plotOn(frame1,Components(ampl1),LineStyle(kDashed));
  pdf.plotOn(frame1,Components(ampl2),LineStyle(kDashed),LineColor(kRed));

  // Make projection on cosa, plot data, pdf and its components
  // Note that components projection may be larger than sum because amplitudes can be negative
  RooPlot* frame2 = cosa.frame();
  data->plotOn(frame2);
  pdf.plotOn(frame2);
  pdf.plotOn(frame2,Components(ampl1),LineStyle(kDashed));
  pdf.plotOn(frame2,Components(ampl2),LineStyle(kDashed),LineColor(kRed));


  regPlot(frame1,"rf704_plot1") ;
  regPlot(frame2,"rf704_plot2") ;
  regTH(hh_cos,"rf704_hh_cos") ;
  regTH(hh_sin,"rf704_hh_sin") ;

  delete data ;
  delete coshGConv ;
  delete sinhGConv ;

  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'SPECIAL PDFS' RooFit tutorial macro #705
//
// Linear interpolation between p.d.f shapes using the 'Alex Read' algorithm
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooIntegralMorph.h"
#include "RooNLLVar.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TH1.h"
using namespace RooFit ;


// Elementary operations on a gaussian PDF
class TestBasic705 : public RooUnitTest
{
public:

  Double_t ctol() { return 5e-2 ; } // very conservative, this is a numerically difficult test

  TestBasic705(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Linear morph operator p.d.f.",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   e n d   p o i n t   p d f   s h a p e s
  // ------------------------------------------------------

  // Observable
  RooRealVar x("x","x",-20,20) ;

  // Lower end point shape: a Gaussian
  RooRealVar g1mean("g1mean","g1mean",-10) ;
  RooGaussian g1("g1","g1",x,g1mean,RooConst(2)) ;

  // Upper end point shape: a Polynomial
  RooPolynomial g2("g2","g2",x,RooArgSet(RooConst(-0.03),RooConst(-0.001))) ;



  // C r e a t e   i n t e r p o l a t i n g   p d f
  // -----------------------------------------------

  // Create interpolation variable
  RooRealVar alpha("alpha","alpha",0,1.0) ;

  // Specify sampling density on observable and interpolation variable
  x.setBins(1000,"cache") ;
  alpha.setBins(50,"cache") ;

  // Construct interpolating pdf in (x,a) represent g1(x) at a=a_min
  // and g2(x) at a=a_max
  RooIntegralMorph lmorph("lmorph","lmorph",g1,g2,x,alpha) ;



  // P l o t   i n t e r p o l a t i n g   p d f   a t   v a r i o u s   a l p h a
  // -----------------------------------------------------------------------------

  // Show end points as blue curves
  RooPlot* frame1 = x.frame() ;
  g1.plotOn(frame1) ;
  g2.plotOn(frame1) ;

  // Show interpolated shapes in red
  alpha.setVal(0.125) ;
  lmorph.plotOn(frame1,LineColor(kRed),Name("alt1")) ;
  alpha.setVal(0.25) ;
  lmorph.plotOn(frame1,LineColor(kRed),Name("alt2")) ;
  alpha.setVal(0.375) ;
  lmorph.plotOn(frame1,LineColor(kRed),Name("alt3")) ;
  alpha.setVal(0.50) ;
  lmorph.plotOn(frame1,LineColor(kRed),Name("alt4")) ;
  alpha.setVal(0.625) ;
  lmorph.plotOn(frame1,LineColor(kRed),Name("alt5")) ;
  alpha.setVal(0.75) ;
  lmorph.plotOn(frame1,LineColor(kRed),Name("alt6")) ;
  alpha.setVal(0.875) ;
  lmorph.plotOn(frame1,LineColor(kRed),Name("alt7")) ;
  alpha.setVal(0.95) ;
  lmorph.plotOn(frame1,LineColor(kRed),Name("alt8")) ;



  // S h o w   2 D   d i s t r i b u t i o n   o f   p d f ( x , a l p h a )
  // -----------------------------------------------------------------------

  // Create 2D histogram
  TH1* hh = lmorph.createHistogram("hh",x,Binning(40),YVar(alpha,Binning(40))) ;
  hh->SetLineColor(kBlue) ;


  // F i t   p d f   t o   d a t a s e t   w i t h   a l p h a = 0 . 8
  // -----------------------------------------------------------------

  // Generate a toy dataset at alpha = 0.8
  alpha=0.8 ;
  RooDataSet* data = lmorph.generate(x,1000) ;

  // Fit pdf to toy data
  lmorph.setCacheAlpha(kTRUE) ;
  lmorph.fitTo(*data) ;

  // Plot fitted pdf and data overlaid
  RooPlot* frame2 = x.frame(Bins(100)) ;
  data->plotOn(frame2) ;
  lmorph.plotOn(frame2) ;


  // S c a n   - l o g ( L )   v s   a l p h a
  // -----------------------------------------

  // Show scan -log(L) of dataset w.r.t alpha
  RooPlot* frame3 = alpha.frame(Bins(100),Range(0.5,0.9)) ;

  // Make 2D pdf of histogram
  RooNLLVar nll("nll","nll",lmorph,*data) ;
  nll.plotOn(frame3,ShiftToZero()) ;

  lmorph.setCacheAlpha(kFALSE) ;


  regPlot(frame1,"rf705_plot1") ;
  regPlot(frame2,"rf705_plot2") ;
  regPlot(frame3,"rf705_plot3") ;
  regTH(hh,"rf705_hh") ;

  delete data ;

  return kTRUE;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'SPECIAL PDFS' RooFit tutorial macro #706
//
// Histogram based p.d.f.s and functions
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;


// Elementary operations on a gaussian PDF
class TestBasic706 : public RooUnitTest
{
public:
  TestBasic706(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Histogram based p.d.f.s",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   p d f   f o r   s a m p l i n g
  // ---------------------------------------------

  RooRealVar x("x","x",0,20) ;
  RooPolynomial p("p","p",x,RooArgList(RooConst(0.01),RooConst(-0.01),RooConst(0.0004))) ;



  // C r e a t e   l o w   s t a t s   h i s t o g r a m
  // ---------------------------------------------------

  // Sample 500 events from p
  x.setBins(20) ;
  RooDataSet* data1 = p.generate(x,500) ;

  // Create a binned dataset with 20 bins and 500 events
  RooDataHist* hist1 = data1->binnedClone() ;

  // Represent data in dh as pdf in x
  RooHistPdf histpdf1("histpdf1","histpdf1",x,*hist1,0) ;

  // Plot unbinned data and histogram pdf overlaid
  RooPlot* frame1 = x.frame(Title("Low statistics histogram pdf"),Bins(100)) ;
  data1->plotOn(frame1) ;
  histpdf1.plotOn(frame1) ;


  // C r e a t e   h i g h   s t a t s   h i s t o g r a m
  // -----------------------------------------------------

  // Sample 100000 events from p
  x.setBins(10) ;
  RooDataSet* data2 = p.generate(x,100000) ;

  // Create a binned dataset with 10 bins and 100K events
  RooDataHist* hist2 = data2->binnedClone() ;

  // Represent data in dh as pdf in x, apply 2nd order interpolation
  RooHistPdf histpdf2("histpdf2","histpdf2",x,*hist2,2) ;

  // Plot unbinned data and histogram pdf overlaid
  RooPlot* frame2 = x.frame(Title("High stats histogram pdf with interpolation"),Bins(100)) ;
  data2->plotOn(frame2) ;
  histpdf2.plotOn(frame2) ;


  regPlot(frame1,"rf607_plot1") ;
  regPlot(frame2,"rf607_plot2") ;

  delete data1 ;
  delete hist1 ;
  delete data2 ;
  delete hist2 ;

  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'SPECIAL PDFS' RooFit tutorial macro #707
//
// Using non-parametric (multi-dimensional) kernel estimation p.d.f.s
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "TH1.h"
#include "RooPlot.h"
using namespace RooFit ;


// Elementary operations on a gaussian PDF
class TestBasic707 : public RooUnitTest
{
public:
  TestBasic707(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Kernel estimation p.d.f.s",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   l o w   s t a t s   1 - D   d a t a s e t
  // -------------------------------------------------------

  // Create a toy pdf for sampling
  RooRealVar x("x","x",0,20) ;
  RooPolynomial p("p","p",x,RooArgList(RooConst(0.01),RooConst(-0.01),RooConst(0.0004))) ;

  // Sample 500 events from p
  RooDataSet* data1 = p.generate(x,200) ;



  // C r e a t e   1 - D   k e r n e l   e s t i m a t i o n   p d f
  // ---------------------------------------------------------------

  // Create adaptive kernel estimation pdf. In this configuration the input data
  // is mirrored over the boundaries to minimize edge effects in distribution
  // that do not fall to zero towards the edges
  RooKeysPdf kest1("kest1","kest1",x,*data1,RooKeysPdf::MirrorBoth) ;

  // An adaptive kernel estimation pdf on the same data without mirroring option
  // for comparison
  RooKeysPdf kest2("kest2","kest2",x,*data1,RooKeysPdf::NoMirror) ;


  // Adaptive kernel estimation pdf with increased bandwidth scale factor
  // (promotes smoothness over detail preservation)
  RooKeysPdf kest3("kest3","kest3",x,*data1,RooKeysPdf::MirrorBoth,2) ;


  // Plot kernel estimation pdfs with and without mirroring over data
  RooPlot* frame = x.frame(Title("Adaptive kernel estimation pdf with and w/o mirroring"),Bins(20)) ;
  data1->plotOn(frame) ;
  kest1.plotOn(frame) ;
  kest2.plotOn(frame,LineStyle(kDashed),LineColor(kRed)) ;


  // Plot kernel estimation pdfs with regular and increased bandwidth
  RooPlot* frame2 = x.frame(Title("Adaptive kernel estimation pdf with regular, increased bandwidth")) ;
  kest1.plotOn(frame2) ;
  kest3.plotOn(frame2,LineColor(kMagenta)) ;



  // C r e a t e   l o w   s t a t s   2 - D   d a t a s e t
  // -------------------------------------------------------

  // Construct a 2D toy pdf for sampleing
  RooRealVar y("y","y",0,20) ;
  RooPolynomial py("py","py",y,RooArgList(RooConst(0.01),RooConst(0.01),RooConst(-0.0004))) ;
  RooProdPdf pxy("pxy","pxy",RooArgSet(p,py)) ;
  RooDataSet* data2 = pxy.generate(RooArgSet(x,y),1000) ;



  // C r e a t e   2 - D   k e r n e l   e s t i m a t i o n   p d f
  // ---------------------------------------------------------------

  // Create 2D adaptive kernel estimation pdf with mirroring
  RooNDKeysPdf kest4("kest4","kest4",RooArgSet(x,y),*data2,"am") ;

  // Create 2D adaptive kernel estimation pdf with mirroring and double bandwidth
  RooNDKeysPdf kest5("kest5","kest5",RooArgSet(x,y),*data2,"am",2) ;

  // Create a histogram of the data
  TH1* hh_data = data2->createHistogram("hh_data",x,Binning(10),YVar(y,Binning(10))) ;

  // Create histogram of the 2d kernel estimation pdfs
  TH1* hh_pdf = kest4.createHistogram("hh_pdf",x,Binning(25),YVar(y,Binning(25))) ;
  TH1* hh_pdf2 = kest5.createHistogram("hh_pdf2",x,Binning(25),YVar(y,Binning(25))) ;
  hh_pdf->SetLineColor(kBlue) ;
  hh_pdf2->SetLineColor(kMagenta) ;

  regPlot(frame,"rf707_plot1") ;
  regPlot(frame2,"rf707_plot2") ;
  regTH(hh_data,"rf707_hhdata") ;
  regTH(hh_pdf,"rf707_hhpdf") ;
  regTH(hh_pdf2,"rf707_hhpdf2") ;

  delete data1 ;
  delete data2 ;

  return kTRUE ;
  }
} ;
//////////////////////////////////////////////////////////////////////////
//
// 'SPECIAL PDFS' RooFit tutorial macro #708
//
// Special decay pdf for B physics with mixing and/or CP violation
//
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCategory.h"
#include "RooBMixDecay.h"
#include "RooBCPEffDecay.h"
#include "RooBDecay.h"
#include "RooFormulaVar.h"
#include "RooTruthModel.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace RooFit ;

// Elementary operations on a gaussian PDF
class TestBasic708 : public RooUnitTest
{
public:
  TestBasic708(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("B Physics p.d.f.s",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  ////////////////////////////////////////////////////
  // B - D e c a y   w i t h   m i x i n g          //
  ////////////////////////////////////////////////////

  // C o n s t r u c t   p d f
  // -------------------------

  // Observable
  RooRealVar dt("dt","dt",-10,10) ;
  dt.setBins(40) ;

  // Parameters
  RooRealVar dm("dm","delta m(B0)",0.472) ;
  RooRealVar tau("tau","tau (B0)",1.547) ;
  RooRealVar w("w","flavour mistag rate",0.1) ;
  RooRealVar dw("dw","delta mistag rate for B0/B0bar",0.1) ;

  RooCategory mixState("mixState","B0/B0bar mixing state") ;
  mixState.defineType("mixed",-1) ;
  mixState.defineType("unmixed",1) ;

  RooCategory tagFlav("tagFlav","Flavour of the tagged B0") ;
  tagFlav.defineType("B0",1) ;
  tagFlav.defineType("B0bar",-1) ;

  // Use delta function resolution model
  RooTruthModel tm("tm","truth model",dt) ;

  // Construct Bdecay with mixing
  RooBMixDecay bmix("bmix","decay",dt,mixState,tagFlav,tau,dm,w,dw,tm,RooBMixDecay::DoubleSided) ;



  // P l o t   p d f   i n   v a r i o u s   s l i c e s
  // ---------------------------------------------------

  // Generate some data
  RooDataSet* data = bmix.generate(RooArgSet(dt,mixState,tagFlav),10000) ;

  // Plot B0 and B0bar tagged data separately
  // For all plots below B0 and B0 tagged data will look somewhat differently
  // if the flavor tagging mistag rate for B0 and B0 is different (i.e. dw!=0)
  RooPlot* frame1 = dt.frame(Title("B decay distribution with mixing (B0/B0bar)")) ;

  data->plotOn(frame1,Cut("tagFlav==tagFlav::B0")) ;
  bmix.plotOn(frame1,Slice(tagFlav,"B0")) ;

  data->plotOn(frame1,Cut("tagFlav==tagFlav::B0bar"),MarkerColor(kCyan)) ;
  bmix.plotOn(frame1,Slice(tagFlav,"B0bar"),LineColor(kCyan),Name("alt")) ;


  // Plot mixed slice for B0 and B0bar tagged data separately
  RooPlot* frame2 = dt.frame(Title("B decay distribution of mixed events (B0/B0bar)")) ;

  data->plotOn(frame2,Cut("mixState==mixState::mixed&&tagFlav==tagFlav::B0")) ;
  bmix.plotOn(frame2,Slice(tagFlav,"B0"),Slice(mixState,"mixed")) ;

  data->plotOn(frame2,Cut("mixState==mixState::mixed&&tagFlav==tagFlav::B0bar"),MarkerColor(kCyan)) ;
  bmix.plotOn(frame2,Slice(tagFlav,"B0bar"),Slice(mixState,"mixed"),LineColor(kCyan),Name("alt")) ;


  // Plot unmixed slice for B0 and B0bar tagged data separately
  RooPlot* frame3 = dt.frame(Title("B decay distribution of unmixed events (B0/B0bar)")) ;

  data->plotOn(frame3,Cut("mixState==mixState::unmixed&&tagFlav==tagFlav::B0")) ;
  bmix.plotOn(frame3,Slice(tagFlav,"B0"),Slice(mixState,"unmixed")) ;

  data->plotOn(frame3,Cut("mixState==mixState::unmixed&&tagFlav==tagFlav::B0bar"),MarkerColor(kCyan)) ;
  bmix.plotOn(frame3,Slice(tagFlav,"B0bar"),Slice(mixState,"unmixed"),LineColor(kCyan),Name("alt")) ;





  ///////////////////////////////////////////////////////
  // B - D e c a y   w i t h   C P   v i o l a t i o n //
  ///////////////////////////////////////////////////////

  // C o n s t r u c t   p d f
  // -------------------------

  // Additional parameters needed for B decay with CPV
  RooRealVar CPeigen("CPeigen","CP eigen value",-1) ;
  RooRealVar absLambda("absLambda","|lambda|",1,0,2) ;
  RooRealVar argLambda("absLambda","|lambda|",0.7,-1,1) ;
  RooRealVar effR("effR","B0/B0bar reco efficiency ratio",1) ;

  // Construct Bdecay with CP violation
  RooBCPEffDecay bcp("bcp","bcp", dt, tagFlav, tau, dm, w, CPeigen, absLambda, argLambda, effR, dw, tm, RooBCPEffDecay::DoubleSided) ;



  // P l o t   s c e n a r i o   1   -   s i n ( 2 b )   =   0 . 7 ,   | l | = 1
  // ---------------------------------------------------------------------------

  // Generate some data
  RooDataSet* data2 = bcp.generate(RooArgSet(dt,tagFlav),10000) ;

  // Plot B0 and B0bar tagged data separately
  RooPlot* frame4 = dt.frame(Title("B decay distribution with CPV(|l|=1,Im(l)=0.7) (B0/B0bar)")) ;

  data2->plotOn(frame4,Cut("tagFlav==tagFlav::B0")) ;
  bcp.plotOn(frame4,Slice(tagFlav,"B0")) ;

  data2->plotOn(frame4,Cut("tagFlav==tagFlav::B0bar"),MarkerColor(kCyan)) ;
  bcp.plotOn(frame4,Slice(tagFlav,"B0bar"),LineColor(kCyan),Name("alt")) ;



  // P l o t   s c e n a r i o   2   -   s i n ( 2 b )   =   0 . 7 ,   | l | = 0 . 7
  // -------------------------------------------------------------------------------

  absLambda=0.7 ;

  // Generate some data
  RooDataSet* data3 = bcp.generate(RooArgSet(dt,tagFlav),10000) ;

  // Plot B0 and B0bar tagged data separately (sin2b = 0.7 plus direct CPV |l|=0.5)
  RooPlot* frame5 = dt.frame(Title("B decay distribution with CPV(|l|=0.7,Im(l)=0.7) (B0/B0bar)")) ;

  data3->plotOn(frame5,Cut("tagFlav==tagFlav::B0")) ;
  bcp.plotOn(frame5,Slice(tagFlav,"B0")) ;

  data3->plotOn(frame5,Cut("tagFlav==tagFlav::B0bar"),MarkerColor(kCyan)) ;
  bcp.plotOn(frame5,Slice(tagFlav,"B0bar"),LineColor(kCyan),Name("alt")) ;



  //////////////////////////////////////////////////////////////////////////////////
  // G e n e r i c   B   d e c a y  w i t h    u s e r   c o e f f i c i e n t s  //
  //////////////////////////////////////////////////////////////////////////////////

  // C o n s t r u c t   p d f
  // -------------------------

  // Model parameters
  RooRealVar DGbG("DGbG","DGamma/GammaAvg",0.5,-1,1);
  RooRealVar Adir("Adir","-[1-abs(l)**2]/[1+abs(l)**2]",0);
  RooRealVar Amix("Amix","2Im(l)/[1+abs(l)**2]",0.7);
  RooRealVar Adel("Adel","2Re(l)/[1+abs(l)**2]",0.7);

  // Derived input parameters for pdf
  RooFormulaVar DG("DG","Delta Gamma","@1/@0",RooArgList(tau,DGbG));

  // Construct coefficient functions for sin,cos,sinh modulations of decay distribution
  RooFormulaVar fsin("fsin","fsin","@0*@1*(1-2*@2)",RooArgList(Amix,tagFlav,w));
  RooFormulaVar fcos("fcos","fcos","@0*@1*(1-2*@2)",RooArgList(Adir,tagFlav,w));
  RooFormulaVar fsinh("fsinh","fsinh","@0",RooArgList(Adel));

  // Construct generic B decay pdf using above user coefficients
  RooBDecay bcpg("bcpg","bcpg",dt,tau,DG,RooConst(1),fsinh,fcos,fsin,dm,tm,RooBDecay::DoubleSided);



  // P l o t   -   I m ( l ) = 0 . 7 ,   R e ( l ) = 0 . 7   | l | = 1,   d G / G = 0 . 5
  // -------------------------------------------------------------------------------------

  // Generate some data
  RooDataSet* data4 = bcpg.generate(RooArgSet(dt,tagFlav),10000) ;

  // Plot B0 and B0bar tagged data separately
  RooPlot* frame6 = dt.frame(Title("B decay distribution with CPV(Im(l)=0.7,Re(l)=0.7,|l|=1,dG/G=0.5) (B0/B0bar)")) ;

  data4->plotOn(frame6,Cut("tagFlav==tagFlav::B0")) ;
  bcpg.plotOn(frame6,Slice(tagFlav,"B0")) ;

  data4->plotOn(frame6,Cut("tagFlav==tagFlav::B0bar"),MarkerColor(kCyan)) ;
  bcpg.plotOn(frame6,Slice(tagFlav,"B0bar"),LineColor(kCyan),Name("alt")) ;


  regPlot(frame1,"rf708_plot1") ;
  regPlot(frame2,"rf708_plot2") ;
  regPlot(frame3,"rf708_plot3") ;
  regPlot(frame4,"rf708_plot4") ;
  regPlot(frame5,"rf708_plot5") ;
  regPlot(frame6,"rf708_plot6") ;

  delete data ;
  delete data2 ;
  delete data3 ;
  delete data4 ;

  return kTRUE ;
  }
} ;
/////////////////////////////////////////////////////////////////////////
//
// 'VALIDATION AND MC STUDIES' RooFit tutorial macro #801
//
// A Toy Monte Carlo study that perform cycles of
// event generation and fittting
//
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH2.h"
#include "RooFitResult.h"
#include "TStyle.h"
#include "TDirectory.h"

using namespace RooFit ;


// Elementary operations on a gaussian PDF
class TestBasic801 : public RooUnitTest
{
public:
  TestBasic801(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("Automated MC studies",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   m o d e l
  // -----------------------

  // Declare observable x
  RooRealVar x("x","x",0,10) ;
  x.setBins(40) ;

  // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
  RooRealVar mean("mean","mean of gaussians",5,0,10) ;
  RooRealVar sigma1("sigma1","width of gaussians",0.5) ;
  RooRealVar sigma2("sigma2","width of gaussians",1) ;

  RooGaussian sig1("sig1","Signal component 1",x,mean,sigma1) ;
  RooGaussian sig2("sig2","Signal component 2",x,mean,sigma2) ;

  // Build Chebychev polynomial p.d.f.
  RooRealVar a0("a0","a0",0.5,0.,1.) ;
  RooRealVar a1("a1","a1",-0.2,-1.,1.) ;
  RooChebychev bkg("bkg","Background",x,RooArgSet(a0,a1)) ;

  // Sum the signal components into a composite signal p.d.f.
  RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
  RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),sig1frac) ;

  // Sum the composite signal and background
  RooRealVar nbkg("nbkg","number of background events,",150,0,1000) ;
  RooRealVar nsig("nsig","number of signal events",150,0,1000) ;
  RooAddPdf  model("model","g1+g2+a",RooArgList(bkg,sig),RooArgList(nbkg,nsig)) ;



  // C r e a t e   m a n a g e r
  // ---------------------------

  // Instantiate RooMCStudy manager on model with x as observable and given choice of fit options
  //
  // The Silence() option kills all messages below the PROGRESS level, leaving only a single message
  // per sample executed, and any error message that occur during fitting
  //
  // The Extended() option has two effects:
  //    1) The extended ML term is included in the likelihood and
  //    2) A poisson fluctuation is introduced on the number of generated events
  //
  // The FitOptions() given here are passed to the fitting stage of each toy experiment.
  // If Save() is specified, the fit result of each experiment is saved by the manager
  //
  // A Binned() option is added in this example to bin the data between generation and fitting
  // to speed up the study at the expemse of some precision

  RooMCStudy* mcstudy = new RooMCStudy(model,x,Binned(kTRUE),Silence(),Extended(),
                                       FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;


  // G e n e r a t e   a n d   f i t   e v e n t s
  // ---------------------------------------------

  // Generate and fit 100 samples of Poisson(nExpected) events
  mcstudy->generateAndFit(100) ;



  // E x p l o r e   r e s u l t s   o f   s t u d y
  // ------------------------------------------------

  // Make plots of the distributions of mean, the error on mean and the pull of mean
  RooPlot* frame1 = mcstudy->plotParam(mean,Bins(40)) ;
  RooPlot* frame2 = mcstudy->plotError(mean,Bins(40)) ;
  RooPlot* frame3 = mcstudy->plotPull(mean,Bins(40),FitGauss(kTRUE)) ;

  // Plot distribution of minimized likelihood
  RooPlot* frame4 = mcstudy->plotNLL(Bins(40)) ;

  regPlot(frame1,"rf801_plot1") ;
  regPlot(frame2,"rf801_plot2") ;
  regPlot(frame3,"rf801_plot3") ;
  regPlot(frame4,"rf801_plot4") ;

  delete mcstudy ;

  return kTRUE ;
  }
} ;
/////////////////////////////////////////////////////////////////////////
//
// 'VALIDATION AND MC STUDIES' RooFit tutorial macro #802
//
// RooMCStudy: using separate fit and generator models, using the chi^2 calculator model
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TDirectory.h"

using namespace RooFit ;


// Elementary operations on a gaussian PDF
class TestBasic802 : public RooUnitTest
{
public:
  TestBasic802(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("MC Study with chi^2 calculator",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   m o d e l
  // -----------------------

  // Observables, parameters
  RooRealVar x("x","x",-10,10) ;
  x.setBins(10) ;
  RooRealVar mean("mean","mean of gaussian",0) ;
  RooRealVar sigma("sigma","width of gaussian",5,1,10) ;

  // Create Gaussian pdf
  RooGaussian gauss("gauss","gaussian PDF",x,mean,sigma) ;



  // C r e a t e   m a n a g e r  w i t h   c h i ^ 2   a d d - o n   m o d u l e
  // ----------------------------------------------------------------------------

  // Create study manager for binned likelihood fits of a Gaussian pdf in 10 bins
  RooMCStudy* mcs = new RooMCStudy(gauss,x,Silence(),Binned()) ;

  // Add chi^2 calculator module to mcs
  RooChi2MCSModule chi2mod ;
  mcs->addModule(chi2mod) ;

  // Generate 200 samples of 1000 events
  mcs->generateAndFit(200,1000) ;

  // Fill histograms with distributions chi2 and prob(chi2,ndf) that
  // are calculated by RooChiMCSModule

  RooRealVar* chi2 = (RooRealVar*) mcs->fitParDataSet().get()->find("chi2") ;
  RooRealVar* prob = (RooRealVar*) mcs->fitParDataSet().get()->find("prob") ;

  TH1* h_chi2  = new TH1F("h_chi2","",40,0,20) ;
  TH1* h_prob  = new TH1F("h_prob","",40,0,1) ;

  mcs->fitParDataSet().fillHistogram(h_chi2,*chi2) ;
  mcs->fitParDataSet().fillHistogram(h_prob,*prob) ;



  // C r e a t e   m a n a g e r  w i t h   s e p a r a t e   f i t   m o d e l
  // ----------------------------------------------------------------------------

  // Create alternate pdf with shifted mean
  RooRealVar mean2("mean2","mean of gaussian 2",0.5) ;
  RooGaussian gauss2("gauss2","gaussian PDF2",x,mean2,sigma) ;

  // Create study manager with separate generation and fit model. This configuration
  // is set up to generate bad fits as the fit and generator model have different means
  // and the mean parameter is not floating in the fit
  RooMCStudy* mcs2 = new RooMCStudy(gauss2,x,FitModel(gauss),Silence(),Binned()) ;

  // Add chi^2 calculator module to mcs
  RooChi2MCSModule chi2mod2 ;
  mcs2->addModule(chi2mod2) ;

  // Generate 200 samples of 1000 events
  mcs2->generateAndFit(200,1000) ;

  // Fill histograms with distributions chi2 and prob(chi2,ndf) that
  // are calculated by RooChiMCSModule

  TH1* h2_chi2  = new TH1F("h2_chi2","",40,0,20) ;
  TH1* h2_prob  = new TH1F("h2_prob","",40,0,1) ;

  mcs2->fitParDataSet().fillHistogram(h2_chi2,*chi2) ;
  mcs2->fitParDataSet().fillHistogram(h2_prob,*prob) ;

  h_chi2->SetLineColor(kRed) ;
  h_prob->SetLineColor(kRed) ;

  regTH(h_chi2,"rf802_hist_chi2") ;
  regTH(h2_chi2,"rf802_hist2_chi2") ;
  regTH(h_prob,"rf802_hist_prob") ;
  regTH(h2_prob,"rf802_hist2_prob") ;

  delete mcs ;
  delete mcs2 ;

  return kTRUE ;
  }
} ;
/////////////////////////////////////////////////////////////////////////
//
// 'VALIDATION AND MC STUDIES' RooFit tutorial macro #803
//
// RooMCStudy: Using the randomizer and profile likelihood add-on models
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooRandomizeParamMCSModule.h"
#include "RooDLLSignificanceMCSModule.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TDirectory.h"

using namespace RooFit ;


class TestBasic803 : public RooUnitTest
{
public:
  TestBasic803(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("MC Study with param rand. and Z calc",refFile,writeRef,verbose) {} ;
  Bool_t testCode() {

  // C r e a t e   m o d e l
  // -----------------------

  // Simulation of signal and background of top quark decaying into
  // 3 jets with background

  // Observable
  RooRealVar mjjj("mjjj","m(3jet) (GeV)",100,85.,350.) ;

  // Signal component (Gaussian)
  RooRealVar mtop("mtop","m(top)",162) ;
  RooRealVar wtop("wtop","m(top) resolution",15.2) ;
  RooGaussian sig("sig","top signal",mjjj,mtop,wtop) ;

  // Background component (Chebychev)
  RooRealVar c0("c0","Chebychev coefficient 0",-0.846,-1.,1.) ;
  RooRealVar c1("c1","Chebychev coefficient 1", 0.112, 0.,1.) ;
  RooRealVar c2("c2","Chebychev coefficient 2", 0.076, 0.,1.) ;
  RooChebychev bkg("bkg","combinatorial background",mjjj,RooArgList(c0,c1,c2)) ;

  // Composite model
  RooRealVar nsig("nsig","number of signal events",53,0,1e3) ;
  RooRealVar nbkg("nbkg","number of background events",103,0,5e3) ;
  RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg)) ;



  // C r e a t e   m a n a g e r
  // ---------------------------

  // Configure manager to perform binned extended likelihood fits (Binned(),Extended()) on data generated
  // with a Poisson fluctuation on Nobs (Extended())
  RooMCStudy* mcs = new RooMCStudy(model,mjjj,Binned(),Silence(),Extended(kTRUE),
                                   FitOptions(Extended(kTRUE),PrintEvalErrors(-1))) ;



  // C u s t o m i z e   m a n a g e r
  // ---------------------------------

  // Add module that randomizes the summed value of nsig+nbkg
  // sampling from a uniform distribution between 0 and 1000
  //
  // In general one can randomize a single parameter, or a
  // sum of N parameters, using either a uniform or a Gaussian
  // distribution. Multiple randomization can be executed
  // by a single randomizer module

  RooRandomizeParamMCSModule randModule ;
  randModule.sampleSumUniform(RooArgSet(nsig,nbkg),50,500) ;
  mcs->addModule(randModule) ;


  // Add profile likelihood calculation of significance. Redo each
  // fit while keeping parameter nsig fixed to zero. For each toy,
  // the difference in -log(L) of both fits is stored, as well
  // a simple significance interpretation of the delta(-logL)
  // using Dnll = 0.5 sigma^2

  RooDLLSignificanceMCSModule sigModule(nsig,0) ;
  mcs->addModule(sigModule) ;



  // R u n   m a n a g e r ,   m a k e   p l o t s
  // ---------------------------------------------

  mcs->generateAndFit(50) ;

  // Make some plots
  RooRealVar* ngen    = (RooRealVar*) mcs->fitParDataSet().get()->find("ngen") ;
  RooRealVar* dll     = (RooRealVar*) mcs->fitParDataSet().get()->find("dll_nullhypo_nsig") ;
  RooRealVar* z       = (RooRealVar*) mcs->fitParDataSet().get()->find("significance_nullhypo_nsig") ;
  RooRealVar* nsigerr = (RooRealVar*) mcs->fitParDataSet().get()->find("nsigerr") ;

  TH1* dll_vs_ngen     = new TH2F("h_dll_vs_ngen"    ,"",40,0,500,40,0,50) ;
  TH1* z_vs_ngen       = new TH2F("h_z_vs_ngen"      ,"",40,0,500,40,0,10) ;
  TH1* errnsig_vs_ngen = new TH2F("h_nsigerr_vs_ngen","",40,0,500,40,0,30) ;
  TH1* errnsig_vs_nsig = new TH2F("h_nsigerr_vs_nsig","",40,0,200,40,0,30) ;

  mcs->fitParDataSet().fillHistogram(dll_vs_ngen,RooArgList(*ngen,*dll)) ;
  mcs->fitParDataSet().fillHistogram(z_vs_ngen,RooArgList(*ngen,*z)) ;
  mcs->fitParDataSet().fillHistogram(errnsig_vs_ngen,RooArgList(*ngen,*nsigerr)) ;
  mcs->fitParDataSet().fillHistogram(errnsig_vs_nsig,RooArgList(nsig,*nsigerr)) ;

  regTH(dll_vs_ngen,"rf803_dll_vs_ngen") ;
  regTH(z_vs_ngen,"rf803_z_vs_ngen") ;
  regTH(errnsig_vs_ngen,"rf803_errnsig_vs_ngen") ;
  regTH(errnsig_vs_nsig,"rf803_errnsig_vs_nsig") ;

  delete mcs ;

  return kTRUE ;

  }
} ;




/////////////////////////////////////////////////////////////////////////
//
// 'VALIDATION AND MC STUDIES' RooFit tutorial macro #804
//
// Using RooMCStudy on models with constrains
//
//
// 07/2008 - Wouter Verkerke
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"

using namespace RooFit ;


class TestBasic804 : public RooUnitTest
{
public:
  TestBasic804(TFile* refFile, Bool_t writeRef, Int_t verbose) : RooUnitTest("MC Studies with aux. obs. constraints",refFile,writeRef,verbose) {} ;

  Double_t htol() { return 0.1 ; } // numerically very difficult test

  Bool_t testCode() {

  // C r e a t e   m o d e l   w i t h   p a r a m e t e r   c o n s t r a i n t
  // ---------------------------------------------------------------------------

  // Observable
  RooRealVar x("x","x",-10,10) ;

  // Signal component
  RooRealVar m("m","m",0,-10,10) ;
  RooRealVar s("s","s",2,0.1,10) ;
  RooGaussian g("g","g",x,m,s) ;

  // Background component
  RooPolynomial p("p","p",x) ;

  // Composite model
  RooRealVar f("f","f",0.4,0.,1.) ;
  RooAddPdf sum("sum","sum",RooArgSet(g,p),f) ;

  // Construct constraint on parameter f
  RooGaussian fconstraint("fconstraint","fconstraint",f,RooConst(0.7),RooConst(0.1)) ;

  // Multiply constraint with p.d.f
  RooProdPdf sumc("sumc","sum with constraint",RooArgSet(sum,fconstraint)) ;



  // S e t u p   t o y   s t u d y   w i t h   m o d e l
  // ---------------------------------------------------

  // Perform toy study with internal constraint on f
  //RooMCStudy mcs(sumc,x,Constrain(f),Silence(),Binned(),FitOptions(PrintLevel(-1))) ;
  RooMCStudy mcs(sumc,x,Constrain(f),Binned()) ;

  // Run 50 toys of 2000 events.
  // Before each toy is generated, a value for the f is sampled from the constraint pdf and
  // that value is used for the generation of that toy.
  mcs.generateAndFit(50,2000) ;

  // Make plot of distribution of generated value of f parameter
  RooRealVar* f_gen = (RooRealVar*) mcs.fitParDataSet().get()->find("f_gen") ;
  TH1* h_f_gen = new TH1F("h_f_gen","",40,0,1) ;
  mcs.fitParDataSet().fillHistogram(h_f_gen,*f_gen) ;

  // Make plot of distribution of fitted value of f parameter
  RooPlot* frame1  = mcs.plotParam(f,Bins(40),Range(0.4,1)) ;
  frame1->SetTitle("Distribution of fitted f values") ;

  // Make plot of pull distribution on f
  RooPlot* frame2 = mcs.plotPull(f,Bins(40),Range(-3,3)) ;
  frame1->SetTitle("Distribution of f pull values") ;

  regTH(h_f_gen,"rf804_h_f_gen") ;
  regPlot(frame1,"rf804_plot1") ;
  regPlot(frame2,"rf804_plot2") ;

  return kTRUE ;
  }
} ;
