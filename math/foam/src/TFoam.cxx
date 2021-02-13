// @(#)root/foam:$Id$
// Author: S. Jadach <mailto:Stanislaw.jadach@ifj.edu.pl>, P.Sawicki <mailto:Pawel.Sawicki@ifj.edu.pl>

/** \class TFoam


TFoam is the main class of the multi-dimensional general purpose
Monte Carlo event generator (integrator) FOAM.

### FOAM  Version 1.02M

\authors
  S. Jadach and P.Sawicki
  Institute of Nuclear Physics, Cracow, Poland
  Stanislaw. Jadach@ifj.edu.pl, Pawel.Sawicki@ifj.edu.pl

### What is FOAM for?

 - Suppose you want to generate randomly points (vectors) according to
   an arbitrary probability distribution  in n dimensions,
   for which you supply your own method. FOAM can do it for you!
   Even if your distributions has quite strong peaks and is discontinuous!
-  FOAM generates random points with weight one or with variable weight.
 - FOAM is capable to integrate using efficient "adaptive" MC method.
   (The distribution does not need to be normalized to one.)

### How does it work?

FOAM is the simplified version of the multi-dimensional general purpose
Monte Carlo event generator (integrator) FOAM.
It creates hyper-rectangular "foam of cells", which is more dense around its peaks.
See the following 2-dim. example of the map of 1000 cells for doubly peaked distribution:

\image html foam_MapCamel1000.png width=400

FOAM is now fully integrated with the ROOT package.
The important bonus of the ROOT use is persistency of the FOAM objects!

For more sophisticated problems full version of FOAM may be more appropriate:
See [full version of FOAM](http://jadach.home.cern.ch/jadach/Foam/Index.html)

### Simple example of the use of FOAM:

Begin_Macro(source)
../../../tutorials/foam/foam_kanwa.C
End_Macro

### Canonical nine steering parameters of FOAM


| Name     | default  | Description                                            |
|----------|----------|--------------------------------------------------------|
| kDim     | 0        | Dimension of the integration space. Must be redefined! |
| nCells   | 1000     | No of allocated number of cells,                       |
| nSampl   | 200      | No. of MC events in the cell MC exploration            |
| nBin     | 8        | No. of bins in edge-histogram in cell exploration      |
| OptRej   | 1        | OptRej = 0, weighted; OptRej=1, wt=1 MC events         |
| OptDrive | 2        | Maximum weight reduction, =1 for variance reduction    |
| EvPerBin | 25       | Maximum number of the effective wt=1 events/bin,       |
|          |          | EvPerBin=0 deactivates this option                     |
| Chat     | 1        | =0,1,2 is the ``chat level'' in the standard output    |
| MaxWtRej | 1.1      | Maximum weight used to get w=1 MC events               |

The above can be redefined before calling Initialize() method,
for instance `FoamObject->SetkDim(15)` sets dimension of the distribution to 15.
Only `kDim` HAS TO BE redefined, the other parameters may be left at their defaults.
`nCell` may be increased up to about million cells for wildly peaked distributions.
Increasing `nSampl` sometimes helps, but it may cost CPU time.
`MaxWtRej` may need to be increased for wild a distribution, while using `OptRej=0`.

Past versions of FOAM: August 2003, v.1.00; September 2003 v.1.01
Adopted starting from FOAM-2.06 by P. Sawicki

Users of FOAM are kindly requested to cite the following work:
S. Jadach, Computer Physics Communications 152 (2003) 55.

*/

#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TFoamMaxwt.h"
#include "TFoamVect.h"
#include "TFoamCell.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TH1.h"
#include "TRefArray.h"
#include "TObjArray.h"
#include "TMethodCall.h"
#include "TRandom.h"
#include "TMath.h"
#include "TInterpreter.h"

ClassImp(TFoam);

//FFFFFF  BoX-FORMATs for nice and flexible outputs
#define BXOPE std::cout<<\
"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"<<std::endl<<\
"F                                                                              F"<<std::endl
#define BXTXT(text) std::cout<<\
"F                   "<<std::setw(40)<<         text           <<"                   F"<<std::endl
#define BX1I(name,numb,text) std::cout<<\
"F "<<std::setw(10)<<name<<" = "<<std::setw(10)<<numb<<" = "          <<std::setw(50)<<text<<" F"<<std::endl
#define BX1F(name,numb,text)     std::cout<<"F "<<std::setw(10)<<name<<\
          " = "<<std::setw(15)<<std::setprecision(8)<<numb<<"   =    "<<std::setw(40)<<text<<" F"<<std::endl
#define BX2F(name,numb,err,text) std::cout<<"F "<<std::setw(10)<<name<<\
" = "<<std::setw(15)<<std::setprecision(8)<<numb<<" +- "<<std::setw(15)<<std::setprecision(8)<<err<< \
                                                      "  = "<<std::setw(25)<<text<<" F"<<std::endl
#define BXCLO std::cout<<\
"F                                                                              F"<<std::endl<<\
"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"<<std::endl
  //FFFFFF  BoX-FORMATs ends here

static const Double_t gHigh= 1.0e150;
static const Double_t gVlow=-1.0e150;

#define SW2 setprecision(7) << std::setw(12)

// class to wrap a global function in a TFoamIntegrand function
class FoamIntegrandFunction : public TFoamIntegrand {

public:

   typedef Double_t (*FunctionPtr)(Int_t, Double_t*);

   FoamIntegrandFunction(FunctionPtr func) : fFunc(func) {}

   virtual ~FoamIntegrandFunction() {}

   // evaluate the density using the provided function pointer
   Double_t Density (Int_t nDim, Double_t * x) {
      return fFunc(nDim,x);
   }

private:

   FunctionPtr fFunc;

};


////////////////////////////////////////////////////////////////////////////////
/// Default constructor for streamer, user should not use it.

TFoam::TFoam() :
   fDim(0), fNCells(0), fRNmax(0),
   fOptDrive(0), fChat(0), fOptRej(0),
   fNBin(0), fNSampl(0), fEvPerBin(0),
   fMaskDiv(0), fInhiDiv(0), fOptPRD(0), fXdivPRD(0),
   fNoAct(0), fLastCe(0), fCells(0),
   fMCMonit(0), fMaxWtRej(0), fPrimAcu(0),
   fHistEdg(0), fHistDbg(0), fHistWt(0),
   fMCvect(0), fMCwt(0), fRvec(0),
   fRho(0), fMethodCall(0), fPseRan(0),
   fNCalls(0), fNEffev(0),
   fSumWt(0), fSumWt2(0),
   fSumOve(0), fNevGen(0),
   fWtMax(0), fWtMin(0),
   fPrime(0), fMCresult(0), fMCerror(0),
   fAlpha(0)
{
}
////////////////////////////////////////////////////////////////////////////////
/// User constructor, to be employed by the user

TFoam::TFoam(const Char_t* Name) :
   fDim(0), fNCells(0), fRNmax(0),
   fOptDrive(0), fChat(0), fOptRej(0),
   fNBin(0), fNSampl(0), fEvPerBin(0),
   fMaskDiv(0), fInhiDiv(0), fOptPRD(0), fXdivPRD(0),
   fNoAct(0), fLastCe(0), fCells(0),
   fMCMonit(0), fMaxWtRej(0), fPrimAcu(0),
   fHistEdg(0), fHistDbg(0), fHistWt(0),
   fMCvect(0), fMCwt(0), fRvec(0),
   fRho(0), fMethodCall(0), fPseRan(0),
   fNCalls(0), fNEffev(0),
   fSumWt(0), fSumWt2(0),
   fSumOve(0), fNevGen(0),
   fWtMax(0), fWtMin(0),
   fPrime(0), fMCresult(0), fMCerror(0),
   fAlpha(0)
{
   if(strlen(Name)  >129) {
      Error("TFoam","Name too long %s \n",Name);
   }
   fName=Name;                                            // Class name
   fDate="  Release date:  2005.04.10";                   // Release date
   fVersion= "1.02M";                                     // Release version
   fMaskDiv  = 0;             // Dynamic Mask for  cell division, h-cubic
   fInhiDiv  = 0;             // Flag allowing to inhibit cell division in certain projection/edge
   fXdivPRD  = 0;             // Lists of division values encoded in one vector per direction
   fCells    = 0;
   fAlpha    = 0;
   fPrimAcu  = 0;
   fHistEdg  = 0;
   fHistWt   = 0;
   fHistDbg  = 0;
   fDim     = 0;                // dimension of hyp-cubical space
   fNCells   = 1000;             // Maximum number of Cells,    is usually re-set
   fNSampl   = 200;              // No of sampling when dividing cell
   fOptPRD   = 0;                // General Option switch for PRedefined Division, for quick check
   fOptDrive = 2;                // type of Drive =1,2 for TrueVol,Sigma,WtMax
   fChat     = 1;                // Chat=0,1,2 chat level in output, Chat=1 normal level
   fOptRej   = 1;                // OptRej=0, wted events; OptRej=1, wt=1 events
   /////////////////////////////////////////////////////////////////////////////

   fNBin     = 8;                // binning of edge-histogram in cell exploration
   fEvPerBin =25;                // maximum no. of EFFECTIVE event per bin, =0 option is inactive
   //------------------------------------------------------
   fNCalls = 0;                  // No of function calls
   fNEffev = 0;                  // Total no of eff. wt=1 events in build=up
   fLastCe =-1;                  // Index of the last cell
   fNoAct  = 0;                  // No of active cells (used in MC generation)
   fWtMin = gHigh;               // Minimal weight
   fWtMax = gVlow;               // Maximal weight
   fMaxWtRej =1.10;              // Maximum weight in rejection for getting wt=1 events
   fPseRan   = 0;                // Initialize private copy of random number generator
   fMCMonit  = 0;                // MC efficiency monitoring
   fRho = 0;                     // pointer to abstract class providing function to integrate
   fMCvect   = 0;
   fRvec     = 0;
   fPseRan   = 0;                // generator of pseudorandom numbers
   fMethodCall=0;                // ROOT's pointer to global distribution function
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor

TFoam::~TFoam()
{
   Int_t i;

   if(fCells!= 0) {
      for(i=0; i<fNCells; i++) delete fCells[i]; // TFoamCell*[]
      delete [] fCells;
   }
   if (fRvec)    delete [] fRvec;    //double[]
   if (fAlpha)   delete [] fAlpha;   //double[]
   if (fMCvect)  delete [] fMCvect;  //double[]
   if (fPrimAcu) delete [] fPrimAcu; //double[]
   if (fMaskDiv) delete [] fMaskDiv; //int[]
   if (fInhiDiv) delete [] fInhiDiv; //int[]

   if( fXdivPRD!= 0) {
      for(i=0; i<fDim; i++) delete fXdivPRD[i]; // TFoamVect*[]
      delete [] fXdivPRD;
   }
   if (fMCMonit) delete fMCMonit;
   if (fHistWt)  delete fHistWt;

   // delete histogram arrays
   if (fHistEdg) {
      fHistEdg->Delete();
      delete fHistEdg;
   }
   if (fHistDbg) {
      fHistDbg->Delete();
      delete fHistDbg;
   }
   // delete function object if it has been created here in SetRhoInt
   if (fRho && dynamic_cast<FoamIntegrandFunction*>(fRho) ) delete fRho;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy Constructor  NOT IMPLEMENTED (NEVER USED)

TFoam::TFoam(const TFoam &From): TObject(From)
{
   Error("TFoam", "COPY CONSTRUCTOR NOT IMPLEMENTED \n");
}

////////////////////////////////////////////////////////////////////////////////
/// ### Basic initialization of FOAM invoked by the user. Mandatory!
///
/// This method starts the process of the cell build-up.
/// User must invoke Initialize with two arguments or Initialize without arguments.
/// This is done BEFORE generating first MC event and AFTER allocating FOAM object
/// and reseting (optionally) its internal parameters/switches.
/// The overall operational scheme of the FOAM is the following:
///
/// \image html foam_schema2.png width=600
///
/// ### This method invokes several other methods:
///
/// InitCells initializes memory storage for cells and begins exploration process
/// from the root cell. The empty cells are allocated/filled using  CellFill.
/// The procedure Grow which loops over cells, picks up the cell with the biggest
/// ``driver integral'', see Comp. Phys. Commun. 152 152 (2003) 55 for explanations,
/// with the help of PeekMax procedure. The chosen cell is split using Divide.
/// Subsequently, the procedure Explore called by the Divide
/// (and by InitCells for the root cell) does the most important
/// job in the FOAM object build-up: it performs a small MC run for each
/// newly allocated daughter cell.
/// Explore calculates how profitable the future split of the cell will be
/// and defines the optimal cell division geometry with the help of Carver or Varedu
/// procedures, for maximum weight or variance optimization respectively.
/// All essential results of the exploration are written into
/// the explored cell object. At the very end of the foam build-up,
/// Finally, MakeActiveList is invoked to create a list of pointers to
/// all active cells, for the purpose of the quick access during the MC generation.
/// The procedure Explore employs MakeAlpha to generate random coordinates
/// inside a given cell with the uniform distribution.
/// The above sequence of the procedure calls is depicted in the following figure:
///
/// \image html foam_Initialize_schema.png width=600

void TFoam::Initialize(TRandom *PseRan, TFoamIntegrand *fun )
{


   SetPseRan(PseRan);
   SetRho(fun);
   Initialize();
}

////////////////////////////////////////////////////////////////////////////////
/// Basic initialization of FOAM invoked by the user.
/// IMPORTANT: Random number generator and the distribution object has to be
/// provided using SetPseRan and SetRho prior to invoking this initialiser!

void TFoam::Initialize()
{
   Bool_t addStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   Int_t i;

   if(fChat>0){
      BXOPE;
      BXTXT("****************************************");
      BXTXT("******      TFoam::Initialize    ******");
      BXTXT("****************************************");
      BXTXT(fName);
      BX1F("  Version",fVersion,  fDate);
      BX1I("     kDim",fDim,     " Dimension of the hyper-cubical space             ");
      BX1I("   nCells",fNCells,   " Requested number of Cells (half of them active)  ");
      BX1I("   nSampl",fNSampl,   " No of MC events in exploration of a cell         ");
      BX1I("     nBin",fNBin,     " No of bins in histograms, MC exploration of cell ");
      BX1I(" EvPerBin",fEvPerBin, " Maximum No effective_events/bin, MC exploration  ");
      BX1I(" OptDrive",fOptDrive, " Type of Driver   =1,2 for Sigma,WtMax            ");
      BX1I("   OptRej",fOptRej,   " MC rejection on/off for OptRej=0,1               ");
      BX1F(" MaxWtRej",fMaxWtRej, " Maximum wt in rejection for wt=1 evts");
      BXCLO;
   }

   if(fPseRan==0) Error("Initialize", "Random number generator not set \n");
   if(fRho==0 && fMethodCall==0 ) Error("Initialize", "Distribution function not set \n");
   if(fDim==0) Error("Initialize", "Zero dimension not allowed \n");

   /////////////////////////////////////////////////////////////////////////
   //                   ALLOCATE SMALL LISTS                              //
   //  it is done globally, not for each cell, to save on allocation time //
   /////////////////////////////////////////////////////////////////////////
   fRNmax= fDim+1;
   fRvec = new Double_t[fRNmax];   // Vector of random numbers
   if(fRvec==0)  Error("Initialize", "Cannot initialize buffer fRvec \n");

   if(fDim>0){
      fAlpha = new Double_t[fDim];    // sum<1 for internal parametrization of the simplex
      if(fAlpha==0)  Error("Initialize", "Cannot initialize buffer fAlpha \n" );
   }
   fMCvect = new Double_t[fDim]; // vector generated in the MC run
   if(fMCvect==0)  Error("Initialize", "Cannot initialize buffer fMCvect  \n" );

   //====== List of directions inhibited for division
   if(fInhiDiv == 0){
      fInhiDiv = new Int_t[fDim];
      for(i=0; i<fDim; i++) fInhiDiv[i]=0;
   }
   //====== Dynamic mask used in Explore for edge determination
   if(fMaskDiv == 0){
      fMaskDiv = new Int_t[fDim];
      for(i=0; i<fDim; i++) fMaskDiv[i]=1;
   }
   //====== List of predefined division values in all directions (initialized as empty)
   if(fXdivPRD == 0){
      fXdivPRD = new TFoamVect*[fDim];
      for(i=0; i<fDim; i++)  fXdivPRD[i]=0; // Artificially extended beyond fDim
   }
   //====== Initialize list of histograms
   fHistWt  = new TH1D("HistWt","Histogram of MC weight",100,0.0, 1.5*fMaxWtRej); // MC weight
   fHistEdg = new TObjArray(fDim);           // Initialize list of histograms
   TString hname;
   TString htitle;
   for(i=0;i<fDim;i++){
      hname=fName+TString("_HistEdge_");
      hname+=i;
      htitle=TString("Edge Histogram No. ");
      htitle+=i;
      //std::cout<<"i= "<<i<<"  hname= "<<hname<<"  htitle= "<<htitle<<std::endl;
      (*fHistEdg)[i] = new TH1D(hname.Data(),htitle.Data(),fNBin,0.0, 1.0); // Initialize histogram for each edge
      ((TH1D*)(*fHistEdg)[i])->Sumw2();
   }
   //======  extra histograms for debug purposes
   fHistDbg = new TObjArray(fDim);         // Initialize list of histograms
   for(i=0;i<fDim;i++){
      hname=fName+TString("_HistDebug_");
      hname+=i;
      htitle=TString("Debug Histogram ");
      htitle+=i;
      (*fHistDbg)[i] = new TH1D(hname.Data(),htitle.Data(),fNBin,0.0, 1.0); // Initialize histogram for each edge
   }

   // ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| //
   //                     BUILD-UP of the FOAM                            //
   // ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| //
   //
   //        Define and explore root cell(s)
   InitCells();
   //        PrintCells(); std::cout<<" ===== after InitCells ====="<<std::endl;
   Grow();
   //        PrintCells(); std::cout<<" ===== after Grow      ====="<<std::endl;

   MakeActiveList(); // Final Preparations for the M.C. generation

   // Preparations for the M.C. generation
   fSumWt  = 0.0;               // M.C. generation sum of Wt
   fSumWt2 = 0.0;               // M.C. generation sum of Wt**2
   fSumOve = 0.0;               // M.C. generation sum of overweighted
   fNevGen = 0.0;               // M.C. generation sum of 1d0
   fWtMax  = gVlow;               // M.C. generation maximum wt
   fWtMin  = gHigh;               // M.C. generation minimum wt
   fMCresult=fCells[0]->GetIntg(); // M.C. Value of INTEGRAL,temporary assignment
   fMCresult=fCells[0]->GetIntg(); // M.C. Value of INTEGRAL,temporary assignment
   fMCerror =fCells[0]->GetIntg(); // M.C. Value of ERROR   ,temporary assignment
   fMCMonit = new TFoamMaxwt(5.0,1000);  // monitoring M.C. efficiency
   //
   if(fChat>0){
      Double_t driver = fCells[0]->GetDriv();
      BXOPE;
      BXTXT("***  TFoam::Initialize FINISHED!!!  ***");
      BX1I("    nCalls",fNCalls,  "Total number of function calls         ");
      BX1F("    XPrime",fPrime,   "Primary total integral                 ");
      BX1F("    XDiver",driver,    "Driver  total integral                 ");
      BX1F("  mcResult",fMCresult,"Estimate of the true MC Integral       ");
      BXCLO;
   }
   if(fChat==2) PrintCells();
   TH1::AddDirectory(addStatus);
} // Initialize

////////////////////////////////////////////////////////////////////////////////
/// Internal method used by Initialize.
/// It initializes "root part" of the FOAM of the tree of cells.

void TFoam::InitCells()
{
   Int_t i;

   fLastCe =-1;                             // Index of the last cell
   if(fCells!= 0) {
      for(i=0; i<fNCells; i++) delete fCells[i];
      delete [] fCells;
   }
   //
   fCells = new TFoamCell*[fNCells];
   for(i=0;i<fNCells;i++){
      fCells[i]= new TFoamCell(fDim); // Allocate BIG list of cells
      fCells[i]->SetSerial(i);
   }
   if(fCells==0) Error("InitCells", "Cannot initialize CELLS \n"  );

   /////////////////////////////////////////////////////////////////////////////
   //              Single Root Hypercube                                      //
   /////////////////////////////////////////////////////////////////////////////
   CellFill(1,   0);  //  0-th cell ACTIVE

   // Exploration of the root cell(s)
   for(Long_t iCell=0; iCell<=fLastCe; iCell++){
      Explore( fCells[iCell] );               // Exploration of root cell(s)
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method used by Initialize.
/// It initializes content of the newly allocated active cell.

Int_t TFoam::CellFill(Int_t Status, TFoamCell *parent)
{
   TFoamCell *cell;
   if (fLastCe==fNCells){
      Error( "CellFill", "Too many cells\n");
   }
   fLastCe++;   // 0-th cell is the first
   if (Status==1) fNoAct++;

   cell = fCells[fLastCe];

   cell->Fill(Status, parent, 0, 0);

   cell->SetBest( -1);         // pointer for planning division of the cell
   cell->SetXdiv(0.5);         // factor for division
   Double_t xInt2,xDri2;
   if(parent!=0){
      xInt2  = 0.5*parent->GetIntg();
      xDri2  = 0.5*parent->GetDriv();
      cell->SetIntg(xInt2);
      cell->SetDriv(xDri2);
   }else{
      cell->SetIntg(0.0);
      cell->SetDriv(0.0);
   }
   return fLastCe;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method used by Initialize.
///
/// It explores newly defined cell with help of special short MC sampling.
/// As a result, estimates of true and drive volume is defined/determined
/// Average and dispersion of the weight distribution will is found along
/// each edge and the best edge (minimum dispersion, best maximum weight)
/// is memorized for future use.
///
/// The optimal division point for eventual future cell division is
/// determined/recorded. Recorded are also minimum and maximum weight etc.
/// The volume estimate in all (inactive) parent cells is updated.
/// Note that links to parents and initial volume = 1/2 parent has to be
/// already defined prior to calling this routine.

void TFoam::Explore(TFoamCell *cell)
{
   Double_t wt, dx, xBest=0, yBest=0;
   Double_t intOld, driOld;

   Long_t iev;
   Double_t nevMC;
   Int_t i, j, k;
   Int_t nProj, kBest;
   Double_t ceSum[5], xproj;

   TFoamVect  cellSize(fDim);
   TFoamVect  cellPosi(fDim);

   cell->GetHcub(cellPosi,cellSize);

   TFoamCell  *parent;

   Double_t *xRand = new Double_t[fDim];

   Double_t *volPart=0;

   cell->CalcVolume();
   dx = cell->GetVolume();
   intOld = cell->GetIntg(); //memorize old values,
   driOld = cell->GetDriv(); //will be needed for correcting parent cells


   /////////////////////////////////////////////////////
   //    Special Short MC sampling to probe cell      //
   /////////////////////////////////////////////////////
   ceSum[0]=0;
   ceSum[1]=0;
   ceSum[2]=0;
   ceSum[3]=gHigh;  //wtmin
   ceSum[4]=gVlow;  //wtmax
   //
   for(i=0;i<fDim;i++) ((TH1D *)(*fHistEdg)[i])->Reset(); // Reset histograms
   fHistWt->Reset();
   //
   // ||||||||||||||||||||||||||BEGIN MC LOOP|||||||||||||||||||||||||||||
   Double_t nevEff=0.;
   for(iev=0;iev<fNSampl;iev++){
      MakeAlpha();               // generate uniformly vector inside hypercube

      if(fDim>0){
      for(j=0; j<fDim; j++)
         xRand[j]= cellPosi[j] +fAlpha[j]*(cellSize[j]);
      }

      wt=dx*Eval(xRand);

      nProj = 0;
      if(fDim>0) {
         for(k=0; k<fDim; k++) {
            xproj =fAlpha[k];
            ((TH1D *)(*fHistEdg)[nProj])->Fill(xproj,wt);
            nProj++;
         }
      }
      //
      fNCalls++;
      ceSum[0] += wt;    // sum of weights
      ceSum[1] += wt*wt; // sum of weights squared
      ceSum[2]++;        // sum of 1
      if (ceSum[3]>wt) ceSum[3]=wt;  // minimum weight;
      if (ceSum[4]<wt) ceSum[4]=wt;  // maximum weight
      // test MC loop exit condition
      nevEff = ceSum[1] == 0. ? 0. : ceSum[0]*ceSum[0]/ceSum[1];
      if( nevEff >= fNBin*fEvPerBin) break;
   }   // ||||||||||||||||||||||||||END MC LOOP|||||||||||||||||||||||||||||
   //------------------------------------------------------------------
   //---  predefine logics of searching for the best division edge ---
   for(k=0; k<fDim;k++){
      fMaskDiv[k] =1;                       // default is all
      if( fInhiDiv[k]==1) fMaskDiv[k] =0; // inhibit some...
   }
   // Note that predefined division below overrule inhibition above
   kBest=-1;
   Double_t rmin,rmax,rdiv;
   if(fOptPRD) {          // quick check
      for(k=0; k<fDim; k++) {
         rmin= cellPosi[k];
         rmax= cellPosi[k] +cellSize[k];
         if( fXdivPRD[k] != 0) {
            Int_t n= (fXdivPRD[k])->GetDim();
            for(j=0; j<n; j++) {
               rdiv=(*fXdivPRD[k])[j];
               // check predefined divisions is available in this cell
               if( (rmin +1e-99 <rdiv) && (rdiv< rmax -1e-99)) {
                  kBest=k;
                  xBest= (rdiv-cellPosi[k])/cellSize[k] ;
                  goto ee05;
               }
            }
         }
      }//k
   }
   ee05:
   /////////////////////////////////////////////////////////////////////////////

   fNEffev += (Long_t)nevEff;
   nevMC          = ceSum[2];
   Double_t intTrue = ceSum[0]/(nevMC+0.000001);
   Double_t intDriv=0.;
   Double_t intPrim=0.;

   switch(fOptDrive){
   case 1:                       // VARIANCE REDUCTION
      if(kBest == -1) Varedu(ceSum,kBest,xBest,yBest); // determine the best edge,
      //intDriv =sqrt( ceSum[1]/nevMC -intTrue*intTrue ); // Older ansatz, numerically not bad
      intDriv =sqrt(ceSum[1]/nevMC) -intTrue; // Foam build-up, sqrt(<w**2>) -<w>
      intPrim =sqrt(ceSum[1]/nevMC);          // MC gen. sqrt(<w**2>) =sqrt(<w>**2 +sigma**2)
      break;
   case 2:                       // WTMAX  REDUCTION
      if(kBest == -1) Carver(kBest,xBest,yBest);  // determine the best edge
      intDriv =ceSum[4] -intTrue; // Foam build-up, wtmax-<w>
      intPrim =ceSum[4];          // MC generation, wtmax!
      break;
   default:
      Error("Explore", "Wrong fOptDrive = \n" );
   }//switch
   //=================================================================================
   //hist_Neff_distrib.Fill( fLastCe/2.0+0.01, nevEff+0.01);  //
   //hist_kBest_distrib.Fill( kBest+0.50, 1.0 ); //  debug
   //hist_xBest_distrib.Fill( xBest+0.01, 1.0 ); //  debug
   //=================================================================================
   cell->SetBest(kBest);
   cell->SetXdiv(xBest);
   cell->SetIntg(intTrue);
   cell->SetDriv(intDriv);
   cell->SetPrim(intPrim);
   // correct/update integrals in all parent cells to the top of the tree
   Double_t  parIntg, parDriv;
   for(parent = cell->GetPare(); parent!=0; parent = parent->GetPare()){
      parIntg = parent->GetIntg();
      parDriv = parent->GetDriv();
      parent->SetIntg( parIntg   +intTrue -intOld );
      parent->SetDriv( parDriv   +intDriv -driOld );
   }
   delete [] volPart;
   delete [] xRand;
   //cell->Print();
} // TFoam::Explore

////////////////////////////////////////////////////////////////////////////////
/// Internal method used by Initialize.
///
/// It determines the best edge candidate and the position of the cell division plane
/// in case of the variance reduction for future cell division,
/// using results of the MC exploration run stored in fHistEdg

void TFoam::Varedu(Double_t ceSum[5], Int_t &kBest, Double_t &xBest, Double_t &yBest)
{
   Double_t nent   = ceSum[2];
   Double_t swAll  = ceSum[0];
   Double_t sswAll = ceSum[1];
   Double_t ssw    = sqrt(sswAll)/sqrt(nent);
   //
   Double_t swIn,swOut,sswIn,sswOut,xLo,xUp;
   kBest =-1;
   xBest =0.5;
   yBest =1.0;
   Double_t maxGain=0.0;
   // Now go over all projections kProj
   for(Int_t kProj=0; kProj<fDim; kProj++) {
      if( fMaskDiv[kProj]) {
         // initialize search over bins
         Double_t sigmIn =0.0; Double_t sigmOut =0.0;
         Double_t sswtBest = gHigh;
         Double_t gain =0.0;
         Double_t xMin=0.0; Double_t xMax=0.0;
         // Double loop over all pairs jLo<jUp
         for(Int_t jLo=1; jLo<=fNBin; jLo++) {
            Double_t aswIn=0;  Double_t asswIn=0;
            for(Int_t jUp=jLo; jUp<=fNBin;jUp++) {
               aswIn  +=     ((TH1D *)(*fHistEdg)[kProj])->GetBinContent(jUp);
               asswIn += Sqr(((TH1D *)(*fHistEdg)[kProj])->GetBinError(  jUp));
               xLo=(jLo-1.0)/fNBin;
               xUp=(jUp*1.0)/fNBin;
               swIn  =        aswIn/nent;
               swOut = (swAll-aswIn)/nent;
               sswIn = sqrt(asswIn)       /sqrt(nent*(xUp-xLo))     *(xUp-xLo);
               sswOut= sqrt(sswAll-asswIn)/sqrt(nent*(1.0-xUp+xLo)) *(1.0-xUp+xLo);
               if( (sswIn+sswOut) < sswtBest) {
                  sswtBest = sswIn+sswOut;
                  gain     = ssw-sswtBest;
                  sigmIn   = sswIn -swIn;  // Debug
                  sigmOut  = sswOut-swOut; // Debug
                  xMin    = xLo;
                  xMax    = xUp;
               }
            }//jUp
         }//jLo
         Int_t iLo = (Int_t) (fNBin*xMin);
         Int_t iUp = (Int_t) (fNBin*xMax);
         //----------DEBUG printout
         //std::cout<<"@@@@@  xMin xMax = "<<xMin   <<" "<<xMax<<"  iLo= "<<iLo<<"  iUp= "<<iUp;
         //std::cout<<"  sswtBest/ssw= "<<sswtBest/ssw<<"  Gain/ssw= "<< Gain/ssw<<std::endl;
         //----------DEBUG auxiliary Plot
         for(Int_t iBin=1;iBin<=fNBin;iBin++) {
            if( ((iBin-0.5)/fNBin > xMin) && ((iBin-0.5)/fNBin < xMax) ){
               ((TH1D *)(*fHistDbg)[kProj])->SetBinContent(iBin,sigmIn/(xMax-xMin));
            } else {
               ((TH1D *)(*fHistDbg)[kProj])->SetBinContent(iBin,sigmOut/(1-xMax+xMin));
            }
         }
         if(gain>=maxGain) {
            maxGain=gain;
            kBest=kProj; // <--- !!!!! The best edge
            xBest=xMin;
            yBest=xMax;
            if(iLo == 0)     xBest=yBest; // The best division point
            if(iUp == fNBin) yBest=xBest; // this is not really used
         }
      }
   } //kProj
   //----------DEBUG printout
   //std::cout<<"@@@@@@@>>>>> kBest= "<<kBest<<"  maxGain/ssw= "<< maxGain/ssw<<std::endl;
   if( (kBest >= fDim) || (kBest<0) ) Error("Varedu", "Something wrong with kBest - kBest = %d dim = %d\n",kBest,fDim);
}          //TFoam::Varedu

////////////////////////////////////////////////////////////////////////////////
/// Internal method used by Initialize.
///
/// Determines the best edge-candidate and the position of the division plane
/// for the future cell division, in the case of the optimization of the maximum weight.
/// It exploits results of the cell MC exploration run stored in fHistEdg.

void TFoam::Carver(Int_t &kBest, Double_t &xBest, Double_t &yBest)
{
   Int_t    kProj,iBin;
   Double_t carve,carvTot,carvMax,carvOne,binMax,binTot;
   Int_t    jLow,jUp,iLow,iUp;
   Double_t theBin;
   // Int_t    jDivi; // TEST
   Int_t j;

   Double_t *bins  = new Double_t[fNBin];      // bins of histogram for single  PROJECTION
   if(bins==0)    Error("Carver", "Cannot initialize buffer Bins \n" );

   kBest =-1;
   xBest =0.5;
   yBest =1.0;
   carvMax = gVlow;
   // primMax = gVlow;
   for(kProj=0; kProj<fDim; kProj++)
      if( fMaskDiv[kProj] ) {
         //if( kProj==1 ){
         //std::cout<<"==================== Carver histogram: kProj ="<<kProj<<"==================="<<std::endl;
         //((TH1D *)(*fHistEdg)[kProj])->Print("all");
         binMax = gVlow;
         for(iBin=0; iBin<fNBin;iBin++){
            bins[iBin]= ((TH1D *)(*fHistEdg)[kProj])->GetBinContent(iBin+1);
            binMax = TMath::Max( binMax, bins[iBin]);       // Maximum content/bin
         }
         if(binMax < 0 ) {       //case of empty cell
            delete [] bins;
            return;
         }
         carvTot = 0.0;
         binTot  = 0.0;
         for(iBin=0;iBin<fNBin;iBin++){
            carvTot = carvTot + (binMax-bins[iBin]);     // Total Carve (more stable)
            binTot  +=bins[iBin];
         }
         // primTot = binMax*fNBin;
         //std::cout <<"Carver:  CarvTot "<<CarvTot<< "    primTot "<<primTot<<std::endl;
         jLow =0;
         jUp  =fNBin-1;
         carvOne = gVlow;
         Double_t yLevel = gVlow;
         for(iBin=0; iBin<fNBin;iBin++) {
            theBin = bins[iBin];
            //-----  walk to the left and find first bin > theBin
            iLow = iBin;
            for(j=iBin; j>-1; j-- ) {
               if(theBin< bins[j]) break;
               iLow = j;
            }
            //iLow = iBin;
            //if(iLow>0)     while( (theBin >= bins[iLow-1])&&(iLow >0) ){iLow--;} // horror!!!
            //------ walk to the right and find first bin > theBin
            iUp  = iBin;
            for(j=iBin; j<fNBin; j++) {
               if(theBin< bins[j]) break;
               iUp = j;
            }
            //iUp  = iBin;
            //if(iUp<fNBin-1) while( (theBin >= bins[iUp+1])&&( iUp<fNBin-1 ) ){iUp++;} // horror!!!
            //
            carve = (iUp-iLow+1)*(binMax-theBin);
            if( carve > carvOne) {
               carvOne = carve;
               jLow = iLow;
               jUp  = iUp;
               yLevel = theBin;
            }
         }//iBin
         if( carvTot > carvMax) {
            carvMax   = carvTot;
            //primMax   = primTot;
            //std::cout <<"Carver:   primMax "<<primMax<<std::endl;
            kBest = kProj;    // Best edge
            xBest = ((Double_t)(jLow))/fNBin;
            yBest = ((Double_t)(jUp+1))/fNBin;
            if(jLow == 0 )       xBest = yBest;
            if(jUp  == fNBin-1) yBest = xBest;
            // division ratio in units of 1/fNBin, testing
            // jDivi = jLow;
            // if(jLow == 0 )     jDivi=jUp+1;
         }
         //======  extra histograms for debug purposes
         //std::cout<<"kProj= "<<kProj<<" jLow= "<<jLow<<" jUp= "<<jUp<<std::endl;
         for(iBin=0;    iBin<fNBin;  iBin++)
            ((TH1D *)(*fHistDbg)[kProj])->SetBinContent(iBin+1,binMax);
         for(iBin=jLow; iBin<jUp+1;   iBin++)
            ((TH1D *)(*fHistDbg)[kProj])->SetBinContent(iBin+1,yLevel);
      }//kProj
   if( (kBest >= fDim) || (kBest<0) ) Error("Carver", "Something wrong with kBest - kBest = %d dim = %d\n",kBest,fDim);
   delete [] bins;
}          //TFoam::Carver

////////////////////////////////////////////////////////////////////////////////
/// Internal method used by Initialize.
/// Provides random vector Alpha  0< Alpha(i) < 1

void TFoam::MakeAlpha()
{
   Int_t k;
   if(fDim<1) return;

   // simply generate and load kDim uniform random numbers
   fPseRan->RndmArray(fDim,fRvec);   // kDim random numbers needed
   for(k=0; k<fDim; k++) fAlpha[k] = fRvec[k];
} //MakeAlpha


////////////////////////////////////////////////////////////////////////////////
/// Internal method used by Initialize.
/// It grow new cells by the binary division process.

void TFoam::Grow()
{
   Long_t iCell;
   TFoamCell* newCell;

   while ( (fLastCe+2) < fNCells ) {  // this condition also checked inside Divide
      iCell   = PeekMax();            // peek up cell with maximum driver integral
      if( (iCell<0) || (iCell>fLastCe) ) {
        Error("Grow", "Wrong iCell \n");
        break;
      }
      newCell = fCells[iCell];

      if(fLastCe !=0) {
         Int_t kEcho=10;
         if(fLastCe>=10000) kEcho=100;
         if( (fLastCe%kEcho)==0) {
            if (fChat>0) {
               if(fDim<10)
                  std::cout<<fDim<<std::flush;
               else
                  std::cout<<"."<<std::flush;
               if( (fLastCe%(100*kEcho))==0)  std::cout<<"|"<<fLastCe<<std::endl<<std::flush;
            }
         }
      }
      if( Divide( newCell )==0) break;  // and divide it into two
   }
   if (fChat>0) {
      std::cout<<std::endl<<std::flush;
   }
   CheckAll(0);   // set arg=1 for more info
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method used by Initialize.
/// It finds cell with maximal driver integral for the purpose of the division.

Long_t  TFoam::PeekMax()
{
   Long_t  i;
   Long_t iCell = -1;
   Double_t  drivMax, driv;

   drivMax = gVlow;
   for(i=0; i<=fLastCe; i++) {//without root
      if( fCells[i]->GetStat() == 1 ) {
         driv =  TMath::Abs( fCells[i]->GetDriv());
         //std::cout<<"PeekMax: Driv = "<<driv<<std::endl;
         if(driv > drivMax) {
            drivMax = driv;
            iCell = i;
         }
      }
   }
   //  std::cout << 'TFoam_PeekMax: iCell=' << iCell << std::endl;
   if (iCell == -1)
      std::cout << "STOP in TFoam::PeekMax: not found iCell=" <<  iCell << std::endl;
   return(iCell);
}                 // TFoam_PeekMax

////////////////////////////////////////////////////////////////////////////////
/// Internal method used by Initialize.
///
/// It divides cell iCell into two daughter cells.
/// The iCell is retained and tagged as inactive, daughter cells are appended
/// at the end of the buffer.
/// New vertex is added to list of vertices.
/// List of active cells is updated, iCell removed, two daughters added
/// and their properties set with help of MC sampling (TFoam_Explore)
/// Returns Code RC=-1 of buffer limit is reached,  fLastCe=fnBuf.

Int_t TFoam::Divide(TFoamCell *cell)
{
   Int_t   kBest;

   if(fLastCe+1 >= fNCells) Error("Divide", "Buffer limit is reached, fLastCe=fnBuf \n");

   cell->SetStat(0); // reset cell as inactive
   fNoAct--;

   kBest = cell->GetBest();
   if( kBest<0 || kBest>=fDim ) Error("Divide", "Wrong kBest \n");

   //////////////////////////////////////////////////////////////////
   //           define two daughter cells (active)                 //
   //////////////////////////////////////////////////////////////////

   Int_t d1 = CellFill(1,   cell);
   Int_t d2 = CellFill(1,   cell);
   cell->SetDau0((fCells[d1]));
   cell->SetDau1((fCells[d2]));
   Explore( (fCells[d1]) );
   Explore( (fCells[d2]) );
   return 1;
} // TFoam_Divide


////////////////////////////////////////////////////////////////////////////////
/// Internal method used by Initialize.
///
/// It finds out number of active cells fNoAct,
/// creates list of active cell fCellsAct and primary cumulative fPrimAcu.
/// They are used during the MC generation to choose randomly an active cell.

void TFoam::MakeActiveList()
{
   Long_t iCell;
   Double_t sum;
   // flush previous result
   if(fPrimAcu  != 0) delete [] fPrimAcu;
   fCellsAct.clear();
   fCellsAct.reserve(fNoAct);

   // Count Active cells and find total Primary
   // Fill-in tables of active cells

   fPrime = 0.0;
   for(iCell=0; iCell<=fLastCe; iCell++) {
      if (fCells[iCell]->GetStat()==1) {
         fPrime += fCells[iCell]->GetPrim();
         fCellsAct.push_back(iCell);
      }
   }

   if(fNoAct != static_cast<Int_t>(fCellsAct.size()))  Error("MakeActiveList", "Wrong fNoAct               \n"  );
   if(fPrime == 0.) Error("MakeActiveList", "Integrand function is zero  \n"  );

   fPrimAcu  = new  Double_t[fNoAct]; // cumulative primary for MC generation
   if( fPrimAcu==0 ) Error("MakeActiveList", "Cant allocate fCellsAct or fPrimAcu \n");

   sum =0.0;
   for(iCell=0; iCell<fNoAct; iCell++) {
      sum = sum + ( fCells[fCellsAct[iCell]] )->GetPrim()/fPrime;
      fPrimAcu[iCell]=sum;
   }

} //MakeActiveList

////////////////////////////////////////////////////////////////////////////////
/// User may optionally reset random number generator using this method.
///
/// Usually it is done when FOAM object is restored from the disk.
/// IMPORTANT: this method deletes existing  random number generator registered in the FOAM object.
/// In particular such an object is created by the streamer during the disk-read operation.

void TFoam::ResetPseRan(TRandom *PseRan)
{
   if(fPseRan) {
      Info("ResetPseRan", "Resetting random number generator  \n");
      delete fPseRan;
   }
   SetPseRan(PseRan);
}

////////////////////////////////////////////////////////////////////////////////
/// User may use this method to set the distribution object

void TFoam::SetRho(TFoamIntegrand *fun)
{
   if (fun)
      fRho=fun;
   else
      Error("SetRho", "Bad function \n" );
}
////////////////////////////////////////////////////////////////////////////////
/// User may use this method to set the distribution object as a global function pointer
/// (and not as an interpreted function).

void TFoam::SetRhoInt(Double_t (*fun)(Int_t, Double_t *) )
{
   // This is needed for both AClic and Cling
   if (fun) {
      // delete function object if it has been created here in SetRho
      if (fRho && dynamic_cast<FoamIntegrandFunction*>(fRho) ) delete fRho;
      fRho= new FoamIntegrandFunction(fun);
   } else
      Error("SetRho", "Bad function \n" );
}

////////////////////////////////////////////////////////////////////////////////
/// User may optionally reset the distribution using this method.
///
/// Usually it is done when FOAM object is restored from the disk.
/// IMPORTANT: this method deletes existing  distribution object registered in the FOAM object.
/// In particular such an object is created by the streamer diring the disk-read operation.
/// This method is used only in very special cases, because the distribution in most cases
/// should be "owned" by the FOAM object and should not be replaced by another one after initialization.

void TFoam::ResetRho(TFoamIntegrand *fun)
{
   if(fRho) {
      Info("ResetRho", "!!! Resetting distribution function  !!!\n");
      delete fRho;
   }
   SetRho(fun);
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method.
/// Evaluates distribution to be generated.

Double_t TFoam::Eval(Double_t *xRand)
{
   Double_t result;

   if(!fRho) {   //interactive mode
      Long_t paramArr[3];
      paramArr[0]=(Long_t)fDim;
      paramArr[1]=(Long_t)xRand;
      fMethodCall->SetParamPtrs(paramArr);
      fMethodCall->Execute(result);
   } else {       //compiled mode
      result=fRho->Density(fDim,xRand);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method.
/// Return randomly chosen active cell with probability equal to its
/// contribution into total driver integral using interpolation search.

void TFoam::GenerCel2(TFoamCell *&pCell)
{
   Long_t  lo, hi, hit;
   Double_t fhit, flo, fhi;
   Double_t random;

   random=fPseRan->Rndm();
   lo  = 0;              hi =fNoAct-1;
   flo = fPrimAcu[lo];  fhi=fPrimAcu[hi];
   while(lo+1<hi) {
      hit = lo + (Int_t)( (hi-lo)*(random-flo)/(fhi-flo)+0.5);
      if (hit<=lo)
         hit = lo+1;
      else if(hit>=hi)
         hit = hi-1;
      fhit=fPrimAcu[hit];
      if (fhit>random) {
         hi = hit;
         fhi = fhit;
      } else {
         lo = hit;
         flo = fhit;
      }
   }
   if (fPrimAcu[lo]>random)
      pCell = fCells[fCellsAct[lo]];
   else
      pCell = fCells[fCellsAct[hi]];
}       // TFoam::GenerCel2


////////////////////////////////////////////////////////////////////////////////
/// User method.
/// It generates randomly point/vector according to user-defined distribution.
/// Prior initialization with help of Initialize() is mandatory.
/// Generated MC point/vector is available using GetMCvect and the MC weight with GetMCwt.
/// MC point is generated with wt=1 or with variable weight, see OptRej switch.

void TFoam::MakeEvent(void)
{
   Int_t      j;
   Double_t   wt,dx,mcwt;
   TFoamCell *rCell;
   //
   //********************** MC LOOP STARS HERE **********************
ee0:
   GenerCel2(rCell);   // choose randomly one cell

   MakeAlpha();

   TFoamVect  cellPosi(fDim); TFoamVect  cellSize(fDim);
   rCell->GetHcub(cellPosi,cellSize);
   for(j=0; j<fDim; j++)
      fMCvect[j]= cellPosi[j] +fAlpha[j]*cellSize[j];
   dx = rCell->GetVolume();      // Cartesian volume of the Cell
   //  weight average normalized to PRIMARY integral over the cell

   wt=dx*Eval(fMCvect);

   mcwt = wt / rCell->GetPrim();  // PRIMARY controls normalization
   fNCalls++;
   fMCwt   =  mcwt;
   // accumulation of statistics for the main MC weight
   fSumWt  += mcwt;           // sum of Wt
   fSumWt2 += mcwt*mcwt;      // sum of Wt**2
   fNevGen++;                 // sum of 1d0
   fWtMax  =  TMath::Max(fWtMax, mcwt);   // maximum wt
   fWtMin  =  TMath::Min(fWtMin, mcwt);   // minimum wt
   fMCMonit->Fill(mcwt);
   fHistWt->Fill(mcwt,1.0);          // histogram
   //*******  Optional rejection ******
   if(fOptRej == 1) {
      Double_t random;
      random=fPseRan->Rndm();
      if( fMaxWtRej*random > fMCwt) goto ee0;  // Wt=1 events, internal rejection
      if( fMCwt<fMaxWtRej ) {
         fMCwt = 1.0;                  // normal Wt=1 event
      } else {
         fMCwt = fMCwt/fMaxWtRej;    // weight for overweighted events! kept for debug
         fSumOve += fMCwt-fMaxWtRej; // contribution of overweighted
      }
   }
   //********************** MC LOOP ENDS HERE **********************
} // MakeEvent

////////////////////////////////////////////////////////////////////////////////
/// User may get generated MC point/vector with help of this method

void TFoam::GetMCvect(Double_t *MCvect)
{
   for ( Int_t k=0 ; k<fDim ; k++) *(MCvect +k) = fMCvect[k];
}

////////////////////////////////////////////////////////////////////////////////
/// User may get weight MC weight using this method

Double_t TFoam::GetMCwt(void)
{
   return(fMCwt);
}
////////////////////////////////////////////////////////////////////////////////
/// User may get weight MC weight using this method

void TFoam::GetMCwt(Double_t &mcwt)
{
   mcwt=fMCwt;
}

////////////////////////////////////////////////////////////////////////////////
/// User method which generates MC event and returns MC weight

Double_t TFoam::MCgenerate(Double_t *MCvect)
{
   MakeEvent();
   GetMCvect(MCvect);
   return(fMCwt);
}

////////////////////////////////////////////////////////////////////////////////
/// User method.
/// It provides the value of the integral calculated from the averages of the MC run
/// May be called after (or during) the MC run.

void TFoam::GetIntegMC(Double_t &mcResult, Double_t &mcError)
{
   Double_t mCerelat;
   mcResult = 0.0;
   mCerelat = 1.0;
   if (fNevGen>0) {
      mcResult = fPrime*fSumWt/fNevGen;
      mCerelat = sqrt( fSumWt2/(fSumWt*fSumWt) - 1/fNevGen);
   }
   mcError = mcResult *mCerelat;
}

////////////////////////////////////////////////////////////////////////////////
/// User method.
/// It returns NORMALIZATION integral to be combined with the average weights
/// and content of the histograms in order to get proper absolute normalization
/// of the integrand and distributions.
/// It can be called after initialization, before or during the MC run.

void  TFoam::GetIntNorm(Double_t& IntNorm, Double_t& Errel )
{
   if(fOptRej == 1) {    // Wt=1 events, internal rejection
      Double_t intMC,errMC;
      GetIntegMC(intMC,errMC);
      IntNorm = intMC;
      Errel   = errMC;
   } else {                // Wted events, NO internal rejection
      IntNorm = fPrime;
      Errel   = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// May be called optionally after the MC run.
/// Returns various parameters of the MC weight for efficiency evaluation

void  TFoam::GetWtParams(Double_t eps, Double_t &aveWt, Double_t &wtMax, Double_t &sigma)
{
   Double_t mCeff, wtLim;
   fMCMonit->GetMCeff(eps, mCeff, wtLim);
   wtMax = wtLim;
   aveWt = fSumWt/fNevGen;
   sigma = sqrt( fSumWt2/fNevGen -aveWt*aveWt );
}

////////////////////////////////////////////////////////////////////////////////
/// May be called optionally by the user after the MC run.
/// It provides normalization and also prints some information/statistics on the MC run.

void TFoam::Finalize(Double_t& IntNorm, Double_t& Errel)
{
   GetIntNorm(IntNorm,Errel);
   Double_t mcResult,mcError;
   GetIntegMC(mcResult,mcError);
   Double_t mCerelat= mcError/mcResult;
   //
   if(fChat>0) {
      Double_t eps = 0.0005;
      Double_t mCeff, mcEf2, wtMax, aveWt, sigma;
      GetWtParams(eps, aveWt, wtMax, sigma);
      mCeff=0;
      if(wtMax>0.0) mCeff=aveWt/wtMax;
      mcEf2 = sigma/aveWt;
      Double_t driver = fCells[0]->GetDriv();
      //
      BXOPE;
      BXTXT("****************************************");
      BXTXT("******     TFoam::Finalize       ******");
      BXTXT("****************************************");
      BX1I("    NevGen",fNevGen, "Number of generated events in the MC generation   ");
      BX1I("    nCalls",fNCalls, "Total number of function calls                    ");
      BXTXT("----------------------------------------");
      BX1F("     AveWt",aveWt,    "Average MC weight                      ");
      BX1F("     WtMin",fWtMin,  "Minimum MC weight (absolute)           ");
      BX1F("     WtMax",fWtMax,  "Maximum MC weight (absolute)           ");
      BXTXT("----------------------------------------");
      BX1F("    XPrime",fPrime,  "Primary total integral, R_prime        ");
      BX1F("    XDiver",driver,   "Driver  total integral, R_loss         ");
      BXTXT("----------------------------------------");
      BX2F("    IntMC", mcResult,  mcError,      "Result of the MC Integral");
      BX1F(" mCerelat", mCerelat,  "Relative error of the MC integral      ");
      BX1F(" <w>/WtMax",mCeff,     "MC efficiency, acceptance rate");
      BX1F(" Sigma/<w>",mcEf2,     "MC efficiency, variance/ave_wt");
      BX1F("     WtMax",wtMax,     "WtMax(esp= 0.0005)            ");
      BX1F("     Sigma",sigma,     "variance of MC weight         ");
      if(fOptRej==1) {
         Double_t avOve=fSumOve/fSumWt;
         BX1F("<OveW>/<W>",avOve,     "Contrib. of events wt>MaxWtRej");
      }
      BXCLO;
   }
}  // Finalize

////////////////////////////////////////////////////////////////////////////////
/// This can be called before Initialize, after setting kDim
/// It defines which variables are excluded in the process of the cell division.
/// For example 'FoamX->SetInhiDiv(1, 1);' inhibits division of y-variable.
/// The resulting map of cells in 2-dim. case will look as follows:
///
/// \image html foam_Map2.png width=400

void  TFoam::SetInhiDiv(Int_t iDim, Int_t InhiDiv)
{


   if(fDim==0) Error("TFoam","SetInhiDiv: fDim=0 \n");
   if(fInhiDiv == 0) {
      fInhiDiv = new Int_t[ fDim ];
      for(Int_t i=0; i<fDim; i++) fInhiDiv[i]=0;
   }
   //
   if( ( 0<=iDim) && (iDim<fDim)) {
      fInhiDiv[iDim] = InhiDiv;
   } else
      Error("SetInhiDiv:","Wrong iDim \n");
}

////////////////////////////////////////////////////////////////////////////////
/// This should be called before Initialize, after setting  kDim
/// It predefines values of the cell division for certain variable iDim.
/// For example setting 3 predefined division lines using:
///
/// ~~~ {.cpp}
///     xDiv[0]=0.30; xDiv[1]=0.40; xDiv[2]=0.65;
///     FoamX->SetXdivPRD(0,3,xDiv);
/// ~~~
///
/// results in the following 2-dim. pattern of the cells:
///
/// \image html foam_Map3.png width=400

void  TFoam::SetXdivPRD(Int_t iDim, Int_t len, Double_t xDiv[])
{
   Int_t i;

   if(fDim<=0)  Error("SetXdivPRD", "fDim=0 \n");
   if(   len<1 )  Error("SetXdivPRD", "len<1 \n");
   // allocate list of pointers, if it was not done before
   if(fXdivPRD == 0) {
      fXdivPRD = new TFoamVect*[fDim];
      for(i=0; i<fDim; i++)  fXdivPRD[i]=0;
   }
  // set division list for direction iDim in H-cubic space!!!
   if( ( 0<=iDim) && (iDim<fDim)) {
      fOptPRD =1;      // !!!!
      if( fXdivPRD[iDim] != 0)
         Error("SetXdivPRD", "Second allocation of XdivPRD not allowed \n");
      fXdivPRD[iDim] = new TFoamVect(len); // allocate list of division points
      for(i=0; i<len; i++) {
         (*fXdivPRD[iDim])[i]=xDiv[i]; // set list of division points
      }
   } else {
      Error("SetXdivPRD", "Wrong iDim  \n");
   }
   // Printing predefined division points
   std::cout<<" SetXdivPRD, idim= "<<iDim<<"  len= "<<len<<"   "<<std::endl;
   for(i=0; i<len; i++) {
      if (iDim < fDim) std::cout<< (*fXdivPRD[iDim])[i] <<"  ";
   }
   std::cout<<std::endl;
   for(i=0; i<len; i++)  std::cout<< xDiv[i] <<"   ";
   std::cout<<std::endl;
   //
}

////////////////////////////////////////////////////////////////////////////////
///  User utility, miscellaneous and debug.
///  Checks all pointers in the tree of cells. This is useful auto diagnostic.
/// - level=0, no printout, failures causes STOP
/// - level=1, printout, failures lead to WARNINGS only

void TFoam::CheckAll(Int_t level)
{
   Int_t errors, warnings;
   TFoamCell *cell;
   Long_t iCell;

   errors = 0; warnings = 0;
   if (level==1) std::cout << "///////////////////////////// FOAM_Checks /////////////////////////////////" << std::endl;
   for(iCell=1; iCell<=fLastCe; iCell++) {
      cell = fCells[iCell];
      //  checking general rules
      if( ((cell->GetDau0()==0) && (cell->GetDau1()!=0) ) ||
         ((cell->GetDau1()==0) && (cell->GetDau0()!=0) ) ) {
         errors++;
         if (level==1) Error("CheckAll","ERROR: Cell's no %ld has only one daughter \n",iCell);
      }
      if( (cell->GetDau0()==0) && (cell->GetDau1()==0) && (cell->GetStat()==0) ) {
         errors++;
         if (level==1) Error("CheckAll","ERROR: Cell's no %ld  has no daughter and is inactive \n",iCell);
      }
      if( (cell->GetDau0()!=0) && (cell->GetDau1()!=0) && (cell->GetStat()==1) ) {
         errors++;
         if (level==1) Error("CheckAll","ERROR: Cell's no %ld has two daughters and is active \n",iCell);
      }

      // checking parents
      if( (cell->GetPare())!=fCells[0] ) { // not child of the root
         if ( (cell != cell->GetPare()->GetDau0()) && (cell != cell->GetPare()->GetDau1()) ) {
            errors++;
            if (level==1) Error("CheckAll","ERROR: Cell's no %ld parent not pointing to this cell\n ",iCell);
         }
      }

      // checking daughters
      if(cell->GetDau0()!=0) {
         if(cell != (cell->GetDau0())->GetPare()) {
            errors++;
            if (level==1)  Error("CheckAll","ERROR: Cell's no %ld daughter 0 not pointing to this cell \n",iCell);
         }
      }
      if(cell->GetDau1()!=0) {
         if(cell != (cell->GetDau1())->GetPare()) {
            errors++;
            if (level==1) Error("CheckAll","ERROR: Cell's no %ld daughter 1 not pointing to this cell \n",iCell);
         }
      }
   }// loop after cells;

   // Check for empty cells
   for(iCell=0; iCell<=fLastCe; iCell++) {
      cell = fCells[iCell];
      if( (cell->GetStat()==1) && (cell->GetDriv()==0) ) {
         warnings++;
         if(level==1) Warning("CheckAll", "Warning: Cell no. %ld is active but empty \n", iCell);
      }
   }
   // summary
   if(level==1){
      Info("CheckAll","Check has found %d errors and %d warnings \n",errors, warnings);
   }
   if(errors>0){
      Info("CheckAll","Check - found total %d  errors \n",errors);
   }
} // Check

////////////////////////////////////////////////////////////////////////////////
/// Prints geometry of ALL cells of the FOAM

void TFoam::PrintCells(void)
{
   Long_t iCell;

   for(iCell=0; iCell<=fLastCe; iCell++) {
      std::cout<<"Cell["<<iCell<<"]={ ";
      //std::cout<<"  "<< fCells[iCell]<<"  ";  // extra DEBUG
      std::cout<<std::endl;
      fCells[iCell]->Print("1");
      std::cout<<"}"<<std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Debugging tool which plots 2-dimensional cells as rectangles
/// in C++ format readable for root

void TFoam::RootPlot2dim(Char_t *filename)
{
   std::ofstream outfile(filename, std::ios::out);
   Double_t   x1,y1,x2,y2,x,y;
   Long_t    iCell;
   Double_t offs =0.1;
   Double_t lpag   =1-2*offs;
   outfile<<"{" << std::endl;
   outfile<<"cMap = new TCanvas(\"Map1\",\" Cell Map \",600,600);"<<std::endl;
   //
   outfile<<"TBox*a=new TBox();"<<std::endl;
   outfile<<"a->SetFillStyle(0);"<<std::endl;  // big frame
   outfile<<"a->SetLineWidth(4);"<<std::endl;
   outfile<<"a->SetLineColor(2);"<<std::endl;
   outfile<<"a->DrawBox("<<offs<<","<<offs<<","<<(offs+lpag)<<","<<(offs+lpag)<<");"<<std::endl;
   //
   outfile<<"TText*t=new TText();"<<std::endl;  // text for numbering
   outfile<<"t->SetTextColor(4);"<<std::endl;
   if(fLastCe<51)
      outfile<<"t->SetTextSize(0.025);"<<std::endl;  // text for numbering
   else if(fLastCe<251)
      outfile<<"t->SetTextSize(0.015);"<<std::endl;
   else
      outfile<<"t->SetTextSize(0.008);"<<std::endl;
   //
   outfile<<"TBox*b=new TBox();"<<std::endl;  // single cell
   outfile <<"b->SetFillStyle(0);"<<std::endl;
   //
   if(fDim==2 && fLastCe<=2000) {
      TFoamVect  cellPosi(fDim); TFoamVect  cellSize(fDim);
      outfile << "// =========== Rectangular cells  ==========="<< std::endl;
      for(iCell=1; iCell<=fLastCe; iCell++) {
         if( fCells[iCell]->GetStat() == 1) {
            fCells[iCell]->GetHcub(cellPosi,cellSize);
            x1 = offs+lpag*(        cellPosi[0]); y1 = offs+lpag*(        cellPosi[1]);
            x2 = offs+lpag*(cellPosi[0]+cellSize[0]); y2 = offs+lpag*(cellPosi[1]+cellSize[1]);
            //     cell rectangle
            if(fLastCe<=2000)
            outfile<<"b->DrawBox("<<x1<<","<<y1<<","<<x2<<","<<y2<<");"<<std::endl;
            //     cell number
            if(fLastCe<=250) {
               x = offs+lpag*(cellPosi[0]+0.5*cellSize[0]); y = offs+lpag*(cellPosi[1]+0.5*cellSize[1]);
               outfile<<"t->DrawText("<<x<<","<<y<<","<<"\""<<iCell<<"\""<<");"<<std::endl;
            }
         }
      }
      outfile<<"// ============== End Rectangles ==========="<< std::endl;
   }//kDim=2
   //
   //
   outfile << "}" << std::endl;
   outfile.close();
}

void TFoam::LinkCells()
{
// Void function for backward compatibility

   Info("LinkCells", "VOID function for backward compatibility \n");
   return;
}

