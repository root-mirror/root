// @(#)root/foam:$Id$
// Author: S. Jadach <mailto:Stanislaw.jadach@ifj.edu.pl>, P.Sawicki <mailto:Pawel.Sawicki@ifj.edu.pl>

#ifndef ROOT_TFoam
#define ROOT_TFoam


#include "TObject.h"
#include "TString.h"

#include <vector>

class TH1D;
class TMethodCall;
class TRandom;
class TFoamIntegrand;
class TFoamMaxwt;
class TFoamVect;
class TFoamCell;

class TFoam : public TObject {
protected:

   TString fName;             ///< Name of a given instance of the FOAM class
   TString fVersion;          ///< Actual version of the FOAM like (1.01m)
   TString fDate;             ///< Release date of FOAM
   Int_t   fDim;              ///< Dimension of the integration/simulation space
   Int_t   fNCells;           ///< Maximum number of cells
   Int_t   fRNmax;            ///< Maximum No. of the rand. numb. requested at once

   Int_t   fOptDrive;         ///< Optimization switch =1,2 for variance or maximum weight optimization
   Int_t   fChat;             ///< Chat=0,1,2 chat level in output, Chat=1 normal level
   Int_t   fOptRej;           ///< Switch =0 for weighted events; =1 for unweighted events in MC

   Int_t   fNBin;             ///< No. of bins in the edge histogram for cell MC exploration
   Int_t   fNSampl;           ///< No. of MC events, when dividing (exploring) cell
   Int_t   fEvPerBin;         ///< Maximum number of effective (wt=1) events per bin

   Int_t  *fMaskDiv;          ///<! [fDim] Dynamic Mask for  cell division
   Int_t  *fInhiDiv;          ///<! [fDim] Flags for inhibiting cell division
   Int_t   fOptPRD;           ///<  Option switch for predefined division, for quick check
   TFoamVect **fXdivPRD;      ///<! Lists of division values encoded in one vector per direction

   Int_t   fNoAct;            ///< Number of active cells
   Int_t   fLastCe;           ///< Index of the last cell
   TFoamCell **fCells;        ///< [fNCells] Array of ALL cells

   TFoamMaxwt   *fMCMonit;    ///< Monitor of the MC weight for measuring MC efficiency
   Double_t   fMaxWtRej;      ///< Maximum weight in rejection for getting wt=1 events
   std::vector<Long_t> fCellsAct; ///< Index of active cells, constructed at the end of foam build-up
   Double_t  *fPrimAcu;       ///< [fNoAct] Array of cumulative probability of all active cells
   TObjArray *fHistEdg;       ///< Histograms of wt, one for each cell edge
   TObjArray *fHistDbg;       ///< Histograms of wt, for debug
   TH1D      *fHistWt;        ///< Histogram of the MC wt

   Double_t *fMCvect;         ///< [fDim] Generated MC vector for the outside user
   Double_t  fMCwt;           ///< MC weight
   Double_t *fRvec;           ///< [fRNmax] random number vector from r.n. generator fDim+1 maximum elements

   TFoamIntegrand *fRho;      ///<! Pointer to the user-defined integrand function/distribution
   TMethodCall *fMethodCall;  ///<! ROOT's pointer to user-defined global distribution function
   TRandom         *fPseRan;  ///< Pointer to user-defined generator of pseudorandom numbers

   Long_t   fNCalls;          ///< Total number of the function calls
   Long_t   fNEffev;          ///< Total number of effective events (wt=1) in the foam buildup
   Double_t fSumWt, fSumWt2;  ///< Total sum of wt and wt^2
   Double_t fSumOve;          ///< Total Sum of overweighted events
   Double_t fNevGen;          ///< Total number of the generated MC events
   Double_t fWtMax, fWtMin;   ///< Maximum/Minimum MC weight
   Double_t fPrime;           ///< Primary integral R' (R=R'<wt>)
   Double_t fMCresult;        ///< True Integral R from MC series
   Double_t fMCerror;         ///< and its error

   Double_t *fAlpha;          ///< [fDim] Internal parameters of the hyper-rectangle

public:
   TFoam();                          // Default constructor (used only by ROOT streamer)
   TFoam(const Char_t*);             // Principal user-defined constructor
   virtual ~TFoam();                 // Default destructor
   TFoam(const TFoam&);              // Copy Constructor  NOT USED
   // Initialization
   virtual void Initialize();                // Initialization of the FOAM (grid, cells, etc), mandatory!
   virtual void Initialize(TRandom *, TFoamIntegrand *); // Alternative initialization method, backward compatibility
   virtual void InitCells();                 // Initializes first cells inside original cube
   virtual Int_t  CellFill(Int_t, TFoamCell*);  // Allocates new empty cell and return its index
   virtual void Explore(TFoamCell *Cell);       // Exploration of the new cell, determine <wt>, wtMax etc.
   virtual void Carver(Int_t&,Double_t&,Double_t&);// Determines the best edge, wt_max reduction
   virtual void Varedu(Double_t [], Int_t&, Double_t&,Double_t&); // Determines the best edge, variance reduction
   virtual void MakeAlpha();                 // Provides random point inside hyper-rectangle
   virtual void Grow();                      // Adds new cells to FOAM object until buffer is full
   virtual Long_t PeekMax();                 // Choose one active cell, used by Grow and also in MC generation
   virtual Int_t  Divide(TFoamCell *);       // Divide iCell into two daughters; iCell retained, taged as inactive
   virtual void MakeActiveList();            // Creates table of active cells
   virtual void GenerCel2(TFoamCell *&);     // Chose an active cell the with probability ~ Primary integral
   // Generation
   virtual Double_t Eval(Double_t *);        // Evaluates value of the distribution function
   virtual void     MakeEvent();             // Makes (generates) single MC event
   virtual void     GetMCvect(Double_t *);   // Provides generated randomly MC vector
   virtual void     GetMCwt(Double_t &);     // Provides generated MC weight
   virtual Double_t GetMCwt();               // Provides generates MC weight
   virtual Double_t MCgenerate(Double_t *MCvect);// All three above function in one
   // Finalization
   virtual void GetIntegMC(Double_t&, Double_t&);// Provides Integrand and abs. error from MC run
   virtual void GetIntNorm(Double_t&, Double_t&);// Provides normalization Inegrand
   virtual void GetWtParams(Double_t, Double_t&, Double_t&, Double_t&);// Provides MC weight parameters
   virtual void Finalize(  Double_t&, Double_t&);  // Prints summary of MC integration
   virtual TFoamIntegrand  *GetRho(){return fRho;} // Gets pointer of the distribut. (after restoring from disk)
   virtual TRandom *GetPseRan() const {return fPseRan;}   // Gets pointer of r.n. generator (after restoring from disk)
   virtual void SetRhoInt(Double_t (*fun)(Int_t, Double_t *));    // Set new integrand distr. in compiled mode
   virtual void SetRho(TFoamIntegrand *Rho);       // Set new integrand distr. in compiled mode
   virtual void ResetRho(TFoamIntegrand *Rho);                // Set new distribution, delete old
   virtual void SetPseRan(TRandom *PseRan){fPseRan=PseRan;}   // Set new r.n. generator
   virtual void ResetPseRan(TRandom *PseRan);                 // Set new r.n.g, delete old
   // Getters and Setters
   virtual void SetkDim(Int_t kDim){fDim = kDim;}            // Sets dimension of cubical space
   virtual void SetnCells(Long_t nCells){fNCells =nCells;}  // Sets maximum number of cells
   virtual void SetnSampl(Long_t nSampl){fNSampl =nSampl;}  // Sets no of MC events in cell exploration
   virtual void SetnBin(Int_t nBin){fNBin = nBin;}          // Sets no of bins in histogs in cell exploration
   virtual void SetChat(Int_t Chat){fChat = Chat;}          // Sets option Chat, chat level
   virtual void SetOptRej(Int_t OptRej){fOptRej =OptRej;}   // Sets option for MC rejection
   virtual void SetOptDrive(Int_t OptDrive){fOptDrive =OptDrive;}  // Sets optimization switch
   virtual void SetEvPerBin(Int_t EvPerBin){fEvPerBin =EvPerBin;}  // Sets max. no. of effective events per bin
   virtual void SetMaxWtRej(Double_t MaxWtRej){fMaxWtRej=MaxWtRej;}  // Sets max. weight for rejection
   virtual void SetInhiDiv(Int_t, Int_t );            // Set inhibition of cell division along certain edge
   virtual void SetXdivPRD(Int_t, Int_t, Double_t[]); // Set predefined division points
   // Getters and Setters
   virtual const char *GetVersion() const {return fVersion.Data();}// Get version of the FOAM
   virtual Int_t    GetTotDim() const { return fDim;}              // Get total dimension
   virtual Double_t GetPrimary() const {return fPrime;}            // Get value of primary integral R'
   virtual void GetPrimary(Double_t &prime) {prime = fPrime;}      // Get value of primary integral R'
   virtual Long_t GetnCalls() const {return fNCalls;}            // Get total no. of the function calls
   virtual Long_t GetnEffev() const {return fNEffev;}            // Get total no. of effective wt=1 events
   // Debug
   virtual void CheckAll(Int_t);     // Checks correctness of the entire data structure in the FOAM object
   virtual void PrintCells();        // Prints content of all cells
   virtual void RootPlot2dim(Char_t*);   // Generates C++ code for drawing foam
   virtual void LinkCells(void);         // Void function for backward compatibility
   // Inline
private:
   Double_t Sqr(Double_t x) const { return x*x;}      // Square function

   ClassDef(TFoam,2);   // General purpose self-adapting Monte Carlo event generator
};

#endif
