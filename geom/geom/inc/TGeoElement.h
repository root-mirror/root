// @(#)root/geom:$Id$
// Author: Andrei Gheata   17/06/04
// Added support for radionuclides: Mihaela Gheata 24/08/2006
/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoElement
#define ROOT_TGeoElement

#include "TNamed.h"

#include "TAttLine.h"

#include "TAttFill.h"

#include "TAttMarker.h"

#include "TObjArray.h"

#include <map>

class TGeoElementTable;
class TGeoIsotope;

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoElement - a chemical element                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoElement : public TNamed
{
protected:
   enum EGeoElement {
      kElemUsed    =   BIT(17),
      kElemDefined =   BIT(18),
      kElementChecked = BIT(19)
   };

   Int_t                    fZ;          // Z of element
   Int_t                    fN;          // Number of nucleons
   Int_t                    fNisotopes;  // Number of isotopes for the element
   Double_t                 fA;          // A of element
   TObjArray               *fIsotopes;   // List of isotopes
   Double_t                *fAbundances; //[fNisotopes] Array of relative isotope abundances
   Double_t                 fCoulomb;    // Coulomb correction factor
   Double_t                 fRadTsai;    // Tsai formula for the radiation length

private:
   TGeoElement(const TGeoElement &other) = delete;
   TGeoElement &operator=(const TGeoElement &other) = delete;

   // Compute the Coulomb correction factor
   void                     ComputeCoulombFactor();
   // Compute the Tsai formula for the radiation length
   void                     ComputeLradTsaiFactor();

public:
   // constructors
   TGeoElement();
   TGeoElement(const char *name, const char *title, Int_t z, Double_t a);
   TGeoElement(const char *name, const char *title, Int_t nisotopes);
   TGeoElement(const char *name, const char *title, Int_t z, Int_t n, Double_t a);
   // destructor
   virtual ~TGeoElement();
   // methods
   virtual Int_t            ENDFCode()    const { return 0;}
   Int_t                    Z() const {return fZ;}
   Int_t                    N() const {return fN;}
   Double_t                 Neff() const;
   Double_t                 A() const {return fA;}
   void                     AddIsotope(TGeoIsotope *isotope, Double_t relativeAbundance);
   Int_t                    GetNisotopes() const {return fNisotopes;}
   TGeoIsotope             *GetIsotope(Int_t i) const;
   Double_t                 GetRelativeAbundance(Int_t i) const;
   // Calculate properties for an atomic number
   void                     ComputeDerivedQuantities();
   // Specific activity (in Bq/gram)
   virtual Double_t         GetSpecificActivity() const {return 0.;}
   Bool_t                   HasIsotopes() const {return (fNisotopes==0)?kFALSE:kTRUE;}
   Bool_t                   IsDefined() const {return TObject::TestBit(kElemDefined);}
   virtual Bool_t           IsRadioNuclide() const {return kFALSE;}
   Bool_t                   IsUsed() const {return TObject::TestBit(kElemUsed);}
   virtual void             Print(Option_t *option = "") const;
   void                     SetDefined(Bool_t flag=kTRUE) {TObject::SetBit(kElemDefined,flag);}
   void                     SetUsed(Bool_t flag=kTRUE) {TObject::SetBit(kElemUsed,flag);}
   static TGeoElementTable *GetElementTable();
   // Coulomb correction factor
   inline Double_t          GetfCoulomb() const {return fCoulomb;}
   // Tsai formula for the radiation length
   inline Double_t          GetfRadTsai() const {return fRadTsai;}

   ClassDef(TGeoElement, 3)              // base element class
};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoIsotope - a isotope defined by the atomic number, number of        //
// nucleons and atomic weight (g/mole)                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoIsotope : public TNamed
{
protected:
   Int_t                    fZ;           // atomic number
   Int_t                    fN;           // number of nucleons
   Double_t                 fA;           // atomic mass (g/mole)

public:
   TGeoIsotope();
   TGeoIsotope(const char *name, Int_t z, Int_t n, Double_t a);
   virtual ~TGeoIsotope() {}

   Int_t                    GetZ() const {return fZ;}
   Int_t                    GetN() const {return fN;}
   Double_t                 GetA() const {return fA;}
   static TGeoIsotope      *FindIsotope(const char *name);
   virtual void             Print(Option_t *option = "") const;

   ClassDef(TGeoIsotope, 1)              // Isotope class defined by Z,N,A
};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoElementRN - a radionuclide.                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoDecayChannel;
class TGeoBatemanSol;

class TGeoElementRN : public TGeoElement
{
protected:
   Int_t                    fENDFcode; // ENDF element code
   Int_t                    fIso;      // Isomer number
   Double_t                 fLevel;    // Isomeric level
   Double_t                 fDeltaM;   // Mass excess
   Double_t                 fHalfLife; // Half life
   Double_t                 fNatAbun;  // Natural Abundance
//   char                     fJP[11];   // Spin-parity
   Double_t                 fTH_F;     // Hynalation toxicity
   Double_t                 fTG_F;     // Ingestion toxicity
   Double_t                 fTH_S;     // Hynalation toxicity
   Double_t                 fTG_S;     // Ingestion toxicity
   Int_t                    fStatus;   // Status code
   TGeoBatemanSol          *fRatio;    // Time evolution of proportion by number

   TObjArray               *fDecays;   // List of decay modes

   void                     MakeName(Int_t a, Int_t z, Int_t iso);

private:
   TGeoElementRN(const TGeoElementRN& elem) = delete;
   TGeoElementRN& operator=(const TGeoElementRN& elem) = delete;

public:
   TGeoElementRN();
   TGeoElementRN(Int_t A, Int_t Z, Int_t iso, Double_t level,
         Double_t deltaM, Double_t halfLife, const char* JP,
         Double_t natAbun, Double_t th_f, Double_t tg_f, Double_t th_s,
         Double_t tg_s, Int_t status);
   virtual ~TGeoElementRN();

   void                     AddDecay(Int_t decay, Int_t diso, Double_t branchingRatio, Double_t qValue);
   void                     AddDecay(TGeoDecayChannel *dc);
   void                     AddRatio(TGeoBatemanSol &ratio);
   void                     ResetRatio();
   static Int_t             ENDF(Int_t a, Int_t z, Int_t iso) {return 10000*z+10*a+iso;}

   // Getters
   virtual Int_t            ENDFCode()    const {return fENDFcode;}
   virtual Double_t         GetSpecificActivity() const;
   virtual Bool_t           IsRadioNuclide() const {return kTRUE;}
   Int_t                    MassNo()      const {return (Int_t)fA;}
   Int_t                    AtomicNo()    const {return fZ;}
   Int_t                    IsoNo()       const {return fIso;}
   Double_t                 Level()       const {return fLevel;}
   Double_t                 MassEx()      const {return fDeltaM;}
   Double_t                 HalfLife()    const {return fHalfLife;}
   Double_t                 NatAbun()     const {return fNatAbun;}
   const char*              PJ()          const {return fTitle.Data();}
   Double_t                 TH_F()        const {return fTH_F;}
   Double_t                 TG_F()        const {return fTG_F;}
   Double_t                 TH_S()        const {return fTH_S;}
   Double_t                 TG_S()        const {return fTG_S;}
   Double_t                 Status()      const {return fStatus;}
   Bool_t                   Stable()      const {return !fDecays;}
   TObjArray               *Decays()      const {return fDecays;}
   Int_t                    GetNdecays()  const;
   TGeoBatemanSol          *Ratio()       const {return fRatio;}

   // Utilities
   Bool_t                   CheckDecays() const;
   Int_t                    DecayResult(TGeoDecayChannel *dc) const;
   void                     FillPopulation(TObjArray *population, Double_t precision=0.001, Double_t factor=1.);
   virtual void             Print(Option_t *option = "") const;
   static TGeoElementRN    *ReadElementRN(const char *record, Int_t &ndecays);
   virtual void             SavePrimitive(std::ostream &out, Option_t *option = "");

   ClassDef(TGeoElementRN, 2)           // radionuclides class
};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoDecayChannel - decay channel utility class.                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoDecayChannel : public TObject
{
private:
   UInt_t                   fDecay;          // Decay mode
   Int_t                    fDiso;           // Delta isomeric number
   Double_t                 fBranchingRatio; // Branching Ratio
   Double_t                 fQvalue;         // Qvalue in GeV
   TGeoElementRN           *fParent;         // Parent element
   TGeoElementRN           *fDaughter;       // Daughter element
public:
   enum ENuclearDecayMode {
      kBitMask32  = 0xffffffff,
      k2BetaMinus   = BIT(0),
      kBetaMinus    = BIT(1),
      kNeutronEm    = BIT(2),
      kProtonEm     = BIT(3),
      kAlpha        = BIT(4),
      kECF          = BIT(5),
      kElecCapt     = BIT(6),
      kIsoTrans     = BIT(7),
      kI            = BIT(8),
      kSpontFiss    = BIT(9),
      k2P           = BIT(10),
      k2N           = BIT(11),
      k2A           = BIT(12),
      kCarbon12     = BIT(13),
      kCarbon14     = BIT(14)
   };
   TGeoDecayChannel() : fDecay(0), fDiso(0), fBranchingRatio(0), fQvalue(0), fParent(0), fDaughter(0) {}
   TGeoDecayChannel(Int_t decay, Int_t diso, Double_t branchingRatio, Double_t qValue)
                  : fDecay(decay), fDiso(diso), fBranchingRatio(branchingRatio), fQvalue(qValue), fParent(0), fDaughter(0) {}
   TGeoDecayChannel(const TGeoDecayChannel &dc) : TObject(dc),fDecay(dc.fDecay),fDiso(dc.fDiso),fBranchingRatio(dc.fBranchingRatio),
                                                  fQvalue(dc.fQvalue),fParent(dc.fParent),fDaughter(dc.fDaughter) {}
   virtual ~TGeoDecayChannel() {}

   TGeoDecayChannel& operator=(const TGeoDecayChannel& dc);

   // Getters
   Int_t                    GetIndex()       const;
   virtual const char      *GetName()        const;
   UInt_t                   Decay()          const {return fDecay;}
   Double_t                 BranchingRatio() const {return fBranchingRatio;}
   Double_t                 Qvalue()         const {return fQvalue;}
   Int_t                    DeltaIso()       const {return fDiso;}
   TGeoElementRN           *Daughter()       const {return fDaughter;}
   TGeoElementRN           *Parent()         const {return fParent;}
   static void              DecayName(UInt_t decay, TString &name);
   // Setters
   void                     SetParent(TGeoElementRN *parent) {fParent = parent;}
   void                     SetDaughter(TGeoElementRN *daughter) {fDaughter = daughter;}
   // Services
   virtual void             Print(Option_t *opt = " ") const;
   static TGeoDecayChannel *ReadDecay(const char *record);
   virtual void             SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void             DecayShift(Int_t &dA, Int_t &dZ, Int_t &dI) const ;

   ClassDef(TGeoDecayChannel,1)    // Decay channel for Elements
};

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// TGeoBatemanSol -Class representing the Bateman solution for a decay branch //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

class TGeoBatemanSol : public TObject, public TAttLine, public TAttFill, public TAttMarker
{
private:
   typedef struct {
      Double_t   cn;     // Concentration for element 'i': Ni/Ntop
      Double_t   lambda; // Decay coef. for element 'i'
   } BtCoef_t;
   TGeoElementRN           *fElem;            // Referred RN element
   TGeoElementRN           *fElemTop;         // Top RN element
   Int_t                    fCsize;           // Size of the array of coefficients
   Int_t                    fNcoeff;          // Number of coefficients
   Double_t                 fFactor;          // Constant factor that applies to all coefficients
   Double_t                 fTmin;            // Minimum value of the time interval
   Double_t                 fTmax;            // Maximum value of the time interval
   BtCoef_t                *fCoeff;           //[fNcoeff] Array of coefficients
public:
   TGeoBatemanSol() : TObject(), TAttLine(), TAttFill(), TAttMarker(), fElem(NULL), fElemTop(NULL), fCsize(0), fNcoeff(0), fFactor(1.), fTmin(0.), fTmax(0), fCoeff(NULL) {}
   TGeoBatemanSol(TGeoElementRN *elem);
   TGeoBatemanSol(const TObjArray *chain);
   TGeoBatemanSol(const TGeoBatemanSol& other);
   ~TGeoBatemanSol();

   TGeoBatemanSol& operator=(const TGeoBatemanSol& other);
   TGeoBatemanSol& operator+=(const TGeoBatemanSol& other);

   Double_t                 Concentration(Double_t time) const;
   virtual void             Draw(Option_t *option="");
   void                     GetCoeff(Int_t i, Double_t &cn, Double_t &lambda) const {cn=fCoeff[i].cn; lambda=fCoeff[i].lambda;}
   void                     GetRange(Double_t &tmin, Double_t &tmax) const {tmin=fTmin; tmax=fTmax;}
   TGeoElementRN           *GetElement()    const {return fElem;}
   TGeoElementRN           *GetTopElement() const {return fElemTop;}
   Int_t                    GetNcoeff()     const  {return fNcoeff;}
   virtual void             Print(Option_t *option = "") const;
   void                     SetRange(Double_t tmin=0., Double_t tmax=0.) {fTmin=tmin; fTmax=tmax;}
   void                     SetFactor(Double_t factor) {fFactor = factor;}
   void                     FindSolution(const TObjArray *array);
   void                     Normalize(Double_t factor);

   ClassDef(TGeoBatemanSol,1)       // Solution for the Bateman equation
};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoElemIter - iterator for decay chains.                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoElemIter
{
private:
   const TGeoElementRN     *fTop;            // Top element of the iteration
   const TGeoElementRN     *fElem;           // Current element
   TObjArray               *fBranch;         // Current branch
   Int_t                    fLevel;          // Current level
   Double_t                 fLimitRatio;     // Minimum cumulative branching ratio
   Double_t                 fRatio;          // Current ratio

protected:
   TGeoElemIter() : fTop(0), fElem(0), fBranch(0), fLevel(0), fLimitRatio(0), fRatio(0) {}
   TGeoElementRN           *Down(Int_t ibranch);
   TGeoElementRN           *Up();

public:
   TGeoElemIter(TGeoElementRN *top, Double_t limit=1.e-4);
   TGeoElemIter(const TGeoElemIter &iter);
   virtual ~TGeoElemIter();

   TGeoElemIter   &operator=(const TGeoElemIter &iter);
   TGeoElementRN  *operator()();
   TGeoElementRN           *Next();

   TObjArray               *GetBranch() const              {return fBranch;}
   const TGeoElementRN     *GetTop() const                 {return fTop;}
   const TGeoElementRN     *GetElement() const             {return fElem;}
   Int_t                    GetLevel() const               {return fLevel;}
   Double_t                 GetRatio() const               {return fRatio;}
   virtual void             Print(Option_t *option="") const;
   void                     SetLimitRatio(Double_t limit)  {fLimitRatio = limit;}

   ClassDef(TGeoElemIter,0)    // Iterator for radionuclide chains.
};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoElementTable - table of elements                                   //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoElementTable : public TObject
{
private:
// data members
   Int_t                    fNelements;    // number of elements
   Int_t                    fNelementsRN;  // number of RN elements
   Int_t                    fNisotopes;    // number of isotopes
   TObjArray               *fList;         // list of elements
   TObjArray               *fListRN;       // list of RN elements
   TObjArray               *fIsotopes;     // list of user-defined isotopes
   // Map of radionuclides
   typedef std::map<Int_t, TGeoElementRN *>   ElementRNMap_t;
   typedef ElementRNMap_t::iterator           ElementRNMapIt_t;
   ElementRNMap_t           fElementsRN; //! map of RN elements with ENDF key

protected:
   TGeoElementTable(const TGeoElementTable&);
   TGeoElementTable& operator=(const TGeoElementTable&);

public:
   // constructors
   TGeoElementTable();
   TGeoElementTable(Int_t nelements);
   // destructor
   virtual ~TGeoElementTable();
   // methods

   enum EGeoETStatus {
      kETDefaultElements = BIT(14),
      kETRNElements      = BIT(15)
   };
   void                     AddElement(const char *name, const char *title, Int_t z, Double_t a);
   void                     AddElement(const char *name, const char *title, Int_t z, Int_t n, Double_t a);
   void                     AddElement(TGeoElement *elem);
   void                     AddElementRN(TGeoElementRN *elem);
   void                     AddIsotope(TGeoIsotope *isotope);
   void                     BuildDefaultElements();
   void                     ImportElementsRN();
   Bool_t                   CheckTable() const;
   TGeoElement             *FindElement(const char *name) const;
   TGeoIsotope             *FindIsotope(const char *name) const;
   TGeoElement             *GetElement(Int_t z) {return (TGeoElement*)fList->At(z);}
   TGeoElementRN           *GetElementRN(Int_t ENDFcode) const;
   TGeoElementRN           *GetElementRN(Int_t a, Int_t z, Int_t iso=0) const;
   TObjArray               *GetElementsRN() const {return fListRN;}
   Bool_t                   HasDefaultElements() const {return TObject::TestBit(kETDefaultElements);}
   Bool_t                   HasRNElements() const {return TObject::TestBit(kETRNElements);}

   Int_t                    GetNelements() const {return fNelements;}
   Int_t                    GetNelementsRN() const {return fNelementsRN;}
   void                     ExportElementsRN(const char *filename="");
   virtual void             Print(Option_t *option = "") const;

   ClassDef(TGeoElementTable,4)              // table of elements
};

#endif

