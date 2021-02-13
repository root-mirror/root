// @(#)root/tree:$Id$
// Author: G Ganis Sep 2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TQueryResult
#define ROOT_TQueryResult


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TQueryResult                                                         //
//                                                                      //
// A container class for the results of a query.                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TDatime.h"
#include "TString.h"

#ifdef R__LESS_INCLUDES
class TMacro;
#else
#include "TMacro.h"
#endif

class TBrowser;
class TTreePlayer;
class TQueryResult;

Bool_t operator==(const TQueryResult &qr1, const TQueryResult &qr2);


class TQueryResult : public TNamed {

friend class TTreePlayer;
friend class TProofPlayerLite;
friend class TProofPlayerRemote;
friend class TProof;
friend class TProofLite;
friend class TProofServ;
friend class TQueryResultManager;

public:
   enum EQueryStatus {
      kAborted = 0, kSubmitted, kRunning, kStopped, kCompleted
   };

protected:
   Int_t           fSeqNum;       ///< query unique sequential number
   Bool_t          fDraw;         ///< true if draw action query
   EQueryStatus    fStatus;       ///< query status
   TDatime         fStart;        ///< time when processing started
   TDatime         fEnd;          ///< time when processing ended
   Float_t         fUsedCPU;      ///< real CPU time used (seconds)
   TString         fOptions;      ///< processing options + aclic mode (< opt >#< aclic_mode >)
   TList          *fInputList;    ///< input list; contains also data sets, entry list, ...
   Long64_t        fEntries;      ///< number of entries processed
   Long64_t        fFirst;        ///< first entry processed
   Long64_t        fBytes;        ///< number of bytes processed
   TMacro         *fLogFile;      ///< file with log messages from the query
   TMacro         *fSelecHdr;     ///< selector header file
   TMacro         *fSelecImp;     ///< selector implementation file
   TString         fLibList;      ///< blank-separated list of libs loaded at fStart
   TString         fParList;      ///< colon-separated list of PAR loaded at fStart
   TList          *fOutputList;   ///< output list
   Bool_t          fFinalized;    ///< whether Terminate has been run
   Bool_t          fArchived;     ///< whether the query has been archived
   TString         fResultFile;   ///< URL of the file where results have been archived
   Float_t         fPrepTime;     ///< Prepare time (seconds) (millisec precision)
   Float_t         fInitTime;     ///< Initialization time (seconds) (millisec precision)
   Float_t         fProcTime;     ///< Processing time (seconds) (millisec precision)
   Float_t         fMergeTime;    ///< Merging time (seconds) (millisec precision)
   Float_t         fRecvTime;     ///< Transfer-to-client time (seconds) (millisec precision)
   Float_t         fTermTime;     ///< Terminate time (seconds) (millisec precision)
   Int_t           fNumWrks;      ///< Number of workers at start
   Int_t           fNumMergers;   ///< Number of submergers

   TQueryResult(Int_t seqnum, const char *opt, TList *inlist,
                Long64_t entries, Long64_t first,
                const char *selec);

   void            AddInput(TObject *obj);
   void            AddLogLine(const char *logline);
   TQueryResult   *CloneInfo();
   virtual void    RecordEnd(EQueryStatus status, TList *outlist = 0);
   void            SaveSelector(const char *selec);
   void            SetArchived(const char *archfile);
   virtual void    SetFinalized() { fFinalized = kTRUE; }
   virtual void    SetInputList(TList *in, Bool_t adopt = kTRUE);
   virtual void    SetOutputList(TList *out, Bool_t adopt = kTRUE);
   virtual void    SetProcessInfo(Long64_t ent, Float_t cpu = 0.,
                                  Long64_t siz = -1,
                                  Float_t inittime = 0., Float_t proctime = 0.);
   void            SetPrepTime(Float_t preptime) { fPrepTime = preptime; }
   void            SetMergeTime(Float_t mergetime) { fMergeTime = mergetime; }
   void            SetRecvTime(Float_t recvtime) { fRecvTime = recvtime; }
   void            SetTermTime(Float_t termtime) { fTermTime = termtime; }
   void            SetNumMergers(Int_t nmergers) { fNumMergers = nmergers; }

public:
   TQueryResult() : fSeqNum(-1), fDraw(0), fStatus(kSubmitted), fUsedCPU(0.),
                    fInputList(0), fEntries(-1), fFirst(-1), fBytes(0),
                    fLogFile(0), fSelecHdr(0), fSelecImp(0),
                    fLibList("-"), fOutputList(0),
                    fFinalized(kFALSE), fArchived(kFALSE), fPrepTime(0.),
                    fInitTime(0.), fProcTime(0.), fMergeTime(0.),
                    fRecvTime(-1), fTermTime(0.), fNumWrks(-1), fNumMergers(-1) { }
   virtual ~TQueryResult();

   void           Browse(TBrowser *b = 0);

   Int_t          GetSeqNum() const { return fSeqNum; }
   EQueryStatus   GetStatus() const { return fStatus; }
   TDatime        GetStartTime() const { return fStart; }
   TDatime        GetEndTime() const { return fEnd; }
   const char    *GetOptions() const { return fOptions; }
   TList         *GetInputList() { return fInputList; }
   TObject       *GetInputObject(const char *classname) const;
   Long64_t       GetEntries() const { return fEntries; }
   Long64_t       GetFirst() const { return fFirst; }
   Long64_t       GetBytes() const { return fBytes; }
   Float_t        GetUsedCPU() const { return fUsedCPU; }
   TMacro        *GetLogFile() const { return fLogFile; }
   TMacro        *GetSelecHdr() const { return fSelecHdr; }
   TMacro        *GetSelecImp() const { return fSelecImp; }
   const char    *GetLibList() const { return fLibList; }
   const char    *GetParList() const { return fParList; }
   TList         *GetOutputList() { return fOutputList; }
   const char    *GetResultFile() const { return fResultFile; }
   Float_t        GetPrepTime() const { return fPrepTime; }
   Float_t        GetInitTime() const { return fInitTime; }
   Float_t        GetProcTime() const { return fProcTime; }
   Float_t        GetMergeTime() const { return fMergeTime; }
   Float_t        GetRecvTime() const { return fRecvTime; }
   Float_t        GetTermTime() const { return fTermTime; }
   Int_t          GetNumWrks() const { return fNumWrks; }
   Int_t          GetNumMergers() const { return fNumMergers; }

   Bool_t         IsArchived() const { return fArchived; }
   virtual Bool_t IsDone() const { return (fStatus > kRunning); }
   Bool_t         IsDraw() const { return fDraw; }
   Bool_t         IsFinalized() const { return fFinalized; }

   Bool_t         Matches(const char *ref);

   void Print(Option_t *opt = "") const;

   ClassDef(TQueryResult,5)  //Class describing a query
};

inline Bool_t operator!=(const TQueryResult &qr1,  const TQueryResult &qr2)
   { return !(qr1 == qr2); }

#endif
