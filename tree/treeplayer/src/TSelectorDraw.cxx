// @(#)root/treeplayer:$Id$
// Author: Rene Brun   08/01/2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TSelectorDraw
A specialized TSelector for TTree::Draw.
*/

#include "TSelectorDraw.h"
#include "TROOT.h"
#include "TH2.h"
#include "TH3.h"
#include "TView.h"
#include "TGraph.h"
#include "TPolyMarker3D.h"
#include "TDirectory.h"
#include "TVirtualPad.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTreeFormulaManager.h"
#include "TEnv.h"
#include "TTree.h"
#include "TCut.h"
#include "TEntryList.h"
#include "TEventList.h"
#include "TEntryListArray.h"
#include "THLimitsFinder.h"
#include "TStyle.h"
#include "TClass.h"
#include "TColor.h"
#include "strlcpy.h"

ClassImp(TSelectorDraw);

const Int_t kCustomHistogram = BIT(17);

////////////////////////////////////////////////////////////////////////////////
/// Default selector constructor.

TSelectorDraw::TSelectorDraw()
{
   fTree           = 0;
   fW              = 0;
   fValSize        = 4;
   fVal            = new Double_t*[fValSize];
   fVmin           = new Double_t[fValSize];
   fVmax           = new Double_t[fValSize];
   fNbins          = new Int_t[fValSize];
   fVarMultiple    = new Bool_t[fValSize];
   fVar            = new TTreeFormula*[fValSize];
   for (Int_t i = 0; i < fValSize; ++i) {
      fVal[i] = 0;
      fVar[i] = 0;
   }
   fManager        = 0;
   fMultiplicity   = 0;
   fSelect         = 0;
   fSelectedRows   = 0;
   fDraw           = 0;
   fObject         = 0;
   fOldHistogram   = 0;
   fObjEval        = kFALSE;
   fSelectMultiple = kFALSE;
   fCleanElist     = kFALSE;
   fTreeElist      = 0;
   fAction         = 0;
   fNfill          = 0;
   fDimension      = 0;
   fOldEstimate    = 0;
   fForceRead      = 0;
   fWeight         = 1;
   fCurrentSubEntry = -1;
   fTreeElistArray  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Selector destructor.

TSelectorDraw::~TSelectorDraw()
{
   ClearFormula();
   delete [] fVar;
   if (fVal) {
      for (Int_t i = 0; i < fValSize; ++i)
         delete [] fVal[i];
      delete [] fVal;
   }
   if (fVmin) delete [] fVmin;
   if (fVmax) delete [] fVmax;
   if (fNbins) delete [] fNbins;
   if (fVarMultiple) delete [] fVarMultiple;
   if (fW)     delete [] fW;
}

////////////////////////////////////////////////////////////////////////////////
/// Called every time a loop on the tree(s) starts.

void TSelectorDraw::Begin(TTree *tree)
{
   SetStatus(0);
   ResetAbort();
   ResetBit(kCustomHistogram);
   fSelectedRows   = 0;
   fTree = tree;
   fDimension = 0;
   fAction = 0;

   TObject *obj = fInput->FindObject("varexp");
   const char *varexp0   = obj ? obj->GetTitle() : "";
   obj = fInput->FindObject("selection");
   const char *selection = obj ? obj->GetTitle() : "";
   const char *option    = GetOption();

   TString  opt, abrt;
   char *hdefault = (char *)"htemp";
   char *varexp = nullptr;
   Int_t i, j, hkeep;
   opt = option;
   opt.ToLower();
   fOldHistogram = 0;
   TEntryList *enlist = 0;
   TEventList *evlist = 0;
   TString htitle;
   Bool_t profile = kFALSE;
   Bool_t optSame = kFALSE;
   Bool_t optEnlist = kFALSE;
   Bool_t optEnlistArray = kFALSE;
   Bool_t optpara = kFALSE;
   Bool_t optcandle = kFALSE;
   Bool_t opt5d = kFALSE;
   if (opt.Contains("same")) {
      optSame = kTRUE;
      opt.ReplaceAll("same", "");
   }
   if (opt.Contains("entrylist")) {
      optEnlist = kTRUE;
      if (opt.Contains("entrylistarray")) {
         optEnlistArray = kTRUE;
         opt.ReplaceAll("entrylistarray", "");
      } else {
         opt.ReplaceAll("entrylist", "");
      }
   }
   if (opt.Contains("para")) {
      optpara = kTRUE;
      opt.ReplaceAll("para", "");
   }
   if (opt.Contains("candle")) {
      optcandle = kTRUE;
      opt.ReplaceAll("candle", "");
   }
   if (opt.Contains("gl5d")) {
      opt5d = kTRUE;
      opt.ReplaceAll("gl5d", "");
   }
   TCut realSelection(selection);
   //input list - only TEntryList
   TEntryList *inElist = fTree->GetEntryList();
   evlist = fTree->GetEventList();
   if (evlist && inElist) {
      //this is needed because the input entry list was created
      //by the fTree from the input TEventList and is owned by the fTree.
      //Calling GetEntryList function changes ownership and here
      //we want fTree to still delete this entry list

      inElist->SetBit(kCanDelete, kTRUE);
   }
   fCleanElist = kFALSE;
   fTreeElist = inElist;

   fTreeElistArray = inElist ? dynamic_cast<TEntryListArray*>(fTreeElist) : 0;


   if (inElist && inElist->GetReapplyCut()) {
      realSelection *= inElist->GetTitle();
   }

   // what each variable should contain:
   //   varexp0   - original expression eg "a:b>>htest"
   //   hname     - name of new or old histogram
   //   hkeep     - flag if to keep new histogram
   //   hnameplus - flag if to add to current histo
   //   i         - length of variable expression stipped of everything after ">>"
   //   varexp    - variable expression stipped of everything after ">>"
   //   fOldHistogram     - pointer to hist hname
   //   elist     - pointer to selection list of hname

   Bool_t canExtend = kTRUE;
   if (optSame) canExtend = kFALSE;

   Int_t nbinsx = 0, nbinsy = 0, nbinsz = 0;
   Double_t xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;

   fObject  = 0;
   char *hname = nullptr;
   TString hnamealloc;
   i = 0;
   if (varexp0 && strlen(varexp0)) {
      for (UInt_t k = strlen(varexp0) - 1; k > 0; k--) {
         if (varexp0[k] == '>' && varexp0[k-1] == '>') {
            i = (int)(&(varexp0[k-1]) - varexp0);    //  length of varexp0 before ">>"
            hnamealloc = &(varexp0[k+1]);   // use TString to copy value of the
            hname = (char *) hnamealloc.Data();
            break;
         }
      }
   }
   //   char *hname = (char*)strstr(varexp0,">>");
   if (hname) {
      hkeep  = 1;
      varexp = new char[i+1];
      varexp[0] = 0; //necessary if i=0
      Bool_t hnameplus = kFALSE;
      while (*hname == ' ') hname++;
      if (*hname == '+') {
         hnameplus = kTRUE;
         hname++;
         while (*hname == ' ') hname++; //skip ' '
      }
      j = strlen(hname) - 1;   // skip ' '  at the end
      while (j) {
         if (hname[j] != ' ') break;
         hname[j] = 0;
         j--;
      }

      if (i) {
         strlcpy(varexp,varexp0,i+1);

         Int_t mustdelete = 0;
         SetBit(kCustomHistogram);

         // parse things that follow the name of the histo between '(' and ')'.
         // At this point hname contains the name of the specified histogram.
         //   Now the syntax is extended to handle an hname of the following format
         //   hname(nBIN [[,[xlow]][,xhigh]],...)
         //   so enclosed in brackets is the binning information, xlow, xhigh, and
         //   the same for the other dimensions

         char *pstart;    // pointer to '('
         char *pend;      // pointer to ')'
         char *cdummy;    // dummy pointer
         int ncomma;      // number of commas between '(' and ')', later number of arguments
         int ncols;       // number of columns in varexpr
         Double_t value;  // parsed value (by sscanf)

         const Int_t maxvalues = 9;

         pstart = strchr(hname, '(');
         pend =  strchr(hname, ')');
         if (pstart != 0) {   // found the bracket

            mustdelete = 1;

            // check that there is only one open and close bracket
            if (pstart == strrchr(hname, '(')  &&  pend == strrchr(hname, ')')) {

               // count number of ',' between '(' and ')'
               ncomma = 0;
               cdummy = pstart;
               cdummy = strchr(&cdummy[1], ',');
               while (cdummy != 0) {
                  cdummy = strchr(&cdummy[1], ',');
                  ncomma++;
               }

               if (ncomma + 1 > maxvalues) {
                  Error("DrawSelect", "ncomma+1>maxvalues, ncomma=%d, maxvalues=%d", ncomma, maxvalues);
                  ncomma = maxvalues - 1;
               }

               ncomma++; // number of arguments
               cdummy = pstart;

               //   number of columns
               ncols  = 1;
               for (j = 0; j < i; j++) {
                  if (varexp[j] == ':'
                      && !((j > 0 && varexp[j-1] == ':') || varexp[j+1] == ':')
                     ) {
                     ncols++;
                  }
               }
               if (ncols > 3) {   // max 3 columns
                  Error("DrawSelect", "ncols > 3, ncols=%d", ncols);
                  ncols = 0;
               }

               // check dimensions before and after ">>"
               if (ncols * 3 < ncomma) {
                  Error("DrawSelect", "ncols*3 < ncomma ncols=%d, ncomma=%d", ncols, ncomma);
                  ncomma = ncols * 3;
               }

               // scan the values one after the other
               for (j = 0; j < ncomma; j++) {
                  cdummy++;  // skip '(' or ','
                  if (sscanf(cdummy, " %lf ", &value) == 1) {
                     cdummy = strchr(&cdummy[1], ',');

                     switch (j) {  // do certain settings depending on position of argument
                        case 0:  // binning x-axis
                           nbinsx = (Int_t)value;
                           if (ncols < 2) {
                              gEnv->SetValue("Hist.Binning.1D.x", nbinsx);
                           } else if (ncols < 3) {
                              gEnv->SetValue("Hist.Binning.2D.x", nbinsx);
                              gEnv->SetValue("Hist.Binning.2D.Prof", nbinsx);
                           } else {
                              gEnv->SetValue("Hist.Binning.3D.x", nbinsx);
                              gEnv->SetValue("Hist.Binning.3D.Profx", nbinsx);
                           }

                           break;
                        case 1:  // lower limit x-axis
                           xmin = value;
                           break;
                        case 2:  // upper limit x-axis
                           xmax = value;
                           break;
                        case 3:  // binning y-axis
                           nbinsy = (Int_t)value;
                           if (ncols < 3) gEnv->SetValue("Hist.Binning.2D.y", nbinsy);
                           else {
                              gEnv->SetValue("Hist.Binning.3D.y", nbinsy);
                              gEnv->SetValue("Hist.Binning.3D.Profy", nbinsy);
                           }
                           break;
                        case 4:  // lower limit y-axis
                           ymin = value;
                           break;
                        case 5:  // upper limit y-axis
                           ymax = value;
                           break;
                        case 6:  // binning z-axis
                           nbinsz = (Int_t)value;
                           gEnv->SetValue("Hist.Binning.3D.z", nbinsz);
                           break;
                        case 7:  // lower limit z-axis
                           zmin = value;
                           break;
                        case 8:  // upper limit z-axis
                           zmax = value;
                           break;
                        default:
                           Error("DrawSelect", "j>8");
                           break;
                     }
                  }  // if sscanf == 1
               } // for j=0;j<ncomma;j++
            } else {
               Error("Begin", "Two open or close brackets found, hname=%s", hname);
            }

            // fix up hname
            pstart[0] = '\0'; // removes things after (and including) '('
         } // if '(' is found

         j = strlen(hname) - 1; // skip ' '  at the end
         while (j > 0) {
            if (hname[j] != ' ') break; // skip ' '  at the end
            hname[j] = 0;
            j--;
         }

         TObject *oldObject = gDirectory->Get(hname);  // if hname contains '(...)' the return values is NULL, which is what we want
         fOldHistogram = oldObject ? dynamic_cast<TH1*>(oldObject) : 0;

         if (!fOldHistogram && oldObject && !oldObject->InheritsFrom(TH1::Class())) {
            abrt.Form("An object of type '%s' has the same name as the requested histo (%s)", oldObject->IsA()->GetName(), hname);
            Abort(abrt);
            delete[] varexp;
            return;
         }
         if (fOldHistogram && !hnameplus) fOldHistogram->Reset();  // reset unless adding is wanted

         if (mustdelete) {
            if (gDebug) {
               Warning("Begin", "Deleting old histogram, since (possibly new) limits and binnings have been given");
            }
            delete fOldHistogram; fOldHistogram=0;
         }

      } else {
         // make selection list (i.e. varexp0 starts with ">>")
         TObject *oldObject = gDirectory->Get(hname);
         if (optEnlist) {
            //write into a TEntryList
            enlist = oldObject ? dynamic_cast<TEntryList*>(oldObject) : 0;

            if (!enlist && oldObject) {
               abrt.Form("An object of type '%s' has the same name as the requested event list (%s)",
                         oldObject->IsA()->GetName(), hname);
               Abort(abrt);
               delete[] varexp;
               return;
            }
            if (!enlist) {
               if (optEnlistArray) {
                  enlist = new TEntryListArray(hname, realSelection.GetTitle());
               } else {
                  enlist = new TEntryList(hname, realSelection.GetTitle());
               }
            }
            if (enlist) {
               if (!hnameplus) {
                  if (enlist == inElist) {
                     // We have been asked to reset the input list!!
                     // Let's set it aside for now ...
                     if (optEnlistArray) {
                        inElist = new TEntryListArray(*enlist);
                     } else {
                        inElist = new TEntryList(*enlist);
                     }
                     fCleanElist = kTRUE;
                     fTree->SetEntryList(inElist);
                  }
                  enlist->Reset();
                  enlist->SetTitle(realSelection.GetTitle());
               } else {
                  TCut old = enlist->GetTitle();
                  TCut upd = old || realSelection.GetTitle();
                  enlist->SetTitle(upd.GetTitle());
               }
            }
         } else {
            //write into a TEventList
            evlist = oldObject ? dynamic_cast<TEventList*>(oldObject) : 0;

            if (!evlist && oldObject) {
               abrt.Form("An object of type '%s' has the same name as the requested event list (%s)",
                         oldObject->IsA()->GetName(), hname);
               Abort(abrt);
               delete[] varexp;
               return;
            }
            if (!evlist) {
               evlist = new TEventList(hname, realSelection.GetTitle(), 1000, 0);
            }
            if (evlist) {
               if (!hnameplus) {
                  if (evlist == fTree->GetEventList()) {
                     // We have been asked to reset the input list!!
                     // Let's set it aside for now ...
                     Abort("Input and output lists are the same!");
                     delete[] varexp;
                     return;
                  }
                  evlist->Reset();
                  evlist->SetTitle(realSelection.GetTitle());
               } else {
                  TCut old = evlist->GetTitle();
                  TCut upd = old || realSelection.GetTitle();
                  evlist->SetTitle(upd.GetTitle());
               }
            }
         }

      }  // if (i)
   } else { // if (hname)
      hname  = hdefault;
      hkeep  = 0;
      const size_t varexpLen = strlen(varexp0) + 1;
      varexp = new char[varexpLen];
      strlcpy(varexp, varexp0, varexpLen);
      if (gDirectory) {
         fOldHistogram = (TH1*)gDirectory->Get(hname);
         if (fOldHistogram) { fOldHistogram->Delete(); fOldHistogram = 0;}
      }
   }

   // Decode varexp and selection
   if (!CompileVariables(varexp, realSelection.GetTitle())) {
      abrt.Form("Variable compilation failed: {%s,%s}", varexp, realSelection.GetTitle());
      Abort(abrt);
      delete[] varexp;
      return;
   }
   if (fDimension > 4 && !(optpara || optcandle || opt5d || opt.Contains("goff"))) {
      Abort("Too many variables. Use the option \"para\", \"gl5d\" or \"candle\" to display more than 4 variables.");
      delete[] varexp;
      return;
   }
   if (fDimension < 2 && (optpara || optcandle)) {
      Abort("The options \"para\" and \"candle\" require at least 2 variables.");
      delete[] varexp;
      return;
   }

   // In case fOldHistogram exists, check dimensionality
   Int_t nsel = strlen(selection);
   if (nsel > 1) {
      htitle.Form("%s {%s}", varexp, selection);
   } else {
      htitle = varexp;
   }
   if (fOldHistogram) {
      Int_t olddim = fOldHistogram->GetDimension();
      Int_t mustdelete = 0;
      if (fOldHistogram->InheritsFrom(TProfile::Class())) {
         profile = kTRUE;
         olddim = 2;
      }
      if (fOldHistogram->InheritsFrom(TProfile2D::Class())) {
         profile = kTRUE;
         olddim = 3;
      }
      if (opt.Contains("prof") && fDimension > 1) {
         // ignore "prof" for 1D.
         if (!profile || olddim != fDimension) mustdelete = 1;
      } else if (opt.Contains("col") && fDimension>2) {
         if (olddim+1 != fDimension) mustdelete = 1;
      } else {
         if (olddim != fDimension) mustdelete = 1;
      }
      if (mustdelete) {
         Warning("Begin", "Deleting old histogram with different dimensions");
         delete fOldHistogram;
         fOldHistogram = 0;
      }
   }

   // Create a default canvas if none exists
   fDraw = 0;
   if (!gPad && !opt.Contains("goff") && fDimension > 0) {
      gROOT->MakeDefCanvas();
      if (!gPad) {
         Abort("Creation of default canvas failed");
         delete[] varexp;
         return;
      }
   }

   // 1-D distribution
   TH1 *hist;
   if (fDimension == 1) {
      fAction = 1;
      if (!fOldHistogram) {
         fNbins[0] = gEnv->GetValue("Hist.Binning.1D.x", 100);
         if (gPad && optSame) {
            TListIter np(gPad->GetListOfPrimitives());
            TObject *op;
            TH1 *oldhtemp = 0;
            while ((op = np()) && !oldhtemp) {
               if (op->InheritsFrom(TH1::Class())) oldhtemp = (TH1 *)op;
            }
            if (oldhtemp) {
               fNbins[0] = oldhtemp->GetXaxis()->GetNbins();
               fVmin[0]  = oldhtemp->GetXaxis()->GetXmin();
               fVmax[0]  = oldhtemp->GetXaxis()->GetXmax();
            } else {
               fVmin[0]  = gPad->GetUxmin();
               fVmax[0]  = gPad->GetUxmax();
            }
         } else {
            fAction   = -1;
            fVmin[0] = xmin;
            fVmax[0] = xmax;
            if (xmin < xmax) canExtend = kFALSE;
         }
      }
      if (fOldHistogram) {
         hist    = fOldHistogram;
         fNbins[0] = hist->GetXaxis()->GetNbins();
      } else {
         TString precision = gEnv->GetValue("Hist.Precision.1D", "float");
         if (precision.Contains("float")) {
            hist = new TH1F(hname, htitle.Data(), fNbins[0], fVmin[0], fVmax[0]);
         } else {
            hist = new TH1D(hname, htitle.Data(), fNbins[0], fVmin[0], fVmax[0]);
         }
         hist->SetLineColor(fTree->GetLineColor());
         hist->SetLineWidth(fTree->GetLineWidth());
         hist->SetLineStyle(fTree->GetLineStyle());
         hist->SetFillColor(fTree->GetFillColor());
         hist->SetFillStyle(fTree->GetFillStyle());
         hist->SetMarkerStyle(fTree->GetMarkerStyle());
         hist->SetMarkerColor(fTree->GetMarkerColor());
         hist->SetMarkerSize(fTree->GetMarkerSize());
         if (canExtend) hist->SetCanExtend(TH1::kAllAxes);
         if (!hkeep) {
            hist->GetXaxis()->SetTitle(fVar[0]->GetTitle());
            hist->SetBit(kCanDelete);
            if (!opt.Contains("goff")) hist->SetDirectory(0);
         }
         if (opt.Length() && opt.Contains("e")) hist->Sumw2();
      }
      fVar[0]->SetAxis(hist->GetXaxis());
      fObject = hist;

      // 2-D distribution
   } else if (fDimension == 2 && !(optpara || optcandle)) {
      fAction = 2;
      if (!fOldHistogram || !optSame) {
         fNbins[0] = gEnv->GetValue("Hist.Binning.2D.y", 40);
         fNbins[1] = gEnv->GetValue("Hist.Binning.2D.x", 40);
         if (opt.Contains("prof")) fNbins[1] = gEnv->GetValue("Hist.Binning.2D.Prof", 100);
         if (optSame) {
            TH1 *oldhtemp = (TH1*)gPad->FindObject(hdefault);
            if (oldhtemp) {
               fNbins[1] = oldhtemp->GetXaxis()->GetNbins();
               fVmin[1]  = oldhtemp->GetXaxis()->GetXmin();
               fVmax[1]  = oldhtemp->GetXaxis()->GetXmax();
               fNbins[0] = oldhtemp->GetYaxis()->GetNbins();
               fVmin[0]  = oldhtemp->GetYaxis()->GetXmin();
               fVmax[0]  = oldhtemp->GetYaxis()->GetXmax();
            } else {
               fNbins[1] = gEnv->GetValue("Hist.Binning.2D.x", 40);
               fVmin[1]  = gPad->GetUxmin();
               fVmax[1]  = gPad->GetUxmax();
               fNbins[0] = gEnv->GetValue("Hist.Binning.2D.y", 40);
               fVmin[0]  = gPad->GetUymin();
               fVmax[0]  = gPad->GetUymax();
            }
         } else {
            if (!fOldHistogram) fAction = -2;
            fVmin[1] = xmin;
            fVmax[1] = xmax;
            fVmin[0] = ymin;
            fVmax[0] = ymax;
            if (xmin < xmax && ymin < ymax) canExtend = kFALSE;
         }
      }
      if (profile || opt.Contains("prof")) {
         TProfile *hp;
         if (fOldHistogram) {
            fAction = 4;
            hp = (TProfile*)fOldHistogram;
         } else {
            if (fAction < 0) {
               fAction = -4;
               fVmin[1] = xmin;
               fVmax[1] = xmax;
               if (xmin < xmax) canExtend = kFALSE;
            }
            if (fAction == 2) {
               //we come here when option = "same prof"
               fAction = -4;
               TH1 *oldhtemp = (TH1*)gPad->FindObject(hdefault);
               if (oldhtemp) {
                  fNbins[1] = oldhtemp->GetXaxis()->GetNbins();
                  fVmin[1]  = oldhtemp->GetXaxis()->GetXmin();
                  fVmax[1]  = oldhtemp->GetXaxis()->GetXmax();
               }
            }
            if (opt.Contains("profs")) {
               hp = new TProfile(hname, htitle.Data(), fNbins[1], fVmin[1], fVmax[1], "s");
            } else if (opt.Contains("profi")) {
               hp = new TProfile(hname, htitle.Data(), fNbins[1], fVmin[1], fVmax[1], "i");
            } else if (opt.Contains("profg")) {
               hp = new TProfile(hname, htitle.Data(), fNbins[1], fVmin[1], fVmax[1], "g");
            } else {
               hp = new TProfile(hname, htitle.Data(), fNbins[1], fVmin[1], fVmax[1], "");
            }
            if (!hkeep) {
               hp->SetBit(kCanDelete);
               if (!opt.Contains("goff")) hp->SetDirectory(0);
            }
            hp->SetLineColor(fTree->GetLineColor());
            hp->SetLineWidth(fTree->GetLineWidth());
            hp->SetLineStyle(fTree->GetLineStyle());
            hp->SetFillColor(fTree->GetFillColor());
            hp->SetFillStyle(fTree->GetFillStyle());
            hp->SetMarkerStyle(fTree->GetMarkerStyle());
            hp->SetMarkerColor(fTree->GetMarkerColor());
            hp->SetMarkerSize(fTree->GetMarkerSize());
            if (canExtend) hp->SetCanExtend(TH1::kAllAxes);
         }
         fVar[1]->SetAxis(hp->GetXaxis());
         fObject = hp;

      } else {
         TH2 *h2;
         if (fOldHistogram) {
            h2 = (TH2F*)fOldHistogram;
         } else {
            TString precision = gEnv->GetValue("Hist.Precision.2D", "float");
            if (precision.Contains("float")) {
               h2 = new TH2F(hname, htitle.Data(), fNbins[1], fVmin[1], fVmax[1], fNbins[0], fVmin[0], fVmax[0]);
            } else {
               h2 = new TH2D(hname, htitle.Data(), fNbins[1], fVmin[1], fVmax[1], fNbins[0], fVmin[0], fVmax[0]);
            }
            h2->SetLineColor(fTree->GetLineColor());
            h2->SetLineWidth(fTree->GetLineWidth());
            h2->SetLineStyle(fTree->GetLineStyle());
            h2->SetFillColor(fTree->GetFillColor());
            h2->SetFillStyle(fTree->GetFillStyle());
            h2->SetMarkerStyle(fTree->GetMarkerStyle());
            h2->SetMarkerColor(fTree->GetMarkerColor());
            h2->SetMarkerSize(fTree->GetMarkerSize());
            if (canExtend) h2->SetCanExtend(TH1::kAllAxes);
            if (!hkeep) {
               h2->GetXaxis()->SetTitle(fVar[1]->GetTitle());
               h2->GetYaxis()->SetTitle(fVar[0]->GetTitle());
               h2->SetBit(TH1::kNoStats);
               h2->SetBit(kCanDelete);
               if (!opt.Contains("goff")) h2->SetDirectory(0);
            }
         }
         fVar[0]->SetAxis(h2->GetYaxis());
         fVar[1]->SetAxis(h2->GetXaxis());
         Bool_t graph = kFALSE;
         Int_t l = opt.Length();
         if (l == 0 || optSame) graph = kTRUE;
         if (opt.Contains("p")     || opt.Contains("*")    || opt.Contains("l"))    graph = kTRUE;
         if (opt.Contains("surf")  || opt.Contains("lego") || opt.Contains("cont")) graph = kFALSE;
         if (opt.Contains("col")   || opt.Contains("hist") || opt.Contains("scat")) graph = kFALSE;
         if (opt.Contains("box"))                                                   graph = kFALSE;
         fObject = h2;
         if (graph) {
            fAction = 12;
            if (!fOldHistogram && !optSame) fAction = -12;
         }
      }

      // 3-D distribution
   } else if ((fDimension == 3 || fDimension == 4) && !(optpara || optcandle)) {
      fAction = 3;
      if (fDimension == 4) fAction = 40;
      if (!fOldHistogram || !optSame) {
         fNbins[0] = gEnv->GetValue("Hist.Binning.3D.z", 20);
         fNbins[1] = gEnv->GetValue("Hist.Binning.3D.y", 20);
         fNbins[2] = gEnv->GetValue("Hist.Binning.3D.x", 20);
         if (fDimension == 3 && opt.Contains("prof")) {
            fNbins[1] = gEnv->GetValue("Hist.Binning.3D.Profy", 20);
            fNbins[2] = gEnv->GetValue("Hist.Binning.3D.Profx", 20);
         } else if (fDimension == 3 && opt.Contains("col")) {
            fNbins[0] = gEnv->GetValue("Hist.Binning.2D.y", 40);
            fNbins[1] = gEnv->GetValue("Hist.Binning.2D.x", 40);
         }
         if (optSame) {
            TH1 *oldhtemp = (TH1*)gPad->FindObject(hdefault);
            if (oldhtemp) {
               fNbins[2] = oldhtemp->GetXaxis()->GetNbins();
               fVmin[2]  = oldhtemp->GetXaxis()->GetXmin();
               fVmax[2]  = oldhtemp->GetXaxis()->GetXmax();
               fNbins[1] = oldhtemp->GetYaxis()->GetNbins();
               fVmin[1]  = oldhtemp->GetYaxis()->GetXmin();
               fVmax[1]  = oldhtemp->GetYaxis()->GetXmax();
               fNbins[0] = oldhtemp->GetZaxis()->GetNbins();
               fVmin[0]  = oldhtemp->GetZaxis()->GetXmin();
               fVmax[0]  = oldhtemp->GetZaxis()->GetXmax();
            } else {
               TView *view = gPad->GetView();
               if (!view) {
                  Error("Begin", "You cannot use option same when no 3D view exists");
                  fVmin[0] = fVmin[1] = fVmin[2] = -1;
                  fVmax[0] = fVmax[1] = fVmax[2] = 1;
                  view = TView::CreateView(1, fVmin, fVmax);
               }
               Double_t *rmin = view->GetRmin();
               Double_t *rmax = view->GetRmax();
               fNbins[2] = gEnv->GetValue("Hist.Binning.3D.z", 20);
               fVmin[2]  = rmin[0];
               fVmax[2]  = rmax[0];
               fNbins[1] = gEnv->GetValue("Hist.Binning.3D.y", 20);
               fVmin[1]  = rmin[1];
               fVmax[1]  = rmax[1];
               fNbins[0] = gEnv->GetValue("Hist.Binning.3D.x", 20);
               fVmin[0]  = rmin[2];
               fVmax[0]  = rmax[2];
            }
         } else {
            if (!fOldHistogram && fDimension == 3) fAction = -3;
            fVmin[2] = xmin;
            fVmax[2] = xmax;
            fVmin[1] = ymin;
            fVmax[1] = ymax;
            fVmin[0] = zmin;
            fVmax[0] = zmax;
            if (xmin < xmax && ymin < ymax && zmin < zmax) canExtend = kFALSE;
         }
      }
      if ((fDimension == 3) && (profile || opt.Contains("prof"))) {
         TProfile2D *hp;
         if (fOldHistogram) {
            fAction = 23;
            hp = (TProfile2D*)fOldHistogram;
         } else {
            if (fAction < 0) {
               fAction = -23;
               fVmin[2] = xmin;
               fVmax[2] = xmax;
               fVmin[1] = ymin;
               fVmax[1] = ymax;
               if (xmin < xmax && ymin < ymax) canExtend = kFALSE;
            }
            if (opt.Contains("profs")) {
               hp = new TProfile2D(hname, htitle.Data(), fNbins[2], fVmin[2], fVmax[2], fNbins[1], fVmin[1], fVmax[1], "s");
            } else if (opt.Contains("profi")) {
               hp = new TProfile2D(hname, htitle.Data(), fNbins[2], fVmin[2], fVmax[2], fNbins[1], fVmin[1], fVmax[1], "i");
            } else if (opt.Contains("profg")) {
               hp = new TProfile2D(hname, htitle.Data(), fNbins[2], fVmin[2], fVmax[2], fNbins[1], fVmin[1], fVmax[1], "g");
            } else {
               hp = new TProfile2D(hname, htitle.Data(), fNbins[2], fVmin[2], fVmax[2], fNbins[1], fVmin[1], fVmax[1], "");
            }
            if (!hkeep) {
               hp->SetBit(kCanDelete);
               if (!opt.Contains("goff")) hp->SetDirectory(0);
            }
            hp->SetLineColor(fTree->GetLineColor());
            hp->SetLineWidth(fTree->GetLineWidth());
            hp->SetLineStyle(fTree->GetLineStyle());
            hp->SetFillColor(fTree->GetFillColor());
            hp->SetFillStyle(fTree->GetFillStyle());
            hp->SetMarkerStyle(fTree->GetMarkerStyle());
            hp->SetMarkerColor(fTree->GetMarkerColor());
            hp->SetMarkerSize(fTree->GetMarkerSize());
            if (canExtend) hp->SetCanExtend(TH1::kAllAxes);
         }
         fVar[1]->SetAxis(hp->GetYaxis());
         fVar[2]->SetAxis(hp->GetXaxis());
         fObject = hp;
      } else if (fDimension == 3 && opt.Contains("col")) {
         TH2F *h2;
         if (fOldHistogram) {
            h2 = (TH2F*)fOldHistogram;
         } else {
            h2 = new TH2F(hname, htitle.Data(), fNbins[1], fVmin[1], fVmax[1], fNbins[0], fVmin[0], fVmax[0]);
            h2->SetLineColor(fTree->GetLineColor());
            h2->SetLineWidth(fTree->GetLineWidth());
            h2->SetLineStyle(fTree->GetLineStyle());
            h2->SetFillColor(fTree->GetFillColor());
            h2->SetFillStyle(fTree->GetFillStyle());
            h2->SetMarkerStyle(fTree->GetMarkerStyle());
            h2->SetMarkerColor(fTree->GetMarkerColor());
            h2->SetMarkerSize(fTree->GetMarkerSize());
            if (canExtend) h2->SetCanExtend(TH1::kAllAxes);
            if (!hkeep) {
               h2->GetXaxis()->SetTitle(fVar[1]->GetTitle());
               h2->GetYaxis()->SetTitle(fVar[0]->GetTitle());
               h2->GetZaxis()->SetTitle(fVar[2]->GetTitle());
               h2->SetBit(TH1::kNoStats);
               h2->SetBit(kCanDelete);
               if (!opt.Contains("goff")) h2->SetDirectory(0);
            }
         }
         fVar[0]->SetAxis(h2->GetYaxis());
         fVar[1]->SetAxis(h2->GetXaxis());
         fObject = h2;
         fAction = 33;
      } else {
         TH3 *h3;
         if (fOldHistogram) {
            h3 = (TH3F*)fOldHistogram;
         } else {
            TString precision = gEnv->GetValue("Hist.Precision.3D", "float");
            if (precision.Contains("float")) {
               h3 = new TH3F(hname, htitle.Data(), fNbins[2], fVmin[2], fVmax[2], fNbins[1], fVmin[1], fVmax[1], fNbins[0], fVmin[0], fVmax[0]);
            } else {
               h3 = new TH3D(hname, htitle.Data(), fNbins[2], fVmin[2], fVmax[2], fNbins[1], fVmin[1], fVmax[1], fNbins[0], fVmin[0], fVmax[0]);
            }
            h3->SetLineColor(fTree->GetLineColor());
            h3->SetLineWidth(fTree->GetLineWidth());
            h3->SetLineStyle(fTree->GetLineStyle());
            h3->SetFillColor(fTree->GetFillColor());
            h3->SetFillStyle(fTree->GetFillStyle());
            h3->SetMarkerStyle(fTree->GetMarkerStyle());
            h3->SetMarkerColor(fTree->GetMarkerColor());
            h3->SetMarkerSize(fTree->GetMarkerSize());
            if (canExtend) h3->SetCanExtend(TH1::kAllAxes);
            if (!hkeep) {
               //small correction for the title offsets in x,y to take into account the angles
               Double_t xoffset = h3->GetXaxis()->GetTitleOffset();
               Double_t yoffset = h3->GetYaxis()->GetTitleOffset();
               h3->GetXaxis()->SetTitleOffset(1.2 * xoffset);
               h3->GetYaxis()->SetTitleOffset(1.2 * yoffset);
               h3->GetXaxis()->SetTitle(fVar[2]->GetTitle());
               h3->GetYaxis()->SetTitle(fVar[1]->GetTitle());
               h3->GetZaxis()->SetTitle(fVar[0]->GetTitle());
               h3->SetBit(kCanDelete);
               h3->SetBit(TH1::kNoStats);
               if (!opt.Contains("goff")) h3->SetDirectory(0);
            }
         }
         fVar[0]->SetAxis(h3->GetZaxis());
         fVar[1]->SetAxis(h3->GetYaxis());
         fVar[2]->SetAxis(h3->GetXaxis());
         fObject = h3;
         Int_t noscat = strlen(option);
         if (optSame) noscat -= 4;
         if (!noscat && fDimension == 3) {
            fAction = 13;
            if (!fOldHistogram && !optSame) fAction = -13;
         }
      }
      // An Event List
   } else if (enlist) {
      fAction = 5;
      fOldEstimate = fTree->GetEstimate();
      fTree->SetEstimate(1);
      fObject = enlist;
   } else if (evlist) {
      fAction = 5;
      fOldEstimate = fTree->GetEstimate();
      fTree->SetEstimate(1);
      fObject = evlist;
   } else if (optcandle || optpara || opt5d) {
      if (optcandle)  fAction = 7;
      else if (opt5d) fAction = 8;
      else            fAction = 6;
   }
   if (varexp) delete[] varexp;
   for (i = 0; i < fValSize; ++i)
      fVarMultiple[i] = kFALSE;
   fSelectMultiple = kFALSE;
   for (i = 0; i < fDimension; ++i) {
      if (fVar[i] && fVar[i]->GetMultiplicity()) fVarMultiple[i] = kTRUE;
   }

   if (fSelect && fSelect->GetMultiplicity()) fSelectMultiple = kTRUE;

   fForceRead = fTree->TestBit(TTree::kForceRead);
   fWeight  = fTree->GetWeight();
   fNfill   = 0;

   for (i = 0; i < fDimension; ++i) {
      if (!fVal[i] && fVar[i]) {
         fVal[i] = new Double_t[(Int_t)fTree->GetEstimate()];
      }
   }

   if (!fW)             fW  = new Double_t[(Int_t)fTree->GetEstimate()];

   for (i = 0; i < fValSize; ++i) {
      fVmin[i] = DBL_MAX;
      fVmax[i] = -DBL_MAX;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete internal buffers.

void TSelectorDraw::ClearFormula()
{
   ResetBit(kWarn);
   for (Int_t i = 0; i < fValSize; ++i) {
      delete fVar[i];
      fVar[i] = 0;
   }
   delete fSelect; fSelect = 0;
   fManager = 0;
   fMultiplicity = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Compile input variables and selection expression.
///
/// varexp is an expression of the general form e1:e2:e3
/// where e1,etc is a formula referencing a combination of the columns
///
/// Example:
///
///     varexp = x  simplest case: draw a 1-Dim distribution of column named x
///            = sqrt(x)         : draw distribution of sqrt(x)
///            = x*y/z
///            = y:sqrt(x) 2-Dim distribution of y versus sqrt(x)
///
/// selection is an expression with a combination of the columns
///
/// Example:
///
///     selection = "x<y && sqrt(z)>3.2"
///
/// in a selection all the C++ operators are authorized
///
/// Return kFALSE if any of the variable is not compilable.

Bool_t TSelectorDraw::CompileVariables(const char *varexp, const char *selection)
{
   Int_t i, nch, ncols;

   // Compile selection expression if there is one
   fDimension = 0;
   ClearFormula();
   fMultiplicity = 0;
   fObjEval = kFALSE;

   if (strlen(selection)) {
      fSelect = new TTreeFormula("Selection", selection, fTree);
      fSelect->SetQuickLoad(kTRUE);
      if (!fSelect->GetNdim()) {
         delete fSelect;
         fSelect = 0;
         return kFALSE;
      }
   }

   // if varexp is empty, take first column by default
   nch = strlen(varexp);
   if (nch == 0) {
      fDimension = 0;
      if (fSelect) {
         fManager = fSelect->GetManager();
      }
      fTree->ResetBit(TTree::kForceRead);

      if (fManager) {
         fManager->Sync();

         if (fManager->GetMultiplicity() == -1) fTree->SetBit(TTree::kForceRead);
         if (fManager->GetMultiplicity() >= 1) fMultiplicity = fManager->GetMultiplicity();
      }

      return kTRUE;
   }

   // otherwise select only the specified columns
   std::vector<TString> varnames;
   ncols  = SplitNames(varexp, varnames);

   InitArrays(ncols);

   fManager = new TTreeFormulaManager();
   if (fSelect) fManager->Add(fSelect);
   fTree->ResetBit(TTree::kForceRead);
   for (i = 0; i < ncols; ++i) {
      fVar[i] = new TTreeFormula(TString::Format("Var%i", i + 1), varnames[i].Data(), fTree);
      fVar[i]->SetQuickLoad(kTRUE);
      if(!fVar[i]->GetNdim()) { ClearFormula(); return kFALSE; }
      fManager->Add(fVar[i]);
   }
   fManager->Sync();

   if (fManager->GetMultiplicity() == -1) fTree->SetBit(TTree::kForceRead);
   if (fManager->GetMultiplicity() >= 1) fMultiplicity = fManager->GetMultiplicity();

   fDimension    = ncols;

   if (ncols == 1) {
      TClass *cl = fVar[0]->EvalClass();
      if (cl) {
         fObjEval = kTRUE;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the last values corresponding to the i-th component
/// of the formula being processed (where the component are ':' separated).
/// The actual number of entries is:
///
///     GetSelectedRows() % tree->GetEstimate()
///
/// Note GetSelectedRows currently returns the actual number of values plotted
/// and thus if the formula contains arrays, this number might be greater than
/// the number of entries in the trees.
///
/// By default TTree::Draw creates the arrays obtained
/// with all GetVal and GetW with a length corresponding to the
/// parameter fEstimate. By default fEstimate=10000 and can be modified
/// via TTree::SetEstimate. A possible recipe is to do
///
///     tree->SetEstimate(tree->GetEntries());
///
/// You must call SetEstimate if the expected number of selected rows
/// is greater than 10000.
///
/// See TTree::Draw for additional details.

Double_t* TSelectorDraw::GetVal(Int_t i) const
{
   if (i < 0 || i >= fDimension)
      return 0;
   else
      return fVal[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Return the TTreeFormula corresponding to the i-th component
/// of the request formula (where the component are ':' separated).

TTreeFormula* TSelectorDraw::GetVar(Int_t i) const
{
   if (i < 0 || i >= fDimension)
      return 0;
   else
      return fVar[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Initialization of the primitive type arrays if the new size is bigger than the available space.

void TSelectorDraw::InitArrays(Int_t newsize)
{
   if (newsize > fValSize) {
      Int_t oldsize = fValSize;
      while (fValSize < newsize)
         fValSize *= 2;          // Double the available space until it matches the new size.
      if (fNbins) delete [] fNbins;
      if (fVmin) delete [] fVmin;
      if (fVmax) delete [] fVmax;
      if (fVarMultiple) delete [] fVarMultiple;

      fNbins = new Int_t[fValSize];
      fVmin = new Double_t[fValSize];
      fVmax = new Double_t[fValSize];
      fVarMultiple = new Bool_t[fValSize];

      for (Int_t i = 0; i < oldsize; ++i)
         delete [] fVal[i];
      delete [] fVal;
      delete [] fVar;
      fVal = new Double_t*[fValSize];
      fVar = new TTreeFormula*[fValSize];
      for (Int_t i = 0; i < fValSize; ++i) {
         fVal[i] = 0;
         fVar[i] = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Build Index array for names in varexp.
/// This will allocated a C style array of TString and Ints

UInt_t TSelectorDraw::SplitNames(const TString &varexp, std::vector<TString> &names)
{
   names.clear();

   Bool_t ternary = kFALSE;
   Int_t prev = 0;
   for (Int_t i = 0; i < varexp.Length(); i++) {
      if (varexp[i] == ':'
          && !((i > 0 && varexp[i-1] == ':') || varexp[i+1] == ':')
         ) {
         if (ternary) {
            ternary = kFALSE;
         } else {
            names.push_back(varexp(prev, i - prev));
            prev = i + 1;
         }
      }
      if (varexp[i] == '?') {
         ternary = kTRUE;
      }
   }
   names.push_back(varexp(prev, varexp.Length() - prev));
   return names.size();
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called at the first entry of a new tree in a chain.

Bool_t TSelectorDraw::Notify()
{
   if (fTree) fWeight  = fTree->GetWeight();
   if (fVar) {
      for (Int_t i = 0; i < fDimension; ++i) {
         if (fVar[i]) fVar[i]->UpdateFormulaLeaves();
      }
   }
   if (fSelect) fSelect->UpdateFormulaLeaves();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Called in the entry loop for all entries accepted by Select.

void TSelectorDraw::ProcessFill(Long64_t entry)
{
   if (fObjEval) {
      ProcessFillObject(entry);
      return;
   }

   if (fMultiplicity) {
      ProcessFillMultiple(entry);
      return;
   }

   // simple case with no multiplicity
   if (fForceRead && fManager->GetNdata() <= 0) return;

   if (fSelect) {
      fW[fNfill] = fWeight * fSelect->EvalInstance(0);
      if (!fW[fNfill]) return;
   } else fW[fNfill] = fWeight;
   if (fVal) {
      for (Int_t i = 0; i < fDimension; ++i) {
         if (fVar[i]) fVal[i][fNfill] = fVar[i]->EvalInstance(0);
      }
   }
   fNfill++;
   if (fNfill >= fTree->GetEstimate()) {
      TakeAction();
      fNfill = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called in the entry loop for all entries accepted by Select.
/// Complex case with multiplicity.

void TSelectorDraw::ProcessFillMultiple(Long64_t entry)
{
   // Grab the array size of the formulas for this entry
   Int_t ndata = fManager->GetNdata();

   // No data at all, let's move on to the next entry.
   if (!ndata) return;

   // If the entry list is a TEntryListArray, get the selected subentries for this entry
   TEntryList *subList = 0;
   if (fTreeElistArray) {
      subList = fTreeElistArray->GetSubListForEntry(entry, fTree->GetTree());
   }

   Int_t nfill0 = fNfill;

   // Calculate the first values
   if (fSelect) {
      // coverity[var_deref_model] fSelectMultiple==kTRUE => fSelect != 0
      fW[fNfill] = fWeight * fSelect->EvalInstance(0);
      if (!fW[fNfill] && !fSelectMultiple) return;
   } else fW[fNfill] = fWeight;

   // Always call EvalInstance(0) to insure the loading
   // of the branches.
   if (fW[fNfill] && (!subList || subList->Contains(0))) {
      if (fDimension == 0 && fSelectMultiple) fCurrentSubEntry = (Long64_t) 0; // to fill TEntryListArray
      for (Int_t i = 0; i < fDimension; ++i) {
         if (fVar[i]) fVal[i][fNfill] = fVar[i]->EvalInstance(0);
      }
      fNfill++;
      if (fNfill >= fTree->GetEstimate()) {
         TakeAction();
         fNfill = 0;
      }
   } else {
      for (Int_t i = 0; i < fDimension; ++i) {
         if (fVar[i]) fVar[i]->ResetLoading();
      }
   }
   Double_t ww = fW[nfill0];

   for (Int_t i = 1; i < ndata; i++) {
      if (subList && !subList->Contains(i)) continue;
      if (fSelectMultiple) {
         // coverity[var_deref_model] fSelectMultiple==kTRUE => fSelect != 0
         ww = fWeight * fSelect->EvalInstance(i);
         if (ww == 0) continue;
         if (fNfill == nfill0) {
            for (Int_t k = 0; k < fDimension; ++k) {
               if (!fVarMultiple[k]) fVal[k][fNfill] = fVar[k]->EvalInstance(0);
            }
         }
         if (fDimension == 0) fCurrentSubEntry = (Long64_t) i; // to fill TEntryListArray
      }
      for (Int_t k = 0; k < fDimension; ++k) {
         if (fVarMultiple[k]) fVal[k][fNfill] = fVar[k]->EvalInstance(i);
         else                 fVal[k][fNfill] = fVal[k][nfill0];
      }
      fW[fNfill] = ww;

      fNfill++;
      if (fNfill >= fTree->GetEstimate()) {
         TakeAction();
         fNfill = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called in the entry loop for all entries accepted by Select.
/// Case where the only variable returns an object (or pointer to).

void TSelectorDraw::ProcessFillObject(Long64_t /*entry*/)
{
   // Complex case with multiplicity.

   // Grab the array size of the formulas for this entry
   Int_t ndata = fManager->GetNdata();

   // No data at all, let's move on to the next entry.
   if (!ndata) return;

   Int_t nfill0 = fNfill;
   Double_t ww = 0;

   for (Int_t i = 0; i < ndata; i++) {
      if (i == 0) {
         if (fSelect) {
            fW[fNfill] = fWeight * fSelect->EvalInstance(0);
            if (!fW[fNfill] && !fSelectMultiple) return;
         } else fW[fNfill] = fWeight;
         ww = fW[nfill0];
      } else if (fSelectMultiple) {
         ww = fWeight * fSelect->EvalInstance(i);
         if (ww == 0) continue;
      }
      if (fDimension >= 1 && fVar[0]) {
         TClass *cl = fVar[0]->EvalClass();
         if (cl == TBits::Class()) {

            void *obj = fVar[0]->EvalObject(i);

            if (obj) {
               TBits *bits = (TBits*)obj;
               Int_t nbits = bits->GetNbits();

               Int_t nextbit = -1;
               while (1) {
                  nextbit = bits->FirstSetBit(nextbit + 1);
                  if (nextbit >= nbits) break;
                  fVal[0][fNfill] = nextbit;
                  fW[fNfill] =  ww;
                  fNfill++;
               }
            }

         } else {

            if (!TestBit(kWarn)) {
               Warning("ProcessFillObject",
                       "Not implemented for %s",
                       cl ? cl->GetName() : "unknown class");
               SetBit(kWarn);
            }

         }
      }
      if (fNfill >= fTree->GetEstimate()) {
         TakeAction();
         fNfill = 0;
      }
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Set number of entries to estimate variable limits.

void TSelectorDraw::SetEstimate(Long64_t)
{
   if (fVal) {
      for (Int_t i = 0; i < fValSize; ++i) {
         delete [] fVal[i];
         fVal[i] = 0;
      }
   }
   delete [] fW;   fW  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action for object obj fNfill times.

void TSelectorDraw::TakeAction()
{
   Int_t i;
   //__________________________1-D histogram_______________________
   if (fAction ==  1)((TH1*)fObject)->FillN(fNfill, fVal[0], fW);
   //__________________________2-D histogram_______________________
   else if (fAction ==  2) {
      TH2 *h2 = (TH2*)fObject;
      for (i = 0; i < fNfill; i++) h2->Fill(fVal[1][i], fVal[0][i], fW[i]);
   }
   //__________________________Profile histogram_______________________
   else if (fAction ==  4)((TProfile*)fObject)->FillN(fNfill, fVal[1], fVal[0], fW);
   //__________________________Event List______________________________
   else if (fAction ==  5) {
      if (fObject->InheritsFrom(TEntryListArray::Class())) {
         TEntryListArray *enlistarray = (TEntryListArray*)fObject;
         Long64_t enumb = fTree->GetTree()->GetReadEntry();
         enlistarray->Enter(enumb, 0, fCurrentSubEntry);
      } else if (fObject->InheritsFrom(TEntryList::Class())) {
         TEntryList *enlist = (TEntryList*)fObject;
         Long64_t enumb = fTree->GetTree()->GetReadEntry();
         enlist->Enter(enumb);
      } else {
         TEventList *evlist = (TEventList*)fObject;
         Long64_t enumb = fTree->GetChainOffset() + fTree->GetTree()->GetReadEntry();
         if (evlist->GetIndex(enumb) < 0) evlist->Enter(enumb);
      }
   }
   //__________________________2D scatter plot_______________________
   else if (fAction == 12) {
      TH2 *h2 = (TH2*)fObject;
      if (h2->CanExtendAllAxes() && h2->TestBit(kCanDelete)) {
         for (i = 0; i < fNfill; i++) {
            if (fVmin[0] > fVal[0][i]) fVmin[0] = fVal[0][i];
            if (fVmax[0] < fVal[0][i]) fVmax[0] = fVal[0][i];
            if (fVmin[1] > fVal[1][i]) fVmin[1] = fVal[1][i];
            if (fVmax[1] < fVal[1][i]) fVmax[1] = fVal[1][i];
         }
         THLimitsFinder::GetLimitsFinder()->FindGoodLimits(h2, fVmin[1], fVmax[1], fVmin[0], fVmax[0]);
      }
      TGraph *pm = new TGraph(fNfill, fVal[1], fVal[0]);
      pm->SetEditable(kFALSE);
      pm->SetBit(kCanDelete);
      pm->SetMarkerStyle(fTree->GetMarkerStyle());
      pm->SetMarkerColor(fTree->GetMarkerColor());
      pm->SetMarkerSize(fTree->GetMarkerSize());
      pm->SetLineColor(fTree->GetLineColor());
      pm->SetLineWidth(fTree->GetLineWidth());
      pm->SetLineStyle(fTree->GetLineStyle());
      pm->SetFillColor(fTree->GetFillColor());
      pm->SetFillStyle(fTree->GetFillStyle());

      if (!fDraw && !strstr(fOption.Data(), "goff")) {
         if (fOption.Length() == 0 || strcasecmp(fOption.Data(), "same") == 0)  pm->Draw("p");
         else                                                                pm->Draw(fOption.Data());
      }
      if (!h2->TestBit(kCanDelete)) {
         for (i = 0; i < fNfill; i++) h2->Fill(fVal[1][i], fVal[0][i], fW[i]);
      }
   }
   //__________________________3D scatter plot_______________________
   else if (fAction ==  3) {
      TH3 *h3 = (TH3*)fObject;
      if (!h3->TestBit(kCanDelete)) {
         for (i = 0; i < fNfill; i++) h3->Fill(fVal[2][i], fVal[1][i], fVal[0][i], fW[i]);
      }
   } else if (fAction == 13) {
      TPolyMarker3D *pm3d = new TPolyMarker3D(fNfill);
      pm3d->SetMarkerStyle(fTree->GetMarkerStyle());
      pm3d->SetMarkerColor(fTree->GetMarkerColor());
      pm3d->SetMarkerSize(fTree->GetMarkerSize());
      for (i = 0; i < fNfill; i++) {
         pm3d->SetPoint(i, fVal[2][i], fVal[1][i], fVal[0][i]);
      }
      pm3d->Draw();
      TH3 *h3 = (TH3*)fObject;
      if (!h3->TestBit(kCanDelete)) {
         for (i = 0; i < fNfill; i++) h3->Fill(fVal[2][i], fVal[1][i], fVal[0][i], fW[i]);
      }
   }
   //__________________________3D scatter plot (3rd variable = col)__
   else if (fAction == 33) {
      TH2 *h2 = (TH2*)fObject;
      TakeEstimate();
      Int_t ncolors  = gStyle->GetNumberOfColors();
      TObjArray *grs = (TObjArray*)h2->GetListOfFunctions()->FindObject("graphs");
      Int_t col;
      TGraph *gr;
      if (!grs) {
         grs = new TObjArray(ncolors);
         grs->SetOwner();
         grs->SetName("graphs");
         h2->GetListOfFunctions()->Add(grs, "P");
         for (col = 0; col < ncolors; col++) {
            gr = new TGraph();
            gr->SetMarkerColor(gStyle->GetColorPalette(col));
            gr->SetMarkerStyle(fTree->GetMarkerStyle());
            gr->SetMarkerSize(fTree->GetMarkerSize());
            grs->AddAt(gr, col);
         }
      }
      h2->SetEntries(fNfill);
      h2->SetMinimum(fVmin[2]);
      h2->SetMaximum(fVmax[2]);
      // Fill the graphs according to the color
      for (i = 0; i < fNfill; i++) {
         col = Int_t(ncolors * ((fVal[2][i] - fVmin[2]) / (fVmax[2] - fVmin[2])));
         if (col < 0) col = 0;
         if (col > ncolors - 1) col = ncolors - 1;
         gr = (TGraph*)grs->UncheckedAt(col);
         if (gr) gr->SetPoint(gr->GetN(), fVal[1][i], fVal[0][i]);
      }
      // Remove potential empty graphs
      for (col = 0; col < ncolors; col++) {
         gr = (TGraph*)grs->At(col);
         if (gr && gr->GetN() <= 0) grs->Remove(gr);
      }
   }
   //__________________________2D Profile Histogram__________________
   else if (fAction == 23) {
      TProfile2D *hp2 = (TProfile2D*)fObject;
      for (i = 0; i < fNfill; i++) hp2->Fill(fVal[2][i], fVal[1][i], fVal[0][i], fW[i]);
   }
   //__________________________4D scatter plot_______________________
   else if (fAction ==  40) {
      TakeEstimate();
      TH3 *h3 = (TH3*)fObject;
      Int_t ncolors  = gStyle->GetNumberOfColors();
      if (ncolors == 0) {
         TColor::InitializeColors();
         ncolors  = gStyle->GetNumberOfColors();
      }
      TObjArray *pms = (TObjArray*)h3->GetListOfFunctions()->FindObject("polymarkers");
      Int_t col;
      TPolyMarker3D *pm3d;
      if (!pms) {
         pms = new TObjArray(ncolors);
         pms->SetOwner();
         pms->SetName("polymarkers");
         h3->GetListOfFunctions()->Add(pms);
         for (col = 0; col < ncolors; col++) {
            pm3d = new TPolyMarker3D();
            pm3d->SetMarkerColor(gStyle->GetColorPalette(col));
            pm3d->SetMarkerStyle(fTree->GetMarkerStyle());
            pm3d->SetMarkerSize(fTree->GetMarkerSize());
            pms->AddAt(pm3d, col);
         }
      }
      h3->SetEntries(fNfill);
      h3->SetMinimum(fVmin[3]);
      h3->SetMaximum(fVmax[3]);
      for (i = 0; i < fNfill; i++) {
         col = Int_t(ncolors * ((fVal[3][i] - fVmin[3]) / (fVmax[3] - fVmin[3])));
         if (col > ncolors-1) col = ncolors-1;
         if (col < 0) col = 0;
         pm3d = (TPolyMarker3D*)pms->UncheckedAt(col);
         pm3d->SetPoint(pm3d->GetLastPoint() + 1, fVal[2][i], fVal[1][i], fVal[0][i]);
      }
   }
   //__________________________Parallel coordinates / candle chart_______________________
   else if (fAction == 6 || fAction == 7) {
      TakeEstimate();
      Bool_t candle = (fAction == 7);
      // Using CINT to avoid a dependency in TParallelCoord
      if (!fOption.Contains("goff"))
         gROOT->ProcessLine(TString::Format("TParallelCoord::BuildParallelCoord((TSelectorDraw*)0x%lx,0x%lx)",
                                (ULong_t)this, (ULong_t)candle));
   } else if (fAction == 8) {
      //gROOT->ProcessLineFast(TString::Format("(new TGL5DDataSet((TTree *)0x%1x))->Draw(\"%s\");", fTree, fOption.Data()));
   }
   //__________________________something else_______________________
   else if (fAction < 0) {
      fAction = -fAction;
      TakeEstimate();
   }

   // Do we need to update screen?
   fSelectedRows += fNfill;
   if (!fTree->GetUpdate()) return;
   if (fSelectedRows > fDraw + fTree->GetUpdate()) {
      if (fDraw) gPad->Modified();
      else       fObject->Draw(fOption.Data());
      gPad->Update();
      fDraw = fSelectedRows;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Estimate limits for 1-D, 2-D or 3-D objects.

void TSelectorDraw::TakeEstimate()
{
   Int_t i;
   Double_t rmin[3], rmax[3];
   Double_t vminOld[4], vmaxOld[4];
   for (i = 0; i < fValSize && i < 4; i++) {
      vminOld[i] = fVmin[i];
      vmaxOld[i] = fVmax[i];
   }
   for (i = 0; i < fValSize; ++i) {
      fVmin[i] = DBL_MAX;
      fVmax[i] = - DBL_MAX;
   }
   //__________________________1-D histogram_______________________
   if (fAction ==  1) {
      TH1 *h1 = (TH1*)fObject;
      if (h1->CanExtendAllAxes()) {
         for (i = 0; i < fNfill; i++) {
            if (fVmin[0] > fVal[0][i]) fVmin[0] = fVal[0][i];
            if (fVmax[0] < fVal[0][i]) fVmax[0] = fVal[0][i];
         }
         THLimitsFinder::GetLimitsFinder()->FindGoodLimits(h1, fVmin[0], fVmax[0]);
      }
      h1->FillN(fNfill, fVal[0], fW);
   //__________________________2-D histogram_______________________
   } else if (fAction ==  2) {
      TH2 *h2 = (TH2*)fObject;
      if (h2->CanExtendAllAxes()) {
         for (i = 0; i < fNfill; i++) {
            if (fVmin[0] > fVal[0][i]) fVmin[0] = fVal[0][i];
            if (fVmax[0] < fVal[0][i]) fVmax[0] = fVal[0][i];
            if (fVmin[1] > fVal[1][i]) fVmin[1] = fVal[1][i];
            if (fVmax[1] < fVal[1][i]) fVmax[1] = fVal[1][i];
         }
         THLimitsFinder::GetLimitsFinder()->FindGoodLimits(h2, fVmin[1], fVmax[1], fVmin[0], fVmax[0]);
      }
      for (i = 0; i < fNfill; i++) h2->Fill(fVal[1][i], fVal[0][i], fW[i]);
   //__________________________Profile histogram_______________________
   } else if (fAction ==  4) {
      TProfile *hp = (TProfile*)fObject;
      if (hp->CanExtendAllAxes()) {
         for (i = 0; i < fNfill; i++) {
            if (fVmin[0] > fVal[0][i]) fVmin[0] = fVal[0][i];
            if (fVmax[0] < fVal[0][i]) fVmax[0] = fVal[0][i];
            if (fVmin[1] > fVal[1][i]) fVmin[1] = fVal[1][i];
            if (fVmax[1] < fVal[1][i]) fVmax[1] = fVal[1][i];
         }
         THLimitsFinder::GetLimitsFinder()->FindGoodLimits(hp, fVmin[1], fVmax[1]);
      }
      hp->FillN(fNfill, fVal[1], fVal[0], fW);
   //__________________________2D scatter plot_______________________
   } else if (fAction == 12) {
      TH2 *h2 = (TH2*)fObject;
      if (h2->CanExtendAllAxes()) {
         for (i = 0; i < fNfill; i++) {
            if (fVmin[0] > fVal[0][i]) fVmin[0] = fVal[0][i];
            if (fVmax[0] < fVal[0][i]) fVmax[0] = fVal[0][i];
            if (fVmin[1] > fVal[1][i]) fVmin[1] = fVal[1][i];
            if (fVmax[1] < fVal[1][i]) fVmax[1] = fVal[1][i];
         }
         THLimitsFinder::GetLimitsFinder()->FindGoodLimits(h2, fVmin[1], fVmax[1], fVmin[0], fVmax[0]);
         // In case the new lower limits of h2 axis are 0, it is better to set them to the minimum of
         // the data set (which should be >0) to avoid data cut when plotting in log scale.
         TAxis *aX = h2->GetXaxis();
         TAxis *aY = h2->GetYaxis();
         Double_t xmin = aX->GetXmin();
         Double_t ymin = aY->GetXmin();
         if (xmin == 0 || ymin == 0) {
            if (aX->GetBinUpEdge(aX->FindFixBin(0.01*aX->GetBinWidth(aX->GetFirst()))) > fVmin[1]) xmin = fVmin[1];
            if (aY->GetBinUpEdge(aY->FindFixBin(0.01*aY->GetBinWidth(aY->GetFirst()))) > fVmin[0]) ymin = fVmin[0];
            h2->SetBins(aX->GetNbins(), xmin, aX->GetXmax(), aY->GetNbins(), ymin, aY->GetXmax());
         }
      }

      if (!strstr(fOption.Data(), "same") && !strstr(fOption.Data(), "goff")) {
         if (!h2->TestBit(kCanDelete)) {
            // case like: T.Draw("y:x>>myhist")
            // we must draw a copy before filling the histogram h2=myhist
            // because h2 will be filled below and we do not want to show
            // the binned scatter-plot, the TGraph being better.
            TH1 *h2c = h2->DrawCopy(fOption.Data(),"");
            if (h2c) h2c->SetStats(kFALSE);
         } else {
            // case like: T.Draw("y:x")
            // h2 is a temporary histogram (htemp). This histogram
            // will be automatically deleted by TPad::Clear
            h2->Draw();
         }
         gPad->Update();
      }
      TGraph *pm = new TGraph(fNfill, fVal[1], fVal[0]);
      pm->SetEditable(kFALSE);
      pm->SetBit(kCanDelete);
      pm->SetMarkerStyle(fTree->GetMarkerStyle());
      pm->SetMarkerColor(fTree->GetMarkerColor());
      pm->SetMarkerSize(fTree->GetMarkerSize());
      pm->SetLineColor(fTree->GetLineColor());
      pm->SetLineWidth(fTree->GetLineWidth());
      pm->SetLineStyle(fTree->GetLineStyle());
      pm->SetFillColor(fTree->GetFillColor());
      pm->SetFillStyle(fTree->GetFillStyle());
      if (!fDraw && !strstr(fOption.Data(),"goff")) {
         if (fOption.Length() == 0 || strcasecmp(fOption.Data(),"same")==0) {
            pm->Draw("p");
         }
         else {
            TString opt = fOption;
            opt.ToLower();
            if (opt.Contains("a")) {
               TString temp(opt);
               temp.ReplaceAll("same","");
               if (temp.Contains("a")) {
                  if (h2->TestBit(kCanDelete)) {
                     // h2 will be deleted, the axis setting is delegated to only
                     // the TGraph.
                     h2 = 0;
                     fObject = pm->GetHistogram();
                  }
               }
            }
            pm->Draw(fOption.Data());
         }
      }
      if (h2 && !h2->TestBit(kCanDelete)) {
         for (i = 0; i < fNfill; i++) h2->Fill(fVal[1][i], fVal[0][i], fW[i]);
      }
   //__________________________3D scatter plot with option col_______________________
   } else if (fAction == 33) {
      TH2 *h2 = (TH2*)fObject;
      Bool_t process2 = kFALSE;
      if (h2->CanExtendAllAxes()) {
         if (vminOld[2] == DBL_MAX)
            process2 = kTRUE;
         for (i = 0; i < fValSize && i < 4; i++) {
            fVmin[i] = vminOld[i];
            fVmax[i] = vmaxOld[i];
         }
         for (i = 0; i < fNfill; i++) {
            if (fVmin[0] > fVal[0][i]) fVmin[0] = fVal[0][i];
            if (fVmax[0] < fVal[0][i]) fVmax[0] = fVal[0][i];
            if (fVmin[1] > fVal[1][i]) fVmin[1] = fVal[1][i];
            if (fVmax[1] < fVal[1][i]) fVmax[1] = fVal[1][i];
            if (process2) {
               if (fVmin[2] > fVal[2][i]) fVmin[2] = fVal[2][i];
               if (fVmax[2] < fVal[2][i]) fVmax[2] = fVal[2][i];
            }
         }
         THLimitsFinder::GetLimitsFinder()->FindGoodLimits(h2, fVmin[1], fVmax[1], fVmin[0], fVmax[0]);
         // In case the new lower limits of h2 axis are 0, it is better to set them to the minimum of
         // the data set (which should be >0) to avoid data cut when plotting in log scale.
         TAxis *aX = h2->GetXaxis();
         TAxis *aY = h2->GetYaxis();
         Double_t xmin = aX->GetXmin();
         Double_t ymin = aY->GetXmin();
         if (xmin == 0 || ymin == 0) {
            if (aX->GetBinUpEdge(aX->FindFixBin(0.01*aX->GetBinWidth(aX->GetFirst()))) > fVmin[1]) xmin = fVmin[1];
            if (aY->GetBinUpEdge(aY->FindFixBin(0.01*aY->GetBinWidth(aY->GetFirst()))) > fVmin[0]) ymin = fVmin[0];
            h2->SetBins(aX->GetNbins(), xmin, aX->GetXmax(), aY->GetNbins(), ymin, aY->GetXmax());
         }
      } else {
         for (i = 0; i < fNfill; i++) {
            if (fVmin[2] > fVal[2][i]) fVmin[2] = fVal[2][i];
            if (fVmax[2] < fVal[2][i]) fVmax[2] = fVal[2][i];
         }
      }
      //__________________________3D scatter plot_______________________
   } else if (fAction == 3 || fAction == 13) {
      TH3 *h3 = (TH3*)fObject;
      if (h3->CanExtendAllAxes()) {
         for (i = 0; i < fNfill; i++) {
            if (fVmin[0] > fVal[0][i]) fVmin[0] = fVal[0][i];
            if (fVmax[0] < fVal[0][i]) fVmax[0] = fVal[0][i];
            if (fVmin[1] > fVal[1][i]) fVmin[1] = fVal[1][i];
            if (fVmax[1] < fVal[1][i]) fVmax[1] = fVal[1][i];
            if (fVmin[2] > fVal[2][i]) fVmin[2] = fVal[2][i];
            if (fVmax[2] < fVal[2][i]) fVmax[2] = fVal[2][i];
         }
         THLimitsFinder::GetLimitsFinder()->FindGoodLimits(h3, fVmin[2], fVmax[2], fVmin[1], fVmax[1], fVmin[0], fVmax[0]);
      }
      if (fAction == 3) {
         for (i = 0; i < fNfill; i++) h3->Fill(fVal[2][i], fVal[1][i], fVal[0][i], fW[i]);
         return;
      }
      if (!strstr(fOption.Data(), "same") && !strstr(fOption.Data(), "goff")) {
         if (!h3->TestBit(kCanDelete)) {
            // case like: T.Draw("y:x>>myhist")
            // we must draw a copy before filling the histogram h3=myhist
            // because h3 will be filled below and we do not want to show
            // the binned scatter-plot, the TGraph being better.
            TH1 *h3c = h3->DrawCopy(fOption.Data(),"");
            if (h3c) h3c->SetStats(kFALSE);
         } else {
            // case like: T.Draw("y:x")
            // h3 is a temporary histogram (htemp). This histogram
            // will be automatically deleted by TPad::Clear
            h3->Draw(fOption.Data());
         }
         gPad->Update();
      } else {
         rmin[0] = fVmin[2]; rmin[1] = fVmin[1]; rmin[2] = fVmin[0];
         rmax[0] = fVmax[2]; rmax[1] = fVmax[1]; rmax[2] = fVmax[0];
         gPad->Clear();
         gPad->Range(-1, -1, 1, 1);
         TView::CreateView(1, rmin, rmax);
      }
      TPolyMarker3D *pm3d = new TPolyMarker3D(fNfill);
      pm3d->SetMarkerStyle(fTree->GetMarkerStyle());
      pm3d->SetMarkerColor(fTree->GetMarkerColor());
      pm3d->SetMarkerSize(fTree->GetMarkerSize());
      for (i = 0; i < fNfill; i++) {
         pm3d->SetPoint(i, fVal[2][i], fVal[1][i], fVal[0][i]);
      }
      if (!fDraw && !strstr(fOption.Data(), "goff")) pm3d->Draw();
      if (!h3->TestBit(kCanDelete)) {
         for (i = 0; i < fNfill; i++) h3->Fill(fVal[2][i], fVal[1][i], fVal[0][i], fW[i]);
      }

   //__________________________2D Profile Histogram__________________
   } else if (fAction == 23) {
      TProfile2D *hp = (TProfile2D*)fObject;
      if (hp->CanExtendAllAxes()) {
         for (i = 0; i < fNfill; i++) {
            if (fVmin[0] > fVal[0][i]) fVmin[0] = fVal[0][i];
            if (fVmax[0] < fVal[0][i]) fVmax[0] = fVal[0][i];
            if (fVmin[1] > fVal[1][i]) fVmin[1] = fVal[1][i];
            if (fVmax[1] < fVal[1][i]) fVmax[1] = fVal[1][i];
            if (fVmin[2] > fVal[2][i]) fVmin[2] = fVal[2][i];
            if (fVmax[2] < fVal[2][i]) fVmax[2] = fVal[2][i];
         }
         THLimitsFinder::GetLimitsFinder()->FindGoodLimits(hp, fVmin[2], fVmax[2], fVmin[1], fVmax[1]);
      }
      for (i = 0; i < fNfill; i++) hp->Fill(fVal[2][i], fVal[1][i], fVal[0][i], fW[i]);
   //__________________________4D scatter plot_______________________
   } else if (fAction == 40) {
      TH3 *h3 = (TH3*)fObject;
      if (h3->CanExtendAllAxes()) {
         for (i = 0; i < fValSize && i < 4; i++) {
            fVmin[i] = vminOld[i];
            fVmax[i] = vmaxOld[i];
         }
         for (i = 0; i < fNfill; i++) {
            if (fVmin[0] > fVal[0][i]) fVmin[0] = fVal[0][i];
            if (fVmax[0] < fVal[0][i]) fVmax[0] = fVal[0][i];
            if (fVmin[1] > fVal[1][i]) fVmin[1] = fVal[1][i];
            if (fVmax[1] < fVal[1][i]) fVmax[1] = fVal[1][i];
            if (fVmin[2] > fVal[2][i]) fVmin[2] = fVal[2][i];
            if (fVmax[2] < fVal[2][i]) fVmax[2] = fVal[2][i];
            if (fVmin[3] > fVal[3][i]) fVmin[3] = fVal[3][i];
            if (fVmax[3] < fVal[3][i]) fVmax[3] = fVal[3][i];
         }
         THLimitsFinder::GetLimitsFinder()->FindGoodLimits(h3, fVmin[2], fVmax[2], fVmin[1], fVmax[1], fVmin[0], fVmax[0]);
      } else {
         for (i = 0; i < fNfill; i++) {
            if (fVmin[3] > fVal[3][i]) fVmin[3] = fVal[3][i];
            if (fVmax[3] < fVal[3][i]) fVmax[3] = fVal[3][i];
         }
      }
   }
   //__________________________Parallel coordinates plot / candle chart_______________________
   else if (fAction == 6 || fAction == 7) {
      for (i = 0; i < fDimension; ++i) {
         for (Long64_t entry = 0; entry < fNfill; entry++) {
            if (fVmin[i] > fVal[i][entry]) fVmin[i] = fVal[i][entry];
            if (fVmax[i] < fVal[i][entry]) fVmax[i] = fVal[i][entry];
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called at the end of a loop on a TTree.

void TSelectorDraw::Terminate()
{
   if (fNfill) TakeAction();

   if ((fSelectedRows == 0) && (TestBit(kCustomHistogram) == 0)) fDraw = 1; // do not draw

   SetStatus(fSelectedRows);
}
