// @(#)root/hist:$Id$
// Author: Christophe.Delaere@cern.ch   21/08/2002

///////////////////////////////////////////////////////////////////////////
//
// TLimit
//
// Class to compute 95% CL limits
//
// adapted from the mclimit code from Tom Junk (CLs method)
// see http://root.cern/root/doc/TomJunk.pdf
// see http://cern.ch/thomasj/searchlimits/ecl.html
// see: Tom Junk,NIM A434, p. 435-443, 1999
//
// see also the following interesting references:
//  Alex Read, "Presentation of search results: the CLs technique"
//  Journal of Physics G: Nucl. Part. Phys. 28 2693-2704 (2002).
//  http://www.iop.org/EJ/abstract/0954-3899/28/10/313/
//
// A nice article is also available in the CERN yellow report with the proceeding
// of the 2000 CERN workshop on confidence intervals.
//
//  Alex Read, "Modified Frequentist Analysis of Search Results (The CLs Method)"
//  CERN 2000-005 (30 May 2000)
//
// see note about: "Should I use TRolke, TFeldmanCousins, TLimit?"
//  in the TRolke class description.
//
///////////////////////////////////////////////////////////////////////////

/** \class TLimit
    \ingroup Hist
 Algorithm to compute 95% C.L. limits using the Likelihood ratio
 semi-bayesian method.
 
 Implemented by C. Delaere from the mclimit code written by Tom Junk [HEP-EX/9902006]. 
 See [http://cern.ch/thomasj/searchlimits/ecl.html](http://cern.ch/thomasj/searchlimits/ecl.html) for more details.

 It takes signal, background and data histograms wrapped in a
 TLimitDataSource as input and runs a set of Monte Carlo experiments in
 order to compute the limits. If needed, inputs are fluctuated according
 to systematics. The output is a TConfidenceLevel.

 The class TLimitDataSource takes the signal, background and data histograms as well as different
 systematics sources to form the TLimit input.

 The class TConfidenceLevel represents the final result of the TLimit algorithm. It is created just after the
 time-consuming part and can be stored in a TFile for further processing.
 It contains light methods to return CLs, CLb and other interesting
 quantities.

 The actual algorithm...

 From an input (TLimitDataSource) it produces an output TConfidenceLevel.
 For this, nmc Monte Carlo experiments are performed.
 As usual, the larger this number, the longer the compute time,
 but the better the result.

 Supposing that there is a plotfile.root file containing 3 histograms
 (signal, background and data), you can imagine doing things like:

~~~{.cpp}
 TFile* infile=new TFile("plotfile.root","READ");
 infile->cd();
 TH1* sh=(TH1*)infile->Get("signal");
 TH1* bh=(TH1*)infile->Get("background");
 TH1* dh=(TH1*)infile->Get("data");
 TLimitDataSource* mydatasource = new TLimitDataSource(sh,bh,dh);
 TConfidenceLevel *myconfidence = TLimit::ComputeLimit(mydatasource,50000);
 std::cout << "  CLs    : " << myconfidence->CLs()  << std::endl;
 std::cout << "  CLsb   : " << myconfidence->CLsb() << std::endl;
 std::cout << "  CLb    : " << myconfidence->CLb()  << std::endl;
 std::cout << "< CLs >  : " << myconfidence->GetExpectedCLs_b()  << std::endl;
 std::cout << "< CLsb > : " << myconfidence->GetExpectedCLsb_b() << std::endl;
 std::cout << "< CLb >  : " << myconfidence->GetExpectedCLb_b()  << std::endl;
 delete myconfidence;
 delete mydatasource;
 infile->Close();
~~~
 More information can still be found on [this page](http://cern.ch/aleph-proj-alphapp/doc/tlimit.html)
 */

#include "TLimit.h"
#include "TArrayD.h"
#include "TVectorD.h"
#include "TOrdCollection.h"
#include "TConfidenceLevel.h"
#include "TLimitDataSource.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TIterator.h"
#include "TObjString.h"
#include "TClassTable.h"
#include "Riostream.h"

ClassImp(TLimit);

TArrayD *TLimit::fgTable = new TArrayD(0);
TOrdCollection *TLimit::fgSystNames = new TOrdCollection();

TConfidenceLevel *TLimit::ComputeLimit(TLimitDataSource * data,
                                       Int_t nmc, bool stat,
                                       TRandom * generator)
{
   // The final object returned...
   TConfidenceLevel *result = new TConfidenceLevel(nmc);
   // The random generator used...
   TRandom *myrandom = generator ? generator : new TRandom3;
   // Compute some total quantities on all the channels
   Int_t maxbins = 0;
   Double_t nsig = 0;
   Double_t nbg  = 0;
   Int_t ncand   = 0;
   Int_t i;
   for (i = 0; i <= data->GetSignal()->GetLast(); i++) {
      maxbins = (((TH1 *) (data->GetSignal()->At(i)))->GetNbinsX() + 2) > maxbins ?
                 (((TH1 *) (data->GetSignal()->At(i)))->GetNbinsX() + 2) : maxbins;
      nsig   +=  ((TH1 *) (data->GetSignal()->At(i)))->Integral();
      nbg    +=  ((TH1 *) (data->GetBackground()->At(i)))->Integral();
      ncand  +=  (Int_t) ((TH1 *) (data->GetCandidates()->At(i)))->Integral();
   }
   result->SetBtot(nbg);
   result->SetStot(nsig);
   result->SetDtot(ncand);
   Double_t buffer = 0;
   fgTable->Set(maxbins * (data->GetSignal()->GetLast() + 1));
   for (Int_t channel = 0; channel <= data->GetSignal()->GetLast(); channel++)
      for (Int_t bin = 0;
           bin <= ((TH1 *) (data->GetSignal()->At(channel)))->GetNbinsX()+1;
           bin++) {
         Double_t s = (Double_t) ((TH1 *) (data->GetSignal()->At(channel)))->GetBinContent(bin);
         Double_t b = (Double_t) ((TH1 *) (data->GetBackground()->At(channel)))->GetBinContent(bin);
         Double_t d = (Double_t) ((TH1 *) (data->GetCandidates()->At(channel)))->GetBinContent(bin);
         // Compute the value of the "-2lnQ" for the actual data
         if ((b == 0) && (s > 0)) {
            std::cout << "WARNING: Ignoring bin " << bin << " of channel "
                 << channel << " which has s=" << s << " but b=" << b << std::endl;
            std::cout << "         Maybe the MC statistic has to be improved..." << std::endl;
         }
         if ((s > 0) && (b > 0))
            buffer += LogLikelihood(s, b, b, d);
         // precompute the log(1+s/b)'s in an array to speed up computation
         // background-free bins are set to have a maximum t.s. value
         // for protection (corresponding to s/b of about 5E8)
         if ((s > 0) && (b > 0))
            fgTable->AddAt(LogLikelihood(s, b, b, 1), (channel * maxbins) + bin);
         else if ((s > 0) && (b == 0))
            fgTable->AddAt(20, (channel * maxbins) + bin);
      }
   result->SetTSD(buffer);
   // accumulate MC experiments.  Hold the test statistic function fixed, but
   // fluctuate s and b within syst. errors for computing probabilities of
   // having that outcome.  (Alex Read's prescription -- errors are on the ensemble,
   // not on the observed test statistic.  This technique does not split outcomes.)
   // keep the tstats as sum log(1+s/b). convert to -2lnQ when preparing the results
   // (reason -- like to keep the < signs right)
   Double_t *tss = new Double_t[nmc];
   Double_t *tsb = new Double_t[nmc];
   Double_t *lrs = new Double_t[nmc];
   Double_t *lrb = new Double_t[nmc];
   TLimitDataSource* tmp_src = (TLimitDataSource*)(data->Clone());
   TLimitDataSource* tmp_src2 = (TLimitDataSource*)(data->Clone());
   for (i = 0; i < nmc; i++) {
      tss[i] = 0;
      tsb[i] = 0;
      lrs[i] = 0;
      lrb[i] = 0;
      // fluctuate signal and background
      TLimitDataSource*  fluctuated = Fluctuate(data, tmp_src, !i, myrandom, stat) ? tmp_src : data;
      // make an independent set of fluctuations for reweighting pseudoexperiments
      // it turns out using the same fluctuated ones for numerator and denominator
      // is biased.
      TLimitDataSource*  fluctuated2 = Fluctuate(data, tmp_src2, false, myrandom, stat) ? tmp_src2 : data;

      for (Int_t channel = 0;
           channel <= fluctuated->GetSignal()->GetLast(); channel++) {
         for (Int_t bin = 0;
              bin <=((TH1 *) (fluctuated->GetSignal()->At(channel)))->GetNbinsX()+1;
              bin++) {
            if ((Double_t) ((TH1 *) (fluctuated->GetSignal()->At(channel)))->GetBinContent(bin) != 0) {
               // s+b hypothesis
               Double_t rate = (Double_t) ((TH1 *) (fluctuated->GetSignal()->At(channel)))->GetBinContent(bin) +
                               (Double_t) ((TH1 *) (fluctuated->GetBackground()->At(channel)))->GetBinContent(bin);
               Double_t rand = myrandom->Poisson(rate);
               tss[i] += rand * fgTable->At((channel * maxbins) + bin);
               Double_t s = (Double_t) ((TH1 *) (fluctuated->GetSignal()->At(channel)))->GetBinContent(bin);
               Double_t s2= (Double_t) ((TH1 *) (fluctuated2->GetSignal()->At(channel)))->GetBinContent(bin);
               Double_t b = (Double_t) ((TH1 *) (fluctuated->GetBackground()->At(channel)))->GetBinContent(bin);
               Double_t b2= (Double_t) ((TH1 *) (fluctuated2->GetBackground()->At(channel)))->GetBinContent(bin);
               if ((s > 0) && (b2 > 0))
                  lrs[i] += LogLikelihood(s, b, b2, rand) - s - b + b2;
               else if ((s > 0) && (b2 == 0))
                  lrs[i] += 20 * rand - s;
               // b hypothesis
               rate = (Double_t) ((TH1 *) (fluctuated->GetBackground()->At(channel)))->GetBinContent(bin);
               rand = myrandom->Poisson(rate);
               tsb[i] += rand * fgTable->At((channel * maxbins) + bin);
               if ((s2 > 0) && (b > 0))
                  lrb[i] += LogLikelihood(s2, b2, b, rand) - s2 - b2 + b;
               else if ((s > 0) && (b == 0))
                  lrb[i] += 20 * rand - s;
            }
         }
      }
      lrs[i] = TMath::Exp(lrs[i] < 710 ? lrs[i] : 710);
      lrb[i] = TMath::Exp(lrb[i] < 710 ? lrb[i] : 710);
   }
   delete tmp_src;
   delete tmp_src2;
   // lrs and lrb are the LR's (no logs) = prob(s+b)/prob(b) for
   // that choice of s and b within syst. errors in the ensemble.  These are
   // the MC experiment weights for relating the s+b and b PDF's of the unsmeared
   // test statistic (in which cas one can use another test statistic if one likes).

   // Now produce the output object.
   // The final quantities are computed on-demand form the arrays tss, tsb, lrs and lrb.
   result->SetTSS(tss);
   result->SetTSB(tsb);
   result->SetLRS(lrs);
   result->SetLRB(lrb);
   if (!generator)
      delete myrandom;
   return result;
}

bool TLimit::Fluctuate(TLimitDataSource * input, TLimitDataSource * output,
                       bool init, TRandom * generator, bool stat)
{
   // initialisation: create a sorted list of all the names of systematics
   if (init) {
      // create a "map" with the systematics names
      TIterator *errornames = input->GetErrorNames()->MakeIterator();
      TObjArray *listofnames = 0;
      delete fgSystNames;
      fgSystNames = new TOrdCollection();
      while ((listofnames = ((TObjArray *) errornames->Next()))) {
         TObjString *name = 0;
         TIterator *loniter = listofnames->MakeIterator();
         while ((name = (TObjString *) (loniter->Next())))
            if ((fgSystNames->IndexOf(name)) < 0)
               fgSystNames->AddLast(name);
      }
      fgSystNames->Sort();
   }
   // if the output is not given, create it from the input
   if (!output)
   output = (TLimitDataSource*)(input->Clone());
   // if there are no systematics, just returns the input as "fluctuated" output
   if ((fgSystNames->GetSize() <= 0)&&(!stat)) {
      return 0;
   }
   // if there are only stat, just fluctuate stats.
   if (fgSystNames->GetSize() <= 0) {
      output->SetOwner();
      for (Int_t channel = 0; channel <= input->GetSignal()->GetLast(); channel++) {
         TH1 *newsignal = (TH1*)(output->GetSignal()->At(channel));
         TH1 *oldsignal = (TH1*)(input->GetSignal()->At(channel));
         if(stat)
            for(int i=1; i<=newsignal->GetNbinsX(); i++) {
               newsignal->SetBinContent(i,oldsignal->GetBinContent(i)+generator->Gaus(0,oldsignal->GetBinError(i)));
            }
         newsignal->SetDirectory(0);
         TH1 *newbackground = (TH1*)(output->GetBackground()->At(channel));
         TH1 *oldbackground = (TH1*)(input->GetBackground()->At(channel));
         if(stat)
            for(int i=1; i<=newbackground->GetNbinsX(); i++)
               newbackground->SetBinContent(i,oldbackground->GetBinContent(i)+generator->Gaus(0,oldbackground->GetBinError(i)));
         newbackground->SetDirectory(0);
      }
      return 1;
   }
   // Find a choice for the random variation and
   // re-toss all random numbers if any background or signal
   // goes negative.  (background = 0 is bad too, so put a little protection
   // around it -- must have at least 10% of the bg estimate).
   Bool_t retoss   = kTRUE;
   Double_t *serrf = new Double_t[(input->GetSignal()->GetLast()) + 1];
   Double_t *berrf = new Double_t[(input->GetSignal()->GetLast()) + 1];
   do {
      Double_t *toss = new Double_t[fgSystNames->GetSize()];
      for (Int_t i = 0; i < fgSystNames->GetSize(); i++)
         toss[i] = generator->Gaus(0, 1);
      retoss = kFALSE;
      for (Int_t channel = 0;
           channel <= input->GetSignal()->GetLast();
           channel++) {
         serrf[channel] = 0;
         berrf[channel] = 0;
         for (Int_t bin = 0;
              bin <((TVectorD *) (input->GetErrorOnSignal()->At(channel)))->GetNrows();
              bin++) {
            serrf[channel] += ((TVectorD *) (input->GetErrorOnSignal()->At(channel)))->operator[](bin) *
                toss[fgSystNames->BinarySearch((TObjString*) (((TObjArray *) (input->GetErrorNames()->At(channel)))->At(bin)))];
            berrf[channel] += ((TVectorD *) (input->GetErrorOnBackground()->At(channel)))->operator[](bin) *
                toss[fgSystNames->BinarySearch((TObjString*) (((TObjArray *) (input->GetErrorNames()->At(channel)))->At(bin)))];
         }
         if ((serrf[channel] < -1.0) || (berrf[channel] < -0.9)) {
            retoss = kTRUE;
            continue;
         }
      }
      delete[]toss;
   } while (retoss);
   // adjust the fluctuated signal and background counts with a legal set
   // of random fluctuations above.
   output->SetOwner();
   for (Int_t channel = 0; channel <= input->GetSignal()->GetLast();
        channel++) {
      TH1 *newsignal = (TH1*)(output->GetSignal()->At(channel));
      TH1 *oldsignal = (TH1*)(input->GetSignal()->At(channel));
      if(stat)
         for(int i=1; i<=newsignal->GetNbinsX(); i++)
            newsignal->SetBinContent(i,oldsignal->GetBinContent(i)+generator->Gaus(0,oldsignal->GetBinError(i)));
      else
         for(int i=1; i<=newsignal->GetNbinsX(); i++)
            newsignal->SetBinContent(i,oldsignal->GetBinContent(i));
      newsignal->Scale(1 + serrf[channel]);
      newsignal->SetDirectory(0);
      TH1 *newbackground = (TH1*)(output->GetBackground()->At(channel));
      TH1 *oldbackground = (TH1*)(input->GetBackground()->At(channel));
      if(stat)
         for(int i=1; i<=newbackground->GetNbinsX(); i++)
            newbackground->SetBinContent(i,oldbackground->GetBinContent(i)+generator->Gaus(0,oldbackground->GetBinError(i)));
      else
         for(int i=1; i<=newbackground->GetNbinsX(); i++)
            newbackground->SetBinContent(i,oldbackground->GetBinContent(i));
      newbackground->Scale(1 + berrf[channel]);
      newbackground->SetDirectory(0);
   }
   delete[] serrf;
   delete[] berrf;
   return 1;
}

TConfidenceLevel *TLimit::ComputeLimit(TH1* s, TH1* b, TH1* d,
                                       Int_t nmc, bool stat,
                                       TRandom * generator)
{
   // Compute limit.

   TLimitDataSource* lds = new TLimitDataSource(s,b,d);
   TConfidenceLevel* out = ComputeLimit(lds,nmc,stat,generator);
   delete lds;
   return out;
}

TConfidenceLevel *TLimit::ComputeLimit(TH1* s, TH1* b, TH1* d,
                                       TVectorD* se, TVectorD* be, TObjArray* l,
                                       Int_t nmc, bool stat,
                                       TRandom * generator)
{
   // Compute limit.

   TLimitDataSource* lds = new TLimitDataSource(s,b,d,se,be,l);
   TConfidenceLevel* out = ComputeLimit(lds,nmc,stat,generator);
   delete lds;
   return out;
}

TConfidenceLevel *TLimit::ComputeLimit(Double_t s, Double_t b, Int_t d,
                                       Int_t nmc,
                                       bool stat,
                                       TRandom * generator)
{
   // Compute limit.

   TH1D* sh = new TH1D("__sh","__sh",1,0,2);
   sh->Fill(1,s);
   TH1D* bh = new TH1D("__bh","__bh",1,0,2);
   bh->Fill(1,b);
   TH1D* dh = new TH1D("__dh","__dh",1,0,2);
   dh->Fill(1,d);
   TLimitDataSource* lds = new TLimitDataSource(sh,bh,dh);
   TConfidenceLevel* out = ComputeLimit(lds,nmc,stat,generator);
   delete lds;
   delete sh;
   delete bh;
   delete dh;
   return out;
}

TConfidenceLevel *TLimit::ComputeLimit(Double_t s, Double_t b, Int_t d,
                                       TVectorD* se, TVectorD* be, TObjArray* l,
                                       Int_t nmc,
                                       bool stat,
                                       TRandom * generator)
{
   // Compute limit.

   TH1D* sh = new TH1D("__sh","__sh",1,0,2);
   sh->Fill(1,s);
   TH1D* bh = new TH1D("__bh","__bh",1,0,2);
   bh->Fill(1,b);
   TH1D* dh = new TH1D("__dh","__dh",1,0,2);
   dh->Fill(1,d);
   TLimitDataSource* lds = new TLimitDataSource(sh,bh,dh,se,be,l);
   TConfidenceLevel* out = ComputeLimit(lds,nmc,stat,generator);
   delete lds;
   delete sh;
   delete bh;
   delete dh;
   return out;
}

Double_t TLimit::LogLikelihood(Double_t s, Double_t b, Double_t b2, Double_t d)
{
   // Compute LogLikelihood (static function)

   return d*TMath::Log((s+b)/b2);
}
